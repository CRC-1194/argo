#!/usr/bin/env python

from argparse import ArgumentParser
from subprocess import run
import os
import sys
import time


app_description = """ Initialize all test cases that match a given pattern.
    Initialization consists of the following steps:
        1. Create a mesh using the prescribed meshing application.\n
        2. Initialize required fields.\n
        3. Optional: perform domain decomposition for parallel simulation.\n
"""

def find_case_directories(pattern):
    """
    Find all directories whose name starts with 'pattern'.
    """
    case_dirs = [case_dir for case_dir in os.listdir(os.curdir) if os.path.isdir(case_dir) and \
                  case_dir.startswith(pattern)] 

    case_dirs.sort()

    return case_dirs

def extract_job_id(sbatch_output):
    # The job id is the last word in the first line of the
    # output given by the sbatch command
    first_line = sbtach_output.split('\n')[0]
    job_id = first_line.rsplit(' ', 1)[1]

    return job_id


def submit_slurm_job(command, log_name="", depends_on=None):
    account = "special00005"
    ntasks = "1"
    mem_per_cpu = "5000"
    time_limit = "00:15:00" # Read: hh:mm:ss, h-> hours, m-> minutes, s-> seconds
    wait_for_job = None

    if depends_on:
        wait_for_job = "--dependency=afterok:" + depends_on
    slurm_command = ("srun --account %s --ntasks %s --mem-per-cpu %s "
                     "--time %s %s "
                     "--job-name=%s --output=log.%s %s >/dev/null 2>&1 &" % 
                     (account, ntasks, mem_per_cpu,
                      time_limit, wait_for_job,
                      log_name, log_name, command))

    sbatch_out = run(slurm_command, shell=True, capture_output=True)
    job_id = extract_job_id(sbatch_out.stdout)

    # Force wait after submitting job to avoid overloading the SLURM manager
    time.sleep(1)

    # Return job id to make initialization job dependent on the meshing job.
    # That ensures that field initialization is done after mesh creation.
    return job_id

def execute_init_step(command, use_slurm, log_name="", depends_on=None):
    job_id = None
    if use_slurm:
        job_id = submit_slurm_job(command, log_name=log_name, depends_on=depends_on)
    else:
        run(command, shell=True)

    return job_id


#---- Command line arguments ----------------------------------------------
parser = ArgumentParser(description=app_description)

parser.add_argument("directory_pattern",
                    help="Pattern matching those directory names which are to be initialized, e.g. myAwesomeStudy_000")
parser.add_argument("-m", "--meshing-application",
                    help="Use this application for mesh creation, e.g. blockMesh. Use 'none' to skip meshing.",
                    required=True,
                    dest="mesh_app")
parser.add_argument("-f", "--field-init-script",
                    help="Name of the script to use for field initialization.\nDefault: None",
                    default="",
                    dest="init_script")
parser.add_argument("-j", "--job-mode",
                    help="Submit meshing and preprocessing as SLURM jobs",
                    action="store_true",
                    default=False,
                    dest="job_mode")
parser.add_argument("-par", "--parallel",
                    help="Use domain decomposition to setup cases for parallel execution.",
                    action="store_true",
                    default=False,
                    dest="parallel")

#--------------------------------------------------------------------------

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    if not args.init_script:
        print("No init script specified via --field-init-script. Not initializing fields.")


    case_dirs = find_case_directories(args.directory_pattern)
    
    case_count = 0
    for case in case_dirs:
        pwd = os.getcwd()
        os.chdir(case)

        case_count = case_count + 1
        print("(%s/%s) Initializing %s ..." % (str(case_count), str(len(case_dirs)), case))

        # Meshing. Can be skipped. No dependency on a previous job.
        mesh_job_id = None
        if args.mesh_app != "none":
            mesh_job_id = execute_init_step(args.mesh_app, args.job_mode, log_name=args.mesh_app)
        else:
            print("Warning: skipping mesh generation. Make sure your init script generates a mesh.")

        # Field initialization
        field_job_id = None
        if args.init_script:
            # Note: the template replacement process invoked in 'argo-create-parameter-study.py'
            #       removes the She-Bang from the init script. So we need to use the file extension
            #       to detect the suitable interpreter (TT)
            script_env = ""
            if args.init_script.endswith('.py'):
                script_env = "python3"
            elif args.init_script.endswith('.sh'):
                script_env = "bash"
            else:
                sys.exit("Error: unsupported script format", args.init_script.rsplit('.')[0], ". Use either .sh for BASH or .py for Python3.")
            execute_init_step(script_env + " " + args.init_script, args.job_mode,
                              log_name=args.init_script, depends_on=mesh_job_id)

        # Domain decomposition
        if args.parallel:
            execute_init_step("decomposePar -force", args.job_mode, log_name="decomposePar",
                                depends_on=field_job_id)

        os.chdir(pwd)
