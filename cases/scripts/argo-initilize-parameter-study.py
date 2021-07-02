#!/usr/bin/env python3

from argparse import ArgumentParser
from subprocess import run
import os
import sys


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

def submit_slurm_job(command, log_name=""):
    # Default SLURM options
    slurm_opts = {}
    # Required arguments
    slurm_opts["--partition"] = "test30m"
    slurm_opts["--account"] = "special000005"
    slurm_opts["--mem-per-cpu"] = "5000"
    slurm_opts["--time"] = "00:15:00"
    # Optional arguments
    slurm_opts["--job-name"] = log_name
    slurm_opts["--output"] = "log." + log_name

    slurm_command = ["srun"]
    
    for key,value in slurm_opts.items():
        slurm_command.append(key)
        slurm_command.append(value)

    for entry in command:
        slurm_command.append(entry)

    run(slurm_command)

    # Force wait after submitting job to avoid overloading the SLURM manager
    time.sleep(1)

def execute_init_step(command, use_slurm, log_name=""):
    if use_slurm:
        submit_slurm_job(command, log_name=log_name)
    else:
        run(command)


#---- Command line arguments ----------------------------------------------
parser = ArgumentParser(description=app_description)

parser.add_argument("directory_pattern",
                    help="Pattern matching those directory names which are to be initialized.")
parser.add_argument("-m", "--meshing-application",
                    help="Use this application for mesh creation, e.g. blockMesh.\nDefault: blockMesh",
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

    for case in case_dirs:
        pwd = os.getcwd()
        os.chdir(case)

        print("Initializing", case, "...")

        # Meshing
        execute_init_step([args.mesh_app], args.job_mode, log_name=args.mesh_app)

        # Field initialization
        if args.init_script:
            execute_init_step(["sh", "./"+args.init_script], args.job_mode, log_name=args.init_script)

        # Domain decomposition
        if args.parallel:
            execute_init_step(["decomposePar", "-force"], args.job_mode, log_name="decomposePar")

        os.chdir(pwd)
