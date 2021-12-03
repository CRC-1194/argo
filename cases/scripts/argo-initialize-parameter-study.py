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
    In each test case a script named 'argo-initialize-case.sh' is created\n
    that contains all requested steps.
"""


slurm_job_config = (
    "#SBATCH --account special00005\n"
    "#SBATCH --ntasks 1\n"
    "#SBATCH --mem-per-cpu 4000\n"
    "#SBATCH --time 00:15:00\n"
    "#SBATCH --job-name argo-case-init\n"
    "#SBATCH --output=log.argo-initialze-case\n"
    )


def find_case_directories(pattern):
    """
    Find all directories whose name starts with 'pattern'.
    """
    case_dirs = [case_dir for case_dir in os.listdir(os.curdir) if os.path.isdir(case_dir) and \
                  case_dir.startswith(pattern)] 
    case_dirs.sort()

    return case_dirs


def meshing_command(mesher):
    return mesher + "\n"


def field_init_command(field_init_script):
    """
    Assemble command for field initialization based on file extension.
    """
    if args.init_script.endswith('.py'):
        return "python3 " + field_init_script + "\n"
    elif args.init_script.endswith('.sh'):
        return "bash " + field_init_script + "\n"
    else:
        sys.exit("Error: unsupported script format", field_init_script.rsplit('.')[0], ". Use either .sh for BASH or .py for Python3.")
  

def decompose_command():
    return "decomposePar -force\n"


def create_init_script(mesher, field_init_script, parallel=False):
    """
    Combine all requested initilization steps into a single script.

    The three initialization steps are:
    - meshing
    - setting fields
    - domain decomposition
    where each step can be omitted.
    """
    script_name = "argo-initialze-case.sh"

    script = "#!/usr/bin/bash\n" + slurm_job_config
    if mesher != "none":
        script = script + meshing_command(mesher)
    if field_init_script:
        script = script + field_init_command(field_init_script)
    if parallel:
        script = script + decompose_command()

    with open(script_name, "w") as script_file:
        script_file.write(script)

    return script_name


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

    if args.mesh_app == "none":
        print("Warning: skipping mesh generation. Make sure your init script generates a mesh.")
    if not args.init_script:
        print("No init script specified via --field-init-script. Not initializing fields.")

    case_dirs = find_case_directories(args.directory_pattern)
    
    case_count = 0
    for case in case_dirs:
        pwd = os.getcwd()
        os.chdir(case)

        case_count = case_count + 1
        print("(%s/%s) Initializing %s ..." % (str(case_count), str(len(case_dirs)), case))

        script_name = create_init_script(args.mesh_app, args.init_script, args.parallel)

        if args.job_mode:
            run("sbatch " + script_name, shell=True)
            # Forced wait time to prevent overloading the SLURM job manager
            time.sleep(1)
        else:
            run("nohup bash " + script_name + " > log.argo-initialze-case", shell=True)

        os.chdir(pwd)
