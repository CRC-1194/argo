#!/usr/bin/env python

from argparse import ArgumentParser
import datetime 
import os
import shlex        # Parse shell arguments for solver, see
                    # https://docs.python.org/3/library/subprocess.html#subprocess.Popen
import subprocess
import sys
import time

import testReportCore as trc

app_description="""Execute a parameter study with a user-specified solver.
The script can execute multiple cases of a study in parallel using the option
'--num-cases-parallel / -nc'.
If you also want to run each case in parallel, set the number of MPI processes with
the option '--num-mpi-processes / -np'. 
The total number of processes created is nc*np.
"""

# The solution to control the number parallel running processes is taken / adapted from
# https://stackoverflow.com/questions/18123325/always-run-a-constant-number-of-subprocesses-in-parallel

nextVariant = 0
maxNumProcesses = 1
nVariants = 1
solver = ""
Processes = []
return_codes = []
proc_to_variant = {}

def start_new(variantList, n_mpi, solver_args):
    """ Start a new subprocess if there is work to do """
    global nextVariant
    global Processes
    global proc_to_variant
    global return_codes

    args = shlex.split(solver_args)

    solver_command = [solver] + args + ["-case", variantList[nextVariant]]
    mpi_prefix = ["mpirun", "-n", str(n_mpi)]

    with open(variantList[nextVariant] + "/log." + solver, 'w') as logfile:
        if n_mpi > 1:
            proc = subprocess.Popen(mpi_prefix + solver_command + ["-parallel"],
                                    stdout=logfile, stderr=logfile)
        else:
            proc = subprocess.Popen(solver_command, stdout=logfile, stderr=logfile)

    print ("Started variant", variantList[nextVariant].rsplit('/',1)[1],"; process ID", proc.pid)
    proc_to_variant[proc.pid] = nextVariant

    nextVariant += 1
    Processes.append(proc)

def check_running(variantList, n_mpi, solver_args):
   """ Check any running processes and start new ones if there are spare slots."""
   global Processes
   global nextVariant
   global proc_to_variant
   global return_codes

   for p in range(len(Processes)-1,-1,-1): # Check the processes in reverse order
      if Processes[p].poll() is not None: # If the process hasn't finished will return None
         variant = proc_to_variant[Processes[p].pid]
         return_codes[variant] =  Processes[p].returncode
         print("Finished variant", variantList[variant].rsplit('/',1)[1])
         del Processes[p] # Remove from list - this is why we needed reverse order

   while (len(Processes) < maxNumProcesses) and (nextVariant < nVariants): # More to do and some spare slots
      start_new(variantList, n_mpi, solver_args)

def report_failed_variants(return_codes, study_directories, solver):
    if all(rc == 0 for rc in return_codes):
        print("\nAll variants finished properly.\n")
        return 0

    failed_variants = [v for v,rc in enumerate(return_codes) if rc != 0]

    recreation_string = "-v"
    print("\n\nPrinting logs from failed variants:")
    print("__________________________________________________________________________\n")
    for fv in failed_variants:
        subprocess.run("cat " + study_directories[fv] + "/log." + solver, shell=True)
        print("__________________________________________________________________________\n\n")
        recreation_string = recreation_string + str(fv) + ","
    recreation_string = recreation_string.rstrip(",")

    print("Number of failed variants:", len(failed_variants))
    print("Failed variants:", failed_variants)
    print("Run 'argo-create-parameter-study.py' with option '", recreation_string,
            "' to re-create failed variants.")

    return 1


def main():
    
    global solver
    global maxNumProcesses
    global nVariants
    global return_codes

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description=app_description)
    parser.add_argument("solver",
                        help="Name of the solver to run the study with.")
    parser.add_argument("--solver-args",
                        help="Additional arguments to pass to the solver. This option is ignored if SLURM is used.\n Default: None",
                        type=str,
                        default="",
                        dest="solver_args")
    parser.add_argument("-d","--dir-pattern",
                        help="Pattern of the directories belonging to the study.",
                        type=str,
                        required=True,
                        dest="dir_pattern")
    parser.add_argument("-nc","--num-cases-parallel",
                        help="The number of cases run in parallel.\n Default: 1",
                        type=int,
                        default=1,
                        dest="numcases")
    parser.add_argument("-j", "--job-mode",
                        help="Submit solver execution as SLURM jobs. Expects to find a SLURM script named solver_name.sbatch.",
                        action="store_true",
                        default=False,
                        dest="job_mode")
    parser.add_argument("-np","--num-mpi-processes",
                        help="Number of MPI processes per case for parallel case execution.\n Has no effect if --job-mode is used.",
                        type=int,
                        default=0,
                        dest="num_mpi_procs")

    args = parser.parse_args()

    # Get list of variation folders
    studyDirectories = trc.pattern_cases(args.dir_pattern + "*")

    if not studyDirectories:
        print("Error: no directories found for study",args.dir_pattern,".")
        print("Exiting.")
        sys.exit()

    # Set the global parameters for parallel run
    solver = args.solver
    nVariants = len(studyDirectories)
    return_codes = [None]*nVariants
    maxNumProcesses = args.numcases

    #--------------------------------------------------------------------------
    #   SLURM batch mode
    #--------------------------------------------------------------------------
    if args.job_mode:
        job_script = solver + ".sbatch";
        if not os.path.isfile(job_script):
            sys.exit("Error: no SLURM script named '" + job_script + "' found.\nExiting.")

        print("Submitting jobs to SLURM using job script", job_script)

        for index in range(len(studyDirectories)):
           pwd = os.getcwd()
           os.chdir(studyDirectories[index])

           print("\nSubmitting job", index+1, "/", len(studyDirectories))
           subprocess.run(["sbatch", "../" + job_script])
           time.sleep(1)

           os.chdir(pwd)

    else:
    #--------------------------------------------------------------------------
    #   Direct execution
    #--------------------------------------------------------------------------
        # User Info
        print("------------ Run Info ----------------")
        print("Running",args.dir_pattern,"using the solver",solver,"with",maxNumProcesses,
                "processes.")
        if args.num_mpi_procs > 1:
            print("Running with MPI using", args.num_mpi_procs, "mpi processes each.")
        print("Start_time: ",datetime.datetime.now().time())

        # Spawn new processes as long as there are variants left and
        # wait until all processes have finished
        while (nextVariant < nVariants) or len(Processes) > 0:
            time.sleep(1) # You may wish to change the time for this
            check_running(studyDirectories, args.num_mpi_procs, args.solver_args)

        print("End_time: ",datetime.datetime.now().time())

        grc = report_failed_variants(return_codes, studyDirectories, solver)
        if grc != 0:
            sys.exit("\n\nError: one or more variants failed to finish properly.\n")

if __name__ == "__main__":
    main()
