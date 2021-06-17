#!/usr/bin/env python3

from argparse import ArgumentParser
import datetime 
import subprocess
import sys
import time

import testReportCore as trc

# The solution to control the number parallel running processes is taken / adapted from
# https://stackoverflow.com/questions/18123325/always-run-a-constant-number-of-subprocesses-in-parallel

nextVariant = 0
maxNumProcesses = 1
nVariants = 1
solver = ""
Processes = []

def start_new(variantList):
    """ Start a new subprocess if there is work to do """
    global nextVariant
    global Processes

    if nextVariant < nVariants:
        # Make the script wait for the last spawned subprocess.
        # This makes it more likely that the end time info
        # will be correct
        if nextVariant == nVariants-1:
            proc = subprocess.Popen([solver, "-case", variantList[nextVariant]]).wait()
        else:
            proc = subprocess.Popen([solver, "-case", variantList[nextVariant]])
        print ("Started to process variant", variantList[nextVariant].rsplit('/',1)[1])
        nextVariant += 1
        Processes.append(proc)

def check_running(variantList):
   """ Check any running processes and start new ones if there are spare slots."""
   global Processes
   global nextVariant

   for p in range(len(Processes)-1,0,-1): # Check the processes in reverse order
      if Processes[p].poll() is not None: # If the process hasn't finished will return None
         del Processes[p] # Remove from list - this is why we needed reverse order

   while (len(Processes) < maxNumProcesses) and (nextVariant < nVariants): # More to do and some spare slots
      start_new(variantList)


def main():
    
    global solver
    global maxNumProcesses
    global nVariants

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Execute a parameter study in parallel using N processes.")
    parser.add_argument("solver",
                        help="Name of the solver to run the study with.")
    parser.add_argument("-d","--dir-pattern",
                        help="Pattern of the directories belonging to the study.",
                        type=str,
                        required=True,
                        dest="dir_pattern")
    parser.add_argument("-n","--num-processes",
                        help="The number of processes running in parallel",
                        type=int,
                        required=True,
                        dest="numprocesses")

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
    maxNumProcesses = args.numprocesses

    # User Info
    print("------------ Run Info ----------------")
    print("Running",args.dir_pattern,"using the solver",solver,"with",maxNumProcesses,
            "processes.")
    print("Start_time: ",datetime.datetime.now().time())

    # Spawn processes
    check_running(studyDirectories) # This will start the max processes running
    while (nextVariant < nVariants): # Some thing still going on.
        time.sleep(0.1) # You may wish to change the time for this
        check_running(studyDirectories)

    print("End_time: ",datetime.datetime.now().time())
    print("Done!")

if __name__ == "__main__":
    main()