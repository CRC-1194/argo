#!/usr/bin/env python3 

from optparse import OptionParser
import os
import sys
import time
from subprocess import call

usage = """A wrapper for pyFoamRunParameterVariation.py that generates the
directory structure for a parameter study. 

Used to prepare for execution on a cluster. Performs meshing and field
initialization.

Usage: ./create-study.py -c templateCase -p paramFile -s studyName"""

parser = OptionParser(usage=usage)

parser.add_option("-c", "--case", dest="casedir",
                  help="Template case directory.", 
                  metavar="CASEDIR")

parser.add_option("-p", "--parameter-file", dest="paramfile", 
                  help="PyFoam parameter file used by pyFoamRunParameterVariation.py.", 
                  metavar="PARAMFILE")

parser.add_option("-s", "--study-name", dest="studyname", 
                  help="Name of the parameter study.", 
                  metavar="STUDYNAME")

parser.add_option("-j", "--job-mode", action="store_true", dest="use_slurm",
                  help="Submit meshing and initialization to SLURM.",
                  metavar="SLURM")

(options, args) = parser.parse_args()

if ((options.casedir == None) or  
    (options.paramfile == None) or 
    (options.studyname == None)): 
    print ("Error: case or parameter option not used. Use --help option for more information.") 
    sys.exit(1)

(options, args) = parser.parse_args()

call(["pyFoamRunParameterVariation.py", "--no-execute-solver", "--no-server-process", 
      "--no-mesh-create", "--no-case-setup", "--cloned-case-prefix=%s" % options.studyname, 
      "--every-variant-one-case-execution",
      "--create-database", options.casedir, options.paramfile])

# Mesh each case and set initial fields
parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                  if os.path.isdir(parameter_dir) \
                  and options.studyname in parameter_dir]
parameter_dirs.sort()

for parameter_dir in parameter_dirs: 
    pwd = os.getcwd()
    os.chdir(parameter_dir)

    print(parameter_dir)

    if (options.use_slurm):
        call(["sbatch", "../caseSetup.sbatch"])
        time.sleep(2)
    else: 
        call("blockMesh")
        call("setAlphaField")
        call("decomposePar", "-force")

    os.chdir(pwd)
