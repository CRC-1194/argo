#!/usr/bin/env python3 

from optparse import OptionParser
import os
import sys
from subprocess import call

usage = """A wrapper for pyFoamRunParameterVariation.py that generates the
directory structure for a parameter study. 

Meshes are not generated and preprocessing is not done. 
Used to prepare for execution on a cluster.

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

    if (args.slurm_run):
        call("sbatch", "../caseSetup.sbatch")
    else: 
        call(args.solver)

    os.chdir(pwd)
