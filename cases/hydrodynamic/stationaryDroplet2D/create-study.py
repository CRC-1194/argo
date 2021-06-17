#!/usr/bin/env python3

from argparse import ArgumentParser
from subprocess import call

usage = """A wrapper for pyFoamRunParameterVariation.py that generates the
directory structure for a parameter study.

Meshes are not generated and preprocessing is not done.
Used to prepare for execution on a cluster.

Usage: ./create-study.py -c templateCase -p paramFile -s studyName"""

parser = ArgumentParser(usage=usage)
parser.add_argument("-c", "--case", dest="casedir",
                    help="Template case directory.",
                    metavar="CASEDIR")

parser.add_argument("-p", "--parameter-file", dest="paramfile",
                    help="PyFoam parameter file used by pyFoamRunParameterVariation.py.",
                    metavar="PARAMFILE")

parser.add_argument("-s", "--study-name", dest="studyname",
                    help="Name of the parameter study.",
                    metavar="STUDYNAME")
args = vars(parser.parse_args())

call(["pyFoamRunParameterVariation.py", "--no-execute-solver", "--no-server-process",
      "--no-mesh-create", "--no-case-setup", "--cloned-case-prefix=%s" % args['studyname'],
      "--every-variant-one-case-execution",
      "--create-database", args['casedir'], args['paramfile']])
