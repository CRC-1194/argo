#!/usr/bin/env bash

# Comparison data 

##  Jones, B. W. S., Malan, A. G., & Ilangakoon, N. A. (2019). The initialisation of volume fractions for unstructured grids using implicit surface definitions. Computers and Fluids, 179, 194–205. https://doi.org/10.1016/j.compfluid.2018.10.021¶

## Table 3

# Level Set Sphere 

## Create the parameter study folder structure
argo-create-parameter-study.py j.compfluid.2018.10.021-table3-levelSet.parameter

## Mesh the solution domain with gmsh using the custom init script 
argo-initilize-parameter-study.py -m none -f j.compfluid.2018.10.021-table3-init-levelSet.sh j.compfluid.2018.10.021-table3-levelSet 

## Run the study  
argo-run-study.py -d  j.compfluid.2018.10.021-table3-levelSet surfaceInitVolumeFraction

## Agglomerate data
argo-agglomerate-study-data.py -p j.compfluid.2018.10.021-table3-levelSet.parameter j.compfluid.2018.10.021-table3-levelSet_00000_templateCase/vof-init-results-SMCA.csv -f j.compfluid.2018.10.021-table3-levelSet
