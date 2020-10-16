#!/sr/bin/env bash

# Create the sphere convergence study.
./create_parameter_study.py --slurm --parameter_file blockMesh.parameter --surface sphere --mesh_generator blockMesh 2>&1 | tee log.create
