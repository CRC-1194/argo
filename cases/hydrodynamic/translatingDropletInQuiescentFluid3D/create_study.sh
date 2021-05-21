#!/usr/bin/env bash
pyFoamRunParameterVariation.py \
    --allow-derived-changes \
    --every-variant-one-case-execution \
    --create-database \
    --no-server-process \
    --no-execute-solver \
    --parameter-file=default.parameter \
    --mesh-create-script=blockMesh.sh \
    --cloned-case-prefix=$3 \
    $1 \
    $2
