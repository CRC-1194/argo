#!/usr/bin/env bash

# Header of interface height function object is broken. Add it manually before
# agglomeration.
for DIR in benchmark-inter*/;
do
    cd $DIR/postProcessing/interfaceHeight1/0
    # Below, commented out, is the original header as generated by the function object.
    # However, for clearity in the test case context, it has been renamed.
    #echo "#time    height_above_boundary       height_above_location" > height_with_header.dat
    echo "#time    height_above_boundary       major_semi_axis_length" > height_with_header.dat
    cat height.dat >> height_with_header.dat
    cd ../../../../
done

# Need to collect data for interFoam and interIsoFoam separately due to different case base names.
argo-agglomerate-study-data.py benchmark-interFoam-benchmarkpaper_00000_templateCase/postProcessing/interfaceHeight1/0/height_with_header.dat -p benchmarkpaper.parameter -f oscillating_droplet_2D_interFoam
argo-agglomerate-study-data.py benchmark-interIsoFoam-benchmarkpaper_00009_templateCase/postProcessing/interfaceHeight1/0/height_with_header.dat -p benchmarkpaper.parameter -f oscillating_droplet_2D_interIsoFoam
