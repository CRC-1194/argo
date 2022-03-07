#!/usr/bin/env bash

# Remove everything except the header from the data files.
# Remove the first whitespace so the header is read correctly by pandas
for DIR in benchmark-inter*/;
do
    cd $DIR/postProcessing/minMaxU/0
    tail -n +2 fieldMinMax.dat > fieldMinMax_fixed_header.dat
    sed -i 's/ //' fieldMinMax_fixed_header.dat
    cd ../../../../
done

# Need to collect data for interFoam and interIsoFoam separately due to different case base names.
argo-agglomerate-study-data.py benchmark-interFoam-benchmarkpaper_00000_templateCase/postProcessing/minMaxU/0/fieldMinMax_fixed_header.dat -p benchmarkpaper.parameter -f translating_droplet_3D_interFoam
argo-agglomerate-study-data.py benchmark-interIsoFoam-benchmarkpaper_00009_templateCase/postProcessing/minMaxU/0/fieldMinMax_fixed_header.dat -p benchmarkpaper.parameter -f translating_droplet_3D_interIsoFoam
