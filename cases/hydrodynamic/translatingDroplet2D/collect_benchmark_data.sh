#!/usr/bin/env bash

# Fuse data from different function objects into a single csv file for
# further agglomeration
for DIR in benchmark-inter*/;
do
    cd $DIR/postProcessing
    python -c "import pandas as pd; \
        dfmax = pd.read_csv('minMaxUError/0/fieldMinMax.dat', delim_whitespace=True, \
                            comment='#', header=None, names=['time', 'min_error_velocity', 'max_error_velocity']); \
        dfl1 = pd.read_csv('l1normUError/0/volFieldValue.dat', delim_whitespace=True, comment='#', header=None); \
        dfl2 = pd.read_csv('l2normUError/0/volFieldValue.dat', delim_whitespace=True, comment='#', header=None); \
        dfmax['mean_absolute_error_velocity'] = dfl1[1]; \
        dfmax['root_mean_square_deviation_velocity'] = dfl2[1]; \
        dfmax.to_csv('velocity_data.csv', index=False)"
    cd ../../
done

# Need to collect data for interFoam and interIsoFoam separately due to different case base names.
argo-agglomerate-study-data.py benchmark-interFoam-benchmarkpaper_00000_templateCase/postProcessing/velocity_data.csv -p benchmarkpaper.parameter -f translating_droplet_2D_interFoam
argo-agglomerate-study-data.py benchmark-interIsoFoam-benchmarkpaper_00009_templateCase/postProcessing/velocity_data.csv -p benchmarkpaper.parameter -f translating_droplet_2D_interIsoFoam