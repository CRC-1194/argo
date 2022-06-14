#!/usr/bin/env bash

# Header of interface height function object is broken. Add it manually before
# agglomeration.
for DIR in benchmark-inter*/;
do
    cd $DIR/postProcessing/interfaceHeight1/0
    # FIXME: ensure that the equilibrium interface height, 0.0013 currently, is set correctly.
    python -c "import pandas as pd; \
        df = pd.read_csv('height.dat', delim_whitespace=True, comment='#', header=None, \
                          names=['time', 'height_above_boundary', 'height_above_location']); \
        df['amplitude_at_center'] = df['height_above_location'] - 0.0013; \
        df.to_csv('interface_data.csv', index=False)"
    cd ../../../../
done

# Need to collect data for interFoam and interIsoFoam separately due to different case base names.
argo-agglomerate-study-data.py benchmark-interFoam-benchmarkpaper_00000_templateCase/postProcessing/interfaceHeight1/0/interface_data.csv -p benchmarkpaper.parameter -f oscillating_wave_2D_interFoam
argo-agglomerate-study-data.py benchmark-interIsoFoam-benchmarkpaper_00009_templateCase/postProcessing/interfaceHeight1/0/interface_data.csv -p benchmarkpaper.parameter -f oscillating_wave_2D_interIsoFoam
argo-agglomerate-study-data.py benchmark-interFlow-benchmarkpaper_00009_templateCase/postProcessing/interfaceHeight1/0/interface_data.csv -p benchmarkpaper.parameter -f oscillating_wave_2D_interFlow
