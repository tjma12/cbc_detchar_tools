#!/bin/bash

N=$1
cache=$2
ifo=$3
outdir=$4

for i in $(seq 1 $N); do
    echo "Generating plot $i"
    python plot_loudest_daily_pycbc_xml_omicron.py --single-ifo-trigs $cache --ifo $ifo --ranking-statistic newsnr --output-file ${outdir}/${ifo}_loudest_daily_pycbc_omicron_${i}.png --plot-window 20 --N $i
done
