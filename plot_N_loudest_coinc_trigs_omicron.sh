#!/bin/bash

N=$1

for i in $(seq 1 $N); do
    echo "Generating plot $i"
    python plot_Nth_loudest_cbc_trig_omicron.py --coinc-trig-file ./1117400416-928800/H1L1-STATMAP_FULL_DATA_FULL_CUMULATIVE_CAT_12H-1117400416-928800.hdf --L1-trigs ./1117400416-928800/L1-HDF_TRIGGER_MERGE_FULL_DATA-1117400416-928800.hdf --H1-trigs ./1117400416-928800/H1-HDF_TRIGGER_MERGE_FULL_DATA-1117400416-928800.hdf --template-file ./1117400416-928800/H1L1-BANK2HDF-1117400416-928800.hdf --output-dir /home/tjmassin/public_html/cbc/DQ/ER7/coinc_trig_omicron_map/1117400416-928800/ --N $i --L1-omicron-dir /home/tjmassin/CBC_DQ/Omicron_trigs/ER7/L1/ --H1-omicron-dir /home/detchar/triggers/ER7/H1/ --omicron-snr-thresh 5
done
