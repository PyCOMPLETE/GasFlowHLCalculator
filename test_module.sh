#!/bin/bash
set -e

cd $(dirname "$0")

python 005_test_qbs_lhc.py --noshow
python valve_LT.py
python Pressure_drop.py

cd '/afs/cern.ch/work/l/lhcscrub/LHC_2016_25ns_beforeTS1'
python ./027_special_instrumented_cells.py 5219
python ./026_recalculate_heat_loads.py 5219
python ./016_heat_load_arcs_vs_SAMs.py 5219
