# GasFlowHLCalculator
Reconstruct the LHC heat loads from valves, gas flows and temperatures

Usage of the most important scripts and modules:

create_csv.m:

- This GNU Octave script reads the data_QBS_LHC.m file and creates a versioned csv file.
- Currently there is no such script for the Helium_proporties.m as I do not expect it to change.

data_qbs.py:

- This module parses the csv file. 
- Some constants such as the cooling pipe diameter are hardcoded.

h5_storage.py:

- The most recent version is defined here.
- This module controls the read and write access to the directory where the offline recomputed heat loads are stored.

qbs_fill.py:

- This module serves as the frontend for other scripts.
- It only needs the fill number as input.
- If the corrseponding data is already stored as h5 it returns it. Otherwise the computation is performed.

compute_QBS_LHC.py, valve_LT.py, Pressure_drop.py

- The main calculations are performed here. 
- These modules are translations of Benjmain Bradu's matlab scripts.

001_store_recalculated.py

- This script computes the qbs for all fills in 2015 and 2016 and stores the results as h5.

004_compare_qbs_versions.py

- This not polished script can compare the results from 2 versions of csv files.

data_qbs_lhc_*.csv

- Tab separated values
- Timber variables
- Valve characteristics
- Cell names 
- and more

variable_list.txt

- Used by the 001f script (not part of this module) that downloads the necessary timber data.
