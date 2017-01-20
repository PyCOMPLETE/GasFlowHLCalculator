# GasFlowHLCalculator
Reconstruct the LHC heat loads from valves, gas flows and temperatures

Usage of the most important scripts and modules:

qbs_fill.py:

- This module serves as the frontend for other scripts.
- It only needs the fill number as input.
- If the corrseponding data is already stored as h5 it returns it. Otherwise the computation is performed.
- Functions:
  - compute_qbs_fill: Obtain all cell heat loads for a given fill as aligned data.
  - compute_qbs_arc_avg: Takes aligned data and returns the arc averages as aligned data.
  - arc_histograms: Takes recomputed cell heat loads and groups them to the arcs in form of a dictionary. This can then be used for histograms.
  - get_fill_dict: This can be used to replace (update) the old fill dicts that are used in the 016\_ and other scripts. Currently only works for the arcs.

create_csv.m:

- This GNU Octave script reads the data_QBS_LHC.m file and creates a versioned csv file.
- Currently there is no such script for the Helium_proporties.m as I do not expect it to change.

data_qbs.py:

- This module parses the csv file. 
- Some constants such as the cooling pipe diameter are hardcoded.

h5_storage.py:

- The most recent version is defined here.
- This module controls the read and write access to the directory where the offline recomputed heat loads are stored.

compute_QBS_LHC.py, valve_LT.py, Pressure_drop.py

- The main calculations are performed here. 
- These modules are translations of Benjmain Bradu's matlab scripts.

001_store_recalculated.py

- This script computes the qbs for all fills in 2015 and 2016 and stores the results as h5.

004_compare_qbs_versions.py

- This not polished script can compare the results from 2 versions of csv files.

data_qbs_lhc\_*.csv

- Tab separated values
- Timber variables
- Valve characteristics
- Cell names 
- and more

variable_list.txt

- Used by the 001f script (not part of this module) that downloads the necessary timber data.

Issues:

- The naming convention of the special instrumented cells can be confusing.
  - Notation in the original MATLAB scripts:
    - 12R4
    - 32R4 (broken sensor)
    - 13L5 (reversed gas flow)
  - Notation of the special cell timber variables:
    - 13L5
    - 33L5 (broken sensor)
    - 13R4 (reversed gas flow)
- There are duplicate cell names in the data_qbs_lhc\_\*.csv
  - These only affect the LSS, not the arcs
  -  ['05L4_947', '05R4_947', '05L6_947', '05R6_947']
- The Quadrupoles are not yet identified, the logged data of those can not yet be replaced.
