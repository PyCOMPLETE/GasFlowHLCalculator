# GasFlowHLCalculator

This is a tool to reconstruct the heat loads on the LHC beam screens based on cryogenics measurements and settings.
It provides tools to:
 - Download the raw data from the LHC logging database
 - Store the raw data for each fill in hdf5 format
 - Compute the corresponding heat loads 
 - Store the recomputed heat-load data for each fill in hdf5 format
 - Provide access to the recomputed data with an interface similar to the one used with the logging database
 
The code relies on the [LHCMeasurementTools](https://github.com/PyCOMPLETE/LHCMeasurementTools) package

Cryogenics calibration info is stored in a separate repository:

https://github.com/PyCOMPLETE/LHCCryoHeatLoadCalibration

