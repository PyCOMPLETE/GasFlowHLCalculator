import numpy as np

from .calibration_config import calibration_config
from .calibration import Calibration, CalibrationManager
from .h5_storage import H5_storage
from . import recalc_multiple_circuits as rmc

filln = 6737
with_P_drop = True

circuit_selection = 'full_lhc'
circuit_selection = 'all_instrumented'

cal_manager = CalibrationManager(calibration_config=calibration_config)
h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')

if circuit_selection == 'full_lhc':
    obraw = h5_storage.load_data_file(filln=filln)
elif circuit_selection == 'all_instrumented':
    obraw = h5_storage.load_special_data_file(filln=filln)
else:
    raise ValueError

calibration = cal_manager.get_calibration(obraw.timestamps[0])

obhl_store, other = rmc.recalc_multiple_circuits(obraw,
        calibration, circuit_selection=circuit_selection,
        with_P_drop=with_P_drop)

h5_storage.store_qbs(filln=filln, qbs_ob=obhl_store, use_dP=with_P_drop)

if circuit_selection == 'full_lhc':
    h5_storage.store_qbs(filln=filln, qbs_ob=obhl_store, use_dP=with_P_drop)
elif circuit_selection == 'all_instrumented':
    h5_storage.store_special_qbs(filln=filln, qbs_ob=obhl_store)
else:
    raise ValueError
