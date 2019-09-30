import numpy as np

from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager
from h5_storage import H5_storage
import recalc_multiple_circuits as rmc

cal_manager = CalibrationManager(calibration_config=calibration_config)
h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')
with_P_drop = True

filln = 6737

obraw = h5_storage.load_data_file(filln=filln)
calibration = cal_manager.get_calibration(obraw.timestamps[0])

obhl_store, other = rmc.recalc_multiple_circuits(obraw,
        calibration, circuit_selection='full_lhc',
        with_P_drop=with_P_drop)

h5_storage.store_qbs(filln=filln, qbs_ob=obhl_store, use_dP=with_P_drop)
