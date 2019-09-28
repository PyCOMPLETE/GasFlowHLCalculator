from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager

cal_manager = CalibrationManager(calibration_config=calibration_config)

calib = Calibration(calibration_csv_file=calibration_config[1]['file'])
