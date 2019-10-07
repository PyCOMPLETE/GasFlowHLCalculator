import pandas

import LHCMeasurementTools.TimestampHelpers as th

class Calibration(object):

    def __init__(self, calibration_csv_file):

        self.calibration_csv_file = calibration_csv_file

        self.caldata = pandas.read_csv(calibration_csv_file)
        self.caldata.set_index('Qbs', inplace=True,
                verify_integrity=True)

    @property
    def circuits(self):
        return list(self.caldata.index)


    def get_circuit(self, name):

        dat = self.caldata.loc[name, :]
        ddd = {
            'QBS': name,
            'P1': dat.P1,
            'P4': dat.P4,
            'T1': dat.T1,
            'T3': dat.T3,
            'T2': dat.T2,
            'CV': dat.CV1,
            'EH': dat.QEH,

            'channel_radius': 3.7e-3/2.,
            'roughness': 1e-5,

            'n_channels_tot': dat.nc,
            'length': float(dat.L),
            'R_calib': float(dat.R),
            'Qs_calib': float(dat.QS),
            'Kv_calib': float(dat.Kvmax)
            }
        return ddd


class CalibrationManager(object):

    def __init__(self, calibration_config):

        self.calibration_config = calibration_config
        self.names = [cc['name'] for cc in calibration_config]
        self.files = [cc['file'] for cc in calibration_config]
        self.calibrations = [Calibration(ff) for ff in self.files]

        self.start_timestamps = [th.localtime2unixstamp(cc['start'])
                for cc in calibration_config]

        self.end_timestamps = [th.localtime2unixstamp(cc['end'])
                for cc in calibration_config]


    def get_calibration(self, t_stamp):
        for i_calib, calib in enumerate(self.calibrations):
            if ((t_stamp > self.start_timestamps[i_calib]) and
                    (t_stamp < self.end_timestamps[i_calib])):
                return calib
        raise ValueError('Calibration not found!')


