import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm


# plug-in replacement of old heat load procedure, the fill dict
def get_fill_dict(filln, h5_storage=None, use_dP=True):

    fname =  h5_storage.get_qbs_file(filln, use_dP=use_dP)
    obhl = ob = mfm.h5_to_obj(fname)

    ms0 = 0.*obhl.timestamps

    dict_out = {}
    for ii, vv in enumerate(obhl.variables):
        tv = tm.timber_variable_list()
        tv.t_stamps = obhl.timestamps
        tv.ms = ms0
        tv.values = obhl.data[:, ii]

        dict_out[vv] = tv

    return dict_out
