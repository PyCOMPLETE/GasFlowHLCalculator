import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm


# plug-in replacement of old heat load procedure, the fill dict
def get_fill_dict(filln, h5_storage=None, use_dP=True):

    obhl_list = []

    obhl = h5_storage.load_qbs(filln, use_dP=use_dP)
    obhl_list.append(obhl)
    try:
        obhl_instrum = h5_storage.load_special_qbs(filln)
        obhl_list.append(obhl_instrum)
    except FileNotFoundError as err:
        print('Instrumented cell heat loads not found: ', err)

    ms0 = 0.*obhl.timestamps

    dict_out = {}
    for ob in obhl_list:
        for ii, vv in enumerate(ob.variables):
            tv = tm.timber_variable_list()
            tv.t_stamps = ob.timestamps
            tv.ms = ms0
            tv.values = ob.data[:, ii]

            dict_out[vv] = tv

    return dict_out
