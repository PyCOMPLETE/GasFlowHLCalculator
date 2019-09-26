from h5_storage import H5_storage

cell_description = {
    'P1': 'QRLAA_29L2_PT961.POSST',
    'P4': 'QRLAA_29L2_PT991.POSST',
    'T1': 'QRLAA_29L2_TT961.POSST',
    'T3': 'QRLAB_31L2_TT943.POSST',

    'circuit_1_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT824.POSST',
        'QBBI_B31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST', 'QBQI_32L2_TT825.POSST'],
    'circuit_2_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT826.POSST',
        'QBBI_B31L2_TT824.POSST', 'QQBI_31L2_TT826.POSST', 'QBQI_32L2_TT825.POSST'],

    'magnet_sequence': ['Q1', 'D2', 'D3', 'D4'],

    'b1_circuit': [1,2,1,2],
    'b2_circuit': [2,1,2,1]

}

cell_calibration = {
 'R': 38.,
 'Qs': 8.,
 'Kv': 0.39
 }

h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')

obraw = h5_storage.load_special_data_file(filln=6737)
