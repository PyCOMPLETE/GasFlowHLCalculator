
instrumented_cells_config = {

    'QRLAB_31L2_QBS943.POSST': {

        'circuit_A_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT824.POSST',
                          'QBBI_B31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST',
                          'QBQI_32L2_TT825.POSST'],
        'circuit_B_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT826.POSST',
                          'QBBI_B31L2_TT824.POSST', 'QQBI_31L2_TT826.POSST',
                          'QBQI_32L2_TT825.POSST'],

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
        'magnet_lengths': [15.7, 15.7, 15.7, 5.9],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_13L5_QBS943.POSST': {

        'circuit_A_sensors': ['QBQI_13L5_TT825.POSST', 'QBBI_A13L5_TT826.POSST',
                        'QBBI_B13L5_TT824.POSST', 'QQBI_13L5_TT826.POSST',
                        'QBQI_14L5_TT825.POSST'],
        'circuit_B_sensors': ['QBQI_13L5_TT825.POSST', 'QBBI_A13L5_TT824.POSST',
                        'QBBI_B13L5_TT826.POSST', 'QQBI_13L5_TT824.POSST',
                        'QBQI_14L5_TT825.POSST'],

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
        'magnet_lengths': [15.7, 15.7, 15.7, 5.9],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_13R4_QBS947.POSST': {

        'circuit_A_sensors': ['QBQI_12R4_TT825.POSST', 'QQBI_12R4_TT824.POSST',
                    'QBBI_A13R4_TT826.POSST', 'QBBI_B13R4_TT824.POSST',
                    'QBQI_13R4_TT825.POSST'],
        'circuit_B_sensors': ['QBQI_12R4_TT825.POSST', 'QQBI_12R4_TT826.POSST',
                    'QBBI_A13R4_TT824.POSST', 'QBBI_B13R4_TT826.POSST',
                    'QBQI_13R4_TT825.POSST'],

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['Q1', 'D2', 'D3', 'D4'],
        'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_33L5_QBS947.POSST': {

        'circuit_A_sensors': ['QBQI_34R4_TT825.POSST', 'QQBI_34L5_TT824.POSST',
                'QBBI_B34L5_TT826.POSST', 'QBBI_A34L5_TT824.POSST',
                'QBQI_34L5_TT825.POSST'],
        'circuit_B_sensors': ['QBQI_34R4_TT825.POSST', 'QQBI_34L5_TT826.POSST',
                'QBBI_B34L5_TT824.POSST', 'QBBI_A34L5_TT826.POSST',
                'QBQI_34L5_TT825.POSST'],

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['Q1', 'D2', 'D3', 'D4'],
        'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },
}
