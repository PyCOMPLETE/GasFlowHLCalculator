
instrumented_cells_config = {

    'QRLAB_31L2_QBS943.POSST': {

        'circuit_A_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT824.POSST',
                          'QBBI_B31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST',
                          'QBQI_32L2_TT824.POSST'],
        'circuit_B_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT826.POSST',
                          'QBBI_B31L2_TT824.POSST', 'QQBI_31L2_TT826.POSST',
                          'QBQI_32L2_TT826.POSST'],

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'C1-FE': 'QBQI_31L2_FE826.POSST',
        'C2-FE': 'QBQI_31L2_FE824.POSST',

        'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
        'magnet_lengths': [15.7, 15.7, 15.7, 5.9],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_13L5_QBS943.POSST': {

        'circuit_A_sensors': ['QBQI_13L5_TT825.POSST', 'QBBI_A13L5_TT826.POSST',
                        'QBBI_B13L5_TT824.POSST', 'QQBI_13L5_TT826.POSST',
                        'QBQI_14L5_TT826.POSST'],
        'circuit_B_sensors': ['QBQI_13L5_TT825.POSST', 'QBBI_A13L5_TT824.POSST',
                        'QBBI_B13L5_TT826.POSST', 'QQBI_13L5_TT824.POSST',
                        'QBQI_14L5_TT824.POSST'],

        'C1-FE': 'QBQI_13L5_FE824.POSST',
        'C2-FE': 'QBQI_13L5_FE826.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
        'magnet_lengths': [15.7, 15.7, 15.7, 5.9],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_13R4_QBS947.POSST': {

        ## in Run 3, TT5 are inverted between beam 1 and beam 2
        'circuit_A_sensors': ['QBQI_12R4_TT825.POSST', 'QQBI_12R4_TT824.POSST',
                    'QBBI_A13R4_TT826.POSST', 'QBBI_B13R4_TT824.POSST',
                    'QBQI_13R4_TT824.POSST'],
        'circuit_B_sensors': ['QBQI_12R4_TT825.POSST', 'QQBI_12R4_TT826.POSST',
                    'QBBI_A13R4_TT824.POSST', 'QBBI_B13R4_TT826.POSST',
                    'QBQI_13R4_TT826.POSST'],

        ## Run2? Theoretical??
        ## 'circuit_A_sensors': ['QBQI_12R4_TT825.POSST', 'QQBI_12R4_TT824.POSST',
        ##             'QBBI_A13R4_TT826.POSST', 'QBBI_B13R4_TT824.POSST',
        ##             'QBQI_13R4_TT826.POSST'],
        ## 'circuit_B_sensors': ['QBQI_12R4_TT825.POSST', 'QQBI_12R4_TT826.POSST',
        ##             'QBBI_A13R4_TT824.POSST', 'QBBI_B13R4_TT826.POSST',
        ##             'QBQI_13R4_TT824.POSST'],

        'C1-FE': 'QBQI_12R4_FE824.POSST',
        'C2-FE': 'QBQI_12R4_FE826.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['Q1', 'D2', 'D3', 'D4'],
        'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_33L5_QBS947.POSST': {

        ## in Run 3, TT5 are inverted between beam 1 and beam 2
        'circuit_A_sensors': ['QBQI_34R4_TT825.POSST', 'QQBI_34L5_TT824.POSST',
                'QBBI_B34L5_TT826.POSST', 'QBBI_A34L5_TT824.POSST',
                'QBQI_34L5_TT824.POSST'],
        'circuit_B_sensors': ['QBQI_34R4_TT825.POSST', 'QQBI_34L5_TT826.POSST',
                'QBBI_B34L5_TT824.POSST', 'QBBI_A34L5_TT826.POSST',
                'QBQI_34L5_TT826.POSST'],

        ## Run2? Theoretical??
        ## 'circuit_A_sensors': ['QBQI_34R4_TT825.POSST', 'QQBI_34L5_TT824.POSST',
        ##         'QBBI_B34L5_TT826.POSST', 'QBBI_A34L5_TT824.POSST',
        ##         'QBQI_34L5_TT826.POSST'],
        ## 'circuit_B_sensors': ['QBQI_34R4_TT825.POSST', 'QQBI_34L5_TT826.POSST',
        ##         'QBBI_B34L5_TT824.POSST', 'QBBI_A34L5_TT826.POSST',
        ##         'QBQI_34L5_TT824.POSST'],

        'C1-FE': 'QBQI_34R4_FE824.POSST',
        'C2-FE': 'QBQI_34R4_FE826.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['Q1', 'D2', 'D3', 'D4'],
        'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAB_15R2_QBS943.POSST': {

        'circuit_A_sensors': ['QBQI_16R2_TT825.POSST', 'QBBI_B16R2_TT826.POSST',
                              'QBBI_A16R2_TT824.POSST', 'QQBI_15R2_TT826.POSST',
                              'QBQI_15R2_TT826.POSST'],
        'circuit_B_sensors': ['QBQI_16R2_TT825.POSST', 'QBBI_B16R2_TT824.POSST',
                              'QBBI_A16R2_TT826.POSST', 'QQBI_15R2_TT824.POSST',
                              'QBQI_15R2_TT824.POSST'],

        'C1-FE': 'QBQI_16R2_FE824.POSST',
        'C2-FE': 'QBQI_16R2_FE826.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
        'magnet_lengths': [15.7, 15.7, 15.7, 5.9],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAA_17L6_QBS943.POSST': {

        'circuit_A_sensors': ['QBQI_17L6_TT825.POSST', 'QBBI_A17L6_TT824.POSST',
                              'QBBI_B17L6_TT826.POSST', 'QQBI_17L6_TT824.POSST',
                              'QBQI_18L6_TT824.POSST'],
        'circuit_B_sensors': ['QBQI_17L6_TT825.POSST', 'QBBI_A17L6_TT826.POSST',
                              'QBBI_B17L6_TT824.POSST', 'QQBI_17L6_TT826.POSST',
                              'QBQI_18L6_TT826.POSST'],

        'C1-FE': 'QBQI_17L6_FE826.POSST',
        'C2-FE': 'QBQI_17L6_FE824.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
        'magnet_lengths': [15.7, 15.7, 15.7, 5.9],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAB_27L8_QBS947.POSST': {

        ## in Run 3, TT5 are inverted between beam 1 and beam 2
        'circuit_A_sensors': ['QBQI_29L8_TT825.POSST', 'QQBI_28L8_TT826.POSST',
                              'QBBI_B28L8_TT824.POSST', 'QBBI_A28L8_TT826.POSST',
                              'QBQI_28L8_TT826.POSST'],
        'circuit_B_sensors': ['QBQI_29L8_TT825.POSST', 'QQBI_28L8_TT824.POSST',
                              'QBBI_B28L8_TT826.POSST', 'QBBI_A28L8_TT824.POSST',
                              'QBQI_28L8_TT824.POSST'],

        ## Didn't exist in Run2, Theoretical??
        ## 'circuit_A_sensors': ['QBQI_29L8_TT825.POSST', 'QQBI_28L8_TT826.POSST',
        ##                       'QBBI_B28L8_TT824.POSST', 'QBBI_A28L8_TT826.POSST',
        ##                       'QBQI_28L8_TT824.POSST'],
        ## 'circuit_B_sensors': ['QBQI_29L8_TT825.POSST', 'QQBI_28L8_TT824.POSST',
        ##                       'QBBI_B28L8_TT826.POSST', 'QBBI_A28L8_TT824.POSST',
        ##                       'QBQI_28L8_TT826.POSST'],

        'C1-FE': 'QBQI_29L8_FE826.POSST',
        'C2-FE': 'QBQI_29L8_FE824.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['Q1', 'D2', 'D3', 'D4'],
        'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },

    'QRLAD_33R2_QBS947.POSST': {

        ## in Run 3, TT5 seem inverted between beam 1 and beam 2
        'circuit_A_sensors': ['QBQI_32R2_TT825.POSST', 'QQBI_32R2_TT824.POSST',
                              'QBBI_A33R2_TT826.POSST', 'QBBI_B33R2_TT824.POSST',
                              'QBQI_33R2_TT824.POSST'],
        'circuit_B_sensors': ['QBQI_32R2_TT825.POSST', 'QQBI_32R2_TT826.POSST',
                              'QBBI_A33R2_TT824.POSST', 'QBBI_B33R2_TT826.POSST',
                              'QBQI_33R2_TT826.POSST'],

        ## Didn't exist in Run2, Theoretical??
        ## 'circuit_A_sensors': ['QBQI_32R2_TT825.POSST', 'QQBI_32R2_TT824.POSST',
        ##                       'QBBI_A33R2_TT826.POSST', 'QBBI_B33R2_TT824.POSST',
        ##                       'QBQI_33R2_TT826.POSST'],
        ## 'circuit_B_sensors': ['QBQI_32R2_TT825.POSST', 'QQBI_32R2_TT826.POSST',
        ##                       'QBBI_A33R2_TT824.POSST', 'QBBI_B33R2_TT826.POSST',
        ##                       'QBQI_33R2_TT824.POSST'],

        'C1-FE': 'QBQI_32R2_FE824.POSST',
        'C2-FE': 'QBQI_32R2_FE826.POSST',

        'n_channels_circuit_A': 2,
        'n_channels_circuit_B': 2,

        'magnet_names': ['Q1', 'D2', 'D3', 'D4'],
        'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

        'circuit_A_beam': [1,2,1,2],
        'circuit_B_beam': [2,1,2,1],

    },
}
