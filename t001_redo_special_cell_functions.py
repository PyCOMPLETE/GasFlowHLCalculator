import numpy as np

from .h5_storage import H5_storage
from . import heatload_recalc as hlr

cell_description = {
    'QBS': 'QRLAB_31L2_QBS943.POSST',

    'P1': 'QRLAA_29L2_PT961.POSST',
    'P4': 'QRLAA_29L2_PT991.POSST',
    'T1': 'QRLAA_29L2_TT961.POSST',
    'T3': 'QRLAB_31L2_TT943.POSST',

    'T2': 'LBARB_31L2_TT843.POSST',

    'CV': 'QRLAB_31L2_CV943.POSST',
    'EH': 'LBARB_31L2_EH843.POSST',

    'channel_radius': 3.7e-3/2.,
    'roughness': 1e-5,

    'n_channels_tot': 4,

    'length': 53.,

    'circuit_A_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT824.POSST',
                          'QBBI_B31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST',
                          'QBQI_32L2_TT825.POSST'],
    'circuit_B_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT826.POSST',
                          'QBBI_B31L2_TT824.POSST', 'QQBI_31L2_TT826.POSST',
                          'QBQI_32L2_TT825.POSST'],

    'n_channels_circuit_A': 2,
    'n_channels_circuit_B': 2,


    'magnet_names': ['D4', 'D3', 'D2', 'Q1'],
    'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

    'circuit_A_beam': [1,2,1,2],
    'circuit_B_beam': [2,1,2,1],
}

cell_calibration = {
 'R': 54.,
 'Qs': 5.,
 'Kv': 0.39
 }

h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')

filln = 6737
#filln = 6966
obraw = h5_storage.load_special_data_file(filln=filln)

T1 = obraw.dictionary[cell_description['T1']]
T3 = obraw.dictionary[cell_description['T3']]
P1 = obraw.dictionary[cell_description['P1']]
P4 = obraw.dictionary[cell_description['P4']]
CV= obraw.dictionary[cell_description['CV']]
EH = obraw.dictionary[cell_description['EH']]

# T2 = obraw.dictionary[cell_description['T2']]

Q_bs, other = hlr.compute_heat_load(P1, T1, T3, P4, CV, EH,
        Qs_calib=cell_calibration['Qs'],
        Kv_calib=cell_calibration['Kv'],
        R_calib=cell_calibration['R'],
        cell_length=cell_description['length'],
        n_channels=cell_description['n_channels_tot'],
        channel_radius=cell_description['channel_radius'],
        channel_roughness=cell_description['roughness'],
        with_P_drop=True, N_iter_max=100, scale_correction=0.3,
        iter_toll=1e-3)




#####################################################################
# Compute share between the two beam screens
N_iter_max = 100
m_L = other['mass_flow']

n_channels_circuits = [cell_description['n_channels_circuit_%s'%cc]
                            for cc in ['A', 'B']]
magnet_lengths_circuits = [cell_description['magnet_lengths']
                                for _ in ['A', 'B']]

out_sensor_names_circuits = [cell_description['circuit_%s_sensors'%cc][1:]
                                for cc in ['A', 'B']]

in_sensor_names_circuits = [cell_description['circuit_%s_sensors'%cc][:-1]
                                for cc in ['A', 'B']]

magnet_beam_circuits = [cell_description['circuit_%s_beam'%cc]
                                for cc in ['A', 'B']]

T_out_magnets_circuits = [[obraw.dictionary[vv] for vv in
        out_sensor_names_circuits[ii]] for ii in [0, 1]]

T_in_magnets_circuits = [[obraw.dictionary[vv] for vv in
        in_sensor_names_circuits[ii]] for ii in [0, 1]]


Qbs_magnets_circuits, other = hlr.compute_heat_loads_instrumented_cell(
        mass_flow = m_L, P1=P1,
        T_in_magnets_circuits=T_in_magnets_circuits,
        T_out_magnets_circuits=T_out_magnets_circuits,
        magnet_lengths_circuits=magnet_lengths_circuits,
        n_channels_circuits=n_channels_circuits,
        channel_radius=cell_description['channel_radius'],
        channel_roughness=cell_description['roughness'],
        dp_toll = 0.001, N_iter_max=200)

QBS_name = cell_description['QBS']

dict_output = {
    QBS_name: Q_bs}
magnet_names = cell_description['magnet_names']
for i_circ in [0, 1]:
    magnets_beam_c = magnet_beam_circuits[i_circ]

    for i_mag, name_mag in enumerate(magnet_names):
        dict_output[QBS_name.split('.POSST')[0]
               +'_%sB%s.POSST'%(name_mag, magnets_beam_c[i_mag])]=\
        Qbs_magnets_circuits[i_circ][i_mag]

for name_mag in magnet_names:
    dict_output[QBS_name.split('.POSST')[0]
        +'_%s.POSST'%(name_mag)] = \
         dict_output[QBS_name.split('.POSST')[0]
               +'_%sB%s.POSST'%(name_mag, 1)] +\
         dict_output[QBS_name.split('.POSST')[0]
               +'_%sB%s.POSST'%(name_mag, 2)]

# Hide last magnet
name_last_magnet = magnet_names[-1]

for kk in list(dict_output.keys()):
    for bb in [1,2]:
        if '_%sB%d'%(name_last_magnet, bb) in kk:
            dict_output[kk] *= 0.

# Some plots

import matplotlib.pyplot as plt
plt.close('all')

for i_mag, name_mag in enumerate(magnet_names):
    fig = plt.figure(i_mag+1)
    ax = fig.add_subplot(111)

    nn = QBS_name.replace('.POSST', '_%s.POSST'%name_mag)
    nnb1 = QBS_name.replace('.POSST', '_%sB1.POSST'%name_mag)
    nnb2 = QBS_name.replace('.POSST', '_%sB2.POSST'%name_mag)

    ax.plot(dict_output[nn], color='k')
    ax.plot(dict_output[nnb1], color='b')
    ax.plot(dict_output[nnb2], color='r')

    fig.suptitle(name_mag)
plt.show()
