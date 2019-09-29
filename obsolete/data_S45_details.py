# Note: only the cell_timber_vars_dict and the cell_list is supposed to be used from this file.
# The reason is that this data structure is better suited in a dictionary, not lists.


# Data extracted from QBS scripts on 31/08/2016.
# New cell added in June 2017
# The name of this file is outdated, as the new cell is not is S45. But let's not be pedantic.

EH84x_list = ['LQATI_12R4_EH847.POSST', 'LQOAN_34R4_EH847.POSST', 'LBALB_13L5_EH843.POSST', 'LBARB_31L2_EH843.POSST']
TT84x_list = ['LQATI_12R4_TT847.POSST', 'LQOAN_34R4_TT847.POSST', 'LBALB_13L5_TT843.POSST', 'LBARB_31L2_TT843.POSST']
CV94x_list = ['QRLAA_13R4_CV947.POSST', 'QRLAA_33L5_CV947.POSST', 'QRLAA_13L5_CV943.POSST', 'QRLAB_31L2_CV943.POSST']
PT961_list = ['QRLAA_13R4_PT961.POSST', 'QRLAA_33L5_PT961.POSST', 'QRLAA_13L5_PT961.POSST', 'QRLAA_29L2_PT961.POSST']
PT991_list = ['QRLAA_13R4_PT991.POSST', 'QRLAA_33L5_PT991.POSST', 'QRLAA_13L5_PT991.POSST', 'QRLAA_29L2_PT991.POSST']
TT94x_list = ['QRLAA_13R4_TT947.POSST', 'QRLAA_33L5_TT947.POSST', 'QRLAA_13L5_TT943.POSST', 'QRLAB_31L2_TT943.POSST']
TT961_list = ['QRLAA_13R4_TT961.POSST', 'QRLAA_33L5_TT961.POSST', 'QRLAA_13L5_TT961.POSST', 'QRLAA_29L2_TT961.POSST']

cell_list = ['13R4', '33L5', '13L5', '31L2']
cell_list_pre_EYETS16 = ['13R4', '33L5', '13L5']

Qs_list = [7,7,7,8]
R_list = [40,40,40,38]
Kv_list = [0.39,0.39,0.39, 0.39]

dd = {}

for index, cell in enumerate(cell_list):
    dd[cell] = cc = {}
    cc['EH84x'] = EH84x_list[index]
    cc['TT84x'] = TT84x_list[index]
    cc['CV94x'] = CV94x_list[index]
    cc['PT961'] = PT961_list[index]
    cc['PT991'] = PT991_list[index]
    cc['TT94x'] = TT94x_list[index]
    cc['TT961'] = TT961_list[index]
    cc['Qs'] = Qs_list[index]
    cc['R'] = R_list[index]
    cc['Kv'] = Kv_list[index]

dd['13R4']['first_element'] = 'Q1'
dd['33L5']['first_element'] = 'Q1'
dd['13L5']['first_element'] = 'D4'
dd['31L2']['first_element'] = 'D4'

#intermedia data 12R4-13R4
Q1_Tin_12R4 = ['QBQI_12R4_TT825.POSST']
Q1_Tout_12R4 = ['QQBI_12R4_TT826.POSST', 'QQBI_12R4_TT824.POSST']
D2_Tin_12R4 = Q1_Tout_12R4
D2_Tout_12R4 = ['QBBI_A13R4_TT826.POSST', 'QBBI_A13R4_TT824.POSST']
D3_Tin_12R4 = D2_Tout_12R4
D3_Tout_12R4 = ['QBBI_B13R4_TT826.POSST', 'QBBI_B13R4_TT824.POSST']
D4_Tin_12R4 = D3_Tout_12R4
D4_Tout_12R4 = ['QBQI_13R4_TT825.POSST']
# Note the swapped naming conventions!
QBS_12R4 = ['QRLAA_13L5_QBS943_Q1.POSST', 'QRLAA_13L5_QBS943_D2.POSST', 'QRLAA_13L5_QBS943_D3.POSST', 'QRLAA_13L5_QBS943_D4.POSST']

dd['13R4'].update({
    'Q1': {
        'Tin': Q1_Tin_12R4,
        'Tout': Q1_Tout_12R4,
    },
    'D2': {
        'Tin': D2_Tin_12R4,
        'Tout': D2_Tout_12R4,
    },
    'D3': {
        'Tin': D3_Tin_12R4,
        'Tout': D3_Tout_12R4,
    },
    'D4': {
        'Tin': D4_Tin_12R4,
        'Tout': D4_Tout_12R4,
    },
})
dd['13R4']['QBS'] = QBS_12R4


#intermedia data 34R4-33L5
#broken sensor
Q1_Tin_32R4 = ['QBQI_34R4_TT825.POSST']
Q1_Tout_32R4 =['QQBI_34L5_TT826.POSST', 'QQBI_34L5_TT824.POSST']
D2_Tin_32R4 = Q1_Tout_32R4
D2_Tout_32R4 = ['QBBI_B34L5_TT826.POSST', 'QBBI_B34L5_TT824.POSST']
D3_Tin_32R4 = D2_Tout_32R4
# Update in 2017: Sensor QBBI_A34L5_TT826.POSST is now working
D3_Tout_32R4 = ['QBBI_A34L5_TT826.POSST', 'QBBI_A34L5_TT824.POSST']
D4_Tin_32R4 = D3_Tout_32R4
D4_Tout_32R4 = ['QBQI_34L5_TT825.POSST']
QBS_32R4 = ['QRLAA_33L5_QBS947_Q1.POSST', 'QRLAA_33L5_QBS947_D2.POSST', 'QRLAA_33L5_QBS947_D3.POSST', 'QRLAA_33L5_QBS947_D4.POSST']

dd['33L5'].update({
    'Q1': {
        'Tin': Q1_Tin_32R4,
        'Tout': Q1_Tout_32R4,
    },
    'D2': {
        'Tin': D2_Tin_32R4,
        'Tout': D2_Tout_32R4,
    },
    'D3': {
        'Tin': D3_Tin_32R4,
        'Tout': D3_Tout_32R4,
    },
    'D4': {
        'Tin': D4_Tin_32R4,
        'Tout': D4_Tout_32R4,
    },
})
dd['33L5']['QBS'] = QBS_32R4


#intermedia data 14L5-13L5
# This is the cell with the reversed gas flow
D4_Tin_13L5 = ['QBQI_13L5_TT825.POSST']
D4_Tout_13L5 =  ['QBBI_A13L5_TT826.POSST', 'QBBI_A13L5_TT824.POSST']
D3_Tin_13L5 = D4_Tout_13L5
D3_Tout_13L5 = ['QBBI_B13L5_TT826.POSST']#, 'QBBI_B13L5_TT824.POSST']
print('WARNING: Wrong for old fills: Functioning sensor QBBI_B13L5_TT824.POSST is ignored!')
D2_Tin_13L5 = D3_Tout_13L5
D2_Tout_13L5 = ['QQBI_13L5_TT826.POSST', 'QQBI_13L5_TT824.POSST']
Q1_Tin_13L5 = D2_Tout_13L5
Q1_Tout_13L5 = ['QBQI_14L5_TT825.POSST']
# Note the swapped naming conventions!
QBS_13L5 = ['QRLAA_13R4_QBS947_D2.POSST', 'QRLAA_13R4_QBS947_D3.POSST', 'QRLAA_13R4_QBS947_D4.POSST', 'QRLAA_13R4_QBS947_Q1.POSST']

dd['13L5'].update({
    'Q1': {
        'Tin': Q1_Tin_13L5,
        'Tout': Q1_Tout_13L5,
    },
    'D2': {
        'Tin': D2_Tin_13L5,
        'Tout': D2_Tout_13L5,
    },
    'D3': {
        'Tin': D3_Tin_13L5,
        'Tout': D3_Tout_13L5,
    },
    'D4': {
        'Tin': D4_Tin_13L5,
        'Tout': D4_Tout_13L5,
    },
})
dd['13L5']['QBS'] = QBS_13L5


#intermediate data 31L2-32L2
# This is the new cell
D4_Tin_32L2 = ['QBQI_31L2_TT825.POSST']
D4_Tout_32L2 =  ['QBBI_A31L2_TT826.POSST', 'QBBI_A31L2_TT824.POSST']
D3_Tin_32L2 = D4_Tout_32L2
D3_Tout_32L2 = ['QBBI_B31L2_TT826.POSST', 'QBBI_B31L2_TT824.POSST']
D2_Tin_32L2 = D3_Tout_32L2
D2_Tout_32L2 = ['QQBI_31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST']
Q1_Tin_32L2 = D2_Tout_32L2
Q1_Tout_32L2 = ['QBQI_32L2_TT825.POSST']
QBS_32L2 = ['QRLAB_31L2_QBS943_D2.POSST', 'QRLAB_31L2_QBS943_D3.POSST', 'QRLAB_31L2_QBS943_D4.POSST', 'QRLAB_31L2_QBS943_Q1.POSST']

dd['31L2'].update({
    'Q1': {
        'Tin': Q1_Tin_32L2,
        'Tout': Q1_Tout_32L2,
    },
    'D2': {
        'Tin': D2_Tin_32L2,
        'Tout': D2_Tout_32L2,
    },
    'D3': {
        'Tin': D3_Tin_32L2,
        'Tout': D3_Tout_32L2,
    },
    'D4': {
        'Tin': D4_Tin_32L2,
        'Tout': D4_Tout_32L2,
    },
})
dd['31L2']['QBS'] = QBS_32L2

cell_timber_vars_dict = dd
del dd

