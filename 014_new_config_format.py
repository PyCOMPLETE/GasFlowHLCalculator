with open('./LHCCryoHeatLoadCalibration/CryoBeamScreenData.csv') as f:
    f_new = f.readlines()
    f_new = [x.split(',') for x in f_new]

with open('./config_qbs_lhc_3.csv') as f:
    f_old = f.readlines()
    f_old = [x.split() for x in f_old]


print('Different naming:\n')
for i in xrange(len(f_new[0])):
        for j in xrange(len(f_old[0])):
            if f_old[1][j] == f_new[1][i]:
                print f_old[0][j], f_new[0][i]


cell_names = set()
print('\nDuplicated cells in new config:\n')
for i, line in enumerate(f_new):
    cell_name = line[3]
    if cell_name in cell_names:
        print cell_name
    else:
        cell_names.add(cell_name)

