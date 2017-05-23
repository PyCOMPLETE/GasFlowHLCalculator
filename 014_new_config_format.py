
with open('./LHCCryoHeatLoadCalibration/CryoBeamScreenData.csv') as f:
    f_new = f.readlines()
    f_new = map(lambda x: x.split(','), f_new)

with open('./config_qbs_lhc_3.csv') as f:
    f_old = f.readlines()
    f_old = map(lambda x: x.split(), f_old)


for i in xrange(len(f_new[0])):
    #print f_new[0][i], f_new[1][i]
        for j in xrange(len(f_old[0])):
            if f_old[1][j] == f_new[1][i]:
                print f_old[0][j], f_new[0][i]


cell_names = set()
print "duplicated"
for i, line in enumerate(f_new):
    cell_name = line[3]
    if cell_name in cell_names:
        print cell_name
    else:
        cell_names.add(cell_name)

