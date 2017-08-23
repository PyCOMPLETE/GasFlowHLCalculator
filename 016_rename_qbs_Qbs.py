import os
import h5_storage
import h5py
import numpy as np

#import LHCMeasurementTools.myfilemanager as mfm
#import LHCMeasurementTools.TimberManager as tm

dir_ = os.path.dirname(h5_storage.get_special_qbs_file(0))
h5_files = filter(lambda x: x.endswith('.h5'), os.listdir(dir_))
h5_files = map(lambda x: dir_ + '/' + x, h5_files)


# Rename 'qbs' to 'Qbs'

#for h5_file in h5_files:
#    filln = int(h5_file[-7:-3])
#    if filln < 4485: continue
#    print filln
#    try:
#        with h5py.File(h5_file, 'r+') as open_file:
#            print open_file.keys()
#            vars_ = map(str, open_file['variables'])
#            if 'qbs' in ''.join(vars_):
#                vars_ = map(lambda x: x.replace('qbs', 'Qbs'),vars_)
#                del open_file['variables']
#                open_file.create_dataset('variables', data=np.array(vars_))
#    except IOError:
#        os.remove(h5_file)
#        print 'Removing ' + h5_file

#    try:
#        qbs_ob = mfm.h5_to_obj(h5_file)
#    except:
#        os.remove(h5_file)
#        print 'Removing ' + h5_file
#        continue
#    if 'qbs' in ''.join(qbs_ob.variables):
#        print 'Fixing ' + h5_file
#        vars_ = map(lambda x: x.replace('qbs', 'Qbs'),qbs_ob.variables)
#        #print qbs_ob.variables, vars_
#
#        atd = tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, vars_)
#        h5_storage.store_special_qbs(filln, atd, special_version=1)



# Rename 05R4 etc

n_tries = 5
for h5_file in h5_files:
    n_try = 0
    while n_try < 5:
        try:
            n_tries += 1
            with h5py.File(h5_file, 'r+') as open_file:
                print open_file.keys()
                vars_ = map(str, open_file['variables'])
                old_vars = vars_[:]
                for ii, var in enumerate(old_vars):
                    if var in ('05R4_947', '05L4_947_2'):
                        vars_[ii] += '_comb'
                    elif var in ('05R4_947_2', '05L4_947'):
                        vars_[ii] += '_quad'

                print(filter(lambda x: x.startswith('05R4') or x.startswith('05L4'), old_vars))
                print(filter(lambda x: x.startswith('05R4') or x.startswith('05L4'), vars_))
                import pdb; pdb.set_trace()
                open_file.create_dataset('variables', data=np.array(vars_))
        except IOError:
            print 'Fail for %s' % h5_file
            pass
    else:
        os.remove(h5_file)
        print 'Removing %s' % h5_file

