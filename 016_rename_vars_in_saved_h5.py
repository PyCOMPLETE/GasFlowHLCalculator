import os
import h5_storage
import h5py
import config_qbs as cq
import numpy as np

#import LHCMeasurementTools.myfilemanager as mfm
#import LHCMeasurementTools.TimberManager as tm

dir_special = os.path.dirname(h5_storage.get_special_qbs_file(0))
h5_files_special = filter(lambda x: x.endswith('.h5'), os.listdir(dir_special))
h5_files_special = map(lambda x: dir_special + '/' + x, h5_files_special)

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

for use_dP in (True, False):
    dir_lhc = os.path.dirname(h5_storage.get_qbs_file(0, use_dP))
    h5_files_lhc = filter(lambda x: x.endswith('.h5'), os.listdir(dir_lhc))
    h5_files_lhc = map(lambda x: dir_lhc + '/' + x, h5_files_lhc)
    n_tries_max = 5
    for h5_file in h5_files_lhc:
        filln = int(h5_file[-7:-3])
        if use_dP and filln < 5955: continue
        print h5_file

        n_try = 0
        success = False
        while not success and n_try < n_tries_max:
            try:
                n_try += 1
                with h5py.File(h5_file, 'r+') as open_file:
                    print open_file.keys()
                    vars_ = map(str, open_file['variables'])
                    old_vars = vars_[:]
                    for ii, var in enumerate(old_vars):
                        if var in ('05R4_947', '05L4_947'):
                            if cq.config_qbs.CV94x_list[ii].startswith('QRLEB'):
                                vars_[ii] += '_comb'
                            elif cq.config_qbs.CV94x_list[ii][:5] in ('QRLFE', 'QRLFF'):
                                vars_[ii] += '_quad'
                            else:
                                raise ValueError

                    print(filter(lambda x: x.startswith('05R4') or x.startswith('05L4'), old_vars))
                    print(filter(lambda x: x.startswith('05R4') or x.startswith('05L4'), vars_))
                    del open_file['variables']
                    open_file.create_dataset('variables', data=np.array(vars_))
                    success = True
            except IOError:
                success = False
                print 'Fail for %s' % h5_file
                pass
            except KeyError as e:
                print(e)
                success = False
                break
        if not success:
            os.remove(h5_file)
            print 'Removing %s' % h5_file

