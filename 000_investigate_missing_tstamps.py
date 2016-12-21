from __future__ import division
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('..')
import LHCMeasurementTools.myfilemanager as mfm


plt.close('all')
filln_list = [5219, 5222, 5223, 5416]
for ctr, filln in enumerate(filln_list):

    root_dir = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/'
    h5_name = 'cryo_data_fill_%i.h5' % filln

    ob = mfm.h5_to_obj(root_dir + h5_name)

    ob.timestamps -= ob.timestamps[0]
    ob.timestamps /= 3600.

    print(len(ob.timestamps), ob.data.shape)


    plt.plot(ob.timestamps, np.ones_like(ob.timestamps)*ctr, marker='o', label=filln)

plt.legend()
plt.show()
