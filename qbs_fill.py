import sys
import os
import h5py
import numpy as np

sys.path.append('..')
import LHCMeasurementTools.myfilemanager as mfm
import h5_storage
from data_QBS_LHC import arc_index, arc_list

version = h5_storage.version
h5_dir = h5_storage.h5_dir

# Load data for one fill
def compute_qbs_fill(filln, use_dP=True, version=version):
    if use_dP:
        h5_file = h5_storage.get_qbs_file(filln, version)
        if os.path.isfile(h5_file):
            return h5_storage.load_qbs(filln, version=version)

    import compute_QBS_LHC as cql
    atd_ob = mfm.h5_to_obj(h5_dir + 'cryo_data_fill_%i.h5' % filln)
    qbs_ob = cql.compute_qbs(atd_ob, use_dP)
    if use_dP:
        h5_storage.store_qbs(filln, qbs_ob, use_dP, version=version)
        print('Stored h5 for fill %i.' % filln)
    return qbs_ob

# Compute average per ARC
def compute_qbs_arc_avg(qbs_ob):
    n_timestamps = qbs_ob.data.shape[0]
    qbs_arc_avg = np.zeros((n_timestamps,8), dtype=float)
    for k in xrange(8):
        first = arc_index[k,0]
        last = arc_index[k,1]
        qbs_arc_avg[:,k] = np.mean(qbs_ob.data[:,first:last+1], axis=1)
    return qbs_arc_avg

# Returns data for histograms, not histograms itself!
def arc_histograms(qbs_ob, avg_time_hrs, avg_pm_hrs):
    qbs_tt = (qbs_ob.timestamps - qbs_ob.timestamps[0])/3600.
    mask_mean = np.abs(qbs_tt - avg_time_hrs) < avg_pm_hrs
    arc_hist_dict = {}
    for ctr, arc_str in enumerate(arc_list):
        arc = arc_str[-2:]
        first, last = arc_index[ctr,:]
        arc_hist_dict[arc] = np.mean(qbs_ob.data[mask_mean,first:last+1], axis=0)
        if ctr == 0:
            arc_hist_total = arc_hist_dict[arc]
        else:
            arc_hist_total = np.append(arc_hist_total, arc_hist_dict[arc])
    return arc_hist_total, arc_hist_dict
