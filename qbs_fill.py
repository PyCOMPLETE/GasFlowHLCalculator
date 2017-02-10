import os
import numpy as np

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.LHC_Heatloads as HL
import h5_storage
from config_qbs import config_qbs, arc_index, arc_list
from compute_QBS_special import compute_qbs_special
import compute_QBS_LHC as cql

version = h5_storage.version

# Load data for one fill
def compute_qbs_fill(filln, use_dP=True, version=version, recompute_if_missing=False):
    """
    Arguments:
        -filln
        -use_dP = True
        -version = h5_storage.version
        -force_recompute = False
    """
    if use_dP:
        h5_file = h5_storage.get_qbs_file(filln, version)
        if os.path.isfile(h5_file):
            return h5_storage.load_qbs(filln, version=version)

    if not recompute_if_missing:
        raise ValueError('File does not exist. Set the recompute_if_missing flag!')

    atd_ob = h5_storage.load_data_file(filln)
    qbs_ob = cql.compute_qbs(atd_ob, use_dP, version=version)
    if use_dP:
        h5_storage.store_qbs(filln, qbs_ob, use_dP, version=version)
        print('Stored h5 for fill %i.' % filln)
    return qbs_ob

def test_compute_qbs(filln, use_dP=True, version=version):
    """
    Never loads or saves recomputed data.
    """
    atd_ob = h5_storage.load_data_file(filln)
    return cql.compute_qbs(atd_ob, use_dP, version=version)

# Special cells
def special_qbs_fill(filln):
    atd_ob = h5_storage.load_special_data_file(filln)
    return compute_qbs_special(atd_ob)

# Compute average per ARC
def compute_qbs_arc_avg(qbs_ob):
    qbs_arc_avg = np.zeros((len(qbs_ob.timestamps),8), dtype=float)
    for k in xrange(8):
        first, last = arc_index[k,:]
        qbs_arc_avg[:,k] = np.nanmean(qbs_ob.data[:,first:last+1], axis=1)
    return tm.AlignedTimberData(qbs_ob.timestamps, qbs_arc_avg, arc_list)

# plug-in replacement of old heat load procedure, the fill dict
def get_fill_dict(filln_or_obj, version=version):
    if isinstance(filln_or_obj, int):
        qbs_ob = compute_qbs_fill(filln_or_obj, version=version)
    else:
        qbs_ob = filln_or_obj
    qbs_arc_avg = compute_qbs_arc_avg(qbs_ob)

    # Fill dict consists of timber_variable_list objects
    def make_timber_variable_list(values):
        tvl = tm.timber_variable_list()
        tvl.t_stamps = qbs_ob.timestamps
        tvl.ms = np.zeros_like(tvl.t_stamps)
        tvl.values = values
        return tvl

    fill_dict = {}
    # Arcs
    for arc_ctr, arc in enumerate(arc_list):
        key = '%s_QBS_AVG_ARC.POSST' % arc
        values = qbs_arc_avg.dictionary[arc]
        fill_dict[key] = make_timber_variable_list(values)

    # Quads
    quad_key_list = []
    for key, list_ in HL.variable_lists_heatloads.iteritems():
        if key[:2] in ('Q6', 'Q5'):
            quad_key_list.extend(list_)

    for ctr, cell in enumerate(qbs_ob.variables):
        if cell[:2] in ('05', '06'):
            cell = cell[:4]
            for quad_key in quad_key_list:
                if quad_key[6:10] == cell:
                    fill_dict[quad_key] = make_timber_variable_list(qbs_ob.data[:,ctr])
                    break

    return fill_dict

def lhc_histograms(qbs_ob, avg_time, avg_pm, in_hrs=True):
    """
    Returns data for histograms, not histograms itself!
    """
    if in_hrs:
        qbs_tt = (qbs_ob.timestamps - qbs_ob.timestamps[0])/3600.
    else:
        qbs_tt = qbs_ob.timestamps
    mask_mean = np.abs(qbs_tt - avg_time) < avg_pm
    if sum(mask_mean) == 0:
        raise ValueError('No valid timestamps')
    lhc_hist_dict = {}
    lhc_hist_dict['arcs'] = arc_hist_dict = {}
    varlist = []
    for ctr, arc in enumerate(arc_list):
        first, last = arc_index[ctr,:]
        cell_names = config_qbs.Cell_list[first:last+1]
        mean = np.nanmean(qbs_ob.data[mask_mean,first:last+1], axis=0)
        mask_nan = np.logical_not(np.isnan(mean))
        arc_hist_dict[arc] = mean[mask_nan]
        varlist.extend(np.array(cell_names)[mask_nan])
        if ctr == 0:
            arc_hist_total = arc_hist_dict[arc]
        else:
            arc_hist_total = np.append(arc_hist_total, arc_hist_dict[arc])
    lhc_hist_dict['total'] = arc_hist_total
    lhc_hist_dict['variables'] = varlist
    return lhc_hist_dict

def lhc_arcs(qbs_ob):
    lhc_hl_dict = {}
    for arc_ctr, arc in enumerate(arc_list):
        first, last = arc_index[arc_ctr,:]
        data = qbs_ob.data[:,first:last+1]
        variables = qbs_ob.variables[first:last+1]
        lhc_hl_dict[arc] = tm.AlignedTimberData(qbs_ob.timestamps, data, variables)
    return lhc_hl_dict
