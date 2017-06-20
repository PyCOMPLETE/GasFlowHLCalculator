import os
import numpy as np
import re

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.LHC_Heatloads as HL
import h5_storage
from config_qbs import config_qbs, arc_index, arc_list
from compute_QBS_special import compute_qbs_special, cell_list
import compute_QBS_LHC as cql

default_version = h5_storage.version

# Load data for one fill
def compute_qbs_fill(filln, use_dP=True, version=default_version, recompute_if_missing=False):
    """
    Arguments:
        -filln
        -use_dP = True
        -version = h5_storage.version
        -recompute_if_missing = False
    """
        # Calib is changed for pre LS1 data
    if filln < 3600:
        version = -1
        print 'Warning in GasflowHLCalculator.qbs_fill: special case for pre LS1 fills. Specified version is ignored.'
    else:
        version = default_version

    h5_file = h5_storage.get_qbs_file(filln, version=version, use_dP=use_dP)
    if os.path.isfile(h5_file):
        return h5_storage.load_qbs(filln, version=version, use_dP=use_dP)

    if not recompute_if_missing:
        raise ValueError("""File %s does not exist.
                         Set the correct flag if you want to recompute!""" % h5_file)



    atd_ob = h5_storage.load_data_file(filln)
    qbs_ob = cql.compute_qbs(atd_ob, use_dP, version=version)
    h5_storage.store_qbs(filln, qbs_ob, use_dP, version=version)
    print('Stored h5 for fill %i.' % filln)
    return qbs_ob

def test_compute_qbs(filln, use_dP=True, version=default_version):
    """
    Never loads or saves recomputed data.
    """
    atd_ob = h5_storage.load_data_file(filln)
    return cql.compute_qbs(atd_ob, use_dP, version=version)

# Special cells
def special_qbs_fill(filln, recompute_if_missing=False):

    h5_file = h5_storage.get_special_qbs_file(filln)

    if os.path.isfile(h5_file):
        qbs_ob = h5_storage.load_special_qbs(filln)
        return aligned_to_dict(qbs_ob)
    elif recompute_if_missing:
        new_cell = filln > 5500
        atd_ob = h5_storage.load_special_data_file(filln)
        qbs_dict = compute_qbs_special(atd_ob, new_cell)
        h5_storage.store_special_qbs(dict_to_aligned(qbs_dict))
        print('Stored h5 for fill %i.' % filln)
        return qbs_dict
    else:
        raise ValueError('Set the correct flag if you want to recompute!')

def special_qbs_fill_aligned(filln, recompute_if_missing=False):
    qbs_dict = special_qbs_fill(filln, recompute_if_missing)
    return dict_to_aligned(qbs_dict)

def dict_to_aligned(dict_):
    timestamps = dict_['timestamps']
    variables = []
    data = []
    for cell in dict_['cells']:
        dd = dict_[cell]
        for key, arr in dd.iteritems():
            main_key = cell + '_' + key
            variables.append(main_key)
            data.append(arr)

    data_arr = np.array(data).T
    return tm.AlignedTimberData(timestamps, data_arr, np.array(variables))

def aligned_to_dict(qbs_ob):
    output = {}
    output['timestamps'] = qbs_ob.timestamps
    output['cells'] = cell_list
    for cell in cell_list:
        dd = {}
        output[cell] = dd
        for key in qbs_ob.variables:
            if cell in key:
                subkey = key.split(cell+'_')[1]
                dd[subkey] = qbs_ob.dictionary[key]
    return output

# Compute average per ARC
def compute_qbs_arc_avg(qbs_ob):
    qbs_arc_avg = np.zeros((len(qbs_ob.timestamps),8), dtype=float)
    for k in xrange(8):
        first, last = arc_index[k,:]
        qbs_arc_avg[:,k] = np.nanmean(qbs_ob.data[:,first:last+1], axis=1)
    return tm.AlignedTimberData(qbs_ob.timestamps, qbs_arc_avg, arc_list)

# plug-in replacement of old heat load procedure, the fill dict
def get_fill_dict(filln, version=default_version, use_dP=True):
    qbs_ob = compute_qbs_fill(filln, version=version, use_dP=use_dP)
    qbs_special = special_qbs_fill_aligned(filln)

    # arcs
    qbs_arc_avg = compute_qbs_arc_avg(qbs_ob)
    output = {}
    for arc_ctr, arc in enumerate(arc_list):
        key = '%s_QBS_AVG_ARC.POSST' % arc
        tvl = tm.timber_variable_list()
        tvl.t_stamps = qbs_arc_avg.timestamps
        tvl.ms = np.zeros_like(tvl.t_stamps)
        tvl.values = qbs_arc_avg.dictionary[arc]
        output[key] = tvl

    #others
    varlist_tmb = []
    for kk in HL.variable_lists_heatloads.keys():
        varlist_tmb+=HL.variable_lists_heatloads[kk]

    varlist_tmb+=HL.arcs_varnames_static

    # Remove new special instrumented cell for older fills
    if filln <= 5456:
        regex = re.compile('^QRLAB_31L2_QBS943_\w\w.POSST$')
        varlist_tmb = filter(lambda x: regex.match(x) == None, varlist_tmb)

    for varname in varlist_tmb:
        tvl = tm.timber_variable_list()
        special_id = varname.split('.POSST')[0][-3:]
        if special_id in('_Q1', '_D2', '_D3', '_D4'):
            cell = varname.split('_')[1]
            try:
                tvl.values = qbs_special.dictionary[cell+special_id]
            except:
                import pdb ; pdb.set_trace()
            tvl.t_stamps = qbs_special.timestamps
        elif '_QBS9' in varname:
            firstp, lastp = tuple(varname.split('_QBS'))
            kkk = firstp.split('_')[-1]+'_'+lastp.split('.')[0]
            tvl.t_stamps = qbs_ob.timestamps
            try:
                tvl.values = qbs_ob.dictionary[kkk]
            except KeyError as err:
                print 'Skipped %s! Got:'%kkk
                print err
                tvl.values = np.zeros_like(tvl.t_stamps)
        elif 'QBS_AVG_ARC' in varname or 'QBS_CALCULATED_ARC' in varname:
            continue
        else:
            print('Variable %s not yet implemented in new fill_dict' % varname)
            continue
        tvl.ms = np.zeros_like(tvl.t_stamps)
        output[varname] = tvl
    return output

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
