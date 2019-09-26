import os
import numpy as np
import re

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.LHC_Heatloads as HL
from config_qbs import arc_index, arc_list
import compute_QBS_special as cqs
import compute_QBS_LHC as cql

# Load data for one fill
def compute_qbs_fill(filln, h5_storage=None, use_dP=True, recompute_if_missing=False):
    """
    Arguments:
        -filln
        -use_dP = True
        -recompute_if_missing = False
    """
    h5_file = h5_storage.get_qbs_file(filln, use_dP=use_dP)
    if os.path.isfile(h5_file):
        return h5_storage.load_qbs(filln, use_dP=use_dP)

    if not recompute_if_missing:
        raise ValueError("""File %s does not exist.
                         Set the correct flag if you want to recompute!""" % h5_file)

    atd_ob = h5_storage.load_data_file(filln)
    qbs_ob = cql.compute_qbs(atd_ob, use_dP)
    h5_storage.store_qbs(filln, qbs_ob, use_dP)
    print('Stored h5 for fill %i.' % filln)
    return qbs_ob

# def test_compute_qbs(filln, use_dP=True):
#     """
#     Never loads or saves recomputed data.
#     """
#     atd_ob = h5_storage.load_data_file(filln)
#     return cql.compute_qbs(atd_ob, use_dP)

# Special cells
def special_qbs_fill(filln, h5_storage=None, recompute_if_missing=False, force_recalc=False, aligned=False):

    if force_recalc:
        print('Force recalculated')
        new_cell = filln > 5500
        atd_ob = h5_storage.load_special_data_file(filln)
        return cqs.compute_qbs_special(atd_ob, new_cell, aligned=aligned)

    h5_file = h5_storage.get_special_qbs_file(filln)

    if os.path.isfile(h5_file):
        qbs_ob = h5_storage.load_special_qbs(filln)
        cqs.aligned_to_dict_separate(qbs_ob)

    elif recompute_if_missing:


        new_cell = filln > 5500
        atd_ob = h5_storage.load_special_data_file(filln)
        qbs_dict = cqs.compute_qbs_special(atd_ob, new_cell, aligned=False)
        qbs_aligned = cqs.dict_to_aligned_separate(qbs_dict)
        h5_storage.store_special_qbs(filln, qbs_aligned)
        print('Stored h5 for fill %i.' % filln)
        if aligned:
            return qbs_aligned
        else:
            return qbs_dict
    else:
            raise ValueError('Set the correct flag if you want to recompute!')

# Compute average per ARC
def compute_qbs_arc_avg(qbs_ob):
    qbs_arc_avg = np.zeros((len(qbs_ob.timestamps),8), dtype=float)
    for k in xrange(8):
        first, last = arc_index[k,:]
        qbs_arc_avg[:,k] = np.nanmean(qbs_ob.data[:,first:last+1], axis=1)
    return tm.AlignedTimberData(qbs_ob.timestamps, qbs_arc_avg, arc_list)

# plug-in replacement of old heat load procedure, the fill dict
def get_fill_dict(filln, h5_storage=None, use_dP=True):
    qbs_ob = compute_qbs_fill(filln, h5_storage=h5_storage, use_dP=use_dP)
    qbs_special = special_qbs_fill(filln, h5_storage=h5_storage, aligned=True)

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
    varlist_tmb+=HL.other_varnames_static

    # Remove new special instrumented cell for older fills
    if filln <= 5456:
        regex = re.compile('^QRLAB_31L2_QBS943_\w\w.POSST$')
        varlist_tmb = filter(lambda x: regex.match(x) is None, varlist_tmb)

    for varname in varlist_tmb:
        tvl = tm.timber_variable_list()
        special_id = varname.split('.POSST')[0][-3:]
        if special_id in('_Q1', '_D2', '_D3', '_D4'):
            cell = varname.split('_')[1]
            tvl.values = qbs_special.dictionary[cell+special_id]
            tvl.t_stamps = qbs_special.timestamps
            for beam in (1,2):
                tvl2 = tm.timber_variable_list()
                tvl2.values = qbs_special.dictionary[cell+special_id+'_%i' % beam]
                tvl2.t_stamps = qbs_special.timestamps
                tvl2.ms = np.zeros_like(tvl2.t_stamps)
                output[varname+'_B%i' % beam] = tvl2

        elif varname.startswith('QRLEB_05L4'):
            temp_name = '05L4_947_comb'
            tvl.t_stamps = qbs_ob.timestamps
            if temp_name in qbs_ob.dictionary:
                tvl.values = qbs_ob.dictionary[temp_name]
            else:
                tvl.values = tvl.t_stamps*0.
                print 'Skipped %s due to key error!'%temp_name
        elif varname.startswith('QRLEB_05R4'):
            temp_name = '05R4_947_comb'
            tvl.t_stamps = qbs_ob.timestamps
            if temp_name in qbs_ob.dictionary:
                tvl.values = qbs_ob.dictionary[temp_name]
            else:
                tvl.values = tvl.t_stamps*0.
                print 'Skipped %s due to key error!'%temp_name
        elif varname.startswith('QRLFF_05L4'):
            temp_name = '05L4_947_quad'
            tvl.t_stamps = qbs_ob.timestamps
            if temp_name in qbs_ob.dictionary:
                tvl.values = qbs_ob.dictionary[temp_name]
            else:
                tvl.values = tvl.t_stamps*0.
                print 'Skipped %s due to key error!'%temp_name
        elif varname.startswith('QRLFF_05R4'):
            temp_name = '05R4_947_quad'
            tvl.t_stamps = qbs_ob.timestamps
            if temp_name in qbs_ob.dictionary:
                tvl.values = qbs_ob.dictionary[temp_name]
            else:
                tvl.values = tvl.t_stamps*0.
                print 'Skipped %s due to key error!'%temp_name
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

# def lhc_histograms(qbs_ob, avg_time, avg_pm, in_hrs=True):
#     """
#     Returns data for histograms, not histograms itself!
#     """
#     if in_hrs:
#         qbs_tt = (qbs_ob.timestamps - qbs_ob.timestamps[0])/3600.
#     else:
#         qbs_tt = qbs_ob.timestamps
#     mask_mean = np.abs(qbs_tt - avg_time) < avg_pm
#     if sum(mask_mean) == 0:
#         raise ValueError('No valid timestamps')
#     lhc_hist_dict = {}
#     lhc_hist_dict['arcs'] = arc_hist_dict = {}
#     varlist = []
#     for ctr, arc in enumerate(arc_list):
#         first, last = arc_index[ctr,:]
#         cell_names = config_qbs.Cell_list[first:last+1]
#         mean = np.nanmean(qbs_ob.data[mask_mean,first:last+1], axis=0)
#         mask_nan = np.logical_not(np.isnan(mean))
#         arc_hist_dict[arc] = mean[mask_nan]
#         varlist.extend(np.array(cell_names)[mask_nan])
#         if ctr == 0:
#             arc_hist_total = arc_hist_dict[arc]
#         else:
#             arc_hist_total = np.append(arc_hist_total, arc_hist_dict[arc])
#     lhc_hist_dict['total'] = arc_hist_total
#     lhc_hist_dict['variables'] = varlist
#     return lhc_hist_dict
# 
# def lhc_arcs(qbs_ob):
#     lhc_hl_dict = {}
#     for arc_ctr, arc in enumerate(arc_list):
#         first, last = arc_index[arc_ctr,:]
#         data = qbs_ob.data[:,first:last+1]
#         variables = qbs_ob.variables[first:last+1]
#         lhc_hl_dict[arc] = tm.AlignedTimberData(qbs_ob.timestamps, data, variables)
#     return lhc_hl_dict

