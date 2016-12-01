import sys
import cPickle
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('..')
import LHCMeasurementTools.LHC_Heatloads as HL
import LHCMeasurementTools.LHC_Energy as Energy
import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.TimberManager as tm
from LHCMeasurementTools.SetOfHomogeneousVariables import SetOfHomogeneousNumericVariables
from LHCMeasurementTools.LHC_BCT import BCT
myfontsz = 16
ms.mystyle_arial(fontsz=myfontsz, dist_tick_lab=8)

filln = 5416
colstr = {1: 'b', 2:'r'}

arc_keys_list = HL.variable_lists_heatloads['AVG_ARC']
beams_list = [1, 2]
with open('../fills_and_bmodes.pkl', 'rb') as fid:
    dict_fill_bmodes = cPickle.load(fid)
t_ref = dict_fill_bmodes[filln]['t_startfill']
tref_string = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t_ref))

fill_dict = {}
fill_dict.update(tm.parse_timber_file('../fill_basic_data_csvs/basic_data_fill_%d.csv' % filln, verbose=False))
fill_dict.update(tm.parse_timber_file('../fill_heatload_data_csvs/heatloads_fill_%d.csv' % filln, verbose=False))
fill_dict.update(tm.parse_timber_file('../fill_bunchbybunch_data_csvs/bunchbybunch_data_fill_%d.csv' % filln, verbose=False))

heatloads = SetOfHomogeneousNumericVariables(variable_list=arc_keys_list, timber_variables=fill_dict)
energy = Energy.energy(fill_dict, beam=1)
bct_bx = {}
for beam_n in beams_list:
    bct_bx[beam_n] = BCT(fill_dict, beam=beam_n)

plt.close('all')
fig = plt.figure()

fig.patch.set_facecolor('w')
fig.canvas.set_window_title('LHC Arcs')
fig.set_size_inches(15., 8.)

plt.suptitle(' Fill. %d started on %s' % (filln, tref_string), fontsize=25)
plt.subplots_adjust(right=0.7, wspace=0.30)

# Intensity and Energy
sptotint = plt.subplot(2, 1, 1)
sptotint.set_ylabel('Total intensity [p+]')
sptotint.grid('on')
for beam_n in beams_list:
    sptotint.plot((bct_bx[beam_n].t_stamps-t_ref)/3600., bct_bx[beam_n].values, '-', color=colstr[beam_n])

spenergy = sptotint.twinx()
spenergy.plot((energy.t_stamps-t_ref)/3600., energy.energy/1e3, c='black', lw=2.)  # alpha=0.1)
spenergy.set_ylabel('Energy [TeV]')
spenergy.set_ylim(0, 7)

# Cell heat loads
sphlcell = plt.subplot(2, 1, 2, sharex=sptotint)
sphlcell.set_ylabel('Heat load [W]')
sphlcell.grid('on')

# Heat loads arcs
arc_ctr = 0
arc_keys_list.sort()
for output_key in arc_keys_list:
    key = output_key
    sp = sphlcell
    color = ms.colorprog(arc_ctr, len(arc_keys_list)+1)
    arc_ctr += 1
    xx_time = (heatloads.timber_variables[key].t_stamps-t_ref)/3600.
    yy_heatloads = (heatloads.timber_variables[key].values)
    
    sp.plot(xx_time, yy_heatloads, '-', lw=2., label=output_key, color=color)

from compute_QBS_LHC import QBS_ARC_AVG, timber_data
t = timber_data.timestamps
for arc_ctr, name in zip(xrange(8),['ARC12','ARC23','ARC34','ARC45','ARC56','ARC67','ARC78','ARC81']):
    color = ms.colorprog(arc_ctr, len(arc_keys_list)+1)
    sp.plot(t/3600,QBS_ARC_AVG[:,arc_ctr], label=name, color=color, ls='--')

# Heat loads model
#    sphlcell.plot(imp_data[:,0], imp_data[:,1]*len_cell, label=imp_label)
#    sphlcell.plot(sr_data[:,0], sr_data[:,1]*len_cell, label=sr_label)
#    sphlcell.plot(imp_data[:,0], (imp_data[:,1]+sr_data[:,1])*len_cell, label='Imp+SR')

sphlcell.legend(prop={'size': myfontsz}, bbox_to_anchor=(1.1, 1),  loc='upper left')


plt.show()
