import sys
if '..' not in sys.path: sys.path.append('..')
import operator

import numpy as np

import qbs_fill as qf
import hl_dicts.dict_utils as du

main_dict = du.load_dict('../hl_dicts/test.pkl')

versions = (6,7)
filln = 5219
moment = 'stable_beams'

qbs_obs = map(lambda v: qf.compute_qbs_fill(5219,version=v), versions)

average_arcs = map(qf.compute_qbs_arc_avg, qbs_obs)

dicts = [aa.dictionary for aa in average_arcs]

diff = du.operate_on_dicts(*dicts, operator=operator.sub)

index = np.argwhere(main_dict['filln'] == filln)
tt = main_dict[moment]['t_stamps'][index]


index_tt = np.argmin(np.abs(average_arcs[0].timestamps - tt))

for average_arc in average_arcs:
    print tt
    print average_arc.timestamps[index_tt]
    for ctr, (arc, arr) in enumerate(sorted(average_arc.dictionary.items())):
        hl_dict = main_dict[moment]['heat_load']['arc_averages'][arc][index] + main_dict['hl_subtracted_offset']['arc_averages'][arc][index]
        print ('hl_dict', arc, hl_dict)
        print arr[index_tt-1], arr[index_tt], arr[index_tt+1]
        best_index = np.argmin(np.abs(arr - hl_dict))
        print(arr[best_index-1], arr[best_index+1], arr[best_index])

    print

