from __future__ import division
import argparse
import matplotlib.pyplot as plt

import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.savefig as sf

import h5_storage
import compute_QBS_special as cqs
#from data_S45_details import cell_timber_vars_dict

#(30-07-17)
#
# Prolems in fill 5786 with only beam 2:
# cell 13L5: beam 2 shows lower heat load for this fill -> sensors wrong?
# cell 31L2: B1 HL higher for Q1
# cell 31L2: 0 HL for D3
# Matching observations when running this script with fill 5783, only beam 1

# (31-07-17):
#
# --> Solved by excluding the last magnet from the calculation, and also adding 31L2 to the normal gasflow list as the temp sensors are reversed (826 instead of 824).
# To be checked with B.B.

parser = argparse.ArgumentParser()
parser.add_argument('filln', help='LHC fill number', type=int)
parser.add_argument('--savefig', help='Save in file.', action='store_true')
parser.add_argument('--noshow', help='Do not call plt.show', action='store_true')
args = parser.parse_args()



plt.close('all')
ms.mystyle(10)

filln = args.filln

blacklist = [
    ('13L5', 'D3'),
    ('13L5', 'D2'),
    ('33L5', 'D3'),
    ('33L5', 'D4'),
]

new_cell = filln > 5700
atd = h5_storage.load_special_data_file(filln)

special_dict_not_separate = cqs.compute_qbs_special(atd, new_cell, False)
special_dict_separate = cqs.compute_qbs_special(atd, new_cell, True)

tt = (atd.timestamps - atd.timestamps[0])/3600


sp = None
for cell_ctr, cell in enumerate(special_dict_separate['cells']):
    fig = ms.figure(cell)
    sp1 = sp = plt.subplot(2,2,2, sharex=sp)
    sp1.set_title('Cell %s - sum of magnet Fill %i' % (cell, filln))
    sp2 = sp = plt.subplot(2,2,1, sharex=sp)
    sp2.set_title('Cell %s - separate beams Fill %i' % (cell, filln))

    for sp in (sp1, sp2):
        sp.set_xlabel('Time [h]')
        sp.set_ylabel('Heat load [W]')


    label_ctr = 0
    for magnet_ctr, magnet_id in enumerate(cqs.magnet_ids):
        if (cell, magnet_id) in blacklist: continue
        #if (cell_timber_vars_dict[cell]['first_element'], magnet_id) in [('D4', 'Q1') or ('Q1', 'D4')]: continue
        label_ctr += 1
        if label_ctr == 1:
            label_old, label_separate = 'Combined', 'Sum of beams'
        else:
            label_old, label_separate = None, None
        color = ms.colorprog(magnet_ctr, cqs.magnet_ids)
        sp1.plot(tt, special_dict_not_separate[cell][magnet_id], color=color, label=label_old, lw=2)
        sum_magnet = 0
        for beam in (1,2):
            label = magnet_id
            if label_ctr == 1:
                label += ' Beam %i' % beam
            elif beam == 2:
                label = None
            hl_beam = special_dict_separate[cell][magnet_id][beam]
            sum_magnet += hl_beam
            ls = ['-',':'][beam-1]
            sp2.plot(tt, hl_beam, ls=ls, color=color, label=label, lw=2)
        sp1.plot(tt, sum_magnet, color=color, ls=':', label=label_separate, lw=2)

    for sp in sp1, sp2:
        sp.grid(True)
        sp.legend(loc=0)
        sp.set_ylim(-5, None)

if args.savefig:
    for num in plt.get_fignums():
        fig = plt.figure(num)
        plt.suptitle('')
    sf.saveall_pdijksta()


if not args.noshow:
    plt.show()

