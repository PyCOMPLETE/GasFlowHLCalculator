from __future__ import division
import matplotlib.pyplot as plt

import LHCMeasurementTools.mystyle as ms

import h5_storage
import compute_QBS_special as cqs

# Prolems in fill 5786 with only beam 2:
# cell 13L5: beam 2 shows lower heat load for this fill -> sensors wrong?
# cell 31L2: B1 HL higher for Q1
# cell 31L2: 0 HL for D3

# Matching observations when running this script with fill 5783, only beam 1

plt.close('all')
ms.mystyle(10)

filln = 5783

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
    sp1 = sp = plt.subplot(2,2,1, sharex=sp)
    sp1.set_title('Sum of magnet Fill %i' % filln)
    sp2 = sp = plt.subplot(2,2,2, sharex=sp)
    sp2.set_title('Separate beams Fill %i' % filln)


    for magnet_ctr, magnet_id in enumerate(cqs.magnet_ids):
        if (cell, magnet_id) in blacklist: continue
        if magnet_ctr == 0:
            label_old, label_separate = 'Combined', 'Sum of beams'
        else:
            label_old, label_separate = None, None
        color = ms.colorprog(magnet_ctr, cqs.magnet_ids)
        sp1.plot(tt, special_dict_not_separate[cell][magnet_id], color=color, label=label_old)
        sum_magnet = 0
        for beam in (1,2):
            label = magnet_id
            if magnet_ctr == 0:
                label += ' Beam %i' % beam
            hl_beam = special_dict_separate[cell][magnet_id][beam]
            sum_magnet += hl_beam
            ls = ['-',':'][beam-1]
            sp2.plot(tt, hl_beam, ls=ls, color=color, label=label)
        sp1.plot(tt, sum_magnet, color=color, ls=':', label=label_separate)

    for sp in sp1, sp2:
        sp.grid(True)
        sp.legend()
        sp.set_ylim(-5, None)

plt.show()


