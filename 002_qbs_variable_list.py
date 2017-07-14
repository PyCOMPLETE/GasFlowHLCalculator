import argparse

import LHCMeasurementTools.myfilemanager as  mfm

#import data_S45_details as dd

parser = argparse.ArgumentParser()
parser.add_argument('--dry', help='Dry run', action='store_true')
args = parser.parse_args()

filln = 5219

h5_filename = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/cryo_data_fill_%i.h5' % filln

atd = mfm.h5_to_obj(h5_filename)

variable_list = atd.variables


variable_list_special = set()

# Do not use this for special cells, as one sensor would be lacking.
# This is the one that has worked last year but not any more this year.
#for attr in dir(dd):
#    if '__' in attr: continue
#    name = getattr(dd, attr)
#    if type(name) is list:
#        for thing in name:
#            if type(thing) is str and 'POSST' in thing:
#                variable_list_special.add(thing)

if not args.dry:
    with open('./variable_list_complete.txt', 'w') as f:
        f.write(','.join(variable_list))
#    with open('./variable_list_special.txt', 'w') as f:
#        f.write(','.join(list(variable_list_special)))

