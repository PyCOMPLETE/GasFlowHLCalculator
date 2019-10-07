import LHCMeasurementTools.lhc_log_db_query as lldb
import LHCMeasurementTools.TimestampHelpers as th
import LHCMeasurementTools.LHC_Fills as Fills

import pickle

periods = [
'2018_05_28 16:00:00!2018_05_29 01:00:00',
]

pkl_name = '/afs/cern.ch/work/e/ecldcode/heat_load_workspace/heat_load_storage/fills_and_bmodes.pkl'

dict_fill_info = {}
for period in periods:

	t_start_string = period.split('!')[0]
	t_stop_string = period.split('!')[1]

	t_start = th.localtime2unixstamp(t_start_string)
	t_stop = th.localtime2unixstamp(t_stop_string)

	csv_name = 'fills_and_bmodes_temp.csv'

	# Get data from database
	varlist = Fills.get_varlist()
	lldb.dbquery(varlist, t_start, t_stop, csv_name)

	# Make pickle
	dict_fill_info.update(Fills.make_dict(csv_name, t_stop))

with open(pkl_name, 'wb') as fid:
      pickle.dump(dict_fill_info, fid)


