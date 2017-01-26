from __future__ import division, print_function
import numpy as np

class VarGetter(object):
    """
    This class can be used to relate the raw timber data for a given cell or all cells.
    Parameters:
        -atd_ob: Timber data
        -dq    : Data_qbs instance
        -strict: Raise error if there are missing variables. Default: True
        -report: Print information on failed sensors.
    """
    def __init__(self, atd_ob, dq, strict=True, report=False):

        self.atd_ob  = atd_ob
        self.dq      = dq
        self.strict  = strict

        self.Ncell = len(dq.Cell_list)
        self.Nvalue = len(atd_ob.timestamps)

        self.problem_cells = {}
        self.missing_variables = []
        self.nan_arr = np.zeros(len(dq.Cell_list), dtype=np.bool)

        # correct: If no data available, copy it from previous cell.
        # negative: This data may be negative.
        self.var_data_dict = {
            'T1': {
                'vars': dq.TT961_list,
                'correct': True,
                'negative': False,
            },
            'T3': {
                'vars': dq.TT94x_list,
                'correct': False,
                'negative': False,
            },
            'CV': {
                'vars': dq.CV94x_list,
                'correct': False,
                'negative': False,
            },
            'EH': {
                'vars': dq.EH84x_list,
                'correct': False,
                'negative': True,
            },
            'P1': {
                'vars': dq.PT961_list,
                'correct': True,
                'negative': False,
            },
            'P4': {
                'vars': dq.PT991_list,
                'correct': True,
                'negative': False,
            },
            'T2': {
                'vars': dq.TT84x_list,
                'correct': False,
                'negative': False,
            },
        }

        self.data_dict = self._store_all_cell_data()
        self.assure()
        if report:
            VarGetter.report(self)

    def _store_all_cell_data(self):
        data_dict = {}
        for key, dd in self.var_data_dict.iteritems():
            arr = np.zeros((self.Nvalue, self.Ncell), dtype=float)
            negative_allowed = dd['negative']

            correct_first = False
            for cell_ctr, var_name in enumerate(dd['vars']):
                try:
                    data = self.atd_ob.dictionary[var_name]
                except KeyError:
                    self.missing_variables.append(var_name)
                    continue

                arr[:,cell_ctr] = data
                if (negative_allowed and np.all(arr[:,cell_ctr] == 0)) or \
                (not negative_allowed and np.all(arr[:,cell_ctr] <= 0)):
                    if dd['correct']:
                        if cell_ctr != 0:
                            arr[:,cell_ctr] = arr[:,cell_ctr-1]
                            self._insert_to_problem_cells(cell_ctr, var_name, 'corrected')
                        else:
                            correct_first = True
                    else:
                        self._insert_to_problem_cells(cell_ctr, var_name, 'no_data')
                        self.nan_arr[cell_ctr] = True
                elif not negative_allowed and np.any(arr[:,cell_ctr] <= 0):
                    self._insert_to_problem_cells(cell_ctr, var_name, 'questionable_data')

            if correct_first:
                arr[:,0] = arr[:,-1]

            data_dict[key] = arr

            if self.missing_variables:
                print('Warning! Some variables are missing!')
                print(self.missing_variables)
                if self.strict:
                    raise ValueError('Missing variables!')
        return data_dict

    def get_single_cell_data(self, cell):
        """
        Returns a dictionary with the data for every variable type
        """
        output_dict = {}
        for index, c in enumerate(self.dq.Cell_list):
            if cell in c:
                break
        else:
            raise ValueError('Cell not found!')
        for key, arr in self.data_dict.iteritems():
            output_dict[key] = arr[:,index]
        return output_dict

    def _insert_to_problem_cells(self, cell_ctr, var, type_):
        problem_cells = self.problem_cells
        cell = self.dq.Cell_list[cell_ctr]
        if type_ not in problem_cells:
            problem_cells[type_] = {}
        if cell not in problem_cells[type_]:
            problem_cells[type_][cell] = {
                'sector': self.dq.Sector_list[cell_ctr],
                'type': self.dq.Type_list[cell_ctr],
                'list': []
            }
        problem_cells[type_][cell]['list'].append(var)

    def assure(self):
        """
        P1 > P4
        """
        P1 = self.data_dict['P1']
        P4 = self.data_dict['P4']
        for cell_ctr, isnan in enumerate(self.nan_arr):
            if not isnan and np.any(P1[:,cell_ctr] < P4[:,cell_ctr]):
                self._insert_to_problem_cells(cell_ctr, 'P1 < P4', 'failed_checks')

    def report(self, details=False):
        """
        Print out the missing and corrected variables as well as the nan cells.
        """
        for type_, dd in self.problem_cells.iteritems():
            print('%i problems of type %s' % (len(dd), type_))
            if details:
                for cell, subdict in dd.iteritems():
                    print('%s in S%s of type %s: %s' % (cell, subdict['sector'], subdict['type'], subdict['list']))

