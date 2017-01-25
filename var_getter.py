from __future__ import division, print_function
import numpy as np

class VarGetter(object):
    """
    This class can be used to relate the raw timber data for a given cell or all cells.
    """
    def __init__(self, atd_ob, dq, strict=True, report=False):
        """
        Parameters:
            -atd_ob: Timber data
            -dq    : Data_qbs instance
            -strict: Raise error if there are missing variables. Default: True
            -report: Print information on failed sensors.
        """

        self.atd_ob  = atd_ob
        self.dq      = dq
        self.strict  = strict

        self.Ncell = len(dq.Cell_list)
        self.Nvalue = len(atd_ob.timestamps)

        self.nan_cells = {}
        self.missing_variables   = []
        self.corrected_variables = []
        self.nan_arr = np.zeros(len(dq.Cell_list), dtype=np.bool)

        self._var_data_dict = {
            'T1': {
                'vars': dq.TT961_list,
                'correct': True,
            },
            'T3': {
                'vars': dq.TT94x_list,
                'correct': False,
            },
            'CV': {
                'vars': dq.CV94x_list,
                'correct': False,
            },
            'EH': {
                'vars': dq.EH84x_list,
                'correct': False,
            },
            'P1': {
                'vars': dq.PT961_list,
                'correct': True,
            },
            'P4': {
                'vars': dq.PT991_list,
                'correct': True,
            },
            'T2': {
                'vars': dq.TT84x_list,
                'correct': False,
            },
        }

        self.data_dict = self._store_all_cell_data()
        if report:
            self.report()

    def _insert_to_nan_cells(self, cell_ctr, var):
        cell = self.dq.Cell_list[cell_ctr]
        if cell not in self.nan_cells:
            self.nan_cells[cell] = {
                'ctr': cell_ctr,
                'varlist': [var],
            }
        else:
            self.nan_cells[cell]['varlist'].append(var)

    def _get_value(self, arr, names, correct_zero):
        """
        Inserts values to array arr.
        Parameters:
            - arr: Empty array of shape timestamps * n_cells
            - names: List of Timber variable names
            - correct_zero: If True then the value of the previous cell is taken.
                This is only useful for selected variable types.
        """
        atd_ob = self.atd_ob
        correct_first = False
        for cell_ctr, name in enumerate(names):
            try:
                data = atd_ob.dictionary[name]
            except KeyError:
                self.missing_variables.append(name)
            else:
                arr[:,cell_ctr] = data
                if np.any(arr[:,cell_ctr] <= 0):
                    if correct_zero:
                        if cell_ctr != 0:
                            arr[:,cell_ctr] = arr[:,cell_ctr-1]
                            self.corrected_variables.append(name)
                        else:
                            correct_first = True
                    else:
                        self._insert_to_nan_cells(cell_ctr, name)
                        self.nan_arr[cell_ctr] = True
        if correct_first:
            arr[:,0] = arr[:,-1]

    def _store_all_cell_data(self):

        data_dict = {}
        for key, dd in self._var_data_dict.iteritems():
            data = np.zeros((self.Nvalue, self.Ncell), dtype=float)
            self._get_value(data, dd['vars'], dd['correct'])
            data_dict[key] = data
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
        for key, dd in self._var_data_dict.iteritems():
            name = dd['vars'][index]
            output_dict[key] = self.atd_ob.dictionary[name]
        return output_dict

    def report(self):
        """
        Print out the missing and corrected variables as well as the nan cells.
        """
        if self.missing_variables:
            print('Missing variables:', self.missing_variables)
            if self.strict:
                raise ValueError('There have been missing variables!')

        if self.corrected_variables:
            print('Corrected variables:', self.corrected_variables)

        if self.nan_cells:
            print('\nFailing variables for %i cells:' % len(self.nan_cells))
        for cell, dd in self.nan_cells.iteritems():
            ctr = dd['ctr']
            varlist = dd['varlist']
            string = '%s in S%s:\n' % (cell, self.dq.Sector_list[ctr])
            for var in varlist:
                string += '\t%s' % var
            print(string)

