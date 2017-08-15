import numpy as np
import qbs_fill as qf
import compute_QBS_special as cqs

filln = 5219

def compare_dicts_recursively(dd1, dd2):
    equal = True
    for key, value1 in dd1.iteritems():
        value2 = dd2[key]
        if type(value1) is np.ndarray:
            equal = np.all(value1 == value2)
            if not equal:
                print(value1[0], value2[0])
                print key
        elif type(value1) is dict:
            equal = compare_dicts_recursively(value1, value2)
        else:
            equal = value1 == value2
        if not equal:
            break
    return equal

qbs_dict = qf.special_qbs_fill(5219)

qbs_ob = cqs.dict_to_aligned(qbs_dict)

qbs_dict_2 = cqs.aligned_to_dict(qbs_ob)
print(compare_dicts_recursively(qbs_dict_2, qbs_dict))


qbs_dict_2['timestamps'] = np.copy(qbs_dict['timestamps'])
qbs_dict_2['timestamps'][0] = 0

print(compare_dicts_recursively(qbs_dict_2, qbs_dict))
