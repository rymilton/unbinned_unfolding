from test_utils import perform_test, get_combination

if __name__ == '__main__':
    parms = {
        'method':['1', '2', '4'],
        'ftrainx': ['7'],
        'doerror':['2'],
        'dosys': ['1'],
        'verbose': ['3'],
        'seed':['42']
    }

    parms_name = list(parms.keys())
    parms_name.sort()
    combined_parm = get_combination(parms, parms_name)

    parms = {
        'method':['2'],
        'ftrainx': ['7'],
        'doerror':['3'],
        'dosys': ['1'],
        'ntoys':['50', '500' ,'5000'],
        'verbose': ['3'],
        'seed':['42']
    }
    
    parms_name = list(parms.keys())
    parms_name.sort()
    combined_parm.extend(get_combination(parms, parms_name))

    ref_file_name = "../ref/test_uncertainty.ref"
    test_name = 'test_uncertainty'
    field_to_compare = ['uncertainty']
    perform_test(combined_parm, ref_file_name, test_name, field_to_compare, is_combined = True)