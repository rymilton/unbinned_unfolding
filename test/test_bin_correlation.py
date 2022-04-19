from test_utils import perform_test

if __name__ == '__main__':
    parms = {
        'bincorr': ['1'],
        'seed':['42'],
        'verbose': ['3']
    }
    ref_file_name = "../ref/test_correlation.ref"
    test_name = 'test_correlation'
    field_to_compare = ['uncertainty']
    perform_test(parms, ref_file_name, test_name, field_to_compare)