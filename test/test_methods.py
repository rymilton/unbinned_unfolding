from test_utils import perform_test
import os

if __name__ == '__main__':
    parms = {
        'method': ['1', '2', '3', '4', '5', '6'],
        'verbose': ['3']
    }
    ref_file_name = "../ref/test_methods.ref"
    test_name = 'test_methods'
    field_to_compare = ['unfold']
    perform_test(parms, ref_file_name, test_name, field_to_compare)

    ## For plotting
    command_str = "../build/RooUnfoldTest ploterrors=2"
    os.system(command_str)
    command_str = "../build/RooUnfoldTest ploterrors=1"
    os.system(command_str)
    command_str = "../build/RooUnfoldTest plotparms=2"
    os.system(command_str)