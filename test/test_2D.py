import os
from test_utils import get_field, compare

def delete_files():
    os.system("rm -f RooUnfoldTest2D.root")
    os.system("rm -f RooUnfoldTest2D.ps")

if __name__ == '__main__':

    ref_file_name = "../ref/test_2D.ref"
    test_name = 'test_2D'
    field_to_compare = ['unfold2D']
    all_output = {}
    command_str = "../build/RooUnfoldTest2D  verbose=3"
    os.system(command_str)
    u = get_field("RooUnfoldTest2D.root", field_to_compare)
    all_output['default'] = u
    delete_files()
    if compare(all_output, ref_file_name, test_name, 1) == 1:
        print("Test failed")
        exit(1)
