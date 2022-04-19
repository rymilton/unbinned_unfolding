import os
from test_utils import get_field, compare

def delete_files():
    os.system("rm -f RooUnfoldTest3D.root")
    os.system("rm -f RooUnfoldTest3D.ps")

if __name__ == '__main__':

    ref_file_name = "../ref/test_3D.ref"
    test_name = 'test_3D'
    field_to_compare = ['unfold3D']

    all_output = {}
    delete_files()
    command_str = "../build/RooUnfoldTest3D  verbose=3"
    os.system(command_str)
    u = get_field("RooUnfoldTest3D.root", field_to_compare)
    all_output['default'] = u
    delete_files()
    if compare(all_output, ref_file_name, test_name, 1) == 1:
        print("Test failed")
        exit(1)