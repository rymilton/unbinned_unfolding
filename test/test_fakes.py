# BEGIN ROOUNFOLD COPYRIGHT
# RooUnfold — Unfolding library for particle-physics inverse problems
#
# Copyright © 2021–2025 CERN and the authors’ respective research institutions
# Please refer to the CONTRIBUTORS file for details.
#
# License: BSD-3-Clause
# SPDX-License-Identifier: BSD-3-Clause
#
# END ROOUNFOLD COPYRIGHT

from test_utils import perform_test

if __name__ == "__main__":
    parms = {"method": ["1", "2", "3", "4", "5", "6"], "addfakes": ["1"], "verbose": ["3"]}

    ref_file_name = "../ref/test_fakes.ref"
    test_name = "test_fakes"
    field_to_compare = ["unfold"]
    perform_test(parms, ref_file_name, test_name, field_to_compare)
