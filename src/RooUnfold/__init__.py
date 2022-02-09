import os
import pathlib

def insert_path(variable, path):
    curr_path = os.environ.get(variable, '')
    if path in curr_path: return
    os.environ[variable] = ':'.join([path, curr_path])

DIR = str(pathlib.Path(__file__).parent.resolve())
# for ROOT to find it via gSystem.Load
insert_path('LD_LIBRARY_PATH', DIR)
# for ROOT to find the header files
insert_path('ROOT_INCLUDE_PATH', DIR)

import sys
try:
    import ROOT
except ImportError:
    print("Must have ROOT installed.", file=sys.stderr)
    raise

from ROOT import RooUnfold
# hide away this __init__.py
sys.modules[__name__] = RooUnfold
