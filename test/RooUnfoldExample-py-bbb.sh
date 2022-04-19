#!/bin/bash
outfile=RooUnfoldExample.py.bbb.ref
python -u examples/RooUnfoldExample.py bbb > $outfile
bash ref/cleanup.sh $outfile
diff $outfile ref/$outfile
