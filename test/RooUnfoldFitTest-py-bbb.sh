#!/bin/bash
outfile=RooUnfoldFitTest.py.bbb.ref
python -u examples/RooUnfoldFitTest.py bbb --mode RooUnfoldSpec > $outfile
bash ref/cleanup.sh $outfile
diff $outfile ref/$outfile
