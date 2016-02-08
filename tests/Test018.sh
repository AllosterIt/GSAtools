#!/bin/bash
TESTN=Test018
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test014/lf_transmat_analyze.out  $DATADIR/test.lf_transmat.out > $OUTDIR/$TESTN/$TESTN.log
