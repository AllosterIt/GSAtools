#!/bin/bash
TESTN=Test027
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test023/lf_nMImat.out $DATADIR/test.lf_nMImat.out > $OUTDIR/$TESTN/$TESTN.log
