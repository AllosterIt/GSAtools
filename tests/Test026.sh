#!/bin/bash
TESTN=Test026
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test023/lf_jHmat.out $DATADIR/test.lf_jHmat.out > $OUTDIR/$TESTN/$TESTN.log
