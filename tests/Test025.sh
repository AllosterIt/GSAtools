#!/bin/bash
TESTN=Test025
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test023/lf_eeMImat.out $DATADIR/test.lf_eeMImat.out > $OUTDIR/$TESTN/$TESTN.log
