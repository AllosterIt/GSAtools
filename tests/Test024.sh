#!/bin/bash
TESTN=Test024
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test023/lf_MImat.out $DATADIR/test.lf_MImat.out > $OUTDIR/$TESTN/$TESTN.log
