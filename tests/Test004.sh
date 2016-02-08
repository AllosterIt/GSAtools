#!/bin/bash
TESTN=Test004
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '#' $OUTDIR/Test001/lf_entropy.xvg $DATADIR/test.lf_entropy.xvg > $OUTDIR/$TESTN/$TESTN.log
