#!/bin/bash
TESTN=Test016
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '#' $OUTDIR/Test014/lf_entropy_analyze.xvg $DATADIR/test.lf_entropy.xvg > $OUTDIR/$TESTN/$TESTN.log
