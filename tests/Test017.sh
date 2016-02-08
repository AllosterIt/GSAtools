#!/bin/bash
TESTN=Test017
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '/*' $OUTDIR/Test014/lf_analyze.xpm $DATADIR/test.lf.xpm > $OUTDIR/$TESTN/$TESTN.log
