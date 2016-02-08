#!/bin/bash
TESTN=Test005
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '/*' $OUTDIR/Test001/lf.xpm $DATADIR/test.lf.xpm > $OUTDIR/$TESTN/$TESTN.log
