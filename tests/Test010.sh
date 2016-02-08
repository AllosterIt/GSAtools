#!/bin/bash
TESTN=Test010
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '/*' $OUTDIR/Test006/gf10.xpm $DATADIR/test.gf.xpm >> $OUTDIR/$TESTN/$TESTN.log
