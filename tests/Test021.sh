#!/bin/bash
TESTN=Test021
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '#' $OUTDIR/Test019/lf_MI.80ns.xvg $DATADIR/test.lf_MI.80ns.xvg > $OUTDIR/$TESTN/$TESTN.log
