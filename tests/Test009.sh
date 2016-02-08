#!/bin/bash
TESTN=Test009
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff -I '#' $OUTDIR/Test006/gf_entropy10.xvg $DATADIR/test.gf_entropy.xvg > $OUTDIR/$TESTN/$TESTN.log
