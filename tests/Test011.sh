#!/bin/bash
TESTN=Test011
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test006/gf10.pdb $DATADIR/test.gf.pdb >> $OUTDIR/$TESTN/$TESTN.log
