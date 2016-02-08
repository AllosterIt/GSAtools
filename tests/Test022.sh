#!/bin/bash
TESTN=Test022
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test019/lf_MI.80ns.log $DATADIR/test.lf_MI.80ns.log > $OUTDIR/$TESTN/$TESTN.log
