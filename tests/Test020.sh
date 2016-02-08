#!/bin/bash
TESTN=Test020
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test019/lf_MI.80ns.out $DATADIR/test.lf_MI.80ns.out > $OUTDIR/$TESTN/$TESTN.log
