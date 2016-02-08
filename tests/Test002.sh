#!/bin/bash
TESTN=Test002
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test001/lf_str.out $DATADIR/test.lf_str.out > $OUTDIR/$TESTN/$TESTN.log
