#!/bin/bash
TESTN=Test007
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test006/gf_str10.out $DATADIR/test.gf_str.out > $OUTDIR/$TESTN/$TESTN.log
