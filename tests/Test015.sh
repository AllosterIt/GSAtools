#!/bin/bash
TESTN=Test015
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test014/lf_prof_analyze.dat $DATADIR/test.lf_prof.dat > $OUTDIR/$TESTN/$TESTN.log
