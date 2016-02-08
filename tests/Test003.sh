#!/bin/bash
TESTN=Test003
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test001/lf_prof.dat $DATADIR/test.lf_prof.dat > $OUTDIR/$TESTN/$TESTN.log
