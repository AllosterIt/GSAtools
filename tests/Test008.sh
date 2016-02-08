#!/bin/bash
TESTN=Test008
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test006/gf_prof10.dat $DATADIR/test.gf_prof.dat > $OUTDIR/$TESTN/$TESTN.log
