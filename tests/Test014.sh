#!/bin/bash
TESTN=Test014
DATADIR=data
OUTDIR=output/$TESTN
SRCDIR=../src

$SRCDIR/g_sa_analyze -sa $DATADIR/test.lf_str.out -xpmout $OUTDIR/lf_analyze.xpm -H $OUTDIR/lf_entropy_analyze.xvg -pro $OUTDIR/lf_prof_analyze.dat -trans $OUTDIR/lf_transmat_analyze.out -xpm -profile -entropy -trmat >& $OUTDIR/$TESTN.log
