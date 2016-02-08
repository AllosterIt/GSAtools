#!/bin/bash
TESTN=Test019
DATADIR=data
OUTDIR=output/$TESTN
SRCDIR=../src

$SRCDIR/g_sa_analyze -sa $DATADIR/test.lf.80ns.out -fvalue $DATADIR/test.80ns.pc1.xvg -n $DATADIR/test.R80.ndx -MIout $OUTDIR/lf_MI.80ns.out -MIxvg $OUTDIR/lf_MI.80ns.xvg -MIlog $OUTDIR/lf_MI.80ns.log >& $OUTDIR/$TESTN.log
