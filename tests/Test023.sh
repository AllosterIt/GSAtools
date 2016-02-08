#!/bin/bash
TESTN=Test023
DATADIR=data
OUTDIR=output/$TESTN
SRCDIR=../src

$SRCDIR/g_sa_analyze -sa $DATADIR/test.lf.80ns.out -MImat $OUTDIR/lf_MImat.out -eeMImat $OUTDIR/lf_eeMImat.out -jHmat $OUTDIR/lf_jHmat.out -nMImat $OUTDIR/lf_nMImat.out -MImatrix >& $OUTDIR/$TESTN.log
