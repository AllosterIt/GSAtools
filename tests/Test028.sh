#!/bin/bash
TESTN=Test028
DATADIR=data
OUTDIR=output/$TESTN
SRCDIR=../src

$SRCDIR/g_sa_analyze -sa $DATADIR/test.lf.80ns.out -MImat $OUTDIR/lf_MImat.out -eeMImat $OUTDIR/lf_eeMImat.out -jHmat $OUTDIR/lf_jHmat.out -nMImat $OUTDIR/lf_nMImat.out -nSample 10 -ZMImat $OUTDIR/lf_ZMImat.out -meanMImat $OUTDIR/lf_meanMImat.out -stdMImat $OUTDIR/lf_stdMImat.out -pvalueMImat $OUTDIR/lf_pvalueMImat.out -MImatrix >& $OUTDIR/$TESTN.log
