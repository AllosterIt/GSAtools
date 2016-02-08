#!/bin/bash
TESTN=Test006
DATADIR=data
OUTDIR=output/$TESTN
SRCDIR=../src

$SRCDIR/g_sa_encode -f $DATADIR/test.xtc -s $DATADIR/test.tpr -strlf $OUTDIR/lf_str10.out -rmsdlf $OUTDIR/lf_rmsd10.xvg -xpmlf $OUTDIR/lf10.xpm -Hlf $OUTDIR/lf_entropy10.xvg -prolf $OUTDIR/lf_prof10.dat -strgf $OUTDIR/gf_str10.out -rmsdgf $OUTDIR/gf_rmsd10.xvg -xpmgf $OUTDIR/gf10.xpm -Hgf $OUTDIR/gf_entropy10.xvg -progf $OUTDIR/gf_prof10.dat -pdbgf $OUTDIR/gf10.pdb -xpm -profile -entropy -fasta -global -heapsize 50 -log $OUTDIR/log10.log -e 10 >& $OUTDIR/$TESTN.log
