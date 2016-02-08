#!/bin/bash
TESTN=Test001
DATADIR=data
OUTDIR=output/$TESTN
SRCDIR=../src

echo "--------------------------------------------------"
echo " To run these tests a GSATOOLSRC script has been  "
echo " sourced to set the correct GROMACS environmental "
echo " variables.                                       "
echo
echo " Please remember to source the GSATOOLSRC before  " 
echo " using the g_sa_encode and/or g_sa_analyze.       " 
echo " You can find the script suitable for your shell  " 
echo " in the scripts directory.                        " 
echo

$SRCDIR/g_sa_encode -f $DATADIR/test.xtc -s $DATADIR/test.tpr -strlf $OUTDIR/lf_str.out -rmsdlf $OUTDIR/lf_rmsd.xvg -xpmlf $OUTDIR/lf.xpm -Hlf $OUTDIR/lf_entropy.xvg -prolf $OUTDIR/lf_prof.dat -xpm -profile -entropy -log $OUTDIR/log.log >& $OUTDIR/$TESTN.log
