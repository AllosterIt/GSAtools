#!/bin/bash
DATADIR=../00_input_data
OUTDIR=.
INPREFIX=T00
OUTPREFIX=T01
SRCDIR=../../src

######################################################################
#
#       check GSATools have been compiled
/bin/bash ../00_input_data/00_check_bin.sh || exit

######################################################################
#
#       check GMX environment variables
if [[ "$GMXLIB" == "" ]]; then
        source ../../scripts/GSATOOLSRC.bash

        echo "--------------------------------------------------"
        echo " To run this script  a GSATOOLSRC script has been "
        echo " sourced to set the correct GROMACS environmental "
        echo " variables.                                       "
        echo
        echo " Please remember to source the GSATOOLSRC before  " 
        echo " using the g_sa_encode and/or g_sa_analyze.       " 
        echo " You can find the script suitable for your shell  " 
        echo " in the scripts directory.                        " 
        echo
fi

echo $SRCDIR/g_sa_encode -f $DATADIR/$INPREFIX.xtc\
                         -s $DATADIR/$INPREFIX.tpr\
                         -strlf $OUTDIR/$OUTPREFIX.lf_str.out\
                         -rmsdlf $OUTDIR/$OUTPREFIX.lf_rmsd.xvg\
                         -xpmlf $OUTDIR/$OUTPREFIX.lf.xpm\
                         -fasta -xpm -log $OUTDIR/$OUTPREFIX.log

$SRCDIR/g_sa_encode -f $DATADIR/$INPREFIX.xtc\
                    -s $DATADIR/$INPREFIX.tpr\
                    -strlf $OUTDIR/$OUTPREFIX.lf_str.out\
                    -rmsdlf $OUTDIR/$OUTPREFIX.lf_rmsd.xvg\
                    -xpmlf $OUTDIR/$OUTPREFIX.lf.xpm\
                    -fasta -xpm -log $OUTDIR/$OUTPREFIX.log >& $OUTDIR/$OUTPREFIX.stdout

echo
echo xpm2ps -f $OUTDIR/$OUTPREFIX.lf.xpm\
            -o $OUTDIR/$OUTPREFIX.lf.eps -legend none\
            -di 01_ps.m2p

xpm2ps -f $OUTDIR/$OUTPREFIX.lf.xpm\
       -o $OUTDIR/$OUTPREFIX.lf.eps -legend none\
       -di 01_ps.m2p >& $OUTDIR/$OUTPREFIX.eps.stdout

PS2PDF=`which ps2pdf`
if [ -n "$PS2PDF" ]
then
        echo
        echo ps2pdf $OUTDIR/$OUTPREFIX.lf.eps
        ps2pdf $OUTDIR/$OUTPREFIX.lf.eps
fi

echo
echo R CMD BATCH 01_plot_RMSD.R
R CMD BATCH 01_plot_RMSD.R
