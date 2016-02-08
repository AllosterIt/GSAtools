#!/bin/bash
DATADIR=.
OUTDIR=.
INPREFIX=T00
OUTPREFIX=T05
SRCDIR=../../src

######################################################################
#
#       check GSATools have been compiled
/bin/bash ../00_input_data/00_check_bin.sh || exit

echo $SRCDIR/g_sa_analyze -sa $DATADIR/$INPREFIX.lf_str.out\
                          -fvalue $DATADIR/$INPREFIX.PC1.xvg\
                          -n $DATADIR/$INPREFIX.F80.ndx\
                          -MIout $OUTDIR/$OUTPREFIX.lf_MI.F80-PC1.out\
                          -MIxvg $OUTDIR/$OUTPREFIX.lf_MI.F80-PC1.xvg\
                          -MIlog $OUTDIR/$OUTPREFIX.lf_MI.F80-PC1.log

$SRCDIR/g_sa_analyze -sa $DATADIR/$INPREFIX.lf_str.out\
                     -fvalue $DATADIR/$INPREFIX.PC1.xvg\
                     -n $DATADIR/$INPREFIX.F80.ndx\
                     -MIout $OUTDIR/$OUTPREFIX.lf_MI.F80-PC1.out\
                     -MIxvg $OUTDIR/$OUTPREFIX.lf_MI.F80-PC1.xvg\
                     -MIlog $OUTDIR/$OUTPREFIX.lf_MI.F80-PC1.log >& $OUTDIR/$OUTPREFIX.stdout

echo
echo R CMD BATCH 05_plot_MI_over_time.R
R CMD BATCH 05_plot_MI_over_time.R
echo
echo R CMD BATCH 05_plot_MI_SA-PC1.R
R CMD BATCH 05_plot_MI_SA-PC1.R
