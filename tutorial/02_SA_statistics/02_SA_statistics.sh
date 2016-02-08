#!/bin/bash
DATADIR=.
OUTDIR=.
INPREFIX=T00
OUTPREFIX=T02
SRCDIR=../../src

######################################################################
#
#       check GSATools have been compiled
/bin/bash ../00_input_data/00_check_bin.sh || exit

echo $SRCDIR/g_sa_analyze -sa $DATADIR/$INPREFIX.lf_str.out\
                          -H $OUTDIR/$OUTPREFIX.lf_entropy.xvg\
                          -pro $OUTDIR/$OUTPREFIX.lf_prof.dat\
                          -trans $OUTDIR/$OUTPREFIX.lf_transmat.out\
                          -profile -entropy -trmat

$SRCDIR/g_sa_analyze -sa $DATADIR/$INPREFIX.lf_str.out\
                     -H $OUTDIR/$OUTPREFIX.lf_entropy.xvg\
                     -pro $OUTDIR/$OUTPREFIX.lf_prof.dat\
                     -trans $OUTDIR/$OUTPREFIX.lf_transmat.out\
                     -profile -entropy -trmat >& $OUTDIR/$OUTPREFIX.stdout

echo
echo R CMD BATCH 02_plot_entropy.R
R CMD BATCH 02_plot_entropy.R
echo
echo R CMD BATCH 02_plot_profile.R
R CMD BATCH 02_plot_profile.R
echo
echo R CMD BATCH 02_plot_transmat.R
R CMD BATCH 02_plot_transmat.R
