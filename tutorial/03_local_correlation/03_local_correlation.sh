#!/bin/bash
DATADIR=../02_SA_statistics
OUTDIR=.
INPREFIX=T00
OUTPREFIX=T03
SRCDIR=../../src

######################################################################
#
#       check GSATools have been compiled
/bin/bash ../00_input_data/00_check_bin.sh || exit

echo $SRCDIR/g_sa_analyze -sa $DATADIR/$INPREFIX.lf_str.out\
                          -MImat $OUTDIR/$OUTPREFIX.lf_MImat.out\
                          -eeMImat $OUTDIR/$OUTPREFIX.lf_eeMImat.out\
                          -jHmat $OUTDIR/$OUTPREFIX.lf_jHmat.out\
                          -nMImat $OUTDIR/$OUTPREFIX.lf_nMImat.out\
                          -MImatrix

echo
echo "Calculation started. Please be patient..."

$SRCDIR/g_sa_analyze -sa $DATADIR/$INPREFIX.lf_str.out\
                     -MImat $OUTDIR/$OUTPREFIX.lf_MImat.out\
                     -eeMImat $OUTDIR/$OUTPREFIX.lf_eeMImat.out\
                     -jHmat $OUTDIR/$OUTPREFIX.lf_jHmat.out\
                     -nMImat $OUTDIR/$OUTPREFIX.lf_nMImat.out\
                     -MImatrix >& $OUTDIR/$OUTPREFIX.stdout

echo
echo R CMD BATCH 03_plot_nMI_matrix.R
R CMD BATCH 03_plot_nMI_matrix.R
