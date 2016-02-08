#!/bin/bash
DATADIR=.
OUTDIR=.
GSA=/home/jkleinj/jk.software/develop/SAsuite/GSAtools/gsatools-4.5.x-1.00/src

source $GSA/../scripts/GSATOOLSRC.bash
#srun -J dfb0cor -p long -n 32 -o cor.o -e cor.e\
#               $GSA/g_sa_analyze -sa $DATADIR/lf_str_short.out\

mpirun -n 8 $GSA/g_sa_analyze -sa $DATADIR/lf_str.out\
                     -MImat $OUTDIR/out.lf_MImat.out\
                     -eeMImat $OUTDIR/out.lf_eeMImat.out\
                     -jHmat $OUTDIR/out.lf_jHmat.out\
                     -nMImat $OUTDIR/out.lf_nMImat.out\
		     -nSample 100\
		     -ZMImat $OUTDIR/out.lf_ZImat.out\
		     -meanMImat $OUTDIR/out.lf_meanMImat.out\
		     -stdMImat $OUTDIR/out.lf_stdMImat.out -pvalueMImat $OUTDIR/out.lf_pvalueMImat.out\
                     -MImatrix
