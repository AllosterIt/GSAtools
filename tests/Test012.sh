#!/bin/bash
TESTN=Test012
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test006/gf_str10.fasta $DATADIR/test.gf_str.fasta >> $OUTDIR/$TESTN/$TESTN.log
