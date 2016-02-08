#!/bin/bash
TESTN=Test013
DATADIR=data
OUTDIR=output
SRCDIR=../src

diff $OUTDIR/Test006/lf_str10.fasta $DATADIR/test.lf_str.fasta >> $OUTDIR/$TESTN/$TESTN.log
