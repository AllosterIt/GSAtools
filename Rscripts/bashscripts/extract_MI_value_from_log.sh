#!/bin/bash
DATADIR=.
OUTDIR=.
OUTPREFIX=example


BC=`which bc`

if [ -n "$BC" ]
then
echo "     MI      jH    nMI filename" > $OUTPREFIX.MI.table
else
echo "     MI      jH filename" > $OUTPREFIX.MI.table
fi

for FILE in `ls -1 *lf_MI.*.log`
do
MI=`grep -h 'Mutual Information:' $FILE | cut -d ':' -f 2`
jH=`grep -h 'joint Entropy:' $FILE | cut -d ':' -f 2`
        
if [ -n "$BC" ]
then
nMI=`echo "scale=5;$MI/$jH" | bc -l`
echo $MI $jH $nMI $FILE>>  $OUTPREFIX.MI.table
else
echo $MI $jH $FILE>>  $OUTPREFIX.MI.table
fi

done
