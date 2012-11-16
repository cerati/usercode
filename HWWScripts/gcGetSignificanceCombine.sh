#!/bin/bash

#this is cut, shape, combine, ...
CARD=$1

IT=1
#for MASS in {110,125,160,250,600}
for MASS in {110,115,120,125,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
  do
  echo 'processing mass='${MASS}
  echo combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS}
  combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS} >& logSig
  SIG=`tail logSig | grep "Significance" | awk '{printf ("%5.3f\n", $2)}'`
  echo combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS} --expectSignal=1 -t -1 -n Expected
  combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS} --expectSignal=1 -t -1 -n Expected >& logSig
  EXS=`tail logSig | grep "Significance" | awk '{printf ("%5.3f\n", $2)}'`
  if [ ${IT} == 1 ] ; then  
      rm logSignif_${CARD}.txt       
  fi
  echo $MASS $SIG $EXS >> logSignif_${CARD}.txt       
  rm logSig
  rm higgsCombineTest.ProfileLikelihood.mH${MASS}.root
  rm higgsCombineExpected.ProfileLikelihood.mH${MASS}.root
  IT=0
done

#./gcGetSignificanceCombine.sh hww_012j_combine



