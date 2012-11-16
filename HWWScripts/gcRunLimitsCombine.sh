#!/bin/bash

#like: hww_012j_combine
CARD=$1

IT=1
for MASS in {110,125,160,250,600}
  do
  echo 'processing mass='${MASS}
  echo combine -M Asymptotic cards/${MASS}/${CARD}.txt  --newExpected -m ${MASS}
  combine -M Asymptotic cards/${MASS}/${CARD}.txt  --newExpected -m ${MASS} >& logLim
  OBS=`tail logLim | grep "Observed Li" | awk '{printf ("%5.2f\n", $5)}'`
  S2D=`tail logLim | grep "Expected  2" | awk '{printf ("%5.2f\n", $5)}'`
  S1D=`tail logLim | grep "Expected 16" | awk '{printf ("%5.2f\n", $5)}'`
  EXP=`tail logLim | grep "Expected 50" | awk '{printf ("%5.2f\n", $5)}'`
  S1U=`tail logLim | grep "Expected 84" | awk '{printf ("%5.2f\n", $5)}'`
  S2U=`tail logLim | grep "Expected 97" | awk '{printf ("%5.2f\n", $5)}'`
  if [ ${IT} == 1 ] ; then  
      rm logLimits_${CARD}.txt       
  fi
  echo $MASS $NEXP $OBS $EXP '['$S1D','$S1U']' '['$S2D','$S2U']' >> logLimits_${CARD}.txt       
  rm logLim
  rm higgsCombineTest.Asymptotic.mH${MASS}.root
  IT=0
done

#./gcRunLimitsCombine.sh hww_012j_combine







