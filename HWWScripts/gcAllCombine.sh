#!/bin/bash

#like: hww_012j_combine
CARD=$1

IT=1
#for MASS in {110,125,160,250,300,600}
for MASS in {110,115,120,125,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
  do
  echo 'processing mass='${MASS}
  #get limits
  echo combine -M Asymptotic cards/${MASS}/${CARD}.txt  --newExpected -m ${MASS}
  combine -M Asymptotic cards/${MASS}/${CARD}.txt  --newExpected -m ${MASS} >& logLim
  OBS=`tail logLim | grep "Observed Li" | awk '{printf ("%5.2f\n", $5)}'`
  S2D=`tail logLim | grep "Expected  2" | awk '{printf ("%5.2f\n", $5)}'`
  S1D=`tail logLim | grep "Expected 16" | awk '{printf ("%5.2f\n", $5)}'`
  EXP=`tail logLim | grep "Expected 50" | awk '{printf ("%5.2f\n", $5)}'`
  S1U=`tail logLim | grep "Expected 84" | awk '{printf ("%5.2f\n", $5)}'`
  S2U=`tail logLim | grep "Expected 97" | awk '{printf ("%5.2f\n", $5)}'`
  #get signal strength
  echo combine -M MaxLikelihoodFit cards/${MASS}/${CARD}.txt -m ${MASS}
  combine -M MaxLikelihoodFit cards/${MASS}/${CARD}.txt -m ${MASS} >& logStr
  sed -i 's/r\:/ /g' logStr
  sed -i 's/\// /g' logStr
  sed -i 's/\+/ /g' logStr
  MU=`tail -15 logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $3)}'`
  EM=`tail -15 logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $4)}'`
  EP=`tail -15 logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $5)}'`
  #get significance
  echo combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS}
  combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS} >& logSig
  SIG=`tail logSig | grep "Significance" | awk '{printf ("%5.3f\n", $2)}'`
  echo combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS} --expectSignal=1 -t -1 -n Expected
  combine cards/${MASS}/${CARD}.txt -M ProfileLikelihood -v 1 --significance -m ${MASS} --expectSignal=1 -t -1 -n Expected >& logSig
  EXS=`tail logSig | grep "Significance" | awk '{printf ("%5.3f\n", $2)}'`
  if [ ${IT} == 1 ] ; then  
      rm logAll_${CARD}.txt       
  fi
  echo 'limit: '$MASS $NEXP $OBS $EXP '['$S1D','$S1U']' '['$S2D','$S2U'] --- strength: '$MU $EM $EP' --- significance: '$SIG $EXS >> logAll_${CARD}.txt       
  rm logLim
  rm higgsCombineTest.Asymptotic.mH${MASS}.root
  rm logStr
  rm higgsCombineTest.MaxLikelihoodFit.mH${MASS}.root
  rm mlfit.root
  rm logSig
  rm higgsCombineTest.ProfileLikelihood.mH${MASS}.root
  rm higgsCombineExpected.ProfileLikelihood.mH${MASS}.root
  rm roostats-*.root
  IT=0
done

#./gcAllCombine.sh hww_012j_combine







