#!/bin/bash

#this is cut, shape, combine, ...
CARD=$1

IT=1
#for MASS in {110,125,160,250,300,600}
for MASS in {110,115,120,125,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
  do
  echo 'processing mass='${MASS}
  echo combine -M MaxLikelihoodFit cards/${MASS}/${CARD}.txt -m ${MASS}
  combine -M MaxLikelihoodFit cards/${MASS}/${CARD}.txt -m ${MASS} >& logStr
  sed -i 's/r\:/ /g' logStr
  sed -i 's/\// /g' logStr
  MU=`cat logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $3)}'`
  EM=`cat logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $4)}'`
  EP=`cat logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $5)}'`
  if [ ${IT} == 1 ] ; then  
      rm logStrength_${CARD}.txt       
  fi
  echo $MASS $MU $EM $EP >> logStrength_${CARD}.txt       
  rm logStr
  rm higgsCombineTest.MaxLikelihoodFit.mH${MASS}.root
  rm mlfit.root
  IT=0
done

#./gcGetStrengthCombine.sh hww_012j_combine



