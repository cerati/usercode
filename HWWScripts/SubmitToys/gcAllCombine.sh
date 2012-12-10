#!/bin/bash

#like: hww_012j_combine
CARD=$1

MASS=$2

IT=1
#get limits
#echo combine -M Asymptotic ${CARD}  --newExpected -m ${MASS}
combine -M Asymptotic ${CARD}  --newExpected -m ${MASS} >& logLim
OBS=`tail logLim | grep "Observed Li" | awk '{printf ("%5.2f\n", $5)}'`
S2D=`tail logLim | grep "Expected  2" | awk '{printf ("%5.2f\n", $5)}'`
S1D=`tail logLim | grep "Expected 16" | awk '{printf ("%5.2f\n", $5)}'`
EXP=`tail logLim | grep "Expected 50" | awk '{printf ("%5.2f\n", $5)}'`
S1U=`tail logLim | grep "Expected 84" | awk '{printf ("%5.2f\n", $5)}'`
S2U=`tail logLim | grep "Expected 97" | awk '{printf ("%5.2f\n", $5)}'`
#get signal strength
#echo combine -M MaxLikelihoodFit ${CARD} -m ${MASS}
combine -M MaxLikelihoodFit ${CARD} -m ${MASS} >& logStr
sed -i 's/r\:/ /g' logStr
sed -i 's/\// /g' logStr
sed -i 's/\+/ /g' logStr
MU=`cat logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $3)}'`
EM=`cat logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $4)}'`
EP=`cat logStr | grep "Best fit" | awk '{printf ("%5.3f\n", $5)}'`
#get significance
#echo combine ${CARD} -M ProfileLikelihood -v 1 --significance -m ${MASS}
combine ${CARD} -M ProfileLikelihood -v 1 --significance -m ${MASS} >& logSig
SIG=`tail logSig | grep "Significance" | awk '{printf ("%5.3f\n", $2)}'`
#echo combine ${CARD} -M ProfileLikelihood -v 1 --significance -m ${MASS} --expectSignal=1 -t -1 -n Expected
combine ${CARD} -M ProfileLikelihood -v 1 --significance -m ${MASS} --expectSignal=1 -t -1 -n Expected >& logSig
EXS=`tail logSig | grep "Significance" | awk '{printf ("%5.3f\n", $2)}'`
echo 'limit: '$MASS $NEXP $OBS $EXP '['$S1D','$S1U']' '['$S2D','$S2U'] --- strength: '$MU $EM $EP' --- significance: '$SIG $EXS
rm logLim
rm higgsCombineTest.Asymptotic.mH${MASS}.root
rm logStr
rm higgsCombineTest.MaxLikelihoodFit.mH${MASS}.root
rm mlfit.root
rm logSig
rm higgsCombineTest.ProfileLikelihood.mH${MASS}.root
rm higgsCombineExpected.ProfileLikelihood.mH${MASS}.root
nrootstats=$(ls roostats-*.root 2> /dev/null | wc -l)
if [ "$nrootstats" != "0" ]; then
    rm roostats-*.root
fi

#./gcAllCombine.sh hww_012j_combine_8TeV







