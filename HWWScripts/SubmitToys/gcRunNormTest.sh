#!/bin/bash

#this is def or 125
INJ=$1

#these specify the range of toy to process
MIN=$2
MAX=$3

# for the batch submission
WRAPPER="process_job.sh"
QUEUE="1nd"

# write the wrapper
cat > ${WRAPPER} << EOF
#!/bin/bash
INJ=\$1
MIN=\$2
MAX=\$3
MASS=\$4
cd /tmp/cerati/
scramv1 p CMSSW CMSSW_5_3_3
cd CMSSW_5_3_3/src
eval \`scramv1 ru -sh\`
addpkg HiggsAnalysis/CombinedLimit V02-02-03 
scramv1 b
perl -p -i -e "s/exp\_/exp\_final\_/g" HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py
cp /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/getPseudoData.C .
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"shape\",0,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"shape\",1,\"of\",\${MIN},\${MAX}\)
echo 'processing mass='\${MASS} >> /tmp/cerati/CMSSW_5_3_3/src/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}
for NEXP in \$(seq \${MIN} \${MAX})
  do
  echo 'toy='\${NEXP} 
  DIRPATH=testcards/cards_\${INJ}_N\${NEXP}/\${MASS}/
  cd \${DIRPATH}
  combine -M MaxLikelihoodFit --saveNormalizations hwwof_0j_shape_8TeV.txt -m \${MASS} >& /dev/null
  python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfit.root >> /tmp/cerati/CMSSW_5_3_3/src/logNorm_\${INJ}_\${MASS}_hwwof_0j_shape_8TeV_\${MIN}_\${MAX}.log
  python ../../../HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a mlfit.root >> /tmp/cerati/CMSSW_5_3_3/src/logNuis_\${INJ}_\${MASS}_hwwof_0j_shape_8TeV_\${MIN}_\${MAX}.log
  rm mlfit.root higgsCombineTest.MaxLikelihoodFit.mH125.root
  combine -M MaxLikelihoodFit --saveNormalizations hwwof_1j_shape_8TeV.txt -m \${MASS} >& /dev/null
  python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfit.root >> /tmp/cerati/CMSSW_5_3_3/src/logNorm_\${INJ}_\${MASS}_hwwof_1j_shape_8TeV_\${MIN}_\${MAX}.log
  python ../../../HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a mlfit.root >> /tmp/cerati/CMSSW_5_3_3/src/logNuis_\${INJ}_\${MASS}_hwwof_1j_shape_8TeV_\${MIN}_\${MAX}.log
  rm mlfit.root higgsCombineTest.MaxLikelihoodFit.mH125.root
  cd /tmp/cerati/CMSSW_5_3_3/src
done
mkdir -p /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/logsNorm/\${MASS}/
cp /tmp/cerati/CMSSW_5_3_3/src/log* /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/logsNorm/\${MASS}/
cd /tmp/cerati/
rm -rf CMSSW_5_3_3
EOF

# submit the jobs
echo 'injection: '$INJ' - processing toys:'$MIN'-'$MAX
for MASS in 125
#for MASS in {110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
  do
  chmod u+x ${WRAPPER}
  bsub -q ${QUEUE} -o joblog_${INJ}_${MASS}_${MIN}_${MAX} "${WRAPPER} ${INJ} ${MIN} ${MAX} ${MASS}"
done

#./gcRunNormTest.sh 125 0 0
# for j in {0,100,200,300,400,500,600,700,800,900}; do ./gcRunNormTest.sh 125 $(($j+1)) $(($j+100)); done


