#!/bin/bash

# this code submits the bias studies to the lsf queue
#--------Test example, submitting 1 job to the queue
# ./gcRunNormTest.sh 125 0 0 hww cards_inj_stat
#--------Realistic example, submit 10 jobs into the queue, which 
#     analyses 100 toys 
# for j in {0,100,200,300,400,500,600,700,800,900}; do ./gcRunNormTest.sh 125 $(($j+1)) $(($j+100)) hww cards_inj_stat; done


#this is def or 125
INJ=$1
#these specify the range of toy to process
MIN=$2
MAX=$3
ANA=$4
CAR=$5

INPUTDIR=/afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/$CAR
OUTPUTDIR=$INPUTDIR
mkdir -p $OUTPUTDIR/logsNorm/


if [ ! $# -eq 5 ]; then
    echo " USAGE: ./gcRunNormTest.sh inj min max ana
    inj - choose from def and 125 (set to 125 for signal injection 
    min - range of the toys to process 
    max - range of the toys to process 
    ana - analysis, choose from hww and xww, xww is for spin 2 "
    exit 1
fi

# for the batch submission
WRAPPER="process_job.sh"
QUEUE="1nd"
WORKDIR=`pwd`
# write the wrapper
cat > ${WRAPPER} << EOF
#!/bin/bash
INJ=\$1
MIN=\$2
MAX=\$3
MASS=\$4
cd /tmp/cerati/
export SCRAM_ARCH=slc5_amd64_gcc472
scramv1 p CMSSW CMSSW_6_1_1
cd CMSSW_6_1_1/src
eval \`scramv1 ru -sh\`
addpkg HiggsAnalysis/CombinedLimit V03-01-00 
scramv1 b -j 4
perl -p -i -e "s/exp\_/exp\_final\_/g" HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py
cp $WORKDIR/diffNuisances.py HiggsAnalysis/CombinedLimit/test/
cp $WORKDIR/getPseudoData.C .
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"shape\",0,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"shape\",1,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
echo 'processing mass='\${MASS} analysis ${ANA} $NJET jet>> /tmp/cerati/CMSSW_6_1_1/src/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}_$ANA_summary.log
for NEXP in \$(seq \${MIN} \${MAX})
  do
  echo 'toy='\${NEXP} 
  DIRPATH=testcards/cards_\${INJ}_N\${NEXP}/\${MASS}/
  cd \${DIRPATH}
  combine -M MaxLikelihoodFit --saveNorm ${ANA}of_0j_shape_8TeV.txt -m \${MASS} #>& /dev/null
  python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfit.root >> /tmp/cerati/CMSSW_6_1_1/src/logNorm_\${INJ}_\${MASS}_${ANA}of_0j_shape_8TeV_\${MIN}_\${MAX}.log
  python ../../../HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a mlfit.root >> /tmp/cerati/CMSSW_6_1_1/src/logNuis_\${INJ}_\${MASS}_${ANA}of_0j_shape_8TeV_\${MIN}_\${MAX}.log
  cp mlfit.root /tmp/cerati/CMSSW_6_1_1/src/mlfit_injm${INJ}_m\${MASS}_${ANA}of_0j_id\${NEXP}.root 
  rm mlfit.root higgsCombineTest.MaxLikelihoodFit.mH125.root
  combine -M MaxLikelihoodFit --saveNorm ${ANA}of_1j_shape_8TeV.txt -m \${MASS} #>& /dev/null
  python ../../../HiggsAnalysis/CombinedLimit/test/mlfitNormsToText.py mlfit.root >> /tmp/cerati/CMSSW_6_1_1/src/logNorm_\${INJ}_\${MASS}_${ANA}of_1j_shape_8TeV_\${MIN}_\${MAX}.log
  python ../../../HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a mlfit.root >> /tmp/cerati/CMSSW_6_1_1/src/logNuis_\${INJ}_\${MASS}_${ANA}of_1j_shape_8TeV_\${MIN}_\${MAX}.log
  cp mlfit.root /tmp/cerati/CMSSW_6_1_1/src/mlfit_injm${INJ}_m\${MASS}_${ANA}of_1j_id\${NEXP}.root 
  rm mlfit.root higgsCombineTest.MaxLikelihoodFit.mH125.root
  cd /tmp/cerati/CMSSW_6_1_1/src
done
mkdir -p $OUTPUTDIR/logsNorm/\${MASS}/
# clear up old files
rm -f $OUTPUTDIR/logsNorm/\${MASS}/logNorm_\${INJ}_\${MASS}_${ANA}of_0j_shape_8TeV_\${MIN}_\${MAX}.log
rm -f $OUTPUTDIR/logsNorm/\${MASS}/logNuis_\${INJ}_\${MASS}_${ANA}of_0j_shape_8TeV_\${MIN}_\${MAX}.log
rm -f $OUTPUTDIR/logsNorm/\${MASS}/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}_$ANA_0j_summary.log
rm -f $OUTPUTDIR/logsNorm/\${MASS}/logNorm_\${INJ}_\${MASS}_${ANA}of_1j_shape_8TeV_\${MIN}_\${MAX}.log
rm -f $OUTPUTDIR/logsNorm/\${MASS}/logNuis_\${INJ}_\${MASS}_${ANA}of_1j_shape_8TeV_\${MIN}_\${MAX}.log
rm -f $OUTPUTDIR/logsNorm/\${MASS}/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}_$ANA_1j_summary.log
cp /tmp/cerati/CMSSW_6_1_1/src/log* $OUTPUTDIR/logsNorm/\${MASS}/
cp /tmp/cerati/CMSSW_6_1_1/src/mlfit_* $OUTPUTDIR/logsNorm/\${MASS}/
cd /tmp/cerati/
rm -rf CMSSW_6_1_1
EOF

# submit the jobs
echo 'injection: '$INJ' - processing toys:'$MIN'-'$MAX
for MASS in 125
  do
  chmod u+x ${WRAPPER}
  rm -f joblog_${CAR}_${INJ}_${MASS}_${MIN}_${MAX}_${ANA}.log
  bsub -q ${QUEUE} -o joblog_${CAR}_${INJ}_${MASS}_${MIN}_${MAX}_${ANA}.log "${WRAPPER} ${INJ} ${MIN} ${MAX} ${MASS}"
done


