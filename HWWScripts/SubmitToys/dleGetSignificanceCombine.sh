#!/bin/bash

#this is def or 125
INJ=$1

#these specify the range of toy to process
MIN=$2
MAX=$3

ANA=$4
CAR=$5

INPUTDIR=/afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/$CAR
OUTPUTDIR=$INPUTDIR
mkdir -p $OUTPUTDIR/logs/

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
scramv1 p CMSSW CMSSW_6_1_1
cd CMSSW_6_1_1/src
eval \`scramv1 ru -sh\`
addpkg HiggsAnalysis/CombinedLimit V03-01-00 
scramv1 b -j 4
cp /afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/getPseudoData.C .
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"cut\",0,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"cut\",0,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"cut\",1,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"cut\",1,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
#root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"cut\",2,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
#root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"cut\",2,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"shape\",0,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"8TeV\",\"shape\",1,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"cut\",0,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"cut\",0,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"cut\",1,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"cut\",1,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
#root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"cut\",2,\"\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"shape\",0,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
#root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"shape\",0,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"shape\",1,\"of\",\"${ANA}\",\${MIN},\${MAX}\)
#root -b -q getPseudoData.C\(\"$INPUTDIR\",\${MASS},\"7TeV\",\"shape\",1,\"sf\",\"${ANA}\",\${MIN},\${MAX}\)
echo 'processing mass='\${MASS} >> /tmp/cerati/CMSSW_6_1_1/src/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}
for NEXP in \$(seq \${MIN} \${MAX})
  do
  echo 'toy='\${NEXP} >> /tmp/cerati/CMSSW_6_1_1/src/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  DIRPATH=testcards/cards_\${INJ}_N\${NEXP}/\${MASS}/
  CARDS_CUT="hwwof_0j_cut_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwof_1j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_0j_cut_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwof_1j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt"
  COMBC_CUT="hww_01j_cut_8TeV.txt"
  #CARDS_SHAPE="hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt"
  #COMBC_SHAPE="hwwof_01j_shape_8TeV.txt"
  #CARDS_COMB="hwwof_0j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_2j_cut_8TeV.txt hwwsf_2j_cut_8TeV.txt"
  #COMBC_COMB="hww_012j_combine_8TeV.txt"
  #CARDS_ALLC="hwwof_0j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_0j_shape_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_1j_cut_7TeV.txt"
  #COMBC_ALLC="hww_012j_allcomb_7p8TeV.txt"
  cd \${DIRPATH}
  combineCards.py -S \${CARDS_CUT}   > \${COMBC_CUT}
  #combineCards.py -S \${CARDS_SHAPE} > \${COMBC_SHAPE}
  #combineCards.py -S \${CARDS_COMB}  > \${COMBC_COMB}
  #combineCards.py -S \${CARDS_ALLC}  > \${COMBC_ALLC}
  cp /afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/gcAllCombine.sh .
  ./gcAllCombine.sh \${COMBC_CUT}   \${MASS}   >> /tmp/cerati/CMSSW_6_1_1/src/log_CUT_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  #./gcAllCombine.sh \${COMBC_SHAPE} \${MASS} >> /tmp/cerati/CMSSW_6_1_1/src/log_SHAPE_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  #./gcAllCombine.sh \${COMBC_COMB}  \${MASS}  >> /tmp/cerati/CMSSW_6_1_1/src/log_COMB_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  #./gcAllCombine.sh \${COMBC_ALLC}  \${MASS}  >> /tmp/cerati/CMSSW_6_1_1/src/log_ALLC_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  cd /tmp/cerati/CMSSW_6_1_1/src
done
mkdir -p /afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/logs/\${MASS}/
cp /tmp/cerati/CMSSW_6_1_1/src/log*_\${INJ}_\${MASS}_\${MIN}_\${MAX} /afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/logs/\${MASS}/
cd /tmp/cerati/
rm -rf CMSSW_6_1_1
EOF

# submit the jobs
echo 'injection: '$INJ' - processing toys:'$MIN'-'$MAX
#for MASS in 160
for MASS in {110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
  do
  chmod u+x ${WRAPPER}
  bsub -q ${QUEUE} -o joblog_${INJ}_${MASS}_${MIN}_${MAX} "${WRAPPER} ${INJ} ${MIN} ${MAX} ${MASS}"
done

#./dleGetSignificanceCombine.sh 125 0 0 hww cards_inj_statsyst78_new
#for j in {0,100,200,300,400,500,600,700,800,900}; do ./dleGetSignificanceCombine.sh 125 $(($j+1)) $(($j+100)) hww cards_inj_statsyst78_new; done
#for j in {0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950}; do ./dleGetSignificanceCombine.sh 125 $(($j+1001)) $(($j+1050)) hww cards_inj_statsyst78_new; done


