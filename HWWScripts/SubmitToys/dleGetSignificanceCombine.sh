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
cp /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/getPseudoData.C .
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"cut\",0,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"cut\",0,\"sf\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"cut\",1,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"cut\",1,\"sf\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"cut\",2,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"cut\",2,\"sf\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"shape\",0,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"8TeV\",\"shape\",1,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"cut\",0,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"cut\",0,\"sf\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"cut\",1,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"cut\",1,\"sf\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"cut\",2,\"\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"shape\",0,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"shape\",0,\"sf\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"shape\",1,\"of\",\${MIN},\${MAX}\)
root -b -q getPseudoData.C\(\${MASS},\"7TeV\",\"shape\",1,\"sf\",\${MIN},\${MAX}\)
echo 'processing mass='\${MASS} >> /tmp/cerati/CMSSW_5_3_3/src/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}
for NEXP in \$(seq \${MIN} \${MAX})
  do
  echo 'toy='\${NEXP} >> /tmp/cerati/CMSSW_5_3_3/src/log_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  DIRPATH=testcards/cards_\${INJ}_N\${NEXP}/\${MASS}/
  CARDS_CUT="hwwof_0j_cut_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwof_1j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_2j_cut_8TeV.txt hwwsf_2j_cut_8TeV.txt"
  COMBC_CUT="hww_012j_cut_8TeV.txt"
  CARDS_SHAPE="hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt"
  COMBC_SHAPE="hwwof_01j_shape_8TeV.txt"
  CARDS_COMB="hwwof_0j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_2j_cut_8TeV.txt hwwsf_2j_cut_8TeV.txt"
  COMBC_COMB="hww_012j_combine_8TeV.txt"
  CARDS_ALLC="hwwof_0j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_2j_cut_8TeV.txt hwwsf_2j_cut_8TeV.txt hwwof_0j_shape_7TeV.txt hwwsf_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_1j_shape_7TeV.txt hww_2j_cut_7TeV.txt"
  COMBC_ALLC="hww_012j_allcomb_7p8TeV.txt"
  cd \${DIRPATH}
  combineCards.py -S \${CARDS_CUT}   > \${COMBC_CUT}
  combineCards.py -S \${CARDS_SHAPE} > \${COMBC_SHAPE}
  combineCards.py -S \${CARDS_COMB}  > \${COMBC_COMB}
  combineCards.py -S \${CARDS_ALLC}  > \${COMBC_ALLC}
  cp /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/gcAllCombine.sh .
  ./gcAllCombine.sh \${COMBC_CUT} \${MASS}   >> /tmp/cerati/CMSSW_5_3_3/src/log_CUT_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  ./gcAllCombine.sh \${COMBC_SHAPE} \${MASS} >> /tmp/cerati/CMSSW_5_3_3/src/log_SHAPE_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  ./gcAllCombine.sh \${COMBC_COMB} \${MASS}  >> /tmp/cerati/CMSSW_5_3_3/src/log_COMB_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  ./gcAllCombine.sh \${COMBC_ALLC} \${MASS}  >> /tmp/cerati/CMSSW_5_3_3/src/log_ALLC_\${INJ}_\${MASS}_\${MIN}_\${MAX}
  cd /tmp/cerati/CMSSW_5_3_3/src
done
mkdir -p /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/logs/\${MASS}/
cp /tmp/cerati/CMSSW_5_3_3/src/log*_\${INJ}_\${MASS}_\${MIN}_\${MAX} /afs/cern.ch/user/c/cerati/scratch0/HCP_Injection/logs/\${MASS}/
cd /tmp/cerati/
rm -rf CMSSW_5_3_3
EOF

# submit the jobs
echo 'injection: '$INJ' - processing toys:'$MIN'-'$MAX
#for MASS in 110
for MASS in {110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600}
  do
  chmod u+x ${WRAPPER}
  bsub -q ${QUEUE} -o joblog_${INJ}_${MASS}_${MIN}_${MAX} "${WRAPPER} ${INJ} ${MIN} ${MAX} ${MASS}"
done

#./dleGetSignificanceCombine.sh 125 0 50; ./dleGetSignificanceCombine.sh 125 51 100; ./dleGetSignificanceCombine.sh 125 101 150; ./dleGetSignificanceCombine.sh 125 151 200


