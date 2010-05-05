#! /bin/sh
limit vmemoryuse unlimited
cd /afs/cern.ch/user/c/cerati/scratch0/zanalisys/CMSSW_3_5_6/src/Tests/ZAnalysis7TeV
eval `scramv1 runtime -sh`
cd /tmp/
cmsRun /afs/cern.ch/user/c/cerati/scratch0/zanalisys/CMSSW_3_5_6/src/Tests/ZAnalysis7TeV/ntuplizer_zanalysis7tev_PAT_MinBias_cfg.py
#rfcp .root /castor/cern.ch/user/c/cerati/ 
#rm .root
cp ntuple_Zanalisys7TeV_mc_MinBias.root /afs/cern.ch/user/c/cerati/scratch0/zanalisys/CMSSW_3_5_6/src/Tests/ZAnalysis7TeV/
rm ntuple_Zanalisys7TeV_mc_MinBias.root
