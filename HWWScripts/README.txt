*** Instructions to check out the code and produce limits ***

* installation *
cvs co -d HWWScripts UserCode/GCerati/HWWScripts
cd HWWScripts
cvs co -d CMS2/NtupleMacros/Tools UserCode/JRibnik/CMS2/NtupleMacros/Tools
cvs co -d Smurf/Core UserCode/Smurf/Core
cvs co -d Smurf/Analysis UserCode/Smurf/Analysis
cvs co -d Smurf/LimitCalc UserCode/Smurf/LimitCalc
cvs co -d LandS UserCode/mschen/LandS
ln -s /smurf/ceballos/tmva/weights/ntuples_160train_0jets_BDTG.class.C BDTG.class.C
chmod u+x shapeMaker.py
cd LandS
export SCRAM_ARCH=slc5_amd64_gcc462 
cmsrel CMSSW_5_2_5
cd CMSSW_5_2_5/src 
eval `scramv1 runtime -sh`
addpkg HiggsAnalysis/CombinedLimit V02-02-03 
scramv1 b -j 4
cd ../..
make
cd ../Smurf/LimitCalc/
ln -s ../../LandS/lands.so
cd ../..

* background estimation *
root -b -q -l dyBg.C++\(4.7\) >& log.dyBg
mv DYBkgScaleFactors.h Smurf/Analysis/HWWlvlv/
root -b -q -l topBg.C++\(4.7\) >& log.topBg
mv TopBkgScaleFactors.h Smurf/Analysis/HWWlvlv/
root -b -q -l wwBkg.C++\(4.7\) >& log.wwBg
mv WWBkgScaleFactors.h Smurf/Analysis/HWWlvlv/

* prepare cards *
root -b -q -l cardMaker.C++\(4.7,\"cut\"\)
root -b -q -l cardMaker.C++\(4.7,\"shape\"\)
./shapeMaker.py 4.7

* run limits *
cd Smurf/LimitCalc/
./fixPath.pl ../../cards/
./GetExpectLimits.pl ../../cards/limits_nj_shape.txt CLs-asymptotic
./postprocess.pl ../../cards/limits_nj_shape.txt "0/1/2-jets"  CLs-asymptotic

* commit only in my area *
cvs tag -l TAGNAME .
cvs tag -l TAGNAME cards/limits_*.txt
