Instructions to check out the code and produce limits:
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
cmsrel CMSSW_4_2_3
cd CMSSW_4_2_3
eval `scramv1 runtime -sh`
cd ..
make
cd ../Smurf/LimitCalc/
ln -s ../../LandS/lands.so
cd ../..
root -b -q -l cardMaker.C++\(4.7,\"cut\"\)
root -b -q -l cardMaker.C++\(4.7,\"shape\"\)
./shapeMaker.py
cd Smurf/LimitCalc/
./fixPath.pl ../../cards/
./GetExpectLimits.pl ../../cards/limits_nj_shape.txt CLs-asymptotic
./postprocess.pl ../../cards/limits_nj_shape.txt "0/1/2-jets"  CLs-asymptotic
