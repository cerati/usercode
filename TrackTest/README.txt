cmsrel CMSSW_6_2_0_pre2
cd CMSSW_6_2_0_pre2/src/
ls -latr
mkdir Test
cd Test/
cvs co -d TrackTest UserCode/GCerati/TrackTest
cd TrackTest/
scramv1 b -j 4
cmsenv
cmsRun tracktest_cfg.py 
