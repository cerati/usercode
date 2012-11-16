#!/bin/bash

root -b -q -l dyBg.C++\(12.1\) >& log.dyBg
#mv DYBkgScaleFactors.h Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h
root -b -q -l topBg.C++\(12.1\) >& log.topBg
#mv TopBkgScaleFactors.h Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h
root -b -q -l wwBkg.C++\(12.1\) >& log.wwBg
#mv WWBkgScaleFactors.h Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h
