#!/bin/bash

root -b -q -l dyBg.C++\(4.7\) >& log.dyBg
mv DYBkgScaleFactors.h Smurf/Analysis/HWWlvlv/
root -b -q -l topBg.C++\(4.7\) >& log.topBg
mv TopBkgScaleFactors.h Smurf/Analysis/HWWlvlv/
root -b -q -l wwBkg.C++\(4.7\) >& log.wwBg
mv WWBkgScaleFactors.h Smurf/Analysis/HWWlvlv/
