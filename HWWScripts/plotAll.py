#! /usr/bin/python

import os, sys

masses = ['115','120','125','130','140','150','160','300']
njs = ['0','1','2']
fss = ['','of','sf']

lumi = '19.467'

masses = ['0']
njs = ['0','1']
fss = ['of','sf']
fss = ['of','sf','mm','ee','em','me']

for mass in masses:
    for nj in njs:
        for fs in fss:
            os.system('root -b -q plotBaby.C+\('+lumi+','+nj+','+mass+',\\"'+fs+'\\"\)')
