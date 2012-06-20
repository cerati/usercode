#! /usr/bin/python

import os, sys

masses = ['115','120','125','130','140','150','160','300']
njs = ['0','1','2']
fss = ['','of','sf']

#masses = ['160']
#njs = ['0']
fss = ['of','sf']

for mass in masses:
    for nj in njs:
        for fs in fss:
            os.system('root -b -q plotBaby.C\(3.553,'+nj+','+mass+',\\"'+fs+'\\"\)')
