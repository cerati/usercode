#! /usr/bin/python

import os, sys

cuts = ['MetGt45ptll45','HWW125','HWW145','HWW150','HWW160','HWW170','HWW180','HWW190','HWW200']
#cuts = ['HWW160']
#fss = ['ll','mm','ee']
fss = ['ll']
njs = ['0','1']

for cut in cuts:
    for fs in fss:
        for nj in njs:
            os.system('./runOFVZSubtrPlots.py '+cut+' '+nj+' 0 '+fs)
