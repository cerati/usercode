#! /usr/bin/env python

import os

masses = [115,120,130,140,150,160,170,180,190,200,250,300]
jetbins = [0,1]
fstates = ['offs','sffs']
#masses = [120]
#jetbins = [0]
#fstates = ['offs']

for mass in masses:
    for njets in jetbins:
        for fs in fstates:
            os.system('rm BDTG.class.C')
            os.system('ln -s /smurf/ceballos/tmva/weights/ntuples_'+str(mass)+'train_'+str(njets)+'jets_BDTG.class.C BDTG.class.C')
            os.system('root -b -l -q shapeMaker.C++\(4.7,'+str(njets)+','+str(mass)+',\\"'+fs+'\\"\)')
            os.system('mkdir -p cards/'+str(mass))
            os.system('mv hww*j.input.root cards/'+str(mass))
            #print('root -b -l -q shapeMaker.C++\(4.7,'+str(njets)+','+str(mass)+',\\"'+fs+'\\"\)')
