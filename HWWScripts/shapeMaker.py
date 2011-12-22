#! /usr/bin/env python

import os

masses = [120,130,140]
jetbins = [0,1]
fstates = ['offs','sffs']

for mass in masses:
    for njets in jetbins:
        for fs in fstates:
            os.system('rm BDTG.class.C')
            os.system('ln -s /smurf/ceballos/tmva/weights/ntuples_'+str(mass)+'train_'+str(njets)+'jets_BDTG.class.C BDTG.class.C')
            os.system('root -b -l -q shapeMaker.C++\(4.7,'+str(njets)+','+str(mass)+',\\"'+fs+'\\"\)')
            os.system('mkdir -p cards/'+str(mass))
            os.system('mv hww*j.input.root cards/'+str(mass))
            #print('root -b -l -q shapeMaker.C++\(4.7,'+str(njets)+','+str(mass)+',\\"'+fs+'\\"\)')
