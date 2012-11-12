#! /usr/bin/env python

import os, sys

masses = [110,115,120,125,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600]
jetbins = [0,1]
fstates = ['offs']
shape = 'mtmll2D'

#masses = [125,200,350]
#jetbins = [0]

if len(sys.argv)!=2: sys.exit('plese specify the lumi as only argument')

for mass in masses:
    for njets in jetbins:
        for fs in fstates:
            os.system('rm BDTG.class.C')
            #os.system('ln -s /home/users/cerati/HWWMVAweights/ntuples_'+str(mass)+'train_'+str(njets)+'jets_BDTG.class.C BDTG.class.C')
            os.system('ln -s /smurf/ceballos/tmva/weights/ntuples2012_'+str(mass)+'train_'+str(njets)+'jets_BDTG.class.C BDTG.class.C')
            os.system('root -b -l -q shapeMaker.C++\('+sys.argv[1]+','+str(njets)+','+str(mass)+',\\"'+fs+'\\",\\"'+shape+'\\"\)')
            os.system('mkdir -p cards/'+str(mass))
            os.system('mv hww*j.input.root cards/'+str(mass))
            #print('root -b -l -q shapeMaker.C++\('+sys.argv[1]+','+str(njets)+','+str(mass)+',\\"'+fs+'\\"\)')
