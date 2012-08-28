#! /usr/bin/python

import os, sys

print 'cut:'+sys.argv[1]
print 'njets:'+sys.argv[2]
print 'norm:'+sys.argv[3]
print 'fs:'+sys.argv[4]

os.system('root -b -q plotOFVZsub.C\(\\"dilep.mass\(\)\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
os.system('root -b -q plotOFVZsub.C\(\\"dilep.pt\(\)\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
os.system('root -b -q plotOFVZsub.C\(\\"pmet\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
os.system('root -b -q plotOFVZsub.C\(\\"pTrackMet\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
os.system('root -b -q plotOFVZsub.C\(\\"mt\\",40,0,300,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
os.system('root -b -q plotOFVZsub.C\(\\"lep1.pt\(\)\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
os.system('root -b -q plotOFVZsub.C\(\\"lep2.pt\(\)\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')

#os.system('root -b -q plotOFVZsub.C\(\\"dPhiJet1MET\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"recoil\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"nvtx\\",40,0,40,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhiDiLepJet1\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"jet1.pt\(\)\\",40,0,400,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhiDiLepMET\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"met/sqrt\(sumet\)\\",40,0,20,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')

#os.system('root -b -q plotOFVZsub.C\(\\"sumet\\",50,0,2000,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"met\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"trackMet\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"metSig\\",40,0,200,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dymva\\",20,-1,1,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhiLep1Jet1\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhiLep2Jet1\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhiLep1MET\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhiLep2MET\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dRLep1Jet1\\",35,0,7,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dRLep2Jet1\\",35,0,7,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
#os.system('root -b -q plotOFVZsub.C\(\\"dPhi\\",32,0,3.2,'+sys.argv[3]+',\\"'+sys.argv[1]+'\\",'+sys.argv[2]+',\\"'+sys.argv[4]+'\\"\)')
