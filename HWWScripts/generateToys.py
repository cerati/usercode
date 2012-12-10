#! /usr/bin/env python

import os, sys

masses = [110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600]
#masses = [120]
#masses = [110,125,160,250,300,600]
jetbins = [0,1,2]
fstates = ['of','sf','']
cmss = ['7TeV','8TeV']
types = ['cut','shape']

#masses = [110]
#jetbins = [0]
#fstates = ['sffs']

if len(sys.argv)!=2: sys.exit('plese specify the card directory')
dir = sys.argv[1]

seed = 1
os.system('rm '+dir+'/*/hww*PseudoData*.root')
os.system('rm '+dir+'/*/toy*.log')
for mass in masses:
    os.chdir(os.getcwd()+'/'+dir+'/'+str(mass)+'/')
    for cms in cmss:
        for njets in jetbins:
            for fs in fstates:
                for type in types:
                    if cms=='8TeV' and type=='shape' and fs!='of': continue
                    if type=='shape' and njets==2: continue
                    if fs=='' and (cms!='7TeV' or type!='cut' or njets!=2): continue
                    if cms=='7TeV' and njets==2 and fs!='': continue
                    card = 'hww'+fs+'_'+str(njets)+'j_'+type+'_'+cms
                    os.system('../../LandS/test/lands.exe -d '+card+'.txt -M Hybrid  -m '+str(mass)+' --minuitSTRATEGY 0  --bWriteToys 1 -n \"'+card+'\" --nToysForCLsb 1000 --nToysForCLb 1 --singlePoint 1 --seed '+str(seed)+'  -rMin 0 -rMax 5 >& toy_'+card+'.log') # --freq
                    os.system('mv '+card+'_PseudoData_sb_seed'+str(seed)+'.root '+card+'_PseudoData_sb.root')
                    os.system('mv '+card+'_PseudoData_b_seed'+str(seed)+'.root '+card+'_PseudoData_b.root')
                    seed=seed+1
    os.chdir(os.getcwd()+'/../../')

#./generateToys.py cards_inj_stat
