import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_13_1.root'
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_10_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_11_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_12_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_13_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_14_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_15_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_16_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_17_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_18_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_19_1.root',
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_1_1.root',#passed = 0
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_20_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_21_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_22_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_23_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_24_1.root',
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_25_1.root',#passed = 0
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_26_1.root',#passed = 0
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_2_1.root',#passed = 0
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_3_1.root',#passed = 0
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_4_1.root',#passed = 0
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_5_1.root',#passed = 0
        #'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_6_2.root',#passed = 0
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_7_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_8_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_GoodColl_9_1.root'
)
)

process.out = cms.OutputModule("PoolOutputModule",                                 
  fileName = cms.untracked.string('dummy.root'),
  outputCommands = cms.untracked.vstring('drop *')
)

# NEW PAT CONFIGURATION
# see: 
# http://cmscvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatExamples/test/patLayer1_fromRECO_7TeV_firstdata_cfg.py?revision=1.1.2.2&view=markup&sortby=date&pathrev=V00-02-16
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookPATExampleData)

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
# add pf met
from PhysicsTools.PatAlgos.tools.metTools import *
removeMCMatching(process, ['All'])
addPfMET(process, 'PF')
# add genMET
#process.metProducer.addGenMET = True
# get the 7 TeV jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Summer09_7TeV_ReReco332")
# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = False,
                 doBTagging   = False,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )

process.demo = cms.EDAnalyzer('NtuplizerZAnalysis7TeVPAT',
                              outfile = cms.string('ntuple_Zanalisys7TeV_data7TeV.root')
)

#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.p = cms.Path(#process.hltLevel1GTSeed*
                     process.patDefaultSequence*
                     #process.patTriggerSequence*
                     process.demo
                     )
