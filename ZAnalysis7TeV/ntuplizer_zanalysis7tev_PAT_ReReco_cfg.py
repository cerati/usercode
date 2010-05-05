import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_1_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_2_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_3_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_5_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_6_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_7_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_8_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_9_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_10_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_11_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_12_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_14_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_15_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_16_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_17_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_18_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_19_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_20_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_22_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_24_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_26_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_27_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_28_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_29_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_30_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_31_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_32_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_33_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_34_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_35_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_36_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_37_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_38_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_39_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_40_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_41_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_42_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_43_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_44_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_45_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_46_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_47_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_48_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_49_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_50_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_51_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_52_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_53_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_54_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_55_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_56_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_57_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_58_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_59_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_60_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_61_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_62_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_63_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_64_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_65_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_66_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_67_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_ReReco_68_1.root'
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
                              outfile = cms.string('ntuple_Zanalisys7TeV_data_ReReco.root')
)

#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options = cms.untracked.PSet(
    fileMode = cms.untracked.string('NOMERGE')
)

process.p = cms.Path(#process.hltLevel1GTSeed*
                     process.patDefaultSequence*
                     #process.patTriggerSequence*
                     process.demo
                     )
