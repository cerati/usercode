import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_1_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_2_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_3_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_4_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_5_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_6_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_7_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_8_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_9_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_10_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_11_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_12_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_13_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_14_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_15_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_16_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_17_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_18_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_19_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_20_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_21_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_22_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_23_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_24_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_25_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_26_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_27_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_28_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_29_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_30_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_31_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_32_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_33_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_34_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_35_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_36_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_37_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_38_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_39_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_40_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_41_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_42_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_43_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_44_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_45_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_46_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_47_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_48_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_49_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_50_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_51_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_52_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_53_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_54_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_55_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_56_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_57_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_58_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_59_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_60_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_61_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_62_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_63_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_64_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_65_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_66_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_67_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_68_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_69_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_70_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_71_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_72_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_73_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_74_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_75_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_76_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_77_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_78_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_79_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_80_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_81_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_82_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_83_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_84_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_85_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_86_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_87_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_88_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_89_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_90_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_91_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_92_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_93_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_94_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_95_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_96_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_97_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_98_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_99_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_100_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_101_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_102_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_103_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_104_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_105_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_106_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_107_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_108_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_109_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_110_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_111_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_112_1.root',
        'rfio:/castor/cern.ch/user/c/cerati/Tmp/elecMuSkim_MinBias_113_1.root'
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
#removeMCMatching(process, ['All'])
addPfMET(process, 'PF')
# add genMET
#process.patMETs.addGenMET = True
#process.patMETsPF.addGenMET = True
#process.patMETs.genMETSource = cms.InputTag("genMetTrue")
#process.patMETsPF.genMETSource = cms.InputTag("genMetTrue")
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
                              outfile = cms.string('ntuple_Zanalisys7TeV_mc_MinBias.root')
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
