import FWCore.ParameterSet.Config as cms

process = cms.Process("elecMuSkim")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_35X_V7A::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/A4DD1FAE-B178-DE11-B608-001D09F24EAC.root',
        '/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/9408B54D-CB78-DE11-9AEB-001D09F2503C.root'
    )
)

#--------------------------------------------------------------------------------
# select electrons and muons
#--------------------------------------------------------------------------------

process.selectedElectrons = cms.EDFilter("GsfElectronSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string("pt > 0.1 & abs(eta) < 2.5"),
    filter = cms.bool(True)
)

process.selectedMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string("pt > 0.1 & abs(eta) < 2.5"),
    filter = cms.bool(True)
)

#process.elecMuPairs = cms.EDProducer("DiCandidatePairProducer",
#    useLeadingTausOnly = cms.bool(False),
#    srcLeg1 = cms.InputTag('selectedElectrons'),
#    srcLeg2 = cms.InputTag('selectedMuons'),
#    dRmin12 = cms.double(0.),
#    srcMET = cms.InputTag(''),
#    recoMode = cms.string(""),
#    scaleFuncImprovedCollinearApprox = cms.string('1'),                                 
#    verbosity = cms.untracked.int32(0)
#)


#--------------------------------------------------------------------------------
# save events passing the electron + muon selection
#--------------------------------------------------------------------------------

elecMuEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('elecMuSkimPath')
    )
)

process.elecMuSkimOutputModule = cms.OutputModule("PoolOutputModule",
    elecMuEventSelection,
    fileName = cms.untracked.string('elecMuSkim_ReReco.root')                                                  
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
process.physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
process.bit40 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('40 AND NOT (36 OR 37 OR 38 OR 39)'))
process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"),
   filter = cms.bool(True)
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

#--------------------------------------------------------------------------------
# keep event in case it passed the muon + electron selection
#--------------------------------------------------------------------------------

process.elecMuSkimPath = cms.Path(process.physDecl*process.bit40*process.noScraping*process.oneGoodVertexFilter*
    (process.selectedElectrons + process.selectedMuons)
)

process.o = cms.EndPath(process.elecMuSkimOutputModule)
