import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/FE48CBD7-437A-E211-AE8D-003048F1CA6E.root'
    )
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START61_V11::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.tracking = cms.Sequence(
    process.siPixelRecHits*
    process.siStripMatchedRecHits*
    process.trackingGlobalReco
)

# this is the seeding presented here: https://indico.cern.ch/getFile.py/access?contribId=6&sessionId=1&resId=0&materialId=slides&confId=238024
process.load("Demo.SeedAndTrackAnalyzer.StripTriplets_cff")

process.demo = cms.EDAnalyzer('SeedAndTrackAnalyzer',
                              seeds  = cms.untracked.InputTag('stripTripletStepSeeds'),
                              tracks = cms.untracked.InputTag('stripTripletStepTracksWithQuality'),
                              TTRHBuilder = cms.string('WithTrackAngle')
)


process.p = cms.Path(process.tracking*process.StripTripletStep*process.demo)
