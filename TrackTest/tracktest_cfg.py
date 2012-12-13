import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START60_V4::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_6_0_0-START60_V4/RelValSingleMuPt10/GEN-SIM-RECO/v1/0000/F2040CF7-6AF3-E111-BD00-003048678DA2.root',
        '/store/relval/CMSSW_6_0_0-START60_V4/RelValSingleMuPt10/GEN-SIM-RECO/v1/0000/D60179F4-6AF3-E111-8B51-00261894395B.root',
        '/store/relval/CMSSW_6_0_0-START60_V4/RelValSingleMuPt10/GEN-SIM-RECO/v1/0000/48B4E0F9-6AF3-E111-A775-0018F3D09624.root'
    )
)

process.demo = cms.EDAnalyzer('TrackTest',
        tracks = cms.untracked.InputTag('TrackRefitter'),
        TTRHBuilder = cms.string('WithAngleAndTemplate')
)


process.p = cms.Path(process.TrackRefitter*process.demo)
