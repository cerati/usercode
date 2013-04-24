import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('SeedAndTrackAnalyzer'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
