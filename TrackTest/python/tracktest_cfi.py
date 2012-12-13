import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('TrackTest'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
