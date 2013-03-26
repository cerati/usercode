import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOALLTRK")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.categories     = ["CkfPattern"]
#process.MessageLogger.debugModules   = [ "*" ]
#process.MessageLogger.cout = cms.untracked.PSet(
#    threshold = cms.untracked.string("DEBUG"),
#    INFO = cms.untracked.PSet( limit = cms.untracked.int32( 0 )),
#    DEBUG = cms.untracked.PSet( limit = cms.untracked.int32( 0 )),
#    CkfPattern = cms.untracked.PSet( limit = cms.untracked.int32( -1 ))
#)


# source
sample = 0 #single mu
sample = 1 #ttbar

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
if sample==0:
    readFiles.extend( ['/store/relval/CMSSW_6_2_0_pre2-START61_V11/RelValSingleMuPt10/GEN-SIM-RECO/v1/00000/52A31EA2-7478-E211-9031-002590489C9E.root'] );
    secFiles.extend( [
        '/store/relval/CMSSW_6_2_0_pre2-START61_V11/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/A8320C0E-6378-E211-8B40-003048FEB91E.root',
        '/store/relval/CMSSW_6_2_0_pre2-START61_V11/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/2465FA67-5D78-E211-A770-003048F0E412.root',
        '/store/relval/CMSSW_6_2_0_pre2-START61_V11/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/1CEFE826-6078-E211-A886-002590494E1E.root'
        ] );
if sample==1:
    readFiles.extend( [
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/FE48CBD7-437A-E211-AE8D-003048F1CA6E.root',
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/F600059F-2E7A-E211-9399-003048F1186A.root',
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/76BF1027-B17A-E211-92CC-003048F00AF8.root',
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/6C0AE730-317A-E211-8416-003048F0E00A.root',
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/22D46130-377A-E211-81F5-003048F0258A.root',
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/0E566081-337A-E211-9E7F-003048CF97B2.root',
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/0A6F8ADE-2D7A-E211-BE3D-0025901D5E20.root' 
        ] );
    secFiles.extend( [
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/2E9776B7-0D7A-E211-AC00-003048CF94A4.root', 
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/F0205BD7-0D7A-E211-A5EA-003048CEA256.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/ECEF5FB8-117A-E211-976D-003048CF677E.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/EA8F7CCA-0D7A-E211-852D-001D09F29849.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/D89A9DB5-0D7A-E211-BFBE-003048F17F04.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/D88286B2-0D7A-E211-9FDE-00304867342C.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/C675E1B7-0D7A-E211-AF1E-0025B3203564.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/C2067D7A-0E7A-E211-847B-001E673983A9.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/BE3983A6-117A-E211-933A-003048F17C90.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/BCEF0576-157A-E211-ADD2-00304894568C.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/B8EB28A5-117A-E211-A156-001E4F3DA517.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/B46F359B-117A-E211-9C97-00304896B92C.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/960A1F9E-117A-E211-8169-003048CFAE0E.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/928191DB-117A-E211-BF99-003048F1E1E4.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/826554DF-0D7A-E211-975B-00304894555E.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/72BC9698-117A-E211-A138-003048C9C960.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/54C28FC5-0D7A-E211-9EF1-001E4F3DA513.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/4E9757B9-1E7A-E211-A011-003048CF4C16.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/2EB873FD-0D7A-E211-B8B3-0025B3244176.root',#
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/287043A8-0E7A-E211-BEE5-0025904B2C2C.root', 
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/26CB2BF9-0D7A-E211-B862-0025904B3354.root', 
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/22E41AC3-0D7A-E211-8919-003048CEFCFA.root', 
        '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/00000/0A74EFA9-117A-E211-B1A9-003048F00BA0.root'  
        ] );
    
process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

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

AllTrkEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_siPixelClusters_*_*',
    'keep *_siStripClusters_*_*',
    'keep *_*_*_RECOALLTRK'),
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize=cms.untracked.int32(5*1024*1024)
)

### here are lions
if sample==0: process.load("RecoTracker.IterativeTracking.StripTriplets_allclust_cff")
else: process.load("RecoTracker.IterativeTracking.StripTriplets_cff")

testReverse = 0
if testReverse:
    process.stripTripletStepSeedLayers.layerList = cms.vstring('TIB2+TIB1+BPix3')
    process.stripTripletStepSeedLayers.BPix = cms.PSet(
        useErrorsFromParam = cms.bool(True),
        hitErrorRZ = cms.double(0.006),
        hitErrorRPhi = cms.double(0.0027),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4MixedTriplets'),
        HitProducer = cms.string('siPixelRecHits')
        #skipClusters = cms.InputTag('stripTripletStepClusters')
        )
    process.stripTripletStepSeeds.SeedCreatorPSet.propagator = cms.string('PropagatorWithMaterialOpposite')

testRPhi = 1
if testRPhi:
    process.stripTripletStepSeedLayers = cms.ESProducer("SeedingLayersESProducer",
        ComponentName = cms.string('stripTripletStepSeedLayers'),
        layerList = cms.vstring('TIB1+TIB2+TIB3'),
        TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        #useSimpleRphiHitsCleaner = cms.bool(False),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")#Unmatched
        #skipClusters = cms.InputTag('stripTripletStepClusters')
        )
        )


import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
process.pixelLessStepTracksWithQuality = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('pixelLessStepTracks')),
    hasSelector=cms.vint32(1),
    indivShareFrac=cms.vdouble(1.0),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("pixelLessStepSelector","pixelLessStep") ),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0), pQual=cms.bool(True) )  ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)
process.stripTripletStepTracksWithQuality = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('stripTripletStepTracks')),
    hasSelector=cms.vint32(1),
    indivShareFrac=cms.vdouble(1.0),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("stripTripletStepSelector","stripTripletStep") ),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0), pQual=cms.bool(True) )  ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)
process.tobTecStepTracksWithQuality = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('tobTecStepTracks')),
    hasSelector=cms.vint32(1),
    indivShareFrac=cms.vdouble(1.0),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("tobTecStepSelector","tobTecStep") ),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0), pQual=cms.bool(True) )  ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)

process.stripPairStepTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('pixelLessStepTracks'),cms.InputTag('tobTecStepTracks')),
    hasSelector=cms.vint32(0,0),
    indivShareFrac=cms.vdouble(1.0,0.09),
    selectedTrackQuals = cms.VInputTag(cms.InputTag(''),cms.InputTag('')),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1), pQual=cms.bool(True) )  ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)

process.stripPairStepTracksWithQuality = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('pixelLessStepTracks'),cms.InputTag('tobTecStepTracks')),
    hasSelector=cms.vint32(1,1),
    indivShareFrac=cms.vdouble(1.0,0.09),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("pixelLessStepSelector","pixelLessStep"),cms.InputTag("tobTecStepSelector","tobTecStep") ),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1), pQual=cms.bool(True) )  ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)

### associators
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
### seed validator
process.load("Validation.RecoTrack.TrackerSeedValidator_cff")
process.trackerSeedValidator.associators = cms.vstring('quickTrackAssociatorByHits')
process.trackerSeedValidator.label = ['pixelLessStepSeeds','stripTripletStepSeeds','tobTecStepSeeds']#these need to be sorted in alphabetical order?!?
process.trackerSeedValidator.outputFile = 'seedValid.root'
if sample==0: process.trackerSeedValidator.label = ['stripTripletStepSeeds']

#replaces for larger TP region
#process.trackerSeedValidator.histoProducerAlgoBlock.generalTpSelector.lip   = cms.double(30.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.generalTpSelector.ptMin = cms.double(0.050)
#process.trackerSeedValidator.histoProducerAlgoBlock.generalTpSelector.tip   = cms.double(70.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.generalTpSelector.minRapidity = cms.double(-0.9)
#process.trackerSeedValidator.histoProducerAlgoBlock.generalTpSelector.maxRapidity = cms.double(0.9)
#
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.lip   = cms.double(30.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.ptMin = cms.double(0.050)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.tip   = cms.double(70.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.minRapidity = cms.double(-0.9)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.maxRapidity = cms.double(0.9)
#
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.lip   = cms.double(30.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.ptMin = cms.double(0.050)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.tip   = cms.double(70.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.minRapidity = cms.double(-0.9)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.maxRapidity = cms.double(0.9)
#
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.lip   = cms.double(30.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.ptMin = cms.double(0.050)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.tip   = cms.double(70.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.minRapidity = cms.double(-0.9)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.maxRapidity = cms.double(0.9)
#
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.lip   = cms.double(30.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.ptMin = cms.double(0.050)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.tip   = cms.double(70.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.minRapidity = cms.double(-0.9)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.maxRapidity = cms.double(0.9)
#
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.lip   = cms.double(30.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.ptMin = cms.double(0.050)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.tip   = cms.double(70.0)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.minRapidity = cms.double(-0.9)
#process.trackerSeedValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.maxRapidity = cms.double(0.9)
#process.trackerSeedValidator.outputFile = 'seedValid_largeTpReg.root'
#replaces for larger TP region

### track validator
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.multiTrackValidator.associators = cms.vstring('quickTrackAssociatorByHits')
process.multiTrackValidator.label = ['pixelLessStepTracks',   'pixelLessStepTracksHP',   'pixelLessStepTracksLoose',
                                     'stripPairStepTracks','stripPairStepTracksHP','stripPairStepTracksLoose',
                                     'stripTripletStepTracks','stripTripletStepTracksHP','stripTripletStepTracksLoose',
                                     'tobTecStepTracks',      'tobTecStepTracksHP',      'tobTecStepTracksLoose']
process.multiTrackValidator.useLogPt=cms.untracked.bool(True)
process.multiTrackValidator.minpT = cms.double(0.1)
process.multiTrackValidator.maxpT = cms.double(3000.0)
process.multiTrackValidator.nintpT = cms.int32(40)
process.multiTrackValidator.UseAssociators = cms.bool(True)
process.multiTrackValidator.runStandalone = cms.bool(True)
process.multiTrackValidator.skipHistoFit=cms.untracked.bool(False)
process.multiTrackValidator.outputFile = 'trackValid.root'
if sample==0: process.multiTrackValidator.label = ['stripTripletStepTracks']


#replaces for larger TP region
#process.multiTrackValidator.histoProducerAlgoBlock.generalTpSelector.lip   = cms.double(30.0)
#process.multiTrackValidator.histoProducerAlgoBlock.generalTpSelector.ptMin = cms.double(0.050)
#process.multiTrackValidator.histoProducerAlgoBlock.generalTpSelector.tip   = cms.double(70.0)
#process.multiTrackValidator.histoProducerAlgoBlock.generalTpSelector.minRapidity = cms.double(-0.9)
#process.multiTrackValidator.histoProducerAlgoBlock.generalTpSelector.maxRapidity = cms.double(0.9)
#
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.lip   = cms.double(30.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.ptMin = cms.double(0.050)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.tip   = cms.double(70.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.minRapidity = cms.double(-0.9)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsEta.maxRapidity = cms.double(0.9)
#
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.lip   = cms.double(30.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.ptMin = cms.double(0.050)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.tip   = cms.double(70.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.minRapidity = cms.double(-0.9)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPhi.maxRapidity = cms.double(0.9)
#
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.lip   = cms.double(30.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.ptMin = cms.double(0.050)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.tip   = cms.double(70.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.minRapidity = cms.double(-0.9)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsPt.maxRapidity = cms.double(0.9)
#
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.lip   = cms.double(30.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.ptMin = cms.double(0.050)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.tip   = cms.double(70.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.minRapidity = cms.double(-0.9)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXR.maxRapidity = cms.double(0.9)
#
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.lip   = cms.double(30.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.ptMin = cms.double(0.050)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.tip   = cms.double(70.0)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.minRapidity = cms.double(-0.9)
#process.multiTrackValidator.histoProducerAlgoBlock.TpSelectorForEfficiencyVsVTXZ.maxRapidity = cms.double(0.9)
#process.multiTrackValidator.outputFile = 'trackValid_largeTpReg.root'
#replaces for larger TP region

process.cutsRecoTracks.quality = cms.vstring('highPurity')
process.pixelLessStepTracksHP    = process.cutsRecoTracks.clone(src = cms.InputTag('pixelLessStepTracksWithQuality'))
process.stripTripletStepTracksHP = process.cutsRecoTracks.clone(src = cms.InputTag('stripTripletStepTracksWithQuality'))
process.tobTecStepTracksHP       = process.cutsRecoTracks.clone(src = cms.InputTag('tobTecStepTracksWithQuality'))
process.stripPairStepTracksHP    = process.cutsRecoTracks.clone(src = cms.InputTag('stripPairStepTracksWithQuality'))

process.pixelLessStepTracksLoose    = process.pixelLessStepTracksHP.clone(quality = cms.vstring('loose'))
process.stripTripletStepTracksLoose = process.stripTripletStepTracksHP.clone(quality = cms.vstring('loose'))
process.tobTecStepTracksLoose       = process.tobTecStepTracksHP.clone(quality = cms.vstring('loose'))
process.stripPairStepTracksLoose    = process.stripPairStepTracksHP.clone(quality = cms.vstring('loose'))

process.out = cms.OutputModule("PoolOutputModule",
    AllTrkEventContent,
    fileName = cms.untracked.string('reco_trk_stripTriplets_ttbar.root')
)
if sample==0: process.out.fileName = cms.untracked.string('reco_trk_stripTriplets_singleMu.root')

process.beamSpot = cms.Sequence(
    process.offlineBeamSpot
)
process.clustToHits = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits
)

process.tracking = cms.Sequence(
    process.trackingGlobalReco
)

process.validation = cms.Sequence(
    process.trackerSeedValidator*
    process.pixelLessStepTracksWithQuality*process.stripTripletStepTracksWithQuality*process.tobTecStepTracksWithQuality*
    process.stripPairStepTracks*process.stripPairStepTracksWithQuality*
    process.pixelLessStepTracksLoose*process.stripTripletStepTracksLoose*process.tobTecStepTracksLoose*process.stripPairStepTracksLoose*
    process.pixelLessStepTracksHP*process.stripTripletStepTracksHP*process.tobTecStepTracksHP*process.stripPairStepTracksHP*
    process.multiTrackValidator
)

if sample==0:
    process.validation.remove(process.pixelLessStepTracksWithQuality)
    #process.validation.remove(process.stripTripletStepTracksWithQuality)
    process.validation.remove(process.tobTecStepTracksWithQuality)
    process.validation.remove(process.stripPairStepTracks)
    process.validation.remove(process.stripPairStepTracksWithQuality)
    process.validation.remove(process.pixelLessStepTracksLoose)
    process.validation.remove(process.stripTripletStepTracksLoose)
    process.validation.remove(process.tobTecStepTracksLoose)
    process.validation.remove(process.stripPairStepTracksLoose)
    process.validation.remove(process.pixelLessStepTracksHP)
    process.validation.remove(process.stripTripletStepTracksHP)
    process.validation.remove(process.tobTecStepTracksHP)
    process.validation.remove(process.stripPairStepTracksHP)

# paths
process.p = cms.Path(
      process.clustToHits
    * process.beamSpot
    * process.tracking
    * process.StripTripletStep
    * process.validation
)
process.o = cms.EndPath(
      process.out
)

process.schedule = cms.Schedule(
      process.p,
      process.o
)


#process.Timing = cms.Service("Timing")
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

