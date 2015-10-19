import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
				     SkipEvent = cms.untracked.vstring('ProductNotFound'))
corrJetsOnTheFly = True
runOnMC = True#False
#****************************************************************************************************#

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   process.GlobalTag.globaltag =  '74X_mcRun2_asymptotic_realisticBS_v1'#'74X_mcRun2_asymptotic_v2'#'MCRUN2_74_V9::All'
elif not(runOnMC):
   process.GlobalTag.globaltag = '74X_dataRun2_v2' #74X_dataRun2_Prompt_v2

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters
hltFiltersProcessName = 'RECO'
if runOnMC:
   hltFiltersProcessName = 'PAT' #'RECO'

process.load("VAJets.PKUCommon.goodMuons_cff")
process.load("VAJets.PKUCommon.goodElectrons_cff")
process.load("VAJets.PKUCommon.goodJets_cff")
process.load("VAJets.PKUCommon.goodPhotons_cff")
process.load("VAJets.PKUCommon.leptonicW_cff")

# If Update
process.goodMuons.src = "slimmedMuons"
process.goodElectrons.src = "slimmedElectrons"
process.goodAK4Jets.src = "slimmedJets"
process.goodPhotons.src = "slimmedPhotons"
process.Wtoenu.MET  = "slimmedMETs"
process.Wtomunu.MET = "slimmedMETs"

process.goodOfflinePrimaryVertex = cms.EDFilter("VertexSelector",
                                       src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                       cut = cms.string("chi2!=0 && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"),
                                       filter = cms.bool(False)
                                       )

WBOSONCUT = "pt > 0.0"

process.leptonicVSelector = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("leptonicV"),
                                       cut = cms.string( WBOSONCUT ), 
                                       filter = cms.bool(False)
                                       )

process.leptonicVFilter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("leptonicV"),
                                       minNumber = cms.uint32(0),
                                       filter = cms.bool(False)
                                       )


process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.eleSequence +
                                      process.leptonicVSequence +
                                      process.leptonicVSelector +
                                      process.leptonicVFilter )

process.jetSequence = cms.Sequence(process.NJetsSequence)


#begin------------JEC on the fly--------
if runOnMC:
   jecLevelsAK4chs = [
                                   'Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt',
                                   'Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt',
                                   'Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt'
     ]
else:
   jecLevelsAK4chs = [
                                   'Summer15_25nsV5_DATA_L1FastJet_AK4PFchs.txt',
                                   'Summer15_25nsV5_DATA_L2Relative_AK4PFchs.txt',
                                   'Summer15_25nsV5_DATA_L3Absolute_AK4PFchs.txt',
                                   'Summer15_25nsV5_DATA_L2L3Residual_AK4PFchs.txt'
     ] 
#end------------JEC on the fly--------

 
   
process.treeDumper = cms.EDAnalyzer("PKUTreeMaker",
                                    originalNEvents = cms.int32(1),
                                    crossSectionPb = cms.double(1),
                                    targetLumiInvPb = cms.double(1.0),
                                    PKUChannel = cms.string("VW_CHANNEL"),
                                    isGen = cms.bool(False),
				    RunOnMC = cms.bool(runOnMC), 
                                    leptonicVSrc = cms.string("leptonicV"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),   
                                    ak4jetsSrc = cms.string("cleanAK4Jets"),      
                                    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
                                    metSrc = cms.string("slimmedMETs"),
                                    electronIDs = cms.InputTag("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                    looseelectronSrc = cms.InputTag("vetoElectrons"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    conversions = cms.InputTag("reducedEgamma","reducedConversions","PAT"),
                                    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
                                    photonSrc = cms.string("goodPhotons"),  
                                    loosemuonSrc = cms.InputTag("looseMuons"),
                                    hltToken    = cms.InputTag("TriggerResults","","HLT"),
                                    elPaths     = cms.vstring("HLT_Ele27*"),
                                    muPaths     = cms.vstring("HLT_Mu45*"), 
				    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
				    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
				    noiseFilterSelection_EarlyRunsHBHENoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"),
				    noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
				    noiseFilterSelection_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
				    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
				    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
				    noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
				    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
				    noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
				    noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
				    # and the sub-filters
				    noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
    				    noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
    				    noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),
    				    # summary
    				    noiseFilterSelection_metFilters = cms.string('Flag_METFilters')
                                    )


process.analysis = cms.Path(
                            process.goodOfflinePrimaryVertex +
                            process.leptonSequence +
                            process.jetSequence +
                            process.photonSequence +
                            process.treeDumper)

### Source
process.load("VAJets.PKUCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
process.source.fileNames = [
#"/store/user/qili/Bulk/test-WAmatchingnew-v3/150724_013843/0000/TOP-RunIISpring15DR74-00001_93.root",
"/store/user/qili/Bulk/test-WAmatchingnew-v3/150724_013843/0000/TOP-RunIISpring15DR74-00001_90.root",
#"/store/data/Run2015D/SingleElectron/MINIAOD/05Oct2015-v1/10000/00991D45-4E6F-E511-932C-0025905A48F2.root"
]
                       
process.maxEvents.input = 24

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.FwkReport.limit = 99999999

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string("treePKU.root")
                                   )
