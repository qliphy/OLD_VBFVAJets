import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
				     SkipEvent = cms.untracked.vstring('ProductNotFound'))
corrJetsOnTheFly = True
runOnMC = True #False
#****************************************************************************************************#
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12' #'74X_mcRun2_asymptotic_v2'#'for version2 miniaod 
elif not(runOnMC):
   process.GlobalTag.globaltag = '76X_dataRun2_v15' #'74X_dataRun2_reMiniAOD_v0' #'74X_dataRun2_Prompt_v4' # for 2015D prompt v4
# 74X_dataRun2_reMiniAOD_v0 for D_05Oct2015

##########					                                                             
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters
hltFiltersProcessName = 'RECO'
if runOnMC:
   hltFiltersProcessName = 'PAT' #'RECO'


process.load("VAJets.PKUCommon.goodMuons_cff")
process.load("VAJets.PKUCommon.goodElectrons_cff")
process.load("VAJets.PKUCommon.goodJets_cff")
process.load("VAJets.PKUCommon.goodPhotons_cff")
process.load("VAJets.PKUCommon.leptonicZ_cff")

# If Update
process.goodMuons.src = "slimmedMuons"
process.goodElectrons.src = "slimmedElectrons"
process.goodAK4Jets.src = "slimmedJets"
process.goodPhotons.src = "slimmedPhotons"

#process.goodOfflinePrimaryVertex = cms.EDFilter("VertexSelector",
#                                       src = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                                       cut = cms.string("chi2!=0 && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"),
#                                       filter = cms.bool(False)
#                                       )

ZBOSONCUT = "pt > 0.0"

process.leptonicVSelector = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("leptonicV"),
                                       cut = cms.string( ZBOSONCUT ), 
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
                                   'Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt',
                                   'Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt',
                                   'Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt'
     ]
else:
   jecLevelsAK4chs = [
                                   'Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt',
                                   'Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt',
                                   'Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt',
                                   'Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt'
     ] 
#end------------JEC on the fly--------

process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
   
process.treeDumper = cms.EDAnalyzer("ZPKUTreeMaker",
                                    originalNEvents = cms.int32(1),
                                    crossSectionPb = cms.double(1),
                                    targetLumiInvPb = cms.double(1.0),
                                    PKUChannel = cms.string("VW_CHANNEL"),
                                    isGen = cms.bool(False),
				    RunOnMC = cms.bool(runOnMC), 
                                    generator =  cms.InputTag("generator"),
                                    pileup  =   cms.InputTag("slimmedAddPileupInfo"),
                                    leptonicVSrc = cms.InputTag("leptonicV"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),   
                                    ak4jetsSrc = cms.InputTag("cleanAK4Jets"),      
#                                    photonSrc = cms.InputTag("goodPhotons"),
                                    photonSrc = cms.InputTag("slimmedPhotons"),
                                    genSrc =  cms.InputTag("prunedGenParticles"),       
                                    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
                                    metSrc = cms.InputTag("slimmedMETs"),
                                    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),  
                                    t1jetSrc = cms.InputTag("slimmedJets"),      
                                    t1muSrc = cms.InputTag("slimmedMuons"),       
                                    electronIDs = cms.InputTag("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                    looseelectronSrc = cms.InputTag("vetoElectrons"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    conversions = cms.InputTag("reducedEgamma","reducedConversions","PAT"),
                                    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
                                    loosemuonSrc = cms.InputTag("looseMuons"),
                                    hltToken    = cms.InputTag("TriggerResults","","HLT"),
                                    elPaths     = cms.vstring("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"),
                                    muPaths     = cms.vstring("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"), 
				    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
				    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
                                    noiseFilterSelection_HBHENoiseIsoFilter = cms.string("Flag_HBHENoiseIsoFilter"),
				    noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHalo2015Filter'),
                                    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
				    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
				    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
                                    full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                                    effAreaChHadFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
                                    effAreaNeuHadFile= cms.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
                                    effAreaPhoFile   = cms.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt")
                                    )


process.analysis = cms.Path(
#                            process.goodOfflinePrimaryVertex +
#			    process.HBHENoiseFilterResultProducer+ #produces HBHE baseline bools
#			    process.ApplyBaselineHBHENoiseFilter+  #reject events based 
#			    process.ApplyBaselineHBHEIsoNoiseFilter+   #reject events based  < 10e-3 mistake rate 
                            process.leptonSequence +
                            process.jetSequence +
#                           process.photonSequence +
                            process.photonIDValueMapProducer*process.treeDumper)

### Source
process.load("VAJets.PKUCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
process.source.fileNames = [
"/store/mc/RunIIFall15MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/12629477-D0B8-E511-929C-0025905C9740.root"
#"/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/FEEA7AEA-12A8-E511-97A6-0025905B860E.root"
]
                       
process.maxEvents.input = 1000

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.MessageLogger.cerr.FwkReport.limit = 99999999

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string("treePKU.root")
                                   )
