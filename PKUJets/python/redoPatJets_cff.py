import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS

chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

ak4PFJetsCHS.src = cms.InputTag('chs') 


patJetCorrFactorsAK4 = cms.EDProducer("JetCorrFactorsProducer",
                                      src = cms.InputTag("ak4PFJetsCHS"),
                                      emf = cms.bool(False),
                                      extraJPTOffset = cms.string('L1FastJet'),
                                      primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                      levels = cms.vstring('L1FastJet',
                                                           'L2Relative',
                                                           'L3Absolute'),
                                      useNPV = cms.bool(True),
                                      rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                      useRho = cms.bool(True),
                                      payload = cms.string('AK4PFchs'),
                                      flavorType = cms.string('J')
                                      )

patJetsAK4 = cms.EDProducer("PATJetProducer",
                            addJetCharge = cms.bool(False),
                            addGenJetMatch = cms.bool(False),
                            embedGenJetMatch = cms.bool(False),
                            addAssociatedTracks = cms.bool(False),
                            addBTagInfo = cms.bool(False),
                            partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
                            addGenPartonMatch = cms.bool(False),
                            JetPartonMapSource = cms.InputTag(""),
                            resolutions = cms.PSet(),
                            genPartonMatch = cms.InputTag(""),
                            addTagInfos = cms.bool(False),
                            addPartonJetMatch = cms.bool(False),
                            embedGenPartonMatch = cms.bool(False),
                            efficiencies = cms.PSet(),
                            genJetMatch = cms.InputTag(""),
                            useLegacyJetMCFlavour = cms.bool(False),
                            jetSource = cms.InputTag("ak4PFJetsCHS"),
                            addEfficiencies = cms.bool(False),
                            discriminatorSources = cms.VInputTag(),
                            trackAssociationSource = cms.InputTag(""),
                            tagInfoSources = cms.VInputTag(cms.InputTag("")),
                            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK4")),
                            embedPFCandidates = cms.bool(False),
                            addJetFlavourInfo = cms.bool(False),
                            addResolutions = cms.bool(False),
                            getJetMCFlavour = cms.bool(False),
                            addDiscriminators = cms.bool(False),
                            jetChargeSource = cms.InputTag(""),
                            JetFlavourInfoSource = cms.InputTag(""),
                            addJetCorrFactors = cms.bool(True),
                            jetIDMap = cms.InputTag(""),
                            addJetID = cms.bool(False)
                            )

selectedPatJetsAK4 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsAK4"),
    cut = cms.string('pt > 20')
)

redoPatJets = cms.Sequence(chs + ak4PFJetsCHS + patJetCorrFactorsAK4 +
                           patJetsAK4 + 
                           selectedPatJetsAK4)
