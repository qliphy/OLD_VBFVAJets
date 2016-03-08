import FWCore.ParameterSet.Config as cms

tightMuIdLabel = "tight"
looseMuIdLabel = "loose"

goodMuons = cms.EDProducer("PATMuonIdSelector",
    src = cms.InputTag( "slimmedMuons" ),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    idLabel = cms.string(tightMuIdLabel)
)   

looseMuons = cms.EDProducer("PATMuonIdSelector",
    src = cms.InputTag( "slimmedMuons" ),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    idLabel = cms.string(looseMuIdLabel)
)   

#goodMuons = cms.EDFilter("PATMuonIdSelector",
#                             src = cms.InputTag("tightMuons"),
#                             cut = cms.string("pt > 20 && abs(eta) < 2.5 ") 
#                         )

#isolationCutString = cms.string("")
#isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.12"

#goodMuons = cms.EDFilter("PATMuonSelector",
#                             src = cms.InputTag("slimmedMuons"),
#                             cut = cms.string("pt > 30 && abs(eta) < 2.4" 
#                                              "&& isGlobalMuon && isPFMuon "
#                                              " && globalTrack().normalizedChi2<10"
#                                              " && globalTrack().hitPattern().numberOfValidMuonHits>0"
#                                              " && numberOfMatchedStations() > 1"
##                                              " && abs(muonBestTrack()->dxy(vertex->position())) < 0.2 "   
#                                              " && dB() < 0.2 "
##                                              " && abs(muonBestTrack()->dz(vertex->position())) < 0.5 "   
#                                              " && globalTrack().hitPattern().numberOfValidPixelHits>0"
#                                              " && numberOfMatchedStations>1"
#                                              " && globalTrack().hitPattern().trackerLayersWithMeasurement>5"
#                                              " && " + isolationCutString
#                                             )
#                             )



muSequence = cms.Sequence(goodMuons+looseMuons)
