import FWCore.ParameterSet.Config as cms

#eleisolationCutString = cms.string("")
#eleisolationCutString = "(pfIsolationVariables().sumChargedHadronPt+max(0.0, pfIsolationVariables().sumNeutralHadronEt+pfIsolationVariables().sumPhotonEt-0.5*pfIsolationVariables().sumPUPt))/pt < 0.10"


tightEleIdLabel = "tight"
mediumEleIdLabel = "medium"
looseEleIdLabel = "loose"
vetoEleIdLabel = "veto"

goodElectrons = cms.EDProducer("PATElectronIdSelector",
    src = cms.InputTag( "slimmedElectrons" ),
    idLabel = cms.string(tightEleIdLabel)
)

vetoElectrons = cms.EDProducer("PATElectronIdSelector",
    src = cms.InputTag( "slimmedElectrons" ),
    idLabel = cms.string(vetoEleIdLabel)
)


#goodElectrons = cms.EDFilter("PATElectronIdSelector",
#                             src = cms.InputTag("tightElectrons"),
#                             cut = cms.string("pt > 20 && abs(eta) < 2.5 "), 
#			     )

#goodElectrons = cms.EDFilter("PATElectronSelector",
#                             src = cms.InputTag("slimmedElectrons"),
#                             cut = cms.string(" pt > 20 && abs(eta) < 2.5 "
#                                              " && ecalDrivenSeed()==1"
#                                              " && abs(1.0/ecalEnergy() - eSuperClusterOverP()/ecalEnergy())<0.05 "
#                                              " && abs(gsfTrack()->dxy())<0.02"
#                                              " && abs(gsfTrack()->dz())<0.1"
#                                              " && gsfTrack()->trackerExpectedHitsInner().numberOfLostHits()==0 "
#                                              " && ( abs(convDist())>0.02 || abs(convDcot())>0.02 ) " 
#                                              " && passConversionVeto()==1 "
#                                              " && ( (isEB() && sigmaIetaIeta()<0.01 && abs(deltaPhiSuperClusterTrackAtVtx())<0.03 && abs(deltaEtaSuperClusterTrackAtVtx())<0.004 && hadronicOverEm()<0.12 ) || " +\
#                                              "      (isEE() && sigmaIetaIeta()<0.03 && abs(deltaPhiSuperClusterTrackAtVtx())<0.02 && abs(deltaEtaSuperClusterTrackAtVtx())<0.005 && hadronicOverEm()<0.10 ))"
#                                              " && " + eleisolationCutString
#                                             )
#                             )

eleSequence = cms.Sequence(goodElectrons+vetoElectrons)
