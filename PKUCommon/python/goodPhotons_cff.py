import FWCore.ParameterSet.Config as cms



#isolationCutString = cms.string("")
#isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.12"

goodPhotons = cms.EDFilter("PATPhotonSelector",
                             src = cms.InputTag("slimmedPhotons"),
                             cut = cms.string("pt > 15 && abs(eta) < 2.5" 
                                             )
                             )



photonSequence = cms.Sequence(goodPhotons)
