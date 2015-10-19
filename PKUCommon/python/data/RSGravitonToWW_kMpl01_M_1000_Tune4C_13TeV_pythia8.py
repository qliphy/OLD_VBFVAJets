import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ['/store/user/qili/Bulk/test-LOWVA-v3/150714_002636/0000/TOP-RunIISpring15DR74-00001_10.root'] );

secFiles.extend( [
               ] )

