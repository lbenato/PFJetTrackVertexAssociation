import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v2/710000/FE472D58-AE73-E811-90E3-24BE05CE1E51.root'
    )
)

process.selectedPrimaryVertices = cms.EDProducer('SelectVertexProducer',
   vertices = cms.InputTag('offlinePrimaryVertices'),
)

process.out = cms.OutputModule("PoolOutputModule",
    #outputCommands = cms.untracked.vstring(['drop *','keep *Vert*_*PrimaryVertices_*_*']),
    outputCommands = cms.untracked.vstring(['keep *_*_*_*']),
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.selectedPrimaryVertices)

process.e = cms.EndPath(process.out)
