import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process("USER")
task = cms.Task()

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

## Messagge logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #MINIAOD
        '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FEDB9DFE-2542-E811-B39F-0CC47A745284.root',
        #AOD
        #'/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v2/710000/FE472D58-AE73-E811-90E3-24BE05CE1E51.root'
    )
)

isData = ('/store/data/' in process.source.fileNames[0])

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag (process.GlobalTag, 'auto:run2_mc')#default option, but we have the best global tag manually

GT = ''
if isData:
    GT = '80X_dataRun2_2016SeptRepro_v7'
    print "data 2017, Miniaod GT"
elif not(isData):
    GT = '94X_mc2017_realistic_v17'#Moriond17 GT

process.GlobalTag = GlobalTag (process.GlobalTag, GT)

#################################################
## Remake jets
#################################################

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
        src = cms.InputTag("packedGenParticles"), 
        cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
)

task.add(process.packedGenParticlesForJetsNoNu)

## Define GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')
task.add(process.ak4GenJetsNoNu)

## Select charged hadron subtracted packed PF candidates
#process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))


from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

## Define PFJetsCHS
#Test
##process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
#process.ak4PFJets = ak4PFJets.clone(src = AAA, doAreaFastjet = True, jetPtMin = cms.double(15.0),)#HERE!!! add jetPtMin = (fatjet_ptmin)
#task.add(process.ak4PFJets)


### New1
##step 1: new PF collection
process.pfNew1CHS = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)")## new selections for PF Candidates
)
task.add(process.pfNew1CHS)

##step2: recluster jets
process.ak4PFJetsNew1CHS = ak4PFJets.clone(src = 'pfNew1CHS',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew1CHS)

### New2
process.pfNew2CHS = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("fromPV()>0 || charge() ==0") ## you can change it here #loosely matched to the PV that it's matched to-->it's matched to 1 PV
)
task.add(process.pfNew2CHS)

process.ak4PFJetsNew2CHS = ak4PFJets.clone(src = 'pfNew2CHS',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew2CHS)


### New3
process.pfNew3CHS = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )") ## you can change it here
)

task.add(process.pfNew3CHS)

process.ak4PFJetsNew3CHS = ak4PFJets.clone(src = 'pfNew3CHS',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew3CHS)


### New4
process.pfNew4CHS = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0 ") ## you can change it here
)
task.add(process.pfNew4CHS)


process.ak4PFJetsNew4CHS = ak4PFJets.clone(src = 'pfNew4CHS',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew4CHS)




#################################################
## Remake PAT jets
#################################################

## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]

from PhysicsTools.PatAlgos.tools.jetTools import *
## Add PAT jet collection based on the above-defined ak4PFJetsCHS
'''addJetCollection(
    process,
    labelName = 'AK4PFCHS',
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.4
)'''
addJetCollection(
    process,
    labelName = 'New1CHS',
    jetSource = cms.InputTag('ak4PFJetsNew1CHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('pfNew1CHS'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'New2CHS',
    jetSource = cms.InputTag('ak4PFJetsNew2CHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('pfNew2CHS'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'New3CHS',
    jetSource = cms.InputTag('ak4PFJetsNew3CHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('pfNew3CHS'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'New4CHS',
    jetSource = cms.InputTag('ak4PFJetsNew4CHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('pfNew4CHS'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.4
)

process.selectedPatJets = cms.EDFilter("PATJetSelector",
    cut = cms.string('pt > 15'),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("slimmedJets")
)
task.add(process.selectedPatJets)

getattr(process,'selectedPatJetsNew1CHS').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew2CHS').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew3CHS').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew4CHS').cut = cms.string('pt > 15')
###cut also slimmedJets!!!



from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False), # while the timing of this is not reliable in unscheduled mode, it still helps understanding what was actually run
        allowUnscheduled = cms.untracked.bool(True)
)


##########################################
### Here: pieces taken from the other config

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('aaaaa.root' if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)



process.ntuple = cms.EDAnalyzer('Ntuplizer',
  #jetSet = cms.PSet(
  #  jets = cms.InputTag('selectedPatJets'),
  #  met = cms.InputTag('slimmedMETs'),
  #  jetid = cms.int32(0),
  #  jet1pt = cms.double(15.),
  #  jet2pt = cms.double(15.),
  #  jeteta = cms.double(5.2),
  #  vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  #  rho = cms.InputTag('fixedGridRhoFastjetAll'),
  #),
  jets = cms.InputTag('selectedPatJets'),
  jetpt = cms.double(15.),
  pileup = cms.InputTag('slimmedAddPileupInfo'),
  vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  jets1 = cms.InputTag('selectedPatJetsNew1CHS'),
  jets2 = cms.InputTag('selectedPatJetsNew2CHS'),
  jets3 = cms.InputTag('selectedPatJetsNew3CHS'),
  jets4 = cms.InputTag('selectedPatJetsNew4CHS'),
  genjets = cms.InputTag('slimmedGenJets'),
  rhosrc = cms.InputTag('fixedGridRhoFastjetAll'),
  verbose = cms.bool(False),
)


#process.seq = cms.Sequence(
#    process.jet
#)


process.seq = cms.Sequence(
    #PF candidates
    process.pfNew1CHS *
    process.pfNew2CHS *
    process.pfNew3CHS *
    process.pfNew4CHS *
    #reco jets
    process.ak4PFJetsNew1CHS *
    process.ak4PFJetsNew2CHS *
    process.ak4PFJetsNew3CHS *
    process.ak4PFJetsNew4CHS *
    #pat jets
    process.patJetsNew1CHS *
    process.patJetsNew2CHS *
    process.patJetsNew3CHS *
    process.patJetsNew4CHS *
    #selected pat jets
    process.selectedPatJets *
    process.selectedPatJetsNew1CHS *
    process.selectedPatJetsNew2CHS *
    process.selectedPatJetsNew3CHS *
    process.selectedPatJetsNew4CHS *
    #analyzer
    process.ntuple
)



process.p = cms.Path(
    process.seq
)

#####
print(process.p)

#####################


## Output file
## if needed
'''
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])#(['drop *','keep patJets_selectedPatJetsAK5PFCHS_*_*'])
)
process.endpath= cms.EndPath(process.OUT)
'''

open('dump_ReclusterJets.py','w').write(process.dumpPython())

process.p.associate(task)
process.p.associate(process.patAlgosToolsTask)
