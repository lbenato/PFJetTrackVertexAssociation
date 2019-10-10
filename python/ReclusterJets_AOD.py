import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process("ReclusterAOD")
task = cms.Task()

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Messagge logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #MINIAOD
        #'/store/mc/RunIIFall17MiniAODv2/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FEDB9DFE-2542-E811-B39F-0CC47A745284.root',
        #AOD
        #'/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v2/710000/FE472D58-AE73-E811-90E3-24BE05CE1E51.root'
        #AOD with selected vertices
        #'file:myOutputFile.root'
        #AOD with general tracks:
        #'file:aodsim_generalTracks.root',
        #AOD with selected tracks:
        'file:aodsim_generalTracks_eta_2p5.root',
    )
)

isData = ('/store/data/' in process.source.fileNames[0])

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag (process.GlobalTag, 'auto:run2_mc')#default option, but we have the best global tag manually

GT = ''
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
if isData:
    GT = '94X_dataRun2_v11'
    print "data 2017,AOD GT"
elif not(isData):
    GT = '94X_mc2017_realistic_v17'#Moriond17 GT

process.GlobalTag = GlobalTag (process.GlobalTag, GT)

#-----------------------#
#       PF tracks       #
#-----------------------#

###Filter tracks
#process.selectTracks = cms.EDFilter("TrackSelector",
#    src = cms.InputTag("generalTracks"),
#    cut = cms.string('abs(eta) < 2.5')
#)
#task.add(process.selectTracks)

###############################################################################################
##Hack PF
#from RecoParticleFlow.PFTracking.pfTrack_cfi import *
#process.pfTrackSelect = pfTrack.clone(TkColList = cms.VInputTag(cms.InputTag("selectTracks")))
#task.add(process.pfTrackSelect)

#from RecoParticleFlow.PFTracking.pfTrackElec_cfi import *
#process.pfTrackElec = pfTrackElec.clone(PFRecTrackLabel = cms.InputTag("pfTrackSelect"),)
#task.add(process.pfTrackElec)
################################################################################################
#-----------------------#
#    VERTEX PRODUCER    #
#-----------------------#

process.selectedPrimaryVerticesFakeRejected = cms.EDProducer('SelectVertexProducer',
   vertices = cms.InputTag('offlinePrimaryVertices'),
)
task.add(process.selectedPrimaryVerticesFakeRejected)

##tentativo random????????????????????????????????????????????????????????????????????????????????
#process.offlinePrimaryVertices = cms.EDProducer('SelectVertexProducer',
#   vertices = cms.InputTag('offlinePrimaryVertices'),
#)
#task.add(process.offlinePrimaryVertices)

#-----------------------#
#    VERTEX FILTER      #
#-----------------------#
#It works, but it's more restrictive. To be reconsidered later.
'''
import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)
task.add(process.primaryVertexFilter)#test!
'''
#-----------------------#
#    PF CANDIDATES      #
#-----------------------#

#Patify
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *

#Puppi PF must be computed
process.load('CommonTools.PileupAlgos.Puppi_cff')
task.add(process.puppi)

process.pfNoLepPUPPI = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"), 
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )
task.add(process.pfNoLepPUPPI)

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
addToProcessAndTask('puppiNoLep', process.puppi.clone(), process, task)

#Vertex association
from CommonTools.RecoAlgos.sortedPFPrimaryVertices_cfi import sortedPFPrimaryVertices
process.primaryVertexAssociation = sortedPFPrimaryVertices.clone(
  jets = cms.InputTag("ak4PFJets"),
  vertices = cms.InputTag("offlinePrimaryVertices"),
  #vertices = cms.InputTag('selectedPrimaryVerticesLisa','','ReclusterAOD'),
  qualityForPrimary = cms.int32(2),
  produceSortedVertices = cms.bool(False),
  producePileUpCollection  = cms.bool(False),  
  produceNoPileUpCollection = cms.bool(False)
)
task.add(process.primaryVertexAssociation)

#Vertex slimmer
process.offlineSlimmedPrimaryVertices = cms.EDProducer("PATVertexSlimmer",
    score = cms.InputTag("primaryVertexAssociation","original"),
    src = cms.InputTag("offlinePrimaryVertices")
    #src = cms.InputTag('selectedPrimaryVerticesLisa','','ReclusterAOD')
)
task.add(process.offlineSlimmedPrimaryVertices)

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)

from PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi import *
from CommonTools.ParticleFlow.pfNoPileUpJME_cff import *
from CommonTools.ParticleFlow.PFBRECO_cff import pfPileUpPFBRECO, pfNoPileUpPFBRECO
from RecoParticleFlow.PFProducer.chargedHadronPFTrackIsolation_cfi import * #not really loading anything

#necessary..?
#process.load("CommonTools.ParticleFlow.pfNoPileUpJME_cff")
#from CommonTools.ParticleFlow.pfPileUp_cfi  import pfPileUp as _pfPileUp
#from CommonTools.ParticleFlow.TopProjectors.pfNoPileUp_cfi import pfNoPileUp as _pfNoPileUp
#from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import *
#process.pfPileUpJME = _pfPileUp.clone(PFCandidates='particleFlowPtrs',
#                              Vertices = 'goodOfflinePrimaryVertices',
#                              checkClosestZVertex = False )
#process.pfNoPileUpJME = _pfNoPileUp.clone(topCollection = 'pfPileUpJME',
#                                  bottomCollection = 'particleFlowPtrs' )
#task.add(process.pfPileUpJME)
#task.add(process.pfNoPileUpJME)

#Prepare packedPFCandidates with generalTracks
process.packedPFCandidates = packedPFCandidates.clone(
    originalVertices = cms.InputTag("offlinePrimaryVertices"),
    #originalVertices = cms.InputTag('selectedPrimaryVerticesLisa','','ReclusterAOD'),
    inputVertices  = cms.InputTag("offlineSlimmedPrimaryVertices"), 
    originalTracks = cms.InputTag("generalTracks"),
    vertexAssociator = cms.InputTag("primaryVertexAssociation","original")
    )
task.add(process.packedPFCandidates)


process.primaryVertexAssociationFakeRejected = sortedPFPrimaryVertices.clone(
  jets = cms.InputTag("ak4PFJets"),
  #vertices = cms.InputTag("offlinePrimaryVertices"),
  vertices = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
  qualityForPrimary = cms.int32(2),
  produceSortedVertices = cms.bool(False),
  producePileUpCollection  = cms.bool(False),  
  produceNoPileUpCollection = cms.bool(False)
)
task.add(process.primaryVertexAssociationFakeRejected)


#Vertex slimmer
process.offlineSlimmedPrimaryVerticesFakeRejected = cms.EDProducer("PATVertexSlimmer",
    score = cms.InputTag("primaryVertexAssociationFakeRejected","original"),
    #src = cms.InputTag("offlinePrimaryVertices")
    src = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD')
)
task.add(process.offlineSlimmedPrimaryVerticesFakeRejected)

process.packedPFCandidatesFakeRejected = packedPFCandidates.clone(
    #originalVertices = cms.InputTag("offlinePrimaryVertices"),
    originalVertices = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    inputVertices  = cms.InputTag("offlineSlimmedPrimaryVerticesFakeRejected"), 
    originalTracks = cms.InputTag("generalTracks"),
    vertexAssociator = cms.InputTag("primaryVertexAssociationFakeRejected","original")
    )
task.add(process.packedPFCandidatesFakeRejected)

#---------------------------------------------#
#  PF candidates input for jet reclustering   #
#---------------------------------------------#

#Default: CHS jets

#New0: PF w/o CHS

### New1
##step 1: new PF collection
process.pfNew1 = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),
      cut = cms.string("(abs(dz()) < cosh(eta())*(0.02+0.01/pt())*5)")
)
task.add(process.pfNew1)

### New2
process.pfNew2 = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("fromPV()>0 || charge() ==0") ## you can change it here #loosely matched to the PV that it's matched to-->it's matched to 1 PV
)
task.add(process.pfNew2)

### New3
process.pfNew3 = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )") ## you can change it here
)

task.add(process.pfNew3)

### New4
process.pfNew4 = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0 ") ## you can change it here
)
task.add(process.pfNew4)


### New5, PUPPI-like
PVUsedInFit = "3"
PVTight = "2"
PVLoose = "1"

process.pfNew5 = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidates"),                                     
      cut = cms.string("fromPV()=="+PVUsedInFit+" || ( (fromPV()=="+PVTight+" || fromPV()=="+PVLoose+") && abs(dz()) < 0.3 )") ## you can change it here
)
task.add(process.pfNew5)

### Vertex association

process.pfFakeRejectedCHS = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidatesFakeRejected"),#%test                                  
      cut = cms.string("fromPV()>0 || charge() ==0")
)
task.add(process.pfFakeRejectedCHS)

process.pfFakeRejectedPuppi = cms.EDFilter("CandPtrSelector",
      src = cms.InputTag("packedPFCandidatesFakeRejected"),#%test                                  
      cut = cms.string("fromPV()=="+PVUsedInFit+" || ( (fromPV()=="+PVTight+" || fromPV()=="+PVLoose+") && abs(dz()) < 0.3 )")
)
task.add(process.pfFakeRejectedPuppi)

### New7: selectTrack w/ non-chs JEC
#so, basically, packedPFCandidates


#-----------------------#
#  Secondary Vertices   #
#-----------------------#

#Secondary vertices needed as input for pat Jets
process.slimmedSecondaryVertices = cms.EDProducer("PATSecondaryVertexSlimmer",
    src = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    packedPFCandidates = cms.InputTag("packedPFCandidates")
)
task.add(process.slimmedSecondaryVertices)

process.slimmedSecondaryVerticesSelect = cms.EDProducer("PATSecondaryVertexSlimmer",
    src = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    packedPFCandidates = cms.InputTag("packedPFCandidates")
)
task.add(process.slimmedSecondaryVerticesSelect)

#-------------------------#
#  Packed Gen Particles   #
#-------------------------#

#Patify packedGenParticles
process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
    select = cms.vstring('drop  *', 
        '++keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15', 
        'drop   status == 2', 
        'keep++ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)', 
        'drop status == 1', 
        'keep+ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)', 
        'keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15', 
        'keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16', 
        '+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())', 
        '+keep abs(pdgId) == 11 && status == 1 && (pt > 3 || isPromptFinalState())', 
        'keep++ abs(pdgId) == 15', 
        'drop  status > 30 && status < 70 ', 
        'drop  pdgId == 21 && pt < 5', 
        'drop   status == 2 && abs(pdgId) == 21', 
        'keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 6 || abs(pdgId) == 37 ', 
        'keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ', 
        '+keep abs(pdgId) == 13 && status == 1', 
        'keep (4 <= abs(pdgId) <= 5)', 
        'keep (1 <= abs(pdgId) <= 3 || pdgId = 21) & (status = 2 || status = 11 || status = 71 || status = 72) && pt>5', 
        'keep+ abs(pdgId) == 333', 
        'keep+ abs(pdgId) == 9920443 || abs(pdgId) == 9042413 || abs(pdgId) == 9000443', 
        'keep+ abs(pdgId) == 443 || abs(pdgId) == 100443 || abs(pdgId) == 10441 || abs(pdgId) == 20443 || abs(pdgId) == 445 || abs(pdgId) == 30443', 
        'keep+ abs(pdgId) == 553 || abs(pdgId) == 100553 || abs(pdgId) == 200553 || abs(pdgId) == 10551 || abs(pdgId) == 20553 || abs(pdgId) == 555', 
        'keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433', 
        'keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543', 
        'keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)', 
        'keep pdgId = 2212', 
        'keep status == 3 || ( 21 <= status <= 29) || ( 11 <= status <= 19)', 
        'keep isHardProcess() || fromHardProcessFinalState() || fromHardProcessDecayed() || fromHardProcessBeforeFSR() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())'),
    src = cms.InputTag("prunedGenParticlesWithStatusOne")
)
task.add(process.prunedGenParticles)

process.prunedGenParticlesWithStatusOne = cms.EDProducer("GenParticlePruner",
    select = cms.vstring('drop  *', 
        '++keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15', 
        'drop   status == 2', 
        'keep++ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)', 
        'drop status == 1', 
        'keep+ (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)', 
        'keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15', 
        'keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16', 
        '+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())', 
        '+keep abs(pdgId) == 11 && status == 1 && (pt > 3 || isPromptFinalState())', 
        'keep++ abs(pdgId) == 15', 
        'drop  status > 30 && status < 70 ', 
        'drop  pdgId == 21 && pt < 5', 
        'drop   status == 2 && abs(pdgId) == 21', 
        'keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 6 || abs(pdgId) == 37 ', 
        'keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ', 
        '+keep abs(pdgId) == 13 && status == 1', 
        'keep (4 <= abs(pdgId) <= 5)', 
        'keep (1 <= abs(pdgId) <= 3 || pdgId = 21) & (status = 2 || status = 11 || status = 71 || status = 72) && pt>5', 
        'keep+ abs(pdgId) == 333', 
        'keep+ abs(pdgId) == 9920443 || abs(pdgId) == 9042413 || abs(pdgId) == 9000443', 
        'keep+ abs(pdgId) == 443 || abs(pdgId) == 100443 || abs(pdgId) == 10441 || abs(pdgId) == 20443 || abs(pdgId) == 445 || abs(pdgId) == 30443', 
        'keep+ abs(pdgId) == 553 || abs(pdgId) == 100553 || abs(pdgId) == 200553 || abs(pdgId) == 10551 || abs(pdgId) == 20553 || abs(pdgId) == 555', 
        'keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433', 
        'keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543', 
        'keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)', 
        'keep pdgId = 2212', 
        'keep status == 3 || ( 21 <= status <= 29) || ( 11 <= status <= 19)', 
        'keep isHardProcess() || fromHardProcessFinalState() || fromHardProcessDecayed() || fromHardProcessBeforeFSR() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())', 
        'keep    status == 1'),
    src = cms.InputTag("genParticles")
)
task.add(process.prunedGenParticlesWithStatusOne)

process.packedGenParticles = cms.EDProducer("PATPackedGenParticleProducer",
    inputCollection = cms.InputTag("prunedGenParticlesWithStatusOne"),
    inputOriginal = cms.InputTag("genParticles"),
    map = cms.InputTag("prunedGenParticles"),
    maxRapidity = cms.double(6)
)
task.add(process.packedGenParticles)


#---------------------#
#      Gen Jets       #
#---------------------#

## Filter out neutrinos from packed GenParticles
#for miniaod
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
        src = cms.InputTag("packedGenParticles"),
        cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
)
task.add(process.packedGenParticlesForJetsNoNu)

#for aod
process.load('RecoJets.Configuration.GenJetParticles_cff')
task.add(process.genParticlesForJetsNoNu)

## Define GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsNoNu = ak4GenJets.clone(
    src = 'genParticlesForJetsNoNu',
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    )#different in miniaod!!!#'packedGenParticlesForJetsNoNu')
task.add(process.ak4GenJetsNoNu)

process.ak4GenJetsNoNuFakeRejected = ak4GenJets.clone(
    src = 'genParticlesForJetsNoNu',
    srcPVs = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    )#different in miniaod!!!#'packedGenParticlesForJetsNoNu')
task.add(process.ak4GenJetsNoNuFakeRejected)

#---------------------#
#     Reco Jets       #
#---------------------#

## recluster reco jets

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJetsNew1 = ak4PFJets.clone(src = 'pfNew1',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      srcPVs = cms.InputTag("offlinePrimaryVertices"),
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew1)

process.ak4PFJetsNew2 = ak4PFJets.clone(src = 'pfNew2',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      srcPVs = cms.InputTag("offlinePrimaryVertices"),
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew2)

process.ak4PFJetsNew3 = ak4PFJets.clone(src = 'pfNew3',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      srcPVs = cms.InputTag("offlinePrimaryVertices"),
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew3)

process.ak4PFJetsNew4 = ak4PFJets.clone(src = 'pfNew4',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      srcPVs = cms.InputTag("offlinePrimaryVertices"),
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew4)

process.ak4PFJetsNew5 = ak4PFJets.clone(src = 'pfNew5',
      #jetPtMin = cms.double(15.0),#Lisa 27.05.2019
      srcPVs = cms.InputTag("offlinePrimaryVertices"),
      doAreaFastjet = True)
task.add(process.ak4PFJetsNew5)

### Fake vertex rejection:

process.ak4PFJetsFakeRejectedCHS = ak4PFJets.clone(src = 'pfFakeRejectedCHS',
      srcPVs = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
      doAreaFastjet = True)
task.add(process.ak4PFJetsFakeRejectedCHS)

process.ak4PFJetsFakeRejectedPuppi = ak4PFJets.clone(src = 'pfFakeRejectedPuppi',
      srcPVs = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
      doAreaFastjet = True)
task.add(process.ak4PFJetsFakeRejectedPuppi)


##needed??? Ask Andreas!
#according to that, it's not relevant: https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/PhysicsTools/PatAlgos/python/tools/jetTools.py#L1242-L1244
#process.jetTracksAssociatorAtVertexSelect = cms.EDProducer("JetTracksAssociatorAtVertex",
#    coneSize = cms.double(0.4),
#    jets = cms.InputTag("ak4PFJetsSelect"),
#    pvSrc = cms.InputTag("offlinePrimaryVertices"),
#    tracks = cms.InputTag("selectTracks"),#HERE?
#    useAssigned = cms.bool(False)
#)
#task.add(process.jetTracksAssociatorAtVertexSelect)


#---------------------#
#      Pat Jets       #
#---------------------#

#patJetPartons needed for pat jets
process.patJetPartons = cms.EDProducer("HadronAndPartonSelector",
    fullChainPhysPartons = cms.bool(True),
    particles = cms.InputTag("genParticles"),
    partonMode = cms.string('Auto'),
    src = cms.InputTag("generator")
)
task.add(process.patJetPartons)

## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]

from PhysicsTools.PatAlgos.tools.jetTools import *
#
#addJetCollection(
#    process,
#    labelName = 'Select',
#    jetSource = cms.InputTag('ak4PFJetsSelect'),
#    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
#    pfCandidates = cms.InputTag('packedPFCandidates'),
#    svSource = cms.InputTag('slimmedSecondaryVertices'),
#    btagDiscriminators = bTagDiscriminators,
#    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
#    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
#    genParticles = cms.InputTag('prunedGenParticles'),
#    #jetTrackAssociation = True,#according to that, it's not relevant: https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/PhysicsTools/PatAlgos/python/tools/jetTools.py#L1242-L1244
#    algo = 'AK',
#    rParam = 0.4
#)
#getattr(process,'selectedPatJetsSelect').cut = cms.string('pt > 15')

addJetCollection(
    process,
    labelName = 'CHS',
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('packedPFCandidates'),#!!is this correct?most likely not!!!!!!!!!!!!check if it changes
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

#PF w/o CHS
addJetCollection(
    process,
    labelName = 'New0',
    jetSource = cms.InputTag('ak4PFJets'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('particleFlow'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)


addJetCollection(
    process,
    labelName = 'New1',
    jetSource = cms.InputTag('ak4PFJetsNew1'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfNew1'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'New2',
    jetSource = cms.InputTag('ak4PFJetsNew2'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfNew2'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'New3',
    jetSource = cms.InputTag('ak4PFJetsNew3'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfNew3'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

#Important! This collection behaves more like PF w/o chs!
addJetCollection(
    process,
    labelName = 'New4',
    jetSource = cms.InputTag('ak4PFJetsNew4'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfNew4'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'New5',
    jetSource = cms.InputTag('ak4PFJetsNew5'),
    pvSource = cms.InputTag('offlinePrimaryVertices'),
    #pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfNew5'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'FakeRejectedCHS',
    jetSource = cms.InputTag('ak4PFJetsFakeRejectedCHS'),
    pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfFakeRejectedCHS'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNuFakeRejected'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'FakeRejectedPuppi',
    jetSource = cms.InputTag('ak4PFJetsFakeRejectedPuppi'),
    pvSource = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
    pfCandidates = cms.InputTag('pfFakeRejectedPuppi'),
    svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNuFakeRejected'),
    genParticles = cms.InputTag('genParticles'),
    algo = 'AK',
    rParam = 0.4
)

#process.patJetCorrFactorsFakeRejected.primaryVertices = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD')

#process.selectedPatJets = cms.EDFilter("PATJetSelector",
#    cut = cms.string('pt > 15'),
#    cutLoose = cms.string(''),
#    nLoose = cms.uint32(0),
#    #src = cms.InputTag("slimmedJets")
#    src = cms.InputTag("updatedPatJetsAK4PFCHS")#slimmedJets with updated JEC
#)
#task.add(process.selectedPatJets)

getattr(process,'selectedPatJetsCHS').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew0').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew1').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew2').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew3').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew4').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsNew5').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsFakeRejectedCHS').cut = cms.string('pt > 15')
#getattr(process,'selectedPatJetsFakeRejected').cut = cms.string('pt > 15')
getattr(process,'selectedPatJetsFakeRejectedPuppi').cut = cms.string('pt > 15')
###cut also slimmedJets!!!

from PhysicsTools.PatAlgos.tools.pfTools import *
## DO NOT USE THAT! ##
## Adapt primary vertex collection
#adaptPVs(process, pvCollection=cms.InputTag('offlinePrimaryVertices'))
#adaptPVs(process, pvCollection=cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'))

#---------------------#
# Slimmed PileUp Info #
#---------------------#

process.slimmedAddPileupInfo = cms.EDProducer("PileupSummaryInfoSlimmer",
    keepDetailedInfoFor = cms.vint32(0),
    src = cms.InputTag("addPileupInfo")
)
task.add(process.slimmedAddPileupInfo)

#2. 'slimmedGenJets'
## FOR SOME REASON THIS DOES NOT WORK
#from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import *
##process.slimmedGenJets = cms.EDProducer("PATGenJetSlimmer",
##    clearDaughters = cms.bool(False),
##    cut = cms.string('pt > 8'),
##    cutLoose = cms.string(''),
##    dropSpecific = cms.bool(False),
##    nLoose = cms.uint32(0),
##    packedGenParticles = cms.InputTag("packedGenParticles"),
##    src = cms.InputTag("ak4GenJetsNoNu")
##)
#process.slimmedGenJets = slimmedGenJets.clone()
#task.add(process.slimmedGenJets)

##process.slimmedGenJetsFlavourInfos = cms.EDProducer("GenJetFlavourInfoPreserver",
##    genJetFlavourInfos = cms.InputTag("ak4GenJetFlavourInfos"),
##    genJets = cms.InputTag("ak4GenJetsNoNu"),
##    slimmedGenJetAssociation = cms.InputTag("slimmedGenJets","slimmedGenJetAssociation"),
##    slimmedGenJets = cms.InputTag("slimmedGenJets")
##)
##task.add(process.slimmedGenJetsFlavourInfos)
##3. 'fixedGridRhoFastjetAll' --> included already


#################################################
## Previous attempts of remaking PF candidates
#################################################

#To be revised!!!
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
####it works but the content is not what I expect....

#from PhysicsTools.PatAlgos.tools.pfTools import *
#addPFCandidates(process, src="particleFlow", patLabel="LLL", cut="")
#This works but PFCandidates do not have the same content as packedPFCandidates, for example no dz --> can't prepare NewCHS1 etc

#from PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi import *
#process.mypackedPFCandidates = packedPFCandidates.clone(PuppiNoLepSrc = cms.InputTag("particleFlow"), PuppiSrc = cms.InputTag("particleFlow"))
#puppi not found!
#task.add(process.mypackedPFCandidates)


'''
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False), # while the timing of this is not reliable in unscheduled mode, it still helps understanding what was actually run
        allowUnscheduled = cms.untracked.bool(True)
)
'''

##########################################
### Here: pieces taken from the other config

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('aaaaa.root' if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)



process.ntuple = cms.EDAnalyzer('Ntuplizer',
  jets = cms.InputTag('selectedPatJetsCHS'),
  jetpt = cms.double(15.),
  pileup = cms.InputTag('slimmedAddPileupInfo'),
  #vertices = cms.InputTag('selectedPrimaryVerticesFakeRejected','','ReclusterAOD'),
  #vertices = cms.InputTag('selectedPrimaryVertices','','OWNPARTICLES'),
  vertices = cms.InputTag('offlinePrimaryVertices'),#!
  pfcandidates = cms.InputTag('packedPFCandidates'),
  pfcandidatesfakerejected = cms.InputTag('packedPFCandidatesFakeRejected'),
  jets0 = cms.InputTag('selectedPatJetsNew0'),
  jets1 = cms.InputTag('selectedPatJetsNew1'),#New1'),#!!
  jets2 = cms.InputTag('selectedPatJetsNew2'),
  jets3 = cms.InputTag('selectedPatJetsNew3'),
  jets4 = cms.InputTag('selectedPatJetsNew4'),
  jets5 = cms.InputTag('selectedPatJetsNew5'),#Puppi
  jets6 = cms.InputTag('selectedPatJetsFakeRejectedCHS'),
  jets7 = cms.InputTag('selectedPatJetsFakeRejectedPuppi'),
  genjets = cms.InputTag('ak4GenJets'),#slimmedGenJets do not work here, for unknown reasons
  rhosrc = cms.InputTag('fixedGridRhoFastjetAll'),
  verbose = cms.bool(False),
)

task.add(process.patAlgosToolsTask)

process.seq = cms.Sequence(
    ##analyzer
    process.selectedPrimaryVerticesFakeRejected *
    process.ntuple
    #process.offlinePrimaryVertices
)


process.lisa = cms.Path(
    process.seq
)

#####
print(process.lisa)

#####################


## Output file
## if needed

'''
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_AOD.root'),
    #outputCommands = cms.untracked.vstring(['keep *'])#(['drop *','keep patJets_selectedPatJetsAK5PFCHS_*_*'])
    outputCommands = cms.untracked.vstring(['drop *','keep recoVertex_*_*_*'])
    #outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)
process.endpath= cms.EndPath(process.OUT, patAlgosToolsTask)
'''

open('dump_ReclusterJets_AOD.py','w').write(process.dumpPython())

process.lisa.associate(task)
