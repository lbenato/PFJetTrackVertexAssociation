// -*- C++ -*-
//
// Package:    PFJetTrackVertexAssociation/Ntuplizer
// Class:      Ntuplizer
// 
/**\class Ntuplizer Ntuplizer.cc Analyzer/PFJetTrackVertexAssociation/plugins/Ntuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Negin Shafiei, Lisa Benato
//         Created:  Wed, 29 May 2019 10:00:00 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "DataFormats/JetReco/interface/JPTJet.h"

#include <vector>
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "Objects.h"
#include "ObjectsFormat.h"
//#include "JetAnalyzer.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Ntuplizer(const edm::ParameterSet&);
      ~Ntuplizer();
      virtual void GenJetAnalyzer(std::vector<pat::Jet>&);
      virtual bool isLooseJet(pat::Jet&);
      virtual bool isTightJet(pat::Jet&);
      virtual double vertex_ptmax2(const reco::Vertex&);
      virtual bool select(const reco::Vertex &, int);

  // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  
  //edm::ParameterSet JetPSet;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  double jetpt;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT< std::vector<reco::Vertex> > PVToken_;
  edm::EDGetTokenT< std::vector<pat::PackedCandidate> > PFCandToken_;
  //edm::EDGetTokenT< std::vector<pat::PackedCandidate> > PFCandFakeRejectedToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken0_;
  edm::EDGetTokenT<pat::JetCollection> jetToken1_;
  edm::EDGetTokenT<pat::JetCollection> jetToken2_;
  edm::EDGetTokenT<pat::JetCollection> jetToken3_;
  edm::EDGetTokenT<pat::JetCollection> jetToken4_;
  edm::EDGetTokenT<pat::JetCollection> jetToken5_;
  edm::EDGetTokenT<pat::JetCollection> jetToken6_;
  edm::EDGetTokenT<pat::JetCollection> jetToken7_;
  //For reco jets:
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken0_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken1_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken2_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken3_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken4_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken5_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken6_;
  //edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken7_;
  edm::EDGetTokenT< std::vector<reco::GenJet> > GenJetToken_;
  edm::EDGetTokenT< double> rhoToken_;
  
  //

  edm::Service<TFileService> fs;
  TTree *tree;
  
  //JetAnalyzer* theJetAnalyzer;


  //const static int njets_max = 100;

  //float ptCHS1[njets_max],etaCHS1[njets_max],phiCHS1[njets_max];
  //float ptCHS2[njets_max],etaCHS2[njets_max],phiCHS2[njets_max];
  //float ptCHS3[njets_max],etaCHS3[njets_max],phiCHS3[njets_max];
  //float ptCHS4[njets_max],etaCHS4[njets_max],phiCHS4[njets_max];
  //float ptPF[njets_max],etaPF[njets_max],phiPF[njets_max],jecPF[njets_max];
  //float ptGen[njets_max],etaGen[njets_max],phiGen[njets_max];
  //float drPFGen[njets_max],ptPFGen[njets_max],jecPFGen[njets_max];
  //float ptGenCHS1[njets_max],drGenCHS1[njets_max];
  //float ptGenCHS2[njets_max],drGenCHS2[njets_max];
  //float ptGenCHS3[njets_max],drGenCHS3[njets_max];
  //float drGenCHS4[njets_max], ptGenCHS4[njets_max];

  long int nGenJets, nJetsNew0, nJetsNew1, nJetsNew2, nJetsNew3, nJetsNew4, nJetsNew5, nJetsNew6, nJetsNew7, nJets, nPUtrue, EventNumber, LumiNumber, RunNumber, nPV, nPVsel_ndof, nPVsel_pt;
  int nPFCandidates, nPFCandidatesTrack, nPFCandidatesHighPurityTrack;
  //int nPFCandidatesFakeRejected, nPFCandidatesFakeRejectedTrack, nPFCandidatesFakeRejectedHighPurityTrack;
  //std::vector<float> PFCand_dz, PFCand_cosh;
  float  rho;
  //int JetId;
  bool isVerbose, isMC;

  //Structures
  //std::vector<PFCandidateType> PFCands;
  //std::vector<PFCandidateType> PFCandsFakeRejected;
  std::vector<JetType> Jets;
  std::vector<JetType> JetsNew0;
  std::vector<JetType> JetsNew1;
  std::vector<JetType> JetsNew2;
  std::vector<JetType> JetsNew3;
  std::vector<JetType> JetsNew4;
  std::vector<JetType> JetsNew5;
  std::vector<JetType> JetsNew6;
  std::vector<JetType> JetsNew7;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
  //JetPSet(iConfig.getParameter<edm::ParameterSet>("jetSet")),
  jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets"))),  
  jetpt(iConfig.getParameter <double>("jetpt")),  
  pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter <edm::InputTag>("pileup"))),
  PVToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter <edm::InputTag>("vertices"))),
  PFCandToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter <edm::InputTag>("pfcandidates"))),
  //PFCandFakeRejectedToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter <edm::InputTag>("pfcandidatesfakerejected"))),
  jetToken0_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets0"))), 
  jetToken1_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets1"))),
  jetToken2_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets2"))),  
  jetToken3_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets3"))),  
  jetToken4_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets4"))),  
  jetToken5_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets5"))),  
  jetToken6_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets6"))),  
  jetToken7_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets7"))),  
  GenJetToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter <edm::InputTag>("genjets"))),  
  rhoToken_(consumes<double>(iConfig.getParameter <edm::InputTag>("rhosrc"))),
  isVerbose(iConfig.getParameter<bool> ("verbose"))

{

  //theJetAnalyzer      = new JetAnalyzer(JetPSet, consumesCollector());
  //jetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("selectedPatJetsNew1CHS"));//("patJetsAK4PF","","USER"));//
  //newJetToken1_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew1CHS"));//("selectedPatJetsak4PFJetsNew1CHS"));//
  //newJetToken2_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew2CHS"));//("selectedPatJetsak4PFJetsNew2CHS"));//
  //newJetToken3_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew3CHS"));//("selectedPatJetsak4PFJetsNew3CHS"));//
  //newJetToken4_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew4CHS"));//("selectedPatJetsak4PFJetsNew4CHS"));//
  //newJetToken5_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew5CHS"));//("selectedPatJetsak4PFJetsNew4CHS"));//
  //GenJetToken_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"));
  //now do what ever initialization is needed
  usesResource("TFileService");
  //edm::Service<TFileService> fs;

 
}


Ntuplizer::~Ntuplizer()
{
 
  //delete theJetAnalyzer;//waiiiiiiiiiit
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;

  //Initialization
  EventNumber = LumiNumber = RunNumber = nPV = nPVsel_ndof = nPVsel_pt = 0;
  nJets = nJetsNew0 = nJetsNew1 = nJetsNew2 = nJetsNew3 = nJetsNew4 = nJetsNew5 = nJetsNew6 = nJetsNew7 = 0;
  nPFCandidates = nPFCandidatesTrack = nPFCandidatesHighPurityTrack = 0;
  //nPFCandidatesFakeRejected = nPFCandidatesFakeRejectedTrack = nPFCandidatesFakeRejectedHighPurityTrack = 0;
  //JetId = 0;
  isMC = false;

  float PtTh(jetpt);

  //Event info
  isMC = !iEvent.isRealData();
  EventNumber=iEvent.id().event();
  LumiNumber=iEvent.id().luminosityBlock();
  RunNumber=iEvent.id().run();
  
  //Pile up info
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(pileupSummaryToken_,PupInfo);
  nPUtrue = PupInfo -> begin()->getTrueNumInteractions();

  //Alternative way
  edm::Handle<reco::VertexCollection> PVCollection;
  iEvent.getByToken(PVToken_, PVCollection);
  nPV = PVCollection->size();

  ////Count good vertices; thanks to Wolfram Erdmann, but works only over AOD
  //int vertex_counter = 0;
  for(std::vector<reco::Vertex>::const_iterator it=PVCollection->begin(); it!=PVCollection->end(); ++it){
    reco::Vertex v=*it;
    //std::cout << "Vertex n." << vertex_counter << std::endl; 
    //vertex_counter++;
   //option 0: usual ndof cut
    if(select(v, 0)){
      nPVsel_ndof++;
    }
    //option 1: Wolfram's selection
    if(select(v, 1)){
      nPVsel_pt++;
    }
  }

  //MET
  //pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho = *rhoH;
  
  edm::Handle<pat::JetCollection> PFJetsCollection;
  iEvent.getByToken(jetToken_,PFJetsCollection);

  //Fill Jet vector
  //manually; otherwise we always need mutiple psets
  std::vector<pat::Jet> JetsVect;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollection->begin(); it!=PFJetsCollection->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());

    /*
    if(jet.genJet()){
      std::vector<edm::Ptr<reco::Candidate>> GenJetConstituentVect = jet.genJet()->getJetConstituents();
      std::cout << "Gen jet.... " << std::endl;
      std::cout << "had energy: " << jet.genJet()->hadEnergy() << std::endl;
      std::cout << "invisible energy: " << jet.genJet()->invisibleEnergy() << std::endl;
      std::cout << "em energy: " << jet.genJet()->emEnergy() << std::endl;

      float emEnergy = 0;
      float eleEnergy = 0;
      float photonEnergy = 0;
      float muEnergy = 0;
      float cHadEnergy = 0;
      float nHadEnergy = 0;
      int chargedMulti = 0;
      int neutralMulti = 0;
      for(unsigned int k = 0; k < GenJetConstituentVect.size(); k++){
	std::cout << "Jet const n. " << k << std::endl;
	std::cout << "pdgId: " << GenJetConstituentVect[k]->pdgId() << std::endl;
	std::cout << "energy: " << GenJetConstituentVect[k]->energy() << std::endl;
	std::cout << "charge: " << GenJetConstituentVect[k]->charge() << std::endl;
	if(abs(GenJetConstituentVect[k]->pdgId())==22) photonEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())==11) eleEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())==13) muEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())>22 && GenJetConstituentVect[k]->charge()==0) nHadEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())>22 && GenJetConstituentVect[k]->charge()!=0) cHadEnergy += GenJetConstituentVect[k]->energy();
	if(GenJetConstituentVect[k]->charge()==0) neutralMulti++;
	if(GenJetConstituentVect[k]->charge()!=0) chargedMulti++;
      //Methods like isMuon, isElectron, isPhoton do not work. Always getting 0.
      }
      std::cout << "from constituents - ele energy: " << eleEnergy << std::endl;
      std::cout << "from constituents - photon energy: " << photonEnergy << std::endl;
      std::cout << "from constituents - nHad energy: " << nHadEnergy << std::endl;
      std::cout << "from constituents - cHad energy: " << cHadEnergy << std::endl;
      std::cout << "from constituents - total Had energy: " << nHadEnergy+cHadEnergy << std::endl;
      std::cout << "from constituents - EM energy: " << eleEnergy+photonEnergy << std::endl;

      jet.addUserFloat("cHadEGen", cHadEnergy);
      jet.addUserFloat("nHadEGen", nHadEnergy);

      jet.addUserFloat("nMultiGen", float(neutralMulti));
      jet.addUserFloat("cMultiGen", float(chargedMulti));

      jet.addUserFloat("emEGen", eleEnergy+photonEnergy);
      jet.addUserFloat("eleEGen", eleEnergy);
      jet.addUserFloat("photonEGen", photonEnergy);
      
    }
    */

    JetsVect.push_back(jet);
  }
  nJets = JetsVect.size();
  Jets.clear();

  Ntuplizer::GenJetAnalyzer(JetsVect);

  ////**********////
  //// JetsNew0 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew0;
  iEvent.getByToken(jetToken0_,PFJetsCollectionNew0);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew0;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew0->begin(); it!=PFJetsCollectionNew0->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew0.push_back(jet);
  }
  nJetsNew0 = JetsVectNew0.size();
  JetsNew0.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew0);

  ////**********////
  //// JetsNew1 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew1;
  iEvent.getByToken(jetToken1_,PFJetsCollectionNew1);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew1;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew1->begin(); it!=PFJetsCollectionNew1->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew1.push_back(jet);
  }
  nJetsNew1 = JetsVectNew1.size();
  JetsNew1.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew1);

  ////**********////
  //// JetsNew2 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew2;
  iEvent.getByToken(jetToken2_,PFJetsCollectionNew2);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew2;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew2->begin(); it!=PFJetsCollectionNew2->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew2.push_back(jet);
  }
  nJetsNew2 = JetsVectNew2.size();
  JetsNew2.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew2);

  ////**********////
  //// JetsNew3 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew3;
  iEvent.getByToken(jetToken3_,PFJetsCollectionNew3);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew3;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew3->begin(); it!=PFJetsCollectionNew3->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew3.push_back(jet);
  }
  nJetsNew3 = JetsVectNew3.size();
  JetsNew3.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew3);


  ////**********////
  //// JetsNew4 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew4;
  iEvent.getByToken(jetToken4_,PFJetsCollectionNew4);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew4;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew4->begin(); it!=PFJetsCollectionNew4->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew4.push_back(jet);

  }
  nJetsNew4 = JetsVectNew4.size();
  JetsNew4.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew4);

  ////**********////
  //// JetsNew5 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew5;
  iEvent.getByToken(jetToken5_,PFJetsCollectionNew5);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew5;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew5->begin(); it!=PFJetsCollectionNew5->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew5.push_back(jet);

  }
  nJetsNew5 = JetsVectNew5.size();
  JetsNew5.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew5);


  ////**********////
  //// JetsNew6 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew6;
  iEvent.getByToken(jetToken6_,PFJetsCollectionNew6);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew6;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew6->begin(); it!=PFJetsCollectionNew6->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew6.push_back(jet);
  }
  nJetsNew6 = JetsVectNew6.size();
  JetsNew6.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew6);

  ////**********////
  //// JetsNew7 ////
  ////**********////

  edm::Handle<pat::JetCollection> PFJetsCollectionNew7;
  iEvent.getByToken(jetToken7_,PFJetsCollectionNew7);

  //Fill Jet vector
  std::vector<pat::Jet> JetsVectNew7;
  for(std::vector<pat::Jet>::const_iterator it=PFJetsCollectionNew7->begin(); it!=PFJetsCollectionNew7->end(); ++it) {
    pat::Jet jet=*it;
    //JetId not from pset!!!
    //if(JetId==1 && !theJetAnalyzer->isLooseJet(jet)) continue;
    //if(JetId==2 && !theJetAnalyzer->isTightJet(jet)) continue;
    //if(JetId==3 && !theJetAnalyzer->isTightLepVetoJet(jet)) continue;
    jet.addUserInt("isLoose", Ntuplizer::isLooseJet(jet) ? 1 : 0);
    jet.addUserInt("isTight", Ntuplizer::isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVectNew7.push_back(jet);
  }
  nJetsNew7 = JetsVectNew7.size();
  JetsNew7.clear();

  Ntuplizer::GenJetAnalyzer(JetsVectNew7);


  ////**************////
  ////PF candidates////
  ////**************////

  edm::Handle<pat::PackedCandidateCollection> PFCandCollection;
  iEvent.getByToken(PFCandToken_, PFCandCollection);

  //std::vector<pat::PackedCandidate> PFCandsVect;

  for(std::vector<pat::PackedCandidate>::const_iterator it=PFCandCollection->begin(); it!=PFCandCollection->end(); ++it) {
    pat::PackedCandidate pfcand=*it;
    nPFCandidates++;
    //PFCandsVect.push_back(pfcand);
    if(pfcand.charge()!=0) nPFCandidatesTrack++;
    if(pfcand.trackHighPurity()) nPFCandidatesHighPurityTrack++;
  }
  //PFCands.clear();

  //edm::Handle<pat::PackedCandidateCollection> PFCandFakeRejectedCollection;
  //iEvent.getByToken(PFCandFakeRejectedToken_, PFCandFakeRejectedCollection);

  //std::vector<pat::PackedCandidate> PFCandsFakeRejectedVect;

  //for(std::vector<pat::PackedCandidate>::const_iterator it=PFCandFakeRejectedCollection->begin(); it!=PFCandFakeRejectedCollection->end(); ++it) {
  //pat::PackedCandidate pfcandFakeRejected=*it;
  //nPFCandidatesFakeRejected++;
  //PFCandsFakeRejectedVect.push_back(pfcandFakeRejected);
  //if(pfcandFakeRejected.charge()!=0) nPFCandidatesFakeRejectedTrack++;
  //if(pfcandFakeRejected.trackHighPurity()) nPFCandidatesFakeRejectedHighPurityTrack++;
  //}
  //PFCandsFakeRejected.clear();

  edm::Handle<std::vector<reco::GenJet> > GenJetsCollection;
  iEvent.getByToken(GenJetToken_,GenJetsCollection);

  //Fill Jet vector
  //manually; otherwise we always need mutiple psets
  std::vector<reco::GenJet> GenJetsVect;
  for(std::vector<reco::GenJet>::const_iterator it=GenJetsCollection->begin(); it!=GenJetsCollection->end(); ++it) {
    reco::GenJet jet=*it;
    if(jet.pt()<PtTh) continue;
    GenJetsVect.push_back(jet);
  }
  nGenJets = GenJetsVect.size();
  //Jets.clear();


  if(isVerbose) {
    std::cout << " --- Event n. " << iEvent.id().event() << ", lumi " << iEvent.luminosityBlock() << ", run " << iEvent.id().run() << std::endl;
    for(unsigned int i = 0; i < GenJetsVect.size(); i++) std::cout << "  Gen AK4 jet  [" << i << "]\tpt: " << GenJetsVect[i].pt() << "\teta: " << GenJetsVect[i].eta() << "\tphi: " << GenJetsVect[i].phi() << "\tmass: " << GenJetsVect[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVect.size(); i++) std::cout << "  CHS AK4 jet  [" << i << "]\tpt: " << JetsVect[i].pt() << "\teta: " << JetsVect[i].eta() << "\tphi: " << JetsVect[i].phi() << "\tmass: " << JetsVect[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew0.size(); i++) std::cout << "  CHS AK4 New0 jet  [" << i << "]\tpt: " << JetsVectNew0[i].pt() << "\teta: " << JetsVectNew0[i].eta() << "\tphi: " << JetsVectNew0[i].phi() << "\tmass: " << JetsVectNew0[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew1.size(); i++) std::cout << "  CHS AK4 New1 jet  [" << i << "]\tpt: " << JetsVectNew1[i].pt() << "\teta: " << JetsVectNew1[i].eta() << "\tphi: " << JetsVectNew1[i].phi() << "\tmass: " << JetsVectNew1[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew2.size(); i++) std::cout << "  CHS AK4 New2 jet  [" << i << "]\tpt: " << JetsVectNew2[i].pt() << "\teta: " << JetsVectNew2[i].eta() << "\tphi: " << JetsVectNew2[i].phi() << "\tmass: " << JetsVectNew2[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew3.size(); i++) std::cout << "  CHS AK4 New3 jet  [" << i << "]\tpt: " << JetsVectNew3[i].pt() << "\teta: " << JetsVectNew3[i].eta() << "\tphi: " << JetsVectNew3[i].phi() << "\tmass: " << JetsVectNew3[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew4.size(); i++) std::cout << "  CHS AK4 New4 jet  [" << i << "]\tpt: " << JetsVectNew4[i].pt() << "\teta: " << JetsVectNew4[i].eta() << "\tphi: " << JetsVectNew4[i].phi() << "\tmass: " << JetsVectNew4[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew5.size(); i++) std::cout << "  CHS AK4 New5 jet  [" << i << "]\tpt: " << JetsVectNew5[i].pt() << "\teta: " << JetsVectNew5[i].eta() << "\tphi: " << JetsVectNew5[i].phi() << "\tmass: " << JetsVectNew5[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew6.size(); i++) std::cout << "  CHS AK4 New6 jet  [" << i << "]\tpt: " << JetsVectNew6[i].pt() << "\teta: " << JetsVectNew6[i].eta() << "\tphi: " << JetsVectNew6[i].phi() << "\tmass: " << JetsVectNew6[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew7.size(); i++) std::cout << "  CHS AK4 New7 jet  [" << i << "]\tpt: " << JetsVectNew7[i].pt() << "\teta: " << JetsVectNew7[i].eta() << "\tphi: " << JetsVectNew7[i].phi() << "\tmass: " << JetsVectNew7[i].mass() << std::endl;

  }

  //for(unsigned int i = 0; i < PFCandsVect.size(); i++) PFCands.push_back( PFCandidateType() );
  //for(unsigned int i = 0; i < PFCandsFakeRejectedVect.size(); i++) PFCandsFakeRejected.push_back( PFCandidateType() );
  for(unsigned int i = 0; i < JetsVect.size(); i++) Jets.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew0.size(); i++) JetsNew0.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew1.size(); i++) JetsNew1.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew2.size(); i++) JetsNew2.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew3.size(); i++) JetsNew3.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew4.size(); i++) JetsNew4.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew5.size(); i++) JetsNew5.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew6.size(); i++) JetsNew6.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew7.size(); i++) JetsNew7.push_back( JetType() );

  //for(unsigned int i = 0; i < PFCandsVect.size(); i++) ObjectsFormat::FillPFCandidateType(PFCands[i], &PFCandsVect[i], isMC);
  //for(unsigned int i = 0; i < PFCandsFakeRejectedVect.size(); i++) ObjectsFormat::FillPFCandidateType(PFCandsFakeRejected[i], &PFCandsFakeRejectedVect[i], isMC);
  for(unsigned int i = 0; i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew0.size(); i++) ObjectsFormat::FillJetType(JetsNew0[i], &JetsVectNew0[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew1.size(); i++) ObjectsFormat::FillJetType(JetsNew1[i], &JetsVectNew1[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew2.size(); i++) ObjectsFormat::FillJetType(JetsNew2[i], &JetsVectNew2[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew3.size(); i++) ObjectsFormat::FillJetType(JetsNew3[i], &JetsVectNew3[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew4.size(); i++) ObjectsFormat::FillJetType(JetsNew4[i], &JetsVectNew4[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew5.size(); i++) ObjectsFormat::FillJetType(JetsNew5[i], &JetsVectNew5[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew6.size(); i++) ObjectsFormat::FillJetType(JetsNew6[i], &JetsVectNew6[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew7.size(); i++) ObjectsFormat::FillJetType(JetsNew7[i], &JetsVectNew7[i], isMC);





  tree->Fill();
}


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif



// ------------ method called once each job just before starting event loop  ------------
void 
Ntuplizer::beginJob()
{
  //Tree Branch
  tree=fs->make<TTree>("tree","tree");
  tree->Branch("isMC" , &isMC, "isMC/O");
  tree->Branch("nJets",&nJets,"nJets/I");
  tree->Branch("nJetsNew0",&nJetsNew0,"nJetsNew0/I");
  tree->Branch("nJetsNew1",&nJetsNew1,"nJetsNew1/I");
  tree->Branch("nJetsNew2",&nJetsNew2,"nJetsNew2/I");
  tree->Branch("nJetsNew3",&nJetsNew3,"nJetsNew3/I");
  tree->Branch("nJetsNew4",&nJetsNew4,"nJetsNew4/I");
  tree->Branch("nJetsNew5",&nJetsNew5,"nJetsNew5/I");
  tree->Branch("nJetsNew6",&nJetsNew6,"nJetsNew6/I");
  tree->Branch("nJetsNew7",&nJetsNew7,"nJetsNew7/I");
  tree->Branch("nGenJets",&nGenJets,"nGenJets/I");
  
  tree-> Branch("nPUtrue",&nPUtrue,"nPUtrue/I");
  tree-> Branch("nPV",&nPV,"nPV/I");
  tree-> Branch("nPVsel_ndof",&nPVsel_ndof,"nPVsel_ndof/I");
  tree-> Branch("nPVsel_pt",&nPVsel_pt,"nPVsel_pt/I");
  tree-> Branch("RunNumber",&RunNumber,"RunNumber/I");
  tree-> Branch("LumiNumber",&LumiNumber,"LumiNumber/I");
  tree-> Branch("EventNumber",&EventNumber,"EventNumber/I");
  tree-> Branch("rho",&rho,"rho/F");

  tree -> Branch("Jets", &Jets);
  tree -> Branch("JetsNew0", &JetsNew0);
  tree -> Branch("JetsNew1", &JetsNew1);
  tree -> Branch("JetsNew2", &JetsNew2);
  tree -> Branch("JetsNew3", &JetsNew3);
  tree -> Branch("JetsNew4", &JetsNew4);
  tree -> Branch("JetsNew5", &JetsNew5);
  tree -> Branch("JetsNew6", &JetsNew6);
  tree -> Branch("JetsNew7", &JetsNew7);
  tree -> Branch("nPFCandidates" , &nPFCandidates, "nPFCandidates/I");
  tree -> Branch("nPFCandidatesTrack", &nPFCandidatesTrack, "nPFCandidatesTrack/I");
  tree -> Branch("nPFCandidatesHighPurityTrack", &nPFCandidatesHighPurityTrack, "nPFCandidatesHighPurityTrack/I");
  //tree -> Branch("PFCands", &PFCands);

  //tree -> Branch("nPFCandidatesFakeRejected" , &nPFCandidatesFakeRejected, "nPFCandidatesFakeRejected/I");
  //tree -> Branch("nPFCandidatesFakeRejectedTrack", &nPFCandidatesFakeRejectedTrack, "nPFCandidatesFakeRejectedTrack/I");
  //tree -> Branch("nPFCandidatesFakeRejectedHighPurityTrack", &nPFCandidatesFakeRejectedHighPurityTrack, "nPFCandidatesFakeRejectedHighPurityTrack/I");
  //tree -> Branch("PFCandsFakeRejected", &PFCandsFakeRejected);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntuplizer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  }*/

void Ntuplizer::GenJetAnalyzer(std::vector<pat::Jet>& Jets)
{

  for(unsigned int i=0; i < Jets.size(); i++) {

    //std::cout << "Jet n. " << i <<std::endl;
    
    if(Jets.at(i).genJet()){
      std::vector<edm::Ptr<reco::Candidate>> GenJetConstituentVect = Jets.at(i).genJet()->getJetConstituents();
      //std::cout << "Gen Jet: " << std::endl;
      //std::cout << "had energy: " << Jets.at(i).genJet()->hadEnergy() << std::endl;
      //std::cout << "invisible energy: " << Jets.at(i).genJet()->invisibleEnergy() << std::endl;
      //std::cout << "em energy: " << Jets.at(i).genJet()->emEnergy() << std::endl;

      float emEnergy = 0;
      float eleEnergy = 0;
      float photonEnergy = 0;
      float muEnergy = 0;
      float cHadEnergy = 0;
      float nHadEnergy = 0;
      int chargedMulti = 0;
      int neutralMulti = 0;
      int multiplicity = 0;
      for(unsigned int k = 0; k < GenJetConstituentVect.size(); k++){
	//std::cout << "Jet const n. " << k << std::endl;
	//std::cout << "pdgId: " << GenJetConstituentVect[k]->pdgId() << std::endl;
	//std::cout << "energy: " << GenJetConstituentVect[k]->energy() << std::endl;
	//std::cout << "charge: " << GenJetConstituentVect[k]->charge() << std::endl;
	if(abs(GenJetConstituentVect[k]->pdgId())==22) photonEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())==11) eleEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())==13) muEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())>22 && GenJetConstituentVect[k]->charge()==0) nHadEnergy += GenJetConstituentVect[k]->energy();
	if(abs(GenJetConstituentVect[k]->pdgId())>22 && GenJetConstituentVect[k]->charge()!=0) cHadEnergy += GenJetConstituentVect[k]->energy();
	if(GenJetConstituentVect[k]->charge()==0) neutralMulti++;
	if(GenJetConstituentVect[k]->charge()!=0) chargedMulti++;
      //Methods like isMuon, isElectron, isPhoton do not work. Always getting 0.
      }
      multiplicity =  GenJetConstituentVect.size();
      //std::cout << "from constituents - ele energy: " << eleEnergy << std::endl;
      //std::cout << "from constituents - photon energy: " << photonEnergy << std::endl;
      //std::cout << "from constituents - nHad energy: " << nHadEnergy << std::endl;
      //std::cout << "from constituents - cHad energy: " << cHadEnergy << std::endl;
      //std::cout << "from constituents - total Had energy: " << nHadEnergy+cHadEnergy << std::endl;
      //std::cout << "from constituents - EM energy: " << eleEnergy+photonEnergy << std::endl;

      //Energies
      Jets.at(i).addUserFloat("cHadEGen", cHadEnergy);
      Jets.at(i).addUserFloat("nHadEGen", nHadEnergy);

      Jets.at(i).addUserFloat("emEGen", eleEnergy+photonEnergy);
      Jets.at(i).addUserFloat("eleEGen", eleEnergy);
      Jets.at(i).addUserFloat("photonEGen", photonEnergy);


      //Multiplicities
      Jets.at(i).addUserFloat("nMultiGen", float(neutralMulti));
      Jets.at(i).addUserFloat("cMultiGen", float(chargedMulti));

      Jets.at(i).addUserFloat("multiGen", float(multiplicity));

    }    

  }

}


// PFJet Quality ID 2015-2016: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
bool Ntuplizer::isLooseJet(pat::Jet& jet)
{
  if(fabs(jet.eta())<=2.7){ /// |eta| < 2.7
    if(jet.neutralHadronEnergyFraction()>=0.99) return false;
    if(jet.neutralEmEnergyFraction()>=0.99) return false;
    if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
    if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
      if(jet.chargedHadronEnergyFraction()<=0.) return false;
      if(jet.chargedMultiplicity()<=0) return false;
      if(jet.chargedEmEnergyFraction()>=0.99) return false;
    }
  }
  else{ /// |eta| > 2.7
    if(jet.neutralEmEnergyFraction()>=0.90) return false;
    if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0
      if(jet.neutralMultiplicity()<=2) return false;
    }
    else{ /// |eta| > 3.0                                                                                                                                                       
      if(jet.neutralMultiplicity()<=10) return false;
    }
  }
  return true;
}

bool Ntuplizer::isTightJet(pat::Jet& jet) {
  if(fabs(jet.eta())<=2.7){ /// |eta| < 2.7
    if(jet.neutralHadronEnergyFraction()>=0.90) return false;
    if(jet.neutralEmEnergyFraction()>=0.90) return false;
    if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
    if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
      if(jet.chargedHadronEnergyFraction()<=0.) return false;
      if(jet.chargedMultiplicity()<=0) return false;
      if(jet.chargedEmEnergyFraction()>=0.99) return false;
    }
  }
  else{ /// |eta| > 2.7
    if(jet.neutralEmEnergyFraction()>=0.90) return false;
    if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0                                                                                                                           
      if(jet.neutralMultiplicity()<=2) return false;
    }
    else{ /// |eta| > 3.0
      if(jet.neutralMultiplicity()<=10) return false;
    }
  }
  return true;
}

//Thanks to Wolfram Erdmann!
double Ntuplizer::vertex_ptmax2(const reco::Vertex & v) {
  double ptmax1 = 0;
  double ptmax2 = 0;
  //std::cout << "vertex_ptmax2 in action!" << std::endl;
  //std::cout << "vertex had ntracks: " << v.nTracks(0.5) << std::endl;
  //std::cout << "vertex has n dof: " << v.ndof() << std::endl;
  //std::cout << "vertex has track size: " << v.tracksSize() << std::endl;
  for(reco::Vertex::trackRef_iterator t = v.tracks_begin(); t != v.tracks_end(); t++){
    //std::cout << "in looop" << std::endl;
    //std::cout << v.trackWeight(*t) << std::endl;
    if (v.trackWeight(*t) > 0.5){
      double pt = t->get()->pt();
      //std::cout << "pt: " << pt << std::endl;
      if (pt > ptmax1){
        ptmax2 = ptmax1;
        ptmax1 = pt;
      }else if(pt > ptmax2){
        ptmax2 = pt;
      }
    }
  }
  return ptmax2; 
}


bool Ntuplizer::select(const reco::Vertex & v, int level){
  /* level
     0  !isFake  && ndof>4
     1  !isFake  && ndof>4 &&  ptmax2 >0.4
  */
  //std::cout << Ntuplizer::vertex_ptmax2(v) << std::endl;
  if( v.isFake() ) return false;
  if( (level == 0) && (v.ndof()>4) ) return true;
  if( (level == 1) && (v.ndof()>4) && (Ntuplizer::vertex_ptmax2(v)>0.4) ) return true;
  return false;
} 

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
