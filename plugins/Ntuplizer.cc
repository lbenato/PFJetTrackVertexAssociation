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
  edm::EDGetTokenT<pat::JetCollection> jetToken1_;
  edm::EDGetTokenT<pat::JetCollection> jetToken2_;
  edm::EDGetTokenT<pat::JetCollection> jetToken3_;
  edm::EDGetTokenT<pat::JetCollection> jetToken4_;
  edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken1_;
  edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken2_;
  edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken3_;
  edm::EDGetTokenT< std::vector<reco::PFJet> > newJetToken4_;
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

  long int nGenJets, nJetsNew1, nJetsNew2, nJetsNew3, nJetsNew4, nJets, nPUtrue, EventNumber, LumiNumber, RunNumber, nPV;
  float  rho;
  //int JetId;
  bool isVerbose, isMC;

  //Structures
  std::vector<JetType> Jets;
  std::vector<JetType> JetsNew1;
  std::vector<JetType> JetsNew2;
  std::vector<JetType> JetsNew3;
  std::vector<JetType> JetsNew4;

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
  jetToken1_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets1"))),  
  jetToken2_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets2"))),  
  jetToken3_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets3"))),  
  jetToken4_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets4"))),  
  GenJetToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter <edm::InputTag>("genjets"))),  
  rhoToken_(consumes<double>(iConfig.getParameter <edm::InputTag>("rhosrc"))),
  isVerbose(iConfig.getParameter<bool> ("verbose"))

{

  //theJetAnalyzer      = new JetAnalyzer(JetPSet, consumesCollector());
  //jetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("selectedPatJetsNew1CHS"));//("patJetsAK4PF","","USER"));//
  newJetToken1_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew1CHS"));//("selectedPatJetsak4PFJetsNew1CHS"));//
  newJetToken2_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew2CHS"));//("selectedPatJetsak4PFJetsNew2CHS"));//
  newJetToken3_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew3CHS"));//("selectedPatJetsak4PFJetsNew3CHS"));//
  newJetToken4_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew4CHS"));//("selectedPatJetsak4PFJetsNew4CHS"));//
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
  EventNumber = LumiNumber = RunNumber = nPV = 0;
  nJets = nJetsNew1 = nJetsNew2 = nJetsNew3 = nJetsNew4 = 0;
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
    //jet.addUserInt("isLoose", theJetAnalyzer->isLooseJet(jet) ? 1 : 0);
    //jet.addUserInt("isTight", theJetAnalyzer->isTightJet(jet) ? 1 : 0);
    //jet.addUserInt("isTightLepVeto", theJetAnalyzer->isTightLepVetoJet(jet) ? 1 : 0);
    jet.addUserFloat("cHadEFrac", jet.chargedHadronEnergyFraction());
    jet.addUserFloat("nHadEFrac", jet.neutralHadronEnergyFraction());
    jet.addUserFloat("nEmEFrac", jet.neutralEmEnergyFraction());
    jet.addUserFloat("cEmEFrac", jet.chargedEmEnergyFraction());
    jet.addUserFloat("cmuEFrac", jet.chargedMuEnergyFraction());
    jet.addUserFloat("muEFrac", jet.muonEnergyFraction());
    jet.addUserFloat("eleEFrac", jet.electronEnergyFraction());
    jet.addUserFloat("photonEFrac", jet.photonEnergyFraction());
    JetsVect.push_back(jet);
  }
  nJets = JetsVect.size();
  Jets.clear();


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
    //jet.addUserInt("isLoose", theJetAnalyzer->isLooseJet(jet) ? 1 : 0);
    //jet.addUserInt("isTight", theJetAnalyzer->isTightJet(jet) ? 1 : 0);
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
    //jet.addUserInt("isLoose", theJetAnalyzer->isLooseJet(jet) ? 1 : 0);
    //jet.addUserInt("isTight", theJetAnalyzer->isTightJet(jet) ? 1 : 0);
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
    //jet.addUserInt("isLoose", theJetAnalyzer->isLooseJet(jet) ? 1 : 0);
    //jet.addUserInt("isTight", theJetAnalyzer->isTightJet(jet) ? 1 : 0);
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
    //jet.addUserInt("isLoose", theJetAnalyzer->isLooseJet(jet) ? 1 : 0);
    //jet.addUserInt("isTight", theJetAnalyzer->isTightJet(jet) ? 1 : 0);
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
    for(unsigned int i = 0; i < JetsVectNew1.size(); i++) std::cout << "  CHS AK4 New1 jet  [" << i << "]\tpt: " << JetsVectNew1[i].pt() << "\teta: " << JetsVectNew1[i].eta() << "\tphi: " << JetsVectNew1[i].phi() << "\tmass: " << JetsVectNew1[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew2.size(); i++) std::cout << "  CHS AK4 New2 jet  [" << i << "]\tpt: " << JetsVectNew2[i].pt() << "\teta: " << JetsVectNew2[i].eta() << "\tphi: " << JetsVectNew2[i].phi() << "\tmass: " << JetsVectNew2[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew3.size(); i++) std::cout << "  CHS AK4 New3 jet  [" << i << "]\tpt: " << JetsVectNew3[i].pt() << "\teta: " << JetsVectNew3[i].eta() << "\tphi: " << JetsVectNew3[i].phi() << "\tmass: " << JetsVectNew3[i].mass() << std::endl;
    for(unsigned int i = 0; i < JetsVectNew4.size(); i++) std::cout << "  CHS AK4 New4 jet  [" << i << "]\tpt: " << JetsVectNew4[i].pt() << "\teta: " << JetsVectNew4[i].eta() << "\tphi: " << JetsVectNew4[i].phi() << "\tmass: " << JetsVectNew4[i].mass() << std::endl;
  }

  for(unsigned int i = 0; i < JetsVect.size(); i++) Jets.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew1.size(); i++) JetsNew1.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew2.size(); i++) JetsNew2.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew3.size(); i++) JetsNew3.push_back( JetType() );
  for(unsigned int i = 0; i < JetsVectNew4.size(); i++) JetsNew4.push_back( JetType() );

  for(unsigned int i = 0; i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew1.size(); i++) ObjectsFormat::FillJetType(JetsNew1[i], &JetsVectNew1[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew2.size(); i++) ObjectsFormat::FillJetType(JetsNew2[i], &JetsVectNew2[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew3.size(); i++) ObjectsFormat::FillJetType(JetsNew3[i], &JetsVectNew3[i], isMC);
  for(unsigned int i = 0; i < JetsVectNew4.size(); i++) ObjectsFormat::FillJetType(JetsNew4[i], &JetsVectNew4[i], isMC);






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
  //tree->Branch("ptPF",ptPF,"ptPF[nJets]/F");
  //tree->Branch("etaPF",etaPF,"etaPF[nJets]/F");
  //tree->Branch("phiPF",phiPF,"phiPF[nJets]/F");
  //tree->Branch("jecPF",jecPF,"jecPF[nJets]/F");
  tree->Branch("nJetsNew1",&nJetsNew1,"nJetsNew1/I");
  tree->Branch("nJetsNew2",&nJetsNew2,"nJetsNew2/I");
  tree->Branch("nJetsNew3",&nJetsNew3,"nJetsNew3/I");
  tree->Branch("nJetsNew4",&nJetsNew4,"nJetsNew4/I");
 
  //tree->Branch("ptCHS1",ptCHS1,"ptCHS1[nJetsNew1]/F");
  //tree->Branch("etaCHS1",etaCHS1,"etaCHS1[nJetsNew1]/F");
  //tree->Branch("phiCHS1",phiCHS1,"phiCHS1[nJetsNew1]/F");
  
  //tree->Branch("ptCHS2",ptCHS2,"ptCHS2[nJetsNew2]/F");
  //tree->Branch("etaCHS2",etaCHS2,"etaCHS2[nJetsNew2]/F");
  //tree->Branch("phiCHS2",phiCHS2,"phiCHS2[nJetsNew2]/F");
  
  //tree->Branch("ptCHS3",ptCHS3,"ptCHS3[nJetsNew3]/F");
  //tree->Branch("etaCHS3",etaCHS3,"etaCHS3[nJetsNew3]/F");
  //tree->Branch("phiCHS3",phiCHS3,"phiCHS3[nJetsNew3]/F");
  
  //tree->Branch("ptCHS4",ptCHS4,"ptCHS4[nJetsNew4]/F");
  //tree->Branch("etaCHS4",etaCHS4,"etaCHS4[nJetsNew4]/F");
  //tree->Branch("phiCHS4",phiCHS4,"phiCHS4[nJetsNew4]/F");

 
  tree->Branch("nGenJets",&nGenJets,"nGenJets/I");
  //tree->Branch("ptGen",ptGen,"ptGen[nGenJets]/F");
  //tree->Branch("ptGen",ptGen,"ptGen[nGenJets]/F");
  //tree->Branch("etaGen",etaGen,"etaGen[nGenJets]/F");
  //tree->Branch("phiGen",phiGen,"phiGen[nGenJets]/F");
  //tree->Branch("ptPFGen",ptPFGen,"ptPFGen[nGenJets]/F");
  //tree->Branch("jecPFGen",jecPFGen,"jecPFGen[nGenJets]/F");
  //tree->Branch("drPFGen",drPFGen,"drPFGen[nGenJets]/F");
  //tree->Branch("ptGenCHS3",ptGenCHS3,"ptGenCHS3[nGenJets]/F");
  //tree->Branch("drGenCHS3",drGenCHS3,"drGenCHS3[nGenJets]/F");

  //tree->Branch("ptGenCHS1",ptGenCHS1,"ptGenCHS1[nGenJets]/F");
  //tree->Branch("drGenCHS1",drGenCHS1,"drGenCHS1[nGenJets]/F");

  //tree->Branch("ptGenCHS2",ptGenCHS2,"ptGenCHS2[nGenJets]/F");
  //tree->Branch("drGenCHS2",drGenCHS2,"drGenCHS2[nGenJets]/F");
 
  //tree->Branch("ptGenCHS4",ptGenCHS4,"ptGenCHS4[nGenJets]/F");
  //tree->Branch("drGenCHS4",drGenCHS4,"drGenCHS4[nGenJets]/F");
  
  tree-> Branch("nPUtrue",&nPUtrue,"nPUtrue/I");
  tree-> Branch("nPV",&nPV,"nPV/I");
  tree-> Branch("RunNumber",&RunNumber,"RunNumber/I");
  tree-> Branch("LumiNumber",&LumiNumber,"LumiNumber/I");
  tree-> Branch("EventNumber",&EventNumber,"EventNumber/I");
  tree-> Branch("rho",&rho,"rho/F");

  tree -> Branch("Jets", &Jets);
  tree -> Branch("JetsNew1", &JetsNew1);
  tree -> Branch("JetsNew2", &JetsNew2);
  tree -> Branch("JetsNew3", &JetsNew3);
  tree -> Branch("JetsNew4", &JetsNew4);

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

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
