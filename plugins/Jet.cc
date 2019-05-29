// -*- C++ -*-
//
// Package:    Analyzer/PFJetTrackVertexAssociation
// Class:      PFJetTrackVertexAssociation
// 
/**\class PFJetTrackVertexAssociation PFJetTrackVertexAssociation.cc Analyzer/PFJetTrackVertexAssociation/plugins/PFJetTrackVertexAssociation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Negin Shafiei, Lisa Benato
//         Created:  Wed, 29 May 2019 16:25:24 GMT
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
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Jet : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Jet(const edm::ParameterSet&);
      ~Jet();

  // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------


  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupSummaryToken_;
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

  TTree *event;
  const static int njets_max = 100;

  float ptCHS1[njets_max],etaCHS1[njets_max],phiCHS1[njets_max];
  float ptCHS2[njets_max],etaCHS2[njets_max],phiCHS2[njets_max];
  float ptCHS3[njets_max],etaCHS3[njets_max],phiCHS3[njets_max];
  float ptCHS4[njets_max],etaCHS4[njets_max],phiCHS4[njets_max];
  float ptPF[njets_max],etaPF[njets_max],phiPF[njets_max],jecPF[njets_max];
  float ptGen[njets_max],etaGen[njets_max],phiGen[njets_max];
  float drPFGen[njets_max],ptPFGen[njets_max],jecPFGen[njets_max];
  float ptGenCHS1[njets_max],drGenCHS1[njets_max];
  float ptGenCHS2[njets_max],drGenCHS2[njets_max];
  float ptGenCHS3[njets_max],drGenCHS3[njets_max];
  float drGenCHS4[njets_max], ptGenCHS4[njets_max];

  long int nGenJets,nNewJets1,nNewJets2,nNewJets3,nNewJets4,nPFJets, nPUtrue, nevt,run;
  float  rho ;
  bool isVerbose;
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
Jet::Jet(const edm::ParameterSet& iConfig):
  jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets"))),  
  pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter <edm::InputTag>("pileup"))),
  jetToken1_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets1"))),  
  jetToken2_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets2"))),  
  jetToken3_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets3"))),  
  jetToken4_(consumes<std::vector<pat::Jet>>(iConfig.getParameter <edm::InputTag>("jets4"))),  
  GenJetToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter <edm::InputTag>("genjets"))),  
  rhoToken_(consumes<double>(iConfig.getParameter <edm::InputTag>("rhosrc"))),
  isVerbose(iConfig.getParameter<bool> ("verbose"))

{
  //jetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("selectedPatJetsNew1CHS"));//("patJetsAK4PF","","USER"));//
 newJetToken1_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew1CHS"));//("selectedPatJetsak4PFJetsNew1CHS"));//
 newJetToken2_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew2CHS"));//("selectedPatJetsak4PFJetsNew2CHS"));//
 newJetToken3_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew3CHS"));//("selectedPatJetsak4PFJetsNew3CHS"));//
 newJetToken4_ = consumes<std::vector<reco::PFJet> >(edm::InputTag("ak4PFJetsNew4CHS"));//("selectedPatJetsak4PFJetsNew4CHS"));//
 //GenJetToken_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"));
   //now do what ever initialization is needed
 usesResource("TFileService");
 edm::Service<TFileService> js;


 event=js->make<TTree>("event","event");
 // add  rho, event and run number 
 event->Branch("nPFJets",&nPFJets,"nPFJets/I");
 event->Branch("ptPF",ptPF,"ptPF[nPFJets]/F");
 event->Branch("etaPF",etaPF,"etaPF[nPFJets]/F");
 event->Branch("phiPF",phiPF,"phiPF[nPFJets]/F");
 event->Branch("jecPF",jecPF,"jecPF[nPFJets]/F");
 event->Branch("nNewJets1",&nNewJets1,"nNewJets1/I");
 event->Branch("nNewJets2",&nNewJets2,"nNewJets2/I");
 event->Branch("nNewJets3",&nNewJets3,"nNewJets3/I");
 event->Branch("nNewJets4",&nNewJets4,"nNewJets4/I");
 
 event->Branch("ptCHS1",ptCHS1,"ptCHS1[nNewJets1]/F");
 event->Branch("etaCHS1",etaCHS1,"etaCHS1[nNewJets1]/F");
 event->Branch("phiCHS1",phiCHS1,"phiCHS1[nNewJets1]/F");

 event->Branch("ptCHS2",ptCHS2,"ptCHS2[nNewJets2]/F");
 event->Branch("etaCHS2",etaCHS2,"etaCHS2[nNewJets2]/F");
 event->Branch("phiCHS2",phiCHS2,"phiCHS2[nNewJets2]/F");

 event->Branch("ptCHS3",ptCHS3,"ptCHS3[nNewJets3]/F");
 event->Branch("etaCHS3",etaCHS3,"etaCHS3[nNewJets3]/F");
 event->Branch("phiCHS3",phiCHS3,"phiCHS3[nNewJets3]/F");

 event->Branch("ptCHS4",ptCHS4,"ptCHS4[nNewJets4]/F");
 event->Branch("etaCHS4",etaCHS4,"etaCHS4[nNewJets4]/F");
 event->Branch("phiCHS4",phiCHS4,"phiCHS4[nNewJets4]/F");

 
 event->Branch("nGenJets",&nGenJets,"nGenJets/I");
 event->Branch("ptGen",ptGen,"ptGen[nGenJets]/F");
 event->Branch("ptGen",ptGen,"ptGen[nGenJets]/F");
 event->Branch("etaGen",etaGen,"etaGen[nGenJets]/F");
 event->Branch("phiGen",phiGen,"phiGen[nGenJets]/F");
 event->Branch("ptPFGen",ptPFGen,"ptPFGen[nGenJets]/F");
 event->Branch("jecPFGen",jecPFGen,"jecPFGen[nGenJets]/F");
 event->Branch("drPFGen",drPFGen,"drPFGen[nGenJets]/F");
 event->Branch("ptGenCHS3",ptGenCHS3,"ptGenCHS3[nGenJets]/F");
 event->Branch("drGenCHS3",drGenCHS3,"drGenCHS3[nGenJets]/F");

 event->Branch("ptGenCHS1",ptGenCHS1,"ptGenCHS1[nGenJets]/F");
 event->Branch("drGenCHS1",drGenCHS1,"drGenCHS1[nGenJets]/F");

 event->Branch("ptGenCHS2",ptGenCHS2,"ptGenCHS2[nGenJets]/F");
 event->Branch("drGenCHS2",drGenCHS2,"drGenCHS2[nGenJets]/F");
 
 event->Branch("ptGenCHS4",ptGenCHS4,"ptGenCHS4[nGenJets]/F");
 event->Branch("drGenCHS4",drGenCHS4,"drGenCHS4[nGenJets]/F");
 
 event-> Branch("nPUtrue",&nPUtrue,"nPUtrue/I");
 event-> Branch("nevt",&nevt,"nevt/I");
 event-> Branch("run",&run,"run/I");
 event-> Branch("rho",&rho,"rho/F");
 

}


Jet::~Jet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Jet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  using namespace std;
  
  using namespace reco;
  using namespace pat;
  
  nevt=iEvent.id().event();
  run=iEvent.id().run();
  
  //std::cout<<nevt<<"number of event"<<std::endl;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(pileupSummaryToken_,PupInfo);
  nPUtrue = PupInfo -> begin()->getTrueNumInteractions();

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho = *rhoH;
  
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_,jets);
  nPFJets=0;
  for (unsigned int i = 0, g=jets->size(); i < g && i<njets_max;++i){
    const pat::Jet &pfjet = (*jets)[i];	 
    ptPF[i]=pfjet.pt();
    etaPF[i]=pfjet.eta();
    phiPF[i]=pfjet.phi();
    jecPF[i]=pfjet.jecFactor("Uncorrected");
    if(isVerbose) printf("PFJet %d pt=%6.1f eta=%5.2f phi=%5.2f nConstituents =%d = %d + %d \n"
    	   , i , ptPF[i], etaPF[i], phiPF[i], pfjet.nConstituents(),pfjet.chargedHadronMultiplicity(),pfjet.neutralHadronMultiplicity());
    nPFJets++;
  }


  
  //moving to pat jets
  edm::Handle<pat::JetCollection> ak4NCHS1;
  iEvent.getByToken(jetToken1_,ak4NCHS1);
  
  nNewJets1 =0;
  //std::cout << "New jet size: " << ak4NCHS1->size() << std::endl;;
  for (unsigned int a = 0, k=ak4NCHS1->size(); a < k && a<njets_max ;++a) {
    const pat::Jet & Sak4NCHS1 = (*ak4NCHS1)[a];
    ptCHS1[a]=Sak4NCHS1.pt();
    etaCHS1[a]=Sak4NCHS1.eta();
    phiCHS1[a]=Sak4NCHS1.phi();
    if(isVerbose) printf("NewJet1 %d pt=%6.1f eta=%5.2f phi=%5.2f nConstituents =%d = %d + %d \n"
  	   , a , ptCHS1[a], etaCHS1[a], phiCHS1[a], Sak4NCHS1.nConstituents(),Sak4NCHS1.chargedHadronMultiplicity(),Sak4NCHS1.neutralHadronMultiplicity());
    nNewJets1++;    
  }

  

  edm::Handle<pat::JetCollection> ak4NCHS2;
  iEvent.getByToken(jetToken2_,ak4NCHS2);
  
  nNewJets2 =0;
  for (unsigned int a = 0, k=ak4NCHS2->size(); a < k && a<njets_max ;++a) {
    const pat::Jet & Sak4NCHS2 = (*ak4NCHS2)[a];
    ptCHS2[a]=Sak4NCHS2.pt();
    etaCHS2[a]=Sak4NCHS2.eta();
    phiCHS2[a]=Sak4NCHS2.phi();
    if(isVerbose) printf("NewJets2 %d pt=%6.1f eta=%5.2f phi=%5.2f nConstituents =%d = %d + %d \n"
	   , a ,ptCHS2[a],  etaCHS2[a], phiCHS2[a], Sak4NCHS2.nConstituents(),Sak4NCHS2.chargedHadronMultiplicity(),Sak4NCHS2.neutralHadronMultiplicity());
    nNewJets2++;    
  }

 
  //..................
  edm::Handle<pat::JetCollection> ak4NCHS3;
  iEvent.getByToken(jetToken3_,ak4NCHS3);

  nNewJets3 =0;
  for (unsigned int a = 0, k=ak4NCHS3->size(); a < k && a<njets_max ;++a) {
    const pat::Jet & Sak4NCHS3 = (*ak4NCHS3)[a];
    ptCHS3[a]=Sak4NCHS3.pt();
    etaCHS3[a]=Sak4NCHS3.eta();
    phiCHS3[a]=Sak4NCHS3.phi();
    if(isVerbose) printf("NewJets3 %d pt=%6.1f eta=%5.2f phi=%5.2f nConstituents =%d = %d + %d \n"
	   , a ,ptCHS3[a],etaCHS3[a],phiCHS3[a], Sak4NCHS3.nConstituents(),Sak4NCHS3.chargedHadronMultiplicity(),Sak4NCHS3.neutralHadronMultiplicity());
    nNewJets3++;
  }


  edm::Handle<pat::JetCollection> ak4NCHS4;
  iEvent.getByToken(jetToken4_,ak4NCHS4);
  
  nNewJets4 =0;
   
  for (unsigned int a = 0, k=ak4NCHS4->size(); a < k && a<njets_max ;++a) {
    const pat::Jet & Sak4NCHS4 = (*ak4NCHS4)[a];
     ptCHS4[a]=Sak4NCHS4.pt();
     etaCHS4[a]=Sak4NCHS4.eta();
     phiCHS4[a]=Sak4NCHS4.phi();
    
     if(isVerbose) printf("NewJets4 %d pt=%3.1f eta=%5.2f phi=%5.2f nConstituents =%d = %d + %d \n"
    	    , a , Sak4NCHS4.pt(), Sak4NCHS4.eta(), Sak4NCHS4.phi(), Sak4NCHS4.nConstituents(),Sak4NCHS4.chargedHadronMultiplicity(),Sak4NCHS4.neutralHadronMultiplicity());
 
     nNewJets4++;
  }
 
  edm::Handle<std::vector<reco::GenJet> > Genj;
  iEvent.getByToken(GenJetToken_,Genj);
  
  /*
  nGenJets =0;
  for (unsigned int j = 0, gg=Genj->size(); j < gg && j<njets_max; ++j){ 
    const reco::GenJet &Gj = (*Genj)[j];
    //    printf("GenJet %d pt=%6.1f eta=%5.2f phi=%5.2f \n"
    //	   , j , Gj.pt(), Gj.eta(), Gj.phi());

    if(Gj.pt()>15){
    ptGen[j]=Gj.pt();
    etaGen[j]=Gj.eta();
    phiGen[j]=Gj.phi();
    nGenJets++;
    //L//printf("GenJet before %d pt=%6.1f eta=%5.2f phi=%5.2f \n"
    //L//	   , j , Gj.pt(), Gj.eta(), Gj.phi());
    //........................
    ptGenCHS1[j]= -1;
    drGenCHS1[j]= 100;
    //   int iMinDRNew1 = -1;
    for (unsigned int a = 0, k=ak4NCHS1->size(); a < k;++a) {
      float  deltaR = reco::deltaR((*Genj)[j].p4(), (*ak4NCHS1)[a].p4());
      //  if((*ak4NCHS1)[a].pt()<15){
     
      if(deltaR< drGenCHS1[j]){
	drGenCHS1[j] = deltaR;
	ptGenCHS1[j]=(*ak4NCHS1)[a].pt();
	//iMinDRNew1 = a;
      }
      // }
      //L//printf("OLd cut %d pt=%6.1f eta=%5.2f phi=%5.2f \n"
      //L//   , a , (*ak4NCHS1)[a].pt(), (*ak4NCHS1)[a].eta(),(*ak4NCHS1)[a].phi());
    }
   
  
    //................................................
    ptGenCHS2[j]= -1;
    drGenCHS2[j]= 100;
    vector<TLorentzVector>  fromPV1, Genjet1;
    int iMinDRNew2 = -1;
    vector<float>  Genjet1_dz, fromPV1_dz;
    vector<float>  Genjet1_dzE, fromPV1_dzE;
    for (unsigned int a = 0, k=ak4NCHS2->size(); a < k;++a) {

       if ( (*ak4NCHS2)[a].pt()>15){
      float  deltaR = reco::deltaR((*Genj)[j].p4(), (*ak4NCHS2)[a].p4());
      if(deltaR< drGenCHS2[j]){
	drGenCHS2[j]= deltaR;
	 ptGenCHS2[j]=(*ak4NCHS2)[a].pt();
	 iMinDRNew2 = a;
      }
       }
    }

  
   if ( drGenCHS2[j]<0.4)
	{
	  const pat::Jet & Sak4NCHS2 = (*ak4NCHS2)[iMinDRNew2];

	  //  printf("matched Both Gen jet %d, pt=%6.2f, eta=%6.2f, phi=%6.2f /n", iMinDRNew2,Sak4NCHS2.pt(),Sak4NCHS2.eta(),Sak4NCHS2.phi());

 
	  const CompositePtrCandidate::daughters& daughters = Sak4NCHS2.daughterPtrVector();
	  
	  for ( CompositePtrCandidate::daughters::const_iterator daughter = daughters.begin();
		
		daughter != daughters.end(); ++daughter ) {
	    
	    TLorentzVector pippo;
	    
	    pippo.SetPtEtaPhiM((*daughter)->pt(),(*daughter)->eta() ,(*daughter)->phi(),(*daughter)->mass());
	    fromPV1.push_back(pippo);
	    
	    // const pat::PackedCandidate * dauPacked = dynamic_cast<const pat::PackedCandidate *>(&(**daughter));
	     
	     
	  }
	  
	 
	  
	  
	  const CompositePtrCandidate::daughters& daughters1 = Gj.daughterPtrVector();
	  
	  for ( CompositePtrCandidate::daughters::const_iterator daughter1 = daughters1.begin();
		daughter1 != daughters1.end(); ++daughter1) {
	    
	    
	    TLorentzVector pippo1;
	    
	    pippo1.SetPtEtaPhiM((*daughter1)->pt(),(*daughter1)->eta() ,(*daughter1)->phi(),(*daughter1)->mass());
	    Genjet1.push_back(pippo1);
	    
	  
	  
	 

	}

	

	  for(int n=Genjet1.size()-1;n>=0;n--)
	    {
	    
	      for (int m=fromPV1.size()-1;m>=0;m--)
		{
		  if (fabs(Genjet1[n].Pt() - fromPV1[m].Pt() ) < 0.01 && fabs(Genjet1[n].Eta() - fromPV1[m].Eta())<0.01 && fabs(Genjet1[n].Phi() - fromPV1[m].Phi() ) < 0.01 ){
		    
		    Genjet1.erase(Genjet1.begin()+n );
		    fromPV1.erase(fromPV1.begin()+m);
		    
		    // Genjet1_dz.erase(Genjet1_dz.begin()+n);
		    //  Genjet1_dzE.erase(Genjet1_dzE.begin()+n);
		    
		    continue;
		   
		  }
		}
	      
	      
	    }
	  
	  fromPV1.clear();  
	  Genjet1.clear();
	  //	Genjet1_dz.clear();  
	  // fromPV1_dz.clear();
	  //	Genjet1_dzE.clear();  
	  // fromPV1_dzE.clear();
	  
	  
	}//if

  
    //......................................................

  
    ptGenCHS3[j]= -1;
    drGenCHS3[j]= 100;
    int iMinDRNew3 = -1;
    vector<TLorentzVector>  BothCut1, Genjet;
     vector<float>  Genjet_dz, BothCut1_dz;
     vector<float>  Genjet_dzE, BothCut1_dzE;
     for (unsigned int a = 0, k=ak4NCHS3->size(); a < k;++a) {
       float  deltaR = reco::deltaR((*Genj)[j].p4(), (*ak4NCHS3)[a].p4());
       // if((*ak4NCHS3)[a].pt()>15){
       if(deltaR< drGenCHS3[j]){
	 drGenCHS3[j] = deltaR;
	 ptGenCHS3[j]=(*ak4NCHS3)[a].pt();
	 iMinDRNew3 = a;
       }
       // }
     }
     
     if(drGenCHS3[j]<0.4){
       const pat::Jet & Sak4NCHS3 = (*ak4NCHS3)[iMinDRNew3];

       // printf("matched Both Gen jet %d, pt=%6.2f, eta=%6.2f, phi=%6.2f /n", iMinDRNew3,Sak4NCHS3.pt(),Sak4NCHS3.eta(),Sak4NCHS3.phi());
       const CompositePtrCandidate::daughters& daughters = Sak4NCHS3.daughterPtrVector();
	
       for ( CompositePtrCandidate::daughters::const_iterator daughter = daughters.begin();
	     daughter != daughters.end(); ++daughter ) {
	 //const pat::PackedCandidate * dauPacked = dynamic_cast<const pat::PackedCandidate *>(&(**daughter));
	 TLorentzVector pippo;
	 
	 pippo.SetPtEtaPhiM((*daughter)->pt(),(*daughter)->eta() ,(*daughter)->phi(),(*daughter)->mass());
	 BothCut1.push_back(pippo);
	 
	 
	 
	 
	 
       }
       const CompositePtrCandidate::daughters& daughters1 = Gj.daughterPtrVector();
       
	for ( CompositePtrCandidate::daughters::const_iterator daughter1 = daughters1.begin();
	      daughter1 != daughters1.end(); ++daughter1) {
	  
	  TLorentzVector pippo1;
	  
	  pippo1.SetPtEtaPhiM((*daughter1)->pt(),(*daughter1)->eta() ,(*daughter1)->phi(),(*daughter1)->mass());
	  Genjet.push_back(pippo1);
	  
	  // const pat::PackedCandidate * dauPacked = dynamic_cast<const pat::PackedCandidate *>(&(**daughter1));
	  
	  

	}
	
	
	for(int n=Genjet.size()-1;n>=0;n--)
	  {
	    
	    for (int m=BothCut1.size()-1;m>=0;m--)
	      {
		if (fabs(Genjet[n].Pt() - BothCut1[m].Pt() ) < 0.01 && fabs(Genjet[n].Eta() - BothCut1[m].Eta())<0.01 && fabs(Genjet[n].Phi() - BothCut1[m].Phi() ) < 0.01 ){
		  
		  Genjet.erase(Genjet.begin()+n );
		  BothCut1.erase(BothCut1.begin()+m);
		  
		  // BothCut1_dz.erase(BothCut1_dz.begin()+n );
		  //  BothCut1_dzE.erase(BothCut1_dzE.begin()+n );
		  
		  // Genjet_dz.erase(Genjet_dz.begin()+m);
		  // Genjet_dzE.erase(Genjet_dzE.begin()+m);
		  continue;
		  
		  
		}
	      }
	  }
	
	
	BothCut1.clear();  
	Genjet.clear();  
	//	Genjet_dz.clear();  
	//	BothCut1_dz.clear();
	//	Genjet_dzE.clear();  
	//	BothCut1_dzE.clear();
	
     }//if


     //.........................Checking the producer

    
     ptGenCHS4[j]= -1;
     drGenCHS4[j]= 100;
     int iMinDRNew4 = -1;
     vector<TLorentzVector>  PR, Genjet2;
     vector<float>  Genjet2_dz, PR_dz;
     vector<float>  Genjet2_dzE, PR_dzE;
     for (unsigned int a = 0, k=ak4NCHS4->size(); a < k && a<njets_max ;++a) {
       //const reco::PFJet & Sak4NCHS4 = (*ak4NCHS4)[a];
      
      float  deltaR = reco::deltaR((*Genj)[j].p4(), (*ak4NCHS4)[a].p4());
       if ((*ak4NCHS4)[a].pt()>15){
      if(deltaR< drGenCHS4[j]){
	drGenCHS4[j] = deltaR;
	ptGenCHS4[j]= (*ak4NCHS4)[a].pt();
	iMinDRNew4 = a;
      }
       }
      
     }
      if( drGenCHS4[j]<0.4){
	const pat::Jet & Sak4NCHS4 = (*ak4NCHS4)[iMinDRNew4];

       //L//printf("matched Pr Gen jet %d, pt=%6.2f, eta=%6.2f, phi=%6.2f /n", iMinDRNew4,Sak4NCHS4.pt(),Sak4NCHS4.eta(),Sak4NCHS4.phi());
       const CompositePtrCandidate::daughters& daughters = Sak4NCHS4.daughterPtrVector();
	
       for ( CompositePtrCandidate::daughters::const_iterator daughter = daughters.begin();
	     daughter != daughters.end(); ++daughter ) {
	
	 TLorentzVector pippo;
	 
	 pippo.SetPtEtaPhiM((*daughter)->pt(),(*daughter)->eta() ,(*daughter)->phi(),(*daughter)->mass());
	 PR.push_back(pippo);
	 const pat::PackedCandidate * dauPacked = dynamic_cast<const pat::PackedCandidate *>(&(**daughter));



	 float dz = -999.;
	 float dzE = -999.;
	 if ( dauPacked->bestTrack() != 0 ){
	   dz = dauPacked->dz();
	   dzE = dauPacked->dzError();
	      }
	 
	      PR_dz.push_back(dz);
	      PR_dzE.push_back(dzE);
		
		
		}
       const CompositePtrCandidate::daughters& daughters1 = Gj.daughterPtrVector();
       
	for ( CompositePtrCandidate::daughters::const_iterator daughter1 = daughters1.begin();
	      daughter1 != daughters1.end(); ++daughter1) {
	  
	  TLorentzVector pippo1;
	  
	  pippo1.SetPtEtaPhiM((*daughter1)->pt(),(*daughter1)->eta() ,(*daughter1)->phi(),(*daughter1)->mass());
	  Genjet2.push_back(pippo1);
	  
	  
	  

	}
	
	for(int n=Genjet2.size()-1;n>=0;n--)
	  {
	    
	    for (int m=PR.size()-1;m>=0;m--)
	      {
		if (fabs(Genjet2[n].Pt() - PR[m].Pt() ) < 0.01 && fabs(Genjet2[n].Eta() - PR[m].Eta())<0.01 && fabs(Genjet2[n].Phi() - PR[m].Phi() ) < 0.01 ){
		  
		  Genjet2.erase(Genjet2.begin()+n );
		  PR.erase(PR.begin()+m);
		  // PR_dz.erase(BothCut1_dz.begin()+n );
		  // PR_dzE.erase(BothCut1_dzE.begin()+n );
		 
		  continue;
		  
		  
		}
	      }
	  }
	
	
	PR.clear();  
	Genjet2.clear();  
	PR_dz.clear();
	PR_dzE.clear();
	
	}//if

     
  
      //.......................  
      
      ptPFGen[j]= -1;
      jecPFGen[j]= -1;
      drPFGen[j] = 100; 
    // int iMinDRPF = -1;
    for (unsigned int i = 0, g=jets->size(); i < g ;++i){
      
      float  deltaR = reco::deltaR((*Genj)[j].p4(), (*jets)[i].p4());
      if ((*jets)[i].pt()*(*jets)[i].jecFactor("Uncorrected")>15  ){
      if(deltaR< drPFGen[j]){
	drPFGen[j] = deltaR;
	ptPFGen[j]= (*jets)[i].pt()*(*jets)[i].jecFactor("Uncorrected");
	
	jecPFGen[j]= (*jets)[i].jecFactor("Uncorrected");
	//	iMinDRPF = i;
      }
      }
    }
    
    }
    //L//printf("GenJet after %d pt=%3.1f eta=%5.2f \n", j ,Gj.pt(),Gj.eta());
  }
  

  //...................
 
  //....................    
  */
  event->Fill();
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
Jet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Jet::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
Jet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  }*/

//define this as a plug-in
DEFINE_FWK_MODULE(Jet);