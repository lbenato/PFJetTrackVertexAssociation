// -*- C++ -*-
//
// Package:    Producer/SelectVertexProducer
// Class:      SelectVertexProducer
// 
/**\class SelectVertexProducer SelectVertexProducer.cc Producer/SelectVertexProducer/plugins/SelectVertexProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lisa Benato
//         Created:  Fri, 04 Oct 2019 11:42:28 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Wrapper.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Ptr.h"

//
// class declaration
//

class SelectVertexProducer : public edm::stream::EDProducer<> {
   public:
      explicit SelectVertexProducer(const edm::ParameterSet&);
      ~SelectVertexProducer();
      virtual double vertex_ptmax2(const reco::Vertex&);
      virtual bool select(const reco::Vertex &, int);

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      edm::EDGetTokenT< std::vector<reco::Vertex> > PVToken_;
      //edm::InputTag PVToken_;
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
SelectVertexProducer::SelectVertexProducer(const edm::ParameterSet& iConfig):
   PVToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter <edm::InputTag>("vertices")))
{
   //register your products
   //PVToken_  = iConfig.getParameter<edm::InputTag>( "src" );
  produces<std::vector<reco::Vertex>>("");//.setBranchAlias("LisaNewVertices");   
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


SelectVertexProducer::~SelectVertexProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SelectVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<reco::VertexCollection> PVCollection;
   iEvent.getByToken(PVToken_, PVCollection);

   std::vector<reco::Vertex> LisaVector;
   //std::unique_ptr<std::vector<reco::Vertex>> LisaNewVertices;
   //int vert_count = 0;

   for(std::vector<reco::Vertex>::const_iterator it=PVCollection->begin(); it!=PVCollection->end(); ++it){
     reco::Vertex v=*it;
     //std::cout << "\n" << std::endl;
     //std::cout << "Vertex n. " << vert_count << std::endl;
     //std::cout << v.p4() << std::endl;

     //////option 0: usual ndof cut
     ////if(select(v, 0)){
     ////  nPVsel_ndof++;
     ////}

     ////option 1: Wolfram's selection
     if(select(v, 1)){
       //std:: cout << "selected as good!" << std::endl;
       LisaVector.push_back(v);
       //nPVsel_pt++;
     }

     //vert_count++;

   }

   //int vert_count_good = 0;
   //for(unsigned int a = 0; a<LisaVector.size(); a++){
     //std::cout << "\n" << std::endl;
     //std::cout << "Vertex good n. " << vert_count_good << std::endl;
     //std::cout << LisaVector.at(a).p4() << std::endl;
     //vert_count_good++;
   //}

   auto LisaCandidates = std::make_unique<std::vector<reco::Vertex>>(LisaVector);
   iEvent.put(std::move(LisaCandidates),"");

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SelectVertexProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SelectVertexProducer::endStream() {
}


// ------------ method called when starting to processes a run  ------------
/*
void
SelectVertexProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SelectVertexProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SelectVertexProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SelectVertexProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SelectVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//Thanks to Wolfram Erdmann!

double SelectVertexProducer::vertex_ptmax2(const reco::Vertex & v) {
  double ptmax1 = 0;
  double ptmax2 = 0;
  //std::cout << "vertex_ptmax2 in action!" << std::endl;
  //std::cout << "vertex had ntracks: " << v.nTracks(0.5) << std::endl;
  //std::cout << "vertex has n dof: " << v.ndof() << std::endl;
  //std::cout << "vertex has track size: " << v.tracksSize() << std::endl;

  for(reco::Vertex::trackRef_iterator t = v.tracks_begin(); t != v.tracks_end(); t++){
    if (v.trackWeight(*t) > 0.5){
      double pt = t->get()->pt();
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

bool SelectVertexProducer::select(const reco::Vertex & v, int level){
  /* level
     0  !isFake  && ndof>4
     1  !isFake  && ndof>4 &&  ptmax2 >0.4
  */
  //std::cout << SelectVertexProducer::vertex_ptmax2(v) << std::endl;

  if( v.isFake() ) return false;
  if( (level == 0) && (v.ndof()>4) ) return true;
  if( (level == 1) && (v.ndof()>4) && (SelectVertexProducer::vertex_ptmax2(v)>0.4) ) return true;
  return false;
}



//define this as a plug-in
DEFINE_FWK_MODULE(SelectVertexProducer);
