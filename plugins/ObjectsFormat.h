#ifndef OBJECTSFORMAT_H
#define OBJECTSFORMAT_H

#include <string>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

//#include "Utilities.h"
#include "Objects.h"

class ObjectsFormat {
    
    public:
        ObjectsFormat() {};
        ~ObjectsFormat() {};
        
        static void FillJetType(JetType&, const pat::Jet*, bool);
        static void ResetJetType(JetType&);
        static std::string ListJetType();
        
    private:
    
};


#endif

