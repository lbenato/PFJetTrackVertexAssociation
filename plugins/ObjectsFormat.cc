
#include "ObjectsFormat.h"


//*******************//
//        Jets       //
//*******************//

void ObjectsFormat::FillJetType(JetType& I, const pat::Jet* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.hasGenJ     = R->genJet() ? true : false;
    if(isMC && R->genJet()) {
      I.ptGenJ    = R->genJet()->pt();
      I.etaGenJ   = R->genJet()->eta();
      I.phiGenJ   = R->genJet()->phi();
      I.massGenJ  = R->genJet()->mass();
      I.dRGenJ    = reco::deltaR(R->eta(),R->phi(),R->genJet()->eta(),R->genJet()->phi());//!!
      I.response  = (reco::deltaR(R->eta(),R->phi(),R->genJet()->eta(),R->genJet()->phi()) < 0.4 && R->genJet()->pt()>0.) ? R->pt()/R->genJet()->pt() : -1.;//!!
      I.multiGen  = R->hasUserFloat("multiGen") ? R->userFloat("multiGen") : -1.;
      I.nMultiGen = R->hasUserFloat("nMultiGen") ? R->userFloat("nMultiGen") : -1.;
      I.cMultiGen = R->hasUserFloat("cMultiGen") ? R->userFloat("cMultiGen") : -1.;
      I.cHadEGen  = R->hasUserFloat("cHadEGen") ? R->userFloat("cHadEGen") : -1.;
      I.nHadEGen  = R->hasUserFloat("nHadEGen") ? R->userFloat("nHadEGen") : -1.;
      I.emEGen    = R->hasUserFloat("emEGen") ? R->userFloat("emEGen") : -1.;
      I.eleEGen   = R->hasUserFloat("eleEGen") ? R->userFloat("eleEGen") : -1.;
      I.photonEGen= R->hasUserFloat("photonEGen") ? R->userFloat("photonEGen") : -1.;

    }
    if(isMC && R->genParton()) {
      I.ptGen     = R->genParton()->pt();
      I.etaGen    = R->genParton()->eta();
      I.phiGen    = R->genParton()->phi();
      I.massGen   = R->genParton()->mass();
      I.pdgIdGen   = R->genParton()->pdgId();
    }
    I.cHadE       = R->chargedHadronEnergy();
    I.nHadE       = R->neutralHadronEnergy();
    I.cHadEFrac   = R->hasUserFloat("cHadEFrac") ? R->userFloat("cHadEFrac") : -1.;//R->chargedHadronEnergyFraction();
    I.nHadEFrac   = R->hasUserFloat("nHadEFrac") ? R->userFloat("nHadEFrac") : -1.;//R->neutralHadronEnergyFraction();
    I.nEmE        = R->neutralEmEnergy();
    I.nEmEFrac    = R->hasUserFloat("nEmEFrac") ? R->userFloat("nEmEFrac") : -1.;//R->neutralEmEnergyFraction();
    I.cEmE        = R->chargedEmEnergy();
    I.cEmEFrac    = R->hasUserFloat("cEmEFrac") ? R->userFloat("cEmEFrac") : -1.;//R->chargedEmEnergyFraction();
    I.cmuE        = R->chargedMuEnergy();
    I.cmuEFrac    = R->hasUserFloat("cmuEFrac") ? R->userFloat("cmuEFrac") : -1.;//R->chargedMuEnergyFraction();
    I.muE         = R->muonEnergy();
    I.muEFrac     = R->hasUserFloat("muEFrac") ? R->userFloat("muEFrac") : -1.;//R->muonEnergyFraction();
    I.muMulti     = R->muonMultiplicity();
    I.muMultiFrac = float(R->muonMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.eleE        = R->electronEnergy();
    I.eleEFrac    = R->hasUserFloat("eleEFrac") ? R->userFloat("eleEFrac") : -1.;//R->electronEnergyFraction();
    I.eleMulti    = R->electronMultiplicity();
    I.eleMultiFrac = float(R->electronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.photonE     = R->photonEnergy();
    I.photonEFrac = R->hasUserFloat("photonEFrac") ? R->userFloat("photonEFrac") : -1.;//R->photonEnergyFraction();
    I.photonMulti = R->photonMultiplicity();
    I.photonMultiFrac = float(R->photonMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.cHadMulti   = R->chargedHadronMultiplicity();
    I.cHadMultiFrac = float(R->chargedHadronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.nHadMulti   = R->neutralHadronMultiplicity();
    I.nHadMultiFrac = float(R->neutralHadronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.npr         = R->chargedMultiplicity() + R->neutralMultiplicity();
    I.cMulti      = R->chargedMultiplicity();
    I.cMultiFrac  = float(R->chargedMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.nMulti      = R->neutralMultiplicity();
    I.nMultiFrac  = float(R->neutralMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.partonFlavour     = R->partonFlavour();
    I.hadronFlavour     = R->hadronFlavour();
    I.mother = 0;
    I.GenRecoMulti      = (R->genJet() && R->hasUserFloat("multiGen") && R->userFloat("multiGen")>0) ? (R->chargedMultiplicity() + R->neutralMultiplicity())/R->userFloat("multiGen")  : -1.;
    I.GenRecoChMulti    =(R->genJet() && R->hasUserFloat("cMultiGen") && R->userFloat("cMultiGen")>0) ? (R->chargedMultiplicity())/R->userFloat("cMultiGen")  : -1.;
    I.GenRecoNeuMulti   =(R->genJet() && R->hasUserFloat("nMultiGen") && R->userFloat("nMultiGen")>0) ? (R->neutralMultiplicity())/R->userFloat("nMultiGen")  : -1.;
    //if(isMC && R->genParton()) I.mother = Utilities::FindMotherId(R->genParton());
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isTightLepVeto     = R->hasUserInt("isTightLepVeto") ? R->userInt("isTightLepVeto") : false;
}

void ObjectsFormat::ResetJetType(JetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.cHadE         = -1.;
    I.nHadE         = -1.;
    I.cHadEFrac     = -1.;
    I.nHadEFrac     = -1.;
    I.nEmE          = -1.;
    I.nEmEFrac      = -1.;
    I.cEmE          = -1.;
    I.cEmEFrac      = -1.;
    I.cmuE          = -1.;
    I.cmuEFrac      = -1.;
    I.muE           = -1.;
    I.muEFrac       = -1.;
    I.muMulti       = -1.;
    I.muMultiFrac   = -1.;
    I.eleE          = -1.;
    I.eleEFrac      = -1.;
    I.eleMulti      = -1.;
    I.eleMultiFrac   = -1.;
    I.photonE       = -1.;
    I.photonEFrac   = -1.;
    I.photonMulti   = -1.;
    I.photonMultiFrac   = -1.;
    I.cHadMulti     = -1.;
    I.cHadMultiFrac = -1.;
    I.nHadMulti     = -1.;
    I.nHadMultiFrac = -1.;
    I.npr           = -1.;
    I.cMulti        = -1.;
    I.cMultiFrac    = -1.;
    I.nMulti        = -1.;
    I.nMultiFrac    = -1.;
    I.hasGenJ       =false;
    I.ptGenJ      = -1.;
    I.etaGenJ     = -9.;
    I.phiGenJ     = -9.;
    I.massGenJ    = -1.;
    I.dRGenJ      = 999.;
    I.response    = -1.;
    I.multiGen    = -1.;
    I.nMultiGen   = -1.;
    I.cMultiGen   = -1.;
    I.cHadEGen    = -1.;
    I.nHadEGen    = -1.;
    I.emEGen      = -1.;
    I.eleEGen     = -1.;
    I.photonEGen  = -1.;
    I.ptGen       = -10.;
    I.etaGen      = -4.;
    I.phiGen      = -4.;
    I.massGen     = -10.;
    I.pdgIdGen     = 0.;
    I.partonFlavour     = 0;
    I.hadronFlavour     = 0;
    I.mother            = false;
    I.GenRecoMulti      = -1.;
    I.GenRecoChMulti    = -1.;
    I.GenRecoNeuMulti   = -1.;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isTightLepVeto     = false;
}

std::string ObjectsFormat::ListJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:cHadE/F:nHadE/F:cHadEFrac/F:nHadEFrac/F:nEmE/F:nEmEFrac/F:cEmE/F:cEmEFrac/F:cmuE/F:cmuEFrac/F:muE/F:muEFrac/F:muMulti/F:muMultiFrac/F:eleE/F:eleEFrac/F:eleMulti/F:eleMultiFrac/F:photonE/F:photonEFrac/F:photonMulti/F:photonMultiFrac/F:cHadMulti/F:cHadMultiFrac/F:nHadMulti/F:nHadMultiFrac/F:npr/F:cMulti/F:cMultiFrac/F:nMulti/F:nMultiFrac/F:hasGenJ/O:ptGenJ/F:etaGenJ/F:phiGenJ/F:massGenJ/F:dRGenJ/F:response/F:multiGen/F:nMultiGen/F:cMultiGen/F:cHadEGen/F:nHadEGen/F:emEGen/F:eleEGen/F:photonEGen/F:ptGen/F:etaGen/F:phiGen/F:massGen/F:pdgIdGen/I:partonFlavour/I:hadronFlavour/I:mother/I:GenRecoMulti/F:GenRecoChMulti/F:GenRecoNeuMulti/F:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O";}


//*******************//
//    PF Candidates  //
//*******************//

void ObjectsFormat::FillPFCandidateType(PFCandidateType& I, const pat::PackedCandidate* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.dz          = abs(R->dz());
    I.cosh        = cosh(R->eta())*(0.02+0.01/R->pt())*5;
    I.isNew1      = abs(R->dz())<cosh(R->eta())*(0.02+0.01/R->pt())*5 ? true : false;
    I.isCHS       = (R->fromPV()>0 || R->charge()==0) ? true : false;
    I.isNew3      = (abs(R->dz()) < cosh(R->eta())*(0.02+0.01/R->pt())*5) ||  (R->hasTrackDetails() && abs(R->dz())<R->pt()*R->dzError()) ? true : false;
    I.isNew4      = (R->vertexRef().key()==0 ||( R->dzAssociatedPV() < cosh(R->eta() )*(0.003+0.01/R->pt())*5) || R->charge()==0) ? true : false;
    I.isPuppi     = (R->fromPV()==3 || ( (R->fromPV()==2 || R->fromPV()==1) && abs(R->dz()) < 0.3 )) ? true : false;
}

void ObjectsFormat::ResetPFCandidateType(PFCandidateType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.dz          = -9.;
    I.cosh        = -9.;
    I.isNew1      = false;
    I.isCHS       = false;
    I.isNew3      = false;
    I.isNew4      = false;
    I.isPuppi     = false;
}

std::string ObjectsFormat::ListPFCandidateType() {return "pt/F:eta/F:phi/F:dz/F:cosh/F:isNew1/O:isCHS/O:isNew3/O:isNew4/O:isPuppi/O";}
