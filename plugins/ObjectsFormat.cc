
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
    I.eleE        = R->electronEnergy();
    I.eleEFrac    = R->hasUserFloat("eleEFrac") ? R->userFloat("eleEFrac") : -1.;//R->electronEnergyFraction();
    I.eleMulti    = R->electronMultiplicity();
    I.photonE     = R->photonEnergy();
    I.photonEFrac = R->hasUserFloat("photonEFrac") ? R->userFloat("photonEFrac") : -1.;//R->photonEnergyFraction();
    I.photonMulti = R->photonMultiplicity();
    I.cHadMulti   = R->chargedHadronMultiplicity();
    I.nHadMulti   = R->neutralHadronMultiplicity();
    I.npr         = R->chargedMultiplicity() + R->neutralMultiplicity();
    I.cMulti      = R->chargedMultiplicity();
    I.nMulti      = R->neutralMultiplicity();
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
    I.eleE          = -1.;
    I.eleEFrac      = -1.;
    I.eleMulti      = -1.;
    I.photonE       = -1.;
    I.photonEFrac   = -1.;
    I.photonMulti   = -1.;
    I.cHadMulti     = -1.;
    I.nHadMulti     = -1.;
    I.npr           = -1.;
    I.cMulti        = -1.;
    I.nMulti        = -1.;
    I.hasGenJ       =false;
    I.ptGenJ      = -10.;
    I.etaGenJ     = -4.;
    I.phiGenJ     = -4.;
    I.massGenJ    = -10.;
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

std::string ObjectsFormat::ListJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:cHadE/F:nHadE/F:cHadEFrac/F:nHadEFrac/F:nEmE/F:nEmEFrac/F:cEmE/F:cEmEFrac/F:cmuE/F:cmuEFrac/F:muE/F:muEFrac/F:eleE/F:eleEFrac/F:eleMulti/F:photonE/F:photonEFrac/F:photonMulti/F:cHadMulti/F:nHadMulti/F:npr/F:cMulti/F:nMulti/F:hasGenJ/O:ptGenJ/F:etaGenJ/F:phiGenJ/F:massGenJ/F:dRGenJ/F:response/F:multiGen/F:nMultiGen/F:cMultiGen/F:cHadEGen/F:nHadEGen/F:emEGen/F:eleEGen/F:photonEGen/F:ptGen/F:etaGen/F:phiGen/F:massGen/F:pdgIdGen/I:partonFlavour/I:hadronFlavour/I:mother/I:GenRecoMulti/F:GenRecoChMulti/F:GenRecoNeuMulti/F:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O";}
