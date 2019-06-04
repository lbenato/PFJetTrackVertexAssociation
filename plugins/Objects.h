#ifndef OBJECTS_H
#define OBJECTS_H

struct JetType {
JetType(): pt(-1.), eta(-9.), phi(-9.), mass(-1.), energy(-1.), cHadE(-1.), nHadE (-1.), cHadEFrac(-1.), nHadEFrac(-1.), nEmE(-1.), nEmEFrac(-1.), cEmE(-1.), cEmEFrac(-1.), cmuE(-1.), cmuEFrac(-1.), muE(-1.), muEFrac(-1.), eleE(-1.), eleEFrac(-1.), eleMulti(-1.), photonE(-1.), photonEFrac(-1.), photonMulti(-1.), cHadMulti(-1.), nHadMulti(-1.), npr(-1.), cMulti(-1.), nMulti(-1.), hasGenJ(false), ptGenJ(-10.), etaGenJ(-4.), phiGenJ(-4.), massGenJ(-10.), dRGenJ(999.), response (-1.), multiGen(-1.), nMultiGen(-1.), cMultiGen(-1.), cHadEGen(-1.), nHadEGen(-1.), emEGen(-1.), eleEGen(-1.), photonEGen(-1.), ptGen(-10.), etaGen(-4.), phiGen(-4.), massGen(-10.), pdgIdGen(0.),partonFlavour(0), hadronFlavour(0), mother(0), GenRecoMulti(-1.), GenRecoChMulti(-1.), GenRecoNeuMulti(-1.), isLoose(false), isMedium(false), isTight(false), isTightLepVeto(false) {}
  float pt;
  float eta;
  float phi;
  float mass;
  float energy;
  float cHadE;
  float nHadE;
  float cHadEFrac;
  float nHadEFrac;
  float nEmE;
  float nEmEFrac;
  float cEmE;
  float cEmEFrac;
  float cmuE;
  float cmuEFrac;
  float muE;
  float muEFrac;
  float eleE;
  float eleEFrac;
  float eleMulti;
  float photonE;
  float photonEFrac;
  float photonMulti;
  float cHadMulti;
  float nHadMulti;
  float npr;
  float cMulti;
  float nMulti;
  bool hasGenJ;
  float ptGenJ;
  float etaGenJ;
  float phiGenJ;
  float massGenJ;
  float dRGenJ;
  float response;
  float multiGen;
  float nMultiGen;
  float cMultiGen;
  float cHadEGen;
  float nHadEGen;
  float emEGen;
  float eleEGen;
  float photonEGen;
  float ptGen;
  float etaGen;
  float phiGen;
  float massGen;
  int pdgIdGen;
  int partonFlavour;
  int hadronFlavour;
  int mother;
  float GenRecoMulti;
  float GenRecoChMulti;
  float GenRecoNeuMulti;
  bool isLoose;
  bool isMedium;
  bool isTight;
  bool isTightLepVeto;
};

#endif