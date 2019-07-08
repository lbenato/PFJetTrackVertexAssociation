variable = {}

var_template = {
    "EventNumber": {
      "title" : "event number",
      "nbins" : 10000000,
      "min" : 0,
      "max" : 1.e7,
      "log" : False,
    },
    "LumiNumber": {
      "title" : "lumisection number",
      "nbins" : 2000,
      "min" : 0,
      "max" : 2000,
      "log" : False,
    },
    "RunNumber": {
      "title" : "run number",
      "nbins" : 7000,
      "min" : 254000,
      "max" : 261000,
      "log" : False,
    },
    "nPV": {
      "title" : "number of reconstructed Primary Vertices",
      "nbins" : 50,
      "min" : -0.5,
      "max" : 49.5,
      "log" : False,
    },
    "nJets": {
      "title" : "number of jets",
      "nbins" : 20,
      "min" : -0.5,
      "max" : 19.5,
      "log" : True,
    },
    "nJetsNew1": {
      "title" : "number of jets",
      "nbins" : 20,
      "min" : -0.5,
      "max" : 19.5,
      "log" : True,
    },
    "nJetsNew2": {
      "title" : "number of jets",
      "nbins" : 20,
      "min" : -0.5,
      "max" : 19.5,
      "log" : True,
    },
    "nJetsNew3": {
      "title" : "number of jets",
      "nbins" : 20,
      "min" : -0.5,
      "max" : 19.5,
      "log" : True,
    },
    "nJetsNew4": {
      "title" : "number of jets",
      "nbins" : 20,
      "min" : -0.5,
      "max" : 19.5,
      "log" : True,
    },
    "nJetsNew5": {
      "title" : "number of jets",
      "nbins" : 20,
      "min" : -0.5,
      "max" : 19.5,
      "log" : True,
    },
    "Jets.pt": {
        "title" : "all jets p_{T} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew1.pt": {
        "title" : "all jets p_{T} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew2.pt": {
        "title" : "all jets p_{T} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew3.pt": {
        "title" : "all jets p_{T} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew4.pt": {
        "title" : "all jets p_{T} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew5.pt": {
        "title" : "all jets p_{T} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },



    "Jets.eta": {
        "title" : "jet #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "JetsNew1.eta": {
        "title" : "jet #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "JetsNew2.eta": {
        "title" : "jet #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "JetsNew3.eta": {
        "title" : "jet #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "JetsNew4.eta": {
        "title" : "jet #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "JetsNew5.eta": {
        "title" : "jet #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },


    "fabs(Jets.eta)": {
        "title" : "all jets #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "fabs(JetsNew1.eta)": {
        "title" : "all jets #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "fabs(JetsNew2.eta)": {
        "title" : "all jets #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "fabs(JetsNew3.eta)": {
        "title" : "all jets #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "fabs(JetsNew4.eta)": {
        "title" : "all jets #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },
    "fabs(JetsNew5.eta)": {
        "title" : "all jets #eta",
        "nbins" : 50,
        "min" : -5.2,
        "max" : 5.2,
        "log" : True,
    },



    "Jets.Deltapt": {
        "title" : ".... p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },


    "Jets.ptGenJ": {
        "title" : "all jets p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew1.ptGenJ": {
        "title" : "all jets p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew2.ptGenJ": {
        "title" : "all jets p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew3.ptGenJ": {
        "title" : "all jets p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew4.ptGenJ": {
        "title" : "all jets p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },
    "JetsNew5.ptGenJ": {
        "title" : "all jets p_{T}^{gen} (GeV)",
        "nbins" : 40,
        "min" : 0,
        "max" : 200,
        "log" : True,
    },



    "Jets.response": {
        "title" : "jets p_{T}^{reco}/p_{T}^{gen}",
        "nbins" : 50,
        "min" : 0,
#        "max" : 5,
        "max" : 3,
        "log" : True,
    },
    "JetsNew1.response": {
        "title" : "jets p_{T}^{reco}/p_{T}^{gen}",
        "nbins" : 50,
        "min" : 0,
#        "max" : 5,
        "max" : 3,
        "log" : True,
    },
    "JetsNew2.response": {
        "title" : "jets p_{T}^{reco}/p_{T}^{gen}",
        "nbins" : 50,
        "min" : 0,
#        "max" : 5,
        "max" : 3,
        "log" : True,
    },
    "JetsNew3.response": {
        "title" : "jets p_{T}^{reco}/p_{T}^{gen}",
        "nbins" : 50,
        "min" : 0,
#        "max" : 5,
        "max" : 3,
        "log" : True,
    },
    "JetsNew4.response": {
        "title" : "jets p_{T}^{reco}/p_{T}^{gen}",
        "nbins" : 50,
        "min" : 0,
#        "max" : 5,
        "max" : 3,
        "log" : True,
    },
    "JetsNew5.response": {
        "title" : "jets p_{T}^{reco}/p_{T}^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },


    "Jets.GenRecoNeuMulti": {
        "title" : "jet neutral multiplicity^{reco}/neutral multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew1.GenRecoNeuMulti": {
        "title" : "jet neutral multiplicity^{reco}/neutral multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew2.GenRecoNeuMulti": {
        "title" : "jet neutral multiplicity^{reco}/neutral multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew3.GenRecoNeuMulti": {
        "title" : "jet neutral multiplicity^{reco}/neutral multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew4.GenRecoNeuMulti": {
        "title" : "jet neutral multiplicity^{reco}/neutral multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew5.GenRecoNeuMulti": {
        "title" : "jet neutral multiplicity^{reco}/neutral multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },


    "Jets.GenRecoChMulti": {
        "title" : "jet charged multiplicity^{reco}/charged multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew1.GenRecoChMulti": {
        "title" : "jet charged multiplicity^{reco}/charged multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew2.GenRecoChMulti": {
        "title" : "jet charged multiplicity^{reco}/charged multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew3.GenRecoChMulti": {
        "title" : "jet charged multiplicity^{reco}/charged multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew4.GenRecoChMulti": {
        "title" : "jet charged multiplicity^{reco}/charged multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },
    "JetsNew5.GenRecoChMulti": {
        "title" : "jet charged multiplicity^{reco}/charged multiplicity^{gen}",
        "nbins" : 50,
        "min" : 0,
        "max" : 5,
        "log" : True,
    },



    "Jets.nHadEFrac": {
        "title" : "jet neutral hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew1.nHadEFrac": {
        "title" : "jet neutral hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew2.nHadEFrac": {
        "title" : "jet neutral hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew3.nHadEFrac": {
        "title" : "jet neutral hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew4.nHadEFrac": {
        "title" : "jet neutral hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew5.nHadEFrac": {
        "title" : "jet neutral hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },

    "Jets.cHadEFrac": {
        "title" : "jet charged hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew1.cHadEFrac": {
        "title" : "jet charged hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew2.cHadEFrac": {
        "title" : "jet charged hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew3.cHadEFrac": {
        "title" : "jet charged hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew4.cHadEFrac": {
        "title" : "jet charged hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew5.cHadEFrac": {
        "title" : "jet charged hadron energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },


    "Jets.nEmEFrac": {
        "title" : "jet neutral electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew1.nEmEFrac": {
        "title" : "jet neutral electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew2.nEmEFrac": {
        "title" : "jet neutral electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew3.nEmEFrac": {
        "title" : "jet neutral electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew4.nEmEFrac": {
        "title" : "jet neutral electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew5.nEmEFrac": {
        "title" : "jet neutral electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },

    "Jets.cEmEFrac": {
        "title" : "jet charged electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew1.cEmEFrac": {
        "title" : "jet charged electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew2.cEmEFrac": {
        "title" : "jet charged electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew3.cEmEFrac": {
        "title" : "jet charged electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew4.cEmEFrac": {
        "title" : "jet charged electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew5.cEmEFrac": {
        "title" : "jet charged electromagnetic energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },


    "Jets.muEFrac": {
        "title" : "jet muon energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew1.muEFrac": {
        "title" : "jet muon energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew2.muEFrac": {
        "title" : "jet muon energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew3.muEFrac": {
        "title" : "jet muon energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew4.muEFrac": {
        "title" : "jet muon energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },
    "JetsNew5.muEFrac": {
        "title" : "jet muon energy fraction",
        "nbins" : 50,
        "min" : 0,
        "max" : 1,
        "log" : True,
    },








    "Jets.hasGenJ": {
        "title" : "has gen jet",
        "nbins" : 2,
        "min" : -0.5,
        "max" : 1.5,
        "log" : False,
    },
    "JetsNew1.hasGenJ": {
        "title" : "has gen jet",
        "nbins" : 2,
        "min" : -0.5,
        "max" : 1.5,
        "log" : False,
    },
    "JetsNew2.hasGenJ": {
        "title" : "has gen jet",
        "nbins" : 2,
        "min" : -0.5,
        "max" : 1.5,
        "log" : False,
    },
    "JetsNew3.hasGenJ": {
        "title" : "has gen jet",
        "nbins" : 2,
        "min" : -0.5,
        "max" : 1.5,
        "log" : False,
    },
    "JetsNew4.hasGenJ": {
        "title" : "has gen jet",
        "nbins" : 2,
        "min" : -0.5,
        "max" : 1.5,
        "log" : False,
    },
    "JetsNew5.hasGenJ": {
        "title" : "has gen jet",
        "nbins" : 2,
        "min" : -0.5,
        "max" : 1.5,
        "log" : False,
    },

}



for n, v in var_template.iteritems():
    if '[N]' in n:
        for i in range(0, 5):
            ni = n.replace('[N]', "%d" % i)
            variable[ni] = v.copy()
            variable[ni]['title'] = variable[ni]['title'].replace('[N]', "%d" % i)
    else:
        variable[n] = v
