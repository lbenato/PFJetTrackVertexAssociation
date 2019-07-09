#! /usr/bin/env python 

sample = {
    'QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1': {
        'xsec' : 246300000.0,#pb
    },
}

lista = []
#for a in range(1,100,1):
#for a in range(3,100,1):
for a in range(1,2,1):
    name = 'output_' + str(a)
    lista.append(name)


samples = {
    'QCD' : {
        'files' : lista,#['output_1'],
        'fillcolor' : 881,
        'fillstyle' : 1001,
        'linecolor' : 881,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : 'QCD',
        'weight': 1.,
        'plot': True,
    },
    'QCD_Pt_15to30' : {
#        'files' : ['QCD_Pt_15to30'],
#        'files' : lista,#['output_1'],
        'files' : ['QCD_Pt_15to30_TuneCP5_13TeV_pythia8-v1'],
        'fillcolor' : 881,
        'fillstyle' : 1001,
        'linecolor' : 881,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : 'QCD Pt 15to30',
        'weight': 1.,
        'plot': True,
    },
    'QCD_Pt_30to50' : {
#        'files' : ['QCD_Pt_30to50'],
#        'files' : lista,#['output_1'],
        'files' : ['QCD_Pt_30to50_TuneCP5_13TeV_pythia8-v1'],
        'fillcolor' : 881,
        'fillstyle' : 1001,
        'linecolor' : 881,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : 'QCD Pt 30to50',
        'weight': 1.,
        'plot': True,
    },
    'QCD_Pt_80to120' : {
#        'files' : ['QCD_Pt_30to50'],
#        'files' : lista,#['output_1'],
        'files' : ['QCD_Pt_80to120_TuneCP5_13TeV_pythia8-v1'],
        'fillcolor' : 881,
        'fillstyle' : 1001,
        'linecolor' : 881,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : 'QCD Pt 80to120',
        'weight': 1.,
        'plot': True,
    },
}
