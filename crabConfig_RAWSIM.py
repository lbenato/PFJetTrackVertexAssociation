from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'QCD_Pt_600to800_TuneCP5_13TeV_pythia8_generalTracks_2p5_10k_RAWSIM'
config.General.workArea = 'QCD_Pt_600to800_generalTracks_2p5_10k_RAWSIM'

config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.Data.inputDataset = '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIIFall17GS-93X_mc2017_realistic_v3_ext1-v2/GEN-SIM'
config.Data.inputDBS = 'global'#'phys03'

config.JobType.psetName = 'step_RAWSIM_cfg.py'
config.JobType.maxMemoryMB = 15900 #(alzare la memoria)
config.JobType.numCores = 8

config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1 #00 #100 events per job seem reasonable (500 events per job are too many)
config.Data.totalUnits = 10000 #10 k as a starting point
config.Data.outLFNDirBase = '/store/user/lbenato/PF_CHS_general_tracks_2p5/10k_events/'

#config.Data.publication = False
#!!!#
config.Data.publication = True
config.Data.outputDatasetTag = 'QCD_Pt_600to800_TuneCP5_13TeV_pythia8_generalTracks_2p5_10k_RAWSIM'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T2_IT_Pisa','T2_FR_IPHC','T3_US_UCR']
#config.Site.whitelist = ['T2_IT_Legnaro']#run the jobs were the minimum bias sample is stored!
