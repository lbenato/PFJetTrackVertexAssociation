from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'QCD_Pt_600to800_TuneCP5_13TeV_pythia8_generalTracks_2p5_10k_AODSIM'
config.General.workArea = 'QCD_Pt_600to800_generalTracks_2p5_10k_AODSIM'

config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.Data.inputDataset = '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/lbenato-QCD_Pt_600to800_TuneCP5_13TeV_pythia8_generalTracks_2p5_10k_RAWSIM-424e4485a07f26f554e82f829d793003/USER'
config.Data.inputDBS = 'phys03'

config.JobType.psetName = 'step_AODSIM_generalTracks_cfg.py'
config.JobType.maxMemoryMB = 15900 #(more memory)
config.JobType.numCores = 8

config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1 #00 #100 events per job seem reasonable (500 events per job are too many)
#config.Data.totalUnits = 10000 #10 k as a starting point
config.Data.outLFNDirBase = '/store/user/lbenato/PF_CHS_general_tracks_2p5/10k_events/'

#config.Data.publication = False
#!!!#
config.Data.publication = True
config.Data.outputDatasetTag = 'QCD_Pt_600to800_TuneCP5_13TeV_pythia8_generalTracks_2p5_10k_AODSIM'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'

config.Site.storageSite = 'T2_DE_DESY'
#config.Site.blacklist = ['T2_IT_Pisa','T2_FR_IPHC','T3_US_UCR']
#config.Site.whitelist = ['T2_IT_Legnaro']#run the jobs were the minimum bias sample is stored!
