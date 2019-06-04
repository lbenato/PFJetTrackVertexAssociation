###Pre-requisites:
###1. cmsenv in a CMSSW area
###2. load CRAB: source /cvmfs/cms.cern.ch/crab3/crab.sh
###3. set a proxy: voms-proxy-init --voms cms   (certificate password will be asked)

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os
config = config()

###Name of the main crab area
config.General.workArea = 'PFCHS_v0_debug_proxy'
###If not planning to transfer on T2, this could probably be commented out...
config.General.transferOutputs = True
config.General.transferLogs = True

###If running on already existing datasets: pluginName must be Analysis
config.JobType.pluginName = 'Analysis'
###Config file for cmsRun
config.JobType.psetName = 'python/ReclusterJets.py'
###If additional files are required (root files for weights, etc..), CRAB must load them
#config.JobType.inputFiles = ['data']

###Name of the working dir
config.General.requestName = 'QCD_Pt_80to120_TuneCP5_13TeV_pythia8-v1'

###Name of the sample in DAS
config.Data.inputDataset =  '/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'

#If private production:
#config.Data.inputDBS = 'phys03'
#If centrally produced
config.Data.inputDBS = 'global'

#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'Automatic'

###If real data, we can specify the json file, CRAB automatically takes care of it:
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
###this number could need to be tuned for larger datasets--CRAB cannot handle too many jobs, or few jobs with too many data (timeout)
config.Data.unitsPerJob = 10000
#### For ttbar: 20000 leads to error 50660 - Application terminated by wrapper because using too much RAM (RSS) 

#NEW: try to see if it stops after 100M events
#config.Data.totalUnits = 100000000

###In case you need to publish your output files in a T2 (they must be edm files)
config.Data.outLFNDirBase = '/store/user/lbenato/PFCHS_v0_debug_proxy'
config.Data.publication = False
#config.Data.outputDatasetTag = 'recluster_ak4Jets_miniaod_29May2018'#only for generation of EDM files
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'#only possible for generation of EDM files

###Storage site in case you want to write on a T2
config.Site.storageSite = 'T2_DE_DESY'

config.Site.blacklist   = ['T2_FR_IPHC']
