#Here: standard crab config file

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ReclusterJets.py'
#config.JobType.inputFiles = ['data']

config.General.requestName = 'QCD_test'

config.Data.inputDataset =  '/VBFH_HToSSTobbbb_MH-125_MS-40_ctauS-0_TuneCUETP8M1_13TeV-powheg-pythia8_PRIVATE-MC/lbenato-RunIISummer16-PU_standard_mixing-Moriond17_80X_mcRun2_2016_MINIAOD-28028af67189b3de7224b79195bd0e1d/USER'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'

config.Data.unitsPerJob = 10000
config.Data.outLFNDirBase = '/store/user/lbenato/choose_a_folder_name'
config.Data.publication = False

config.Site.storageSite = 'T2_DE_DESY'
config.Site.blacklist   = ['T2_FR_IPHC']
            

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    # Selection of samples via python lists
    import os
    from Analyzer.PFJetTrackVertexAssociation.samples import sample, samples

    list_of_samples = ["QCD","QCD_Pt_15to30","QCD_Pt_30to50","QCD_Pt_80to120"]#,"data_obs"
    print "Possible subgroups of samples:"
    for a in list_of_samples:
        print a
    print "---------------"

    ########parser#######
    import optparse
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("-a", "--crabaction", action="store", type="string", dest="crabaction", default="test")
    parser.add_option("-l", "--lists", action="store", type="string", dest="lists", default="")
    parser.add_option("-g", "--groupofsamples", action="store", type="string", dest="groupofsamples", default="")
    (options, args) = parser.parse_args()


    ####################
    folder = ''
    pset = ''
    workarea = ''
    from Analyzer.PFJetTrackVertexAssociation.crab_requests_lists import * #This list is fine for us!
    #crabConfig = 'crabConfig.py'
    pset = "ReclusterJets.py"
    folder = "test" #CHANGE here your crab folder name
    outLFNDirBase = "/store/user/lbenato/"+folder #CHANGE here according to your username!
    workarea = "/nfs/dust/cms/user/lbenato/" + folder #CHANGE here according to your username!


    ####################!!!!!####!#!#!#!#!#!#!

    selected_requests = {}
    if options.groupofsamples not in list_of_samples:
        print "Invalid subgroup of samples, aborting!"
        exit()

    for b, k in enumerate(requests.keys()):
        if k in samples[options.groupofsamples]["files"]:
            print k
            selected_requests[k] = requests[k]
    
    for a, j in enumerate(selected_requests):

        # submission of the python config
        if options.crabaction=="submit":
            config.Data.inputDBS = "global"

            os.system('echo submitting this config...\n')
            #modify parameters here
            config.General.requestName = j
            config.Data.inputDataset = selected_requests[j]
            config.JobType.psetName = "python/" + pset
            config.Data.outLFNDirBase = outLFNDirBase
            config.General.workArea= workarea
            #if isData:
            #    config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
            #    #config.Data.splitting = 'Automatic'
            #    config.Data.unitsPerJob = 100000
            #config.JobType.pyCfgParams = [string_runLocal, string_isData, string_isREHLT, string_isReReco, string_isReMiniAod, string_isPromptReco,string_noLHEinfo, string_isbbH, string_isSignal, string_GT, string_JECstring, string_jsonName, string_triggerTag, string_filterString]
            print config
            submit(config)

        elif options.crabaction=="status":
            os.system('echo status -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab status -d ' + workarea + '/crab_'+j+'\n')
            os.system('echo ----------------------------------------------------\n') 
        elif options.crabaction=="resubmit":
            os.system('echo resubmit -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab resubmit -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="getoutput":
            os.system('echo getoutput -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab getoutput -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="kill":
            os.system('echo kill -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab kill -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="report":
            os.system('echo report -d ' + workarea + '/crab_'+j+'\n')
            os.system('crab report -d ' + workarea + '/crab_'+j+'\n')
        elif options.crabaction=="test":
            if "VBFH_HToSS" in j:
                #automatic implementation of the choice bewteen inputDBS global/phys03
                config.Data.inputDBS = "phys03"
            else:
                config.Data.inputDBS = "global"

            os.system('echo submitting this config...\n')
            #modify parameters here
            config.General.requestName = j
            config.Data.inputDataset = selected_requests[j]
            config.JobType.psetName = "python/" + pset
            config.Data.outLFNDirBase = outLFNDirBase
            config.General.workArea= workarea
            #if isData:
            #    config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
            #    #config.Data.splitting = 'Automatic'
            #    config.Data.unitsPerJob = 100000
            #config.JobType.pyCfgParams = [string_runLocal, string_isData, string_isREHLT, string_isReReco, string_isReMiniAod, string_isPromptReco,string_noLHEinfo, string_isbbH, string_isSignal, string_GT, string_JECstring, string_jsonName, string_triggerTag, string_filterString]
            print config
        else:
            print "Invalid crab action. Please type: -a submit/status/resubmit/getoutput/kill"
            exit()
    os.system('echo -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-\n') 





