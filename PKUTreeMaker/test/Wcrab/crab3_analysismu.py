from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'SMu15D-vA1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt','Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt','Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt','Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'analysis_data.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 14
config.Data.lumiMask = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
#config.Data.runRange = '246908-258750'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'SMu15D-vA1'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'



