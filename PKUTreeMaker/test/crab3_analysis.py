from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'WAmatchingnew_50ns_10_12_1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt','Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt','Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt','Summer15_25nsV5_DATA_L1FastJet_AK4PFchs.txt','Summer15_25nsV5_DATA_L2Relative_AK4PFchs.txt','Summer15_25nsV5_DATA_L3Absolute_AK4PFchs.txt','Summer15_25nsV5_DATA_L2L3Residual_AK4PFchs.txt']
config.JobType.psetName    = 'analysis.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 7000

config.section_("Data")
config.Data.inputDataset = '/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =1
config.Data.totalUnits = -1
config.Data.publication = False

# This string is used to construct the output dataset name
config.Data.publishDataName = 'WAmatchingnew_25ns'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
