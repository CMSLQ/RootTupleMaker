
import FWCore.ParameterSet.Config as cms

process = cms.Process("treeCreator")

process.load("FWCore.MessageService.MessageLogger_cfi")

############## IMPORTANT ########################################
# if you run over many samples ans you save the log remember to reduce
# the size of the output by prescaling the report of the event number
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.default.limit = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#################################################################

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:/pcuscms/pcuscms46/data/MCData/AOD/PYTHIA6_Exotica_LQ_eejj_400_cff_py_GEN_FASTSIM.root'
    'file:/pcuscms/pcuscms46/data/MCData/AOD/PYTHIA6_Exotica_LQ_enujj_400_cff_py_GEN_FASTSIM.root'    
    #   'file:/home/lockner/Data/Summer08_Bkgnd/QCD_AOD.root'     # AOD (QCD)
    #'file:/home/santanas/Data/C81A2D83-ED9A-DD11-98F1-0015C5E59E7F.root'  #FULLSIM RECO (QCD)
    )
)

process.treeCreator = cms.EDAnalyzer('RootTupleMaker'
)

#process.treeCreator.rootfile        = cms.untracked.string("TTree_PYTHIA6_Exotica_LQ_enujj_400_cff_py_GEN_FASTSIM.root")
process.treeCreator.rootfile        = cms.untracked.string("TTree_test_enujj.root")
process.treeCreator.maxgenparticles = cms.untracked.int32(50)
process.treeCreator.maxgenjets      = cms.untracked.int32(10)
process.treeCreator.maxelectrons    = cms.untracked.int32(10)
process.treeCreator.maxcalojets     = cms.untracked.int32(10)
process.treeCreator.maxmuons        = cms.untracked.int32(10)
process.treeCreator.aodsim          = cms.untracked.bool(True)
process.treeCreator.fastSim         = cms.untracked.bool(True)
process.treeCreator.PAT             = cms.untracked.bool(False)
process.treeCreator.debug           = cms.untracked.bool(False)
# overall luminosity normalization  (in pb-1) 	
process.treeCreator.luminosity      =  cms.untracked.double(100)
process.treeCreator.numEvents       = cms.untracked.int32(10)
process.treeCreator.saveTrigger     = cms.untracked.bool(True)

process.treeCreator.useSkim1st2ndGenLQ = cms.untracked.bool(False)
process.treeCreator.useSkim1st2ndGenLQenujj = cms.untracked.bool(False)
process.treeCreator.usePDFweight       = cms.untracked.bool(False)
process.treeCreator.PDFSet             = cms.untracked.string("/cteq61.LHgrid")
process.treeCreator.skim1st2ndGenLQpTEle  =  cms.untracked.double(20)
process.treeCreator.skim1st2ndGenLQpTMu  =  cms.untracked.double(20)
process.treeCreator.skim1st2ndGenLQpTJet  =  cms.untracked.double(10)
process.treeCreator.skim1st2ndGenLQDeltaRJetEle  =  cms.untracked.double(0.1)
##only for Skim1st2ndGenLQenujj
process.treeCreator.skim1st2ndGenLQMET = cms.untracked.double(40)
##

######## electron isolation  ########
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("EgammaAnalysis.EgammaIsolationProducers.egammaElectronTkIsolation_cfi")
process.load("EgammaAnalysis.EgammaIsolationProducers.egammaElectronTkNumIsolation_cfi")
process.load("EgammaAnalysis.EgammaIsolationProducers.egammaEcalRecHitIsolation_cfi") #only in RECO
process.load("EgammaAnalysis.EgammaIsolationProducers.egammaHcalIsolation_cfi") #only in RECO
process.load("EgammaAnalysis.EgammaIsolationProducers.egammaTowerIsolation_cfi") #for AOD

process.egammaEcalRecHitIsolation.extRadius = cms.double(0.3)
process.egammaEcalRecHitIsolation.etMin = cms.double(0.)

#process.egammaElectronTkIsolation.trackProducer = cms.InputTag("gsWithMaterialTracks")
process.egammaElectronTkIsolation.ptMin = cms.double(1.5)
process.egammaElectronTkIsolation.intRadius = cms.double(0.02)
process.egammaElectronTkIsolation.extRadius = cms.double(0.2)
process.egammaElectronTkIsolation.maxVtxDist = cms.double(0.1)

process.egammaElectronTkNumIsolation.ptMin = cms.double(1.5)
process.egammaElectronTkNumIsolation.intRadius = cms.double(0.02)
process.egammaElectronTkNumIsolation.extRadius = cms.double(0.2)
process.egammaElectronTkNumIsolation.maxVtxDist = cms.double(0.1)

process.reducedEcalRecHitIsolation = cms.EDProducer("EgammaEcalRecHitIsolationProducer",
    absolut = cms.bool(True),
    ecalBarrelRecHitProducer = cms.InputTag("reducedEcalRecHitsEB"),
    ecalEndcapRecHitCollection = cms.InputTag(""),
    intRadius = cms.double(0.0),
    ecalEndcapRecHitProducer = cms.InputTag("reducedEcalRecHitsEE"),
    extRadius = cms.double(0.3),
    useIsolEt = cms.bool(True),
    ecalBarrelRecHitCollection = cms.InputTag(""),
    etMin = cms.double(0.0),
    emObjectProducer = cms.InputTag("pixelMatchGsfElectrons")
)

#############   Define the L2 correction service #####
process.L2RelativeJetCorrector = cms.ESSource("L2RelativeCorrectionService", 
    tagName = cms.string('Summer08_L2Relative_IC5Calo'),
    label = cms.string('L2RelativeJetCorrector')
)
#############   Define the L3 correction service #####
process.L3AbsoluteJetCorrector = cms.ESSource("L3AbsoluteCorrectionService", 
    tagName = cms.string('Summer08_L3Absolute_IC5Calo'),
    label = cms.string('L3AbsoluteJetCorrector')
)
#############   Define the L5 correction service #####
process.L5JetCorrector = cms.ESSource("L5FlavorCorrectionService",
#    section = cms.string('b'), 
#    tagName = cms.string('L5Flavor_fromTTbar_iterativeCone5'),
    section = cms.string('uds'), 
    tagName = cms.string('L5Flavor_fromQCD_iterativeCone5'),
    label = cms.string('L5FlavorJetCorrector')
)
#############   Define the chain corrector service ###
process.L2L3JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L2RelativeJetCorrector','L3AbsoluteJetCorrector'),
    label = cms.string('L2L3JetCorrector')
)
#############   Define the chain corrector module ####
process.L2L3CorJet = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("iterativeCone5CaloJets"),
    correctors = cms.vstring('L2L3JetCorrector')
)
## #############   Define the chain corrector service ###
process.L2L3L5JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L2RelativeJetCorrector','L3AbsoluteJetCorrector','L5FlavorJetCorrector'),
    label = cms.string('L2L3L5JetCorrector')
)
#############   Define the chain corrector module ####
process.L2L3L5CorJet = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("iterativeCone5CaloJets"),
    correctors = cms.vstring('L2L3L5JetCorrector')
)
# set the record's IOV. Must be defined once. Choose ANY correction service. #
##process.prefer("L2L3JetCorrector") 
##process.prefer("L2L3L5JetCorrector") 

##process.p = cms.Path(process.L2L3CorJet * process.treeCreator)
#process.p = cms.Path(process.L2L3L5CorJet * process.L2L3CorJet * process.egammaIsolationSequence * process.treeCreator)
process.p = cms.Path(process.L2L3L5CorJet * process.L2L3CorJet
                     * process.egammaElectronTkIsolation * process.egammaElectronTkNumIsolation
#                     * process.egammaEcalRecHitIsolation #for RECO
                     * process.reducedEcalRecHitIsolation #for both
#                     * process.egammaHcalIsolation #for RECO
                     * process.egammaTowerIsolation #for AOD
                     * process.treeCreator)

