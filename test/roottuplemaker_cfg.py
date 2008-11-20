
import FWCore.ParameterSet.Config as cms

process = cms.Process("treeCreator")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'file:/afs/cern.ch/user/l/lockner/scratch0/CMSSW_2_1_8/src/data/LQ300_HLT.root'
       # 'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_0.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_1.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_2.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_3.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_4.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_5.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_6.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_7.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_8.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_9.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_10.root'
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_11.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_12.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_13.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_14.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_15.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_16.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_17.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_18.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_19.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_20.root'
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_21.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_22.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_23.root',
       #'file:/data/groups/LQ/FastSim/Leptoquarks/2_1_6/300/LQ_300_HLT_080916_59_24.root'
        #'file:/data/users/twedt/Leptoquark/PAT/LQ_test_full.root'
    )
)

process.treeCreator = cms.EDAnalyzer('RootTupleMaker'
)

process.treeCreator.rootfile        = cms.untracked.string("TTree.root")
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
process.treeCreator.numEvents       = cms.untracked.int32(100)
process.treeCreator.saveTrigger     = cms.untracked.bool(True)
process.treeCreator.UDSQuarksCorrector = cms.string("L5FlavorJetCorrectorUds")

######## electron isolation  ########
process.load( "RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff") ## import *
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

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
    section = cms.string('b'), 
    tagName = cms.string('L5Flavor_fromTTbar_iterativeCone5'),
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
process.p = cms.Path(process.L2L3L5CorJet * process.L2L3CorJet * process.treeCreator)
