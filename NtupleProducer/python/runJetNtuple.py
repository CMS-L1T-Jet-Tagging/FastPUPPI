import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
import os

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1
inputMC = ['file:/eos/cms/store/cmst3/group/l1tr/FastPUPPI/15_1_X/fpinputs_140X/v0/GluGluHHTo2B2Tau_PU200/inputs140X_1-1.root']
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*inputMC),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*",
            "drop l1tKMTFTracks_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.pfMet_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun4_realistic_v1', '')

# NOTE: we need this to avoid saving the stubs
process.l1tTrackSelectionProducer.processSimulatedTracks = False

from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import l1tPhase2L1CaloEGammaEmulator
process.l1tPhase2L1CaloEGammaEmulator = l1tPhase2L1CaloEGammaEmulator.clone()

from L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi import l1tPhase2CaloPFClusterEmulator
process.l1tPhase2CaloPFClusterEmulator = l1tPhase2CaloPFClusterEmulator.clone()

from L1Trigger.L1CaloTrigger.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator_cfi import l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator = l1tPhase2GCTBarrelToCorrelatorLayer1Emulator.clone()

from L1Trigger.Phase2L1ParticleFlow.l1tMETPFProducer_cfi import l1tMETPFProducer
process.l1tMETPFProducer = l1tMETPFProducer.clone()

process.extraPFStuff = cms.Task(
        process.l1tPhase2L1CaloEGammaEmulator,
        process.l1tPhase2CaloPFClusterEmulator,
        process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator,
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tTrackSelectionProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TLayer2EGTask,
        process.l1tMETPFProducer)

def addJetNTuple(trktype = "extended", nparam = 5, tagged = True):
    # create new jet tupler
    jetColl = "l1tSC4PFL1PuppiExtendedEmulator"
    jetCollCorr = "l1tSC4PFL1PuppiExtendedEmulator"
    if trktype == "baseline":
        jetColl = "l1tSC4PFL1PuppiEmulator"
        jetCollCorr = "l1tSC4PFL1PuppiCorrectedEmulator"
        process.l1tPFTracksFromL1Tracks.nParam = cms.uint32(nparam)
    else:
        process.l1tPFTracksFromL1TracksExtended.nParam = cms.uint32(nparam)
    if tagged:
        jetColl = ("l1tSC4NGJetProducer","l1tSC4NGJets")
        jetCollCorr = "l1tSC4PFL1PuppiCorrectedEmulator"

    process.outnano = cms.EDAnalyzer("JetNTuplizer",
        genJets = cms.InputTag("ak4GenJetsNoNu"),
        genParticles = cms.InputTag("genParticles"),
        scPuppiJets = cms.InputTag(jetColl),
        scPuppiJetsCorr = cms.InputTag(jetCollCorr),
        nnTaus = cms.InputTag("l1tNNTauProducerPuppi","L1PFTausNN"),
        genJetsFlavour = cms.InputTag("genFlavourInfo"),
        vtx = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
        bjetIDs = cms.InputTag("l1tBJetProducerPuppiCorrectedEmulator", "L1PFBJets"),
        electrons = cms.InputTag("l1tLayer2EG","L1CtTkElectron"),
        muons = cms.InputTag("l1tSAMuonsGmt","promptSAMuons"),
    )
    process.endTuple = cms.EndPath(process.outnano)
    outName = "jetTuple_"+trktype+"_"+str(nparam)+".root"
    process.TFileService = cms.Service("TFileService", fileName = cms.string(outName))

# to check available tags:
process.p = cms.Path()
process.p.associate(process.extraPFStuff)
process.p.associate(process.L1TPFJetsExtendedTask)
process.p.associate(process.L1TBJetsTask)
#process.p.associate(process.l1tSC4NGJetTask)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetTuple.root"))

def addCalib():
    process.load("L1Trigger.Phase2L1ParticleFlow.l1tPFClustersFromHGC3DClustersEM_cfi")


    process.l1tLayer1BarrelRaw = process.l1tLayer1Barrel.clone(
        gctEmInputConversionParameters = process.l1tLayer1Barrel.gctEmInputConversionParameters.clone(
            gctEmCorrector = cms.string("")
        )
    )


    process.l1tLayer1HGCalRaw = process.l1tLayer1HGCal.clone(
        hgcalInputConversionParameters = process.l1tLayer1HGCal.hgcalInputConversionParameters.clone(
            corrector= cms.string("")
            )
        )
    process.l1tLayer1HGCalNoTKRaw = process.l1tLayer1HGCalNoTK.clone(
        hgcalInputConversionParameters = process.l1tLayer1HGCalNoTK.hgcalInputConversionParameters.clone(
            corrector= cms.string("")
            )
        )

    process.extraPFStuff.add(
        process.l1tLayer1BarrelRaw,
        process.l1tLayer1HGCalRaw,
        process.l1tLayer1HGCalNoTKRaw
    )

    #uncalibrated
    process.ntuple.objects.L1RawBarrelEcal   = cms.VInputTag('l1tLayer1BarrelRaw:DecodedEmClusters')
    process.ntuple.objects.L1RawBarrelCalo   = cms.VInputTag('l1tPFClustersFromCombinedCaloHCal:uncalibrated')
    process.ntuple.objects.L1RawHGCal   = cms.VInputTag('l1tLayer1HGCalRaw:DecodedHadClusters', 'l1tLayer1HGCalNoTKRaw:DecodedHadClusters')#use only this, try to understand if you have to use emf or emf_tot
    process.ntuple.objects.L1RawHGCalEM = cms.VInputTag('l1tLayer1HGCalRaw:DecodedEmClusters', 'l1tLayer1HGCalNoTKRaw:DecodedEmClusters')
    process.ntuple.objects.L1RawHFCalo  = cms.VInputTag('l1tPFClustersFromCombinedCaloHF:uncalibrated')

    #calibrated
    #process.ntuple.objects.L1BarrelEcal = cms.VInputTag('l1tPFClustersFromL1EGClusters' ) #Change IT. calibrated decodedEmClu in barrel not available
    #process.ntuple.objects.L1BarrelCalo = cms.VInputTag('l1tPFClustersFromCombinedCaloHCal:calibrated')
    #process.ntuple.objects.L1HGCal   = cms.VInputTag('l1tPFClustersFromHGC3DClusters')
    #process.ntuple.objects.L1HFCalo  = cms.VInputTag('l1tPFClustersFromCombinedCaloHF:calibrated')
    #process.ntuple.objects.L1HGCalEM = cms.VInputTag('l1tPFClustersFromHGC3DClustersEM', )


def addNNPuppiTaus():
    process.load("L1Trigger.Phase2L1ParticleFlow.L1NNTauProducer_cff")
    process.l1tNNTauProducerPuppi.maxtaus = cms.int32(500)
    process.extraPFStuff.add(process.l1tNNTauProducerPuppi)

def addSeededConeJets():
    process.extraPFStuff.add(process.L1TPFJetsTask)
    process.extraPFStuff.add(process.L1TPFJetsExtendedTask)

def addMultitagging(trktype = "extended"):
    if trktype == "extended":
        process.l1tSC4NGJetProducer.jets = cms.InputTag("l1tSC4PFL1PuppiExtendedEmulator")
    else:
        process.l1tSC4NGJetProducer.jets = cms.InputTag("l1tSC4PFL1PuppiEmulator")
    process.l1tSC4NGJetProducer.maxJets = cms.int32(500)
    process.l1tSC4NGJetProducer.l1tSC4NGJetModelPath = cms.string(os.environ['CMSSW_BASE']+"/src/L1TSC4NGJetModel/L1TSC4NGJetModel_v0")
    process.extraPFStuff.add(process.l1tSC4NGJetProducer)
    process.l1tSC4NGJetProducer.doJEC = cms.bool(True)
    process.l1tSC4NGJetProducer.correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root")
    process.l1tSC4NGJetProducer.correctorDir = cms.string("L1PuppiSC4EmuJets")

def addBtagging(jetColl): #extended TRK
    process.load("L1Trigger.Phase2L1ParticleFlow.L1BJetProducer_cff")
    process.l1tBJetProducerPuppiCorrectedEmulator.jets = cms.InputTag(jetColl)
    process.l1tBJetProducerPuppiCorrectedEmulator.maxJets = cms.int32(500)
    process.l1tBJetProducerPuppiCorrectedEmulator.useRawPt = cms.bool(True)
    process.extraPFStuff.add(process.L1TBJetsTask)
    #process.l1pfjetTable.jets.scPuppiBJet = cms.InputTag('l1tBJetProducerPuppiCorrectedEmulator')  

def addGenJetFlavourTable():
    process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.selectedHadronsAndPartons.partonMode = cms.string("Pythia8")
    process.genFlavourInfo = process.ak4JetFlavourInfos.clone(jets = "ak4GenJetsNoNu")
    process.p += process.selectedHadronsAndPartons
    process.p += process.genFlavourInfo

def goMT(nthreads=1):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)

if True:
    process.source.fileNames  = cms.untracked.vstring(*inputMC)
    goMT()
    trktype = "extended"
    nparam = 5
    addSeededConeJets()
    addMultitagging(trktype = trktype)
    addBtagging(("l1tSC4NGJetProducer","l1tSC4NGJets"))
    addNNPuppiTaus()
    addGenJetFlavourTable()
    addJetNTuple(trktype = trktype, nparam = nparam)
    if False:
        open("debug_dump_runJetNTuple.py", "w").write(process.dumpPython())