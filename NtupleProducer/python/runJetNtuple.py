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

# for offline btagging
from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
from RecoBTag.FeatureTools.ParticleNetFeatureEvaluator_cfi import ParticleNetFeatureEvaluator
from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4DiscriminatorsJetTags_cfi import *
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

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
    # run NGJets as default
    if tagged:
        jetColl = ("l1tSC4NGJetProducer","l1tSC4NGJets")
        jetCollCorr = "l1tSC4PFL1PuppiCorrectedEmulator"
        if trktype == "baseline":
            process.l1tSC4NGJetProducer.jets = cms.InputTag("l1tSC4PFL1PuppiEmulator")
        elif trktype == "extended":
            process.l1tSC4NGJetProducer.jets = cms.InputTag("l1tSC4PFL1PuppiExtendedEmulator")

    process.l1tSC4NGJetProducer.returnRawPt = cms.bool(True)

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
        offlineJets = cms.InputTag("slimmedJetsPuppi"),
        offlinePVs = cms.InputTag("offlineSlimmedPrimaryVertices"),
    )

    if hasattr(process, "slimmedJetsUpdated"):
        process.outnano.offlineJets = cms.InputTag("slimmedJetsUpdated")

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
    #process.l1tSC4NGJetProducer.l1tSC4NGJetModelPath = cms.string(os.environ['CMSSW_BASE']+"/src/L1TSC4NGJetModel/L1TSC4NGJetModel_v0")
    process.extraPFStuff.add(process.l1tSC4NGJetProducer)
    process.l1tSC4NGJetProducer.doJEC = cms.bool(True)
    process.l1tSC4NGJetProducer.correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root")
    process.l1tSC4NGJetProducer.correctorDir = cms.string("L1PuppiSC4EmuJets")
    process.l1tSC4NGJetProducer.minPt = cms.double(0.)
    process.l1tSC4NGJetProducer.maxEta = cms.double(99.)
    process.l1tSC4NGJetProducer.maxJets = cms.int32(999)

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

def addOfflineBTagging():
    process.pfParticleNetFromMiniAODAK4PuppiCentralTagInfos = ParticleNetFeatureEvaluator.clone(
        jets = "slimmedJetsPuppi",
        jet_radius = 0.4,
        min_jet_pt = 0.,
        min_jet_eta = 0.,
        max_jet_eta = 2.5,
    )

    process.pfParticleNetFromMiniAODAK4PuppiForwardTagInfos = ParticleNetFeatureEvaluator.clone(
        jets = "slimmedJetsPuppi",
        jet_radius = 0.4,
        min_jet_pt = 0.,
        min_jet_eta = 2.5,
        max_jet_eta = 5.0,
    )

    process.pfParticleNetFromMiniAODAK4PuppiCentralTagInfos.puppi_value_map = cms.InputTag("") # edm::ValueMap<float> "puppi" doesn't exist in our miniAOD, let's try without it
    process.pfParticleNetFromMiniAODAK4PuppiForwardTagInfos.puppi_value_map = cms.InputTag("") # edm::ValueMap<float> "puppi" doesn't exist in our miniAOD, let's try without it
    process.pfParticleNetFromMiniAODAK4PuppiCentralTagInfos.fallback_puppi_weight = cms.bool(True)
    process.pfParticleNetFromMiniAODAK4PuppiForwardTagInfos.fallback_puppi_weight = cms.bool(True)

    process.pfParticleNetFromMiniAODAK4PuppiCentralJetTags = boostedJetONNXJetTagsProducer.clone(
        src = 'pfParticleNetFromMiniAODAK4PuppiCentralTagInfos',
        preprocess_json = 'RecoBTag/Combined/data/ParticleNetFromMiniAODAK4/PUPPI/Central/preprocess.json',
        model_path = 'RecoBTag/Combined/data/ParticleNetFromMiniAODAK4/PUPPI/Central/modelfile/model.onnx',
        flav_names = ['probmu','probele','probtaup1h0p','probtaup1h1p','probtaup1h2p','probtaup3h0p','probtaup3h1p','probtaum1h0p','probtaum1h1p','probtaum1h2p','probtaum3h0p','probtaum3h1p','probb','probc','probuds','probg','ptcorr','ptreshigh','ptreslow','ptnu'],
    )
    process.pfParticleNetFromMiniAODAK4PuppiForwardJetTags = boostedJetONNXJetTagsProducer.clone(
        src = 'pfParticleNetFromMiniAODAK4PuppiForwardTagInfos',
        preprocess_json = 'RecoBTag/Combined/data/ParticleNetFromMiniAODAK4/PUPPI/Forward/preprocess.json',
        model_path = 'RecoBTag/Combined/data/ParticleNetFromMiniAODAK4/PUPPI/Forward/modelfile/model.onnx',
        flav_names = ['probq','probg','ptcorr','ptreshigh','ptreslow','ptnu'],
    )


    from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4DiscriminatorsJetTags_cfi import pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags
    from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4DiscriminatorsJetTags_cfi import pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags

    process.pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags = pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags.clone()
    process.pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags = pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags.clone()

    process.pfParticleNetFromMiniAODAK4PuppiTask = cms.Task(process.pfParticleNetFromMiniAODAK4PuppiCentralTagInfos,
                                                    process.pfParticleNetFromMiniAODAK4PuppiForwardTagInfos,
                                                    process.pfParticleNetFromMiniAODAK4PuppiCentralJetTags,
                                                    process.pfParticleNetFromMiniAODAK4PuppiForwardJetTags,
                                                    process.pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags,
                                                    process.pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags
                                                    )

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
    process.slimmedJetsUpdated = updatedPatJets.clone(
        jetSource = "slimmedJetsPuppi",
        addJetCorrFactors = False,
        discriminatorSources = cms.VInputTag(
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probb"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probc"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probuds"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probg"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probmu"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probele"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h0p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h1p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h2p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup3h0p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup3h1p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h0p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h1p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h2p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum3h0p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum3h1p"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptcorr"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptnu"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreshigh"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreslow"),

            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsAll"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsL"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:CvsB"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:TauVsJet"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:TauVsEle"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:TauVsMu"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:QvsG"),

            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:probq"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:probg"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptcorr"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptnu"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptreshigh"),
            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardJetTags:ptreslow"),

            cms.InputTag("pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags:QvsG"),

        )
    )
    process.slimmedJetsUpdatedTask = cms.Task(process.slimmedJetsUpdated)
    process.p.associate(process.pfParticleNetFromMiniAODAK4PuppiTask)
    process.p.associate(process.slimmedJetsUpdatedTask)

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
    addOfflineBTagging()
    addJetNTuple(trktype = trktype, nparam = nparam)
    if False:
        open("debug_dump_runJetNTuple.py", "w").write(process.dumpPython())