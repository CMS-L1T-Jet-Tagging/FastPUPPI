// -*- C++ -*-
//
// Package:    Giovanni/NTuplizer
// Class:      NTuplizer
// 
/**\class NTuplizer NTuplizer.cc Giovanni/NTuplizer/plugins/NTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giovanni Petrucciani
//         Created:  Thu, 01 Sep 2016 11:30:38 GMT
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TParticleFlow/interface/PFTau.h"
#include "DataFormats/L1TParticleFlow/interface/jets.h"

#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkElectronFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkEm.h"
#include "DataFormats/L1TCorrelator/interface/TkEmFwd.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/L1TMuonPhase2/interface/SAMuon.h"

#include <cstdint>
#include <TTree.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/JetId.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TSC4NGJetID.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TParticleFlow/interface/datatypes.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/jetmet/L1SeedConePFJetEmulator.h"

// some tools to fix inputs or calculate them
namespace jettools{

    bool useoffsets = true;
    const float& catchInfs(const float& in,const float& replace_value){
            if(in==in){
                if(std::isinf(in))
                    return replace_value;
                else if(in < -1e32 || in > 1e32)
                    return replace_value;
                return in;
            }
            return replace_value;
    }


    float catchInfsAndBound(const float& in,const float& replace_value,
            const float& lowerbound, const float& upperbound,const float offset=0){
        float withoutinfs=catchInfs(in,replace_value);
        if(withoutinfs+offset<lowerbound) return lowerbound;
        if(withoutinfs+offset>upperbound) return upperbound;
        if(useoffsets)
            withoutinfs+=offset;
        return withoutinfs;
    }


    double etaRel(const math::XYZVector& dir, const math::XYZVector& track) {
        double momPar = dir.Dot(track);
        // double energy = std::sqrt(track.Mag2() + ROOT::Math::Square(reco::ParticleMasses::piPlus));
        double energy = std::sqrt(track.Mag2() + ROOT::Math::Square(0.13957));

        return 0.5 * std::log((energy + momPar) / (energy - momPar));
    }


    // Sorters to order object collections in decreasing order of pT
    template<typename T> 
    class PatPtSorter {
    public:
        bool operator()(const T& i, const T& j) const {
        return (i.pt() > j.pt());
        }
    };

    PatPtSorter<l1t::PFCandidate>  l1PFCandidateSorter;

    template<typename T> 
    class PatRefPtSorter {
    public:
        bool operator()(const T& i, const T& j) const {
        return (i->pt() > j->pt());
        }
    };
    PatRefPtSorter<l1t::PFJetRef>    jetRefSorter;
    PatRefPtSorter<l1t::PFTauRef>    tauRefSorter;
    PatRefPtSorter<l1t::TkElectronRef>    eleRefSorter;
    PatRefPtSorter<l1t::SAMuonRef>    muonRefSorter;
    PatRefPtSorter<reco::GenJetRef>  genJetRefSorter;

}

using namespace jettools;

class JetNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
    public:
        explicit JetNTuplizer(const edm::ParameterSet&);
        ~JetNTuplizer();

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {
            // bZ_ = 3.8112; // avoid loading the event setup
        }
        virtual void endRun(edm::Run const&, edm::EventSetup const& iSetup) override { } // framework wants this to be implemented

        void fill_genParticles(const edm::Event& iEvent);

        template<unsigned int bits=10>
        static float zip(float f) {
            return MiniFloatConverter::reduceMantissaToNbitsRounding<bits>(f);
        }

        edm::EDGetTokenT<std::vector<reco::GenJet>> genjets_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
        edm::EDGetTokenT<std::vector<l1t::PFJet>> scjets_;
        edm::EDGetTokenT<std::vector<l1t::PFJet>> scngjets_;
        edm::EDGetTokenT<std::vector<l1t::PFJet>> scjetsCorr_;
        edm::EDGetTokenT<std::vector<l1t::PFTau>> nntaus_;
        edm::EDGetTokenT<std::vector<l1t::TkElectron>> electrons_;
        edm::EDGetTokenT<std::vector<l1t::SAMuon>> muons_;
        edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetsFlavour_;
        edm::EDGetTokenT<std::vector<l1t::VertexWord>> const fVtxEmu_;
        edm::EDGetTokenT<edm::ValueMap<float>> const bjetids_;
        // const edm::InputTag pileupInfoTag_;
        TTree *tree_;
        uint32_t run_, lumi_; uint64_t event_;
         //   float bZ_;

        // float dRJetGenMatch_ = 0.2;
        float dRJetGenMatch_ = 0.4;
        bool  isMC_;
        float jetPtMin_ = 10;
        float jetEtaMin_ = -2.4;
        float jetEtaMax_ = 2.4;
        float jetPFCandidatePtMin_ = 0.;

        static constexpr size_t max_pfcand_ = 16;
        bool applyJetIDFlag = false;
        /// thresholds for matching
        float dRCone        = 0.2;
        float ptGenLeptonMin = 8;
        float ptGenTauVisibleMin = 15;



    // --------------------
    std::vector<reco::GenParticle> gToBB_;
    std::vector<reco::GenParticle> gToCC_;
    std::vector<reco::GenParticle> neutrinosLepB_;
    std::vector<reco::GenParticle> neutrinosLepB_C_;
    std::vector<reco::GenParticle> alltaus_;
    std::vector<reco::GenParticle> Bhadron_;
    std::vector<reco::GenParticle> Bhadron_daughter_;

    // --------------------
    // from PNET
    // Generator-level information (GEN particles)
    std::vector<float>  gen_particle_pt;
    std::vector<float>  gen_particle_eta;
    std::vector<float>  gen_particle_phi;
    std::vector<float>  gen_particle_mass;
    std::vector<int>    gen_particle_id;
    std::vector<unsigned int>  gen_particle_status;
    std::vector<int>    gen_particle_daughters_id;
    std::vector<unsigned int> gen_particle_daughters_igen;
    std::vector<unsigned int> gen_particle_daughters_status;
    std::vector<float>  gen_particle_daughters_pt;
    std::vector<float>  gen_particle_daughters_eta;
    std::vector<float>  gen_particle_daughters_phi;
    std::vector<float>  gen_particle_daughters_mass;
    std::vector<int>    gen_particle_daughters_charge;

    // Gen leptons from resonance decay 
    std::vector<TLorentzVector> genLepFromResonance4V_;
    std::vector<TLorentzVector> genMuonsFromResonance4V_;
    std::vector<TLorentzVector> genElectronsFromResonance4V_;
    std::vector<TLorentzVector> tau_gen_visible_;
    std::vector<TLorentzVector> tau_gen_;
    std::vector<int> tau_gen_charge_;
    std::vector<unsigned int> tau_gen_nch_;
    std::vector<unsigned int> tau_gen_np0_;
    std::vector<unsigned int> tau_gen_nnh_;

    int jet_tauflav_, jet_muflav_, jet_elflav_, jet_taudecaymode_, jet_lepflav_;
    int jet_taucharge_;
    float jet_genmatch_lep_pt_;
    float jet_genmatch_lep_vis_pt_;
    float jet_genmatch_lep_dR_;
    // --------------------
    bool jet_reject_;
    float jet_eta_;
    float jet_eta_phys_;
    float jet_phi_;
    float jet_phi_phys_;
    float jet_pt_;
    float jet_pt_phys_;
    float jet_pt_raw_;
    float jet_pt_corr_;
    float jet_mass_;
    float jet_energy_;

    // float jet_px_;
    // float jet_py_;
    // float jet_pz_;

    float jet_SC4NGJet_score_light_;
    float jet_SC4NGJet_score_b_;
    float jet_SC4NGJet_score_taup_;
    float jet_SC4NGJet_score_taum_;
    float jet_SC4NGJet_score_gluon_;
    float jet_SC4NGJet_score_charm_;
    float jet_SC4NGJet_score_electron_;
    float jet_SC4NGJet_score_muon_;
    float jet_SC4NGJet_score_regression_;

    float jet_bjetscore_;
    float jet_tauscore_;
    float jet_eletkiso_;
    float jet_elepfiso_;
    float jet_elepuppiiso_;
    float jet_muoniso_;
    float jet_taupt_;

    float jet_taumatch_dR_;
    float jet_elematch_dR_;
    float jet_muonmatch_dR_;

    float jet_jecmatch_dR_;

    unsigned int jet_hflav_;
    int jet_pflav_;


    // --------------------
    float jet_genmatch_pt_;
    float jet_genmatch_eta_;
    float jet_genmatch_phi_;
    float jet_genmatch_mass_;
    float jet_genmatch_dR_;
    unsigned int jet_genmatch_hflav_;
    int jet_genmatch_pflav_;


    // jet pf candidates
    unsigned int njet_pfcand_;
    std::vector<float> jet_pfcand_pt;
    std::vector<float> jet_pfcand_pt_phys;
    std::vector<float> jet_pfcand_pt_rel;
    std::vector<float> jet_pfcand_pt_rel_phys;
    std::vector<bool> jet_pfcand_isfilled;
    // std::vector<float> jet_pfcand_px;
    // std::vector<float> jet_pfcand_py;
    // std::vector<float> jet_pfcand_pz;
    std::vector<float> jet_pfcand_pt_log;
    std::vector<float> jet_pfcand_pt_rel_log;
    std::vector<float> jet_pfcand_eta;
    std::vector<float> jet_pfcand_eta_phys;
    std::vector<float> jet_pfcand_phi;
    std::vector<float> jet_pfcand_phi_phys;
    std::vector<float> jet_pfcand_mass;
    std::vector<float> jet_pfcand_energy;
    std::vector<float> jet_pfcand_energy_log;
    std::vector<float> jet_pfcand_puppiweight;
    std::vector<float> jet_pfcand_emid;
    std::vector<float> jet_pfcand_quality;
    std::vector<float> jet_pfcand_tkquality;
    std::vector<float> jet_pfcand_z0;
    std::vector<float> jet_pfcand_dxy;
    std::vector<float> jet_pfcand_dxy_custom;
    std::vector<unsigned int> jet_pfcand_id;
    std::vector<bool> jet_pfcand_isPhoton;
    std::vector<bool> jet_pfcand_isElectron;
    std::vector<bool> jet_pfcand_isElectronPlus;
    std::vector<bool> jet_pfcand_isElectronMinus;
    std::vector<bool> jet_pfcand_isMuon;
    std::vector<bool> jet_pfcand_isMuonPlus;
    std::vector<bool> jet_pfcand_isMuonMinus;
    std::vector<bool> jet_pfcand_isNeutralHadron;
    std::vector<bool> jet_pfcand_isChargedHadron;
    std::vector<bool> jet_pfcand_isChargedHadronPlus;
    std::vector<bool> jet_pfcand_isChargedHadronMinus;
    std::vector<int> jet_pfcand_charge;
    std::vector<float> jet_pfcand_pperp_ratio;
    std::vector<float> jet_pfcand_ppara_ratio;
    std::vector<float> jet_pfcand_deta;
    std::vector<float> jet_pfcand_deta_phys;
    std::vector<float> jet_pfcand_dphi;
    std::vector<float> jet_pfcand_dphi_phys;
    std::vector<float> jet_pfcand_etarel;

    std::vector<float> jet_pfcand_track_valid;
    std::vector<bool> jet_pfcand_track_isfilled;
    std::vector<float> jet_pfcand_track_rinv;
    std::vector<float> jet_pfcand_track_phizero;
    std::vector<float> jet_pfcand_track_tanl;
    std::vector<float> jet_pfcand_track_z0;
    std::vector<float> jet_pfcand_track_d0;
    std::vector<float> jet_pfcand_track_chi2rphi;
    std::vector<float> jet_pfcand_track_chi2rz;
    std::vector<float> jet_pfcand_track_bendchi2;
    std::vector<float> jet_pfcand_track_hitpattern;
    std::vector<float> jet_pfcand_track_nstubs;
    std::vector<float> jet_pfcand_track_mvaquality;
    std::vector<float> jet_pfcand_track_mvaother;

    std::vector<float> jet_pfcand_track_chi2;
    std::vector<float> jet_pfcand_track_chi2norm;
    std::vector<float> jet_pfcand_track_qual;
    std::vector<float> jet_pfcand_track_npar;
    // std::vector<float> jet_pfcand_track_nstubs;
    std::vector<float> jet_pfcand_track_vx;
    std::vector<float> jet_pfcand_track_vy;
    std::vector<float> jet_pfcand_track_vz;
    std::vector<float> jet_pfcand_track_pterror;

    std::vector<bool> jet_pfcand_cluster_isfilled;
    std::vector<float> jet_pfcand_cluster_hovere;
    std::vector<float> jet_pfcand_cluster_sigmarr;
    std::vector<float> jet_pfcand_cluster_abszbarycenter;
    std::vector<float> jet_pfcand_cluster_emet;
    std::vector<float> jet_pfcand_cluster_egvspion;
    std::vector<float> jet_pfcand_cluster_egvspu;





       
};

JetNTuplizer::JetNTuplizer(const edm::ParameterSet& iConfig) :
    genjets_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    scjets_(consumes<std::vector<l1t::PFJet>>(iConfig.getParameter<edm::InputTag>("scPuppiJets"))), // l1tSCPFL1PuppiEmulator
    scngjets_(consumes<std::vector<l1t::PFJet>>(iConfig.getParameter<edm::InputTag>("scNGPuppiJets"))), // l1tSC4NGJetProducerPuppi
    scjetsCorr_(consumes<std::vector<l1t::PFJet>>(iConfig.getParameter<edm::InputTag>("scPuppiJetsCorr"))), // l1tSCPFL1PuppiEmulator
    nntaus_(consumes<std::vector<l1t::PFTau>>(iConfig.getParameter<edm::InputTag>("nnTaus"))), 
    electrons_(consumes<std::vector<l1t::TkElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))), 
    muons_(consumes<std::vector<l1t::SAMuon>>(iConfig.getParameter<edm::InputTag>("muons"))), 
    genJetsFlavour_   (consumes<reco::JetFlavourInfoMatchingCollection >    (iConfig.getParameter<edm::InputTag>("genJetsFlavour"))),
    fVtxEmu_(consumes<std::vector<l1t::VertexWord>>(iConfig.getParameter<edm::InputTag>("vtx"))),
    bjetids_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("bjetIDs")))
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("Jets","Jets");
    // tree_->Branch("run",  &run_,  "run/i");
    // tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    tree_->Branch("jet_eta", &jet_eta_);
    tree_->Branch("jet_eta_phys", &jet_eta_phys_);
    tree_->Branch("jet_phi", &jet_phi_);
    tree_->Branch("jet_phi_phys", &jet_phi_phys_);
    tree_->Branch("jet_pt", &jet_pt_);
    tree_->Branch("jet_pt_phys", &jet_pt_phys_);
    tree_->Branch("jet_pt_raw", &jet_pt_raw_);
    tree_->Branch("jet_mass", &jet_mass_);
    tree_->Branch("jet_energy", &jet_energy_);
    // tree_->Branch("jet_px", &jet_px_);
    // tree_->Branch("jet_py", &jet_py_);
    // tree_->Branch("jet_pz", &jet_pz_);

    tree_->Branch("jet_SC4NGJet_score_b", &jet_SC4NGJet_score_b_);
    tree_->Branch("jet_SC4NGJet_score_light", &jet_SC4NGJet_score_light_);
    tree_->Branch("jet_SC4NGJet_score_gluon", &jet_SC4NGJet_score_gluon_);
    tree_->Branch("jet_SC4NGJet_score_charm", &jet_SC4NGJet_score_charm_);
    tree_->Branch("jet_SC4NGJet_score_taup", &jet_SC4NGJet_score_taup_);
    tree_->Branch("jet_SC4NGJet_score_taum", &jet_SC4NGJet_score_taum_);
    tree_->Branch("jet_SC4NGJet_score_muon", &jet_SC4NGJet_score_muon_);
    tree_->Branch("jet_SC4NGJet_score_electron", &jet_SC4NGJet_score_electron_);
    tree_->Branch("jet_SC4NGJet_score_regression", &jet_SC4NGJet_score_regression_);

    tree_->Branch("jet_bjetscore", &jet_bjetscore_);
    tree_->Branch("jet_tauscore", &jet_tauscore_);
    tree_->Branch("jet_eletkiso", &jet_eletkiso_);
    tree_->Branch("jet_elepfiso", &jet_elepfiso_);
    tree_->Branch("jet_elepuppiiso", &jet_elepuppiiso_);
    tree_->Branch("jet_muoniso", &jet_muoniso_);
    tree_->Branch("jet_taupt", &jet_taupt_);
    
    tree_->Branch("jet_reject", &jet_reject_);
    tree_->Branch("jet_tauflav", &jet_tauflav_);
    tree_->Branch("jet_muflav", &jet_muflav_);
    tree_->Branch("jet_elflav", &jet_elflav_);
    tree_->Branch("jet_taudecaymode", &jet_taudecaymode_);
    tree_->Branch("jet_lepflav", &jet_lepflav_);
    tree_->Branch("jet_taucharge", &jet_taucharge_);

    tree_->Branch("jet_genmatch_lep_pt", &jet_genmatch_lep_pt_);
    tree_->Branch("jet_genmatch_lep_vis_pt", &jet_genmatch_lep_vis_pt_);
    tree_->Branch("jet_genmatch_lep_dR", &jet_genmatch_lep_dR_);

    tree_->Branch("jet_genmatch_pt", &jet_genmatch_pt_);
    tree_->Branch("jet_genmatch_eta", &jet_genmatch_eta_);
    tree_->Branch("jet_genmatch_phi", &jet_genmatch_phi_);
    tree_->Branch("jet_genmatch_mass", &jet_genmatch_mass_);
    tree_->Branch("jet_genmatch_dR", &jet_genmatch_dR_);
    tree_->Branch("jet_taumatch_dR", &jet_taumatch_dR_);
    tree_->Branch("jet_elematch_dR", &jet_elematch_dR_);
    tree_->Branch("jet_muonmatch_dR", &jet_muonmatch_dR_);
    tree_->Branch("jet_genmatch_hflav", &jet_genmatch_hflav_);
    tree_->Branch("jet_genmatch_pflav", &jet_genmatch_pflav_);

    tree_->Branch("jet_jecmatch_dR", &jet_jecmatch_dR_);
    tree_->Branch("jet_pt_corr", &jet_pt_corr_);

    tree_->Branch("jet_npfcand", &njet_pfcand_);
    tree_->Branch("jet_pfcand_isfilled", &jet_pfcand_isfilled, njet_pfcand_);
    tree_->Branch("jet_pfcand_pt", &jet_pfcand_pt, njet_pfcand_);
    tree_->Branch("jet_pfcand_pt_phys", &jet_pfcand_pt_phys, njet_pfcand_);
    tree_->Branch("jet_pfcand_pt_rel", &jet_pfcand_pt_rel, njet_pfcand_);
    tree_->Branch("jet_pfcand_pt_rel_phys", &jet_pfcand_pt_rel_phys, njet_pfcand_);
    // tree_->Branch("jet_pfcand_px", &jet_pfcand_px, njet_pfcand_);
    // tree_->Branch("jet_pfcand_py", &jet_pfcand_py, njet_pfcand_);
    // tree_->Branch("jet_pfcand_pz", &jet_pfcand_pz, njet_pfcand_);
    tree_->Branch("jet_pfcand_pt_log", &jet_pfcand_pt_log, njet_pfcand_);
    tree_->Branch("jet_pfcand_pt_rel_log", &jet_pfcand_pt_rel_log, njet_pfcand_);
    tree_->Branch("jet_pfcand_eta", &jet_pfcand_eta, njet_pfcand_);
    tree_->Branch("jet_pfcand_eta_phys", &jet_pfcand_eta_phys, njet_pfcand_);
    tree_->Branch("jet_pfcand_phi", &jet_pfcand_phi, njet_pfcand_);
    tree_->Branch("jet_pfcand_phi_phys", &jet_pfcand_phi_phys, njet_pfcand_);
    tree_->Branch("jet_pfcand_mass", &jet_pfcand_mass, njet_pfcand_);
    tree_->Branch("jet_pfcand_energy", &jet_pfcand_energy, njet_pfcand_);
    tree_->Branch("jet_pfcand_energy_log", &jet_pfcand_energy_log, njet_pfcand_);
    tree_->Branch("jet_pfcand_puppiweight", &jet_pfcand_puppiweight, njet_pfcand_);
    tree_->Branch("jet_pfcand_emid", &jet_pfcand_emid, njet_pfcand_);
    tree_->Branch("jet_pfcand_quality", &jet_pfcand_quality, njet_pfcand_);
    tree_->Branch("jet_pfcand_tkquality", &jet_pfcand_tkquality, njet_pfcand_);
    tree_->Branch("jet_pfcand_z0", &jet_pfcand_z0, njet_pfcand_);
    tree_->Branch("jet_pfcand_dxy", &jet_pfcand_dxy, njet_pfcand_);
    tree_->Branch("jet_pfcand_dxy_custom", &jet_pfcand_dxy_custom, njet_pfcand_);
    tree_->Branch("jet_pfcand_id", &jet_pfcand_id, njet_pfcand_);
    tree_->Branch("jet_pfcand_isPhoton", &jet_pfcand_isPhoton, njet_pfcand_);
    tree_->Branch("jet_pfcand_isElectron", &jet_pfcand_isElectron, njet_pfcand_);
    tree_->Branch("jet_pfcand_isElectronPlus", &jet_pfcand_isElectronPlus, njet_pfcand_);
    tree_->Branch("jet_pfcand_isElectronMinus", &jet_pfcand_isElectronMinus, njet_pfcand_);
    tree_->Branch("jet_pfcand_isMuon", &jet_pfcand_isMuon, njet_pfcand_);
    tree_->Branch("jet_pfcand_isMuonPlus", &jet_pfcand_isMuonPlus, njet_pfcand_);
    tree_->Branch("jet_pfcand_isMuonMinus", &jet_pfcand_isMuonMinus, njet_pfcand_);
    tree_->Branch("jet_pfcand_isNeutralHadron", &jet_pfcand_isNeutralHadron, njet_pfcand_);
    tree_->Branch("jet_pfcand_isChargedHadron", &jet_pfcand_isChargedHadron, njet_pfcand_);
    tree_->Branch("jet_pfcand_isChargedHadronPlus", &jet_pfcand_isChargedHadronPlus, njet_pfcand_);
    tree_->Branch("jet_pfcand_isChargedHadronMinus", &jet_pfcand_isChargedHadronMinus, njet_pfcand_);

    tree_->Branch("jet_pfcand_charge", &jet_pfcand_charge, njet_pfcand_);
    tree_->Branch("jet_pfcand_pperp_ratio", &jet_pfcand_pperp_ratio, njet_pfcand_);
    tree_->Branch("jet_pfcand_ppara_ratio", &jet_pfcand_ppara_ratio, njet_pfcand_);
    tree_->Branch("jet_pfcand_deta", &jet_pfcand_deta, njet_pfcand_);
    tree_->Branch("jet_pfcand_deta_phys", &jet_pfcand_deta_phys, njet_pfcand_);
    tree_->Branch("jet_pfcand_dphi", &jet_pfcand_dphi, njet_pfcand_);
    tree_->Branch("jet_pfcand_dphi_phys", &jet_pfcand_dphi_phys, njet_pfcand_);
    tree_->Branch("jet_pfcand_etarel", &jet_pfcand_etarel, njet_pfcand_);

    tree_->Branch("jet_pfcand_track_valid", &jet_pfcand_track_valid, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_isfilled", &jet_pfcand_track_isfilled, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_rinv", &jet_pfcand_track_rinv, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_phizero", &jet_pfcand_track_phizero, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_tanl", &jet_pfcand_track_tanl, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_z0", &jet_pfcand_track_z0, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_d0", &jet_pfcand_track_d0, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_chi2rphi", &jet_pfcand_track_chi2rphi, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_chi2rz", &jet_pfcand_track_chi2rz, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_bendchi2", &jet_pfcand_track_bendchi2, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_hitpattern", &jet_pfcand_track_hitpattern, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_nstubs", &jet_pfcand_track_nstubs, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_mvaquality", &jet_pfcand_track_mvaquality, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_mvaother", &jet_pfcand_track_mvaother, njet_pfcand_);

    tree_->Branch("jet_pfcand_track_chi2", &jet_pfcand_track_chi2, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_chi2norm", &jet_pfcand_track_chi2norm, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_qual", &jet_pfcand_track_qual, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_npar", &jet_pfcand_track_npar, njet_pfcand_);
    // tree_->Branch("jet_pfcand_track_nstubs", &jet_pfcand_track_nstubs, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_vx", &jet_pfcand_track_vx, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_vy", &jet_pfcand_track_vy, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_vz", &jet_pfcand_track_vz, njet_pfcand_);
    tree_->Branch("jet_pfcand_track_pterror", &jet_pfcand_track_pterror, njet_pfcand_);


    tree_->Branch("jet_pfcand_cluster_isfilled", &jet_pfcand_cluster_isfilled, njet_pfcand_);
    tree_->Branch("jet_pfcand_cluster_hovere", &jet_pfcand_cluster_hovere, njet_pfcand_);
    tree_->Branch("jet_pfcand_cluster_sigmarr", &jet_pfcand_cluster_sigmarr, njet_pfcand_);
    tree_->Branch("jet_pfcand_cluster_abszbarycenter", &jet_pfcand_cluster_abszbarycenter, njet_pfcand_);
    tree_->Branch("jet_pfcand_cluster_emet", &jet_pfcand_cluster_emet, njet_pfcand_);
    tree_->Branch("jet_pfcand_cluster_egvspion", &jet_pfcand_cluster_egvspion, njet_pfcand_);
    tree_->Branch("jet_pfcand_cluster_egvspu", &jet_pfcand_cluster_egvspu, njet_pfcand_);


    // -------------------------------------
    // settings for output TFile and TTree
    fs->file().SetCompressionAlgorithm(ROOT::ECompressionAlgorithm::kLZ4);
    fs->file().SetCompressionLevel(4);
    for (int idx = 0; idx < tree_->GetListOfBranches()->GetEntries(); ++idx) {
        TBranch* br = dynamic_cast<TBranch*>(tree_->GetListOfBranches()->At(idx));
        if (br) {
            br->SetBasketSize(1024 * 1024);
        }
    }
    if (tree_->GetListOfBranches()->GetEntries() > 0) {
        tree_->SetAutoFlush(-1024 * 1024 * tree_->GetListOfBranches()->GetEntries());
    }
    // -------------------------------------
}

JetNTuplizer::~JetNTuplizer() { }

// ------------ method called once each job just before starting event loop  ------------
void 
JetNTuplizer::beginJob()
{

}


// ------------ method called for each event  ------------
void
JetNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();

    edm::Handle<std::vector<reco::GenJet>> genjets;
    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    edm::Handle<std::vector<l1t::PFJet>> scjets;
    edm::Handle<std::vector<l1t::PFJet>> scngjets;
    edm::Handle<std::vector<l1t::PFJet>> scjetsCorr;
    edm::Handle<std::vector<l1t::PFTau>> nntaus;
    edm::Handle<std::vector<l1t::TkElectron>> electrons;
    edm::Handle<std::vector<l1t::SAMuon>> muons;

    edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetsFlavour;
    
    if (iEvent.isRealData()){
        isMC_ = false;
    }else{
        isMC_ = true;
    }
    if (isMC_){
        iEvent.getByToken(genjets_, genjets);
        iEvent.getByToken(genparticles_, genparticles);
        // fill genparticle categories once per event
        fill_genParticles(iEvent);

        iEvent.getByToken(genJetsFlavour_, genJetsFlavour);
        
    }
    iEvent.getByToken(scjets_, scjets);
    iEvent.getByToken(scngjets_, scngjets);
    iEvent.getByToken(scjetsCorr_, scjetsCorr);
    iEvent.getByToken(nntaus_, nntaus);
    iEvent.getByToken(electrons_, electrons);
    iEvent.getByToken(muons_, muons);

    float vz = 0.;
    double ptsum = 0;
    edm::Handle<std::vector<l1t::VertexWord>> vtxEmuHandle;
    iEvent.getByToken(fVtxEmu_, vtxEmuHandle);
    for (const auto& vtx : *vtxEmuHandle) {
        if (ptsum == 0 || vtx.pt() > ptsum) {
            ptsum = vtx.pt();
            vz = vtx.z0();
        }
    }

    edm::Handle<edm::ValueMap<float>> bjetIDhandle;
    iEvent.getByToken(bjetids_, bjetIDhandle);

    // gen jets
    std::vector<reco::GenJetRef> jetv_gen;  
    if(genjets.isValid()){
        for (auto jets_iter = genjets->begin(); jets_iter != genjets->end(); ++jets_iter) {                                                                                                   
        reco::GenJetRef jref (genjets, jets_iter - genjets->begin());                                                                                                                      
        jetv_gen.push_back(jref);                                                                                                                                                              
        }
        sort(jetv_gen.begin(), jetv_gen.end(), genJetRefSorter);
    }

    // reco jets
    std::vector<l1t::PFJetRef> jetv_l1;
    for (auto jets_iter = scjets->begin(); jets_iter != scjets->end(); ++jets_iter) {                                                                                                   
        l1t::PFJetRef jref(scjets, jets_iter - scjets->begin());                                                                                                                
        if (jref->pt() < jetPtMin_) continue;
        if (fabs(jref->eta()) > jetEtaMax_) continue;                 
        if (fabs(jref->eta()) < jetEtaMin_) continue;                 
        jetv_l1.push_back(jref);                                                                                                                                                              
    }
    sort(jetv_l1.begin(), jetv_l1.end(), jetRefSorter);

        // reco jets
    std::vector<l1t::PFJetRef> ngjetv_l1;
    for (auto ngjets_iter = scngjets->begin(); ngjets_iter != scngjets->end(); ++ngjets_iter) {                                                                                                   
        l1t::PFJetRef ngjref(scngjets, ngjets_iter - scngjets->begin());                                                                                                                
        if (ngjref->pt() < jetPtMin_) continue;
        if (fabs(ngjref->eta()) > jetEtaMax_) continue;                 
        if (fabs(ngjref->eta()) < jetEtaMin_) continue;                 
        ngjetv_l1.push_back(ngjref);                                                                                                                                                              
    }
    sort(ngjetv_l1.begin(), ngjetv_l1.end(), jetRefSorter);

    // reco jets, with jec
    std::vector<l1t::PFJetRef> jetv_l1_corr;
    for (auto jets_iter_corr = scjetsCorr->begin(); jets_iter_corr != scjetsCorr->end(); ++jets_iter_corr) {                                                                                                   
        l1t::PFJetRef jrefcorr(scjetsCorr, jets_iter_corr - scjetsCorr->begin());              
        jetv_l1_corr.push_back(jrefcorr);                                                                                                                                                              
    }
    sort(jetv_l1_corr.begin(), jetv_l1_corr.end(), jetRefSorter);

    // reco taus
    std::vector<l1t::PFTauRef> tauv_l1;
    for (auto tau_iter = nntaus->begin(); tau_iter != nntaus->end(); ++tau_iter) {                                                                                                   
        l1t::PFTauRef tref(nntaus, tau_iter - nntaus->begin());                                                                                                                             
        if (tref->pt() < 5.) continue;            
        tauv_l1.push_back(tref);                                                                                                                                                              
    }
    sort(tauv_l1.begin(), tauv_l1.end(), tauRefSorter);

    // reco electron
    std::vector<l1t::TkElectronRef> elev_l1;
    for (auto ele_iter = electrons->begin(); ele_iter != electrons->end(); ++ele_iter) {                                                                                                   
        l1t::TkElectronRef eref(electrons, ele_iter - electrons->begin());                                                                                                                             
        if (eref->pt() < 5.) continue;            
        elev_l1.push_back(eref);                                                                                                                                                              
    }
    sort(elev_l1.begin(), elev_l1.end(), eleRefSorter);

    // reco muons
    std::vector<l1t::SAMuonRef> muonv_l1;
    for (auto muon_iter = muons->begin(); muon_iter != muons->end(); ++muon_iter) {                                                                                                   
        l1t::SAMuonRef mref(muons, muon_iter - muons->begin());                                                                                                                             
        if (mref->pt() < 5.) continue;            
        muonv_l1.push_back(mref);                                                                                                                                                              
    }
    sort(muonv_l1.begin(), muonv_l1.end(), muonRefSorter);

    for (size_t i = 0; i < ngjetv_l1.size(); i++) {\
        l1ct::Jet ctngJet = l1ct::Jet::unpack(ngjetv_l1[i]->getHWJetCT());
        if (ngjetv_l1[i]->getTagScores().size() != 0) {
            jet_SC4NGJet_score_light_ = float(ctngJet.hwTagScores[2] );
            jet_SC4NGJet_score_b_ = float(ctngJet.hwTagScores[0] );
            jet_SC4NGJet_score_taup_ = float(ctngJet.hwTagScores[4] );
            jet_SC4NGJet_score_taum_ = float(ctngJet.hwTagScores[5] );
            jet_SC4NGJet_score_gluon_ = float(ctngJet.hwTagScores[3] );
            jet_SC4NGJet_score_charm_ = float(ctngJet.hwTagScores[1] );
            jet_SC4NGJet_score_muon_ = float(ctngJet.hwTagScores[6] );
            jet_SC4NGJet_score_electron_ = float(ctngJet.hwTagScores[7] );
            jet_SC4NGJet_score_regression_ = ngjetv_l1[i]->getPtCorrection();
        }
    }

    // loop over reco jets
    for (size_t i = 0; i < jetv_l1.size(); i++) {
        
        // get hardware value jet
        l1ct::Jet ctJet = l1ct::Jet::unpack(jetv_l1[i]->getHWJetCT());
        jet_pt_ = float(ctJet.hwPt);
        jet_pt_phys_ = jetv_l1[i]->pt();
        jet_eta_ = float(ctJet.hwEta);
        jet_eta_phys_ = jetv_l1[i]->eta();
        jet_phi_ = float(ctJet.hwPhi);
        jet_phi_phys_ = jetv_l1[i]->phi();

        // jet_pt_ = jetv_l1[i]->pt();
        // jet_eta_ = jetv_l1[i]->eta();
        // jet_phi_ = jetv_l1[i]->phi();
        jet_mass_ = jetv_l1[i]->mass();
        jet_pt_raw_ = jetv_l1[i]->rawPt();
        jet_energy_ = jetv_l1[i]->energy();
        // jet_px_ = jetv_l1[i]->px();
        // jet_py_ = jetv_l1[i]->py();
        // jet_pz_ = jetv_l1[i]->pz();

        jet_bjetscore_ = (*bjetIDhandle)[jetv_l1[i]];
        
        // match to GEN
        int   pos_matched = -1;
        float minDR = dRJetGenMatch_;
        for(size_t igen = 0; igen < jetv_gen.size(); igen++){
            if(jetv_l1[i]->pt() <= 0.1 * jetv_gen[igen]->pt()) continue;
            if(reco::deltaR(jetv_gen[igen]->p4(),jetv_l1[i]->p4()) < minDR){
                pos_matched = igen;
                minDR = reco::deltaR(jetv_gen[igen]->p4(),jetv_l1[i]->p4());
            }
        }
        
        if(pos_matched >= 0){
            jet_genmatch_pt_ = jetv_gen[pos_matched]->pt();
            jet_genmatch_eta_ = jetv_gen[pos_matched]->eta();
            jet_genmatch_phi_ = jetv_gen[pos_matched]->phi();
            jet_genmatch_mass_ = jetv_gen[pos_matched]->mass();
            jet_genmatch_dR_ = minDR;
            jet_genmatch_hflav_ = (*genJetsFlavour)[edm::RefToBase<reco::Jet>(jetv_gen[pos_matched])].getHadronFlavour();
            jet_genmatch_pflav_ = (*genJetsFlavour)[edm::RefToBase<reco::Jet>(jetv_gen[pos_matched])].getPartonFlavour();      
        }
        else{
            jet_genmatch_pt_ = 0;
            jet_genmatch_eta_ = 0;
            jet_genmatch_phi_ = 0;
            jet_genmatch_mass_ = 0;
            jet_genmatch_dR_ = 0;
            jet_genmatch_hflav_ = 0;
            jet_genmatch_pflav_ = 0;
        }


        // matching with gen-leptons (muons/electrons/hadronic taus)
        minDR = 1000;
        int nlep_in_cone  = 0;
        int pos_matched_genmu = -1;
        int pos_matched_genele = -1;
        int pos_matched_tauh = -1;
        int gentau_decaymode = -1; 	  
        TLorentzVector genLepton4V;
        TLorentzVector genLeptonVis4V;
        TLorentzVector jet4V;
        jet4V.SetPtEtaPhiM(jetv_l1[i]->pt(), jetv_l1[i]->eta(), jetv_l1[i]->phi(), jetv_l1[i]->mass());
        for(size_t igen = 0; igen < genMuonsFromResonance4V_.size(); igen++){
            float dR = jet4V.DeltaR(genMuonsFromResonance4V_.at(igen));	      
            if(dR < dRCone) nlep_in_cone++;
            if(dR < dRCone and dR < minDR){
            pos_matched_genmu = igen;
            minDR = dR;
            genLepton4V = genMuonsFromResonance4V_.at(igen);
            genLeptonVis4V = genMuonsFromResonance4V_.at(igen);
            }
        }

        for(size_t igen = 0; igen < genElectronsFromResonance4V_.size(); igen++){
            float dR = jet4V.DeltaR(genElectronsFromResonance4V_.at(igen));	      
            if(dR < dRCone) nlep_in_cone++;
            if(dR < dRCone and dR < minDR){
            pos_matched_genmu  = -1;
            pos_matched_genele = igen;
            minDR = dR;
            genLepton4V = genElectronsFromResonance4V_.at(igen);
            genLeptonVis4V = genElectronsFromResonance4V_.at(igen);
            }
        }

        for(size_t itau = 0; itau < tau_gen_visible_.size(); itau++){
            float dR = tau_gen_visible_.at(itau).DeltaR(jet4V); 
            if(dR < dRCone) nlep_in_cone++;
            if(dR < dRCone and dR < minDR){
            pos_matched_genmu  = -1;
            pos_matched_genele = -1;
            pos_matched_tauh = itau;
            minDR = dR;
            gentau_decaymode = 5*(tau_gen_nch_.at(itau)-1)+tau_gen_np0_.at(itau);
            genLepton4V = tau_gen_.at(itau);
            genLeptonVis4V = tau_gen_visible_.at(itau);
            }
        }

        // some basic selection crityeria
        jet_reject_  = false;
        // exclude, when a jet is matched with a lepton, those for which the matched lepton is below the chosen pt threshold
        // Jet id applied only to jets not overlapping with gen-leptons
        if(applyJetIDFlag && pos_matched_genmu  == -1 && pos_matched_genele == -1 && pos_matched_tauh == -1){
            jet_reject_ = true;
        }
        if(pos_matched_genmu != -1 and genLeptonVis4V.Pt() < ptGenLeptonMin){
            jet_reject_ = true;
        }
        if(pos_matched_genele != -1 and genLeptonVis4V.Pt() < ptGenLeptonMin){
            jet_reject_ = true;
        }
        if(pos_matched_tauh != -1 and genLeptonVis4V.Pt() < ptGenTauVisibleMin){
            jet_reject_ = true;
        } 

        // flavor info
        if(pos_matched_genmu >= 0){
            jet_muflav_  = 1;
        }
        else{
            jet_muflav_  = 0;
        }

        if(pos_matched_genele >= 0){
            jet_elflav_  = 1;
        }
        else{
            jet_elflav_  = 0;
        }

        if(pos_matched_tauh >= 0){
            jet_tauflav_ = 1;
            jet_taudecaymode_ = gentau_decaymode;
            jet_taucharge_ = tau_gen_charge_.at(pos_matched_tauh);
        }
        else{
            jet_tauflav_ = 0;
            jet_taudecaymode_ = -1;
            jet_taucharge_ = 0;
        }

        jet_lepflav_ = nlep_in_cone;
        jet_genmatch_lep_pt_ = genLepton4V.Pt();
        jet_genmatch_lep_vis_pt_ = genLeptonVis4V.Pt();
        jet_genmatch_lep_dR_ = minDR;

        // match to TAUS
        int   pos_matched_tau = -1;
        float minDR_tau = 0.4;
        for(size_t itau = 0; itau < tauv_l1.size(); itau++){
            if(reco::deltaR(tauv_l1[itau]->p4(),jetv_l1[i]->p4()) < minDR_tau){
                pos_matched_tau = itau;
                minDR_tau = reco::deltaR(tauv_l1[itau]->p4(),jetv_l1[i]->p4());
            }
        }
        if (pos_matched_tau > -1){
            jet_taumatch_dR_ = minDR_tau;
            jet_tauscore_ = tauv_l1[pos_matched_tau]->chargedIso();
            jet_taupt_ = tauv_l1[pos_matched_tau]->pt();
            // l1gt::Tau hwtau = tauv_l1[pos_matched_tau]->getHWTauGT();
            // jet_tauscoreHW_ = hwtau->hwRawId;
        }else{
            jet_taumatch_dR_ = 999.;
            jet_tauscore_ = -1.;
            jet_taupt_ = -1.;
            // jet_tauscoreHW_ = -1.;
        }

        // match to ELECTRONS
        int   pos_matched_ele = -1;
        float minDR_ele = 0.4;
        for(size_t iele = 0; iele < elev_l1.size(); iele++){
            if(reco::deltaR(elev_l1[iele]->p4(),jetv_l1[i]->p4()) < minDR_ele){
                pos_matched_ele = iele;
                minDR_ele = reco::deltaR(elev_l1[iele]->p4(), jetv_l1[i]->p4());
            }
        }
        if (pos_matched_ele > -1){
            jet_elematch_dR_ = minDR_ele;
            jet_eletkiso_ = elev_l1[pos_matched_ele]->trkIsol();
            jet_elepfiso_ = elev_l1[pos_matched_ele]->pfIsol();
            jet_elepuppiiso_ = elev_l1[pos_matched_ele]->puppiIsol();
        }else{
            jet_elematch_dR_ = 999.;
            jet_eletkiso_ = -1.;
            jet_elepfiso_ = -1.;
            jet_elepuppiiso_ = -1.;
        }

        // match to MUONS
        int   pos_matched_muon = -1;
        float minDR_muon = 0.4;
        for(size_t imuon = 0; imuon < muonv_l1.size(); imuon++){
            if(reco::deltaR(muonv_l1[imuon]->p4(),jetv_l1[i]->p4()) < minDR_muon){
                pos_matched_muon = imuon;
                minDR_muon = reco::deltaR(muonv_l1[imuon]->p4(), jetv_l1[i]->p4());
            }
        }
        if (pos_matched_muon > -1){
            jet_muonmatch_dR_ = minDR_muon;
            jet_muoniso_ = muonv_l1[pos_matched_muon]->hwIso();
        }else{
            jet_muonmatch_dR_ = 999.;
            jet_muoniso_ = -1.;
        }

        // match to jec corrected jets
        int   pos_matched_jec = -1;
        float minDR_jec = 0.4;
        for(size_t ijec = 0; ijec < jetv_l1_corr.size(); ijec++){
            if(reco::deltaR(jetv_l1_corr[ijec]->p4(),jetv_l1[i]->p4()) < minDR_jec){
                pos_matched_jec = ijec;
                minDR_jec = reco::deltaR(jetv_l1_corr[ijec]->p4(),jetv_l1[i]->p4());
            }
        }
        if (pos_matched_jec > -1){
            jet_jecmatch_dR_ = minDR_jec;
            jet_pt_corr_ = jetv_l1_corr[pos_matched_jec]->pt();
        }else{
            jet_jecmatch_dR_ = 999.;
            jet_pt_corr_ = -1.;
        }
        
        // pf candidates
        std::vector<l1t::PFCandidate> vectorOfConstituents;
        for(unsigned ipart = 0; ipart < jetv_l1[i]->numberOfDaughters(); ipart++){
            const l1t::PFCandidate* pfPart = dynamic_cast<const l1t::PFCandidate*> (jetv_l1[i]->daughter(ipart));      
            vectorOfConstituents.push_back(*pfPart);
        }

        TVector3 jet_direction (jetv_l1[i]->momentum().Unit().x(),jetv_l1[i]->momentum().Unit().y(),jetv_l1[i]->momentum().Unit().z());
        GlobalVector jet_global_vec (jetv_l1[i]->px(),jetv_l1[i]->py(),jetv_l1[i]->pz());
        math::XYZVector jetDir = jetv_l1[i]->momentum().Unit();
        GlobalVector jetRefTrackDir(jetv_l1[i]->px(),jetv_l1[i]->py(),jetv_l1[i]->pz());


        // from PNET ntupler
        std::sort(vectorOfConstituents.begin(),vectorOfConstituents.end(),l1PFCandidateSorter);
    

        jet_pfcand_pt.clear();
        jet_pfcand_pt_phys.clear();
        jet_pfcand_pt_rel.clear();
        jet_pfcand_pt_rel_phys.clear();
        jet_pfcand_isfilled.clear();
        // jet_pfcand_px.clear();
        // jet_pfcand_py.clear();
        // jet_pfcand_pz.clear();
        jet_pfcand_pt_log.clear();
        jet_pfcand_pt_rel_log.clear();
        jet_pfcand_eta.clear();
        jet_pfcand_eta_phys.clear();
        jet_pfcand_phi.clear();
        jet_pfcand_phi_phys.clear();
        jet_pfcand_mass.clear();
        jet_pfcand_energy.clear();
        jet_pfcand_energy_log.clear();
        jet_pfcand_puppiweight.clear();
        jet_pfcand_emid.clear();
        jet_pfcand_quality.clear();
        jet_pfcand_tkquality.clear();
        jet_pfcand_z0.clear();
        jet_pfcand_dxy.clear();
        jet_pfcand_dxy_custom.clear();
        jet_pfcand_id.clear();
        jet_pfcand_isPhoton.clear();
        jet_pfcand_isElectron.clear();
        jet_pfcand_isElectronPlus.clear();
        jet_pfcand_isElectronMinus.clear();
        jet_pfcand_isMuon.clear();
        jet_pfcand_isMuonPlus.clear();
        jet_pfcand_isMuonMinus.clear();
        jet_pfcand_isNeutralHadron.clear();
        jet_pfcand_isChargedHadron.clear();
        jet_pfcand_isChargedHadronPlus.clear();
        jet_pfcand_isChargedHadronMinus.clear();

        jet_pfcand_charge.clear();
        jet_pfcand_pperp_ratio.clear();
        jet_pfcand_ppara_ratio.clear();
        jet_pfcand_deta.clear();
        jet_pfcand_deta_phys.clear();
        jet_pfcand_dphi.clear();
        jet_pfcand_dphi_phys.clear();
        jet_pfcand_etarel.clear();

        jet_pfcand_track_valid.clear();
        jet_pfcand_track_isfilled.clear();
        jet_pfcand_track_rinv.clear();
        jet_pfcand_track_phizero.clear();
        jet_pfcand_track_tanl.clear();
        jet_pfcand_track_z0.clear();
        jet_pfcand_track_d0.clear();
        jet_pfcand_track_chi2rphi.clear();
        jet_pfcand_track_chi2rz.clear();
        jet_pfcand_track_bendchi2.clear();
        jet_pfcand_track_hitpattern.clear();
        jet_pfcand_track_nstubs.clear();
        jet_pfcand_track_mvaquality.clear();
        jet_pfcand_track_mvaother.clear();

        jet_pfcand_track_chi2.clear();
        jet_pfcand_track_chi2norm.clear();
        jet_pfcand_track_qual.clear();
        jet_pfcand_track_npar.clear();
        // jet_pfcand_track_nstubs.clear();
        jet_pfcand_track_vx.clear();
        jet_pfcand_track_vy.clear();
        jet_pfcand_track_vz.clear();
        jet_pfcand_track_pterror.clear();

        jet_pfcand_cluster_isfilled.clear();
        jet_pfcand_cluster_hovere.clear();
        jet_pfcand_cluster_sigmarr.clear();
        jet_pfcand_cluster_abszbarycenter.clear();
        jet_pfcand_cluster_emet.clear();
        jet_pfcand_cluster_egvspion.clear();
        jet_pfcand_cluster_egvspu.clear();

        njet_pfcand_ = 0;
        // for(auto const & pfcand : vectorOfConstituents){
        // for(l1t::PFCandidate const & pfcand : vectorOfConstituents){
        for(size_t ipfcand = 0; ipfcand < vectorOfConstituents.size(); ipfcand++){
        
            l1t::PFCandidate pfcand = vectorOfConstituents.at(ipfcand);
            if(pfcand.pt() < jetPFCandidatePtMin_) continue;
            njet_pfcand_++;

            // jet_pfcand_pt.push_back(pfcand.pt());
            jet_pfcand_pt.push_back(float(pfcand.hwPt()));
            jet_pfcand_pt_phys.push_back(pfcand.pt());
            jet_pfcand_pt_rel_phys.push_back(pfcand.pt()/jetv_l1[i]->pt());
            jet_pfcand_pt_rel.push_back(float(pfcand.hwPt())/jet_pt_);
            jet_pfcand_isfilled.push_back(1);
            // jet_pfcand_px.push_back(pfcand.px());
            // jet_pfcand_py.push_back(pfcand.py());
            // jet_pfcand_pz.push_back(pfcand.py());
            // jet_pfcand_pt_log.push_back(std::log(pfcand.pt()));
            jet_pfcand_pt_log.push_back(std::log(float(pfcand.hwPt())));
            // jet_pfcand_pt_rel_log.push_back(std::log(pfcand.pt()/jetv_l1[i]->pt()));
            jet_pfcand_pt_rel_log.push_back(std::log(float(pfcand.hwPt())/jet_pt_));
            // jet_pfcand_eta.push_back(pfcand.eta());
            jet_pfcand_eta.push_back(float(pfcand.hwEta()));
            jet_pfcand_eta_phys.push_back(pfcand.eta());
            // jet_pfcand_phi.push_back(pfcand.phi());
            jet_pfcand_phi.push_back(float(pfcand.hwPhi()));
            jet_pfcand_phi_phys.push_back(pfcand.phi());

            float massCand = 0.13f;
            if (abs(pfcand.charge())) {
                if ((pfcand.id() == l1t::PFCandidate::Muon)) {
                    // type = l1t::PFCandidate::Muon;
                    massCand = 0.105;
                } else if ((pfcand.id() == l1t::PFCandidate::Electron)) {
                    // type = l1t::PFCandidate::Electron;
                    massCand = 0.005;
                } //else
                    // type = l1t::PFCandidate::ChargedHadron;
            } else {
            // type = pfcand.hwId.isPhoton() ? l1t::PFCandidate::Photon : l1t::PFCandidate::NeutralHadron;
                massCand = pfcand.id() == l1t::PFCandidate::Photon ? 0.0 : 0.5;
            }

            // jet_pfcand_mass.push_back(pfcand.mass());
            jet_pfcand_mass.push_back(massCand);
            jet_pfcand_energy.push_back(pfcand.energy());
            jet_pfcand_energy_log.push_back(std::log(pfcand.energy()));
            // jet_pfcand_puppiweight.push_back(pfcand.puppiWeight());
            jet_pfcand_puppiweight.push_back(float(pfcand.hwPuppiWeight()));

            jet_pfcand_emid.push_back(float(pfcand.hwEmID()));
            jet_pfcand_quality.push_back(float(pfcand.hwTkQuality()));
            
            jet_pfcand_id.push_back(abs(pfcand.pdgId()));
            // std::cout<<pfcand.hwQual()<<" "<<pfcand.id()<<" "<<pfcand.pdgId()<<" "<<pfcand.charge()<<std::endl;
            jet_pfcand_isPhoton.push_back(pfcand.id() == l1t::PFCandidate::Photon);
            jet_pfcand_isElectron.push_back(pfcand.id() == l1t::PFCandidate::Electron);
            jet_pfcand_isElectronPlus.push_back((pfcand.id() == l1t::PFCandidate::Electron) && (pfcand.charge() > 0.));
            jet_pfcand_isElectronMinus.push_back((pfcand.id() == l1t::PFCandidate::Electron) && (pfcand.charge() < 0.));
            jet_pfcand_isMuon.push_back(pfcand.id() == l1t::PFCandidate::Muon);
            jet_pfcand_isMuonPlus.push_back((pfcand.id() == l1t::PFCandidate::Muon) && (pfcand.charge() > 0.));
            jet_pfcand_isMuonMinus.push_back((pfcand.id() == l1t::PFCandidate::Muon) && (pfcand.charge() < 0.));
            jet_pfcand_isNeutralHadron.push_back(pfcand.id() == l1t::PFCandidate::NeutralHadron);
            jet_pfcand_isChargedHadron.push_back(pfcand.id() == l1t::PFCandidate::ChargedHadron);
            jet_pfcand_isChargedHadronPlus.push_back((pfcand.id() == l1t::PFCandidate::ChargedHadron) && (pfcand.charge() > 0.));
            jet_pfcand_isChargedHadronMinus.push_back((pfcand.id() == l1t::PFCandidate::ChargedHadron) && (pfcand.charge() < 0.));
            // int sumID = jet_pfcand_isPhoton.at(ipfcand)+jet_pfcand_isElectronPlus.at(ipfcand)+jet_pfcand_isElectronMinus.at(ipfcand)+jet_pfcand_isMuonMinus.at(ipfcand)+jet_pfcand_isMuonPlus.at(ipfcand)+jet_pfcand_isNeutralHadron.at(ipfcand)+jet_pfcand_isChargedHadronMinus.at(ipfcand)+jet_pfcand_isChargedHadronPlus.at(ipfcand);
            // std::cout<<sumID<<std::endl;

            jet_pfcand_charge.push_back(pfcand.charge());
            
            // jet_pfcand_z0.push_back(pfcand.z0());
            jet_pfcand_z0.push_back(float(pfcand.hwZ0()));
            // jet_pfcand_dxy.push_back(pfcand.dxy());
            jet_pfcand_dxy.push_back(float(pfcand.hwDxy()));

            jet_pfcand_tkquality.push_back(float(pfcand.hwTkQuality()));
            
            // Track related
            
            // TVector3 pfcand_momentum (pfcand.momentum().x(),pfcand.momentum().y(),pfcand.momentum().z());
            TVector3 pfcand_momentum (pfcand.momentum().x(),pfcand.momentum().y(),pfcand.momentum().z());
            jet_pfcand_pperp_ratio.push_back(jet_direction.Perp(pfcand_momentum)/pfcand_momentum.Mag());
            jet_pfcand_ppara_ratio.push_back(jet_direction.Dot(pfcand_momentum)/pfcand_momentum.Mag());

            jet_pfcand_dphi_phys.push_back(jet_direction.DeltaPhi(pfcand_momentum));
            // jet_pfcand_dphi.push_back(deltaPhi(float(pfcand.hwPhi()), jet_phi_));
            // https://github.com/cms-sw/cmssw/blob/master/L1Trigger/Phase2L1ParticleFlow/src/jetmet/L1SeedConePFJetEmulator.cc#L11

            L1SCJetEmu::detaphi_t dphi(pfcand.hwPhi() - jet_phi_);
            // phi wrap
            L1SCJetEmu::detaphi_t dphi0 = dphi > L1SCJetEmu::detaphi_t(l1ct::Scales::INTPHI_PI) ? L1SCJetEmu::detaphi_t(l1ct::Scales::INTPHI_TWOPI - dphi) : L1SCJetEmu::detaphi_t(dphi);
            L1SCJetEmu::detaphi_t dphi1 = dphi < L1SCJetEmu::detaphi_t(-l1ct::Scales::INTPHI_PI) ? L1SCJetEmu::detaphi_t(l1ct::Scales::INTPHI_TWOPI + dphi) : L1SCJetEmu::detaphi_t(dphi);
            L1SCJetEmu::detaphi_t dphiw = dphi > L1SCJetEmu::detaphi_t(0) ? dphi0 : dphi1;
            jet_pfcand_dphi.push_back(float(dphiw));

            jet_pfcand_deta_phys.push_back(jet_direction.Eta()-pfcand_momentum.Eta());
            jet_pfcand_deta.push_back(jet_eta_-float(pfcand.hwEta()));
            // jet_pfcand_etarel.push_back(etaRel(jetDir,pfcand.momentum()));
            jet_pfcand_etarel.push_back(etaRel(jetDir,pfcand.momentum()));
            
            const l1t::PFTrackRef track = pfcand.pfTrack();      
        
            if(track.isNonnull()){ // need valid track object
            
                const TTTrack_TrackWord trackWord = track->trackWord();    

                jet_pfcand_track_isfilled.push_back(1);
                jet_pfcand_track_valid.push_back(trackWord.getValid());
                jet_pfcand_track_rinv.push_back(trackWord.getRinv());
                jet_pfcand_track_phizero.push_back(trackWord.getPhi());
                jet_pfcand_track_tanl.push_back(trackWord.getTanl());
                jet_pfcand_track_z0.push_back(trackWord.getZ0());
                jet_pfcand_track_d0.push_back(trackWord.getD0());
                jet_pfcand_track_chi2rphi.push_back(trackWord.getChi2RPhi());
                jet_pfcand_track_chi2rz.push_back(trackWord.getChi2RZ());
                jet_pfcand_track_bendchi2.push_back(trackWord.getBendChi2());
                jet_pfcand_track_hitpattern.push_back(trackWord.getHitPattern());
                jet_pfcand_track_nstubs.push_back(trackWord.getNStubs());
                jet_pfcand_track_mvaquality.push_back(trackWord.getMVAQuality());
                jet_pfcand_track_mvaother.push_back(trackWord.getMVAOther());

                jet_pfcand_track_chi2.push_back(track->chi2());
                jet_pfcand_track_chi2norm.push_back(track->normalizedChi2());
                jet_pfcand_track_qual.push_back(track->quality());
                jet_pfcand_track_npar.push_back(track->nPar());
                // jet_pfcand_track_nstubs.push_back(track->nStubs());
                jet_pfcand_track_vx.push_back(track->vx());
                jet_pfcand_track_vy.push_back(track->vy());
                jet_pfcand_track_vz.push_back(track->vz()-vz);
                jet_pfcand_track_pterror.push_back(track->trkPtError());
                jet_pfcand_dxy_custom.push_back(-track->vx() * sin(track->phi()) + track->vy() * cos(track->phi()));  
            }else{

                jet_pfcand_track_isfilled.push_back(0);
                jet_pfcand_track_valid.push_back(0);
                jet_pfcand_track_rinv.push_back(0);
                jet_pfcand_track_phizero.push_back(0);
                jet_pfcand_track_tanl.push_back(0);
                jet_pfcand_track_z0.push_back(0);
                jet_pfcand_track_d0.push_back(0);
                jet_pfcand_track_chi2rphi.push_back(0);
                jet_pfcand_track_chi2rz.push_back(0);
                jet_pfcand_track_bendchi2.push_back(0);
                jet_pfcand_track_hitpattern.push_back(0);
                jet_pfcand_track_nstubs.push_back(0);
                jet_pfcand_track_mvaquality.push_back(0);
                jet_pfcand_track_mvaother.push_back(0);

                jet_pfcand_track_chi2.push_back(0);
                jet_pfcand_track_chi2norm.push_back(0);
                jet_pfcand_track_qual.push_back(0);
                jet_pfcand_track_npar.push_back(0);
                // jet_pfcand_track_nstubs.push_back(0);
                jet_pfcand_track_vx.push_back(0);
                jet_pfcand_track_vy.push_back(0);
                jet_pfcand_track_vz.push_back(0);
                jet_pfcand_track_pterror.push_back(0);
                jet_pfcand_dxy_custom.push_back(0);
            }	  

            const l1t::PFClusterRef cluster = pfcand.pfCluster();    

            if(cluster.isNonnull()){ // need valid cluster object

                jet_pfcand_cluster_isfilled.push_back(1);
                jet_pfcand_cluster_hovere.push_back(cluster->hOverE());
                jet_pfcand_cluster_sigmarr.push_back(cluster->sigmaRR());
                jet_pfcand_cluster_abszbarycenter.push_back(cluster->absZBarycenter());
                jet_pfcand_cluster_emet.push_back(cluster->emEt());
                jet_pfcand_cluster_egvspion.push_back(cluster->egVsPionMVAOut());
                jet_pfcand_cluster_egvspu.push_back(cluster->egVsPUMVAOut());

            }else{

                jet_pfcand_cluster_isfilled.push_back(0);
                jet_pfcand_cluster_hovere.push_back(0);
                jet_pfcand_cluster_sigmarr.push_back(0);
                jet_pfcand_cluster_abszbarycenter.push_back(0);
                jet_pfcand_cluster_emet.push_back(0);
                jet_pfcand_cluster_egvspion.push_back(0);
                jet_pfcand_cluster_egvspu.push_back(0);

            }
        }

        tree_->Fill();
    }


}





void JetNTuplizer::fill_genParticles(const edm::Event& iEvent)
{
    gToBB_.clear();
    gToCC_.clear();
    neutrinosLepB_.clear();
    neutrinosLepB_C_.clear();
    alltaus_.clear();
    Bhadron_.clear();
    Bhadron_daughter_.clear();

    // Generator-level information (GEN particles)
    gen_particle_pt.clear();
    gen_particle_eta.clear();
    gen_particle_phi.clear();
    gen_particle_mass.clear();
    gen_particle_id.clear();
    gen_particle_status.clear();
    gen_particle_daughters_id.clear();
    gen_particle_daughters_igen.clear();
    gen_particle_daughters_status.clear();
    gen_particle_daughters_pt.clear();
    gen_particle_daughters_eta.clear();
    gen_particle_daughters_phi.clear();
    gen_particle_daughters_mass.clear();
    gen_particle_daughters_charge.clear();

    genLepFromResonance4V_.clear();
    genMuonsFromResonance4V_.clear();
    genElectronsFromResonance4V_.clear();
    tau_gen_visible_.clear();
    tau_gen_.clear();
    tau_gen_charge_.clear();
    tau_gen_nch_.clear();
    tau_gen_np0_.clear();
    tau_gen_nnh_.clear();


    if(!iEvent.isRealData())
    {
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByToken(genparticles_, genParticles);

        for (const reco::Candidate &genC : *genParticles){
            const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
            if((abs(gen.pdgId())>500&&abs(gen.pdgId())<600)||(abs(gen.pdgId())>5000&&abs(gen.pdgId())<6000)){
                Bhadron_.push_back(gen);
                if(gen.numberOfDaughters()>0){
                    if( (abs(gen.daughter(0)->pdgId())>500&&abs(gen.daughter(0)->pdgId())<600)||(abs(gen.daughter(0)->pdgId())>5000&&abs(gen.daughter(0)->pdgId())<6000)){
                        if(gen.daughter(0)->numberOfDaughters()>0){
                            const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*(gen.daughter(0)->daughter(0)));
                            if(daughter_.vx()!=gen.vx()){
                                Bhadron_daughter_.push_back(daughter_);
                            }else Bhadron_daughter_.push_back(gen);
                        }else  Bhadron_daughter_.push_back(gen);
                    }else{
                        const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*gen.daughter(0));
                        Bhadron_daughter_.push_back(daughter_);
                    }
                }else {
                    Bhadron_daughter_.push_back(gen);
                }
            }
        }

        for (const reco::Candidate &genC : *genParticles){
            const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
            if(abs(gen.pdgId())==12||abs(gen.pdgId())==14||abs(gen.pdgId())==16){
                const reco::GenParticle* mother =  static_cast< const reco::GenParticle*> (gen.mother());
                if(mother!=NULL){
                    if((abs(mother->pdgId())>500&&abs(mother->pdgId())<600)||(abs(mother->pdgId())>5000&&abs(mother->pdgId())<6000)){
                        neutrinosLepB_.emplace_back(gen);
                    }
                    if((abs(mother->pdgId())>400&&abs(mother->pdgId())<500)||(abs(mother->pdgId())>4000&&abs(mother->pdgId())<5000)){
                        neutrinosLepB_C_.emplace_back(gen);
                    }
                }else {
                    std::cout << "No mother" << std::endl;
                }
            }
            int id(std::abs(gen.pdgId()));
            int status(gen.status());
            if (id == 21 && status >= 21 && status <= 59){ //// Pythia8 hard scatter, ISR, or FSR
                if ( gen.numberOfDaughters() == 2 ){
                    const reco::Candidate* d0 = gen.daughter(0);
                    const reco::Candidate* d1 = gen.daughter(1);
                    if ( std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5 && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4){
                        gToBB_.push_back(gen);
                    }
                    if ( std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4 && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4){
                        gToCC_.push_back(gen);
                    }
                }
            }
            if(id == 15 && false){
                alltaus_.push_back(gen);
            }
        }

        // ---------------------------------------
        // from PNET

  // GEN particle informatio
    if(genParticles.isValid()){
        unsigned int igen = 0;
        for (auto gens_iter = genParticles->begin(); gens_iter != genParticles->end(); ++gens_iter) {      

            if((abs(gens_iter->pdgId()) == 25 or abs(gens_iter->pdgId()) == 24 or abs(gens_iter->pdgId()) == 23) and
            gens_iter->isLastCopy() and 
            gens_iter->statusFlags().fromHardProcess()){ 

                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());

                for(size_t idau = 0; idau < gens_iter->numberOfDaughters(); idau++){
                    gen_particle_daughters_id.push_back(gens_iter->daughter(idau)->pdgId());
                    gen_particle_daughters_igen.push_back(igen);
                    gen_particle_daughters_pt.push_back(gens_iter->daughter(idau)->pt());
                    gen_particle_daughters_eta.push_back(gens_iter->daughter(idau)->eta());
                    gen_particle_daughters_phi.push_back(gens_iter->daughter(idau)->phi());
                    gen_particle_daughters_mass.push_back(gens_iter->daughter(idau)->mass());
                    gen_particle_daughters_status.push_back(gens_iter->daughter(idau)->status());
                    gen_particle_daughters_charge.push_back(gens_iter->daughter(idau)->charge());
                }
                igen++;
            }

            // Final states Leptons (e,mu) and Neutrinos --> exclude taus. They need to be prompt or from Tau decay      
            if (abs(gens_iter->pdgId()) > 10 and abs(gens_iter->pdgId()) < 17 and abs(gens_iter->pdgId()) != 15  and 
            (gens_iter->isPromptFinalState() or 
            gens_iter->isDirectPromptTauDecayProductFinalState())) { 

                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());

                // No need to save daughters here
                igen++;
            }

            // Final state quarks or gluons from the hard process before the shower --> partons in which H/Z/W/top decay into
            if (((abs(gens_iter->pdgId()) >= 1 and abs(gens_iter->pdgId()) <= 5) or abs(gens_iter->pdgId()) == 21) and 
            gens_iter->statusFlags().fromHardProcess() and 
            gens_iter->statusFlags().isFirstCopy()){
                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());
                igen++;
                // no need to save daughters
            }

            // Special case of taus: last-copy, from hard process and, prompt and decayed
            if(abs(gens_iter->pdgId()) == 15 and 
            gens_iter->isLastCopy() and
            gens_iter->statusFlags().fromHardProcess() and
            gens_iter->isPromptDecayed()){ // hadronic taus

                gen_particle_pt.push_back(gens_iter->pt());
                gen_particle_eta.push_back(gens_iter->eta());
                gen_particle_phi.push_back(gens_iter->phi());
                gen_particle_mass.push_back(gens_iter->mass());
                gen_particle_id.push_back(gens_iter->pdgId());
                gen_particle_status.push_back(gens_iter->status());

                // only store the final decay particles
                for(size_t idau = 0; idau < gens_iter->numberOfDaughters(); idau++){
                    if(not dynamic_cast<const reco::GenParticle*>(gens_iter->daughter(idau))->statusFlags().isPromptTauDecayProduct()) continue;
                    gen_particle_daughters_id.push_back(gens_iter->daughter(idau)->pdgId());
                    gen_particle_daughters_igen.push_back(igen);
                    gen_particle_daughters_pt.push_back(gens_iter->daughter(idau)->pt());
                    gen_particle_daughters_eta.push_back(gens_iter->daughter(idau)->eta());
                    gen_particle_daughters_phi.push_back(gens_iter->daughter(idau)->phi());
                    gen_particle_daughters_mass.push_back(gens_iter->daughter(idau)->mass());    
                    gen_particle_daughters_status.push_back(gens_iter->daughter(idau)->status());
                    gen_particle_daughters_charge.push_back(gens_iter->daughter(idau)->charge());
                }
                igen++;
            }  
        }
    }



        for(size_t igen = 0; igen < gen_particle_pt.size(); igen++){
            // select resonances like Higgs, W, Z, taus
            if(abs(gen_particle_id.at(igen)) == 25 or
            abs(gen_particle_id.at(igen)) == 23 or
            abs(gen_particle_id.at(igen)) == 24 or
            abs(gen_particle_id.at(igen)) == 15){	
            for(size_t idau = 0; idau < gen_particle_daughters_id.size(); idau++){
                // select electrons or muons from the resonance / tau decay
                if(gen_particle_daughters_igen.at(idau) == igen and
                (abs(gen_particle_daughters_id.at(idau)) == 11 or
                abs(gen_particle_daughters_id.at(idau)) == 13)){
                    TLorentzVector gen4V;
                    gen4V.SetPtEtaPhiM(gen_particle_daughters_pt.at(idau),gen_particle_daughters_eta.at(idau),gen_particle_daughters_phi.at(idau),gen_particle_daughters_mass.at(idau));
                    if(std::find(genLepFromResonance4V_.begin(),genLepFromResonance4V_.end(),gen4V) == genLepFromResonance4V_.end())
                    genLepFromResonance4V_.push_back(gen4V);
                    if(abs(gen_particle_daughters_id.at(idau)) == 13 and 
                    std::find(genMuonsFromResonance4V_.begin(),genMuonsFromResonance4V_.end(),gen4V) == genMuonsFromResonance4V_.end()){
                    genMuonsFromResonance4V_.push_back(gen4V);
                    }
                    if(abs(gen_particle_daughters_id.at(idau)) == 11 and 
                    std::find(genElectronsFromResonance4V_.begin(),genElectronsFromResonance4V_.end(),gen4V) == genElectronsFromResonance4V_.end()){
                        genElectronsFromResonance4V_.push_back(gen4V);		
                    }
                }
            }
            }
        }   

        
        // Gen hadronic taus	
        for(size_t igen = 0; igen < gen_particle_pt.size(); igen++){
            if(abs(gen_particle_id.at(igen)) == 15){ // hadronic or leptonic tau
                TLorentzVector tau_gen_tmp;
                unsigned int tau_gen_nch_tmp(0);
                unsigned int tau_gen_np0_tmp(0);
                unsigned int tau_gen_nnh_tmp(0);
                for(size_t idau = 0; idau < gen_particle_daughters_pt.size(); idau++){
                    if(gen_particle_daughters_igen.at(idau) == igen and
                    abs(gen_particle_daughters_id.at(idau)) != 11 and // no mu
                    abs(gen_particle_daughters_id.at(idau)) != 13 and // no el
                    abs(gen_particle_daughters_id.at(idau)) != 12 and // no neutrinos
                    abs(gen_particle_daughters_id.at(idau)) != 14 and
                    abs(gen_particle_daughters_id.at(idau)) != 16){
                    TLorentzVector tmp4V; 
                    tmp4V.SetPtEtaPhiM(gen_particle_daughters_pt.at(idau),gen_particle_daughters_eta.at(idau),gen_particle_daughters_phi.at(idau),gen_particle_daughters_mass.at(idau));
                    tau_gen_tmp += tmp4V;
                    if (gen_particle_daughters_charge.at(idau) != 0 and gen_particle_daughters_status.at(idau) == 1) tau_gen_nch_tmp ++; // charged particles
                    else if(gen_particle_daughters_charge.at(idau) == 0 and gen_particle_daughters_id.at(idau) == 111) tau_gen_np0_tmp++;
                    else if(gen_particle_daughters_charge.at(idau) == 0 and gen_particle_daughters_id.at(idau) != 111) tau_gen_nnh_tmp++;
                    }
                }	
                if(tau_gen_tmp.Pt() > 0){ // good hadronic tau
                    tau_gen_visible_.push_back(tau_gen_tmp);
                    tau_gen_tmp.SetPtEtaPhiM(gen_particle_pt.at(igen),gen_particle_eta.at(igen),gen_particle_phi.at(igen),gen_particle_mass.at(igen));
                    tau_gen_charge_.push_back((gen_particle_id.at(igen) > 0) ? -1 : 1);
                    // std::cout<<((gen_particle_id.at(igen) > 0) ? -1 : 1)<<std::endl;
                    tau_gen_.push_back(tau_gen_tmp);
                    tau_gen_nch_.push_back(tau_gen_nch_tmp);
                    tau_gen_np0_.push_back(tau_gen_np0_tmp);
                    tau_gen_nnh_.push_back(tau_gen_nnh_tmp);
                }
            }
        }       
    }
}










//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetNTuplizer);
