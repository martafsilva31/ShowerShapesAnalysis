///////////////////////////////////////////////////////////////////////////////
// fill_histograms_layer1.C
//
// Two-pass pipeline for LAYER-1 (EM strip) cell-energy reweighting.
//
// Works on the full 56 × 2 = 112-cell strip cluster. Cell corrections are
// accumulated in coordinates relative to the hot Lr1 strip (eta_rel = i -
// hotEtaIdx) so the profiles are translation-invariant, and re-applied in
// absolute coordinates in pass 2 to recompute shower shapes from the same
// full cluster that xAOD uses. Five shower shapes:
//   weta1, w_s,tot, f_side, ΔE, E_ratio.
//
// Usage:
//   root -l -b -q 'fill_histograms_layer1.C("llgamma","unconverted",
//                      "../../output/Layer_1/eta_loose")'
//   root -l -b -q 'fill_histograms_layer1.C("llgamma","unconverted",
//                      "../../output/Layer_1/eta_pt_tight","eta_pt","tight")'
///////////////////////////////////////////////////////////////////////////////

#include "config_layer1.h"

#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <cmath>

using namespace config_l1;  // includes all of config:: too

// ============================================================
// Shower-shape channel indices and metadata
// ============================================================
const int kNBins = 100;

enum SSL1 { kWeta1 = 0, kWstot = 1, kFside = 2, kDeltaE = 3, kEratio = 4, kNSS = 5 };

static const char* ssName[kNSS]  = {"weta1", "wstot", "fside", "deltae", "eratio"};
static const char* ssTitle[kNSS] = {"w_{#eta 1}", "w_{s,tot}", "f_{side}",
                                    "#DeltaE [MeV]", "E_{ratio}"};

static SSRange ssRange(int s) {
    switch (s) {
        case kWeta1:  return rangeWeta1();
        case kWstot:  return rangeWstot();
        case kFside:  return rangeFside();
        case kDeltaE: return rangeDeltaE();
        case kEratio: return rangeEratio();
    }
    return {0, 1};
}

TH1D* makeSS(int iss, const char* tag, const char* suf = "") {
    TString n = Form("h_%s_%s%s", ssName[iss], tag, suf);
    TString t = Form("%s (%s%s);%s;Events", ssTitle[iss], tag, suf, ssTitle[iss]);
    SSRange r = ssRange(iss);
    return new TH1D(n, t, kNBins, r.lo, r.hi);
}

// ============================================================
// Correction storage — per (eta, pT, rel_cell) with rel_cell
// indexed by kr = (eta_rel + kRelEtaHalf) * kRelPhiSize + phi_rel.
// ============================================================
struct Corrections {
    double delta  [kNEtaBins][kNPtBins][kRelGridSize] = {};  // M1
    double shift  [kNEtaBins][kNPtBins][kRelGridSize] = {};  // M2
    double stretch[kNEtaBins][kNPtBins][kRelGridSize] = {};  // M2
};

// Map an absolute Layer-1 cell (iEta, iPhi) to the relative-grid index
// kr = (iEta - hotEtaIdx + kRelEtaHalf) * kRelPhiSize + iPhi.
static inline int absToRel(int iEta, int iPhi, int hotEtaIdx) {
    int er = iEta - hotEtaIdx;            // in [-kRelEtaHalf, kRelEtaHalf]
    return (er + kRelEtaHalf) * kRelPhiSize + iPhi;
}

// ============================================================
// Main
// ============================================================
void fill_histograms_layer1(const char* channel   = "llgamma",
                            const char* scenario  = "unconverted",
                            const char* baseDir   = "../../output/Layer_1/eta_loose",
                            const char* binning   = "eta",
                            const char* isolation = "loose") {

    TH1::SetDefaultSumw2(true);

    // --------------------------------------------------------
    // Binning mode
    // --------------------------------------------------------
    TString binMode(binning);
    bool usePtBins = (binMode == "eta_pt");
    int nPtUsed = usePtBins ? kNPtBins : 1;

    // --------------------------------------------------------
    // Isolation override
    // --------------------------------------------------------
    Selection sel = getSelection(scenario);
    TString isoMode(isolation);
    if (isoMode == "tight") {
        sel.requireTightIso = true;
        sel.requireLooseIso = false;
    }

    // --------------------------------------------------------
    // Inputs
    // --------------------------------------------------------
    TString inputDir = "/dcache/atlas/mfernand/qt_ntuples/data24/";
    TString dataF, mcF;
    TString ch(channel);
    if      (ch == "eegamma")   { dataF = inputDir + "egam3.root";
                                  mcF   = inputDir + "mc_eegamma.root"; }
    else if (ch == "mumugamma") { dataF = inputDir + "egam4.root";
                                  mcF   = inputDir + "mc_mumugamma.root"; }
    else if (ch == "llgamma")   { dataF = inputDir + "egam3.root," + inputDir + "egam4.root";
                                  mcF   = inputDir + "mc_eegamma.root," + inputDir + "mc_mumugamma.root"; }
    else {
        std::cerr << "ERROR: unknown channel '" << channel << "'\n";
        return;
    }

    // --------------------------------------------------------
    // Output
    // --------------------------------------------------------
    TString outPath = Form("%s/%s/%s", baseDir, channel, scenario);
    gSystem->mkdir(outPath, true);
    TString outFile = outPath + "/histograms.root";
    TFile* fout = TFile::Open(outFile, "RECREATE");

    std::cout << "=== fill_histograms_layer1 ===\n"
              << "  channel   = " << channel   << "\n"
              << "  scenario  = " << scenario  << "\n"
              << "  binning   = " << binning   << " (" << nPtUsed << " pT bins)\n"
              << "  isolation = " << isolation << "\n"
              << "  data      = " << dataF     << "\n"
              << "  mc        = " << mcF       << "\n"
              << "  output    = " << outFile   << "\n";

    // ============================================================
    //  Histograms
    // ============================================================

    // Integrated
    TH1D *h_d_comp[kNSS], *h_d_stor[kNSS];
    TH1D *h_m_comp[kNSS], *h_m_fudg[kNSS], *h_m_unf[kNSS];
    TH1D *h_m_M1[kNSS],   *h_m_M2[kNSS];
    TH1D *h_m_M3[kNSS],   *h_m_M4[kNSS];
    for (int s = 0; s < kNSS; ++s) {
        h_d_comp[s] = makeSS(s, "data_computed");
        h_d_stor[s] = makeSS(s, "data_stored");
        h_m_comp[s] = makeSS(s, "mc_computed");
        h_m_fudg[s] = makeSS(s, "mc_fudged");
        h_m_unf[s]  = makeSS(s, "mc_unfudged");
        h_m_M1[s]   = makeSS(s, "mc_M1");
        h_m_M2[s]   = makeSS(s, "mc_M2");
        h_m_M3[s]   = makeSS(s, "mc_M3");
        h_m_M4[s]   = makeSS(s, "mc_M4");
    }

    // Per-eta
    TH1D *h_d_eta[kNSS][kNEtaBins];
    TH1D *h_m_eta[kNSS][kNEtaBins];
    TH1D *h_m_M1_eta[kNSS][kNEtaBins];
    TH1D *h_m_M2_eta[kNSS][kNEtaBins];
    TH1D *h_m_M3_eta[kNSS][kNEtaBins];
    TH1D *h_m_M4_eta[kNSS][kNEtaBins];
    TH1D *h_d_stor_eta[kNSS][kNEtaBins];
    TH1D *h_m_fudg_eta[kNSS][kNEtaBins];
    TH1D *h_m_unf_eta[kNSS][kNEtaBins];
    fout->cd();
    for (int s = 0; s < kNSS; ++s)
        for (int e = 0; e < kNEtaBins; ++e) {
            TString suf = Form("_eta%02d", e);
            h_d_eta[s][e]      = makeSS(s, "data",         suf.Data());
            h_m_eta[s][e]      = makeSS(s, "mc",           suf.Data());
            h_m_M1_eta[s][e]   = makeSS(s, "mc_M1",        suf.Data());
            h_m_M2_eta[s][e]   = makeSS(s, "mc_M2",        suf.Data());
            h_m_M3_eta[s][e]   = makeSS(s, "mc_M3",        suf.Data());
            h_m_M4_eta[s][e]   = makeSS(s, "mc_M4",        suf.Data());
            h_d_stor_eta[s][e] = makeSS(s, "data_stored",  suf.Data());
            h_m_fudg_eta[s][e] = makeSS(s, "mc_fudged",    suf.Data());
            h_m_unf_eta[s][e]  = makeSS(s, "mc_unfudged",  suf.Data());
        }

    // Per-(eta,pT)
    TH1D *h_d_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_M1_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_M2_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_M3_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_M4_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_d_stor_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_fudg_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_unf_etapt[kNSS][kNEtaBins][kNPtBins];
    if (usePtBins) {
        fout->cd();
        for (int s = 0; s < kNSS; ++s)
            for (int e = 0; e < kNEtaBins; ++e)
                for (int p = 0; p < kNPtBins; ++p) {
                    TString suf = Form("_eta%02d_pt%02d", e, p);
                    h_d_etapt[s][e][p]      = makeSS(s, "data",        suf.Data());
                    h_m_etapt[s][e][p]      = makeSS(s, "mc",          suf.Data());
                    h_m_M1_etapt[s][e][p]   = makeSS(s, "mc_M1",       suf.Data());
                    h_m_M2_etapt[s][e][p]   = makeSS(s, "mc_M2",       suf.Data());
                    h_m_M3_etapt[s][e][p]   = makeSS(s, "mc_M3",       suf.Data());
                    h_m_M4_etapt[s][e][p]   = makeSS(s, "mc_M4",       suf.Data());
                    h_d_stor_etapt[s][e][p] = makeSS(s, "data_stored", suf.Data());
                    h_m_fudg_etapt[s][e][p] = makeSS(s, "mc_fudged",   suf.Data());
                    h_m_unf_etapt[s][e][p]  = makeSS(s, "mc_unfudged", suf.Data());
                }
    }

    // M4 rank-2 displacement distribution (per (eta, pT) bin) — captures
    // the position of the 2nd-highest-energy cell relative to the hot
    // strip.  Used only for M4 (stochastic reshuffle).  Axes: eta_rel in
    // [-kRelEtaHalf, kRelEtaHalf], phi_rel in {0, 1}.  Filled from data in
    // pass 1; also filled from MC for diagnostics.
    TH2D *h_dp2_d[kNEtaBins][kNPtBins] = {};
    TH2D *h_dp2_m[kNEtaBins][kNPtBins] = {};
    for (int e = 0; e < kNEtaBins; ++e) {
        for (int p = 0; p < nPtUsed; ++p) {
            TString suf = usePtBins ? Form("_eta%02d_pt%02d", e, p)
                                    : Form("_eta%02d", e);
            h_dp2_d[e][p] = new TH2D(
                Form("h_dp2_data%s",  suf.Data()),
                Form("Rank-2 #Delta p (data) %s;#Delta#eta (strips);#Delta#phi",
                     suf.Data()),
                kRelEtaSize, -kRelEtaHalf - 0.5, kRelEtaHalf + 0.5,
                kRelPhiSize, -0.5, kRelPhiSize - 0.5);
            h_dp2_m[e][p] = new TH2D(
                Form("h_dp2_mc%s",    suf.Data()),
                Form("Rank-2 #Delta p (MC) %s;#Delta#eta (strips);#Delta#phi",
                     suf.Data()),
                kRelEtaSize, -kRelEtaHalf - 0.5, kRelEtaHalf + 0.5,
                kRelPhiSize, -0.5, kRelPhiSize - 0.5);
        }
    }

    // Event counts per eta bin
    TH1D* h_cnt_d = new TH1D("h_counts_data",
        "Passing events per #eta bin (data);#eta bin;events",
        kNEtaBins, 0, kNEtaBins);
    TH1D* h_cnt_m = new TH1D("h_counts_mc",
        "Passing events per #eta bin (MC, weighted);#eta bin;events",
        kNEtaBins, 0, kNEtaBins);

    // pT / eta distributions
    TH1D* h_pt_d       = new TH1D("h_pt_data",      "Photon p_{T} (data);p_{T} [GeV];Events / 2 GeV", 50, 0, 100);
    TH1D* h_pt_m       = new TH1D("h_pt_mc",        "Photon p_{T} (MC);p_{T} [GeV];Events / 2 GeV",   50, 0, 100);
    TH1D* h_eta_d      = new TH1D("h_eta_data",     "Photon |#eta| (data);|#eta|;Events",              48, 0, 2.4);
    TH1D* h_eta_m      = new TH1D("h_eta_mc",       "Photon |#eta| (MC);|#eta|;Events",                48, 0, 2.4);
    TH1D* h_pt_m_nosf  = new TH1D("h_pt_mc_nosf",   "Photon p_{T} (MC, no SFs);p_{T} [GeV];Events / 2 GeV", 50, 0, 100);
    TH1D* h_eta_m_nosf = new TH1D("h_eta_mc_nosf",  "Photon |#eta| (MC, no SFs);|#eta|;Events",             48, 0, 2.4);

    // Cluster-health diagnostics
    Long64_t nL1Fail_d = 0, nL1Try_d = 0;
    Long64_t nL1Fail_m = 0, nL1Try_m = 0;

    // ============================================================
    //  Accumulators (Pass 1)  — per (eta, pT, rel_cell)
    //
    //  Per-cell TH1Ds are used for sigma (M2) because their explicit
    //  per-cell range windows act as outlier clipping via overflow
    //  exclusion (TH1::GetRMS() excludes under/overflow by default).
    //  Direct Σwf² accumulators are numerically unstable at small
    //  variance (cancellation in E[X²]-μ²) — don't replace without
    //  re-introducing equivalent clipping.
    // ============================================================
    double sf_d [kNEtaBins][kNPtBins][kRelGridSize] = {};
    double cnt_d[kNEtaBins][kNPtBins]               = {};
    double swf_m [kNEtaBins][kNPtBins][kRelGridSize] = {};
    double sw_m  [kNEtaBins][kNPtBins]               = {};

    // Per-cell TH1D for sigma (M2)
    const int kCellBins = 100;
    TH1D* h_frac_d[kNEtaBins][kNPtBins][kRelGridSize];
    TH1D* h_frac_m[kNEtaBins][kNPtBins][kRelGridSize];
    for (int e = 0; e < kNEtaBins; ++e) {
        for (int p = 0; p < nPtUsed; ++p) {
            for (int k = 0; k < kRelGridSize; ++k) {
                // Rough ranges: wider for hot strip row, tighter for far strips.
                int er = k / kRelPhiSize;              // 0..kRelEtaSize-1
                int eta_rel = er - kRelEtaHalf;        // -55..55
                double fLo = -0.05, fHi = 0.25;
                if (eta_rel == 0)                    { fLo = -0.05; fHi = 0.80; }
                else if (std::abs(eta_rel) <= 2)     { fLo = -0.05; fHi = 0.50; }
                else if (std::abs(eta_rel) <= 5)     { fLo = -0.05; fHi = 0.30; }
                else if (std::abs(eta_rel) <= 10)    { fLo = -0.05; fHi = 0.15; }
                else                                 { fLo = -0.05; fHi = 0.05; }

                h_frac_d[e][p][k] = new TH1D(
                    Form("h_frac_d_eta%02d_pt%02d_cell%03d", e, p, k), "",
                    kCellBins, fLo, fHi);
                h_frac_m[e][p][k] = new TH1D(
                    Form("h_frac_m_eta%02d_pt%02d_cell%03d", e, p, k), "",
                    kCellBins, fLo, fHi);
                h_frac_d[e][p][k]->SetDirectory(nullptr);
                h_frac_m[e][p][k]->SetDirectory(nullptr);
            }
        }
    }

    // Post-M2 accumulator for cell profiles
    double swf_m2[kNEtaBins][kNPtBins][kRelGridSize] = {};
    double sw_m2 [kNEtaBins][kNPtBins]               = {};

    // ============================================================
    //  Branches
    // ============================================================
    Float_t phPt, phEta, phPhi;
    Float_t clE, clEta2, clEtaMax2;
    Double_t vmll, vmllg, vdrl1, vdrl2, vdrll;
    Bool_t  visconv, vtight, visoL, visoT, vtruth;
    Int_t   vcellSz;  // Lr2 cluster size (used by passSelection as quality cut)
    Float_t vmcw, vmuw;
    UInt_t  vmcch;
    Double_t vl1sf, vl2sf, vphIsosf;

    // Layer-1 stored shower shapes
    Float_t vWeta1, vWstot1, vFside, vDeltaE, vEratio;
    Float_t vWeta1U, vWstot1U, vFsideU, vDeltaEU, vEratioU;

    // Layer-1 cell arrays
    std::vector<double>* vL1E   = nullptr;
    std::vector<double>* vL1Eta = nullptr;
    std::vector<double>* vL1Phi = nullptr;
    Int_t vL1Size;

    auto setBrCommon = [&](TTree* t) {
        vL1E = vL1Eta = vL1Phi = nullptr;
        t->SetBranchAddress(kPhotonPtBranch,       &phPt);
        t->SetBranchAddress(kPhotonEtaBranch,      &phEta);
        t->SetBranchAddress(kPhotonPhiBranch,      &phPhi);
        t->SetBranchAddress(kPhotonEBranch,        &clE);
        t->SetBranchAddress(kClusterEta2Branch,    &clEta2);
        t->SetBranchAddress(kClusterEtamax2Branch, &clEtaMax2);
        t->SetBranchAddress(kMllBranch,            &vmll);
        t->SetBranchAddress(kMllgBranch,           &vmllg);
        t->SetBranchAddress(kDRLepton1Branch,      &vdrl1);
        t->SetBranchAddress(kDRLepton2Branch,      &vdrl2);
        t->SetBranchAddress(kDRllBranch,           &vdrll);
        t->SetBranchAddress(kIsConvBranch,         &visconv);
        t->SetBranchAddress(kTightIDBranch,        &vtight);
        t->SetBranchAddress(kIsoLooseBranch,       &visoL);
        t->SetBranchAddress(kIsoTightBranch,       &visoT);
        t->SetBranchAddress(kCellSizeBranch,       &vcellSz);    // Lr2 size
        // Layer-1 cells
        t->SetBranchAddress(kL1CellBranch,         &vL1E);
        t->SetBranchAddress(kL1CellEtaBranch,      &vL1Eta);
        t->SetBranchAddress(kL1CellPhiBranch,      &vL1Phi);
        t->SetBranchAddress(kL1CellSizeBranch,     &vL1Size);
        // Layer-1 stored shower shapes
        t->SetBranchAddress(kWeta1Branch,          &vWeta1);
        t->SetBranchAddress(kWstot1Branch,         &vWstot1);
        t->SetBranchAddress(kFsideBranch,          &vFside);
        t->SetBranchAddress(kDeltaEBranch,         &vDeltaE);
        t->SetBranchAddress(kEratioBranch,         &vEratio);
    };

    auto setBrMC = [&](TTree* t) {
        setBrCommon(t);
        t->SetBranchAddress(kMCWeightBranch,    &vmcw);
        t->SetBranchAddress(kPUWeightBranch,    &vmuw);
        t->SetBranchAddress("mcChannelNumber",  &vmcch);
        t->SetBranchAddress(kTruthMatchBranch,  &vtruth);
        t->SetBranchAddress(kWeta1UnfBranch,    &vWeta1U);
        t->SetBranchAddress(kWstot1UnfBranch,   &vWstot1U);
        t->SetBranchAddress(kFsideUnfBranch,    &vFsideU);
        t->SetBranchAddress(kDeltaEUnfBranch,   &vDeltaEU);
        t->SetBranchAddress(kEratioUnfBranch,   &vEratioU);
        t->SetBranchAddress(kLepton1SFBranch,   &vl1sf);
        t->SetBranchAddress(kLepton2SFBranch,   &vl2sf);
        t->SetBranchAddress(kPhotonIsoSFBranch, &vphIsosf);
    };

    auto getPtBin = [&](double pt) -> int {
        if (!usePtBins) return 0;
        return findPtBin(pt);
    };

    // Helper: compute all 5 shower shapes from the 112-cell full grid.
    auto computeAll = [&](const std::vector<double>& full, int hotE,
                          const std::vector<double>& cEta, double etaMaxL2,
                          double clusterEta1, double hotEta1,
                          double out[kNSS]) {
        out[kWeta1]  = calcWeta1(full, hotE, cEta, clusterEta1, hotEta1);
        out[kWstot]  = calcWsTot(full, hotE, cEta, etaMaxL2);
        out[kFside]  = calcFside(full, hotE, cEta);
        out[kDeltaE] = calcDeltaE(full, cEta);
        out[kEratio] = calcEratio(full, cEta);
    };

    // ============================================================
    //  PASS 1 — DATA
    // ============================================================
    {
        TChain* td = makeChain(dataF);
        setBrCommon(td);
        Long64_t N = td->GetEntries();
        Long64_t nPass = 0;
        std::cout << "\n--- Pass 1: Data (" << N << " entries) ---\n";

        for (Long64_t i = 0; i < N; ++i) {
            td->GetEntry(i);
            if (i % 500000 == 0) std::cout << "  " << i << " / " << N << "\n";

            if (!passSelection(sel, phPt, clEta2,
                               vmll, vmllg, vdrl1, vdrl2, vdrll, vcellSz,
                               visconv, vtight, visoL, visoT,
                               false, false))
                continue;

            ++nL1Try_d;
            std::vector<double> cE  (vL1E->begin(),   vL1E->end());
            std::vector<double> cEta(vL1Eta->begin(), vL1Eta->end());
            std::vector<double> cPhi(vL1Phi->begin(), vL1Phi->end());

            std::vector<double> full;
            int hotE = -1, hotP = -1;
            if (!buildFullGrid(cE, cEta, cPhi, full, hotE, hotP)) {
                ++nL1Fail_d;
                continue;
            }

            double Etot = 0;
            for (int k = 0; k < kL1GridSize; ++k) Etot += full[k];
            if (Etot <= 0) continue;

            h_pt_d ->Fill(phPt);
            h_eta_d->Fill(std::fabs(clEta2));

            double clusterEta1 = 0.0, sumE_l1 = 0.0;
            for (int kk = 0; kk < kL1GridSize; ++kk) {
                if (isPadded(cEta[kk], cPhi[kk])) continue;
                clusterEta1 += full[kk] * cEta[kk];
                sumE_l1     += full[kk];
            }
            clusterEta1 = (sumE_l1 > 0) ? clusterEta1 / sumE_l1 : 0.0;
            double hotEta1 = cEta[hotE * kL1PhiSize];

            double ss_c[kNSS];
            computeAll(full, hotE, cEta, clEtaMax2, clusterEta1, hotEta1, ss_c);
            double ss_s[kNSS] = {vWeta1, vWstot1, vFside, vDeltaE, vEratio};

            for (int s = 0; s < kNSS; ++s) {
                h_d_comp[s]->Fill(ss_c[s]);
                h_d_stor[s]->Fill(ss_s[s]);
            }

            int eb = findEtaBin(std::fabs(clEta2));
            int pb = getPtBin(phPt);
            if (eb < 0 || pb < 0) continue;

            for (int s = 0; s < kNSS; ++s) {
                h_d_eta     [s][eb]->Fill(ss_c[s]);
                h_d_stor_eta[s][eb]->Fill(ss_s[s]);
                if (usePtBins) {
                    h_d_etapt     [s][eb][pb]->Fill(ss_c[s]);
                    h_d_stor_etapt[s][eb][pb]->Fill(ss_s[s]);
                }
            }

            // Accumulate cell fractions in relative coordinates.
            for (int i = 0; i < kL1EtaSize; ++i) {
                for (int j = 0; j < kL1PhiSize; ++j) {
                    int kAbs = i * kL1PhiSize + j;
                    int kRel = absToRel(i, j, hotE);
                    double fk = full[kAbs] / Etot;
                    sf_d[eb][pb][kRel] += fk;
                    h_frac_d[eb][pb][kRel]->Fill(fk);
                }
            }
            cnt_d[eb][pb] += 1.0;

            // Rank-2 position (relative to hot strip) — fed to M4.
            {
                int kM1 = -1, kM2 = -1;
                if (findTop2Abs(full, hotE, hotP, kM1, kM2)) {
                    int i2 = kM2 / kL1PhiSize;
                    int j2 = kM2 % kL1PhiSize;
                    h_dp2_d[eb][pb]->Fill(i2 - hotE, j2);
                }
            }

            h_cnt_d->Fill(eb);
            ++nPass;
        }
        std::cout << "  Data passing: " << nPass
                  << "  (Lr1 unhealthy " << nL1Fail_d << " / " << nL1Try_d << ")\n";
        delete td;
    }

    // ============================================================
    //  PASS 1 — MC
    // ============================================================
    auto sumWmap = computeSumWeightsFromFiles(mcF);
    {
        TChain* tm = makeChain(mcF);
        setBrMC(tm);
        Long64_t N = tm->GetEntries();
        Long64_t nPass = 0;
        std::cout << "\n--- Pass 1: MC (" << N << " entries) ---\n";

        for (Long64_t i = 0; i < N; ++i) {
            tm->GetEntry(i);
            if (i % 500000 == 0) std::cout << "  " << i << " / " << N << "\n";

            if (!passSelection(sel, phPt, clEta2,
                               vmll, vmllg, vdrl1, vdrl2, vdrll, vcellSz,
                               visconv, vtight, visoL, visoT,
                               (bool)vtruth, true))
                continue;

            ++nL1Try_m;
            std::vector<double> cE  (vL1E->begin(),   vL1E->end());
            std::vector<double> cEta(vL1Eta->begin(), vL1Eta->end());
            std::vector<double> cPhi(vL1Phi->begin(), vL1Phi->end());

            std::vector<double> full;
            int hotE = -1, hotP = -1;
            if (!buildFullGrid(cE, cEta, cPhi, full, hotE, hotP)) {
                ++nL1Fail_m;
                continue;
            }

            double nf = mcNormFactor(vmcch, sumWmap);
            double w  = mcWeight(vmcw, vmuw, nf, vl1sf, vl2sf, vphIsosf);
            double w_nosf = mcWeight(vmcw, vmuw, nf, 1.0, 1.0, 1.0);
            double Etot = 0;
            for (int k = 0; k < kL1GridSize; ++k) Etot += full[k];
            if (Etot <= 0) continue;

            h_pt_m      ->Fill(phPt,              w);
            h_eta_m     ->Fill(std::fabs(clEta2), w);
            h_pt_m_nosf ->Fill(phPt,              w_nosf);
            h_eta_m_nosf->Fill(std::fabs(clEta2), w_nosf);

            double clusterEta1 = 0.0, sumE_l1 = 0.0;
            for (int kk = 0; kk < kL1GridSize; ++kk) {
                if (isPadded(cEta[kk], cPhi[kk])) continue;
                clusterEta1 += full[kk] * cEta[kk];
                sumE_l1     += full[kk];
            }
            clusterEta1 = (sumE_l1 > 0) ? clusterEta1 / sumE_l1 : 0.0;
            double hotEta1 = cEta[hotE * kL1PhiSize];

            double ss_c[kNSS];
            computeAll(full, hotE, cEta, clEtaMax2, clusterEta1, hotEta1, ss_c);
            double ss_f[kNSS] = {vWeta1,  vWstot1,  vFside,  vDeltaE,  vEratio};
            double ss_u[kNSS] = {vWeta1U, vWstot1U, vFsideU, vDeltaEU, vEratioU};

            for (int s = 0; s < kNSS; ++s) {
                h_m_comp[s]->Fill(ss_c[s], w);
                h_m_fudg[s]->Fill(ss_f[s], w);
                h_m_unf [s]->Fill(ss_u[s], w);
            }

            int eb = findEtaBin(std::fabs(clEta2));
            int pb = getPtBin(phPt);
            if (eb < 0 || pb < 0) continue;

            for (int s = 0; s < kNSS; ++s) {
                h_m_eta     [s][eb]->Fill(ss_c[s], w);
                h_m_fudg_eta[s][eb]->Fill(ss_f[s], w);
                h_m_unf_eta [s][eb]->Fill(ss_u[s], w);
                if (usePtBins) {
                    h_m_etapt     [s][eb][pb]->Fill(ss_c[s], w);
                    h_m_fudg_etapt[s][eb][pb]->Fill(ss_f[s], w);
                    h_m_unf_etapt [s][eb][pb]->Fill(ss_u[s], w);
                }
            }

            for (int i = 0; i < kL1EtaSize; ++i) {
                for (int j = 0; j < kL1PhiSize; ++j) {
                    int kAbs = i * kL1PhiSize + j;
                    int kRel = absToRel(i, j, hotE);
                    double fk = full[kAbs] / Etot;
                    swf_m[eb][pb][kRel] += w * fk;
                    h_frac_m[eb][pb][kRel]->Fill(fk, w);
                }
            }
            sw_m[eb][pb] += w;

            // Rank-2 position (MC, for diagnostics only — M4 samples from data).
            {
                int kM1 = -1, kM2 = -1;
                if (findTop2Abs(full, hotE, hotP, kM1, kM2)) {
                    int i2 = kM2 / kL1PhiSize;
                    int j2 = kM2 % kL1PhiSize;
                    h_dp2_m[eb][pb]->Fill(i2 - hotE, j2, w);
                }
            }

            h_cnt_m->Fill(eb, w);
            ++nPass;
        }
        std::cout << "  MC passing: " << nPass
                  << "  (Lr1 unhealthy " << nL1Fail_m << " / " << nL1Try_m << ")\n";
        delete tm;
    }

    std::cout << "\n--- Layer-1 cluster-health summary ---\n";
    std::cout << Form("  Data: %lld / %lld fail Lr1 health (%.2f%%)\n",
                      nL1Fail_d, nL1Try_d,
                      nL1Try_d > 0 ? 100.0 * nL1Fail_d / nL1Try_d : 0.0);
    std::cout << Form("  MC:   %lld / %lld fail Lr1 health (%.2f%%)\n",
                      nL1Fail_m, nL1Try_m,
                      nL1Try_m > 0 ? 100.0 * nL1Fail_m / nL1Try_m : 0.0);

    // ============================================================
    //  COMPUTE CORRECTIONS
    // ============================================================
    Corrections corr;
    fout->mkdir("corrections")->cd();

    std::cout << "\n--- Computing corrections ---\n";
    for (int e = 0; e < kNEtaBins; ++e) {
        for (int p = 0; p < nPtUsed; ++p) {
            TString binLabel;
            if (usePtBins)
                binLabel = Form("  eta %2d [%.2f,%.2f) pT %d [%.0f,%.0f)  data=%6.0f  MC_wt=%.1f",
                    e, kEtaLimits[e], kEtaLimits[e+1],
                    p, kPtLimits[p], kPtLimits[p+1],
                    cnt_d[e][p], sw_m[e][p]);
            else
                binLabel = Form("  eta %2d [%.2f,%.2f)  data=%6.0f  MC_wt=%.1f",
                    e, kEtaLimits[e], kEtaLimits[e+1],
                    cnt_d[e][p], sw_m[e][p]);
            std::cout << binLabel << "\n";

            if (cnt_d[e][p] < 10 || sw_m[e][p] < 10) {
                std::cout << "    -> too few events, skipping corrections\n";
                continue;
            }

            TString suffix = usePtBins ? Form("_eta%02d_pt%02d", e, p)
                                       : Form("_eta%02d", e);

            TH1D* hd = new TH1D(Form("h_delta%s",   suffix.Data()),
                Form("M1 #Delta  %s;cell;#Delta", suffix.Data()),
                kRelGridSize, 0, kRelGridSize);
            TH1D* hs = new TH1D(Form("h_shift%s",   suffix.Data()),
                Form("M2 shift  %s;cell;shift",   suffix.Data()),
                kRelGridSize, 0, kRelGridSize);
            TH1D* ht = new TH1D(Form("h_stretch%s", suffix.Data()),
                Form("M2 stretch  %s;cell;stretch", suffix.Data()),
                kRelGridSize, 0, kRelGridSize);

            for (int k = 0; k < kRelGridSize; ++k) {
                double mu_d  = sf_d[e][p][k]  / cnt_d[e][p];
                double mu_m  = swf_m[e][p][k] / sw_m[e][p];
                double sig_d = h_frac_d[e][p][k]->GetRMS();
                double sig_m = h_frac_m[e][p][k]->GetRMS();

                corr.delta[e][p][k] = mu_d - mu_m;
                double str = (sig_m > 1e-12) ? sig_d / sig_m : 1.0;
                corr.stretch[e][p][k] = str;
                corr.shift[e][p][k]   = mu_d - str * mu_m;

                hd->SetBinContent(k + 1, corr.delta[e][p][k]);
                hs->SetBinContent(k + 1, corr.shift[e][p][k]);
                ht->SetBinContent(k + 1, corr.stretch[e][p][k]);
            }
        }
    }

    // ============================================================
    //  PASS 2 — MC  (apply corrections)
    // ============================================================
    fout->cd();
    // M4 uses a stochastic per-event draw from the data rank-2 position
    // distribution.  Seed gRandom with a fixed value so runs are exactly
    // reproducible (any change in MC ordering will change the output,
    // but re-running on the same input is deterministic).
    delete gRandom;
    gRandom = new TRandom3(42u);
    {
        TChain* tm = makeChain(mcF);
        setBrMC(tm);
        Long64_t N = tm->GetEntries();
        Long64_t nPass = 0;
        std::cout << "\n--- Pass 2: MC (" << N << " entries) ---\n";

        for (Long64_t i = 0; i < N; ++i) {
            tm->GetEntry(i);
            if (i % 500000 == 0) std::cout << "  " << i << " / " << N << "\n";

            if (!passSelection(sel, phPt, clEta2,
                               vmll, vmllg, vdrl1, vdrl2, vdrll, vcellSz,
                               visconv, vtight, visoL, visoT,
                               (bool)vtruth, true))
                continue;

            std::vector<double> cE  (vL1E->begin(),   vL1E->end());
            std::vector<double> cEta(vL1Eta->begin(), vL1Eta->end());
            std::vector<double> cPhi(vL1Phi->begin(), vL1Phi->end());

            std::vector<double> full;
            int hotE = -1, hotP = -1;
            if (!buildFullGrid(cE, cEta, cPhi, full, hotE, hotP)) continue;

            double nf = mcNormFactor(vmcch, sumWmap);
            double w  = mcWeight(vmcw, vmuw, nf, vl1sf, vl2sf, vphIsosf);
            double Etot = 0;
            for (int k = 0; k < kL1GridSize; ++k) Etot += full[k];
            if (Etot <= 0) continue;

            double clusterEta1 = 0.0, sumE_l1 = 0.0;
            for (int kk = 0; kk < kL1GridSize; ++kk) {
                if (isPadded(cEta[kk], cPhi[kk])) continue;
                clusterEta1 += full[kk] * cEta[kk];
                sumE_l1     += full[kk];
            }
            clusterEta1 = (sumE_l1 > 0) ? clusterEta1 / sumE_l1 : 0.0;
            double hotEta1 = cEta[hotE * kL1PhiSize];

            int eb = findEtaBin(std::fabs(clEta2));
            int pb = getPtBin(phPt);
            if (eb < 0 || pb < 0) continue;

            // M1: additive delta in fraction space, applied per absolute cell
            // via its relative-grid index.
            std::vector<double> c1(kL1GridSize, 0.0);
            for (int i = 0; i < kL1EtaSize; ++i) {
                for (int j = 0; j < kL1PhiSize; ++j) {
                    int kAbs = i * kL1PhiSize + j;
                    int kRel = absToRel(i, j, hotE);
                    c1[kAbs] = full[kAbs] + corr.delta[eb][pb][kRel] * Etot;
                }
            }

            double ss1[kNSS];
            computeAll(c1, hotE, cEta, clEtaMax2, clusterEta1, hotEta1, ss1);
            for (int s = 0; s < kNSS; ++s) {
                h_m_M1[s]         ->Fill(ss1[s], w);
                h_m_M1_eta[s][eb] ->Fill(ss1[s], w);
                if (usePtBins) h_m_M1_etapt[s][eb][pb]->Fill(ss1[s], w);
            }

            // M2: shift + stretch, rescale to preserve Etot.
            std::vector<double> c2(kL1GridSize, 0.0);
            for (int i = 0; i < kL1EtaSize; ++i) {
                for (int j = 0; j < kL1PhiSize; ++j) {
                    int kAbs = i * kL1PhiSize + j;
                    int kRel = absToRel(i, j, hotE);
                    c2[kAbs] = Etot * corr.shift[eb][pb][kRel]
                             + corr.stretch[eb][pb][kRel] * full[kAbs];
                }
            }

            double Ecorr = 0;
            for (int k = 0; k < kL1GridSize; ++k) Ecorr += c2[k];
            if (Ecorr > 0) {
                double sc = Etot / Ecorr;
                for (int k = 0; k < kL1GridSize; ++k) c2[k] *= sc;
            }

            double ss2[kNSS];
            computeAll(c2, hotE, cEta, clEtaMax2, clusterEta1, hotEta1, ss2);
            for (int s = 0; s < kNSS; ++s) {
                h_m_M2[s]         ->Fill(ss2[s], w);
                h_m_M2_eta[s][eb] ->Fill(ss2[s], w);
                if (usePtBins) h_m_M2_etapt[s][eb][pb]->Fill(ss2[s], w);
            }

            // Post-M2 cell profile in relative coordinates.
            for (int i = 0; i < kL1EtaSize; ++i) {
                for (int j = 0; j < kL1PhiSize; ++j) {
                    int kAbs = i * kL1PhiSize + j;
                    int kRel = absToRel(i, j, hotE);
                    double fk = c2[kAbs] / Etot;
                    swf_m2[eb][pb][kRel] += w * fk;
                }
            }
            sw_m2[eb][pb] += w;

            // ========================================================
            //  M3 — per-cell quantile transport.
            //
            //  Map each MC cell fraction through F_data^{-1}(F_mc(f_k)).
            //  Falls back to identity when either CDF has too few counts
            //  (see applyQuantileTransport in config_layer1.h).  After
            //  transport, rescale the cluster so Etot is preserved.
            // ========================================================
            std::vector<double> c3(kL1GridSize, 0.0);
            for (int i = 0; i < kL1EtaSize; ++i) {
                for (int j = 0; j < kL1PhiSize; ++j) {
                    int kAbs = i * kL1PhiSize + j;
                    int kRel = absToRel(i, j, hotE);
                    double fk  = full[kAbs] / Etot;
                    double fk3 = applyQuantileTransport(
                        h_frac_m[eb][pb][kRel],
                        h_frac_d[eb][pb][kRel], fk);
                    // Clip non-negative (inverse CDF can undershoot into
                    // the [qLo,qHi] clip window if data has a <0 bin).
                    if (fk3 < 0) fk3 = 0;
                    c3[kAbs] = fk3 * Etot;
                }
            }
            double Ec3 = 0;
            for (int k = 0; k < kL1GridSize; ++k) Ec3 += c3[k];
            if (Ec3 > 0) {
                double sc3 = Etot / Ec3;
                for (int k = 0; k < kL1GridSize; ++k) c3[k] *= sc3;
            }

            double ss3[kNSS];
            computeAll(c3, hotE, cEta, clEtaMax2, clusterEta1, hotEta1, ss3);
            for (int s = 0; s < kNSS; ++s) {
                h_m_M3[s]         ->Fill(ss3[s], w);
                h_m_M3_eta[s][eb] ->Fill(ss3[s], w);
                if (usePtBins) h_m_M3_etapt[s][eb][pb]->Fill(ss3[s], w);
            }

            // ========================================================
            //  M4 — rank-2 position reshuffle on top of M3.
            //
            //  The current rank-2 cell (after M3) is found on the full
            //  grid.  A new (Δη, Δφ) offset is drawn from the data Δp
            //  distribution in this (eta, pT) bin.  If the draw moves
            //  rank-2 to a new position, swap energies between the old
            //  and new positions.  The hot cell is excluded from the
            //  draw to keep the rank-1 position fixed at (0, hotPhi).
            // ========================================================
            std::vector<double> c4(c3);  // starts from M3
            if (h_dp2_d[eb][pb] && h_dp2_d[eb][pb]->GetEntries() >= kMinCDFCounts) {
                int kM1 = -1, kM2 = -1;
                if (findTop2Abs(c4, hotE, hotP, kM1, kM2)) {
                    // Current rank-2 relative position.
                    int i2cur = kM2 / kL1PhiSize;
                    int j2cur = kM2 % kL1PhiSize;
                    int dEtaCur = i2cur - hotE;

                    // Draw a new (Δη, Δφ) from the data Δp distribution.
                    // Retry a few times if the draw lands on the hot
                    // cell (Δη=0, Δφ=hotP); fall back to the current
                    // position if the draw keeps failing.
                    double dEtaD = dEtaCur, dPhiD = j2cur;
                    bool drawn = false;
                    for (int attempt = 0; attempt < 5; ++attempt) {
                        h_dp2_d[eb][pb]->GetRandom2(dEtaD, dPhiD);
                        int dEtaI = (int)std::lround(dEtaD);
                        int dPhiI = (int)std::lround(dPhiD);
                        if (dEtaI == 0 && dPhiI == hotP) continue;  // hot cell
                        int iTgt = hotE + dEtaI;
                        int jTgt = dPhiI;
                        if (iTgt < 0 || iTgt >= kL1EtaSize) continue;
                        if (jTgt < 0 || jTgt >= kL1PhiSize) continue;
                        int kTgt = iTgt * kL1PhiSize + jTgt;
                        if (kTgt != kM2) {
                            // Swap energies between current rank-2 and
                            // the new target position.
                            std::swap(c4[kM2], c4[kTgt]);
                        }
                        drawn = true;
                        break;
                    }
                    (void)drawn;  // no-op if draw kept failing
                }
            }

            double ss4[kNSS];
            computeAll(c4, hotE, cEta, clEtaMax2, clusterEta1, hotEta1, ss4);
            for (int s = 0; s < kNSS; ++s) {
                h_m_M4[s]         ->Fill(ss4[s], w);
                h_m_M4_eta[s][eb] ->Fill(ss4[s], w);
                if (usePtBins) h_m_M4_etapt[s][eb][pb]->Fill(ss4[s], w);
            }

            ++nPass;
        }
        std::cout << "  MC pass 2 passing: " << nPass << "\n";
        delete tm;
    }

    // ============================================================
    //  CELL PROFILES — TH2D heatmaps on 21×2 grid
    // ============================================================
    fout->mkdir("cell_profiles")->cd();
    for (int e = 0; e < kNEtaBins; ++e) {
        double tot_cnt_d = 0, tot_sw_m = 0, tot_sw_m2 = 0;
        double tot_sf_d[kRelGridSize]   = {};
        double tot_swf_m[kRelGridSize]  = {};
        double tot_swf_m2[kRelGridSize] = {};
        for (int p = 0; p < nPtUsed; ++p) {
            tot_cnt_d += cnt_d[e][p];
            tot_sw_m  += sw_m[e][p];
            tot_sw_m2 += sw_m2[e][p];
            for (int k = 0; k < kRelGridSize; ++k) {
                tot_sf_d[k]   += sf_d[e][p][k];
                tot_swf_m[k]  += swf_m[e][p][k];
                tot_swf_m2[k] += swf_m2[e][p][k];
            }
        }

        if (tot_cnt_d < 1 && tot_sw_m < 1) continue;

        auto mk2D = [&](const char* tag, const char* ptSuf = "") {
            return new TH2D(
                Form("h_frac_%s_eta%02d%s", tag, e, ptSuf),
                Form("Cell frac %s  #eta [%.2f,%.2f)%s;#eta idx (rel);#phi idx",
                     tag, kEtaLimits[e], kEtaLimits[e + 1], ptSuf),
                kRelEtaSize, 0, kRelEtaSize, kRelPhiSize, 0, kRelPhiSize);
        };

        TH2D* h2_mean_d  = mk2D("mean_data");
        TH2D* h2_mean_m  = mk2D("mean_mc");
        TH2D* h2_mean_m2 = mk2D("mean_mc_M2");
        TH2D* h2_delta   = mk2D("delta");
        TH2D* h2_rms_d   = mk2D("rms_data");
        TH2D* h2_rms_m   = mk2D("rms_mc");

        for (int k = 0; k < kRelGridSize; ++k) {
            int ie = k / kRelPhiSize;
            int ip = k % kRelPhiSize;

            double md  = (tot_cnt_d > 0) ? tot_sf_d[k]   / tot_cnt_d : 0;
            double mm  = (tot_sw_m  > 0) ? tot_swf_m[k]  / tot_sw_m  : 0;
            double mm2 = (tot_sw_m2 > 0) ? tot_swf_m2[k] / tot_sw_m2 : 0;

            double rms_d_sum2 = 0, rms_m_sum2 = 0, rms_d_w = 0, rms_m_w = 0;
            for (int p = 0; p < nPtUsed; ++p) {
                double rd = h_frac_d[e][p][k]->GetRMS();
                double rm = h_frac_m[e][p][k]->GetRMS();
                double nd = cnt_d[e][p];
                double nm = sw_m[e][p];
                rms_d_sum2 += nd * rd * rd;
                rms_m_sum2 += nm * rm * rm;
                rms_d_w += nd;
                rms_m_w += nm;
            }
            double rd = (rms_d_w > 0) ? std::sqrt(rms_d_sum2 / rms_d_w) : 0;
            double rm = (rms_m_w > 0) ? std::sqrt(rms_m_sum2 / rms_m_w) : 0;

            h2_mean_d ->SetBinContent(ie + 1, ip + 1, md);
            h2_mean_m ->SetBinContent(ie + 1, ip + 1, mm);
            h2_mean_m2->SetBinContent(ie + 1, ip + 1, mm2);
            h2_delta  ->SetBinContent(ie + 1, ip + 1, md - mm);
            h2_rms_d  ->SetBinContent(ie + 1, ip + 1, rd);
            h2_rms_m  ->SetBinContent(ie + 1, ip + 1, rm);
        }

        if (usePtBins) {
            for (int p = 0; p < kNPtBins; ++p) {
                if (cnt_d[e][p] < 1 && sw_m[e][p] < 1) continue;
                TString ptSuf = Form("_pt%02d", p);

                TH2D* h2p_mean_d  = mk2D("mean_data",  ptSuf.Data());
                TH2D* h2p_mean_m  = mk2D("mean_mc",    ptSuf.Data());
                TH2D* h2p_mean_m2 = mk2D("mean_mc_M2", ptSuf.Data());
                TH2D* h2p_delta   = mk2D("delta",      ptSuf.Data());
                TH2D* h2p_rms_d   = mk2D("rms_data",   ptSuf.Data());
                TH2D* h2p_rms_m   = mk2D("rms_mc",     ptSuf.Data());

                for (int k = 0; k < kRelGridSize; ++k) {
                    int ie = k / kRelPhiSize;
                    int ip = k % kRelPhiSize;

                    double md  = (cnt_d[e][p] > 0) ? sf_d[e][p][k]   / cnt_d[e][p] : 0;
                    double mm  = (sw_m[e][p]  > 0) ? swf_m[e][p][k]  / sw_m[e][p]  : 0;
                    double mm2 = (sw_m2[e][p] > 0) ? swf_m2[e][p][k] / sw_m2[e][p] : 0;
                    double rd  = h_frac_d[e][p][k]->GetRMS();
                    double rm  = h_frac_m[e][p][k]->GetRMS();

                    h2p_mean_d ->SetBinContent(ie + 1, ip + 1, md);
                    h2p_mean_m ->SetBinContent(ie + 1, ip + 1, mm);
                    h2p_mean_m2->SetBinContent(ie + 1, ip + 1, mm2);
                    h2p_delta  ->SetBinContent(ie + 1, ip + 1, md - mm);
                    h2p_rms_d  ->SetBinContent(ie + 1, ip + 1, rd);
                    h2p_rms_m  ->SetBinContent(ie + 1, ip + 1, rm);
                }
            }
        }
    }

    // ============================================================
    //  WRITE & CLOSE
    // ============================================================
    fout->cd();
    fout->Write();
    fout->Close();
    std::cout << "\n=== Done: " << outFile << " ===\n";
}

#ifndef __CLING__
#include <TROOT.h>
#include <TApplication.h>
int main(int argc, char** argv) {
    TApplication app("app", nullptr, nullptr);
    gROOT->SetBatch(true);
    const char* ch  = (argc > 1) ? argv[1] : "llgamma";
    const char* sc  = (argc > 2) ? argv[2] : "unconverted";
    const char* bd  = (argc > 3) ? argv[3] : "../../output/Layer_1/eta_loose";
    const char* bn  = (argc > 4) ? argv[4] : "eta";
    const char* iso = (argc > 5) ? argv[5] : "loose";
    fill_histograms_layer1(ch, sc, bd, bn, iso);
    return 0;
}
#endif
