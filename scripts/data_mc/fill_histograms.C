///////////////////////////////////////////////////////////////////////////////
// fill_histograms.C
//
// Two-pass pipeline for cell-energy reweighting corrections.
//
//   Pass 1: Loop data + MC -> accumulate per-cell fraction statistics,
//           fill uncorrected shower-shape histograms.
//   Pass 2: Loop MC only   -> apply M1 (flat shift) and M2 (shift+stretch)
//           corrections, fill corrected shower-shape histograms.
//
// Parameters:
//   channel   "eegamma", "mumugamma", or "llgamma"
//   scenario  conversion scenario: "unconverted", "converted", etc.
//   baseDir   output base directory
//   binning   "eta" (14 eta bins, default) or "eta_pt" (14 eta x 6 pT bins)
//   isolation "loose" (default) or "tight"
//
// Usage:
//   root -l -b -q 'fill_histograms.C("llgamma", "unconverted", "../../output/Layer_2/eta_loose")'
//   root -l -b -q 'fill_histograms.C("llgamma", "unconverted", "../../output/Layer_2/eta_pt_tight", "eta_pt", "tight")'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <cmath>

using namespace config;

// ============================================================
// Shower-shape histogram helpers
// ============================================================
const int kNBins = 100;
const bool kExcludeZeroPadFromCorr = false;
enum SS { kReta = 0, kRphi = 1, kWeta2 = 2, kNSS = 3 };
static const char* ssName[kNSS]  = {"reta", "rphi", "weta2"};
static const char* ssTitle[kNSS] = {"R_{#eta}", "R_{#phi}", "w_{#eta 2}"};
static const double ssLo[kNSS]   = {0.80, 0.50, 0.005};
static const double ssHi[kNSS]   = {1.05, 1.10, 0.020};

TH1D* makeSS(int iss, const char* tag, const char* suf = "") {
    TString n = Form("h_%s_%s%s", ssName[iss], tag, suf);
    TString t = Form("%s (%s%s);%s;Events", ssTitle[iss], tag, suf, ssTitle[iss]);
    return new TH1D(n, t, kNBins, ssLo[iss], ssHi[iss]);
}

// ============================================================
// Correction storage — supports 2D (eta, pT) binning
// ============================================================
struct Corrections {
    double delta  [kNEtaBins][kNPtBins][kClusterSize] = {};   // M1
    double shift  [kNEtaBins][kNPtBins][kClusterSize] = {};   // M2
    double stretch[kNEtaBins][kNPtBins][kClusterSize] = {};   // M2
};

// ============================================================
// Main
// ============================================================
void fill_histograms(const char* channel   = "eegamma",
                     const char* scenario  = "unconverted",
                     const char* baseDir   = "../../output/Layer_2/eta_loose",
                     const char* binning   = "eta",
                     const char* isolation = "loose") {

    TH1::SetDefaultSumw2(true);

    // --------------------------------------------------------
    // Determine binning mode
    // --------------------------------------------------------
    TString binMode(binning);
    bool usePtBins = (binMode == "eta_pt");
    int nPtUsed = usePtBins ? kNPtBins : 1;  // 1 = eta-only mode

    // --------------------------------------------------------
    // Apply isolation override to selection
    // --------------------------------------------------------
    Selection sel = getSelection(scenario);
    TString isoMode(isolation);
    if (isoMode == "tight") {
        sel.requireTightIso = true;
        sel.requireLooseIso = false;
    }
    // else: loose is the default from getSelection

    // --------------------------------------------------------
    // Input files
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

    std::cout << "=== fill_histograms ===\n"
              << "  channel   = " << channel   << "\n"
              << "  scenario  = " << scenario  << "\n"
              << "  binning   = " << binning   << " (" << nPtUsed << " pT bins)\n"
              << "  isolation = " << isolation  << "\n"
              << "  data      = " << dataF     << "\n"
              << "  mc        = " << mcF       << "\n"
              << "  output    = " << outFile   << "\n";

    // ============================================================
    //  CREATE HISTOGRAMS
    // ============================================================

    // --- Integrated ---
    TH1D *h_d_comp[kNSS], *h_d_stor[kNSS];
    TH1D *h_m_comp[kNSS], *h_m_fudg[kNSS], *h_m_unf[kNSS];
    TH1D *h_m_M1[kNSS],   *h_m_M2[kNSS];
    for (int s = 0; s < kNSS; ++s) {
        h_d_comp[s] = makeSS(s, "data_computed");
        h_d_stor[s] = makeSS(s, "data_stored");
        h_m_comp[s] = makeSS(s, "mc_computed");
        h_m_fudg[s] = makeSS(s, "mc_fudged");
        h_m_unf[s]  = makeSS(s, "mc_unfudged");
        h_m_M1[s]   = makeSS(s, "mc_M1");
        h_m_M2[s]   = makeSS(s, "mc_M2");
    }

    // --- Per-eta bin ---
    TH1D *h_d_eta[kNSS][kNEtaBins];
    TH1D *h_m_eta[kNSS][kNEtaBins];
    TH1D *h_m_M1_eta[kNSS][kNEtaBins];
    TH1D *h_m_M2_eta[kNSS][kNEtaBins];
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
            h_d_stor_eta[s][e] = makeSS(s, "data_stored",  suf.Data());
            h_m_fudg_eta[s][e] = makeSS(s, "mc_fudged",    suf.Data());
            h_m_unf_eta[s][e]  = makeSS(s, "mc_unfudged",  suf.Data());
        }

    // --- Per-(eta, pT) bin (only created in eta_pt mode) ---
    TH1D *h_d_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_M1_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_M2_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_d_stor_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_fudg_etapt[kNSS][kNEtaBins][kNPtBins];
    TH1D *h_m_unf_etapt[kNSS][kNEtaBins][kNPtBins];
    if (usePtBins) {
        fout->cd();
        for (int s = 0; s < kNSS; ++s)
            for (int e = 0; e < kNEtaBins; ++e)
                for (int p = 0; p < kNPtBins; ++p) {
                    TString suf = Form("_eta%02d_pt%02d", e, p);
                    h_d_etapt[s][e][p]      = makeSS(s, "data",         suf.Data());
                    h_m_etapt[s][e][p]      = makeSS(s, "mc",           suf.Data());
                    h_m_M1_etapt[s][e][p]   = makeSS(s, "mc_M1",        suf.Data());
                    h_m_M2_etapt[s][e][p]   = makeSS(s, "mc_M2",        suf.Data());
                    h_d_stor_etapt[s][e][p] = makeSS(s, "data_stored",  suf.Data());
                    h_m_fudg_etapt[s][e][p] = makeSS(s, "mc_fudged",    suf.Data());
                    h_m_unf_etapt[s][e][p]  = makeSS(s, "mc_unfudged",  suf.Data());
                }
    }

    // --- Event counts per eta bin ---
    TH1D* h_cnt_d = new TH1D("h_counts_data",
        "Passing events per #eta bin (data);#eta bin;events",
        kNEtaBins, 0, kNEtaBins);
    TH1D* h_cnt_m = new TH1D("h_counts_mc",
        "Passing events per #eta bin (MC, weighted);#eta bin;events",
        kNEtaBins, 0, kNEtaBins);

    // --- pT and eta distributions ---
    TH1D* h_pt_d = new TH1D("h_pt_data",
        "Photon p_{T} (data);p_{T} [GeV];Events / 2 GeV", 50, 0, 100);
    TH1D* h_pt_m = new TH1D("h_pt_mc",
        "Photon p_{T} (MC);p_{T} [GeV];Events / 2 GeV", 50, 0, 100);
    TH1D* h_eta_d = new TH1D("h_eta_data",
        "Photon |#eta| (data);|#eta|;Events", 48, 0, 2.4);
    TH1D* h_eta_m = new TH1D("h_eta_mc",
        "Photon |#eta| (MC);|#eta|;Events", 48, 0, 2.4);
    TH1D* h_pt_m_nosf = new TH1D("h_pt_mc_nosf",
        "Photon p_{T} (MC, no SFs);p_{T} [GeV];Events / 2 GeV", 50, 0, 100);
    TH1D* h_eta_m_nosf = new TH1D("h_eta_mc_nosf",
        "Photon |#eta| (MC, no SFs);|#eta|;Events", 48, 0, 2.4);

    // --- Zero-padding diagnostics ---
    Long64_t nZeroPad_d = 0, nTotal_d = 0;
    Long64_t nZeroPad_m = 0, nTotal_m = 0;
    TH1D* h_nzero_d = new TH1D("h_nzero_data",
        "Cells with E=0 per event (data);N_{zero};Events", 78, 0, 78);
    TH1D* h_nzero_m = new TH1D("h_nzero_mc",
        "Cells with E=0 per event (MC);N_{zero};Events", 78, 0, 78);

    // ============================================================
    //  ACCUMULATION ARRAYS  (Pass 1)
    //  Dimensions: [eta][pT][cell]   (pT dim = 1 for eta-only)
    // ============================================================
    double sf_d [kNEtaBins][kNPtBins][kClusterSize] = {};
    double cnt_d[kNEtaBins][kNPtBins]               = {};
    double swf_m [kNEtaBins][kNPtBins][kClusterSize] = {};
    double sw_m  [kNEtaBins][kNPtBins]               = {};

    // M2: per-cell TH1D for histogram-based sigma
    const int kCellBins = 100;
    TH1D* h_frac_d[kNEtaBins][kNPtBins][kClusterSize];
    TH1D* h_frac_m[kNEtaBins][kNPtBins][kClusterSize];
    for (int e = 0; e < kNEtaBins; ++e) {
        for (int p = 0; p < nPtUsed; ++p) {
            for (int k = 0; k < kClusterSize; ++k) {
                double fLo = -0.1, fHi = 0.2;
                if (k == kCentralCell)                                    fHi = 1.0;
                else if (k == 27 || k == 37 || k == 39 || k == 49)       fHi = 0.5;
                else if (k == 26 || k == 28 || k == 48 || k == 50)       fHi = 0.3;

                h_frac_d[e][p][k] = new TH1D(
                    Form("h_frac_d_eta%02d_pt%02d_cell%02d", e, p, k), "",
                    kCellBins, fLo, fHi);
                h_frac_m[e][p][k] = new TH1D(
                    Form("h_frac_m_eta%02d_pt%02d_cell%02d", e, p, k), "",
                    kCellBins, fLo, fHi);
                h_frac_d[e][p][k]->SetDirectory(nullptr);
                h_frac_m[e][p][k]->SetDirectory(nullptr);
            }
        }
    }

    double swf_m2[kNEtaBins][kNPtBins][kClusterSize] = {};
    double sw_m2 [kNEtaBins][kNPtBins]               = {};

    // ============================================================
    //  BRANCH VARIABLES
    // ============================================================
    Float_t phPt, phEta, phPhi;
    Float_t clE, clEta2, clEtaMax2;
    Double_t vmll, vmllg, vdrl1, vdrl2, vdrll;
    Bool_t  visconv, vtight, visoL, visoT, vtruth;
    Int_t   vcellSz;
    Float_t vmcw, vmuw;
    UInt_t  vmcch;
    Double_t vl1sf, vl2sf, vphIsosf;
    Float_t vreta, vrphi, vweta2;
    Float_t vretaU, vrphiU, vweta2U;
    std::vector<double>* vcellE = nullptr;

    auto setBrCommon = [&](TTree* t) {
        vcellE = nullptr;
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
        t->SetBranchAddress(kCellSizeBranch,       &vcellSz);
        t->SetBranchAddress(kCellBranch,           &vcellE);
        t->SetBranchAddress(kRetaBranch,           &vreta);
        t->SetBranchAddress(kRphiBranch,           &vrphi);
        t->SetBranchAddress(kWeta2Branch,          &vweta2);
    };

    auto setBrMC = [&](TTree* t) {
        setBrCommon(t);
        t->SetBranchAddress(kMCWeightBranch,    &vmcw);
        t->SetBranchAddress(kPUWeightBranch,    &vmuw);
        t->SetBranchAddress("mcChannelNumber",  &vmcch);
        t->SetBranchAddress(kTruthMatchBranch,  &vtruth);
        t->SetBranchAddress(kRetaUnfBranch,     &vretaU);
        t->SetBranchAddress(kRphiUnfBranch,     &vrphiU);
        t->SetBranchAddress(kWeta2UnfBranch,    &vweta2U);
        t->SetBranchAddress(kLepton1SFBranch,   &vl1sf);
        t->SetBranchAddress(kLepton2SFBranch,   &vl2sf);
        t->SetBranchAddress(kPhotonIsoSFBranch, &vphIsosf);
    };

    // Helper: get pT bin index (0 for eta-only mode)
    auto getPtBin = [&](double pt) -> int {
        if (!usePtBins) return 0;
        return findPtBin(pt);
    };

    // ============================================================
    //  PASS 1 -- DATA
    // ============================================================
    {
        TChain* td = makeChain(dataF);
        setBrCommon(td);
        Long64_t N = td->GetEntries();
        Long64_t nPass = 0;
        std::cout << "\n--- Pass 1: Data (" << N << " entries) ---\n";

        for (Long64_t i = 0; i < N; ++i) {
            td->GetEntry(i);
            if (i % 500000 == 0)
                std::cout << "  " << i << " / " << N << "\n";

            if (!passSelection(sel, phPt, clEta2,
                               vmll, vmllg,
                               vdrl1, vdrl2, vdrll, vcellSz,
                               visconv, vtight, visoL, visoT,
                               false, false))
                continue;

            std::vector<double> cells(vcellE->begin(), vcellE->end());
            if (!isHealthyCluster(cells)) continue;

            int nzero = 0;
            for (int k = 0; k < kClusterSize; ++k)
                if (k != kCentralCell && cells[k] == 0.0) ++nzero;
            h_nzero_d->Fill(nzero);
            ++nTotal_d;
            if (nzero > 0) ++nZeroPad_d;

            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cells[k];

            h_pt_d->Fill(phPt);
            h_eta_d->Fill(std::fabs(clEta2));

            double rc = calcReta(cells);
            double pc = calcRphi(cells);
            double wc = calcWeta2(cells, clEta2, clEtaMax2);

            h_d_comp[kReta]->Fill(rc);
            h_d_comp[kRphi]->Fill(pc);
            h_d_comp[kWeta2]->Fill(wc);
            h_d_stor[kReta]->Fill(vreta);
            h_d_stor[kRphi]->Fill(vrphi);
            h_d_stor[kWeta2]->Fill(vweta2);

            int eb = findEtaBin(std::fabs(clEta2));
            int pb = getPtBin(phPt);
            if (eb < 0 || pb < 0) continue;

            h_d_eta[kReta][eb]->Fill(rc);
            h_d_eta[kRphi][eb]->Fill(pc);
            h_d_eta[kWeta2][eb]->Fill(wc);
            h_d_stor_eta[kReta][eb]->Fill(vreta);
            h_d_stor_eta[kRphi][eb]->Fill(vrphi);
            h_d_stor_eta[kWeta2][eb]->Fill(vweta2);

            if (usePtBins) {
                h_d_etapt[kReta][eb][pb]->Fill(rc);
                h_d_etapt[kRphi][eb][pb]->Fill(pc);
                h_d_etapt[kWeta2][eb][pb]->Fill(wc);
                h_d_stor_etapt[kReta][eb][pb]->Fill(vreta);
                h_d_stor_etapt[kRphi][eb][pb]->Fill(vrphi);
                h_d_stor_etapt[kWeta2][eb][pb]->Fill(vweta2);
            }

            if (!kExcludeZeroPadFromCorr || nzero == 0) {
                for (int k = 0; k < kClusterSize; ++k) {
                    double fk = cells[k] / Etot;
                    sf_d [eb][pb][k] += fk;
                    h_frac_d[eb][pb][k]->Fill(fk);
                }
                cnt_d[eb][pb] += 1.0;
            }

            h_cnt_d->Fill(eb);
            ++nPass;
        }
        std::cout << "  Data passing: " << nPass << "\n";
        delete td;
    }

    // ============================================================
    //  PASS 1 -- MC
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
            if (i % 500000 == 0)
                std::cout << "  " << i << " / " << N << "\n";

            if (!passSelection(sel, phPt, clEta2,
                               vmll, vmllg,
                               vdrl1, vdrl2, vdrll, vcellSz,
                               visconv, vtight, visoL, visoT,
                               (bool)vtruth, true))
                continue;

            std::vector<double> cells(vcellE->begin(), vcellE->end());
            if (!isHealthyCluster(cells)) continue;

            int nzero = 0;
            for (int k = 0; k < kClusterSize; ++k)
                if (k != kCentralCell && cells[k] == 0.0) ++nzero;
            h_nzero_m->Fill(nzero, 1.0);
            ++nTotal_m;
            if (nzero > 0) ++nZeroPad_m;

            double nf = mcNormFactor(vmcch, sumWmap);
            double w = mcWeight(vmcw, vmuw, nf, vl1sf, vl2sf, vphIsosf);
            double w_nosf = mcWeight(vmcw, vmuw, nf, 1.0, 1.0, 1.0);
            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cells[k];

            h_pt_m->Fill(phPt, w);
            h_eta_m->Fill(std::fabs(clEta2), w);
            h_pt_m_nosf->Fill(phPt, w_nosf);
            h_eta_m_nosf->Fill(std::fabs(clEta2), w_nosf);

            double rc = calcReta(cells);
            double pc = calcRphi(cells);
            double wc = calcWeta2(cells, clEta2, clEtaMax2);

            h_m_comp[kReta]->Fill(rc, w);
            h_m_comp[kRphi]->Fill(pc, w);
            h_m_comp[kWeta2]->Fill(wc, w);
            h_m_fudg[kReta]->Fill(vreta, w);
            h_m_fudg[kRphi]->Fill(vrphi, w);
            h_m_fudg[kWeta2]->Fill(vweta2, w);
            h_m_unf[kReta]->Fill(vretaU, w);
            h_m_unf[kRphi]->Fill(vrphiU, w);
            h_m_unf[kWeta2]->Fill(vweta2U, w);

            int eb = findEtaBin(std::fabs(clEta2));
            int pb = getPtBin(phPt);
            if (eb < 0 || pb < 0) continue;

            h_m_eta[kReta][eb]->Fill(rc, w);
            h_m_eta[kRphi][eb]->Fill(pc, w);
            h_m_eta[kWeta2][eb]->Fill(wc, w);
            h_m_fudg_eta[kReta][eb]->Fill(vreta, w);
            h_m_fudg_eta[kRphi][eb]->Fill(vrphi, w);
            h_m_fudg_eta[kWeta2][eb]->Fill(vweta2, w);
            h_m_unf_eta[kReta][eb]->Fill(vretaU, w);
            h_m_unf_eta[kRphi][eb]->Fill(vrphiU, w);
            h_m_unf_eta[kWeta2][eb]->Fill(vweta2U, w);

            if (usePtBins) {
                h_m_etapt[kReta][eb][pb]->Fill(rc, w);
                h_m_etapt[kRphi][eb][pb]->Fill(pc, w);
                h_m_etapt[kWeta2][eb][pb]->Fill(wc, w);
                h_m_fudg_etapt[kReta][eb][pb]->Fill(vreta, w);
                h_m_fudg_etapt[kRphi][eb][pb]->Fill(vrphi, w);
                h_m_fudg_etapt[kWeta2][eb][pb]->Fill(vweta2, w);
                h_m_unf_etapt[kReta][eb][pb]->Fill(vretaU, w);
                h_m_unf_etapt[kRphi][eb][pb]->Fill(vrphiU, w);
                h_m_unf_etapt[kWeta2][eb][pb]->Fill(vweta2U, w);
            }

            if (!kExcludeZeroPadFromCorr || nzero == 0) {
                for (int k = 0; k < kClusterSize; ++k) {
                    double fk = cells[k] / Etot;
                    swf_m [eb][pb][k] += w * fk;
                    h_frac_m[eb][pb][k]->Fill(fk, w);
                }
                sw_m[eb][pb] += w;
            }

            h_cnt_m->Fill(eb, w);
            ++nPass;
        }
        std::cout << "  MC passing: " << nPass << "\n";
        delete tm;
    }

    // --- Zero-padding summary ---
    std::cout << "\n--- Zero-padding summary ---\n";
    std::cout << Form("  Data: %lld / %lld events have >=1 zero cell (%.2f%%)\n",
                      nZeroPad_d, nTotal_d,
                      nTotal_d > 0 ? 100.0 * nZeroPad_d / nTotal_d : 0.0);
    std::cout << Form("  MC:   %lld / %lld events have >=1 zero cell (%.2f%%)\n",
                      nZeroPad_m, nTotal_m,
                      nTotal_m > 0 ? 100.0 * nZeroPad_m / nTotal_m : 0.0);

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

            TH1D* hd = new TH1D(Form("h_delta%s", suffix.Data()),
                Form("M1 #Delta  %s;cell;#Delta", suffix.Data()),
                kClusterSize, 0, kClusterSize);
            TH1D* hs = new TH1D(Form("h_shift%s", suffix.Data()),
                Form("M2 shift  %s;cell;shift", suffix.Data()),
                kClusterSize, 0, kClusterSize);
            TH1D* ht = new TH1D(Form("h_stretch%s", suffix.Data()),
                Form("M2 stretch  %s;cell;stretch", suffix.Data()),
                kClusterSize, 0, kClusterSize);

            for (int k = 0; k < kClusterSize; ++k) {
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
    //  PASS 2 -- MC  (apply corrections)
    // ============================================================
    fout->cd();
    {
        TChain* tm = makeChain(mcF);
        setBrMC(tm);
        Long64_t N = tm->GetEntries();
        Long64_t nPass = 0;
        std::cout << "\n--- Pass 2: MC (" << N << " entries) ---\n";

        for (Long64_t i = 0; i < N; ++i) {
            tm->GetEntry(i);
            if (i % 500000 == 0)
                std::cout << "  " << i << " / " << N << "\n";

            if (!passSelection(sel, phPt, clEta2,
                               vmll, vmllg,
                               vdrl1, vdrl2, vdrll, vcellSz,
                               visconv, vtight, visoL, visoT,
                               (bool)vtruth, true))
                continue;

            std::vector<double> cells(vcellE->begin(), vcellE->end());
            if (!isHealthyCluster(cells)) continue;

            double nf = mcNormFactor(vmcch, sumWmap);
            double w = mcWeight(vmcw, vmuw, nf, vl1sf, vl2sf, vphIsosf);
            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cells[k];

            int eb = findEtaBin(std::fabs(clEta2));
            int pb = getPtBin(phPt);
            if (eb < 0 || pb < 0) continue;

            // ---- M1: flat shift ----
            std::vector<double> c1(kClusterSize);
            for (int k = 0; k < kClusterSize; ++k)
                c1[k] = cells[k] + corr.delta[eb][pb][k] * Etot;

            double r1 = calcReta(c1);
            double p1 = calcRphi(c1);
            double w1 = calcWeta2(c1, clEta2, clEtaMax2);

            h_m_M1[kReta]->Fill(r1, w);
            h_m_M1[kRphi]->Fill(p1, w);
            h_m_M1[kWeta2]->Fill(w1, w);
            h_m_M1_eta[kReta][eb]->Fill(r1, w);
            h_m_M1_eta[kRphi][eb]->Fill(p1, w);
            h_m_M1_eta[kWeta2][eb]->Fill(w1, w);

            if (usePtBins) {
                h_m_M1_etapt[kReta][eb][pb]->Fill(r1, w);
                h_m_M1_etapt[kRphi][eb][pb]->Fill(p1, w);
                h_m_M1_etapt[kWeta2][eb][pb]->Fill(w1, w);
            }

            // ---- M2: shift + stretch ----
            std::vector<double> c2(kClusterSize);
            for (int k = 0; k < kClusterSize; ++k)
                c2[k] = Etot * corr.shift[eb][pb][k]
                       + corr.stretch[eb][pb][k] * cells[k];

            double Ecorr = 0;
            for (int k = 0; k < kClusterSize; ++k) Ecorr += c2[k];
            if (Ecorr > 0) {
                double sc = Etot / Ecorr;
                for (int k = 0; k < kClusterSize; ++k) c2[k] *= sc;
            }

            double r2 = calcReta(c2);
            double p2 = calcRphi(c2);
            double w2 = calcWeta2(c2, clEta2, clEtaMax2);

            h_m_M2[kReta]->Fill(r2, w);
            h_m_M2[kRphi]->Fill(p2, w);
            h_m_M2[kWeta2]->Fill(w2, w);
            h_m_M2_eta[kReta][eb]->Fill(r2, w);
            h_m_M2_eta[kRphi][eb]->Fill(p2, w);
            h_m_M2_eta[kWeta2][eb]->Fill(w2, w);

            if (usePtBins) {
                h_m_M2_etapt[kReta][eb][pb]->Fill(r2, w);
                h_m_M2_etapt[kRphi][eb][pb]->Fill(p2, w);
                h_m_M2_etapt[kWeta2][eb][pb]->Fill(w2, w);
            }

            // Accumulate M2-corrected cell fractions for cell profile map
            for (int k = 0; k < kClusterSize; ++k) {
                double fk = c2[k] / Etot;
                swf_m2[eb][pb][k] += w * fk;
            }
            sw_m2[eb][pb] += w;

            ++nPass;
        }
        std::cout << "  MC pass 2 passing: " << nPass << "\n";
        delete tm;
    }

    // ============================================================
    //  CELL PROFILES -- TH2D heatmaps
    // ============================================================
    fout->mkdir("cell_profiles")->cd();
    for (int e = 0; e < kNEtaBins; ++e) {
        // Eta-integrated-over-pT profiles (sum across pT bins)
        double tot_cnt_d = 0, tot_sw_m = 0, tot_sw_m2 = 0;
        double tot_sf_d[kClusterSize] = {}, tot_swf_m[kClusterSize] = {}, tot_swf_m2[kClusterSize] = {};
        for (int p = 0; p < nPtUsed; ++p) {
            tot_cnt_d += cnt_d[e][p];
            tot_sw_m  += sw_m[e][p];
            tot_sw_m2 += sw_m2[e][p];
            for (int k = 0; k < kClusterSize; ++k) {
                tot_sf_d[k]   += sf_d[e][p][k];
                tot_swf_m[k]  += swf_m[e][p][k];
                tot_swf_m2[k] += swf_m2[e][p][k];
            }
        }

        if (tot_cnt_d < 1 && tot_sw_m < 1) continue;

        auto mk2D = [&](const char* tag, const char* ptSuf = "") {
            return new TH2D(
                Form("h_frac_%s_eta%02d%s", tag, e, ptSuf),
                Form("Cell frac %s  #eta [%.2f,%.2f)%s;#eta idx;#phi idx",
                     tag, kEtaLimits[e], kEtaLimits[e + 1], ptSuf),
                kEtaSize, 0, kEtaSize, kPhiSize, 0, kPhiSize);
        };

        // Eta-only profiles (always created)
        TH2D* h2_mean_d  = mk2D("mean_data");
        TH2D* h2_mean_m  = mk2D("mean_mc");
        TH2D* h2_mean_m2 = mk2D("mean_mc_M2");
        TH2D* h2_delta   = mk2D("delta");
        TH2D* h2_rms_d   = mk2D("rms_data");
        TH2D* h2_rms_m   = mk2D("rms_mc");

        for (int k = 0; k < kClusterSize; ++k) {
            int ie = k / kPhiSize;
            int ip = k % kPhiSize;

            double md  = (tot_cnt_d > 0) ? tot_sf_d[k]   / tot_cnt_d : 0;
            double mm  = (tot_sw_m  > 0) ? tot_swf_m[k]  / tot_sw_m  : 0;
            double mm2 = (tot_sw_m2 > 0) ? tot_swf_m2[k] / tot_sw_m2 : 0;

            // RMS: sum over pT bins for eta-integrated RMS
            double rms_d_sum2 = 0, rms_m_sum2 = 0;
            double rms_d_w = 0, rms_m_w = 0;
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

        // Per-pT-bin profiles (only in eta_pt mode)
        if (usePtBins) {
            for (int p = 0; p < kNPtBins; ++p) {
                if (cnt_d[e][p] < 1 && sw_m[e][p] < 1) continue;
                TString ptSuf = Form("_pt%02d", p);

                TH2D* h2p_mean_d  = mk2D("mean_data", ptSuf.Data());
                TH2D* h2p_mean_m  = mk2D("mean_mc",   ptSuf.Data());
                TH2D* h2p_mean_m2 = mk2D("mean_mc_M2", ptSuf.Data());
                TH2D* h2p_delta   = mk2D("delta",     ptSuf.Data());
                TH2D* h2p_rms_d   = mk2D("rms_data",  ptSuf.Data());
                TH2D* h2p_rms_m   = mk2D("rms_mc",    ptSuf.Data());

                for (int k = 0; k < kClusterSize; ++k) {
                    int ie = k / kPhiSize;
                    int ip = k % kPhiSize;

                    double md  = (cnt_d[e][p]  > 0) ? sf_d[e][p][k]   / cnt_d[e][p]  : 0;
                    double mm  = (sw_m[e][p]   > 0) ? swf_m[e][p][k]  / sw_m[e][p]   : 0;
                    double mm2 = (sw_m2[e][p]  > 0) ? swf_m2[e][p][k] / sw_m2[e][p]  : 0;
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
    const char* ch  = (argc > 1) ? argv[1] : "eegamma";
    const char* sc  = (argc > 2) ? argv[2] : "unconverted";
    const char* bd  = (argc > 3) ? argv[3] : "../../output/Layer_2/eta_loose";
    const char* bn  = (argc > 4) ? argv[4] : "eta";
    const char* iso = (argc > 5) ? argv[5] : "loose";
    fill_histograms(ch, sc, bd, bn, iso);
    return 0;
}
#endif
