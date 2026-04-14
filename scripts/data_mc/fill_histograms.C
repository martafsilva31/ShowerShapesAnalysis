///////////////////////////////////////////////////////////////////////////////
// fill_histograms.C
//
// Two-pass pipeline for cell-energy reweighting corrections.
//
//   Pass 1: Loop data + MC → accumulate per-cell fraction statistics,
//           fill uncorrected shower-shape histograms.
//   Pass 2: Loop MC only   → apply M1 (flat shift) and M2 (shift+stretch)
//           corrections, fill corrected shower-shape histograms.
//
// Usage:
//   root -l -b -q 'fill_histograms.C("eegamma", "baseline")'
//   root -l -b -q 'fill_histograms.C("mumugamma", "tight_id")'
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
const bool kExcludeZeroPadFromCorr = true;  // skip zero-padded events in correction derivation
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
// Correction storage
// ============================================================
struct Corrections {
    double delta  [kNEtaBins][kClusterSize] = {};   // M1
    double shift  [kNEtaBins][kClusterSize] = {};   // M2
    double stretch[kNEtaBins][kClusterSize] = {};   // M2
};

// ============================================================
// Main
// ============================================================
void fill_histograms(const char* channel  = "eegamma",
                     const char* scenario = "baseline",
                     const char* baseDir  =
                         "../../output/cell_energy_reweighting_Francisco_method/data24") {

    TH1::SetDefaultSumw2(true);

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

    Selection sel = getSelection(scenario);

    // --------------------------------------------------------
    // Output
    // --------------------------------------------------------
    TString outPath = Form("%s/%s/%s", baseDir, channel, scenario);
    gSystem->mkdir(outPath, true);
    TString outFile = outPath + "/histograms.root";
    TFile* fout = TFile::Open(outFile, "RECREATE");

    std::cout << "=== fill_histograms ===\n"
              << "  channel  = " << channel  << "\n"
              << "  scenario = " << scenario << "\n"
              << "  data     = " << dataF    << "\n"
              << "  mc       = " << mcF      << "\n"
              << "  output   = " << outFile  << "\n";

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
    TH1D *h_d_eta[kNSS][kNEtaBins];   // cell-computed data
    TH1D *h_m_eta[kNSS][kNEtaBins];   // cell-computed MC
    TH1D *h_m_M1_eta[kNSS][kNEtaBins];
    TH1D *h_m_M2_eta[kNSS][kNEtaBins];
    TH1D *h_d_stor_eta[kNSS][kNEtaBins];  // branch-stored data (per-eta)
    TH1D *h_m_fudg_eta[kNSS][kNEtaBins];  // fudged MC (per-eta)
    for (int s = 0; s < kNSS; ++s)
        for (int e = 0; e < kNEtaBins; ++e) {
            TString suf = Form("_eta%02d", e);
            h_d_eta[s][e]      = makeSS(s, "data",         suf.Data());
            h_m_eta[s][e]      = makeSS(s, "mc",           suf.Data());
            h_m_M1_eta[s][e]   = makeSS(s, "mc_M1",        suf.Data());
            h_m_M2_eta[s][e]   = makeSS(s, "mc_M2",        suf.Data());
            h_d_stor_eta[s][e] = makeSS(s, "data_stored",  suf.Data());
            h_m_fudg_eta[s][e] = makeSS(s, "mc_fudged",    suf.Data());
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
        "Photon p_{T} (data);p_{T} [GeV];Events / 2 GeV",
        50, 0, 100);
    TH1D* h_pt_m = new TH1D("h_pt_mc",
        "Photon p_{T} (MC);p_{T} [GeV];Events / 2 GeV",
        50, 0, 100);
    TH1D* h_eta_d = new TH1D("h_eta_data",
        "Photon |#eta| (data);|#eta|;Events",
        48, 0, 2.4);
    TH1D* h_eta_m = new TH1D("h_eta_mc",
        "Photon |#eta| (MC);|#eta|;Events",
        48, 0, 2.4);
    TH1D* h_pt_m_nosf = new TH1D("h_pt_mc_nosf",
        "Photon p_{T} (MC, no SFs);p_{T} [GeV];Events / 2 GeV",
        50, 0, 100);
    TH1D* h_eta_m_nosf = new TH1D("h_eta_mc_nosf",
        "Photon |#eta| (MC, no SFs);|#eta|;Events",
        48, 0, 2.4);

    // --- Zero-padding diagnostics ---
    Long64_t nZeroPad_d = 0, nTotal_d = 0;
    Long64_t nZeroPad_m = 0, nTotal_m = 0;
    TH1D* h_nzero_d = new TH1D("h_nzero_data",
        "Cells with E=0 per event (data);N_{zero};Events",
        78, 0, 78);
    TH1D* h_nzero_m = new TH1D("h_nzero_mc",
        "Cells with E=0 per event (MC);N_{zero};Events",
        78, 0, 78);

    // ============================================================
    //  ACCUMULATION ARRAYS  (Pass 1)
    // ============================================================
    double sf_d [kNEtaBins][kClusterSize] = {};   // sum f_k      (data)
    double sf2_d[kNEtaBins][kClusterSize] = {};   // sum f_k^2    (data)
    double cnt_d[kNEtaBins]               = {};   // count        (data)

    double swf_m [kNEtaBins][kClusterSize] = {};  // sum w*f_k    (MC)
    double swf2_m[kNEtaBins][kClusterSize] = {};  // sum w*f_k^2  (MC)
    double sw_m  [kNEtaBins]               = {};  // sum w        (MC)

    double swf_m2[kNEtaBins][kClusterSize] = {};  // sum w*f_k    (MC M2-corrected)
    double sw_m2 [kNEtaBins]               = {};  // sum w        (MC M2-corrected)

    // ============================================================
    //  BRANCH VARIABLES  (reused across chains)
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
        t->SetBranchAddress(kPhotonIsoSFBranch,  &vphIsosf);
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

            // Zero-padding diagnostic
            int nzero = 0;
            for (int k = 0; k < kClusterSize; ++k)
                if (k != kCentralCell && cells[k] == 0.0) ++nzero;
            h_nzero_d->Fill(nzero);
            ++nTotal_d;
            if (nzero > 0) ++nZeroPad_d;

            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cells[k];

            // pT and eta distributions
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
            if (eb >= 0) {
                h_d_eta[kReta][eb]->Fill(rc);
                h_d_eta[kRphi][eb]->Fill(pc);
                h_d_eta[kWeta2][eb]->Fill(wc);
                h_d_stor_eta[kReta][eb]->Fill(vreta);
                h_d_stor_eta[kRphi][eb]->Fill(vrphi);
                h_d_stor_eta[kWeta2][eb]->Fill(vweta2);

                if (!kExcludeZeroPadFromCorr || nzero == 0) {
                    for (int k = 0; k < kClusterSize; ++k) {
                        double fk = cells[k] / Etot;
                        sf_d [eb][k] += fk;
                        sf2_d[eb][k] += fk * fk;
                    }
                    cnt_d[eb] += 1.0;
                }
            }
            ++nPass;
        }
        std::cout << "  Data passing: " << nPass << "\n";
        delete td;
    }

    // ============================================================
    //  PASS 1 — MC
    // ============================================================
    // Sum-of-weights from h_sumW histogram (bin 2) in the MC ntuple files
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

            // Zero-padding diagnostic
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

            // pT and eta distributions
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
            if (eb >= 0) {
                h_m_eta[kReta][eb]->Fill(rc, w);
                h_m_eta[kRphi][eb]->Fill(pc, w);
                h_m_eta[kWeta2][eb]->Fill(wc, w);
                h_m_fudg_eta[kReta][eb]->Fill(vreta, w);
                h_m_fudg_eta[kRphi][eb]->Fill(vrphi, w);
                h_m_fudg_eta[kWeta2][eb]->Fill(vweta2, w);

                if (!kExcludeZeroPadFromCorr || nzero == 0) {
                    for (int k = 0; k < kClusterSize; ++k) {
                        double fk = cells[k] / Etot;
                        swf_m [eb][k] += w * fk;
                        swf2_m[eb][k] += w * fk * fk;
                    }
                    sw_m[eb] += w;
                }
            }
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
        std::cout << Form("  eta bin %2d  [%.2f, %.2f)  data=%6.0f  MC_wt=%.1f",
                          e, kEtaLimits[e], kEtaLimits[e + 1],
                          cnt_d[e], sw_m[e]) << "\n";

        h_cnt_d->SetBinContent(e + 1, cnt_d[e]);
        h_cnt_m->SetBinContent(e + 1, sw_m[e]);

        if (cnt_d[e] < 10 || sw_m[e] < 10) {
            std::cout << "    -> too few events, skipping corrections\n";
            continue;
        }

        TH1D* hd = new TH1D(Form("h_delta_eta%02d", e),
            Form("M1 #Delta  #eta bin %d;cell;#Delta", e),
            kClusterSize, 0, kClusterSize);
        TH1D* hs = new TH1D(Form("h_shift_eta%02d", e),
            Form("M2 shift  #eta bin %d;cell;shift", e),
            kClusterSize, 0, kClusterSize);
        TH1D* ht = new TH1D(Form("h_stretch_eta%02d", e),
            Form("M2 stretch  #eta bin %d;cell;stretch", e),
            kClusterSize, 0, kClusterSize);

        for (int k = 0; k < kClusterSize; ++k) {
            double mu_d  = sf_d[e][k]  / cnt_d[e];
            double mu_m  = swf_m[e][k] / sw_m[e];
            double var_d = sf2_d[e][k] / cnt_d[e] - mu_d * mu_d;
            double var_m = swf2_m[e][k] / sw_m[e] - mu_m * mu_m;
            double sig_d = (var_d > 0) ? std::sqrt(var_d) : 0;
            double sig_m = (var_m > 0) ? std::sqrt(var_m) : 0;

            // M1: flat shift
            corr.delta[e][k] = mu_d - mu_m;

            // M2: shift + stretch
            double str = (sig_m > 1e-12) ? sig_d / sig_m : 1.0;
            corr.stretch[e][k] = str;
            corr.shift[e][k]   = mu_d - str * mu_m;

            hd->SetBinContent(k + 1, corr.delta[e][k]);
            hs->SetBinContent(k + 1, corr.shift[e][k]);
            ht->SetBinContent(k + 1, corr.stretch[e][k]);
        }
    }

    // ============================================================
    //  PASS 2 — MC  (apply corrections)
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
            if (eb < 0) continue;

            // ---- M1: flat shift ----
            std::vector<double> c1(kClusterSize);
            for (int k = 0; k < kClusterSize; ++k)
                c1[k] = cells[k] + corr.delta[eb][k] * Etot;

            double r1 = calcReta(c1);
            double p1 = calcRphi(c1);
            double w1 = calcWeta2(c1, clEta2, clEtaMax2);

            h_m_M1[kReta]->Fill(r1, w);
            h_m_M1[kRphi]->Fill(p1, w);
            h_m_M1[kWeta2]->Fill(w1, w);
            h_m_M1_eta[kReta][eb]->Fill(r1, w);
            h_m_M1_eta[kRphi][eb]->Fill(p1, w);
            h_m_M1_eta[kWeta2][eb]->Fill(w1, w);

            // ---- M2: shift + stretch ----
            std::vector<double> c2(kClusterSize);
            for (int k = 0; k < kClusterSize; ++k)
                c2[k] = Etot * corr.shift[eb][k]
                       + corr.stretch[eb][k] * cells[k];

            // Rescale to conserve total energy
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

            // Accumulate M2-corrected cell fractions for cell profile map
            for (int k = 0; k < kClusterSize; ++k) {
                double fk = c2[k] / Etot;   // Etot preserved by rescaling
                swf_m2[eb][k] += w * fk;
            }
            sw_m2[eb] += w;

            ++nPass;
        }
        std::cout << "  MC pass 2 passing: " << nPass << "\n";
        delete tm;
    }

    // ============================================================
    //  CELL PROFILES — TH2D heatmaps
    // ============================================================
    fout->mkdir("cell_profiles")->cd();
    for (int e = 0; e < kNEtaBins; ++e) {
        if (cnt_d[e] < 1 && sw_m[e] < 1) continue;

        auto mk2D = [&](const char* tag) {
            return new TH2D(
                Form("h_frac_%s_eta%02d", tag, e),
                Form("Cell frac %s  #eta [%.2f,%.2f);#eta idx;#phi idx",
                     tag, kEtaLimits[e], kEtaLimits[e + 1]),
                kEtaSize, 0, kEtaSize, kPhiSize, 0, kPhiSize);
        };

        TH2D* h2_mean_d  = mk2D("mean_data");
        TH2D* h2_mean_m  = mk2D("mean_mc");
        TH2D* h2_mean_m2 = mk2D("mean_mc_M2");
        TH2D* h2_delta   = mk2D("delta");
        TH2D* h2_rms_d   = mk2D("rms_data");
        TH2D* h2_rms_m   = mk2D("rms_mc");

        for (int k = 0; k < kClusterSize; ++k) {
            int ie = k / kPhiSize;
            int ip = k % kPhiSize;

            double md  = (cnt_d[e]  > 0) ? sf_d[e][k]  / cnt_d[e]  : 0;
            double mm  = (sw_m[e]   > 0) ? swf_m[e][k] / sw_m[e]   : 0;
            double mm2 = (sw_m2[e]  > 0) ? swf_m2[e][k]/ sw_m2[e]  : 0;
            double vd  = (cnt_d[e]  > 0) ? sf2_d[e][k] / cnt_d[e]  - md  * md  : 0;
            double vm  = (sw_m[e]   > 0) ? swf2_m[e][k]/ sw_m[e]   - mm  * mm  : 0;

            h2_mean_d ->SetBinContent(ie + 1, ip + 1, md);
            h2_mean_m ->SetBinContent(ie + 1, ip + 1, mm);
            h2_mean_m2->SetBinContent(ie + 1, ip + 1, mm2);
            h2_delta  ->SetBinContent(ie + 1, ip + 1, md - mm);
            h2_rms_d  ->SetBinContent(ie + 1, ip + 1, (vd > 0) ? std::sqrt(vd) : 0);
            h2_rms_m  ->SetBinContent(ie + 1, ip + 1, (vm > 0) ? std::sqrt(vm) : 0);
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
    const char* ch = (argc > 1) ? argv[1] : "eegamma";
    const char* sc = (argc > 2) ? argv[2] : "baseline";
    const char* bd = (argc > 3) ? argv[3]
                     : "../../output/cell_energy_reweighting_Francisco_method/data24";
    fill_histograms(ch, sc, bd);
    return 0;
}
#endif
