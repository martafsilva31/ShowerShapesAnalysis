///////////////////////////////////////////////////////////////////////////////
// validate_data_mc.C
//
// Applies cell-energy corrections and produces comparison histograms.
// Reads corrections from derive_corrections.C output.
//
// Two SEPARATE histogram sets for fair comparison:
//
//   SET A — Branch-based shower shapes (ATLAS-reconstructed truth):
//     data_br:  Data branch values (photon.reta, .rphi, .weta2)
//     mc_unf:   MC unfudged branch values (photon.unfudged_*)
//     mc_fud:   MC fudged branch values (photon.reta, .rphi, .weta2)
//
//   SET B — Cell-computed shower shapes (position-based):
//     data_cell: Data from calcRetaPos etc
//     mc_cell:   MC from calcRetaPos etc
//     m1:        MC with Method 1 (flat shift from CalCoef.C)
//     m2:        MC with Method 2 (CDF quantile mapping)
//
// Using separate sets avoids the ~0.006-0.008 offset between branch and
// cell-computed values that invalidated previous comparisons.
//
// M1 application (document §5.1, CalCoef.C + NTUP.C):
//   E_new_k = E_old_k + Delta_k * E_total
//   Delta_k = <e_data>_k - <e_MC>_k (population mean fraction shift)
//   sum(Delta_k) = 0 by construction => total energy preserved, NO renorm
//
// M2 application (DeltaR-matched TProfile, Francisco's tree.C):
//   alpha = m2_prof->Interpolate(MC_frac_k)
//   E_new_k = E_old_k + alpha * E_total
//   Then rescale: E_new_k *= sum(E_old)/sum(E_new)  [preserve total E]
//   The TProfile stores alpha(f_MC) = <f_data - f_MC> from matched pairs.
//
// Usage:
//   root -l -b -q 'validate_data_mc.C("data.root","mc.root","corr.root","histos.root")'
//   root -l -b -q 'validate_data_mc.C("data.root","mc.root","corr.root","histos.root",500000)'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace datamc;


// ======================================================================
// Chi2 between two area-normalised histograms
// ======================================================================
double computeChi2(TH1D* h1, TH1D* h2, int& ndf) {
    double chi2 = 0;
    ndf = 0;
    for (int b = 1; b <= h1->GetNbinsX(); ++b) {
        double v1 = h1->GetBinContent(b);
        double v2 = h2->GetBinContent(b);
        if (v2 > 0) {
            chi2 += (v1 - v2) * (v1 - v2) / v2;
            ndf++;
        }
    }
    return chi2;
}


// ======================================================================
// Main
// ======================================================================
int validate_data_mc(const char* dataFile,
                     const char* mcFile,
                     const char* corrFile,
                     const char* histoFile,
                     Long64_t maxEvents = -1) {

    // ==================================================================
    // Load M1 corrections: cellsDelta TH2D per eta bin
    // Delta_k = <e_data>_k - <e_MC>_k (fraction space)
    // ==================================================================
    TFile* fcorr = TFile::Open(corrFile, "READ");
    if (!fcorr || fcorr->IsZombie()) {
        std::cerr << "ERROR: Cannot open corrections: " << corrFile << std::endl;
        return 1;
    }

    // M1: flat shift per cell per eta bin
    double m1_delta[kClusterSize][kNEtaBins] = {};

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* h_delta = (TH2D*)fcorr->Get("cellsDelta_" + suffix);
        if (!h_delta) {
            std::cerr << "WARNING: Missing cellsDelta_" << suffix
                      << " — M1 will be zero for this eta bin" << std::endl;
            continue;
        }

        for (int phi = 1; phi <= kPhiSize; ++phi) {
            for (int eta = 1; eta <= kEtaSize; ++eta) {
                int k = (phi - 1) + kPhiSize * (eta - 1);
                m1_delta[k][n] = h_delta->GetBinContent(phi, eta);
            }
        }
    }

    // M2: DeltaR-matched TProfile per cell per eta bin
    TProfile* m2_prof[kClusterSize][kNEtaBins] = {};
    bool hasM2 = true;

    for (int k = 0; k < kClusterSize; ++k) {
        for (int n = 0; n < kNEtaBins; ++n) {
            TString name = Form("m2_alpha_Cell_%d_Eta_%1.2f_%1.2f",
                                k + 1, kEtaLimits[n], kEtaLimits[n + 1]);
            m2_prof[k][n] = (TProfile*)fcorr->Get(name);
            if (!m2_prof[k][n] && k == 0 && n == 0) {
                std::cerr << "WARNING: M2 TProfiles not found — skipping M2"
                          << std::endl;
                hasM2 = false;
            }
        }
    }

    std::cout << "Corrections loaded from: " << corrFile << std::endl;
    std::cout << "  M1: cellsDelta loaded" << std::endl;
    std::cout << "  M2: TProfile " << (hasM2 ? "loaded" : "NOT found") << std::endl;

    // M3: shift+stretch parameters per cell per eta bin
    double m3_meanData[kClusterSize][kNEtaBins] = {};
    double m3_meanMC[kClusterSize][kNEtaBins] = {};
    double m3_sigRatio[kClusterSize][kNEtaBins] = {};
    bool hasM3 = true;

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);
        TH2D* h_md = (TH2D*)fcorr->Get("cellsMeanData_" + suffix);
        TH2D* h_mm = (TH2D*)fcorr->Get("cellsMeanMC_" + suffix);
        TH2D* h_sr = (TH2D*)fcorr->Get("cellsSigmaRatio_" + suffix);

        if (!h_md || !h_mm || !h_sr) {
            if (n == 0) {
                std::cerr << "WARNING: M3 histograms not found — skipping M3"
                          << std::endl;
                hasM3 = false;
            }
            continue;
        }

        for (int phi = 1; phi <= kPhiSize; ++phi) {
            for (int eta = 1; eta <= kEtaSize; ++eta) {
                int k = (phi - 1) + kPhiSize * (eta - 1);
                m3_meanData[k][n] = h_md->GetBinContent(phi, eta);
                m3_meanMC[k][n]   = h_mm->GetBinContent(phi, eta);
                m3_sigRatio[k][n] = h_sr->GetBinContent(phi, eta);
                // Defensive cap on sigma ratio
                if (m3_sigRatio[k][n] < kSigmaRatioMin)
                    m3_sigRatio[k][n] = kSigmaRatioMin;
                if (m3_sigRatio[k][n] > kSigmaRatioMax)
                    m3_sigRatio[k][n] = kSigmaRatioMax;
            }
        }
    }
    std::cout << "  M3: shift+stretch " << (hasM3 ? "loaded" : "NOT found")
              << std::endl;


    // ==================================================================
    // Create histograms
    //
    // SET A — Branch-based: 3 samples × 3 vars × 14 eta bins
    // SET B — Cell-computed: 4 samples × 3 vars × 14 eta bins
    // ==================================================================

    // Set A: branch-based
    TH1D* h_reta_data_br[kNEtaBins], *h_reta_mc_unf[kNEtaBins],
         *h_reta_mc_fud[kNEtaBins];
    TH1D* h_rphi_data_br[kNEtaBins], *h_rphi_mc_unf[kNEtaBins],
         *h_rphi_mc_fud[kNEtaBins];
    TH1D* h_weta2_data_br[kNEtaBins], *h_weta2_mc_unf[kNEtaBins],
         *h_weta2_mc_fud[kNEtaBins];

    // Set B: cell-computed
    TH1D* h_reta_data_cell[kNEtaBins], *h_reta_mc_cell[kNEtaBins],
         *h_reta_m1[kNEtaBins], *h_reta_m2[kNEtaBins],
         *h_reta_m3[kNEtaBins];
    TH1D* h_rphi_data_cell[kNEtaBins], *h_rphi_mc_cell[kNEtaBins],
         *h_rphi_m1[kNEtaBins], *h_rphi_m2[kNEtaBins],
         *h_rphi_m3[kNEtaBins];
    TH1D* h_weta2_data_cell[kNEtaBins], *h_weta2_mc_cell[kNEtaBins],
         *h_weta2_m1[kNEtaBins], *h_weta2_m2[kNEtaBins],
         *h_weta2_m3[kNEtaBins];

    for (int n = 0; n < kNEtaBins; ++n) {
        TString s = Form("_%d", n);

        // Set A
        h_reta_data_br[n] = new TH1D("reta_data_br" + s, "", 25, 0.90, 1.0);
        h_reta_mc_unf[n]  = new TH1D("reta_mc_unf" + s,  "", 25, 0.90, 1.0);
        h_reta_mc_fud[n]  = new TH1D("reta_mc_fud" + s,  "", 25, 0.90, 1.0);
        h_rphi_data_br[n] = new TH1D("rphi_data_br" + s, "", 25, 0.85, 1.0);
        h_rphi_mc_unf[n]  = new TH1D("rphi_mc_unf" + s,  "", 25, 0.85, 1.0);
        h_rphi_mc_fud[n]  = new TH1D("rphi_mc_fud" + s,  "", 25, 0.85, 1.0);
        h_weta2_data_br[n] = new TH1D("weta2_data_br" + s, "", 25, 0.008, 0.016);
        h_weta2_mc_unf[n]  = new TH1D("weta2_mc_unf" + s,  "", 25, 0.008, 0.016);
        h_weta2_mc_fud[n]  = new TH1D("weta2_mc_fud" + s,  "", 25, 0.008, 0.016);

        // Set B
        h_reta_data_cell[n] = new TH1D("reta_data_cell" + s, "", 25, 0.90, 1.0);
        h_reta_mc_cell[n]   = new TH1D("reta_mc_cell" + s,   "", 25, 0.90, 1.0);
        h_reta_m1[n]        = new TH1D("reta_m1" + s,        "", 25, 0.90, 1.0);
        h_reta_m2[n]        = new TH1D("reta_m2" + s,        "", 25, 0.90, 1.0);
        h_reta_m3[n]        = new TH1D("reta_m3" + s,        "", 25, 0.90, 1.0);
        h_rphi_data_cell[n] = new TH1D("rphi_data_cell" + s, "", 25, 0.85, 1.0);
        h_rphi_mc_cell[n]   = new TH1D("rphi_mc_cell" + s,   "", 25, 0.85, 1.0);
        h_rphi_m1[n]        = new TH1D("rphi_m1" + s,        "", 25, 0.85, 1.0);
        h_rphi_m2[n]        = new TH1D("rphi_m2" + s,        "", 25, 0.85, 1.0);
        h_rphi_m3[n]        = new TH1D("rphi_m3" + s,        "", 25, 0.85, 1.0);
        h_weta2_data_cell[n] = new TH1D("weta2_data_cell" + s, "", 25, 0.008, 0.016);
        h_weta2_mc_cell[n]   = new TH1D("weta2_mc_cell" + s,   "", 25, 0.008, 0.016);
        h_weta2_m1[n]        = new TH1D("weta2_m1" + s,        "", 25, 0.008, 0.016);
        h_weta2_m2[n]        = new TH1D("weta2_m2" + s,        "", 25, 0.008, 0.016);
        h_weta2_m3[n]        = new TH1D("weta2_m3" + s,        "", 25, 0.008, 0.016);

        // Sumw2 for proper error propagation with MC weights
        h_reta_data_br[n]->Sumw2(); h_reta_mc_unf[n]->Sumw2();
        h_reta_mc_fud[n]->Sumw2();
        h_rphi_data_br[n]->Sumw2(); h_rphi_mc_unf[n]->Sumw2();
        h_rphi_mc_fud[n]->Sumw2();
        h_weta2_data_br[n]->Sumw2(); h_weta2_mc_unf[n]->Sumw2();
        h_weta2_mc_fud[n]->Sumw2();

        h_reta_data_cell[n]->Sumw2(); h_reta_mc_cell[n]->Sumw2();
        h_reta_m1[n]->Sumw2(); h_reta_m2[n]->Sumw2(); h_reta_m3[n]->Sumw2();
        h_rphi_data_cell[n]->Sumw2(); h_rphi_mc_cell[n]->Sumw2();
        h_rphi_m1[n]->Sumw2(); h_rphi_m2[n]->Sumw2(); h_rphi_m3[n]->Sumw2();
        h_weta2_data_cell[n]->Sumw2(); h_weta2_mc_cell[n]->Sumw2();
        h_weta2_m1[n]->Sumw2(); h_weta2_m2[n]->Sumw2(); h_weta2_m3[n]->Sumw2();
    }

    // ==================================================================
    // Helper: enable a branch safely
    // ==================================================================
    auto enableBranch = [](TTree* t, const char* name, void* addr) {
        if (t->GetBranch(name)) {
            t->SetBranchStatus(name, 1);
            t->SetBranchAddress(name, addr);
        } else {
            std::cerr << "WARNING: Branch '" << name << "' not found" << std::endl;
        }
    };


    // ==================================================================
    // Cell fraction accumulators for profile maps (written as TH2D)
    // ==================================================================
    double prof_data[kClusterSize][kNEtaBins] = {};
    double prof_data_w[kNEtaBins] = {};
    double prof_mc[kClusterSize][kNEtaBins] = {};
    double prof_m1[kClusterSize][kNEtaBins] = {};
    double prof_m2[kClusterSize][kNEtaBins] = {};
    double prof_m3[kClusterSize][kNEtaBins] = {};
    double prof_mc_w[kNEtaBins] = {};

    // ==================================================================
    // Loop 1: Data
    //   - Set A: branch shower shapes (data_br)
    //   - Set B: cell-computed shower shapes (data_cell)
    // ==================================================================
    std::cout << "\n=== Processing Data ===" << std::endl;
    {
        TChain* t = makeChain(dataFile);
        if (!t || t->GetEntries() == 0) {
            std::cerr << "ERROR: Cannot open data " << dataFile << std::endl;
            delete t; return 1;
        }

        // Cell energies + positions
        std::vector<double>* cellE = nullptr;
        std::vector<double>* cellEta = nullptr;
        std::vector<double>* cellPhi = nullptr;
        Int_t cellSize = 0;
        Float_t eta2 = 0, photonPt = 0, etamax2 = 0;
        Double_t mll = 0, mllg = 0;
        Double_t dR_l1_ph = 0, dR_l2_ph = 0, dR_ll = 0;
        Int_t tightID = 0;
        Bool_t isConv = false;

        // Branch shower shapes (for Set A)
        Float_t br_reta = 0, br_rphi = 0, br_weta2 = 0;

        t->SetBranchStatus("*", 0);
        enableBranch(t, kCellBranch, &cellE);
        enableBranch(t, kCellEtaBranch, &cellEta);
        enableBranch(t, kCellPhiBranch, &cellPhi);
        enableBranch(t, kCellSizeBranch, &cellSize);
        enableBranch(t, kEtaBranch, &eta2);
        enableBranch(t, kEtamax2Branch, &etamax2);
        enableBranch(t, kPhotonPtBranch, &photonPt);
        enableBranch(t, kMllBranch, &mll);
        enableBranch(t, kMllgBranch, &mllg);
        enableBranch(t, kDRLepton1Branch, &dR_l1_ph);
        enableBranch(t, kDRLepton2Branch, &dR_l2_ph);
        enableBranch(t, kDRllBranch, &dR_ll);
        enableBranch(t, kTightIDBranch, &tightID);
        enableBranch(t, kIsConvBranch, &isConv);
        // Data uses fudged branch names (photon.reta etc), same as mc fudged
        enableBranch(t, kRetaFudBranch, &br_reta);
        enableBranch(t, kRphiFudBranch, &br_rphi);
        enableBranch(t, kWeta2FudBranch, &br_weta2);

        Long64_t nEntries = t->GetEntries();
        if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
        std::cout << "Data: " << dataFile << " (" << nEntries
                  << " entries)" << std::endl;

        Long64_t nUsed = 0;
        for (Long64_t i = 0; i < nEntries; ++i) {
            t->GetEntry(i);
            if (i > 0 && i % 2000000 == 0)
                std::cout << "  " << i << " / " << nEntries << std::endl;

            if (!passSelection(photonPt, eta2, mll, mllg,
                               dR_l1_ph, dR_l2_ph, dR_ll,
                               cellSize, true, false,
                               tightID, isConv))
                continue;
            if (!cellE || (int)cellE->size() != kClusterSize) continue;
            if (!isHealthyCluster(*cellE)) continue;

            int etaBin = findEtaBin(std::fabs(eta2));
            if (etaBin < 0) continue;

            // Set A: branch values
            h_reta_data_br[etaBin]->Fill(br_reta);
            h_rphi_data_br[etaBin]->Fill(br_rphi);
            h_weta2_data_br[etaBin]->Fill(br_weta2);

            // Set B: cell-computed values (position-based)
            if (cellEta && cellPhi) {
                h_reta_data_cell[etaBin]->Fill(
                    calcRetaPos(*cellE, *cellEta, *cellPhi));
                h_rphi_data_cell[etaBin]->Fill(
                    calcRphiPos(*cellE, *cellEta, *cellPhi));
                h_weta2_data_cell[etaBin]->Fill(
                    calcWeta2Pos(*cellE, *cellEta, *cellPhi, eta2, etamax2));
            }

            // Cell fraction accumulator for profile maps
            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cellE->at(k);
            if (Etot > 0) {
                for (int k = 0; k < kClusterSize; ++k)
                    prof_data[k][etaBin] += cellE->at(k) / Etot;
                prof_data_w[etaBin] += 1.0;
            }

            nUsed++;
        }
        std::cout << "  Selected: " << nUsed << " events" << std::endl;
        delete t;
    }


    // ==================================================================
    // Loop 2: MC
    //   - Set A: branch MC unfudged (mc_unf) + MC fudged (mc_fud)
    //   - Set B: cell-computed MC (mc_cell), M1, M2
    // ==================================================================
    std::cout << "\n=== Processing MC ===" << std::endl;
    {
        TChain* t = makeChain(mcFile);
        if (!t || t->GetEntries() == 0) {
            std::cerr << "ERROR: Cannot open MC " << mcFile << std::endl;
            delete t; return 1;
        }

        // Cell energies + positions
        std::vector<double>* cellE = nullptr;
        std::vector<double>* cellEta = nullptr;
        std::vector<double>* cellPhi = nullptr;
        Int_t cellSize = 0;
        Float_t eta2 = 0, photonPt = 0, etamax2 = 0;
        Double_t mll = 0, mllg = 0;
        Double_t dR_l1_ph = 0, dR_l2_ph = 0, dR_ll = 0;
        Bool_t truthMatch = false;
        Int_t tightID = 0;
        Bool_t isConv = false;

        // MC weights
        Float_t mcwgt = 1, muwgt = 1;

        // Branch shower shapes
        Float_t reta_fud = 0, rphi_fud = 0, weta2_fud = 0;
        Float_t reta_unf = 0, rphi_unf = 0, weta2_unf = 0;

        t->SetBranchStatus("*", 0);
        enableBranch(t, kCellBranch, &cellE);
        enableBranch(t, kCellEtaBranch, &cellEta);
        enableBranch(t, kCellPhiBranch, &cellPhi);
        enableBranch(t, kCellSizeBranch, &cellSize);
        enableBranch(t, kEtaBranch, &eta2);
        enableBranch(t, kEtamax2Branch, &etamax2);
        enableBranch(t, kPhotonPtBranch, &photonPt);
        enableBranch(t, kMllBranch, &mll);
        enableBranch(t, kMllgBranch, &mllg);
        enableBranch(t, kDRLepton1Branch, &dR_l1_ph);
        enableBranch(t, kDRLepton2Branch, &dR_l2_ph);
        enableBranch(t, kDRllBranch, &dR_ll);
        enableBranch(t, kTruthMatchBranch, &truthMatch);
        enableBranch(t, kTightIDBranch, &tightID);
        enableBranch(t, kIsConvBranch, &isConv);
        enableBranch(t, kMCWeightBranch, &mcwgt);
        enableBranch(t, kPUWeightBranch, &muwgt);
        enableBranch(t, kRetaFudBranch, &reta_fud);
        enableBranch(t, kRphiFudBranch, &rphi_fud);
        enableBranch(t, kWeta2FudBranch, &weta2_fud);
        enableBranch(t, kRetaUnfBranch, &reta_unf);
        enableBranch(t, kRphiUnfBranch, &rphi_unf);
        enableBranch(t, kWeta2UnfBranch, &weta2_unf);

        Long64_t nEntries = t->GetEntries();
        if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
        std::cout << "MC: " << mcFile << " (" << nEntries
                  << " entries)" << std::endl;

        Long64_t nUsed = 0, nM2Applied = 0;
        // Renormalization diagnostic: track scale factors per eta bin
        double m2_scaleSum[kNEtaBins] = {0};
        double m3_scaleSum[kNEtaBins] = {0};
        int    m2_scaleN[kNEtaBins]   = {0};
        int    m3_scaleN[kNEtaBins]   = {0};
        for (Long64_t i = 0; i < nEntries; ++i) {
            t->GetEntry(i);
            if (i > 0 && i % 2000000 == 0)
                std::cout << "  " << i << " / " << nEntries << std::endl;

            if (!passSelection(photonPt, eta2, mll, mllg,
                               dR_l1_ph, dR_l2_ph, dR_ll,
                               cellSize, truthMatch, true,
                               tightID, isConv))
                continue;
            if (!cellE || (int)cellE->size() != kClusterSize) continue;
            if (!isHealthyCluster(*cellE)) continue;
            if (!cellEta || !cellPhi) continue;

            int etaBin = findEtaBin(std::fabs(eta2));
            if (etaBin < 0) continue;

            // MC event weight: pu_wgt * mc_wgt
            // (matches Francisco NTUP.C line 797)
            double w = muwgt * mcwgt;

            // Total cell energy (no clamping, same as Francisco)
            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cellE->at(k);
            if (Etot <= 0) continue;

            // --- Set A: branch values ---
            h_reta_mc_unf[etaBin]->Fill(reta_unf, w);
            h_rphi_mc_unf[etaBin]->Fill(rphi_unf, w);
            h_weta2_mc_unf[etaBin]->Fill(weta2_unf, w);
            h_reta_mc_fud[etaBin]->Fill(reta_fud, w);
            h_rphi_mc_fud[etaBin]->Fill(rphi_fud, w);
            h_weta2_mc_fud[etaBin]->Fill(weta2_fud, w);

            // --- Set B: cell-computed baseline ---
            h_reta_mc_cell[etaBin]->Fill(
                calcRetaPos(*cellE, *cellEta, *cellPhi), w);
            h_rphi_mc_cell[etaBin]->Fill(
                calcRphiPos(*cellE, *cellEta, *cellPhi), w);
            h_weta2_mc_cell[etaBin]->Fill(
                calcWeta2Pos(*cellE, *cellEta, *cellPhi, eta2, etamax2), w);


            // =============================================================
            // Method 1: flat shift (CalCoef.C + NTUP.C, document §5.1)
            //
            // E_new_k = E_old_k + Delta_k * E_total
            //
            // Since sum(Delta_k) = 0, sum(E_new) = sum(E_old) exactly.
            // No clamping, no renormalization. (This is the key insight
            // from CalCoef.C that our previous version got wrong.)
            // =============================================================
            {
                std::vector<double> m1_cells(kClusterSize);
                for (int k = 0; k < kClusterSize; ++k) {
                    m1_cells[k] = cellE->at(k) + m1_delta[k][etaBin] * Etot;
                }
                h_reta_m1[etaBin]->Fill(
                    calcRetaPos(m1_cells, *cellEta, *cellPhi), w);
                h_rphi_m1[etaBin]->Fill(
                    calcRphiPos(m1_cells, *cellEta, *cellPhi), w);
                h_weta2_m1[etaBin]->Fill(
                    calcWeta2Pos(m1_cells, *cellEta, *cellPhi, eta2, etamax2), w);
            }


            // =============================================================
            // Method 2: DeltaR-matched TProfile (Francisco's tree.C)
            //
            // For each cell k:
            //   alpha = m2_prof->Interpolate(mc_frac)
            //   E_new_k = E_old_k + alpha * E_total
            //
            // The TProfile stores alpha(f_MC) = <f_data - f_MC> from
            // DeltaR-matched data/MC photon pairs.
            // Rescale afterwards to preserve total energy.
            // Falls back to M1 shift if TProfile has < 10 entries in bin.
            // =============================================================
            if (hasM2) {
                std::vector<double> m2_cells(kClusterSize);
                for (int k = 0; k < kClusterSize; ++k) {
                    double mc_frac = cellE->at(k) / Etot;
                    double alpha = 0;
                    if (m2_prof[k][etaBin]) {
                        int bin = m2_prof[k][etaBin]->FindBin(mc_frac);
                        if (m2_prof[k][etaBin]->GetBinEntries(bin) >= 10)
                            alpha = m2_prof[k][etaBin]->Interpolate(mc_frac);
                        else
                            alpha = m1_delta[k][etaBin];  // M1 fallback
                    } else {
                        alpha = m1_delta[k][etaBin];  // M1 fallback
                    }
                    m2_cells[k] = cellE->at(k) + alpha * Etot;
                }

                // Rescale to preserve total energy
                // (tree.C line 219: NewMC->at(k) *= sumNew/sumE_MC)
                double sumNew = 0;
                for (int k = 0; k < kClusterSize; ++k) sumNew += m2_cells[k];
                if (sumNew > 0) {
                    double scale = Etot / sumNew;
                    m2_scaleSum[etaBin] += scale;
                    m2_scaleN[etaBin]++;
                    for (int k = 0; k < kClusterSize; ++k) m2_cells[k] *= scale;
                }

                h_reta_m2[etaBin]->Fill(
                    calcRetaPos(m2_cells, *cellEta, *cellPhi), w);
                h_rphi_m2[etaBin]->Fill(
                    calcRphiPos(m2_cells, *cellEta, *cellPhi), w);
                h_weta2_m2[etaBin]->Fill(
                    calcWeta2Pos(m2_cells, *cellEta, *cellPhi, eta2, etamax2), w);
                nM2Applied++;
            }


            // =============================================================
            // Method 3: shift+stretch (document §5.2)
            //
            // f'_k = (sigma_data / sigma_MC) * (f_MC - mu_MC) + mu_data
            // E'_k = f'_k * E_total
            // Then rescale to preserve total energy.
            // =============================================================
            if (hasM3) {
                std::vector<double> m3_cells(kClusterSize);
                for (int k = 0; k < kClusterSize; ++k) {
                    double mc_frac = cellE->at(k) / Etot;
                    double f_new = m3_sigRatio[k][etaBin]
                                 * (mc_frac - m3_meanMC[k][etaBin])
                                 + m3_meanData[k][etaBin];
                    m3_cells[k] = f_new * Etot;
                }

                // Rescale to preserve total energy
                double sumNew = 0;
                for (int k = 0; k < kClusterSize; ++k) sumNew += m3_cells[k];
                if (sumNew > 0) {
                    double scale = Etot / sumNew;
                    m3_scaleSum[etaBin] += scale;
                    m3_scaleN[etaBin]++;
                    for (int k = 0; k < kClusterSize; ++k) m3_cells[k] *= scale;
                }

                h_reta_m3[etaBin]->Fill(
                    calcRetaPos(m3_cells, *cellEta, *cellPhi), w);
                h_rphi_m3[etaBin]->Fill(
                    calcRphiPos(m3_cells, *cellEta, *cellPhi), w);
                h_weta2_m3[etaBin]->Fill(
                    calcWeta2Pos(m3_cells, *cellEta, *cellPhi, eta2, etamax2), w);
            }

            // Cell fraction accumulators (MC baseline + each method)
            for (int k = 0; k < kClusterSize; ++k)
                prof_mc[k][etaBin] += w * (cellE->at(k) / Etot);
            prof_mc_w[etaBin] += w;

            nUsed++;
        }
        std::cout << "  Selected: " << nUsed << " events" << std::endl;
        if (hasM2) std::cout << "  M2 applied: " << nM2Applied << std::endl;

        // Print renormalization diagnostics
        std::cout << "\n--- Renormalization Scale Factors (avg Etot/sumNew) ---" << std::endl;
        std::cout << "EtaBin |  M2 <scale>  |  M3 <scale>" << std::endl;
        std::cout << "-------|--------------|------------" << std::endl;
        for (int n = 0; n < kNEtaBins; ++n) {
            if (n == 8) continue;
            double avgM2 = (m2_scaleN[n] > 0) ? m2_scaleSum[n] / m2_scaleN[n] : 0;
            double avgM3 = (m3_scaleN[n] > 0) ? m3_scaleSum[n] / m3_scaleN[n] : 0;
            if (m2_scaleN[n] > 0 || m3_scaleN[n] > 0)
                std::cout << Form("%2d     |  %10.6f  |  %10.6f", n, avgM2, avgM3)
                          << std::endl;
        }

        delete t;
    }


    // ==================================================================
    // Normalise and compute chi2
    // ==================================================================
    auto safeNorm = [](TH1D* h) {
        double integral = h->Integral();
        if (integral > 0) h->Scale(1.0 / integral);
    };

    std::cout << "\n=== Chi2: Set A (branch-based) — vs data_br ===" << std::endl;
    std::cout << "EtaBin |  R_eta: unf    fud  |"
              << "  R_phi: unf    fud  |"
              << " w_eta2: unf    fud" << std::endl;
    std::cout << "-------|---------------------|"
              << "---------------------|"
              << "--------------------" << std::endl;

    for (int n = 0; n < kNEtaBins; ++n) {
        if (n == 8) continue;  // crack region — always empty
        if (h_reta_data_br[n]->GetEntries() < 100) continue;

        safeNorm(h_reta_data_br[n]); safeNorm(h_reta_mc_unf[n]);
        safeNorm(h_reta_mc_fud[n]);
        safeNorm(h_rphi_data_br[n]); safeNorm(h_rphi_mc_unf[n]);
        safeNorm(h_rphi_mc_fud[n]);
        safeNorm(h_weta2_data_br[n]); safeNorm(h_weta2_mc_unf[n]);
        safeNorm(h_weta2_mc_fud[n]);

        int ndf;
        std::cout << Form(
            "%6d | %6.2f %6.2f |"
            " %6.2f %6.2f |"
            " %6.2f %6.2f",
            n,
            computeChi2(h_reta_mc_unf[n], h_reta_data_br[n], ndf),
            computeChi2(h_reta_mc_fud[n], h_reta_data_br[n], ndf),
            computeChi2(h_rphi_mc_unf[n], h_rphi_data_br[n], ndf),
            computeChi2(h_rphi_mc_fud[n], h_rphi_data_br[n], ndf),
            computeChi2(h_weta2_mc_unf[n], h_weta2_data_br[n], ndf),
            computeChi2(h_weta2_mc_fud[n], h_weta2_data_br[n], ndf)
            ) << std::endl;
    }

    std::cout << "\n=== Chi2: Set B (cell-computed) — vs data_cell ===" << std::endl;
    std::cout << "EtaBin |  R_eta:  MC    M1    M2    M3  |"
              << "  R_phi:  MC    M1    M2    M3  |"
              << " w_eta2:  MC    M1    M2    M3" << std::endl;
    std::cout << "-------|----------------------------|"
              << "----------------------------|"
              << "----------------------------" << std::endl;

    for (int n = 0; n < kNEtaBins; ++n) {
        if (n == 8) continue;  // crack region — always empty
        if (h_reta_data_cell[n]->GetEntries() < 100) continue;

        safeNorm(h_reta_data_cell[n]); safeNorm(h_reta_mc_cell[n]);
        safeNorm(h_reta_m1[n]); safeNorm(h_reta_m2[n]); safeNorm(h_reta_m3[n]);
        safeNorm(h_rphi_data_cell[n]); safeNorm(h_rphi_mc_cell[n]);
        safeNorm(h_rphi_m1[n]); safeNorm(h_rphi_m2[n]); safeNorm(h_rphi_m3[n]);
        safeNorm(h_weta2_data_cell[n]); safeNorm(h_weta2_mc_cell[n]);
        safeNorm(h_weta2_m1[n]); safeNorm(h_weta2_m2[n]); safeNorm(h_weta2_m3[n]);

        int ndf;
        std::cout << Form(
            "%6d | %6.2f %5.2f %5.2f %5.2f |"
            " %6.2f %5.2f %5.2f %5.2f |"
            " %6.2f %5.2f %5.2f %5.2f",
            n,
            computeChi2(h_reta_mc_cell[n], h_reta_data_cell[n], ndf),
            computeChi2(h_reta_m1[n],      h_reta_data_cell[n], ndf),
            computeChi2(h_reta_m2[n],      h_reta_data_cell[n], ndf),
            computeChi2(h_reta_m3[n],      h_reta_data_cell[n], ndf),
            computeChi2(h_rphi_mc_cell[n], h_rphi_data_cell[n], ndf),
            computeChi2(h_rphi_m1[n],      h_rphi_data_cell[n], ndf),
            computeChi2(h_rphi_m2[n],      h_rphi_data_cell[n], ndf),
            computeChi2(h_rphi_m3[n],      h_rphi_data_cell[n], ndf),
            computeChi2(h_weta2_mc_cell[n], h_weta2_data_cell[n], ndf),
            computeChi2(h_weta2_m1[n],      h_weta2_data_cell[n], ndf),
            computeChi2(h_weta2_m2[n],      h_weta2_data_cell[n], ndf),
            computeChi2(h_weta2_m3[n],      h_weta2_data_cell[n], ndf)
            ) << std::endl;
    }

    // ==================================================================
    // Write output
    // ==================================================================
    TFile* fout = TFile::Open(histoFile, "RECREATE");
    for (int n = 0; n < kNEtaBins; ++n) {
        // Set A
        h_reta_data_br[n]->Write(); h_reta_mc_unf[n]->Write();
        h_reta_mc_fud[n]->Write();
        h_rphi_data_br[n]->Write(); h_rphi_mc_unf[n]->Write();
        h_rphi_mc_fud[n]->Write();
        h_weta2_data_br[n]->Write(); h_weta2_mc_unf[n]->Write();
        h_weta2_mc_fud[n]->Write();
        // Set B
        h_reta_data_cell[n]->Write(); h_reta_mc_cell[n]->Write();
        h_reta_m1[n]->Write(); h_reta_m2[n]->Write(); h_reta_m3[n]->Write();
        h_rphi_data_cell[n]->Write(); h_rphi_mc_cell[n]->Write();
        h_rphi_m1[n]->Write(); h_rphi_m2[n]->Write(); h_rphi_m3[n]->Write();
        h_weta2_data_cell[n]->Write(); h_weta2_mc_cell[n]->Write();
        h_weta2_m1[n]->Write(); h_weta2_m2[n]->Write(); h_weta2_m3[n]->Write();
    }

    // Cell profile maps: mean fraction per cell from validation loops
    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* h_profData = new TH2D(
            "profData_" + suffix, "<f_{data}> (validate);phi;eta",
            kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);
        TH2D* h_profMC = new TH2D(
            "profMC_" + suffix, "<f_{MC}> (validate);phi;eta",
            kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);

        for (int phi = 1; phi <= kPhiSize; ++phi) {
            for (int eta = 1; eta <= kEtaSize; ++eta) {
                int k = (phi - 1) + kPhiSize * (eta - 1);
                double dw = prof_data_w[n];
                double mw = prof_mc_w[n];
                h_profData->SetBinContent(phi, eta,
                    dw > 0 ? prof_data[k][n] / dw : 0);
                h_profMC->SetBinContent(phi, eta,
                    mw > 0 ? prof_mc[k][n] / mw : 0);
            }
        }
        h_profData->Write();
        h_profMC->Write();
    }

    fout->Close(); delete fout;

    std::cout << "\nHistograms written to: " << histoFile << std::endl;
    fcorr->Close(); delete fcorr;
    return 0;
}
