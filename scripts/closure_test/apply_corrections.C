///////////////////////////////////////////////////////////////////////////////
// apply_corrections.C
//
// Applies cell-energy reweighting corrections to MC, producing reweighted MC
// that should match the pseudo-data distributions.
//
// Two methods available:
//
// METHOD 1 — Flat shift (Shift Only):
//   f_corr_k = f_MC_k + shift_k
//   where shift_k = <f_pseudo>_k - <f_MC>_k (flat per cell per eta bin).
//
// METHOD 2 — TProfile correction (Shift + Stretch):
//   f_corr_k = f_MC_k + alpha_k(f_MC_k)
//   where alpha_k is interpolated from the TProfile of (f_pseudo - f_MC) vs f_MC.
//   This captures the energy-dependent linear relationship.
//
// Then cluster energy is conserved by rescaling:
//   e_corr_k = f_corr_k / sum(f_corr_j) * E_total
//
// Reference: ATL-COM-PHYS-2021-640, Section 5
//
// Usage:
//   root -l -b -q 'apply_corrections.C("mc.root","corrections.root","rew_m1.root",-1,1)'  # Shift
//   root -l -b -q 'apply_corrections.C("mc.root","corrections.root","rew_m2.root",-1,2)'  # TProfile
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace closure;

int apply_corrections(const char* inputFile,
                      const char* corrFile,
                      const char* outputFile,
                      Long64_t maxEvents = -1,
                      int method = 2) {

    std::cout << "=== Apply Corrections ===" << std::endl;
    std::cout << "Method: " << method;
    if (method == 1) std::cout << " (Flat shift, Shift Only)";
    else if (method == 2) std::cout << " (TProfile, Shift + Stretch)";
    std::cout << std::endl;

    // ======================================================================
    // Load corrections (skip for oracle method 0)
    // ======================================================================
    TFile* fcorr = nullptr;
    TProfile* CellCorrection[kClusterSize][kNEtaBins] = {};
    double flatShift[kClusterSize][kNEtaBins] = {};

    if (method > 0) {
        fcorr = TFile::Open(corrFile, "READ");
        if (!fcorr || fcorr->IsZombie()) {
            std::cerr << "ERROR: Cannot open corrections file: " << corrFile << std::endl;
            return 1;
        }

        for (int n = 0; n < kNEtaBins; ++n) {
            TString suffix = Form("Eta_%1.2f_%1.2f", kEtaLimits[n], kEtaLimits[n + 1]);

            if (method == 1) {
                // Load flat shifts
                TH2D* h_shift = (TH2D*)fcorr->Get("shift_" + suffix);
                if (!h_shift) {
                    std::cerr << "ERROR: Missing shift histogram for eta bin " << n << std::endl;
                    return 1;
                }
                for (int k = 0; k < kClusterSize; ++k) {
                    int phi_bin = k % kPhiSize;
                    int eta_bin = k / kPhiSize;
                    flatShift[k][n] = h_shift->GetBinContent(phi_bin + 1, eta_bin + 1);
                }
            }
            else if (method == 2) {
                // Load TProfiles
                for (int k = 0; k < kClusterSize; ++k) {
                    CellCorrection[k][n] = (TProfile*)fcorr->Get(
                        Form("CellCorrection_%d_%s", k, suffix.Data()));
                    if (!CellCorrection[k][n]) {
                        std::cerr << "ERROR: Missing TProfile for cell " << k
                                  << ", eta bin " << n << std::endl;
                        return 1;
                    }
                }
            }
        }
        std::cout << "Corrections loaded from: " << corrFile << std::endl;
    }

    // ======================================================================
    // Open input (MC)
    // ======================================================================
    TFile* fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file: " << inputFile << std::endl;
        return 1;
    }
    TTree* tin = (TTree*)fin->Get(kTreeName);
    if (!tin) {
        std::cerr << "ERROR: No tree '" << kTreeName << "' in " << inputFile << std::endl;
        return 1;
    }

    Long64_t nEntries = tin->GetEntries();
    if (maxEvents > 0 && maxEvents < nEntries)
        nEntries = maxEvents;
    std::cout << "Input: " << inputFile << " (" << nEntries << " entries)" << std::endl;

    std::vector<double>* cellE = nullptr;
    Int_t cellSize = 0;
    Float_t eta2 = 0;

    tin->SetBranchAddress(kCellBranch, &cellE);
    tin->SetBranchAddress(kCellSizeBranch, &cellSize);
    tin->SetBranchAddress(kEtaBranch, &eta2);

    // Only read the branches we need (huge speedup)
    tin->SetBranchStatus("*", 0);
    tin->SetBranchStatus(kCellBranch, 1);
    tin->SetBranchStatus(kCellSizeBranch, 1);
    tin->SetBranchStatus(kEtaBranch, 1);

    // ======================================================================
    // Create output
    // ======================================================================
    TFile* fout = TFile::Open(outputFile, "RECREATE");
    TTree* tout = tin->CloneTree(0);
    tout->SetDirectory(fout);

    Long64_t nProcessed = 0, nCorrected = 0, nPassthrough = 0;

    // ======================================================================
    // Event loop
    // ======================================================================
    for (Long64_t i = 0; i < nEntries; ++i) {
        tin->GetEntry(i);
        if (i % 1000000 == 0)
            std::cout << "  " << i << " / " << nEntries << std::endl;
        nProcessed++;

        // Passthrough for bad clusters
        if (cellSize != kClusterSize || !cellE ||
            (int)cellE->size() != kClusterSize) {
            nPassthrough++;
            tout->Fill();
            continue;
        }

        int etaBin = findEtaBin(std::fabs(eta2));
        if (etaBin < 0) {
            nPassthrough++;
            tout->Fill();
            continue;
        }

        double Etot = 0;
        for (int k = 0; k < kClusterSize; ++k)
            Etot += cellE->at(k);
        if (Etot <= 0) {
            nPassthrough++;
            tout->Fill();
            continue;
        }

        // Apply correction (FORWARD: MC → pseudo-data direction)
        std::vector<double> fCorr(kClusterSize);
        double sumCorr = 0;

        for (int k = 0; k < kClusterSize; ++k) {
            double fMC = cellE->at(k) / Etot;

            if (method == 1) {
                // Flat shift: f_corr = f_MC + shift
                fCorr[k] = fMC + flatShift[k][etaBin];
            }
            else if (method == 2) {
                // TProfile: f_corr = f_MC + alpha(f_MC)
                // Use Interpolate() exactly like tree.C line 197
                double alpha = CellCorrection[k][etaBin]->Interpolate(fMC);
                fCorr[k] = fMC + alpha;
            }

            sumCorr += fCorr[k];
        }

        // Rescale to conserve cluster energy
        if (std::fabs(sumCorr) > 1e-12) {
            for (int k = 0; k < kClusterSize; ++k)
                cellE->at(k) = fCorr[k] / sumCorr * Etot;
        }

        nCorrected++;
        tout->Fill();
    }

    fout->cd();
    tout->Write();

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total processed:  " << nProcessed << std::endl;
    std::cout << "Corrected:        " << nCorrected << std::endl;
    std::cout << "Passthrough:      " << nPassthrough << std::endl;
    std::cout << "Output written to: " << outputFile << std::endl;

    fout->Close();
    if (fcorr) fcorr->Close();
    fin->Close();
    delete fout;
    if (fcorr) delete fcorr;
    delete fin;

    return 0;
}
