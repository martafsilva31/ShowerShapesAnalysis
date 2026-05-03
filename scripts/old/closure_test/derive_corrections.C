///////////////////////////////////////////////////////////////////////////////
// derive_corrections.C
//
// Derives cell-energy reweighting corrections by comparing MC and pseudo-data.
// Implements TWO methods:
//
// METHOD 1 (§5.1) — Flat shift per cell:
//   shift_k = <f_pseudo>_k - <f_MC>_k
//   Stored in TH2D (11 x 7 bins for phi x eta).
//
// METHOD 2 (§5.2) — Energy-dependent correction via TProfile:
//   For each cell k, fill TProfile with:
//     x = f_MC (MC cell fraction)
//     y = f_pseudo - f_MC (difference)
//   This captures the linear relationship exactly as in supervisor's tree.C.
//   When applied, we interpolate: f_corr = f_MC + alpha(f_MC)
//
// Reference: ATL-COM-PHYS-2021-640, Section 5.2
//
// Usage:
//   root -l -b -q 'derive_corrections.C("mc.root", "pseudo.root", "corrections.root")'
//   root -l -b -q 'derive_corrections.C("mc.root", "pseudo.root", "corrections.root", 200000)'
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

using namespace closure;

int derive_corrections(const char* mcFile,
                       const char* pseudoFile,
                       const char* outputFile,
                       Long64_t maxEvents = -1) {

    // ======================================================================
    // Open MC file
    // ======================================================================
    TFile* fmc = TFile::Open(mcFile, "READ");
    if (!fmc || fmc->IsZombie()) {
        std::cerr << "ERROR: Cannot open MC file: " << mcFile << std::endl;
        return 1;
    }
    TTree* tmc = (TTree*)fmc->Get(kTreeName);
    if (!tmc) {
        std::cerr << "ERROR: No tree '" << kTreeName << "' in " << mcFile << std::endl;
        return 1;
    }

    // ======================================================================
    // Open pseudo-data file
    // ======================================================================
    TFile* fps = TFile::Open(pseudoFile, "READ");
    if (!fps || fps->IsZombie()) {
        std::cerr << "ERROR: Cannot open pseudo-data file: " << pseudoFile << std::endl;
        return 1;
    }
    TTree* tps = (TTree*)fps->Get(kTreeName);
    if (!tps) {
        std::cerr << "ERROR: No tree '" << kTreeName << "' in " << pseudoFile << std::endl;
        return 1;
    }

    // Use the minimum number of entries between MC and pseudo-data
    Long64_t nEntriesMC = tmc->GetEntries();
    Long64_t nEntriesPS = tps->GetEntries();
    Long64_t nEntries = std::min(nEntriesMC, nEntriesPS);

    if (maxEvents > 0 && maxEvents < nEntries)
        nEntries = maxEvents;

    std::cout << "MC file:     " << mcFile << " (" << nEntriesMC << " entries)" << std::endl;
    std::cout << "Pseudo file: " << pseudoFile << " (" << nEntriesPS << " entries)" << std::endl;
    std::cout << "Processing:  " << nEntries << " entries" << std::endl;

    // ======================================================================
    // Set up branches
    // ======================================================================
    std::vector<double>* cellMC = nullptr;
    std::vector<double>* cellPS = nullptr;
    Int_t sizeMC = 0, sizePS = 0;
    Float_t eta2MC = 0, eta2PS = 0;

    tmc->SetBranchAddress(kCellBranch, &cellMC);
    tmc->SetBranchAddress(kCellSizeBranch, &sizeMC);
    tmc->SetBranchAddress(kEtaBranch, &eta2MC);

    // Only read the branches we need (huge speedup)
    tmc->SetBranchStatus("*", 0);
    tmc->SetBranchStatus(kCellBranch, 1);
    tmc->SetBranchStatus(kCellSizeBranch, 1);
    tmc->SetBranchStatus(kEtaBranch, 1);

    tps->SetBranchAddress(kCellBranch, &cellPS);
    tps->SetBranchAddress(kCellSizeBranch, &sizePS);
    tps->SetBranchAddress(kEtaBranch, &eta2PS);

    tps->SetBranchStatus("*", 0);
    tps->SetBranchStatus(kCellBranch, 1);
    tps->SetBranchStatus(kCellSizeBranch, 1);
    tps->SetBranchStatus(kEtaBranch, 1);

    // ======================================================================
    // Create output histograms
    // ======================================================================
    TFile* fout = TFile::Open(outputFile, "RECREATE");

    // Method 2: TProfiles (energy-dependent correction, §5.2)
    // CellCorrection[k][n]: x = f_MC, y = f_pseudo - f_MC
    TProfile* CellCorrection[kClusterSize][kNEtaBins];

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f", kEtaLimits[n], kEtaLimits[n + 1]);
        for (int k = 0; k < kClusterSize; ++k) {
            CellCorrection[k][n] = new TProfile(
                Form("CellCorrection_%d_%s", k, suffix.Data()),
                Form("Cell %d correction;f_{MC};f_{pseudo} - f_{MC}", k),
                kNFracBins, kFracBins);
        }
    }

    // Method 1: Flat shifts (§5.1)
    // We'll compute mean fractions and store shifts in TH2D
    TH2D* h_meanMC[kNEtaBins];
    TH2D* h_meanPS[kNEtaBins];
    TH2D* h_shift[kNEtaBins];

    // Accumulators for online mean calculation
    double sumFracMC[kClusterSize][kNEtaBins] = {};
    double sumFracPS[kClusterSize][kNEtaBins] = {};
    Long64_t counts[kNEtaBins] = {};

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f", kEtaLimits[n], kEtaLimits[n + 1]);
        h_meanMC[n] = new TH2D("meanMC_" + suffix, "Mean fraction (MC);phi;eta",
                               kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);
        h_meanPS[n] = new TH2D("meanPS_" + suffix, "Mean fraction (pseudo);phi;eta",
                               kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);
        h_shift[n]  = new TH2D("shift_" + suffix, "Shift (pseudo - MC);phi;eta",
                               kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);
    }

    // ======================================================================
    // Event loop: fill TProfiles and accumulate means
    // ======================================================================
    std::cout << "\n=== Processing events ===" << std::endl;

    Long64_t nProcessed = 0, nUsed = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tmc->GetEntry(i);
        tps->GetEntry(i);

        if (i % 1000000 == 0)
            std::cout << "  " << i << " / " << nEntries << std::endl;
        nProcessed++;

        // Skip if either cluster is bad
        if (sizeMC != kClusterSize || sizePS != kClusterSize) continue;
        if (!cellMC || !cellPS) continue;
        if ((int)cellMC->size() != kClusterSize || (int)cellPS->size() != kClusterSize) continue;

        int etaBin = findEtaBin(std::fabs(eta2MC));
        if (etaBin < 0) continue;

        // Clamp negative cell energies and check cluster quality
        clampCellEnergies(*cellMC);
        clampCellEnergies(*cellPS);
        if (!isHealthyCluster(*cellMC) || !isHealthyCluster(*cellPS)) continue;

        // Compute total energies
        double EtotMC = 0, EtotPS = 0;
        for (int k = 0; k < kClusterSize; ++k) {
            EtotMC += cellMC->at(k);
            EtotPS += cellPS->at(k);
        }
        if (EtotMC <= 0 || EtotPS <= 0) continue;

        // Fill TProfiles and accumulate means
        for (int k = 0; k < kClusterSize; ++k) {
            double fMC = cellMC->at(k) / EtotMC;
            double fPS = cellPS->at(k) / EtotPS;

            // TProfile: x = f_MC, y = f_pseudo - f_MC (exactly like tree.C)
            CellCorrection[k][etaBin]->Fill(fMC, fPS - fMC);

            // Accumulate for flat means
            sumFracMC[k][etaBin] += fMC;
            sumFracPS[k][etaBin] += fPS;
        }
        counts[etaBin]++;
        nUsed++;
    }

    // ======================================================================
    // Compute flat shifts (Method 1)
    // ======================================================================
    std::cout << "\n=== Computing flat shifts ===" << std::endl;

    for (int n = 0; n < kNEtaBins; ++n) {
        if (counts[n] == 0) continue;

        for (int k = 0; k < kClusterSize; ++k) {
            int phi_bin = k % kPhiSize;
            int eta_bin = k / kPhiSize;

            double meanMC = sumFracMC[k][n] / counts[n];
            double meanPS = sumFracPS[k][n] / counts[n];
            double shift  = meanPS - meanMC;

            h_meanMC[n]->SetBinContent(phi_bin + 1, eta_bin + 1, meanMC);
            h_meanPS[n]->SetBinContent(phi_bin + 1, eta_bin + 1, meanPS);
            h_shift[n]->SetBinContent(phi_bin + 1, eta_bin + 1, shift);
        }

        // Print central cell info
        double centralShift = h_shift[n]->GetBinContent(6, 4);  // phi=5, eta=3 (1-indexed)
        std::cout << Form("  Eta bin %2d [%.2f, %.2f): %8lld events, central shift = %+.6f",
                          n, kEtaLimits[n], kEtaLimits[n + 1], counts[n], centralShift) << std::endl;
    }

    // ======================================================================
    // Write output
    // ======================================================================
    fout->cd();

    for (int n = 0; n < kNEtaBins; ++n) {
        for (int k = 0; k < kClusterSize; ++k) {
            CellCorrection[k][n]->Write();
        }
        h_meanMC[n]->Write();
        h_meanPS[n]->Write();
        h_shift[n]->Write();
    }

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total processed: " << nProcessed << std::endl;
    std::cout << "Used for corr:   " << nUsed << std::endl;
    std::cout << "Output written to: " << outputFile << std::endl;

    fout->Close();
    fps->Close();
    fmc->Close();
    delete fout;
    delete fps;
    delete fmc;

    return 0;
}
