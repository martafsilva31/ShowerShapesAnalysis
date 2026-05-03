///////////////////////////////////////////////////////////////////////////////
// create_pseudodata.C
//
// Creates pseudo-data by applying a systematic per-cell bias plus Gaussian
// random noise to MC23e cell energies. The systematic component (fixed across
// events) gives the reweighting method a signal to correct; the random noise
// (independent per event per cell) adds realistic scatter.
//
// For each healthy event, for each cell k:
//   e'_k = e_k * (1 + distortionLevel * biasPattern_k + Gaus(0, noiseSigma))
//   e'_k *= E_total / sum(e'_j)   (conserve cluster energy)
//
// biasPattern_k = -1 (central), -0.4 (ring1), +0.4 (ring2), +1.0 (outer)
//   => broader showers, mimicking the standard ATLAS data-MC difference
//      (MC has narrower showers than data; fudge factors broaden MC).
//
// distortionLevel  controls the systematic bias magnitude.
// noiseSigma       controls random event-by-event scatter (default: level/2).
//
// Usage:
//   root -l -b -q 'create_pseudodata.C("in.root","ps.root",500000,0.05)'
//   root -l -b -q 'create_pseudodata.C("in.root","ps.root",500000,0.05,0.02)'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace closure;

int create_pseudodata(const char* inputFile,
                      const char* outputFile,
                      Long64_t maxEvents = -1,
                      double distortionLevel = 0.05,
                      double noiseSigma = -1) {

    // Default noise sigma = half the systematic level
    if (noiseSigma < 0) noiseSigma = distortionLevel / 2.0;

    // ======================================================================
    // Open input
    // ======================================================================
    TFile* fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file: " << inputFile << std::endl;
        return 1;
    }
    TTree* tin = (TTree*)fin->Get(kTreeName);
    if (!tin) {
        std::cerr << "ERROR: Cannot find tree '" << kTreeName << "' in " << inputFile << std::endl;
        return 1;
    }

    Long64_t nEntries = tin->GetEntries();
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Total entries: " << nEntries << std::endl;
    if (maxEvents > 0 && maxEvents < nEntries) {
        nEntries = maxEvents;
        std::cout << "Processing only " << nEntries << " entries" << std::endl;
    }

    // ======================================================================
    // Set up branches
    // ======================================================================
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

    // Random number generator for Gaussian noise
    TRandom3 rng(42);

    std::cout << "\n--- Systematic + noise distortion model ---" << std::endl;
    std::cout << "  Systematic bias level: " << distortionLevel << std::endl;
    std::cout << "  Random noise sigma:    " << noiseSigma << std::endl;
    std::cout << "  Formula: e' = e * (1 + " << distortionLevel << "*bias_k + Gaus(0," << noiseSigma << "))" << std::endl;
    std::cout << "  Bias pattern: center=-1, ring1=-0.4, ring2=+0.4, outer=+1" << std::endl;

    // ======================================================================
    // Create output
    // ======================================================================
    TFile* fout = TFile::Open(outputFile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: Cannot create output file: " << outputFile << std::endl;
        return 1;
    }

    TTree* tout = tin->CloneTree(0);
    tout->SetDirectory(fout);

    // Monitoring histograms
    TH1D* h_reta_orig   = new TH1D("h_reta_orig",   ";R_{#eta};Events", 100, 0.85, 1.0);
    TH1D* h_reta_pseudo = new TH1D("h_reta_pseudo", ";R_{#eta};Events", 100, 0.85, 1.0);
    TH1D* h_rphi_orig   = new TH1D("h_rphi_orig",   ";R_{#phi};Events", 100, 0.80, 1.0);
    TH1D* h_rphi_pseudo = new TH1D("h_rphi_pseudo", ";R_{#phi};Events", 100, 0.80, 1.0);
    TH1D* h_weta2_orig  = new TH1D("h_weta2_orig",  ";w_{#eta 2};Events", 100, 0.005, 0.020);
    TH1D* h_weta2_pseudo= new TH1D("h_weta2_pseudo",";w_{#eta 2};Events", 100, 0.005, 0.020);

    Long64_t nProcessed = 0, nDistorted = 0, nPassthrough = 0;

    // ======================================================================
    // Event loop: apply per-cell scaling
    // ======================================================================
    std::cout << "\n=== Processing events ===" << std::endl;

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

        // Clamp negative cell energies and check cluster quality
        clampCellEnergies(*cellE);
        if (!isHealthyCluster(*cellE)) {
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

        // Fill original shower shape histos
        h_reta_orig->Fill(calcReta(*cellE));
        h_rphi_orig->Fill(calcRphi(*cellE));
        h_weta2_orig->Fill(calcWeta2(*cellE));

        // Apply systematic bias + random noise per cell
        double sumScaled = 0;
        for (int k = 0; k < kClusterSize; ++k) {
            double bias = distortionLevel * getBiasPattern(k);
            double noise = rng.Gaus(0, noiseSigma);
            cellE->at(k) *= (1.0 + bias + noise);
            sumScaled += cellE->at(k);
        }

        // Rescale to preserve total cluster energy
        if (sumScaled > 0) {
            double rescale = Etot / sumScaled;
            for (int k = 0; k < kClusterSize; ++k)
                cellE->at(k) *= rescale;
        }

        // Fill pseudo-data shower shape histos
        h_reta_pseudo->Fill(calcReta(*cellE));
        h_rphi_pseudo->Fill(calcRphi(*cellE));
        h_weta2_pseudo->Fill(calcWeta2(*cellE));

        nDistorted++;
        tout->Fill();
    }

    // ======================================================================
    // Write output
    // ======================================================================
    fout->cd();
    tout->Write();
    h_reta_orig->Write();
    h_reta_pseudo->Write();
    h_rphi_orig->Write();
    h_rphi_pseudo->Write();
    h_weta2_orig->Write();
    h_weta2_pseudo->Write();

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Distortion level: " << distortionLevel << std::endl;
    std::cout << "Total processed:  " << nProcessed << std::endl;
    std::cout << "Distorted:        " << nDistorted << std::endl;
    std::cout << "Passthrough:      " << nPassthrough << std::endl;
    std::cout << "R_eta MC mean:    " << h_reta_orig->GetMean() << std::endl;
    std::cout << "R_eta PS mean:    " << h_reta_pseudo->GetMean() << std::endl;
    std::cout << "Output written to: " << outputFile << std::endl;

    fout->Close();
    fin->Close();
    delete fout;
    delete fin;

    return 0;
}
