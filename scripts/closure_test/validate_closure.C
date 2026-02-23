///////////////////////////////////////////////////////////////////////////////
// validate_closure.C
//
// Fills shower shape histograms from all files and computes closure metrics.
// Creates histograms per eta bin for:
//   - Original MC
//   - Pseudo-data (the target)
//   - Reweighted MC method 1 (flat shift, §5.1)
//   - Reweighted MC method 2 (TProfile, §5.2)
//   - (Optional) Oracle method 0
//
// Also computes and stores ratio histograms (reweighted / pseudo-data).
//
// Usage:
//   root -l -b -q 'validate_closure.C("mc.root","out/","out/histos.root")'
//   root -l -b -q 'validate_closure.C("mc.root","out/","out/histos.root",200000)'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace closure;

// ======================================================================
// Fill shower shape histograms from a file
// ======================================================================
void fillHistos(const char* filename, Long64_t maxEvents,
                TH1D* h_reta[], TH1D* h_rphi[], TH1D* h_weta2[],
                const char* tag) {

    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << filename << std::endl;
        return;
    }
    TTree* t = (TTree*)f->Get(kTreeName);
    if (!t) {
        std::cerr << "ERROR: No tree in " << filename << std::endl;
        f->Close();
        return;
    }

    std::vector<double>* cellE = nullptr;
    Int_t cellSize = 0;
    Float_t eta2 = 0;

    t->SetBranchAddress(kCellBranch, &cellE);
    t->SetBranchAddress(kCellSizeBranch, &cellSize);
    t->SetBranchAddress(kEtaBranch, &eta2);

    // Only read the branches we need (huge speedup)
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus(kCellBranch, 1);
    t->SetBranchStatus(kCellSizeBranch, 1);
    t->SetBranchStatus(kEtaBranch, 1);

    Long64_t nEntries = t->GetEntries();
    if (maxEvents > 0 && maxEvents < nEntries)
        nEntries = maxEvents;

    std::cout << "  " << tag << ": " << filename << " (" << nEntries << " entries)" << std::endl;

    Long64_t nUsed = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);

        if (cellSize != kClusterSize) continue;
        if (!cellE || (int)cellE->size() != kClusterSize) continue;

        int etaBin = findEtaBin(std::fabs(eta2));
        if (etaBin < 0) continue;

        double Etot = 0;
        for (int k = 0; k < kClusterSize; ++k)
            Etot += cellE->at(k);
        if (Etot <= 0) continue;

        h_reta[etaBin]->Fill(calcReta(*cellE));
        h_rphi[etaBin]->Fill(calcRphi(*cellE));
        h_weta2[etaBin]->Fill(calcWeta2(*cellE));
        nUsed++;
    }

    std::cout << "    Used: " << nUsed << " events" << std::endl;

    f->Close();
    delete f;
}

// ======================================================================
// Compute chi2/ndf between two normalised histograms
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
int validate_closure(const char* mcFile,
                     const char* outputDir,
                     const char* histoFile,
                     Long64_t maxEvents = -1) {

    // Build file paths from output directory
    TString dir(outputDir);
    if (!dir.EndsWith("/")) dir += "/";
    TString pseudoFile = dir + "pseudo.root";
    TString m1File     = dir + "reweighted_m1.root";
    TString m2File     = dir + "reweighted_m2.root";

    // ======================================================================
    // Create histograms: 4 samples x 3 variables x 14 eta bins
    // ======================================================================
    TH1D* h_reta_mc[kNEtaBins],    *h_reta_ps[kNEtaBins],
         *h_reta_m1[kNEtaBins],    *h_reta_m2[kNEtaBins];
    TH1D* h_rphi_mc[kNEtaBins],    *h_rphi_ps[kNEtaBins],
         *h_rphi_m1[kNEtaBins],    *h_rphi_m2[kNEtaBins];
    TH1D* h_weta2_mc[kNEtaBins],   *h_weta2_ps[kNEtaBins],
         *h_weta2_m1[kNEtaBins],   *h_weta2_m2[kNEtaBins];

    for (int n = 0; n < kNEtaBins; ++n) {
        TString s = Form("_%d", n);

        h_reta_mc[n]  = new TH1D("reta_mc" + s,  "", 25, 0.90, 1.0);
        h_reta_ps[n]  = new TH1D("reta_ps" + s,  "", 25, 0.90, 1.0);
        h_reta_m1[n]  = new TH1D("reta_m1" + s,  "", 25, 0.90, 1.0);
        h_reta_m2[n]  = new TH1D("reta_m2" + s,  "", 25, 0.90, 1.0);

        h_rphi_mc[n]  = new TH1D("rphi_mc" + s,  "", 25, 0.85, 1.0);
        h_rphi_ps[n]  = new TH1D("rphi_ps" + s,  "", 25, 0.85, 1.0);
        h_rphi_m1[n]  = new TH1D("rphi_m1" + s,  "", 25, 0.85, 1.0);
        h_rphi_m2[n]  = new TH1D("rphi_m2" + s,  "", 25, 0.85, 1.0);

        h_weta2_mc[n] = new TH1D("weta2_mc" + s, "", 25, 0.008, 0.016);
        h_weta2_ps[n] = new TH1D("weta2_ps" + s, "", 25, 0.008, 0.016);
        h_weta2_m1[n] = new TH1D("weta2_m1" + s, "", 25, 0.008, 0.016);
        h_weta2_m2[n] = new TH1D("weta2_m2" + s, "", 25, 0.008, 0.016);
    }

    // ======================================================================
    // Fill histograms from all files
    // ======================================================================
    std::cout << "=== Filling histograms ===" << std::endl;
    fillHistos(mcFile,            maxEvents, h_reta_mc, h_rphi_mc, h_weta2_mc, "MC");
    fillHistos(pseudoFile.Data(), maxEvents, h_reta_ps, h_rphi_ps, h_weta2_ps, "Pseudo");
    fillHistos(m1File.Data(),     maxEvents, h_reta_m1, h_rphi_m1, h_weta2_m1, "Method1");
    fillHistos(m2File.Data(),     maxEvents, h_reta_m2, h_rphi_m2, h_weta2_m2, "Method2");

    // ======================================================================
    // Normalise and compute chi2
    // ======================================================================
    std::cout << "\n=== Closure metrics (chi2/ndf) ===" << std::endl;
    std::cout << "EtaBin |   R_eta M1   R_eta M2 |   R_phi M1   R_phi M2 |  w_eta2 M1  w_eta2 M2" << std::endl;
    std::cout << "-------|----------------------|----------------------|----------------------" << std::endl;

    for (int n = 0; n < kNEtaBins; ++n) {
        if (h_reta_ps[n]->GetEntries() < 100) continue;

        // Normalise all histograms
        h_reta_mc[n]->Scale(1.0 / h_reta_mc[n]->Integral());
        h_reta_ps[n]->Scale(1.0 / h_reta_ps[n]->Integral());
        h_reta_m1[n]->Scale(1.0 / h_reta_m1[n]->Integral());
        h_reta_m2[n]->Scale(1.0 / h_reta_m2[n]->Integral());

        h_rphi_mc[n]->Scale(1.0 / h_rphi_mc[n]->Integral());
        h_rphi_ps[n]->Scale(1.0 / h_rphi_ps[n]->Integral());
        h_rphi_m1[n]->Scale(1.0 / h_rphi_m1[n]->Integral());
        h_rphi_m2[n]->Scale(1.0 / h_rphi_m2[n]->Integral());

        h_weta2_mc[n]->Scale(1.0 / h_weta2_mc[n]->Integral());
        h_weta2_ps[n]->Scale(1.0 / h_weta2_ps[n]->Integral());
        h_weta2_m1[n]->Scale(1.0 / h_weta2_m1[n]->Integral());
        h_weta2_m2[n]->Scale(1.0 / h_weta2_m2[n]->Integral());

        // Compare reweighted MC to pseudo-data (not MC)
        int ndf;
        double chi2_reta_m1 = computeChi2(h_reta_m1[n], h_reta_ps[n], ndf);
        double chi2_reta_m2 = computeChi2(h_reta_m2[n], h_reta_ps[n], ndf);
        double chi2_rphi_m1 = computeChi2(h_rphi_m1[n], h_rphi_ps[n], ndf);
        double chi2_rphi_m2 = computeChi2(h_rphi_m2[n], h_rphi_ps[n], ndf);
        double chi2_weta2_m1 = computeChi2(h_weta2_m1[n], h_weta2_ps[n], ndf);
        double chi2_weta2_m2 = computeChi2(h_weta2_m2[n], h_weta2_ps[n], ndf);

        std::cout << Form("%6d | %10.4f %10.4f | %10.4f %10.4f | %10.4f %10.4f",
                          n, chi2_reta_m1, chi2_reta_m2,
                          chi2_rphi_m1, chi2_rphi_m2,
                          chi2_weta2_m1, chi2_weta2_m2) << std::endl;
    }

    // ======================================================================
    // Write histograms
    // ======================================================================
    TFile* fout = TFile::Open(histoFile, "RECREATE");

    for (int n = 0; n < kNEtaBins; ++n) {
        h_reta_mc[n]->Write();   h_reta_ps[n]->Write();
        h_reta_m1[n]->Write();   h_reta_m2[n]->Write();
        h_rphi_mc[n]->Write();   h_rphi_ps[n]->Write();
        h_rphi_m1[n]->Write();   h_rphi_m2[n]->Write();
        h_weta2_mc[n]->Write();  h_weta2_ps[n]->Write();
        h_weta2_m1[n]->Write();  h_weta2_m2[n]->Write();
    }

    std::cout << "\nHistograms written to: " << histoFile << std::endl;

    fout->Close();
    delete fout;

    return 0;
}
