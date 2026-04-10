///////////////////////////////////////////////////////////////////////////////
// compare_computed_vs_branch.C
//
// Diagnostic: overlays cell-computed shower shapes (from 7x11 cluster) vs
// stored branch values (official ATLAS computation), for both data and MC.
//
// Produces 4 histograms per variable per eta bin:
//   - Data branch (unfudged = fudged for data)
//   - Data cell-computed (no clamping)
//   - MC branch unfudged
//   - MC cell-computed (no clamping)
//
// Usage:
//   root -l -b -q 'compare_computed_vs_branch.C("data.root","mc.root","outdir/")'
//   root -l -b -q 'compare_computed_vs_branch.C("data.root","mc.root","outdir/",100000)'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace datamc;

int compare_computed_vs_branch(const char* dataFile,
                               const char* mcFile,
                               const char* outDir,
                               Long64_t maxEvents = -1) {

    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);

    TString dir(outDir);
    if (!dir.EndsWith("/")) dir += "/";
    gSystem->mkdir(dir, true);

    // ==================================================================
    // Create histograms: 4 samples × 3 vars × 14 eta bins
    //   data_br, data_cell, mc_br, mc_cell
    // ==================================================================
    TH1D* h_reta[4][kNEtaBins];
    TH1D* h_rphi[4][kNEtaBins];
    TH1D* h_weta2[4][kNEtaBins];

    const char* tags[4] = {"data_br", "data_cell", "mc_br", "mc_cell"};

    for (int s = 0; s < 4; ++s) {
        for (int n = 0; n < kNEtaBins; ++n) {
            TString suf = Form("_%s_%d", tags[s], n);
            h_reta[s][n]  = new TH1D("reta" + suf,  "", 40, 0.85, 1.0);
            h_rphi[s][n]  = new TH1D("rphi" + suf,  "", 40, 0.80, 1.0);
            h_weta2[s][n] = new TH1D("weta2" + suf, "", 40, 0.006, 0.018);
            h_reta[s][n]->Sumw2();
            h_rphi[s][n]->Sumw2();
            h_weta2[s][n]->Sumw2();
        }
    }

    // ==================================================================
    // Helper
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
    // Loop Data
    // ==================================================================
    std::cout << "=== Processing Data ===" << std::endl;
    {
        TFile* f = TFile::Open(dataFile, "READ");
        TTree* t = (TTree*)f->Get(kTreeName);

        std::vector<double>* cellE = nullptr;
        std::vector<double>* cellEta = nullptr;
        std::vector<double>* cellPhi = nullptr;
        Int_t cellSize = 0;
        Float_t eta2 = 0, photonPt = 0, etamax2 = 0;
        Double_t mll = 0, mllg = 0, dR1 = 0, dR2 = 0, dRll = 0;
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
        enableBranch(t, kDRLepton1Branch, &dR1);
        enableBranch(t, kDRLepton2Branch, &dR2);
        enableBranch(t, kDRllBranch, &dRll);
        enableBranch(t, kRetaUnfBranch, &br_reta);
        enableBranch(t, kRphiUnfBranch, &br_rphi);
        enableBranch(t, kWeta2UnfBranch, &br_weta2);

        Long64_t nEntries = t->GetEntries();
        if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
        Long64_t nUsed = 0;

        for (Long64_t i = 0; i < nEntries; ++i) {
            t->GetEntry(i);
            if (i > 0 && i % 2000000 == 0)
                std::cout << "  " << i << " / " << nEntries << std::endl;

            if (!passSelection(photonPt, eta2, mll, mllg,
                               dR1, dR2, dRll, cellSize, true, false))
                continue;
            if (!cellE || (int)cellE->size() != kClusterSize) continue;
            if (!cellEta || !cellPhi) continue;

            int etaBin = findEtaBin(std::fabs(eta2));
            if (etaBin < 0) continue;

            { std::vector<double> cl(*cellE); clampCellEnergies(cl);
              if (!isHealthyCluster(cl)) continue; }

            // Branch values
            h_reta[0][etaBin]->Fill(br_reta);
            h_rphi[0][etaBin]->Fill(br_rphi);
            h_weta2[0][etaBin]->Fill(br_weta2);

            // Cell-computed (position-based, no clamping)
            h_reta[1][etaBin]->Fill(calcRetaPos(*cellE, *cellEta, *cellPhi));
            h_rphi[1][etaBin]->Fill(calcRphiPos(*cellE, *cellEta, *cellPhi));
            h_weta2[1][etaBin]->Fill(calcWeta2Pos(*cellE, *cellEta, *cellPhi, eta2, etamax2));
            nUsed++;
        }
        std::cout << "  Data selected: " << nUsed << std::endl;
        f->Close(); delete f;
    }

    // ==================================================================
    // Loop MC
    // ==================================================================
    std::cout << "=== Processing MC ===" << std::endl;
    {
        TFile* f = TFile::Open(mcFile, "READ");
        TTree* t = (TTree*)f->Get(kTreeName);

        std::vector<double>* cellE = nullptr;
        std::vector<double>* cellEta = nullptr;
        std::vector<double>* cellPhi = nullptr;
        Int_t cellSize = 0;
        Float_t eta2 = 0, photonPt = 0, etamax2 = 0;
        Double_t mll = 0, mllg = 0, dR1 = 0, dR2 = 0, dRll = 0;
        Bool_t truthMatch = false;
        Float_t mcwgt = 1, muwgt = 1, bswgt = 1;
        Double_t lep1SF = 1, lep2SF = 1;
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
        enableBranch(t, kDRLepton1Branch, &dR1);
        enableBranch(t, kDRLepton2Branch, &dR2);
        enableBranch(t, kDRllBranch, &dRll);
        enableBranch(t, kTruthMatchBranch, &truthMatch);
        enableBranch(t, kMCWeightBranch, &mcwgt);
        enableBranch(t, kPUWeightBranch, &muwgt);
        enableBranch(t, kBSWeightBranch, &bswgt);
        enableBranch(t, kLepton1SFBranch, &lep1SF);
        enableBranch(t, kLepton2SFBranch, &lep2SF);
        enableBranch(t, kRetaUnfBranch, &br_reta);
        enableBranch(t, kRphiUnfBranch, &br_rphi);
        enableBranch(t, kWeta2UnfBranch, &br_weta2);

        Long64_t nEntries = t->GetEntries();
        if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
        Long64_t nUsed = 0;

        for (Long64_t i = 0; i < nEntries; ++i) {
            t->GetEntry(i);
            if (i > 0 && i % 2000000 == 0)
                std::cout << "  " << i << " / " << nEntries << std::endl;

            if (!passSelection(photonPt, eta2, mll, mllg,
                               dR1, dR2, dRll, cellSize, truthMatch, true))
                continue;
            if (!cellE || (int)cellE->size() != kClusterSize) continue;
            if (!cellEta || !cellPhi) continue;

            int etaBin = findEtaBin(std::fabs(eta2));
            if (etaBin < 0) continue;

            { std::vector<double> cl(*cellE); clampCellEnergies(cl);
              if (!isHealthyCluster(cl)) continue; }

            double w = mcwgt * muwgt * bswgt * lep1SF * lep2SF;

            // Branch values (unfudged)
            h_reta[2][etaBin]->Fill(br_reta, w);
            h_rphi[2][etaBin]->Fill(br_rphi, w);
            h_weta2[2][etaBin]->Fill(br_weta2, w);

            // Cell-computed (position-based, no clamping)
            h_reta[3][etaBin]->Fill(calcRetaPos(*cellE, *cellEta, *cellPhi), w);
            h_rphi[3][etaBin]->Fill(calcRphiPos(*cellE, *cellEta, *cellPhi), w);
            h_weta2[3][etaBin]->Fill(calcWeta2Pos(*cellE, *cellEta, *cellPhi, eta2, etamax2), w);
            nUsed++;
        }
        std::cout << "  MC selected: " << nUsed << std::endl;
        f->Close(); delete f;
    }

    // ==================================================================
    // Normalise
    // ==================================================================
    for (int s = 0; s < 4; ++s)
        for (int n = 0; n < kNEtaBins; ++n) {
            double I;
            I = h_reta[s][n]->Integral();  if (I > 0) h_reta[s][n]->Scale(1./I);
            I = h_rphi[s][n]->Integral();  if (I > 0) h_rphi[s][n]->Scale(1./I);
            I = h_weta2[s][n]->Integral(); if (I > 0) h_weta2[s][n]->Scale(1./I);
        }

    // ==================================================================
    // Plot: one page per eta bin, one PDF per variable
    // ==================================================================
    struct VarDef { const char* name; const char* title; const char* pdf; TH1D* (*arr)[kNEtaBins]; };

    // Can't take address of 2D array element portably, so use a lambda
    auto drawVar = [&](TH1D* hists[4][kNEtaBins], const char* varTitle, const char* pdfName) {
        TCanvas* c = new TCanvas("c_diag", "", 800, 700);
        TString pdfPath = dir + pdfName;

        for (int n = 0; n < kNEtaBins; ++n) {
            if (hists[0][n]->GetEntries() < 50 && hists[2][n]->GetEntries() < 50) continue;

            c->Clear();

            // Upper pad: overlay
            TPad* pad1 = new TPad("p1", "", 0, 0.30, 1, 1.0);
            pad1->SetBottomMargin(0.02);
            pad1->SetLeftMargin(0.14);
            pad1->Draw(); pad1->cd();

            // Style
            hists[0][n]->SetMarkerStyle(20); hists[0][n]->SetMarkerSize(0.8);
            hists[0][n]->SetMarkerColor(kBlack); hists[0][n]->SetLineColor(kBlack);
            hists[0][n]->SetLineWidth(1);

            hists[1][n]->SetLineColor(kBlue);
            hists[1][n]->SetLineWidth(2); hists[1][n]->SetLineStyle(1);

            hists[2][n]->SetMarkerStyle(24); hists[2][n]->SetMarkerSize(0.8);
            hists[2][n]->SetMarkerColor(kRed); hists[2][n]->SetLineColor(kRed);
            hists[2][n]->SetLineWidth(1);

            hists[3][n]->SetLineColor(kGreen+2);
            hists[3][n]->SetLineWidth(2); hists[3][n]->SetLineStyle(1);

            double ymax = 1.6 * std::max({hists[0][n]->GetMaximum(),
                                           hists[1][n]->GetMaximum(),
                                           hists[2][n]->GetMaximum(),
                                           hists[3][n]->GetMaximum()});

            hists[0][n]->SetMaximum(ymax); hists[0][n]->SetMinimum(0);
            hists[0][n]->GetXaxis()->SetLabelSize(0);
            hists[0][n]->GetYaxis()->SetTitle("Normalised");
            hists[0][n]->GetYaxis()->SetTitleSize(0.06);
            hists[0][n]->GetYaxis()->SetLabelSize(0.05);
            hists[0][n]->GetYaxis()->SetTitleOffset(0.95);
            hists[0][n]->SetTitle(""); hists[0][n]->SetStats(0);

            hists[0][n]->Draw("P");
            hists[1][n]->Draw("HIST SAME");
            hists[2][n]->Draw("P SAME");
            hists[3][n]->Draw("HIST SAME");

            TLegend* leg = new TLegend(0.48, 0.62, 0.89, 0.89);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
            leg->AddEntry(hists[0][n], "Data (branch)", "lp");
            leg->AddEntry(hists[1][n], "Data (cell-computed)", "l");
            leg->AddEntry(hists[2][n], "MC unfudged (branch)", "lp");
            leg->AddEntry(hists[3][n], "MC (cell-computed)", "l");
            leg->Draw();

            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.050); lat.SetTextFont(72);
            lat.DrawLatex(0.18, 0.86, "ATLAS");
            lat.SetTextFont(42);
            lat.DrawLatex(0.30, 0.86, "Internal");
            lat.SetTextSize(0.036);
            lat.DrawLatex(0.18, 0.80,
                Form("%s, |#eta| #in [%.2f, %.2f)",
                     varTitle, kEtaLimits[n], kEtaLimits[n+1]));

            // Lower pad: ratio (cell-computed / branch) for data and MC
            c->cd();
            TPad* pad2 = new TPad("p2", "", 0, 0.0, 1, 0.30);
            pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.35);
            pad2->SetLeftMargin(0.14);
            pad2->Draw(); pad2->cd();

            TH1D* r_data = (TH1D*)hists[1][n]->Clone(Form("rd_%d", n));
            TH1D* r_mc   = (TH1D*)hists[3][n]->Clone(Form("rm_%d", n));
            r_data->Divide(hists[0][n]);
            r_mc->Divide(hists[2][n]);

            r_data->SetMinimum(0.6); r_data->SetMaximum(1.4);
            r_data->GetYaxis()->SetTitle("Cell / Branch");
            r_data->GetYaxis()->SetTitleSize(0.12);
            r_data->GetYaxis()->SetTitleOffset(0.45);
            r_data->GetYaxis()->SetLabelSize(0.10);
            r_data->GetYaxis()->SetNdivisions(505);
            r_data->GetXaxis()->SetTitle(varTitle);
            r_data->GetXaxis()->SetTitleSize(0.14);
            r_data->GetXaxis()->SetLabelSize(0.10);
            r_data->SetStats(0);

            r_data->Draw("HIST");
            r_mc->Draw("HIST SAME");

            TLine* line = new TLine(r_data->GetXaxis()->GetXmin(), 1.0,
                                    r_data->GetXaxis()->GetXmax(), 1.0);
            line->SetLineColor(kGray+2); line->SetLineStyle(2); line->Draw();

            c->cd();
            if (n == 0 || (n == 1 && hists[0][0]->GetEntries() < 50))
                c->Print(pdfPath + "(");
            else
                c->Print(pdfPath);
        }
        // Close the PDF
        c->Clear();
        c->Print(pdfPath + ")");
        std::cout << "  Written: " << pdfPath << std::endl;
        delete c;
    };

    std::cout << "\n=== Creating diagnostic plots ===" << std::endl;
    drawVar(h_reta,  "R_{#eta}",    "diag_reta.pdf");
    drawVar(h_rphi,  "R_{#phi}",    "diag_rphi.pdf");
    drawVar(h_weta2, "w_{#eta 2}",  "diag_weta2.pdf");

    std::cout << "\n=== Done ===" << std::endl;
    return 0;
}
