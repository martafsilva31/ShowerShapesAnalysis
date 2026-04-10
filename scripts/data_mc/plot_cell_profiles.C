///////////////////////////////////////////////////////////////////////////////
// plot_cell_profiles.C
//
// Six separate PDFs, one page per eta bin (14 bins), portrait canvas:
//   cell_data.pdf      — Data mean cell energy fractions
//   cell_mc.pdf        — MC mean cell energy fractions (original)
//   cell_mc_m1.pdf     — MC after M1 (shift) reweighting
//   cell_mc_m2.pdf     — MC after M2 (shift+stretch) reweighting
//   cell_shift.pdf     — Per-cell shift correction (M1 delta_k)
//   cell_stretch.pdf   — Per-cell stretch correction (M2 sigma_d/sigma_MC)
//
// Header layout matches Francisco's style (3-column, 2 rows above the plot):
//   ATLAS Work in Progress  |  channel label  |  map type
//   sqrt(s) = 13.6 TeV      |  eta range / bin number  |  scenario
//
// Usage:
//   root -l -b -q 'plot_cell_profiles.C("eegamma", "baseline")'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <cmath>

using namespace config;

// ======================================================================
// Reshape a 1D correction histogram (kClusterSize bins) into a 7x11 TH2D.
// Cell index k -> eta = k/kPhiSize, phi = k%kPhiSize
// ======================================================================
TH2D* corrTo2D(TH1D* h1, const char* name) {
    TH2D* h2 = new TH2D(name, "",
                        kEtaSize, 0.5, kEtaSize + 0.5,
                        kPhiSize, 0.5, kPhiSize + 0.5);
    for (int k = 0; k < kClusterSize; ++k)
        h2->SetBinContent(k / kPhiSize + 1, k % kPhiSize + 1,
                          h1->GetBinContent(k + 1));
    return h2;
}


// ======================================================================
// Draw one cell map on canvas c with Francisco-style 3-column header.
// Single-canvas layout (no sub-pads): top margin holds labels,
// margins set directly on c. Canvas should be 800x600.
//
// Header rows (NDC, within top margin):
//   Row 1: ATLAS WIP (left) | channel (centre) | map type (right)
//   Row 2: sqrt(s)  (left)  | eta bin  (centre) | scenario  (right)
// ======================================================================
void drawCellMap(TH2D* hIn,
                 const char* typeLabel,   // "Data", "Shift correction", etc.
                 const char* chLabel,     // "Z#rightarrowee#gamma"
                 const char* etaLabel,    // "Bin 3:  0.60 < |#eta| < 0.80"
                 const char* scenLabel,   // scenario description
                 TCanvas* c) {
    c->cd();
    c->Clear();

    // Single-canvas margins (match Francisco's style)
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.18);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.16);

    // Style and draw the histogram directly (no clone needed)
    hIn->SetStats(0);
    hIn->SetTitle("");

    hIn->GetXaxis()->SetTitle("Cell (#eta direction)");
    hIn->GetYaxis()->SetTitle("Cell (#phi direction)");
    hIn->GetXaxis()->SetTitleSize(0.050);
    hIn->GetYaxis()->SetTitleSize(0.050);
    hIn->GetXaxis()->SetLabelSize(0.050);
    hIn->GetYaxis()->SetLabelSize(0.050);
    hIn->GetXaxis()->SetTitleOffset(1.0);
    hIn->GetYaxis()->SetTitleOffset(1.0);
    hIn->GetYaxis()->SetLabelOffset(0.01);
    hIn->GetZaxis()->SetLabelSize(0.030);
    hIn->GetZaxis()->SetTitleSize(0.045);

    // Per-cell bin labels (1..7 on X, 1..11 on Y) — matches Francisco's style
    for (int i = 1; i <= kEtaSize; ++i) hIn->GetXaxis()->SetBinLabel(i, Form("%d", i));
    for (int i = 1; i <= kPhiSize; ++i) hIn->GetYaxis()->SetBinLabel(i, Form("%d", i));
    // Ticks at every bin boundary (0,1,...,7 and 0,1,...,11)
    hIn->GetYaxis()->SetNdivisions(kPhiSize, false);
    hIn->GetXaxis()->SetNdivisions(kEtaSize, false);

    hIn->Draw("COLZ");

    // Cell value annotations drawn at bin centres (using GetBinCenter to
    // handle any axis range correctly: [0,7], [0.5,7.5], etc.)
    TLatex txt;
    txt.SetTextAlign(22);
    txt.SetTextFont(42);
    txt.SetTextSize(0.022);
    txt.SetTextColor(kBlack);

    for (int ei = 1; ei <= kEtaSize; ++ei) {
        double xc = hIn->GetXaxis()->GetBinCenter(ei);
        for (int pi = 1; pi <= kPhiSize; ++pi) {
            double yc  = hIn->GetYaxis()->GetBinCenter(pi);
            double val  = hIn->GetBinContent(ei, pi);
            double aval = std::fabs(val);
            TString label;
            if      (aval < 1e-9)  label = "0";
            else if (aval < 0.01)  label = Form("%.1e", val);
            else                   label = Form("%.3f", val);
            txt.DrawLatex(xc, yc, label);
        }
    }

    // ---- 3-column header in top margin (NDC) ----
    // With topMargin=0.16, top of plot area is at y_ndc=0.84.
    // Row 1 baseline y=0.955, Row 2 baseline y=0.895 — both within margin.
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);

    // Row 1
    lat.SetTextSize(0.038);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.03, 0.955, "#bf{#it{ATLAS}} Work in Progress");
    lat.DrawLatex(0.44, 0.955, chLabel);
    lat.SetTextAlign(31);                              // right-justified
    lat.DrawLatex(0.97, 0.955, typeLabel);

    // Row 2
    lat.SetTextSize(0.033);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.03, 0.895, "#sqrt{s} = 13.6 TeV");
    lat.DrawLatex(0.44, 0.895, etaLabel);
    lat.SetTextSize(0.028);
    lat.SetTextAlign(31);                              // right-justified
    lat.DrawLatex(0.97, 0.895, scenLabel);
}


// ======================================================================
// Main
// ======================================================================
int plot_cell_profiles(const char* channel  = "eegamma",
                       const char* scenario = "baseline",
                       const char* baseDir  =
                           "../../output/cell_energy_reweighting_Francisco_method/data24") {

    const char* chLabel   = channelLabel(channel);
    const char* scenLabel = scenarioLabel(scenario);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(88);  // kLightTerrain — Francisco's palette

    gROOT->SetBatch(true);


    TString outPath = Form("%s/%s/%s", baseDir, channel, scenario);
    TString corrFile = outPath + "/histograms.root";
    TString plotDir  = outPath + "/plots/";

    TFile* f = TFile::Open(corrFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << corrFile << std::endl;
        return 1;
    }

    TString dir(plotDir);
    if (!dir.EndsWith("/")) dir += "/";
    gSystem->mkdir(dir, true);

    // Canvas: 800x600 landscape (matches Francisco's style).
    // Plot area ~560x432px for 7x11 cells → ~80x39 px/cell.
    TCanvas* c = new TCanvas("c", "Cell Profiles", 800, 600);

    // Build eta label with bin number
    auto etaLbl = [](int n) -> TString {
        return Form("Bin %d:  %.2f < |#eta| < %.2f",
                    n, kEtaLimits[n], kEtaLimits[n + 1]);
    };

    // Check if a TH2D has any non-zero content
    auto nonEmpty = [](TH2D* h) -> bool {
        if (!h) return false;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            for (int j = 1; j <= h->GetNbinsY(); ++j)
                if (h->GetBinContent(i, j) != 0) return true;
        return false;
    };

    // ================================================================
    // 1. cell_data.pdf — Data mean cell energy fractions
    // ================================================================
    {
        TString pdf = dir + "cell_data.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH2D* h = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_data_eta%02d", n));
            if (!nonEmpty(h)) continue;
            drawCellMap(h, "Data", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
        }
        c->Print(pdf + "]");
    }

    // ================================================================
    // 2. cell_mc.pdf — MC mean cell energy fractions (original)
    // ================================================================
    {
        TString pdf = dir + "cell_mc.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH2D* h = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc_eta%02d", n));
            if (!nonEmpty(h)) continue;
            drawCellMap(h, "Original MC", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
        }
        c->Print(pdf + "]");
    }

    // ================================================================
    // 3. cell_mc_m1.pdf — MC after M1 (shift) reweighting
    //    M1-corrected = MC + delta  (equals data by construction)
    // ================================================================
    {
        TString pdf = dir + "cell_mc_m1.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH2D* hMC    = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc_eta%02d",    n));
            TH2D* hDelta = (TH2D*)f->Get(Form("cell_profiles/h_frac_delta_eta%02d",      n));
            if (!nonEmpty(hMC) || !hDelta) continue;
            TH2D* hM1 = (TH2D*)hMC->Clone(Form("cellsM1_%d", n));
            hM1->Add(hDelta);
            drawCellMap(hM1, "MC after M1 reweighting", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
            delete hM1;
        }
        c->Print(pdf + "]");
    }

    // ================================================================
    // 4. cell_mc_m2.pdf — MC after M2 (shift+stretch) reweighting
    // ================================================================
    {
        TString pdf = dir + "cell_mc_m2.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH2D* h = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc_M2_eta%02d", n));
            if (!nonEmpty(h)) continue;
            drawCellMap(h, "MC after M2 reweighting", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
        }
        c->Print(pdf + "]");
    }

    // ================================================================
    // 5. cell_shift.pdf — M1 shift correction (delta_k per cell)
    // ================================================================
    {
        TString pdf = dir + "cell_shift.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH1D* h1 = (TH1D*)f->Get(Form("corrections/h_delta_eta%02d", n));
            if (!h1) continue;
            bool ok = false;
            for (int k = 1; k <= kClusterSize; ++k)
                if (h1->GetBinContent(k) != 0) { ok = true; break; }
            if (!ok) continue;
            TH2D* h2 = corrTo2D(h1, Form("shift2d_%d", n));
            drawCellMap(h2, "Shift correction", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
            delete h2;
        }
        c->Print(pdf + "]");
    }

    // ================================================================
    // 6. cell_stretch.pdf — M2 stretch correction (sigma_d/sigma_MC)
    // ================================================================
    {
        TString pdf = dir + "cell_stretch.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH1D* h1 = (TH1D*)f->Get(Form("corrections/h_stretch_eta%02d", n));
            if (!h1) continue;
            bool ok = false;
            for (int k = 1; k <= kClusterSize; ++k)
                if (h1->GetBinContent(k) != 0) { ok = true; break; }
            if (!ok) continue;
            TH2D* h2 = corrTo2D(h1, Form("stretch2d_%d", n));
            drawCellMap(h2, "Stretch correction", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
            delete h2;
        }
        c->Print(pdf + "]");
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== All cell profile plots created in " << dir << " ===" << std::endl;
    return 0;
}
