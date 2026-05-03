///////////////////////////////////////////////////////////////////////////////
// plot_cell_profiles.C
//
// Six separate PDFs, one page per eta bin (14 bins), portrait canvas:
//   cell_data.pdf      -- Data mean cell energy fractions
//   cell_mc.pdf        -- MC mean cell energy fractions (original)
//   cell_mc_m1.pdf     -- MC after M1 (shift) reweighting
//   cell_mc_m2.pdf     -- MC after M2 (shift+stretch) reweighting
//   cell_shift.pdf     -- Per-cell shift correction (M1 delta_k)
//   cell_stretch.pdf   -- Per-cell stretch correction (M2 sigma_d/sigma_MC)
//
// When binning="eta_pt", also produces per-pT-bin PDFs:
//   cell_<type>_pt<PP>.pdf  (one PDF per pT bin, 14 eta-bin pages each)
//
// Parameters:
//   channel   "eegamma", "mumugamma", or "llgamma"
//   scenario  conversion scenario
//   baseDir   output base directory
//   binning   "eta" or "eta_pt"
//   isolation "loose" or "tight"
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
// Draw one cell map with 3-column header.
// Optional ptLabel drawn in row 3 (below eta label).
// ======================================================================
void drawCellMap(TH2D* hIn,
                 const char* typeLabel,
                 const char* chLabel,
                 const char* etaLabel,
                 const char* scenLabel,
                 TCanvas* c,
                 const char* ptLabel = nullptr) {
    c->cd();
    c->Clear();

    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.18);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(ptLabel ? 0.17 : 0.16);

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

    for (int i = 1; i <= kEtaSize; ++i) hIn->GetXaxis()->SetBinLabel(i, Form("%d", i));
    for (int i = 1; i <= kPhiSize; ++i) hIn->GetYaxis()->SetBinLabel(i, Form("%d", i));
    hIn->GetYaxis()->SetNdivisions(kPhiSize, false);
    hIn->GetXaxis()->SetNdivisions(kEtaSize, false);

    hIn->Draw("COLZ");

    // Cell value annotations
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
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);

    // Col 1 (left 0.03): ATLAS WiP / sqrt(s)
    // Col 2 (center, aligned with plot content): etaLabel / ptLabel
    //   left margin=0.12, right margin=0.18 → plot centre at NDC x≈0.47
    // Col 3 (right 0.97): chLabel / scenLabel
    double row1Y = ptLabel ? 0.955 : 0.955;
    double row2Y = ptLabel ? 0.900 : 0.895;

    // Row 1
    lat.SetTextSize(0.038);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.03, row1Y, "#bf{#it{ATLAS}} Work in Progress");
    lat.SetTextSize(0.033);
    lat.SetTextAlign(22);
    lat.DrawLatex(0.47, row1Y, etaLabel);
    lat.SetTextAlign(31);
    lat.DrawLatex(0.97, row1Y, chLabel);

    // Row 2
    lat.SetTextSize(0.033);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.03, row2Y, "#sqrt{s} = 13.6 TeV");
    if (ptLabel) {
        lat.SetTextAlign(22);
        lat.DrawLatex(0.47, row2Y, ptLabel);
    }
    lat.SetTextSize(0.028);
    lat.SetTextAlign(31);
    lat.DrawLatex(0.97, row2Y, scenLabel);
}


// ======================================================================
// Main
// ======================================================================
int plot_cell_profiles(const char* channel   = "eegamma",
                       const char* scenario  = "unconverted",
                       const char* baseDir   = "../../output/Layer_2/eta_loose",
                       const char* binning   = "eta",
                       const char* isolation = "loose") {

    const char* chLabel   = channelLabel(channel);
    const char* scenLabel = scenarioLabel(scenario, isolation);

    TString binMode(binning);
    bool usePtBins  = (binMode == "eta_pt");
    bool useMuBins  = (binMode == "eta_mu");
    bool useSecBins = usePtBins || useMuBins;
    int  nSecBins   = usePtBins ? kNPtBins : useMuBins ? kNMuBins : 0;
    TString secTag  = usePtBins ? "pt" : useMuBins ? "mu" : "";

    gStyle->SetOptStat(0);
    gStyle->SetPalette(88);
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

    TCanvas* c = new TCanvas("c", "Cell Profiles", 800, 600);

    auto etaLbl = [](int n) -> TString {
        return Form("%.2f < |#eta| < %.2f",
                    kEtaLimits[n], kEtaLimits[n + 1]);
    };

    auto nonEmpty = [](TH2D* h) -> bool {
        if (!h) return false;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            for (int j = 1; j <= h->GetNbinsY(); ++j)
                if (h->GetBinContent(i, j) != 0) return true;
        return false;
    };

    // Correction histogram suffix: eta-only or eta_sec
    auto corrSuf = [&](int e, int p = -1) -> TString {
        if (useSecBins && p >= 0)
            return Form("_eta%02d_%s%02d", e, secTag.Data(), p);
        return Form("_eta%02d", e);
    };

    // ================================================================
    // Per-eta cell profile PDFs (always produced)
    // ================================================================
    struct MapDef {
        const char* fileName;
        const char* typeLabel;
        const char* histPat;  // pattern in cell_profiles/
        bool isCorrHist;      // true = from corrections/ (1D -> corrTo2D)
        bool needsM1Sum;      // true = MC + delta
    };

    MapDef maps[] = {
        {"cell_data",    "Data",                      "h_frac_mean_data",   false, false},
        {"cell_mc",      "Original MC",               "h_frac_mean_mc",     false, false},
        {"cell_mc_m2",   "MC after M2 reweighting",   "h_frac_mean_mc_M2", false, false},
        {"cell_shift",   "Shift correction",          "h_delta",            true,  false},
        {"cell_stretch", "Stretch correction",         "h_stretch",         true,  false},
    };

    // ---- eta-only profiles ----
    for (auto& m : maps) {
        if (m.isCorrHist && useSecBins) continue;  // shift/stretch only meaningful in eta-only mode
        TString pdf = dir + Form("%s.pdf", m.fileName);
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TString hname;
            if (m.isCorrHist)
                hname = Form("corrections/%s%s", m.histPat, corrSuf(n).Data());
            else
                hname = Form("cell_profiles/%s%s", m.histPat, corrSuf(n).Data());

            if (m.isCorrHist) {
                TH1D* h1 = (TH1D*)f->Get(hname);
                if (!h1) continue;
                bool ok = false;
                for (int k = 1; k <= kClusterSize; ++k)
                    if (h1->GetBinContent(k) != 0) { ok = true; break; }
                if (!ok) continue;
                TH2D* h2 = corrTo2D(h1, Form("%s_2d_%d", m.fileName, n));
                drawCellMap(h2, m.typeLabel, chLabel, etaLbl(n), scenLabel, c);
                c->Print(pdf);
                delete h2;
            } else {
                TH2D* h = (TH2D*)f->Get(hname);
                if (!nonEmpty(h)) continue;
                drawCellMap(h, m.typeLabel, chLabel, etaLbl(n), scenLabel, c);
                c->Print(pdf);
            }
        }
        c->Print(pdf + "]");
    }

    // ---- M1 profile (MC + delta, computed on the fly) ----
    {
        TString pdf = dir + "cell_mc_m1.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH2D* hMC    = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc%s",
                                               corrSuf(n).Data()));
            TH2D* hDelta = (TH2D*)f->Get(Form("cell_profiles/h_frac_delta%s",
                                               corrSuf(n).Data()));
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
    // Per-sec cell profile PDFs (in eta_pt and eta_mu modes)
    // eta_pt produces: cell_<type>_pt<PP>.pdf
    // eta_mu produces: cell_<type>_mu<MM>.pdf
    // ================================================================
    if (useSecBins) {
        for (int p = 0; p < nSecBins; ++p) {
            const char* secLbl = usePtBins ? ptBinLabel(p) : muBinLabel(p);

            for (auto& m : maps) {
                TString pdf = dir + Form("%s_%s%02d.pdf", m.fileName, secTag.Data(), p);
                std::cout << "Creating: " << pdf << std::endl;
                c->Print(pdf + "[");
                for (int n = 0; n < kNEtaBins; ++n) {
                    TString hname;
                    TString suf = corrSuf(n, p);
                    if (m.isCorrHist)
                        hname = Form("corrections/%s%s", m.histPat, suf.Data());
                    else
                        hname = Form("cell_profiles/%s%s", m.histPat, suf.Data());

                    if (m.isCorrHist) {
                        TH1D* h1 = (TH1D*)f->Get(hname);
                        if (!h1) continue;
                        bool ok = false;
                        for (int k = 1; k <= kClusterSize; ++k)
                            if (h1->GetBinContent(k) != 0) { ok = true; break; }
                        if (!ok) continue;
                        TH2D* h2 = corrTo2D(h1, Form("%s_%s%d_2d_%d", m.fileName, secTag.Data(), p, n));
                        drawCellMap(h2, m.typeLabel, chLabel, etaLbl(n), scenLabel, c, secLbl);
                        c->Print(pdf);
                        delete h2;
                    } else {
                        TH2D* h = (TH2D*)f->Get(hname);
                        if (!nonEmpty(h)) continue;
                        drawCellMap(h, m.typeLabel, chLabel, etaLbl(n), scenLabel, c, secLbl);
                        c->Print(pdf);
                    }
                }
                c->Print(pdf + "]");
            }

            // M1 per-sec
            {
                TString pdf = dir + Form("cell_mc_m1_%s%02d.pdf", secTag.Data(), p);
                std::cout << "Creating: " << pdf << std::endl;
                c->Print(pdf + "[");
                for (int n = 0; n < kNEtaBins; ++n) {
                    TString suf = corrSuf(n, p);
                    TH2D* hMC    = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc%s", suf.Data()));
                    TH2D* hDelta = (TH2D*)f->Get(Form("cell_profiles/h_frac_delta%s", suf.Data()));
                    if (!nonEmpty(hMC) || !hDelta) continue;
                    TH2D* hM1 = (TH2D*)hMC->Clone(Form("cellsM1_%s%d_%d", secTag.Data(), p, n));
                    hM1->Add(hDelta);
                    drawCellMap(hM1, "MC after M1 reweighting", chLabel, etaLbl(n), scenLabel, c, secLbl);
                    c->Print(pdf);
                    delete hM1;
                }
                c->Print(pdf + "]");
            }
        }
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== All cell profile plots created in " << dir << " ===" << std::endl;
    return 0;
}
