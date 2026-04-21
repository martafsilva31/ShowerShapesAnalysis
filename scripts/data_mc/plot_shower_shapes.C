///////////////////////////////////////////////////////////////////////////////
// plot_shower_shapes.C
//
// Creates comparison plots with ratio panels from histograms produced by
// fill_histograms.C.  Generates PDFs for cell-computed shower shapes:
//
//   rew_<var>.pdf  (2-4 curves, per-eta):
//     Data (cell-computed), MC (cell-computed), M1 (flat shift), M2 (Francisco)
//
// When binning="eta_pt", also produces per-pT-bin PDFs:
//   rew_<var>_pt<PP>.pdf  (one PDF per pT bin, 14 eta-bin pages each)
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
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace config;

// ======================================================================
// Generic panel: overlay N histograms with ratio to first one
// histos[0] is drawn as data points, rest as lines.
// ======================================================================
void drawPanel(std::vector<TH1D*>& histos,
               std::vector<TString>& labels,
               std::vector<int>& colors,
               std::vector<int>& styles,
               const char* varTitle,
               int etaBin,
               const char* chLabel,
               const char* scenLabel,
               TCanvas* c,
               const char* overrideEtaLabel = nullptr,
               const char* ptLabel = nullptr,
               double legX1 = 0.59, double legY1 = 0.65,
               double legX2 = 0.96, double legY2 = 0.92) {

    c->Clear();

    if (histos.empty()) return;
    TH1D* hRef = histos[0];  // reference for ratio

    // --- Upper pad ---
    TPad* pad1 = new TPad("pad1", "", 0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->Draw();
    pad1->cd();

    // Style: first histogram as data markers, rest as lines
    hRef->SetMarkerStyle(20);
    hRef->SetMarkerSize(0.8);
    hRef->SetMarkerColor(colors[0]);
    hRef->SetLineColor(colors[0]);
    hRef->SetLineWidth(1);

    for (size_t i = 1; i < histos.size(); ++i) {
        histos[i]->SetLineColor(colors[i]);
        histos[i]->SetLineWidth(2);
        histos[i]->SetLineStyle(styles[i]);
    }

    // Y-axis range
    double ymax = 0;
    for (auto& h : histos) ymax = std::max(ymax, h->GetMaximum());
    ymax *= 1.6;

    // Auto-zoom X: find first/last bins with content > 0.5% of peak.
    double peak = 0;
    for (auto& h : histos) peak = std::max(peak, h->GetMaximum());
    double threshold = 0.005 * peak;
    int axFirst = hRef->GetXaxis()->GetFirst();
    int axLast  = hRef->GetXaxis()->GetLast();
    int firstBin = axFirst, lastBin = axLast;
    for (int b = axFirst; b <= axLast; ++b) {
        bool hasContent = false;
        for (auto& h : histos)
            if (h->GetBinContent(b) > threshold) { hasContent = true; break; }
        if (hasContent) { firstBin = b; break; }
    }
    for (int b = axLast; b >= axFirst; --b) {
        bool hasContent = false;
        for (auto& h : histos)
            if (h->GetBinContent(b) > threshold) { hasContent = true; break; }
        if (hasContent) { lastBin = b; break; }
    }
    if (firstBin > axFirst) firstBin--;
    if (lastBin < axLast)   lastBin++;
    double xlo = hRef->GetXaxis()->GetBinLowEdge(firstBin);
    double xhi = hRef->GetXaxis()->GetBinUpEdge(lastBin);
    hRef->GetXaxis()->SetRangeUser(xlo, xhi);

    hRef->SetMaximum(ymax);
    hRef->SetMinimum(0);
    hRef->GetXaxis()->SetLabelSize(0);
    hRef->GetYaxis()->SetTitle("Normalised");
    hRef->GetYaxis()->SetTitleSize(0.06);
    hRef->GetYaxis()->SetLabelSize(0.05);
    hRef->GetYaxis()->SetTitleOffset(0.95);
    hRef->SetTitle("");
    hRef->SetStats(0);

    hRef->Draw("E P X0");
    for (size_t i = 1; i < histos.size(); ++i) {
        TH1D* hBand = (TH1D*)histos[i]->Clone(
            Form("%s_band_%d", histos[i]->GetName(), etaBin));
        hBand->SetFillColorAlpha(colors[i], 0.15);
        hBand->SetLineWidth(0);
        hBand->Draw("E2 SAME");
        histos[i]->Draw("HIST SAME");
    }
    hRef->Draw("E P SAME X0");

    // Legend
    TLegend* leg = new TLegend(legX1, legY1, legX2, legY2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    for (size_t i = 0; i < histos.size(); ++i)
        leg->AddEntry(histos[i], labels[i], (i == 0) ? "ep" : "lf");
    leg->Draw();

    // ATLAS label
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050);
    lat.SetTextFont(72);
    lat.DrawLatex(0.18, 0.85, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(0.285, 0.85, "Work in Progress");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.18, 0.79, "#sqrt{s} = 13.6 TeV");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.18, 0.73, chLabel);
    lat.SetTextSize(0.033);
    lat.DrawLatex(0.18, 0.67, scenLabel);
    lat.SetTextSize(0.033);
    if (overrideEtaLabel)
        lat.DrawLatex(0.18, 0.61, overrideEtaLabel);
    else
        lat.DrawLatex(0.18, 0.61,
            Form("%.2f < |#eta| < %.2f",
                 kEtaLimits[etaBin], kEtaLimits[etaBin + 1]));

    // pT label (drawn below eta label if provided)
    if (ptLabel) {
        lat.SetTextSize(0.033);
        lat.DrawLatex(0.18, 0.55, ptLabel);
    }

    // --- Lower pad: ratio to data ---
    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    bool firstRatio = true;
    for (size_t i = 1; i < histos.size(); ++i) {
        TH1D* ratio = (TH1D*)histos[i]->Clone(
            Form("ratio_%zu_%d", i, etaBin));
        ratio->Divide(hRef);

        ratio->SetMinimum(0.50);
        ratio->SetMaximum(1.50);
        ratio->GetXaxis()->SetRangeUser(xlo, xhi);
        ratio->GetYaxis()->SetTitle("MC/Data");
        ratio->GetYaxis()->SetTitleSize(0.12);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitle(varTitle);
        ratio->GetXaxis()->SetTitleSize(0.14);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->SetTitle("");
        ratio->SetStats(0);

        TH1D* rBand = (TH1D*)ratio->Clone(
            Form("ratio_band_%zu_%d", i, etaBin));
        rBand->SetFillColorAlpha(colors[i], 0.15);
        rBand->SetLineWidth(0);
        if (firstRatio) {
            rBand->Draw("E2");
            ratio->Draw("HIST SAME");
        } else {
            rBand->Draw("E2 SAME");
            ratio->Draw("HIST SAME");
        }
        firstRatio = false;
    }

    TLine* line = new TLine(xlo, 1.0, xhi, 1.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
}


// ======================================================================
// Axis-zoom helper
// ======================================================================
void computeDataRange(TH1D* href, double& xLo, double& xHi,
                      double pctLo = 0.005, double pctHi = 0.995) {
    if (!href) { xLo = 0; xHi = 1; return; }
    TAxis* ax = href->GetXaxis();
    if (href->GetEntries() < 10) {
        xLo = ax->GetXmin(); xHi = ax->GetXmax(); return;
    }
    double p[2] = {pctLo, pctHi}, q[2];
    href->GetQuantiles(2, q, p);
    int bLo = std::max(1, ax->FindFixBin(q[0]));
    int bHi = std::min(ax->GetNbins(), ax->FindFixBin(q[1]));
    xLo = ax->GetBinLowEdge(bLo);
    xHi = ax->GetBinUpEdge(bHi);
}

// ======================================================================
// Main
// ======================================================================
int plot_shower_shapes(const char* channel   = "eegamma",
                       const char* scenario  = "unconverted",
                       const char* baseDir   = "../../output/Layer_2/eta_loose",
                       const char* binning   = "eta",
                       const char* isolation = "loose") {

    const char* chLabel   = channelLabel(channel);
    const char* scenLabel = scenarioLabel(scenario, isolation);

    TString binMode(binning);
    bool usePtBins = (binMode == "eta_pt");

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->SetBatch(true);

    TString outPath = Form("%s/%s/%s", baseDir, channel, scenario);
    TString histoFile = outPath + "/histograms.root";
    TString plotDir   = outPath + "/plots/";

    TFile* f = TFile::Open(histoFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << histoFile << std::endl;
        return 1;
    }

    const int nZoomBins = 50;
    int origNBins = 100;
    {
        TH1D* hProbe = (TH1D*)f->Get("h_reta_data_computed");
        if (hProbe) origNBins = hProbe->GetNbinsX();
    }
    const int rebinGroups = std::max(1, origNBins / nZoomBins);

    TString dir(plotDir);
    if (!dir.EndsWith("/")) dir += "/";
    gSystem->mkdir(dir, true);

    TCanvas* c = new TCanvas("c", "Data-MC Comparison", 800, 700);

    struct VarDef { const char* name; const char* title; };
    VarDef vars[] = {
        {"reta",  "R_{#eta}"},
        {"rphi",  "R_{#phi}"},
        {"weta2", "w_{#eta_{2}}"}
    };

    auto zoomNorm = [&](std::vector<TH1D*>& hv, TH1D* hRef) {
        for (auto*& h : hv) {
            h = (TH1D*)h->Clone(Form("%s_zm", h->GetName()));
            h->Rebin(rebinGroups);
        }
        double xLo, xHi;
        computeDataRange(hRef, xLo, xHi);
        for (auto* h : hv)
            h->GetXaxis()->SetRangeUser(xLo, xHi);
        for (auto* h : hv)
            if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
    };

    // ================================================================
    // SET B -- Cell-computed per-eta: data, mc, mc_M1, mc_M2
    // Produces: rew_reta.pdf, rew_rphi.pdf, rew_weta2.pdf
    // ================================================================
    for (auto& v : vars) {
        TString pdfPath = dir + Form("rew_%s.pdf", v.name);
        std::cout << "Set B " << v.title << ": " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (int n = 0; n < kNEtaBins; ++n) {
            TString s = Form("_eta%02d", n);
            TH1D* hData = (TH1D*)f->Get(Form("h_%s_data%s",    v.name, s.Data()));
            TH1D* hMC   = (TH1D*)f->Get(Form("h_%s_mc%s",      v.name, s.Data()));
            TH1D* hM1   = (TH1D*)f->Get(Form("h_%s_mc_M1%s",   v.name, s.Data()));
            TH1D* hM2   = (TH1D*)f->Get(Form("h_%s_mc_M2%s",   v.name, s.Data()));
            if (!hData || !hMC || !hM1 || !hM2) continue;
            if (hData->GetEntries() < 100) continue;

            std::vector<TH1D*> hv = {hData, hMC, hM1, hM2};
            zoomNorm(hv, hv[0]);

            std::vector<TString> labels = {
                "Data", "Original MC",
                "MC reweighted (shift only)", "MC reweighted (shift+stretch)"
            };
            std::vector<int> colors = {kBlack, kRed, kGreen + 2, kBlue};
            std::vector<int> styles = {1, 1, 1, 1};

            drawPanel(hv, labels, colors, styles,
                      v.title, n, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.56, 0.67, 0.93, 0.90);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // ================================================================
    // SET B per-(eta,pT) — only in eta_pt mode
    // Produces: rew_<var>_pt<PP>.pdf  (one PDF per pT bin, 14 pages)
    // ================================================================
    if (usePtBins) {
        for (int p = 0; p < kNPtBins; ++p) {
            for (auto& v : vars) {
                TString pdfPath = dir + Form("rew_%s_pt%02d.pdf", v.name, p);
                std::cout << "Set B per-pT " << v.title
                          << " pT bin " << p << ": " << pdfPath << std::endl;
                c->Print(pdfPath + "[");

                for (int n = 0; n < kNEtaBins; ++n) {
                    TString s = Form("_eta%02d_pt%02d", n, p);
                    TH1D* hData = (TH1D*)f->Get(Form("h_%s_data%s",  v.name, s.Data()));
                    TH1D* hMC   = (TH1D*)f->Get(Form("h_%s_mc%s",    v.name, s.Data()));
                    TH1D* hM1   = (TH1D*)f->Get(Form("h_%s_mc_M1%s", v.name, s.Data()));
                    TH1D* hM2   = (TH1D*)f->Get(Form("h_%s_mc_M2%s", v.name, s.Data()));
                    if (!hData || !hMC || !hM1 || !hM2) continue;
                    if (hData->GetEntries() < 50) continue;

                    std::vector<TH1D*> hv = {hData, hMC, hM1, hM2};
                    zoomNorm(hv, hv[0]);

                    std::vector<TString> labels = {
                        "Data", "Original MC",
                        "MC reweighted (shift only)", "MC reweighted (shift+stretch)"
                    };
                    std::vector<int> colors = {kBlack, kRed, kGreen + 2, kBlue};
                    std::vector<int> styles = {1, 1, 1, 1};

                    drawPanel(hv, labels, colors, styles,
                              v.title, n, chLabel, scenLabel, c,
                              nullptr, ptBinLabel(p), 0.56, 0.67, 0.93, 0.90);
                    c->Print(pdfPath);
                }
                c->Print(pdfPath + "]");
            }
        }
    }

    // ================================================================
    // SET B integrated -- Data, MC, M1, M2 (all eta combined)
    // Produces: rew_integrated.pdf  (3 pages: reta, rphi, weta2)
    // ================================================================
    {
        TString pdfPath = dir + "rew_integrated.pdf";
        std::cout << "Set B integrated: " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (auto& v : vars) {
            TH1D* hData = (TH1D*)f->Get(Form("h_%s_data_computed", v.name));
            TH1D* hMC   = (TH1D*)f->Get(Form("h_%s_mc_computed",  v.name));
            TH1D* hM1   = (TH1D*)f->Get(Form("h_%s_mc_M1",        v.name));
            TH1D* hM2   = (TH1D*)f->Get(Form("h_%s_mc_M2",        v.name));
            if (!hData || !hMC || !hM1 || !hM2) continue;

            std::vector<TH1D*> hv = {hData, hMC, hM1, hM2};
            zoomNorm(hv, hv[0]);

            std::vector<TString> labels = {
                "Data", "Original MC",
                "MC reweighted (shift only)", "MC reweighted (shift+stretch)"
            };
            std::vector<int> colors = {kBlack, kRed, kGreen + 2, kBlue};
            std::vector<int> styles = {1, 1, 1, 1};

            drawPanel(hv, labels, colors, styles,
                      v.title, 0, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.56, 0.67, 0.93, 0.90);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // ================================================================
    // Computed vs Stored -- integrated
    // Produces: computed_vs_stored.pdf
    // ================================================================
    {
        TString pdfPath = dir + "computed_vs_stored.pdf";
        std::cout << "Computed vs stored (integrated): " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (auto& v : vars) {
            TH1D* hComp = (TH1D*)f->Get(Form("h_%s_data_computed", v.name));
            TH1D* hStor = (TH1D*)f->Get(Form("h_%s_data_stored",   v.name));
            if (!hComp || !hStor) continue;

            std::vector<TH1D*> hv = {hComp, hStor};
            zoomNorm(hv, hv[0]);

            std::vector<TString> labels = {"Data (computed)", "Data (stored)"};
            std::vector<int> colors = {kBlack, kRed};
            std::vector<int> styles = {1, 1};

            drawPanel(hv, labels, colors, styles,
                      v.title, 0, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.65, 0.70, 0.96, 0.92);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // ================================================================
    // Computed vs Stored -- per-eta
    // Produces: computed_vs_stored_eta.pdf
    // ================================================================
    {
        TString pdfPath = dir + "computed_vs_stored_eta.pdf";
        std::cout << "Computed vs stored (per-eta): " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (auto& v : vars) {
            for (int n = 0; n < kNEtaBins; ++n) {
                TString s = Form("_eta%02d", n);
                TH1D* hComp = (TH1D*)f->Get(Form("h_%s_data%s",        v.name, s.Data()));
                TH1D* hStor = (TH1D*)f->Get(Form("h_%s_data_stored%s", v.name, s.Data()));
                if (!hComp || !hStor) continue;
                if (hComp->GetEntries() < 100) continue;

                std::vector<TH1D*> hv = {hComp, hStor};
                zoomNorm(hv, hv[0]);

                std::vector<TString> labels = {"Data (computed)", "Data (stored)"};
                std::vector<int> colors = {kBlack, kRed};
                std::vector<int> styles = {1, 1};

                drawPanel(hv, labels, colors, styles,
                          v.title, n, chLabel, scenLabel, c,
                          nullptr, nullptr, 0.65, 0.70, 0.96, 0.92);
                c->Print(pdfPath);
            }
        }
        c->Print(pdfPath + "]");
    }

    // ================================================================
    // Fudge factors -- integrated
    // Produces: fudge_factors.pdf
    // ================================================================
    {
        TString pdfPath = dir + "fudge_factors.pdf";
        std::cout << "Fudge factors (integrated): " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (auto& v : vars) {
            TH1D* hStor = (TH1D*)f->Get(Form("h_%s_data_stored",   v.name));
            TH1D* hFud  = (TH1D*)f->Get(Form("h_%s_mc_fudged",     v.name));
            TH1D* hUnf  = (TH1D*)f->Get(Form("h_%s_mc_unfudged",   v.name));
            if (!hStor || !hFud || !hUnf) continue;

            std::vector<TH1D*> hv = {hStor, hFud, hUnf};
            zoomNorm(hv, hv[0]);

            std::vector<TString> labels = {
                "Data", "MC (fudged)", "MC (unfudged)"
            };
            std::vector<int> colors = {kBlack, kBlue, kRed};
            std::vector<int> styles = {1, 1, 2};

            drawPanel(hv, labels, colors, styles,
                      v.title, 0, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.65, 0.62, 0.96, 0.88);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // ================================================================
    // Fudge factors -- per-eta
    // Produces: fudge_factors_eta.pdf
    // ================================================================
    {
        TString pdfPath = dir + "fudge_factors_eta.pdf";
        std::cout << "Fudge factors (per-eta): " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (auto& v : vars) {
            for (int n = 0; n < kNEtaBins; ++n) {
                TString s = Form("_eta%02d", n);
                TH1D* hStor = (TH1D*)f->Get(Form("h_%s_data_stored%s",  v.name, s.Data()));
                TH1D* hFud  = (TH1D*)f->Get(Form("h_%s_mc_fudged%s",   v.name, s.Data()));
                TH1D* hUnf  = (TH1D*)f->Get(Form("h_%s_mc_unfudged%s", v.name, s.Data()));
                if (!hStor || !hFud || !hUnf) continue;
                if (hStor->GetEntries() < 100) continue;

                std::vector<TH1D*> hv = {hStor, hFud, hUnf};
                zoomNorm(hv, hv[0]);

                std::vector<TString> labels = {
                    "Data", "MC (fudged)", "MC (unfudged)"
                };
                std::vector<int> colors = {kBlack, kBlue, kRed};
                std::vector<int> styles = {1, 1, 2};

                drawPanel(hv, labels, colors, styles,
                          v.title, n, chLabel, scenLabel, c,
                          nullptr, nullptr, 0.65, 0.62, 0.96, 0.88);
                c->Print(pdfPath);
            }
        }
        c->Print(pdfPath + "]");
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== All plots created in " << dir << " ===" << std::endl;
    return 0;
}
