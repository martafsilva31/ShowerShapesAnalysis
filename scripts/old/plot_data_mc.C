///////////////////////////////////////////////////////////////////////////////
// plot_data_mc.C
//
// Creates comparison plots with ratio panels from histograms produced by
// validate_data_mc.C.  Generates SEPARATE PDFs for the two histogram sets:
//
//   SET A — branch_<var>.pdf  (3 curves):
//     Data (branch), MC unfudged (branch), MC fudged (branch)
//     Purpose: Show effect of ATLAS fudge factors.
//
//   SET B — rew_<var>.pdf  (5 curves):
//     Data (cell), MC (cell), M1 (flat shift), M2 (profiled TProfile),
//     M3 (shift+stretch)
//     Purpose: Show effect of cell reweighting methods.
//
// Each PDF has one page per eta bin (14 pages).
//
// Usage:
//   root -l -b -q 'plot_data_mc.C("histos.root", "plots/")'
//   root -l -b -q 'plot_data_mc.C("histos.root", "plots/", "Z#rightarrowee#gamma")'
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

using namespace datamc;

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
               const char* channelLabel,
               TCanvas* c) {

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

    // Auto-zoom X: find first/last bins with content > 0.5% of peak (display only)
    double peak = 0;
    for (auto& h : histos) peak = std::max(peak, h->GetMaximum());
    double threshold = 0.005 * peak;
    int firstBin = 1, lastBin = hRef->GetNbinsX();
    for (int b = 1; b <= hRef->GetNbinsX(); ++b) {
        bool hasContent = false;
        for (auto& h : histos)
            if (h->GetBinContent(b) > threshold) { hasContent = true; break; }
        if (hasContent) { firstBin = b; break; }
    }
    for (int b = hRef->GetNbinsX(); b >= 1; --b) {
        bool hasContent = false;
        for (auto& h : histos)
            if (h->GetBinContent(b) > threshold) { hasContent = true; break; }
        if (hasContent) { lastBin = b; break; }
    }
    // Add 1-bin padding on each side
    if (firstBin > 1) firstBin--;
    if (lastBin < hRef->GetNbinsX()) lastBin++;
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

    hRef->Draw("P");
    for (size_t i = 1; i < histos.size(); ++i)
        histos[i]->Draw("HIST SAME");
    hRef->Draw("P SAME");  // redraw on top

    // Legend
    TLegend* leg = new TLegend(0.52, 0.60, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    for (size_t i = 0; i < histos.size(); ++i)
        leg->AddEntry(histos[i], labels[i], (i == 0) ? "lp" : "l");
    leg->Draw();

    // ATLAS label
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050);
    lat.SetTextFont(72);
    lat.DrawLatex(0.18, 0.86, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(0.30, 0.86, "Internal");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.18, 0.80, channelLabel);
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.18, 0.74,
        Form("%s, |#eta| #in [%.2f, %.2f)",
             varTitle, kEtaLimits[etaBin], kEtaLimits[etaBin + 1]));

    // --- Lower pad: ratio to data ---
    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->Draw();
    pad2->cd();

    bool firstRatio = true;
    for (size_t i = 1; i < histos.size(); ++i) {
        TH1D* ratio = (TH1D*)histos[i]->Clone(
            Form("ratio_%zu_%d", i, etaBin));
        ratio->Divide(hRef);

        ratio->SetMinimum(0.85);
        ratio->SetMaximum(1.15);
        ratio->GetXaxis()->SetRangeUser(xlo, xhi);
        ratio->GetYaxis()->SetTitle("MC / Data");
        ratio->GetYaxis()->SetTitleSize(0.12);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitle(varTitle);
        ratio->GetXaxis()->SetTitleSize(0.14);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->SetStats(0);

        ratio->Draw(firstRatio ? "HIST" : "HIST SAME");
        firstRatio = false;
    }

    // Reference line at 1 (use zoomed range)
    TLine* line = new TLine(xlo, 1.0, xhi, 1.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
}


// ======================================================================
// Main
// ======================================================================
int plot_data_mc(const char* histoFile,
                 const char* plotDir,
                 const char* channelLabel = "Z#rightarrowll#gamma") {

    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);

    TFile* f = TFile::Open(histoFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << histoFile << std::endl;
        return 1;
    }

    TString dir(plotDir);
    if (!dir.EndsWith("/")) dir += "/";
    gSystem->mkdir(dir, true);

    TCanvas* c = new TCanvas("c", "Data-MC Comparison", 800, 700);

    // Variable definitions
    struct VarDef {
        const char* name; const char* title;
    };
    VarDef vars[] = {
        {"reta",  "R_{#eta}"},
        {"rphi",  "R_{#phi}"},
        {"weta2", "w_{#eta 2}"}
    };

    // ================================================================
    // SET A — Branch-based: data_br, mc_unf, mc_fud
    // Produces: branch_reta.pdf, branch_rphi.pdf, branch_weta2.pdf
    // ================================================================
    for (auto& v : vars) {
        TString pdfPath = dir + Form("branch_%s.pdf", v.name);
        std::cout << "Set A " << v.title << ": " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (int n = 0; n < kNEtaBins; ++n) {
            TString s = Form("_%d", n);
            TH1D* hData = (TH1D*)f->Get(Form("%s_data_br%s",  v.name, s.Data()));
            TH1D* hUnf  = (TH1D*)f->Get(Form("%s_mc_unf%s",   v.name, s.Data()));
            TH1D* hFud  = (TH1D*)f->Get(Form("%s_mc_fud%s",   v.name, s.Data()));

            if (!hData || !hUnf || !hFud) continue;
            if (hData->GetEntries() < 100) continue;

            std::vector<TH1D*> histos = {hData, hUnf, hFud};
            std::vector<TString> labels = {"Data", "MC unfudged", "MC fudged"};
            std::vector<int> colors = {kBlack, kRed, kRed};
            std::vector<int> styles = {1, 1, 2};

            drawPanel(histos, labels, colors, styles,
                      v.title, n, channelLabel, c);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // ================================================================
    // SET B — Cell-computed: data_cell, mc_cell, m1, m2, m3
    // Produces: rew_reta.pdf, rew_rphi.pdf, rew_weta2.pdf
    // ================================================================
    for (auto& v : vars) {
        TString pdfPath = dir + Form("rew_%s.pdf", v.name);
        std::cout << "Set B " << v.title << ": " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (int n = 0; n < kNEtaBins; ++n) {
            TString s = Form("_%d", n);
            TH1D* hData = (TH1D*)f->Get(Form("%s_data_cell%s", v.name, s.Data()));
            TH1D* hMC   = (TH1D*)f->Get(Form("%s_mc_cell%s",   v.name, s.Data()));
            TH1D* hM1   = (TH1D*)f->Get(Form("%s_m1%s",        v.name, s.Data()));
            TH1D* hM2   = (TH1D*)f->Get(Form("%s_m2%s",        v.name, s.Data()));
            TH1D* hM3   = (TH1D*)f->Get(Form("%s_m3%s",        v.name, s.Data()));

            if (!hData || !hMC || !hM1 || !hM2) continue;
            if (hData->GetEntries() < 100) continue;

            std::vector<TH1D*> histos = {hData, hMC, hM1, hM2};
            std::vector<TString> labels = {
                "Data (cell)", "MC (cell)",
                "Rew. M1 (flat shift)", "Rew. M2 (profiled)"
            };
            std::vector<int> colors = {kBlack, kRed, kGreen + 2, kBlue};
            std::vector<int> styles = {1, 1, 1, 1};

            if (hM3) {
                histos.push_back(hM3);
                labels.push_back("Rew. M3 (shift+stretch)");
                colors.push_back(kOrange + 1);
                styles.push_back(1);
            }

            drawPanel(histos, labels, colors, styles,
                      v.title, n, channelLabel, c);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== All plots created ===" << std::endl;
    return 0;
}
