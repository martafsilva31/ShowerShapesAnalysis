///////////////////////////////////////////////////////////////////////////////
// plot_stored_vs_computed.C
//
// Validates that cell-level shower shape computation matches the stored
// ntuple branch values.  Uses the same drawPanel layout as plot_data_mc.C
// (800x700 canvas, upper 70% overlay + lower 30% ratio, ATLAS labels).
//
// Produces:
//   stored_vs_computed_reta.pdf   — R_eta: branch vs cell-computed
//   stored_vs_computed_rphi.pdf   — R_phi: branch vs cell-computed
//   stored_vs_computed_weta2.pdf  — w_eta2: branch vs cell-computed
//
// Each PDF has one page per eta bin (14 pages).
// Per page: 4 curves — Data(branch), Data(cell), MC(branch), MC(cell)
// Ratio panel: cell/branch for data and MC separately.
//
// Usage:
//   root -l -b -q 'plot_stored_vs_computed.C("histos.root", "plots/")'
//   root -l -b -q 'plot_stored_vs_computed.C("histos.root", "plots/", "Z#rightarrowee#gamma")'
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
// drawPanel — overlay histograms with ratio to reference
// Same layout as plot_data_mc.C::drawPanel.
// histos[0] is data(branch), drawn as points.
// histos[1..] are lines.  Ratio to histos[0].
// ======================================================================
void drawPanel(std::vector<TH1D*>& histos,
               std::vector<TString>& labels,
               std::vector<int>& colors,
               std::vector<int>& styles,
               std::vector<int>& markers,
               const char* varTitle,
               int etaBin,
               const char* channelLabel,
               TCanvas* c) {

    c->Clear();
    if (histos.empty()) return;
    TH1D* hRef = histos[0];

    // --- Upper pad ---
    TPad* pad1 = new TPad("pad1", "", 0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->Draw();
    pad1->cd();

    // Styles
    for (size_t i = 0; i < histos.size(); ++i) {
        histos[i]->SetLineColor(colors[i]);
        histos[i]->SetLineWidth(2);
        histos[i]->SetLineStyle(styles[i]);
        if (markers[i] > 0) {
            histos[i]->SetMarkerStyle(markers[i]);
            histos[i]->SetMarkerSize(0.8);
            histos[i]->SetMarkerColor(colors[i]);
        }
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

    hRef->Draw(markers[0] > 0 ? "P" : "HIST");
    for (size_t i = 1; i < histos.size(); ++i)
        histos[i]->Draw(markers[i] > 0 ? "P SAME" : "HIST SAME");
    hRef->Draw(markers[0] > 0 ? "P SAME" : "HIST SAME");

    // Legend
    TLegend* leg = new TLegend(0.52, 0.60, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    for (size_t i = 0; i < histos.size(); ++i)
        leg->AddEntry(histos[i], labels[i], (markers[i] > 0) ? "lp" : "l");
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

    // --- Lower pad: ratio cell/branch ---
    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->Draw();
    pad2->cd();

    // Ratio: data_cell / data_branch  and  mc_cell / mc_branch
    // histos = [data_br, data_cell, mc_br, mc_cell]
    // Ratios: [1]/[0] and [3]/[2]
    TH1D* rData = (TH1D*)histos[1]->Clone(Form("rData_%d", etaBin));
    rData->Divide(histos[0]);
    rData->SetLineColor(colors[1]);
    rData->SetLineWidth(2);
    rData->SetLineStyle(styles[1]);
    rData->SetMarkerStyle(0);

    TH1D* rMC = (TH1D*)histos[3]->Clone(Form("rMC_%d", etaBin));
    rMC->Divide(histos[2]);
    rMC->SetLineColor(colors[3]);
    rMC->SetLineWidth(2);
    rMC->SetLineStyle(styles[3]);
    rMC->SetMarkerStyle(0);

    rData->SetMinimum(0.85);
    rData->SetMaximum(1.15);
    rData->GetXaxis()->SetRangeUser(xlo, xhi);
    rData->GetYaxis()->SetTitle("Cell / Branch");
    rData->GetYaxis()->SetTitleSize(0.12);
    rData->GetYaxis()->SetTitleOffset(0.45);
    rData->GetYaxis()->SetLabelSize(0.10);
    rData->GetYaxis()->SetNdivisions(505);
    rData->GetXaxis()->SetTitle(varTitle);
    rData->GetXaxis()->SetTitleSize(0.14);
    rData->GetXaxis()->SetLabelSize(0.10);
    rData->SetStats(0);

    rData->Draw("HIST");
    rMC->Draw("HIST SAME");

    // Reference line at 1 (use zoomed range)
    TLine* line = new TLine(xlo, 1.0, xhi, 1.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
}

// ======================================================================
// Main entry point
// ======================================================================
int plot_stored_vs_computed(const char* histoFile,
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

    TCanvas* c = new TCanvas("c", "Stored vs Computed", 800, 700);

    struct VarDef {
        const char* name; const char* title;
    };
    VarDef vars[] = {
        {"reta",  "R_{#eta}"},
        {"rphi",  "R_{#phi}"},
        {"weta2", "w_{#eta 2}"}
    };

    for (auto& v : vars) {
        TString pdfPath = dir + Form("stored_vs_computed_%s.pdf", v.name);
        std::cout << "Stored vs Computed " << v.title << ": " << pdfPath << std::endl;
        c->Print(pdfPath + "[");

        for (int n = 0; n < kNEtaBins; ++n) {
            TString s = Form("_%d", n);

            TH1D* hDataBr   = (TH1D*)f->Get(Form("%s_data_br%s",   v.name, s.Data()));
            TH1D* hDataCell = (TH1D*)f->Get(Form("%s_data_cell%s",  v.name, s.Data()));
            TH1D* hMCBr     = (TH1D*)f->Get(Form("%s_mc_unf%s",    v.name, s.Data()));
            TH1D* hMCCell   = (TH1D*)f->Get(Form("%s_mc_cell%s",   v.name, s.Data()));

            if (!hDataBr || !hDataCell || !hMCBr || !hMCCell) continue;
            if (hDataBr->GetEntries() < 100) continue;

            std::vector<TH1D*> histos = {hDataBr, hDataCell, hMCBr, hMCCell};
            std::vector<TString> labels = {
                "Data (branch)", "Data (cell)",
                "MC (branch)", "MC (cell)"
            };
            std::vector<int> colors  = {kBlack, kBlack, kRed, kRed};
            std::vector<int> styles  = {1, 2, 1, 2};
            std::vector<int> markers = {20, 0, 0, 0};

            drawPanel(histos, labels, colors, styles, markers,
                      v.title, n, channelLabel, c);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== Stored vs Computed plots created ===" << std::endl;
    return 0;
}
