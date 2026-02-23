///////////////////////////////////////////////////////////////////////////////
// plot_closure.C
//
// Creates comparison plots for the closure test with ratio panels.
// Style follows supervisor's write3HistosRatio.cxx and existing plot_reta.C.
//
// Upper pad: Overlay of pseudo-data (markers), MC (red), M1 (green), M2 (blue)
// Lower pad: Ratio to pseudo-data
//
// Generates multi-page PDFs with one page per eta bin.
//
// Usage:
//   root -l -b -q 'plot_closure.C("histos.root", "plots/")'
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
#include <iostream>

using namespace closure;

// ======================================================================
// Draw one variable for one eta bin
// ======================================================================
void drawVariable(TH1D* hMC, TH1D* hPS, TH1D* hM1, TH1D* hM2,
                  const char* varName, const char* varTitle, int etaBin,
                  TCanvas* c) {

    c->Clear();

    // --- Upper pad: overlay ---
    TPad* pad1 = new TPad("pad1", "", 0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->Draw();
    pad1->cd();

    // Style
    hPS->SetMarkerStyle(20);
    hPS->SetMarkerSize(0.8);
    hPS->SetMarkerColor(kBlack);
    hPS->SetLineColor(kBlack);
    hPS->SetLineWidth(1);

    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(2);
    hMC->SetLineStyle(1);

    hM1->SetLineColor(kGreen + 2);
    hM1->SetLineWidth(2);
    hM1->SetLineStyle(2);

    hM2->SetLineColor(kBlue);
    hM2->SetLineWidth(2);
    hM2->SetLineStyle(1);

    // Find maximum
    double ymax = 1.4 * std::max({hMC->GetMaximum(), hPS->GetMaximum(),
                                  hM1->GetMaximum(), hM2->GetMaximum()});
    hMC->SetMaximum(ymax);
    hMC->SetMinimum(0);
    hMC->GetXaxis()->SetLabelSize(0);
    hMC->GetYaxis()->SetTitle("Normalised");
    hMC->GetYaxis()->SetTitleSize(0.06);
    hMC->GetYaxis()->SetLabelSize(0.05);
    hMC->GetYaxis()->SetTitleOffset(0.95);
    hMC->SetTitle("");
    hMC->SetStats(0);

    hMC->Draw("HIST");
    hM1->Draw("HIST SAME");
    hM2->Draw("HIST SAME");
    hPS->Draw("P SAME");

    // Legend — position adjusted for w_eta2 (distribution peaks left)
    bool isWeta2 = (TString(varName) == "weta2");
    double legX1 = isWeta2 ? 0.16 : 0.55;
    double legX2 = isWeta2 ? 0.52 : 0.88;
    TLegend* leg = new TLegend(legX1, 0.58, legX2, 0.78);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    leg->AddEntry(hPS, "Pseudo-data", "lp");
    leg->AddEntry(hMC, "Original MC", "l");
    leg->AddEntry(hM1, "Reweighted (Shift Only)", "l");
    leg->AddEntry(hM2, "Reweighted (Shift + Stretch)", "l");
    leg->Draw();

    // ATLAS Internal + eta bin label
    double atlasX = isWeta2 ? 0.55 : 0.18;
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050);
    lat.SetTextFont(72);
    lat.DrawLatex(atlasX, 0.86, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(atlasX + 0.12, 0.86, "Internal");
    lat.SetTextSize(0.040);
    lat.DrawLatex(atlasX, 0.80, Form("%s, |#eta| #in [%.2f, %.2f)",
                                   varTitle, kEtaLimits[etaBin], kEtaLimits[etaBin + 1]));

    // --- Lower pad: ratio to pseudo-data ---
    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->Draw();
    pad2->cd();

    // Clone and divide for ratios
    TH1D* ratio_mc = (TH1D*)hMC->Clone(Form("ratio_mc_%s_%d", varName, etaBin));
    TH1D* ratio_m1 = (TH1D*)hM1->Clone(Form("ratio_m1_%s_%d", varName, etaBin));
    TH1D* ratio_m2 = (TH1D*)hM2->Clone(Form("ratio_m2_%s_%d", varName, etaBin));

    ratio_mc->Divide(hPS);
    ratio_m1->Divide(hPS);
    ratio_m2->Divide(hPS);

    ratio_mc->SetMinimum(0.90);
    ratio_mc->SetMaximum(1.10);
    ratio_mc->GetYaxis()->SetTitle("Ratio to PD");
    ratio_mc->GetYaxis()->SetTitleSize(0.12);
    ratio_mc->GetYaxis()->SetTitleOffset(0.45);
    ratio_mc->GetYaxis()->SetLabelSize(0.10);
    ratio_mc->GetYaxis()->SetNdivisions(505);
    ratio_mc->GetXaxis()->SetTitle(varTitle);
    ratio_mc->GetXaxis()->SetTitleSize(0.14);
    ratio_mc->GetXaxis()->SetLabelSize(0.10);
    ratio_mc->SetStats(0);

    ratio_mc->Draw("HIST");
    ratio_m1->Draw("HIST SAME");
    ratio_m2->Draw("HIST SAME");

    // Reference line at 1
    double xmin = ratio_mc->GetXaxis()->GetXmin();
    double xmax = ratio_mc->GetXaxis()->GetXmax();
    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
}

// ======================================================================
// Main
// ======================================================================
int plot_closure(const char* histoFile,
                 const char* plotDir) {

    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);

    TFile* f = TFile::Open(histoFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << histoFile << std::endl;
        return 1;
    }

    TString dir(plotDir);
    if (!dir.EndsWith("/")) dir += "/";

    // Output PDF files
    TString pdfReta  = dir + "closure_reta.pdf";
    TString pdfRphi  = dir + "closure_rphi.pdf";
    TString pdfWeta2 = dir + "closure_weta2.pdf";

    TCanvas* c = new TCanvas("c", "Closure Test", 800, 700);

    // ======================================================================
    // R_eta plots
    // ======================================================================
    std::cout << "Creating R_eta plots..." << std::endl;
    c->Print(pdfReta + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString s = Form("_%d", n);
        TH1D* hMC = (TH1D*)f->Get("reta_mc" + s);
        TH1D* hPS = (TH1D*)f->Get("reta_ps" + s);
        TH1D* hM1 = (TH1D*)f->Get("reta_m1" + s);
        TH1D* hM2 = (TH1D*)f->Get("reta_m2" + s);

        if (!hMC || !hPS || !hM1 || !hM2) continue;
        if (hPS->GetEntries() < 100) continue;

        drawVariable(hMC, hPS, hM1, hM2, "reta", "R_{#eta}", n, c);
        c->Print(pdfReta);
    }

    c->Print(pdfReta + "]");
    std::cout << "  Written: " << pdfReta << std::endl;

    // ======================================================================
    // R_phi plots
    // ======================================================================
    std::cout << "Creating R_phi plots..." << std::endl;
    c->Print(pdfRphi + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString s = Form("_%d", n);
        TH1D* hMC = (TH1D*)f->Get("rphi_mc" + s);
        TH1D* hPS = (TH1D*)f->Get("rphi_ps" + s);
        TH1D* hM1 = (TH1D*)f->Get("rphi_m1" + s);
        TH1D* hM2 = (TH1D*)f->Get("rphi_m2" + s);

        if (!hMC || !hPS || !hM1 || !hM2) continue;
        if (hPS->GetEntries() < 100) continue;

        drawVariable(hMC, hPS, hM1, hM2, "rphi", "R_{#phi}", n, c);
        c->Print(pdfRphi);
    }

    c->Print(pdfRphi + "]");
    std::cout << "  Written: " << pdfRphi << std::endl;

    // ======================================================================
    // w_eta2 plots
    // ======================================================================
    std::cout << "Creating w_eta2 plots..." << std::endl;
    c->Print(pdfWeta2 + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString s = Form("_%d", n);
        TH1D* hMC = (TH1D*)f->Get("weta2_mc" + s);
        TH1D* hPS = (TH1D*)f->Get("weta2_ps" + s);
        TH1D* hM1 = (TH1D*)f->Get("weta2_m1" + s);
        TH1D* hM2 = (TH1D*)f->Get("weta2_m2" + s);

        if (!hMC || !hPS || !hM1 || !hM2) continue;
        if (hPS->GetEntries() < 100) continue;

        drawVariable(hMC, hPS, hM1, hM2, "weta2", "w_{#eta 2}", n, c);
        c->Print(pdfWeta2);
    }

    c->Print(pdfWeta2 + "]");
    std::cout << "  Written: " << pdfWeta2 << std::endl;

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== All plots created ===" << std::endl;

    return 0;
}
