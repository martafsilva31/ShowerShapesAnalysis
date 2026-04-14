///////////////////////////////////////////////////////////////////////////////
// plot_kinematics.C
//
// Plots photon pT and |eta| distributions (data vs MC) and the
// zero-padding diagnostic histogram from fill_histograms.C output.
//
// Produces:
//   kinematics_pt.pdf      — photon pT
//   kinematics_eta.pdf     — photon |eta|
//   zero_padding.pdf       — number of zero cells per event
//
// Usage:
//   root -l -b -q 'plot_kinematics.C("eegamma", "baseline")'
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
#include <cmath>

using namespace config;

// ======================================================================
// Draw a single data-vs-MC comparison with ratio panel
// ======================================================================
void drawKinPanel(TH1D* hData, TH1D* hMC, TH1D* hMCnosf,
                  const char* xTitle, const char* yTitle,
                  const char* chLabel, const char* scenLabel,
                  TCanvas* c, bool logy = false) {

    c->Clear();

    // --- Upper pad ---
    TPad* pad1 = new TPad("pad1", "", 0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    if (logy) pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.8);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(1);

    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(2);

    TH1D* hMCband = (TH1D*)hMC->Clone("hMCband");
    hMCband->SetFillColorAlpha(kRed, 0.15);
    hMCband->SetLineWidth(0);

    double ymaxAll = std::max(hData->GetMaximum(), hMC->GetMaximum());
    if (hMCnosf) ymaxAll = std::max(ymaxAll, hMCnosf->GetMaximum());
    hData->SetMaximum(1.6 * ymaxAll);
    hData->SetMinimum(logy ? 1e-5 : 0);
    hData->GetXaxis()->SetLabelSize(0);
    hData->GetYaxis()->SetTitle(yTitle);
    hData->GetYaxis()->SetTitleSize(0.06);
    hData->GetYaxis()->SetLabelSize(0.05);
    hData->GetYaxis()->SetTitleOffset(0.95);
    hData->SetTitle("");
    hData->SetStats(0);

    hData->Draw("E P X0");
    hMCband->Draw("E2 SAME");
    hMC->Draw("HIST SAME");
    if (hMCnosf) {
        hMCnosf->SetLineColor(kBlue + 1);
        hMCnosf->SetLineWidth(2);
        hMCnosf->SetLineStyle(2);
        hMCnosf->Draw("HIST SAME");
    }
    hData->Draw("E P SAME X0");

    TLegend* leg = new TLegend(0.48, 0.66, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(hData, "Data", "ep");
    leg->AddEntry(hMC,   "MC (with SFs)", "lf");
    if (hMCnosf)
        leg->AddEntry(hMCnosf, "MC (no SFs)", "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050);
    lat.SetTextFont(72);
    lat.DrawLatex(0.18, 0.85, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(0.285, 0.85, "Work in Progress");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.18, 0.79, "#sqrt{s} = 13.6 TeV");
    lat.DrawLatex(0.18, 0.73, chLabel);
    lat.SetTextSize(0.033);
    lat.DrawLatex(0.18, 0.67, scenLabel);

    // --- Ratio pad ---
    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1D* ratio = (TH1D*)hMC->Clone("ratio_kin");
    ratio->Divide(hData);
    ratio->SetMinimum(0.50);
    ratio->SetMaximum(1.50);
    ratio->GetYaxis()->SetTitle("MC/Data");
    ratio->GetYaxis()->SetTitleSize(0.12);
    ratio->GetYaxis()->SetTitleOffset(0.45);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetTitle(xTitle);
    ratio->GetXaxis()->SetTitleSize(0.14);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->SetTitle("");
    ratio->SetStats(0);
    ratio->SetFillColorAlpha(kRed, 0.15);
    ratio->Draw("E2");
    ratio->Draw("HIST SAME");

    if (hMCnosf) {
        TH1D* ratio_nosf = (TH1D*)hMCnosf->Clone("ratio_kin_nosf");
        ratio_nosf->Divide(hData);
        ratio_nosf->SetLineColor(kBlue + 1);
        ratio_nosf->SetLineWidth(2);
        ratio_nosf->SetLineStyle(2);
        ratio_nosf->SetFillStyle(0);
        ratio_nosf->Draw("HIST SAME");
    }

    TLine* line = new TLine(
        ratio->GetXaxis()->GetXmin(), 1.0,
        ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->Draw();
}

// ======================================================================
// Main
// ======================================================================
int plot_kinematics(const char* channel  = "eegamma",
                    const char* scenario = "baseline",
                    const char* baseDir  =
                        "../../output/cell_energy_reweighting_Francisco_method/data24") {

    const char* chLabel   = channelLabel(channel);
    const char* scenLabel = scenarioLabel(scenario);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gROOT->SetBatch(true);

    TString outPath = Form("%s/%s/%s", baseDir, channel, scenario);
    TString histoFile = outPath + "/histograms.root";
    TString plotDir   = outPath + "/plots/";
    gSystem->mkdir(plotDir, true);

    TFile* f = TFile::Open(histoFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << histoFile << std::endl;
        return 1;
    }

    TCanvas* c = new TCanvas("c", "Kinematics", 800, 700);

    // --- pT distribution ---
    {
        TH1D* hd    = (TH1D*)f->Get("h_pt_data");
        TH1D* hm    = (TH1D*)f->Get("h_pt_mc");
        TH1D* hmnsf = (TH1D*)f->Get("h_pt_mc_nosf");
        if (hd && hm) {
            TString pdf = plotDir + "kinematics_pt.pdf";
            drawKinPanel(hd, hm, hmnsf, "Photon p_{T} [GeV]", "Events",
                         chLabel, scenLabel, c);
            c->SaveAs(pdf);
            std::cout << "Created: " << pdf << std::endl;
        } else {
            std::cerr << "WARNING: pt histograms not found in " << histoFile << std::endl;
        }
    }

    // --- |eta| distribution ---
    {
        TH1D* hd    = (TH1D*)f->Get("h_eta_data");
        TH1D* hm    = (TH1D*)f->Get("h_eta_mc");
        TH1D* hmnsf = (TH1D*)f->Get("h_eta_mc_nosf");
        if (hd && hm) {
            TString pdf = plotDir + "kinematics_eta.pdf";
            drawKinPanel(hd, hm, hmnsf, "Photon |#eta|", "Events",
                         chLabel, scenLabel, c);
            c->SaveAs(pdf);
            std::cout << "Created: " << pdf << std::endl;
        } else {
            std::cerr << "WARNING: eta histograms not found in " << histoFile << std::endl;
        }
    }

    // --- Zero-padding diagnostic ---
    {
        TH1D* hd = (TH1D*)f->Get("h_nzero_data");
        TH1D* hm = (TH1D*)f->Get("h_nzero_mc");
        if (hd && hm) {
            TString pdf = plotDir + "zero_padding.pdf";

            c->Clear();
            c->SetLogy(true);
            c->SetLeftMargin(0.14);

            // Normalise
            if (hd->Integral() > 0) hd->Scale(1.0 / hd->Integral());
            if (hm->Integral() > 0) hm->Scale(1.0 / hm->Integral());

            hd->SetMarkerStyle(20);
            hd->SetMarkerSize(0.8);
            hd->SetMarkerColor(kBlack);
            hd->SetLineColor(kBlack);
            hm->SetLineColor(kRed);
            hm->SetLineWidth(2);

            double ymax = 2.0 * std::max(hd->GetMaximum(), hm->GetMaximum());
            hd->SetMaximum(ymax);
            hd->SetMinimum(1e-6);
            hd->GetXaxis()->SetTitle("Number of cells with E = 0 (excl. central)");
            hd->GetYaxis()->SetTitle("Fraction of events");
            hd->GetXaxis()->SetRangeUser(-0.5, 20.5);
            hd->SetTitle("");
            hd->SetStats(0);

            hd->Draw("E P X0");
            hm->Draw("HIST SAME");
            hd->Draw("E P SAME X0");

            TLegend* leg = new TLegend(0.55, 0.75, 0.89, 0.89);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hd, "Data", "ep");
            leg->AddEntry(hm, "MC", "l");
            leg->Draw();

            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.045);
            lat.SetTextFont(72);
            lat.DrawLatex(0.18, 0.85, "ATLAS");
            lat.SetTextFont(42);
            lat.DrawLatex(0.275, 0.85, "Work in Progress");
            lat.SetTextSize(0.033);
            lat.DrawLatex(0.18, 0.79, chLabel);
            lat.DrawLatex(0.18, 0.73, scenLabel);

            c->SaveAs(pdf);
            c->SetLogy(false);
            std::cout << "Created: " << pdf << std::endl;
        } else {
            std::cerr << "WARNING: zero-padding histograms not found in " << histoFile << std::endl;
        }
    }

    f->Close();
    delete f;
    delete c;
    std::cout << "\n=== Kinematics plots done ===" << std::endl;
    return 0;
}
