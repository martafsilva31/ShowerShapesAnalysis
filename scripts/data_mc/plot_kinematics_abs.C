///////////////////////////////////////////////////////////////////////////////
// plot_kinematics_abs.C
//
// Plots photon pT and |eta| distributions (Data vs MC, absolute event counts,
// MC scaled by sigma*L/sumW using AMI sumW from config.h) for every
// channel × scenario combination.
//
// Writes NEW pdf files — does NOT overwrite kinematics_pt.pdf / kinematics_eta.pdf:
//   kinematics_pt_abs.pdf
//   kinematics_eta_abs.pdf
//
// Usage:
//   root -l -b -q plot_kinematics_abs.C
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

using namespace config;

// ======================================================================
// Draw a data-vs-MC comparison with ratio panel (absolute event counts)
// ======================================================================
static void drawAbsPanel(TH1D* hData, TH1D* hMC, TH1D* hMCnosf,
                         const char* xTitle,
                         const char* chLabel, const char* scenLabel,
                         TCanvas* c, bool logy = true)
{
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
    hMC->SetFillStyle(0);

    double ymaxAll = std::max(hData->GetMaximum(), hMC->GetMaximum());
    if (hMCnosf) ymaxAll = std::max(ymaxAll, hMCnosf->GetMaximum());
    hData->SetMaximum(1.6 * ymaxAll);
    hData->SetMinimum(logy ? 1e-2 : 0);
    hData->GetXaxis()->SetLabelSize(0);
    hData->GetYaxis()->SetTitle("Events");
    hData->GetYaxis()->SetTitleSize(0.06);
    hData->GetYaxis()->SetLabelSize(0.05);
    hData->GetYaxis()->SetTitleOffset(0.95);
    hData->SetTitle("");
    hData->SetStats(0);

    hData->Draw("E P X0");
    hMC->Draw("HIST SAME");
    if (hMCnosf) {
        hMCnosf->SetLineColor(kBlue + 1);
        hMCnosf->SetLineWidth(2);
        hMCnosf->SetLineStyle(2);
        hMCnosf->SetFillStyle(0);
        hMCnosf->Draw("HIST SAME");
    }
    hData->Draw("E P SAME X0");

    TLegend* leg = new TLegend(0.59, 0.70, 0.96, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(hData, "Data", "ep");
    leg->AddEntry(hMC,   "MC (with SFs)", "l");
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
    lat.DrawLatex(0.18, 0.79, "#sqrt{s} = 13.6 TeV,  108 fb^{-1}");
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

    TH1D* ratio = (TH1D*)hMC->Clone("ratio_abs");
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
    ratio->SetFillStyle(0);
    ratio->Draw("HIST");

    if (hMCnosf) {
        TH1D* ratio_nosf = (TH1D*)hMCnosf->Clone("ratio_abs_nosf");
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
// Draw a data-vs-MC comparison, area-normalised to unity (shape only)
// ======================================================================
static void drawNormPanel(TH1D* hData, TH1D* hMC, TH1D* hMCnosf,
                          const char* xTitle,
                          const char* chLabel, const char* scenLabel,
                          TCanvas* c, bool logy = false)
{
    c->Clear();

    // Clone and normalise to unity
    TH1D* hd = (TH1D*)hData->Clone("hd_norm");
    TH1D* hm = (TH1D*)hMC->Clone("hm_norm");
    if (hd->Integral() > 0) hd->Scale(1.0 / hd->Integral());
    if (hm->Integral() > 0) hm->Scale(1.0 / hm->Integral());

    TH1D* hmn = nullptr;
    if (hMCnosf) {
        hmn = (TH1D*)hMCnosf->Clone("hmn_norm");
        if (hmn->Integral() > 0) hmn->Scale(1.0 / hmn->Integral());
    }

    // --- Upper pad ---
    TPad* pad1 = new TPad("pad1n", "", 0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    if (logy) pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(0.8);
    hd->SetMarkerColor(kBlack);
    hd->SetLineColor(kBlack);
    hd->SetLineWidth(1);

    hm->SetLineColor(kRed);
    hm->SetLineWidth(2);
    hm->SetFillStyle(0);

    double ymaxAll = std::max(hd->GetMaximum(), hm->GetMaximum());
    if (hmn) ymaxAll = std::max(ymaxAll, hmn->GetMaximum());
    hd->SetMaximum(1.6 * ymaxAll);
    hd->SetMinimum(logy ? 1e-6 : 0);
    hd->GetXaxis()->SetLabelSize(0);
    hd->GetYaxis()->SetTitle("Normalised to unity");
    hd->GetYaxis()->SetTitleSize(0.06);
    hd->GetYaxis()->SetLabelSize(0.05);
    hd->GetYaxis()->SetTitleOffset(0.95);
    hd->SetTitle("");
    hd->SetStats(0);

    hd->Draw("E P X0");
    hm->Draw("HIST SAME");
    if (hmn) {
        hmn->SetLineColor(kBlue + 1);
        hmn->SetLineWidth(2);
        hmn->SetLineStyle(2);
        hmn->SetFillStyle(0);
        hmn->Draw("HIST SAME");
    }
    hd->Draw("E P SAME X0");

    TLegend* leg = new TLegend(0.59, 0.70, 0.96, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(hd, "Data", "ep");
    leg->AddEntry(hm, "MC (with SFs)", "l");
    if (hmn)
        leg->AddEntry(hmn, "MC (no SFs)", "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050);
    lat.SetTextFont(72);
    lat.DrawLatex(0.18, 0.85, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(0.285, 0.85, "Work in Progress");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.18, 0.79, "#sqrt{s} = 13.6 TeV,  108 fb^{-1}");
    lat.DrawLatex(0.18, 0.73, chLabel);
    lat.SetTextSize(0.033);
    lat.DrawLatex(0.18, 0.67, scenLabel);

    // --- Ratio pad ---
    c->cd();
    TPad* pad2 = new TPad("pad2n", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1D* ratio = (TH1D*)hm->Clone("ratio_norm");
    ratio->Divide(hd);
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
    ratio->SetFillStyle(0);
    ratio->Draw("HIST");

    if (hmn) {
        TH1D* ratio_nosf = (TH1D*)hmn->Clone("ratio_norm_nosf");
        ratio_nosf->Divide(hd);
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
// Main: loop over all channels × scenarios
// ======================================================================
int plot_kinematics_abs(
    const char* baseDir =
        "../../output/cell_energy_reweighting_Francisco_method/data24")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->SetBatch(true);

    const char* channels[]  = {"eegamma", "mumugamma", "llgamma"};
    const char* scenarios[] = {"baseline", "iso_tight", "converted", "all_conv"};
    const int nCh  = 3;
    const int nSc  = 4;

    int nDone = 0, nSkip = 0;

    TCanvas* c = new TCanvas("c_abs", "Kinematics abs", 800, 700);

    for (int ic = 0; ic < nCh; ic++) {
        for (int is = 0; is < nSc; is++) {
            const char* ch   = channels[ic];
            const char* sc   = scenarios[is];

            TString histoFile = Form("%s/%s/%s/histograms.root", baseDir, ch, sc);
            TString plotDir   = Form("%s/%s/%s/plots/", baseDir, ch, sc);

            TFile* f = TFile::Open(histoFile, "READ");
            if (!f || f->IsZombie()) {
                std::cerr << "SKIP: " << histoFile << " (not found)\n";
                if (f) delete f;
                nSkip++;
                continue;
            }

            gSystem->mkdir(plotDir, true);

            const char* chLabel   = channelLabel(ch);
            const char* scenLabel = scenarioLabel(sc);

            // --- pT ---
            TH1D* hd_pt    = (TH1D*)f->Get("h_pt_data");
            TH1D* hm_pt    = (TH1D*)f->Get("h_pt_mc");
            TH1D* hmn_pt   = (TH1D*)f->Get("h_pt_mc_nosf");
            if (hd_pt && hm_pt) {
                TString pdf = plotDir + "kinematics_pt_abs.pdf";
                drawAbsPanel(hd_pt, hm_pt, hmn_pt,
                             "Photon p_{T} [GeV]", chLabel, scenLabel, c);
                c->SaveAs(pdf);
                std::cout << "Created: " << pdf << "\n";
                nDone++;
                if (TString(sc) == "iso_tight") {
                    TString pdfn = plotDir + "kinematics_pt_norm.pdf";
                    drawNormPanel(hd_pt, hm_pt, hmn_pt,
                                  "Photon p_{T} [GeV]", chLabel, scenLabel, c);
                    c->SaveAs(pdfn);
                    std::cout << "Created: " << pdfn << "\n";
                    nDone++;
                }
            } else {
                std::cerr << "WARNING: h_pt histograms missing in " << histoFile << "\n";
            }

            // --- |eta| ---
            TH1D* hd_eta   = (TH1D*)f->Get("h_eta_data");
            TH1D* hm_eta   = (TH1D*)f->Get("h_eta_mc");
            TH1D* hmn_eta  = (TH1D*)f->Get("h_eta_mc_nosf");
            if (hd_eta && hm_eta) {
                TString pdf = plotDir + "kinematics_eta_abs.pdf";
                drawAbsPanel(hd_eta, hm_eta, hmn_eta,
                             "Photon |#eta|", chLabel, scenLabel, c, false);
                c->SaveAs(pdf);
                std::cout << "Created: " << pdf << "\n";
                nDone++;
                if (TString(sc) == "iso_tight") {
                    TString pdfn = plotDir + "kinematics_eta_norm.pdf";
                    drawNormPanel(hd_eta, hm_eta, hmn_eta,
                                  "Photon |#eta|", chLabel, scenLabel, c, false);
                    c->SaveAs(pdfn);
                    std::cout << "Created: " << pdfn << "\n";
                    nDone++;
                }
            } else {
                std::cerr << "WARNING: h_eta histograms missing in " << histoFile << "\n";
            }

            f->Close();
            delete f;
        }
    }

    delete c;
    std::cout << "\n=== Done: " << nDone << " plots created, "
              << nSkip << " scenarios skipped ===\n";
    return 0;
}
