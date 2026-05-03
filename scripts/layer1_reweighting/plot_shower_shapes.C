///////////////////////////////////////////////////////////////////////////////
// plot_shower_shapes_layer1.C
//
// Layer-1 version of plot_shower_shapes.C. Produces rew_<var>.pdf,
// rew_integrated.pdf, computed_vs_stored(_eta).pdf, fudge_factors(_eta).pdf
// and optional per-pT versions for the five strip shower shapes:
//   weta1, wstot, fside, deltae, eratio.
///////////////////////////////////////////////////////////////////////////////

#include "config_layer1.h"

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

using namespace config_l1;

// ======================================================================
// Generic panel (identical logic to Layer 2)
// ======================================================================
static void drawPanel(std::vector<TH1D*>& histos,
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
    TH1D* hRef = histos[0];

    TPad* pad1 = new TPad("pad1", "", 0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->Draw(); pad1->cd();

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

    double ymax = 0;
    for (auto& h : histos) ymax = std::max(ymax, h->GetMaximum());
    ymax *= 1.6;

    double peak = 0;
    for (auto& h : histos) peak = std::max(peak, h->GetMaximum());
    double threshold = 0.005 * peak;
    int axFirst = hRef->GetXaxis()->GetFirst();
    int axLast  = hRef->GetXaxis()->GetLast();
    int firstBin = axFirst, lastBin = axLast;
    for (int b = axFirst; b <= axLast; ++b) {
        bool any = false;
        for (auto& h : histos) if (h->GetBinContent(b) > threshold) { any = true; break; }
        if (any) { firstBin = b; break; }
    }
    for (int b = axLast; b >= axFirst; --b) {
        bool any = false;
        for (auto& h : histos) if (h->GetBinContent(b) > threshold) { any = true; break; }
        if (any) { lastBin = b; break; }
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
    hRef->SetTitle(""); hRef->SetStats(0);

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

    TLegend* leg = new TLegend(legX1, legY1, legX2, legY2);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
    for (size_t i = 0; i < histos.size(); ++i)
        leg->AddEntry(histos[i], labels[i], (i == 0) ? "ep" : "lf");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.050); lat.SetTextFont(72);
    lat.DrawLatex(0.18, 0.85, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(0.285, 0.85, "Work in Progress");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.18, 0.79, "#sqrt{s} = 13.6 TeV");
    lat.DrawLatex(0.18, 0.73, chLabel);
    lat.SetTextSize(0.033);
    lat.DrawLatex(0.18, 0.67, scenLabel);
    if (overrideEtaLabel)
        lat.DrawLatex(0.18, 0.61, overrideEtaLabel);
    else
        lat.DrawLatex(0.18, 0.61,
            Form("%.2f < |#eta| < %.2f",
                 kEtaLimits[etaBin], kEtaLimits[etaBin + 1]));
    if (ptLabel) { lat.SetTextSize(0.033); lat.DrawLatex(0.18, 0.55, ptLabel); }

    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.30);
    pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14); pad2->SetGridy();
    pad2->Draw(); pad2->cd();

    bool firstRatio = true;
    for (size_t i = 1; i < histos.size(); ++i) {
        TH1D* ratio = (TH1D*)histos[i]->Clone(Form("ratio_%zu_%d", i, etaBin));
        ratio->Divide(hRef);
        ratio->SetMinimum(0.50); ratio->SetMaximum(1.50);
        ratio->GetXaxis()->SetRangeUser(xlo, xhi);
        ratio->GetYaxis()->SetTitle("MC/Data");
        ratio->GetYaxis()->SetTitleSize(0.12);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitle(varTitle);
        ratio->GetXaxis()->SetTitleSize(0.14);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->SetTitle(""); ratio->SetStats(0);

        TH1D* rBand = (TH1D*)ratio->Clone(Form("ratio_band_%zu_%d", i, etaBin));
        rBand->SetFillColorAlpha(colors[i], 0.15);
        rBand->SetLineWidth(0);
        if (firstRatio) { rBand->Draw("E2"); ratio->Draw("HIST SAME"); }
        else            { rBand->Draw("E2 SAME"); ratio->Draw("HIST SAME"); }
        firstRatio = false;
    }
    TLine* line = new TLine(xlo, 1.0, xhi, 1.0);
    line->SetLineColor(kGray + 2); line->SetLineStyle(2); line->SetLineWidth(1);
    line->Draw();
}

static void computeDataRange(TH1D* href, double& xLo, double& xHi,
                             double pctLo = 0.005, double pctHi = 0.995) {
    if (!href) { xLo = 0; xHi = 1; return; }
    TAxis* ax = href->GetXaxis();
    if (href->GetEntries() < 10) { xLo = ax->GetXmin(); xHi = ax->GetXmax(); return; }
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
int plot_shower_shapes_layer1(const char* channel   = "llgamma",
                              const char* scenario  = "unconverted",
                              const char* baseDir   = "../../output/Layer_1/eta_loose",
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
        TH1D* hProbe = (TH1D*)f->Get("h_weta1_data_computed");
        if (hProbe) origNBins = hProbe->GetNbinsX();
    }
    const int rebinGroups = std::max(1, origNBins / nZoomBins);

    TString dir(plotDir);
    if (!dir.EndsWith("/")) dir += "/";
    gSystem->mkdir(dir, true);

    TCanvas* c = new TCanvas("c", "Data-MC Comparison", 800, 700);

    struct VarDef { const char* name; const char* title; };
    VarDef vars[] = {
        {"weta1",  "w_{#eta 1}"},
        {"wstot",  "w_{s,tot}"},
        {"fside",  "f_{side}"},
        {"deltae", "#DeltaE [MeV]"},
        {"eratio", "E_{ratio}"}
    };

    auto zoomNorm = [&](std::vector<TH1D*>& hv, TH1D* hRef) {
        for (auto*& h : hv) {
            h = (TH1D*)h->Clone(Form("%s_zm", h->GetName()));
            h->Rebin(rebinGroups);
        }
        double xLo, xHi;
        computeDataRange(hRef, xLo, xHi);
        for (auto* h : hv) h->GetXaxis()->SetRangeUser(xLo, xHi);
        for (auto* h : hv) if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
    };

    // --- Set B per-eta ---
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
            TH1D* hM3   = (TH1D*)f->Get(Form("h_%s_mc_M3%s",   v.name, s.Data()));
            TH1D* hM4   = (TH1D*)f->Get(Form("h_%s_mc_M4%s",   v.name, s.Data()));
            if (!hData || !hMC || !hM1 || !hM2) continue;
            if (hData->GetEntries() < 100) continue;

            std::vector<TH1D*>    hv     = {hData, hMC, hM1, hM2};
            std::vector<TString>  labels = {"Data", "Original MC",
                                             "MC reweighted (shift only)",
                                             "MC reweighted (shift+stretch)"};
            std::vector<int>      colors = {kBlack, kRed, kGreen + 2, kBlue};
            std::vector<int>      styles = {1, 1, 1, 1};
            if (hM3) { hv.push_back(hM3); labels.push_back("MC M3 (quantile)");
                       colors.push_back(kOrange + 7); styles.push_back(2); }
            if (hM4) { hv.push_back(hM4); labels.push_back("MC M4 (quant.+reshuffle)");
                       colors.push_back(kMagenta + 1); styles.push_back(1); }
            zoomNorm(hv, hv[0]);
            drawPanel(hv, labels, colors, styles,
                      v.title, n, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.52, 0.60, 0.89, 0.87);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // --- Set B per-(eta,pT) ---
    if (usePtBins) {
        for (int p = 0; p < kNPtBins; ++p) {
            for (auto& v : vars) {
                TString pdfPath = dir + Form("rew_%s_pt%02d.pdf", v.name, p);
                std::cout << "Set B per-pT " << v.title << " pT " << p
                          << ": " << pdfPath << std::endl;
                c->Print(pdfPath + "[");
                for (int n = 0; n < kNEtaBins; ++n) {
                    TString s = Form("_eta%02d_pt%02d", n, p);
                    TH1D* hData = (TH1D*)f->Get(Form("h_%s_data%s",  v.name, s.Data()));
                    TH1D* hMC   = (TH1D*)f->Get(Form("h_%s_mc%s",    v.name, s.Data()));
                    TH1D* hM1   = (TH1D*)f->Get(Form("h_%s_mc_M1%s", v.name, s.Data()));
                    TH1D* hM2   = (TH1D*)f->Get(Form("h_%s_mc_M2%s", v.name, s.Data()));
                    TH1D* hM3   = (TH1D*)f->Get(Form("h_%s_mc_M3%s", v.name, s.Data()));
                    TH1D* hM4   = (TH1D*)f->Get(Form("h_%s_mc_M4%s", v.name, s.Data()));
                    if (!hData || !hMC || !hM1 || !hM2) continue;
                    if (hData->GetEntries() < 50) continue;

                    std::vector<TH1D*>    hv     = {hData, hMC, hM1, hM2};
                    std::vector<TString>  labels = {"Data", "Original MC",
                                                    "MC reweighted (shift only)",
                                                    "MC reweighted (shift+stretch)"};
                    std::vector<int>      colors = {kBlack, kRed, kGreen + 2, kBlue};
                    std::vector<int>      styles = {1, 1, 1, 1};
                    if (hM3) { hv.push_back(hM3); labels.push_back("MC M3 (quantile)");
                               colors.push_back(kOrange + 7); styles.push_back(2); }
                    if (hM4) { hv.push_back(hM4); labels.push_back("MC M4 (quant.+reshuffle)");
                               colors.push_back(kMagenta + 1); styles.push_back(1); }
                    zoomNorm(hv, hv[0]);
                    drawPanel(hv, labels, colors, styles,
                              v.title, n, chLabel, scenLabel, c,
                              nullptr, ptBinLabel(p), 0.52, 0.60, 0.89, 0.87);
                    c->Print(pdfPath);
                }
                c->Print(pdfPath + "]");
            }
        }
    }

    // --- Set B integrated ---
    {
        TString pdfPath = dir + "rew_integrated.pdf";
        std::cout << "Set B integrated: " << pdfPath << std::endl;
        c->Print(pdfPath + "[");
        for (auto& v : vars) {
            TH1D* hData = (TH1D*)f->Get(Form("h_%s_data_computed", v.name));
            TH1D* hMC   = (TH1D*)f->Get(Form("h_%s_mc_computed",   v.name));
            TH1D* hM1   = (TH1D*)f->Get(Form("h_%s_mc_M1",         v.name));
            TH1D* hM2   = (TH1D*)f->Get(Form("h_%s_mc_M2",         v.name));
            TH1D* hM3   = (TH1D*)f->Get(Form("h_%s_mc_M3",         v.name));
            TH1D* hM4   = (TH1D*)f->Get(Form("h_%s_mc_M4",         v.name));
            if (!hData || !hMC || !hM1 || !hM2) continue;

            std::vector<TH1D*>    hv     = {hData, hMC, hM1, hM2};
            std::vector<TString>  labels = {"Data", "Original MC",
                                             "MC reweighted (shift only)",
                                             "MC reweighted (shift+stretch)"};
            std::vector<int>      colors = {kBlack, kRed, kGreen + 2, kBlue};
            std::vector<int>      styles = {1, 1, 1, 1};
            if (hM3) { hv.push_back(hM3); labels.push_back("MC M3 (quantile)");
                       colors.push_back(kOrange + 7); styles.push_back(2); }
            if (hM4) { hv.push_back(hM4); labels.push_back("MC M4 (quant.+reshuffle)");
                       colors.push_back(kMagenta + 1); styles.push_back(1); }
            zoomNorm(hv, hv[0]);
            drawPanel(hv, labels, colors, styles,
                      v.title, 0, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.52, 0.60, 0.89, 0.87);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // --- Computed vs stored (integrated) ---
    {
        TString pdfPath = dir + "computed_vs_stored.pdf";
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
                      nullptr, nullptr, 0.62, 0.70, 0.96, 0.92);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // --- Computed vs stored (per-eta) ---
    {
        TString pdfPath = dir + "computed_vs_stored_eta.pdf";
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
                          nullptr, nullptr, 0.62, 0.70, 0.96, 0.92);
                c->Print(pdfPath);
            }
        }
        c->Print(pdfPath + "]");
    }

    // --- Fudge factors (integrated) ---
    {
        TString pdfPath = dir + "fudge_factors.pdf";
        c->Print(pdfPath + "[");
        for (auto& v : vars) {
            TH1D* hStor = (TH1D*)f->Get(Form("h_%s_data_stored", v.name));
            TH1D* hFud  = (TH1D*)f->Get(Form("h_%s_mc_fudged",   v.name));
            TH1D* hUnf  = (TH1D*)f->Get(Form("h_%s_mc_unfudged", v.name));
            if (!hStor || !hFud || !hUnf) continue;
            std::vector<TH1D*> hv = {hStor, hFud, hUnf};
            zoomNorm(hv, hv[0]);
            std::vector<TString> labels = {"Data", "MC (fudged)", "MC (unfudged)"};
            std::vector<int> colors = {kBlack, kBlue, kRed};
            std::vector<int> styles = {1, 1, 2};
            drawPanel(hv, labels, colors, styles,
                      v.title, 0, chLabel, scenLabel, c,
                      nullptr, nullptr, 0.59, 0.62, 0.96, 0.88);
            c->Print(pdfPath);
        }
        c->Print(pdfPath + "]");
    }

    // --- Fudge factors (per-eta) ---
    {
        TString pdfPath = dir + "fudge_factors_eta.pdf";
        c->Print(pdfPath + "[");
        for (auto& v : vars) {
            for (int n = 0; n < kNEtaBins; ++n) {
                TString s = Form("_eta%02d", n);
                TH1D* hStor = (TH1D*)f->Get(Form("h_%s_data_stored%s", v.name, s.Data()));
                TH1D* hFud  = (TH1D*)f->Get(Form("h_%s_mc_fudged%s",   v.name, s.Data()));
                TH1D* hUnf  = (TH1D*)f->Get(Form("h_%s_mc_unfudged%s", v.name, s.Data()));
                if (!hStor || !hFud || !hUnf) continue;
                if (hStor->GetEntries() < 100) continue;
                std::vector<TH1D*> hv = {hStor, hFud, hUnf};
                zoomNorm(hv, hv[0]);
                std::vector<TString> labels = {"Data", "MC (fudged)", "MC (unfudged)"};
                std::vector<int> colors = {kBlack, kBlue, kRed};
                std::vector<int> styles = {1, 1, 2};
                drawPanel(hv, labels, colors, styles,
                          v.title, n, chLabel, scenLabel, c,
                          nullptr, nullptr, 0.59, 0.62, 0.96, 0.88);
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

#ifndef __CLING__
#include <TApplication.h>
int main(int argc, char** argv) {
    TApplication app("app", nullptr, nullptr);
    gROOT->SetBatch(true);
    const char* ch  = (argc > 1) ? argv[1] : "llgamma";
    const char* sc  = (argc > 2) ? argv[2] : "unconverted";
    const char* bd  = (argc > 3) ? argv[3] : "../../output/Layer_1/eta_loose";
    const char* bn  = (argc > 4) ? argv[4] : "eta";
    const char* iso = (argc > 5) ? argv[5] : "loose";
    return plot_shower_shapes_layer1(ch, sc, bd, bn, iso);
}
#endif
