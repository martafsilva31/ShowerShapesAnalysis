///////////////////////////////////////////////////////////////////////////////
// plot_mu_scaling.C
//
// Pileup-scaling summary for the eta_mu_loose variant.
//
// Reads histograms.root from the eta_mu_loose output and produces:
//
//   plots/mu_chi2_scaling.pdf
//     Three panels (one per shower-shape variable), each showing chi2/ndf
//     vs <mu> bin centre for three correction methods:
//       - Data vs MC (no correction)
//       - Data vs MC + M1 (flat shift)
//       - Data vs MC + M2 (shift + stretch)
//     One page per eta bin.
//
//   plots/mu_m1_shift_scaling.pdf
//     Three panels (one per variable) showing the mean absolute M1 shift
//     summed over cells vs <mu> bin centre.
//     One page per eta bin.
//
// Parameters:
//   baseDir   root of the eta_mu_loose output, e.g.
//             "../../output/Layer_2/eta_mu_loose"
//   channel   "llgamma" (default) — highest statistics
//   scenario  "inclusive" (default)
//   outDir    directory where PDFs are written; defaults to
//             "<baseDir>/<channel>/<scenario>/plots"
//
// Usage:
//   cd scripts/data_mc
//   root -l -b -q 'plot_mu_scaling.C()'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>

using namespace config;

// ── helpers ──────────────────────────────────────────────────────────────────

static double chi2ndf_local(TH1D* htest, TH1D* href) {
    if (!htest || !href) return -1.;
    double chi2 = 0.; int ndf = 0;
    for (int b = 1; b <= href->GetNbinsX(); ++b) {
        double d  = href->GetBinContent(b);
        double dd = href->GetBinError(b);
        double m  = htest->GetBinContent(b);
        double dm = htest->GetBinError(b);
        double s2 = dd * dd + dm * dm;
        if (s2 <= 0.) continue;
        chi2 += (d - m) * (d - m) / s2;
        ++ndf;
    }
    return ndf > 0 ? chi2 / ndf : -1.;
}

static TH1D* normClone_local(TH1D* h, const char* name) {
    if (!h || h->Integral() <= 0.) return nullptr;
    TH1D* c = (TH1D*)h->Clone(name);
    c->Scale(1. / c->Integral());
    return c;
}

// Draw one page: 3 panels side-by-side showing chi2/ndf vs mu bin for
// each of the 3 shower-shape variables.
static void drawChi2Page(TCanvas* c,
                         double chi2[][kNMuBins][3],  // [var][mu][method 0=MC,1=M1,2=M2]
                         bool   valid[][kNMuBins],    // [var][mu]
                         int    etaBin,
                         const char* scenLabel) {

    const char* varTitles[3]    = {"R_{#eta}", "R_{#phi}", "w_{#eta 2}"};
    const int   methodColor[3]  = {kBlue + 1, kGreen + 2, kRed + 1};
    const int   methodStyle[3]  = {20, 21, 22};
    const char* methodLabel[3]  = {"MC (uncorr.)", "M1 (flat shift)", "M2 (shift+stretch)"};

    c->Clear();
    c->Divide(3, 1, 0.01, 0.01);

    for (int iv = 0; iv < 3; ++iv) {
        c->cd(iv + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.12);
        gPad->SetRightMargin(0.04);

        // Build x-axis: mu bin centres
        double x[kNMuBins], xe[kNMuBins];
        for (int p = 0; p < kNMuBins; ++p) {
            x[p]  = 0.5 * (kMuLimits[p] + kMuLimits[p + 1]);
            xe[p] = 0.5 * (kMuLimits[p + 1] - kMuLimits[p]);
        }

        TMultiGraph* mg = new TMultiGraph();

        for (int im = 0; im < 3; ++im) {
            double y[kNMuBins], ye[kNMuBins];
            int np = 0;
            double xp[kNMuBins], xep[kNMuBins];
            for (int p = 0; p < kNMuBins; ++p) {
                if (!valid[iv][p] || chi2[iv][p][im] < 0.) continue;
                xp[np]  = x[p]; xep[np] = xe[p];
                y[np]   = chi2[iv][p][im];
                ye[np]  = 0.;
                ++np;
            }
            if (np == 0) continue;
            TGraphErrors* gr = new TGraphErrors(np, xp, y, xep, ye);
            gr->SetMarkerStyle(methodStyle[im]);
            gr->SetMarkerColor(methodColor[im]);
            gr->SetLineColor(methodColor[im]);
            gr->SetMarkerSize(1.1);
            gr->SetLineWidth(2);
            mg->Add(gr, "pl");
        }

        mg->SetTitle(Form(";#LT#mu#GT;#chi^{2}/n_{df}  (%s)", varTitles[iv]));
        mg->Draw("a");
        mg->GetXaxis()->SetRangeUser(0., 120.);
        mg->GetXaxis()->SetTitleSize(0.055);
        mg->GetYaxis()->SetTitleSize(0.055);
        mg->GetXaxis()->SetLabelSize(0.050);
        mg->GetYaxis()->SetLabelSize(0.050);

        // Reference line at chi2/ndf = 1
        TLine* l1 = new TLine(0., 1., 120., 1.);
        l1->SetLineStyle(2);
        l1->SetLineColor(kGray + 1);
        l1->Draw();

        // Eta label
        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.048);
        lat.DrawLatex(0.18, 0.90,
                      Form("|#eta| #in [%.2f, %.2f)", kEtaLimits[etaBin], kEtaLimits[etaBin + 1]));
        lat.SetTextSize(0.040);
        lat.DrawLatex(0.18, 0.84, scenLabel);

        // Legend only in first panel
        if (iv == 0) {
            TLegend* leg = new TLegend(0.19, 0.62, 0.96, 0.82);
            leg->SetTextSize(0.044);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            // Rebuild graphs for legend entries
            for (int im = 0; im < 3; ++im) {
                TGraphErrors* gleg = new TGraphErrors(1);
                gleg->SetMarkerStyle(methodStyle[im]);
                gleg->SetMarkerColor(methodColor[im]);
                gleg->SetLineColor(methodColor[im]);
                gleg->SetMarkerSize(1.1);
                leg->AddEntry(gleg, methodLabel[im], "pl");
            }
            leg->Draw();
        }
    }
}

// Draw one page: mean |M1 shift| per cell vs mu bin centre.
static void drawShiftPage(TCanvas* c,
                          double shift[3][kNMuBins],  // [var][mu]
                          bool   valid[3][kNMuBins],  // [var][mu]
                          int    etaBin,
                          const char* scenLabel) {

    const char* varTitles[3]   = {"R_{#eta}", "R_{#phi}", "w_{#eta 2}"};
    const int   varColor[3]    = {kBlue + 1, kRed + 1, kGreen + 2};
    const int   varStyle[3]    = {20, 21, 22};

    c->Clear();
    c->Divide(3, 1, 0.01, 0.01);

    for (int iv = 0; iv < 3; ++iv) {
        c->cd(iv + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.12);
        gPad->SetRightMargin(0.04);

        double x[kNMuBins], xe[kNMuBins], y[kNMuBins], ye[kNMuBins];
        int np = 0;
        for (int p = 0; p < kNMuBins; ++p) {
            if (!valid[iv][p]) continue;
            x[np]  = 0.5 * (kMuLimits[p] + kMuLimits[p + 1]);
            xe[np] = 0.5 * (kMuLimits[p + 1] - kMuLimits[p]);
            y[np]  = shift[iv][p];
            ye[np] = 0.;
            ++np;
        }

        TGraphErrors* gr = new TGraphErrors(np, x, y, xe, ye);
        gr->SetMarkerStyle(varStyle[iv]);
        gr->SetMarkerColor(varColor[iv]);
        gr->SetLineColor(varColor[iv]);
        gr->SetMarkerSize(1.1);
        gr->SetLineWidth(2);
        gr->SetTitle(Form(";#LT#mu#GT;Mean |M1 shift|  (%s)", varTitles[iv]));
        gr->Draw("apl");
        gr->GetXaxis()->SetRangeUser(0., 120.);
        gr->GetXaxis()->SetTitleSize(0.055);
        gr->GetYaxis()->SetTitleSize(0.055);
        gr->GetXaxis()->SetLabelSize(0.050);
        gr->GetYaxis()->SetLabelSize(0.050);

        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.048);
        lat.DrawLatex(0.18, 0.90,
                      Form("|#eta| #in [%.2f, %.2f)", kEtaLimits[etaBin], kEtaLimits[etaBin + 1]));
        lat.SetTextSize(0.040);
        lat.DrawLatex(0.18, 0.84, scenLabel);
    }
}

// ── main function ─────────────────────────────────────────────────────────────

void plot_mu_scaling(
    const char* baseDir  = "../../output/Layer_2/eta_mu_loose",
    const char* channel  = "llgamma",
    const char* scenario = "inclusive",
    const char* outDir   = nullptr)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Resolve output directory
    TString plotsDir;
    if (outDir && strlen(outDir) > 0)
        plotsDir = outDir;
    else
        plotsDir = Form("%s/%s/%s/plots", baseDir, channel, scenario);
    gSystem->mkdir(plotsDir, kTRUE);

    // Open histograms.root
    TString histFile = Form("%s/%s/%s/histograms.root", baseDir, channel, scenario);
    TFile* f = TFile::Open(histFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << histFile << std::endl;
        return;
    }
    std::cout << "Opened: " << histFile << std::endl;

    const char* vars[3]   = {"reta", "rphi", "weta2"};
    const char* scLabel   = scenario;

    // ── chi2/ndf vs mu ─────────────────────────────────────────────────────

    TString chi2Pdf = plotsDir + "/mu_chi2_scaling.pdf";
    TCanvas* cChi2  = new TCanvas("cChi2", "", 1400, 500);
    cChi2->Print(chi2Pdf + "[");  // open PDF

    for (int ie = 0; ie < kNEtaBins; ++ie) {
        double chi2[3][kNMuBins][3] = {};
        bool   valid[3][kNMuBins]   = {};

        for (int iv = 0; iv < 3; ++iv) {
            for (int p = 0; p < kNMuBins; ++p) {
                TString suf = Form("_eta%02d_mu%02d", ie, p);
                TH1D* hd  = (TH1D*)f->Get(Form("h_%s_data%s",  vars[iv], suf.Data()));
                TH1D* hmc = (TH1D*)f->Get(Form("h_%s_mc%s",    vars[iv], suf.Data()));
                TH1D* hm1 = (TH1D*)f->Get(Form("h_%s_mc_M1%s", vars[iv], suf.Data()));
                TH1D* hm2 = (TH1D*)f->Get(Form("h_%s_mc_M2%s", vars[iv], suf.Data()));

                if (!hd || hd->GetEntries() < 50) continue;
                valid[iv][p] = true;

                TH1D* nd  = normClone_local(hd,  "sn_d");
                TH1D* nmc = normClone_local(hmc, "sn_mc");
                TH1D* nm1 = normClone_local(hm1, "sn_m1");
                TH1D* nm2 = normClone_local(hm2, "sn_m2");

                chi2[iv][p][0] = (nd && nmc) ? chi2ndf_local(nmc, nd) : -1.;
                chi2[iv][p][1] = (nd && nm1) ? chi2ndf_local(nm1, nd) : -1.;
                chi2[iv][p][2] = (nd && nm2) ? chi2ndf_local(nm2, nd) : -1.;

                delete nd; delete nmc; delete nm1; delete nm2;
            }
        }

        drawChi2Page(cChi2, chi2, valid, ie, scLabel);
        cChi2->Print(chi2Pdf);
    }

    cChi2->Print(chi2Pdf + "]");  // close PDF
    std::cout << "Written: " << chi2Pdf << std::endl;

    // ── Mean |M1 shift| vs mu ──────────────────────────────────────────────

    TString shiftPdf = plotsDir + "/mu_m1_shift_scaling.pdf";
    TCanvas* cShift  = new TCanvas("cShift", "", 1400, 500);
    cShift->Print(shiftPdf + "[");

    for (int ie = 0; ie < kNEtaBins; ++ie) {
        double shiftMean[3][kNMuBins] = {};
        bool   shiftValid[3][kNMuBins] = {};

        for (int iv = 0; iv < 3; ++iv) {
            for (int p = 0; p < kNMuBins; ++p) {
                // The M1 shift is stored per cell in h_shift_<var>_eta<NN>_mu<PP>
                TString suf = Form("_eta%02d_mu%02d", ie, p);
                TH2D* hShift = (TH2D*)f->Get(Form("h_shift_%s%s", vars[iv], suf.Data()));
                if (!hShift || hShift->GetEntries() < 1) continue;

                // Compute mean |content| over all cells (non-zero entries)
                double sumAbs = 0.; int nCells = 0;
                for (int bx = 1; bx <= hShift->GetNbinsX(); ++bx) {
                    for (int by = 1; by <= hShift->GetNbinsY(); ++by) {
                        double v = hShift->GetBinContent(bx, by);
                        if (v != 0.) { sumAbs += std::fabs(v); ++nCells; }
                    }
                }
                if (nCells > 0) {
                    shiftMean[iv][p]  = sumAbs / nCells;
                    shiftValid[iv][p] = true;
                }
            }
        }

        drawShiftPage(cShift, shiftMean, shiftValid, ie, scLabel);
        cShift->Print(shiftPdf);
    }

    cShift->Print(shiftPdf + "]");
    std::cout << "Written: " << shiftPdf << std::endl;

    f->Close();
    delete f;
    delete cChi2;
    delete cShift;

    std::cout << "\n=== plot_mu_scaling complete ===" << std::endl;
}
