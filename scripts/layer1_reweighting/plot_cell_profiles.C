///////////////////////////////////////////////////////////////////////////////
// plot_cell_profiles_layer1.C
//
// Layer-1 version of plot_cell_profiles.C. The cell accumulators span the
// full 56×2 cluster in coordinates relative to the hot strip (111 × 2 = 222
// cells, eta_rel ∈ [-55, +55]). For display we crop to the central ±10
// strips (21 × 2 = 42 cells) where the shower energy is concentrated.
// Produces the same set of PDFs as Layer 2: cell_data, cell_mc, cell_mc_m2,
// cell_shift, cell_stretch and on-the-fly cell_mc_m1 (MC + delta), plus
// per-pT versions in eta_pt.
///////////////////////////////////////////////////////////////////////////////

#include "config_layer1.h"

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

using namespace config_l1;

// Display window: ±kDispEtaHalf strips around the hot strip.
static const int kDispEtaHalf = 10;
static const int kDispEtaSize = 2 * kDispEtaHalf + 1;  // 21

// Reshape a 1D correction histogram (kRelGridSize = 222 bins) into a
// kDispEtaSize × kRelPhiSize TH2D centered on the hot strip.
static TH2D* corrTo2D(TH1D* h1, const char* name) {
    TH2D* h2 = new TH2D(name, "",
                        kDispEtaSize, 0.5, kDispEtaSize + 0.5,
                        kRelPhiSize,  0.5, kRelPhiSize  + 0.5);
    for (int b = 0; b < kDispEtaSize; ++b) {
        int srcEta = kRelCenterEta + (b - kDispEtaHalf);  // 0-based in 111-strip axis
        for (int p = 0; p < kRelPhiSize; ++p) {
            int kr = srcEta * kRelPhiSize + p;
            h2->SetBinContent(b + 1, p + 1, h1->GetBinContent(kr + 1));
        }
    }
    return h2;
}

// Crop a full-size (kRelEtaSize × kRelPhiSize = 111 × 2) cell-profile TH2D
// down to the central kDispEtaSize × kRelPhiSize display window.
static TH2D* cropCentral(TH2D* hIn, const char* name) {
    TH2D* h2 = new TH2D(name, "",
                        kDispEtaSize, 0.5, kDispEtaSize + 0.5,
                        kRelPhiSize,  0.5, kRelPhiSize  + 0.5);
    int cen1 = kRelCenterEta + 1;  // 1-based center bin in source axis
    for (int b = 0; b < kDispEtaSize; ++b) {
        int srcBin = cen1 + (b - kDispEtaHalf);
        for (int p = 1; p <= kRelPhiSize; ++p)
            h2->SetBinContent(b + 1, p, hIn->GetBinContent(srcBin, p));
    }
    return h2;
}

static void drawCellMap(TH2D* hIn,
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
    c->SetTopMargin(ptLabel ? 0.20 : 0.16);

    hIn->SetStats(0); hIn->SetTitle("");
    hIn->GetXaxis()->SetTitle("Strip index  (#eta direction, 0 = hot strip)");
    hIn->GetYaxis()->SetTitle("#phi cell");
    hIn->GetXaxis()->SetTitleSize(0.050);
    hIn->GetYaxis()->SetTitleSize(0.050);
    hIn->GetXaxis()->SetLabelSize(0.045);
    hIn->GetYaxis()->SetLabelSize(0.050);
    hIn->GetXaxis()->SetTitleOffset(1.0);
    hIn->GetYaxis()->SetTitleOffset(1.0);
    hIn->GetZaxis()->SetLabelSize(0.030);
    hIn->GetZaxis()->SetTitleSize(0.045);

    // Signed strip offsets relative to hot strip as x-axis labels.
    for (int i = 1; i <= kDispEtaSize; ++i) {
        int off = (i - 1) - kDispEtaHalf;
        hIn->GetXaxis()->SetBinLabel(i, Form("%d", off));
    }
    for (int j = 1; j <= kRelPhiSize; ++j)
        hIn->GetYaxis()->SetBinLabel(j, Form("%d", j - 1));
    hIn->GetXaxis()->SetNdivisions(kDispEtaSize, false);
    hIn->GetYaxis()->SetNdivisions(kRelPhiSize, false);

    hIn->Draw("COLZ");

    TLatex txt;
    txt.SetTextAlign(22);
    txt.SetTextFont(42);
    txt.SetTextSize(0.022);
    txt.SetTextColor(kBlack);
    for (int ei = 1; ei <= kDispEtaSize; ++ei) {
        double xc = hIn->GetXaxis()->GetBinCenter(ei);
        for (int pi = 1; pi <= kRelPhiSize; ++pi) {
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

    // ---- 3-column header in top margin (NDC), mirroring Layer 2 ----
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);

    double row1Y = ptLabel ? 0.965 : 0.955;
    double row2Y = ptLabel ? 0.920 : 0.895;
    double row3Y = 0.875;

    // Row 1
    lat.SetTextSize(0.038);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.03, row1Y, "#bf{#it{ATLAS}} Work in Progress");
    lat.DrawLatex(0.36, row1Y, chLabel);
    lat.SetTextAlign(31);
    lat.DrawLatex(0.97, row1Y, typeLabel);

    // Row 2
    lat.SetTextSize(0.033);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.03, row2Y, "#sqrt{s} = 13.6 TeV");
    lat.DrawLatex(0.36, row2Y, etaLabel);
    lat.SetTextSize(0.028);
    lat.SetTextAlign(31);
    lat.DrawLatex(0.97, row2Y, scenLabel);

    // Row 3 (pT label, optional)
    if (ptLabel) {
        lat.SetTextSize(0.033);
        lat.SetTextAlign(11);
        lat.DrawLatex(0.36, row3Y, ptLabel);
    }
}

int plot_cell_profiles_layer1(const char* channel   = "llgamma",
                              const char* scenario  = "unconverted",
                              const char* baseDir   = "../../output/Layer_1/eta_loose",
                              const char* binning   = "eta",
                              const char* isolation = "loose") {
    const char* chLabel   = channelLabel(channel);
    const char* scenLabel = scenarioLabel(scenario, isolation);

    TString binMode(binning);
    bool usePtBins = (binMode == "eta_pt");

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

    // Wider canvas than Layer 2 to accommodate the 21-strip eta axis.
    TCanvas* c = new TCanvas("c", "Cell Profiles L1", 1200, 500);

    auto etaLbl = [](int n) -> TString {
        return Form("Bin %d:  %.2f < |#eta| < %.2f",
                    n, kEtaLimits[n], kEtaLimits[n + 1]);
    };
    auto nonEmpty = [](TH2D* h) -> bool {
        if (!h) return false;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
            for (int j = 1; j <= h->GetNbinsY(); ++j)
                if (h->GetBinContent(i, j) != 0) return true;
        return false;
    };
    auto corrSuf = [&](int e, int p = -1) -> TString {
        if (usePtBins && p >= 0) return Form("_eta%02d_pt%02d", e, p);
        return Form("_eta%02d", e);
    };

    struct MapDef {
        const char* fileName;
        const char* typeLabel;
        const char* histPat;
        bool isCorrHist;
    };
    MapDef maps[] = {
        {"cell_data",    "Data",                    "h_frac_mean_data",  false},
        {"cell_mc",      "Original MC",             "h_frac_mean_mc",    false},
        {"cell_mc_m2",   "MC after M2 reweighting", "h_frac_mean_mc_M2", false},
        {"cell_shift",   "Shift correction",        "h_delta",           true },
        {"cell_stretch", "Stretch correction",      "h_stretch",         true },
    };

    // eta-only profiles
    for (auto& m : maps) {
        // In eta_pt mode the correction maps are keyed (eta,pt), so the
        // eta-only correction histograms (h_delta_eta%02d, h_stretch_eta%02d)
        // don't exist → skip to avoid empty/blank PDFs.
        if (usePtBins && m.isCorrHist) continue;

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
                for (int k = 1; k <= kRelGridSize; ++k)
                    if (h1->GetBinContent(k) != 0) { ok = true; break; }
                if (!ok) continue;
                TH2D* h2 = corrTo2D(h1, Form("%s_2d_%d", m.fileName, n));
                drawCellMap(h2, m.typeLabel, chLabel, etaLbl(n), scenLabel, c);
                c->Print(pdf);
                delete h2;
            } else {
                TH2D* h = (TH2D*)f->Get(hname);
                if (!nonEmpty(h)) continue;
                TH2D* hDisp = cropCentral(h, Form("%s_disp_%d", m.fileName, n));
                drawCellMap(hDisp, m.typeLabel, chLabel, etaLbl(n), scenLabel, c);
                c->Print(pdf);
                delete hDisp;
            }
        }
        c->Print(pdf + "]");
    }

    // M1 (MC + delta), computed on the fly from profile + reshape of 1D delta.
    {
        TString pdf = dir + "cell_mc_m1.pdf";
        std::cout << "Creating: " << pdf << std::endl;
        c->Print(pdf + "[");
        for (int n = 0; n < kNEtaBins; ++n) {
            TH2D* hMC    = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc%s", corrSuf(n).Data()));
            TH1D* hDelta = (TH1D*)f->Get(Form("corrections/h_delta%s",          corrSuf(n).Data()));
            if (!nonEmpty(hMC) || !hDelta) continue;
            TH2D* hMCdisp = cropCentral(hMC, Form("cellsM1mc_%d", n));
            TH2D* hDdisp  = corrTo2D(hDelta, Form("cellsM1d_%d", n));
            hMCdisp->Add(hDdisp);
            drawCellMap(hMCdisp, "MC after M1 reweighting", chLabel, etaLbl(n), scenLabel, c);
            c->Print(pdf);
            delete hMCdisp;
            delete hDdisp;
        }
        c->Print(pdf + "]");
    }

    // Per-pT
    if (usePtBins) {
        for (int p = 0; p < kNPtBins; ++p) {
            const char* ptLbl = ptBinLabel(p);
            for (auto& m : maps) {
                TString pdf = dir + Form("%s_pt%02d.pdf", m.fileName, p);
                std::cout << "Creating: " << pdf << std::endl;
                c->Print(pdf + "[");
                for (int n = 0; n < kNEtaBins; ++n) {
                    TString suf = corrSuf(n, p);
                    TString hname;
                    if (m.isCorrHist)
                        hname = Form("corrections/%s%s", m.histPat, suf.Data());
                    else
                        hname = Form("cell_profiles/%s%s", m.histPat, suf.Data());
                    if (m.isCorrHist) {
                        TH1D* h1 = (TH1D*)f->Get(hname);
                        if (!h1) continue;
                        bool ok = false;
                        for (int k = 1; k <= kRelGridSize; ++k)
                            if (h1->GetBinContent(k) != 0) { ok = true; break; }
                        if (!ok) continue;
                        TH2D* h2 = corrTo2D(h1, Form("%s_pt%d_2d_%d", m.fileName, p, n));
                        drawCellMap(h2, m.typeLabel, chLabel, etaLbl(n), scenLabel, c, ptLbl);
                        c->Print(pdf);
                        delete h2;
                    } else {
                        TH2D* h = (TH2D*)f->Get(hname);
                        if (!nonEmpty(h)) continue;
                        TH2D* hDisp = cropCentral(h, Form("%s_pt%d_disp_%d", m.fileName, p, n));
                        drawCellMap(hDisp, m.typeLabel, chLabel, etaLbl(n), scenLabel, c, ptLbl);
                        c->Print(pdf);
                        delete hDisp;
                    }
                }
                c->Print(pdf + "]");
            }

            // M1 per-pT
            TString pdf = dir + Form("cell_mc_m1_pt%02d.pdf", p);
            std::cout << "Creating: " << pdf << std::endl;
            c->Print(pdf + "[");
            for (int n = 0; n < kNEtaBins; ++n) {
                TString suf = corrSuf(n, p);
                TH2D* hMC    = (TH2D*)f->Get(Form("cell_profiles/h_frac_mean_mc%s", suf.Data()));
                TH1D* hDelta = (TH1D*)f->Get(Form("corrections/h_delta%s",          suf.Data()));
                if (!nonEmpty(hMC) || !hDelta) continue;
                TH2D* hMCdisp = cropCentral(hMC, Form("cellsM1mc_pt%d_%d", p, n));
                TH2D* hDdisp  = corrTo2D(hDelta, Form("cellsM1d_pt%d_%d", p, n));
                hMCdisp->Add(hDdisp);
                drawCellMap(hMCdisp, "MC after M1 reweighting", chLabel, etaLbl(n), scenLabel, c, ptLbl);
                c->Print(pdf);
                delete hMCdisp;
                delete hDdisp;
            }
            c->Print(pdf + "]");
        }
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== All cell profile plots created in " << dir << " ===" << std::endl;
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
    return plot_cell_profiles_layer1(ch, sc, bd, bn, iso);
}
#endif
