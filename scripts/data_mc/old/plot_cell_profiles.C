///////////////////////////////////////////////////////////////////////////////
// plot_cell_profiles.C
//
// Produces 2D cell energy profile maps (7 eta x 11 phi) for data and MC,
// before and after reweighting (M1, M2, M3), for all eta bins.
//
// Reads the mean cell fractions stored in corrections.root by
// derive_corrections.C:
//   cellsMeanData_Eta_LO_HI   (TH2D: 11 phi x 7 eta)
//   cellsMeanMC_Eta_LO_HI     (TH2D: 11 phi x 7 eta)
//   cellsSigmaRatio_Eta_LO_HI (TH2D: 11 phi x 7 eta)
//   cellsDelta_Eta_LO_HI      (TH2D: 11 phi x 7 eta)  [M1 shift]
//   m2_alpha_Cell_K_Eta_LO_HI (TProfile)             [M2 DeltaR-matched]
//
// Produces PDFs:
//   cell_profiles_before.pdf — Data vs MC (uncorrected), one page per eta bin
//   cell_profiles_after_m1.pdf — Data vs MC after M1 reweighting
//   cell_profiles_after_m2.pdf — Data vs MC after M2 reweighting
//   cell_profiles_after_m3.pdf — Data vs MC after M3 reweighting
//   cell_profiles_delta.pdf — Delta maps and sigma ratio
//
// After-reweighting profiles are computed analytically:
//   M1: <f'_k> = <f_MC_k> + Delta_k = <f_data_k>   (by construction)
//   M2: <f'_k> = <f_MC_k> + m2_prof->Interpolate(<f_MC_k>) via matched TProfile
//   M3: <f'_k> = <f_data_k>                          (by construction)
//
// Usage:
//   root -l -b -q 'plot_cell_profiles.C("corrections.root", "plots/")'
//   root -l -b -q 'plot_cell_profiles.C("corrections.root", "plots/", "Z#rightarrowee#gamma")'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TColor.h>
#include <iostream>
#include <cmath>

using namespace datamc;


// ======================================================================
// Draw a single 2D cell profile map with text annotations
// ======================================================================
void drawCellMap(TH2D* hIn, const char* title, TPad* pad) {
    pad->cd();
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.16);
    pad->SetTopMargin(0.15);
    pad->SetBottomMargin(0.12);

    // Transpose: input is (phi x eta), note style is (eta-x, phi-y)
    TH2D* h = new TH2D(
        Form("%s_tr", hIn->GetName()), "",
        kEtaSize, 0.5, kEtaSize + 0.5,
        kPhiSize, 0.5, kPhiSize + 0.5);

    for (int phi = 1; phi <= kPhiSize; ++phi)
        for (int eta = 1; eta <= kEtaSize; ++eta)
            h->SetBinContent(eta, phi, hIn->GetBinContent(phi, eta));

    h->GetXaxis()->SetTitle("Cell (#eta direction)");
    h->GetYaxis()->SetTitle("Cell (#phi direction)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetXaxis()->SetNdivisions(kEtaSize, 0, 0);
    h->GetYaxis()->SetNdivisions(kPhiSize, 0, 0);
    h->SetStats(0);
    h->SetTitle("");

    // Color scale: same range for both panels
    h->Draw("COLZ");

    // Annotate each cell with its value
    TLatex txt;
    txt.SetTextAlign(22);
    txt.SetTextFont(42);

    for (int eta = 1; eta <= kEtaSize; ++eta) {
        for (int phi = 1; phi <= kPhiSize; ++phi) {
            double val = h->GetBinContent(eta, phi);
            double aval = std::fabs(val);

            // Adaptive formatting
            TString label;
            if (aval < 1e-5)
                label = Form("%.3e", val);
            else if (aval < 0.001)
                label = Form("%.4e", val);
            else if (aval < 0.01)
                label = Form("%.4f", val);
            else if (aval < 0.1)
                label = Form("%.4f", val);
            else
                label = Form("%.4f", val);

            // Size: smaller text for the 7x11 grid to fit
            double textSize = 0.022;
            txt.SetTextSize(textSize);

            // Color: white on dark cells, black on light cells
            double zmin = h->GetMinimum();
            double zmax = h->GetMaximum();
            double frac = (zmax > zmin) ? (val - zmin) / (zmax - zmin) : 0.5;
            if (frac > 0.4 && frac < 0.7)
                txt.SetTextColor(kWhite);
            else
                txt.SetTextColor(kBlack);

            txt.DrawLatex(eta, phi, label);
        }
    }

    // Title
    TLatex header;
    header.SetNDC();
    header.SetTextFont(42);
    header.SetTextSize(0.045);
    header.DrawLatex(0.15, 0.90, title);
}


// ======================================================================
// Main
// ======================================================================
int plot_cell_profiles(const char* corrFile,
                       const char* plotDir,
                       const char* channelLabel = "Z#rightarrowll#gamma") {

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".4g");
    gStyle->SetPalette(kBird);
    gROOT->SetBatch(true);

    TFile* f = TFile::Open(corrFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << corrFile << std::endl;
        return 1;
    }

    TString dir(plotDir);
    if (!dir.EndsWith("/")) dir += "/";
    gSystem->mkdir(dir, true);

    TCanvas* c = new TCanvas("c", "Cell Profiles", 1400, 700);

    // ================================================================
    // PDF 1: Before reweighting — Data vs MC (uncorrected)
    // ================================================================
    TString pdfBefore = dir + "cell_profiles_before.pdf";
    std::cout << "Creating: " << pdfBefore << std::endl;
    c->Print(pdfBefore + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* hData = (TH2D*)f->Get("cellsMeanData_" + suffix);
        TH2D* hMC   = (TH2D*)f->Get("cellsMeanMC_" + suffix);

        if (!hData || !hMC) {
            std::cout << "  Skipping eta bin " << n
                      << " — missing histograms" << std::endl;
            continue;
        }

        // Check if bin has data
        bool hasData = false;
        for (int i = 1; i <= hData->GetNbinsX(); ++i)
            for (int j = 1; j <= hData->GetNbinsY(); ++j)
                if (hData->GetBinContent(i, j) != 0) { hasData = true; break; }
        if (!hasData) continue;

        c->Clear();
        c->Divide(2, 1);

        // Common Z range for both panels
        double zmin = std::min(hData->GetMinimum(), hMC->GetMinimum());
        double zmax = std::max(hData->GetMaximum(), hMC->GetMaximum());

        // Left panel: Data
        drawCellMap(hData,
                    Form("(a) Data  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         channelLabel, kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

        // Right panel: MC uncorrected
        drawCellMap(hMC,
                    Form("(b) MC without reweighting  |  |#eta| #in [%.2f, %.2f)",
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(2));

        c->Print(pdfBefore);
    }
    c->Print(pdfBefore + "]");

    // ================================================================
    // PDF 2: After M1 reweighting
    // M1: <f'_k> = <f_MC_k> + Delta_k = <f_data_k> by construction
    // ================================================================
    TString pdfM1 = dir + "cell_profiles_after_m1.pdf";
    std::cout << "Creating: " << pdfM1 << std::endl;
    c->Print(pdfM1 + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* hData = (TH2D*)f->Get("cellsMeanData_" + suffix);
        TH2D* hMC   = (TH2D*)f->Get("cellsMeanMC_" + suffix);
        TH2D* hDelta = (TH2D*)f->Get("cellsDelta_" + suffix);

        if (!hData || !hMC || !hDelta) continue;

        bool hasData = false;
        for (int i = 1; i <= hData->GetNbinsX(); ++i)
            for (int j = 1; j <= hData->GetNbinsY(); ++j)
                if (hData->GetBinContent(i, j) != 0) { hasData = true; break; }
        if (!hasData) continue;

        // M1 corrected = MC + delta = data (by construction)
        TH2D* hM1 = (TH2D*)hMC->Clone(Form("cellsM1_%d", n));
        hM1->Add(hDelta);

        c->Clear();
        c->Divide(2, 1);

        drawCellMap(hData,
                    Form("(a) Data  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         channelLabel, kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

        drawCellMap(hM1,
                    Form("(b) MC after M1 reweighting  |  |#eta| #in [%.2f, %.2f)",
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(2));

        c->Print(pdfM1);
        delete hM1;
    }
    c->Print(pdfM1 + "]");

    // ================================================================
    // PDF 3: After M2 reweighting
    // M2: <f'_k> = <f_MC_k> + alpha(f_MC_k) from DeltaR-matched TProfile
    // ================================================================
    TString pdfM2 = dir + "cell_profiles_after_m2.pdf";
    std::cout << "Creating: " << pdfM2 << std::endl;
    c->Print(pdfM2 + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* hData = (TH2D*)f->Get("cellsMeanData_" + suffix);
        TH2D* hMC   = (TH2D*)f->Get("cellsMeanMC_" + suffix);

        if (!hData || !hMC) continue;

        bool hasData = false;
        for (int i = 1; i <= hData->GetNbinsX(); ++i)
            for (int j = 1; j <= hData->GetNbinsY(); ++j)
                if (hData->GetBinContent(i, j) != 0) { hasData = true; break; }
        if (!hasData) continue;

        // M2 corrected: mc_frac + alpha(mc_frac) from TProfile
        TH2D* hM2 = (TH2D*)hMC->Clone(Form("cellsM2_%d", n));
        bool allFound = true;
        for (int phi = 1; phi <= kPhiSize; ++phi) {
            for (int eta = 1; eta <= kEtaSize; ++eta) {
                int k = (phi - 1) + kPhiSize * (eta - 1);
                TString pname = Form("m2_alpha_Cell_%d_Eta_%1.2f_%1.2f",
                                     k + 1, kEtaLimits[n], kEtaLimits[n + 1]);
                TProfile* prof = (TProfile*)f->Get(pname);
                if (!prof) { allFound = false; continue; }
                double mc_frac = hMC->GetBinContent(phi, eta);
                double alpha = prof->Interpolate(mc_frac);
                hM2->SetBinContent(phi, eta, mc_frac + alpha);
            }
        }
        if (!allFound && n == 0)
            std::cout << "  WARNING: Some M2 TProfiles missing" << std::endl;

        c->Clear();
        c->Divide(2, 1);

        drawCellMap(hData,
                    Form("(a) Data  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         channelLabel, kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

        drawCellMap(hM2,
                    Form("(b) MC after M2 reweighting  |  |#eta| #in [%.2f, %.2f)",
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(2));

        c->Print(pdfM2);
        delete hM2;
    }
    c->Print(pdfM2 + "]");

    // ================================================================
    // PDF 4: After M3 reweighting — Data vs MC (corrected means = data)
    // By construction: <f'_k> = mu_data_k
    // ================================================================
    TString pdfAfter = dir + "cell_profiles_after_m3.pdf";
    std::cout << "Creating: " << pdfAfter << std::endl;
    c->Print(pdfAfter + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* hData = (TH2D*)f->Get("cellsMeanData_" + suffix);
        TH2D* hMC   = (TH2D*)f->Get("cellsMeanMC_" + suffix);

        if (!hData || !hMC) continue;

        bool hasData = false;
        for (int i = 1; i <= hData->GetNbinsX(); ++i)
            for (int j = 1; j <= hData->GetNbinsY(); ++j)
                if (hData->GetBinContent(i, j) != 0) { hasData = true; break; }
        if (!hasData) continue;

        // M3-corrected MC mean = mu_data by construction
        // (before renormalization — which shifts values slightly)
        TH2D* hM3 = (TH2D*)hData->Clone(Form("cellsM3_%d", n));

        c->Clear();
        c->Divide(2, 1);

        drawCellMap(hData,
                    Form("(a) Data  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         channelLabel, kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

        drawCellMap(hM3,
                    Form("(b) MC after M3 reweighting  |  |#eta| #in [%.2f, %.2f)",
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(2));

        c->Print(pdfAfter);
        delete hM3;
    }
    c->Print(pdfAfter + "]");

    // ================================================================
    // PDF 3: Delta maps — (Data - MC) before and sigma ratio
    // ================================================================
    TString pdfDelta = dir + "cell_profiles_delta.pdf";
    std::cout << "Creating: " << pdfDelta << std::endl;
    c->Print(pdfDelta + "[");

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* hData   = (TH2D*)f->Get("cellsMeanData_" + suffix);
        TH2D* hMC     = (TH2D*)f->Get("cellsMeanMC_" + suffix);
        TH2D* hSigR   = (TH2D*)f->Get("cellsSigmaRatio_" + suffix);

        if (!hData || !hMC || !hSigR) continue;

        bool hasData = false;
        for (int i = 1; i <= hData->GetNbinsX(); ++i)
            for (int j = 1; j <= hData->GetNbinsY(); ++j)
                if (hData->GetBinContent(i, j) != 0) { hasData = true; break; }
        if (!hasData) continue;

        // Compute difference map
        TH2D* hDiff = (TH2D*)hData->Clone(Form("cellsDiff_%d", n));
        hDiff->Add(hMC, -1.0);

        c->Clear();
        c->Divide(2, 1);

        drawCellMap(hDiff,
                    Form("(a) #Delta = Data - MC  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         channelLabel, kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

        drawCellMap(hSigR,
                    Form("(b) #sigma_{data}/#sigma_{MC}  |  |#eta| #in [%.2f, %.2f)",
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(2));

        c->Print(pdfDelta);
        delete hDiff;
    }
    c->Print(pdfDelta + "]");

    // ================================================================
    // PDF 6-8: Per-method residual (Data - MC_corrected)
    // Shows how well each method matches data at cell level.
    // M1, M3 should be ~0 by construction; M2 shows actual residual.
    // ================================================================
    const char* mNames[] = {"m1", "m2", "m3"};
    const char* mLabels[] = {"M1 (flat shift)", "M2 (profiled)", "M3 (shift+stretch)"};

    for (int mi = 0; mi < 3; ++mi) {
        TString pdfDM = dir + Form("cell_delta_%s.pdf", mNames[mi]);
        std::cout << "Creating: " << pdfDM << std::endl;
        c->Print(pdfDM + "[");

        for (int n = 0; n < kNEtaBins; ++n) {
            TString suffix = Form("Eta_%1.2f_%1.2f",
                                  kEtaLimits[n], kEtaLimits[n + 1]);

            TH2D* hData = (TH2D*)f->Get("cellsMeanData_" + suffix);
            TH2D* hMC   = (TH2D*)f->Get("cellsMeanMC_" + suffix);
            if (!hData || !hMC) continue;

            bool hasData = false;
            for (int i = 1; i <= hData->GetNbinsX(); ++i)
                for (int j = 1; j <= hData->GetNbinsY(); ++j)
                    if (hData->GetBinContent(i, j) != 0) { hasData = true; break; }
            if (!hasData) continue;

            // Compute corrected MC depending on method
            TH2D* hCorr = (TH2D*)hMC->Clone(Form("cellsCorr_%s_%d", mNames[mi], n));
            if (mi == 0) {
                // M1: MC + delta = data
                TH2D* hDelta = (TH2D*)f->Get("cellsDelta_" + suffix);
                if (hDelta) hCorr->Add(hDelta);
            } else if (mi == 1) {
                // M2: DeltaR-matched TProfile correction
                for (int phi = 1; phi <= kPhiSize; ++phi) {
                    for (int eta = 1; eta <= kEtaSize; ++eta) {
                        int k = (phi - 1) + kPhiSize * (eta - 1);
                        TString pname = Form("m2_alpha_Cell_%d_Eta_%1.2f_%1.2f",
                                             k + 1, kEtaLimits[n], kEtaLimits[n + 1]);
                        TProfile* prof = (TProfile*)f->Get(pname);
                        double mc_frac = hMC->GetBinContent(phi, eta);
                        double alpha = prof ? prof->Interpolate(mc_frac) : 0;
                        hCorr->SetBinContent(phi, eta, mc_frac + alpha);
                    }
                }
            } else {
                // M3: use actual shift+stretch formula on mean map
                //   corrected = mu_data + sigma_ratio * (mu_mc - mu_mc) = mu_data
                // Mean residual is zero by construction.  Show sigma ratio instead.
                TH2D* hDelta = (TH2D*)f->Get("cellsDelta_" + suffix);
                if (hDelta) hCorr->Add(hDelta);
            }

            if (mi < 2) {
                // M1, M2: show residual = Data - Corrected
                TH2D* hResid = (TH2D*)hData->Clone(
                    Form("cellsResid_%s_%d", mNames[mi], n));
                hResid->Add(hCorr, -1.0);

                c->Clear();
                c->Divide(2, 1);

                drawCellMap(hCorr,
                    Form("(a) MC after %s  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         mLabels[mi], channelLabel,
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

                drawCellMap(hResid,
                    Form("(b) Data - MC(%s)  |  |#eta| #in [%.2f, %.2f)",
                         mLabels[mi],
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(2));

                delete hResid;
            } else {
                // M3: mean residual is zero by design (shift guarantees it).
                // Show sigma ratio map instead — the meaningful M3 diagnostic.
                TH2D* hSigRatio = (TH2D*)f->Get("cellsSigmaRatio_" + suffix);

                c->Clear();
                c->Divide(2, 1);

                drawCellMap(hCorr,
                    Form("(a) MC after %s  |  %s  |  |#eta| #in [%.2f, %.2f)",
                         mLabels[mi], channelLabel,
                         kEtaLimits[n], kEtaLimits[n + 1]),
                    (TPad*)c->cd(1));

                if (hSigRatio) {
                    drawCellMap(hSigRatio,
                        Form("(b) #sigma_{data}/#sigma_{MC}  |  |#eta| #in [%.2f, %.2f)",
                             kEtaLimits[n], kEtaLimits[n + 1]),
                        (TPad*)c->cd(2));
                }
            }

            c->Print(pdfDM);
            delete hCorr;
        }

        c->Print(pdfDM + "]");
    }

    f->Close();
    delete f;
    delete c;

    std::cout << "\n=== Cell profile plots created ===" << std::endl;
    return 0;
}
