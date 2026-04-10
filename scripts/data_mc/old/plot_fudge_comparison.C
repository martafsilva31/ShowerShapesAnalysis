// plot_fudge_comparison.C
// Inclusive (all eta bins summed) comparison:
//   Data (branch) vs MC fudged (branch) vs MC unfudged (branch)
// for both Z→eeγ and Z→μμγ channels, three shower shape variables.

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLine.h"
#include <iostream>

void plot_fudge_comparison(
    TString zeegFile   = "../../output/data_mc_reweighting/Zeeg/histos.root",
    TString zmumugFile = "../../output/data_mc_reweighting/Zmumug/histos.root",
    TString outDir     = "../../output/data_mc_reweighting/"
) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const int nEta = 14;
    const int skipBin = 8; // crack region

    // Variable definitions: name prefix, axis label, nbins, xmin, xmax
    struct VarDef {
        const char* prefix;
        const char* label;
        int nbins;
        double xmin, xmax;
    };
    VarDef vars[3] = {
        {"reta",  "R_{#eta}",    25, 0.90, 1.0},
        {"rphi",  "R_{#phi}",    25, 0.85, 1.0},
        {"weta2", "w_{#eta_{2}}", 25, 0.008, 0.016}
    };

    // Channel definitions
    struct Channel {
        const char* name;
        const char* label;
        TString file;
    };
    Channel channels[2] = {
        {"Zeeg",   "Z#rightarrowee#gamma",       zeegFile},
        {"Zmumug", "Z#rightarrow#mu#mu#gamma",   zmumugFile}
    };

    // Lambda: sum histograms across eta bins
    auto sumEtaBins = [&](TFile* f, const char* pattern) -> TH1D* {
        TH1D* hsum = nullptr;
        for (int i = 0; i < nEta; ++i) {
            if (i == skipBin) continue;
            TH1D* h = (TH1D*)f->Get(Form("%s_%d", pattern, i));
            if (!h) continue;
            if (!hsum) {
                hsum = (TH1D*)h->Clone(Form("%s_inclusive", pattern));
                hsum->SetDirectory(0);
            } else {
                hsum->Add(h);
            }
        }
        return hsum;
    };

    // Normalize to unit area
    auto normalize = [](TH1D* h) {
        if (!h) return;
        double integ = h->Integral();
        if (integ > 0) h->Scale(1.0 / integ);
    };

    // Create one PDF per variable with both channels side by side
    for (int iv = 0; iv < 3; ++iv) {
        TCanvas* c = new TCanvas("c", "", 1400, 600);
        c->Divide(2, 1);

        for (int ic = 0; ic < 2; ++ic) {
            TFile* f = TFile::Open(channels[ic].file);
            if (!f || f->IsZombie()) {
                std::cerr << "Cannot open " << channels[ic].file << std::endl;
                continue;
            }

            // Sum across eta bins
            TH1D* h_data = sumEtaBins(f, Form("%s_data_br", vars[iv].prefix));
            TH1D* h_fud  = sumEtaBins(f, Form("%s_mc_fud",  vars[iv].prefix));
            TH1D* h_unf  = sumEtaBins(f, Form("%s_mc_unf",  vars[iv].prefix));

            if (!h_data || !h_fud || !h_unf) {
                std::cerr << "Missing histograms for " << vars[iv].prefix
                          << " in " << channels[ic].name << std::endl;
                f->Close();
                continue;
            }

            normalize(h_data);
            normalize(h_fud);
            normalize(h_unf);

            // Style
            h_data->SetMarkerStyle(20);
            h_data->SetMarkerSize(0.8);
            h_data->SetMarkerColor(kBlack);
            h_data->SetLineColor(kBlack);

            h_fud->SetLineColor(kRed);
            h_fud->SetLineWidth(2);
            h_fud->SetFillColorAlpha(kRed, 0.15);

            h_unf->SetLineColor(kBlue);
            h_unf->SetLineWidth(2);
            h_unf->SetLineStyle(2);

            // Pad with upper plot + ratio
            c->cd(ic + 1);
            TPad* pUp = new TPad(Form("pUp_%d_%d", iv, ic), "", 0, 0.3, 1, 1);
            pUp->SetBottomMargin(0.02);
            pUp->SetLeftMargin(0.14);
            pUp->SetRightMargin(0.04);
            pUp->Draw();
            pUp->cd();

            // Y-axis range
            double ymax = 1.5 * std::max({h_data->GetMaximum(), h_fud->GetMaximum(), h_unf->GetMaximum()});
            h_data->GetYaxis()->SetRangeUser(0, ymax);
            h_data->GetYaxis()->SetTitle("Normalized");
            h_data->GetYaxis()->SetTitleSize(0.06);
            h_data->GetYaxis()->SetLabelSize(0.05);
            h_data->GetXaxis()->SetLabelSize(0);
            h_data->GetXaxis()->SetTitleSize(0);

            h_data->Draw("EP");
            h_fud->Draw("HIST SAME");
            h_unf->Draw("HIST SAME");
            h_data->Draw("EP SAME");

            // Legend
            TLegend* leg = new TLegend(0.16, 0.62, 0.50, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.05);
            leg->AddEntry(h_data, "Data", "lep");
            leg->AddEntry(h_fud,  "MC (fudged)", "lf");
            leg->AddEntry(h_unf,  "MC (unfudged)", "l");
            leg->Draw();

            // Channel label
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.055);
            latex.DrawLatex(0.16, 0.55, channels[ic].label);
            latex.SetTextSize(0.045);
            latex.DrawLatex(0.16, 0.48, "Inclusive #eta");

            // Ratio pad
            c->cd(ic + 1);
            TPad* pLo = new TPad(Form("pLo_%d_%d", iv, ic), "", 0, 0, 1, 0.3);
            pLo->SetTopMargin(0.02);
            pLo->SetBottomMargin(0.35);
            pLo->SetLeftMargin(0.14);
            pLo->SetRightMargin(0.04);
            pLo->Draw();
            pLo->cd();

            // Ratio: MC/Data
            TH1D* r_fud = (TH1D*)h_fud->Clone(Form("r_fud_%d_%d", iv, ic));
            TH1D* r_unf = (TH1D*)h_unf->Clone(Form("r_unf_%d_%d", iv, ic));
            r_fud->Divide(h_data);
            r_unf->Divide(h_data);

            r_fud->SetLineColor(kRed);
            r_fud->SetLineWidth(2);
            r_fud->SetFillColorAlpha(kRed, 0.15);
            r_unf->SetLineColor(kBlue);
            r_unf->SetLineWidth(2);
            r_unf->SetLineStyle(2);

            r_fud->GetYaxis()->SetRangeUser(0.8, 1.2);
            r_fud->GetYaxis()->SetTitle("MC / Data");
            r_fud->GetYaxis()->SetTitleSize(0.12);
            r_fud->GetYaxis()->SetTitleOffset(0.5);
            r_fud->GetYaxis()->SetLabelSize(0.10);
            r_fud->GetYaxis()->SetNdivisions(505);
            r_fud->GetXaxis()->SetTitle(vars[iv].label);
            r_fud->GetXaxis()->SetTitleSize(0.14);
            r_fud->GetXaxis()->SetLabelSize(0.10);

            r_fud->Draw("HIST");
            r_unf->Draw("HIST SAME");

            // Unity line
            TLine* line = new TLine(vars[iv].xmin, 1.0, vars[iv].xmax, 1.0);
            line->SetLineStyle(3);
            line->SetLineColor(kGray + 2);
            line->Draw();

            f->Close();
        }

        TString outName = Form("%sfudge_comparison_%s.pdf", outDir.Data(), vars[iv].prefix);
        c->SaveAs(outName);
        std::cout << "Saved: " << outName << std::endl;
        delete c;
    }
}
