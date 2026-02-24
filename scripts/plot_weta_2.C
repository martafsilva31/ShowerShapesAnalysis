/// \file plot_weta_2.C
/// \brief Overlay fudged / unfudged / computed w_eta_2 with ratio panel.
///
/// Usage:
///   root -l -b -q 'plot_weta_2.C("../output/weta2_Zeeg.root","../output/plots/weta2_Zeeg")'

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>

void plot_weta_2(const char* infile  = "../output/weta2_Zeeg.root",
                 const char* outbase = "../output/plots/weta2_Zeeg") {

    TFile* f = TFile::Open(infile);
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << infile << "\n"; return; }

    TH1F* h_unfudged = (TH1F*)f->Get("h_weta2_unfudged");
    TH1F* h_fudged   = (TH1F*)f->Get("h_weta2_fudged");
    TH1F* h_computed = (TH1F*)f->Get("h_weta2_computed");

    h_unfudged->SetTitle("");
    h_fudged->SetTitle("");
    h_computed->SetTitle("");

    h_unfudged->Scale(1.0 / h_unfudged->Integral());
    h_fudged->Scale(1.0 / h_fudged->Integral());
    h_computed->Scale(1.0 / h_computed->Integral());

    TString fname(infile);
    TString sample = fname(fname.Last('/') + 1, fname.Length());
    sample.ReplaceAll("weta2_", "");
    sample.ReplaceAll(".root", "");

    gStyle->SetOptStat(0);

    h_unfudged->SetLineColor(kBlue);    h_unfudged->SetLineWidth(2);
    h_fudged->SetLineColor(kRed);       h_fudged->SetLineWidth(2);
    h_computed->SetLineColor(kGreen+2); h_computed->SetLineWidth(2);
    h_computed->SetLineStyle(2);

    TCanvas* c = new TCanvas("c", "", 700, 600);

    TPad* pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    h_unfudged->GetXaxis()->SetLabelSize(0);
    h_unfudged->GetXaxis()->SetTitleSize(0);
    h_unfudged->GetYaxis()->SetTitle("Normalised");
    h_unfudged->GetYaxis()->SetTitleSize(0.06);
    h_unfudged->GetYaxis()->SetLabelSize(0.05);
    h_unfudged->GetYaxis()->SetTitleOffset(1.0);
    h_unfudged->Draw("HIST");
    h_fudged->Draw("HIST SAME");
    h_computed->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.16, 0.05, 0.50, 0.25);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(h_unfudged, "Unfudged (stored)", "l");
    leg->AddEntry(h_fudged,   "Fudged (stored)",   "l");
    leg->AddEntry(h_computed, "Computed (cells)",   "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.055);
    lat.SetTextFont(72);
    lat.DrawLatex(0.55, 0.84, "ATLAS");
    lat.SetTextFont(42);
    lat.DrawLatex(0.67, 0.84, "Internal");
    lat.SetTextSize(0.040);
    TString sampleLatex = sample;
    if (sample == "Zeeg")   sampleLatex = "Z#rightarrowee#gamma";
    if (sample == "Zmumug") sampleLatex = "Z#rightarrow#mu#mu#gamma";
    lat.DrawLatex(0.55, 0.78, Form("MC23e %s, #sqrt{s} = 13.6 TeV", sampleLatex.Data()));

    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->Draw();
    pad2->cd();

    TH1F* r_fudged  = (TH1F*)h_fudged->Clone("r_fudged");
    TH1F* r_computed = (TH1F*)h_computed->Clone("r_computed");
    r_fudged->SetTitle("");
    r_computed->SetTitle("");
    r_fudged->Divide(h_unfudged);
    r_computed->Divide(h_unfudged);

    r_fudged->SetMinimum(0.8);
    r_fudged->SetMaximum(1.2);
    r_fudged->GetXaxis()->SetTitle("w_{#eta_{2}}");
    r_fudged->GetXaxis()->SetTitleSize(0.14);
    r_fudged->GetXaxis()->SetLabelSize(0.12);
    r_fudged->GetXaxis()->SetTitleOffset(1.0);
    r_fudged->GetYaxis()->SetTitle("Ratio");
    r_fudged->GetYaxis()->SetTitleSize(0.12);
    r_fudged->GetYaxis()->SetLabelSize(0.10);
    r_fudged->GetYaxis()->SetTitleOffset(0.50);
    r_fudged->GetYaxis()->SetNdivisions(505);
    r_fudged->Draw("HIST");
    r_computed->Draw("HIST SAME");

    TLine line(r_fudged->GetXaxis()->GetXmin(), 1.0,
               r_fudged->GetXaxis()->GetXmax(), 1.0);
    line.SetLineColor(kGray+1);
    line.SetLineStyle(2);
    line.Draw();

    c->SaveAs(Form("%s.pdf", outbase));
    std::cout << "Saved: " << outbase << ".pdf\n";
    f->Close();
}
