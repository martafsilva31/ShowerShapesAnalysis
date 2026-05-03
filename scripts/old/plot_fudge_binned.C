// plot_fudge_binned.C
// Comprehensive fudge factor validation: data vs MC unfudged vs MC fudged
// Uses TUNE27AD η×pT binning (matching the actual fudge factor derivation).
// Three views: inclusive, η-binned (summed over pT), full η×pT grid.
// Three channels: Z→eeγ, Z→μμγ, Combined.
// Also updates the old fudge_comparison_*.pdf files.

#include "config.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLine.h"
#include "TSystem.h"
#include <iostream>
#include <algorithm>
#include <vector>

using namespace datamc;

// ===================================================================
// TUNE27AD binning (from EGammaVariableCorrection/TUNE27AD/*.conf)
// ===================================================================
namespace tune27 {
    const int nEta = 8;
    const double etaEdge[9] = {0., 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37};
    const int nPt = 11;
    const double ptEdge[12] = {
        0., 15e3, 20e3, 30e3, 40e3, 60e3, 80e3, 100e3, 150e3, 300e3, 600e3, 3e6
    }; // MeV (ntuple units)
    const int skipEta = 4; // crack: 1.37–1.52

    int findEta(double abseta) {
        for (int i = 0; i < nEta; ++i)
            if (abseta >= etaEdge[i] && abseta < etaEdge[i + 1]) return i;
        return -1;
    }
    int findPt(double ptMeV) {
        for (int i = 0; i < nPt; ++i)
            if (ptMeV >= ptEdge[i] && ptMeV < ptEdge[i + 1]) return i;
        return -1;
    }
    TString etaLabel(int i) {
        return Form("%.2f < |#eta| < %.2f", etaEdge[i], etaEdge[i + 1]);
    }
    TString ptLabel(int i) {
        double lo = ptEdge[i] / 1e3, hi = ptEdge[i + 1] / 1e3;
        if (hi > 1000) return Form("p_{T} > %.0f GeV", lo);
        return Form("%.0f < p_{T} < %.0f GeV", lo, hi);
    }
}

// ===================================================================
// Variable definitions
// ===================================================================
struct VarDef { const char* name; const char* label; int nb; double lo, hi; };
static const int nVar = 3;
static const VarDef kVars[nVar] = {
    {"reta",  "R_{#eta}",      25, 0.90,  1.0},
    {"rphi",  "R_{#phi}",      25, 0.85,  1.0},
    {"weta2", "w_{#eta_{2}}",  25, 0.008, 0.016}
};

// ===================================================================
// Histogram container — [var][etaBin][ptBin]
// ===================================================================
struct Histos {
    TH1D* h[nVar][tune27::nEta][tune27::nPt];

    void create(const char* tag) {
        for (int v = 0; v < nVar; ++v)
            for (int e = 0; e < tune27::nEta; ++e)
                for (int p = 0; p < tune27::nPt; ++p) {
                    h[v][e][p] = new TH1D(
                        Form("%s_%s_e%d_p%d", tag, kVars[v].name, e, p),
                        "", kVars[v].nb, kVars[v].lo, kVars[v].hi);
                    h[v][e][p]->Sumw2();
                    h[v][e][p]->SetDirectory(0);
                }
    }

    void fill(int v, int e, int p, double val, double w = 1.) {
        if (v >= 0 && v < nVar && e >= 0 && e < tune27::nEta
            && p >= 0 && p < tune27::nPt)
            h[v][e][p]->Fill(val, w);
    }

    // Sum over all pT bins for a given (var, eta)
    TH1D* sumPt(int v, int e) const {
        TH1D* s = nullptr;
        for (int p = 0; p < tune27::nPt; ++p) {
            if (!s) {
                s = (TH1D*)h[v][e][p]->Clone(Form("spt_%d_%d_%p", v, e, (void*)this));
                s->SetDirectory(0);
            } else {
                s->Add(h[v][e][p]);
            }
        }
        return s;
    }

    // Sum over all η (skip crack) and all pT for a given var
    TH1D* sumAll(int v) const {
        TH1D* s = nullptr;
        for (int e = 0; e < tune27::nEta; ++e) {
            if (e == tune27::skipEta) continue;
            for (int p = 0; p < tune27::nPt; ++p) {
                if (!s) {
                    s = (TH1D*)h[v][e][p]->Clone(Form("sall_%d_%p", v, (void*)this));
                    s->SetDirectory(0);
                } else {
                    s->Add(h[v][e][p]);
                }
            }
        }
        return s;
    }

    // Add another Histos set
    void add(const Histos& o) {
        for (int v = 0; v < nVar; ++v)
            for (int e = 0; e < tune27::nEta; ++e)
                for (int p = 0; p < tune27::nPt; ++p)
                    h[v][e][p]->Add(o.h[v][e][p]);
    }
};

// ===================================================================
// Helpers
// ===================================================================
static void norm(TH1D* h) {
    double s = h->Integral();
    if (s > 0) h->Scale(1.0 / s);
}

static double chi2ndf(TH1D* h1, TH1D* h2) {
    double chi2 = 0;
    int ndf = 0;
    for (int b = 1; b <= h1->GetNbinsX(); ++b) {
        double v1 = h1->GetBinContent(b), v2 = h2->GetBinContent(b);
        double e1 = h1->GetBinError(b), e2 = h2->GetBinError(b);
        double den = e1 * e1 + e2 * e2;
        if (den > 0) { chi2 += (v1 - v2) * (v1 - v2) / den; ndf++; }
    }
    return ndf > 0 ? chi2 / ndf : -1;
}

// ===================================================================
// Draw one comparison panel (upper pad + ratio)
// ===================================================================
static int gPadCtr = 0;  // unique pad names

void drawPanel(TVirtualPad* mother, int padIdx,
               TH1D* hd_raw, TH1D* hu_raw, TH1D* hf_raw,
               const char* varLabel, const char* binLabel, bool showLegend) {
    ++gPadCtr;
    mother->cd(padIdx);

    // Clone and normalize
    TH1D* hd = (TH1D*)hd_raw->Clone(Form("hd_%d", gPadCtr)); hd->SetDirectory(0);
    TH1D* hu = (TH1D*)hu_raw->Clone(Form("hu_%d", gPadCtr)); hu->SetDirectory(0);
    TH1D* hf = (TH1D*)hf_raw->Clone(Form("hf_%d", gPadCtr)); hf->SetDirectory(0);
    norm(hd); norm(hu); norm(hf);

    // Style: data = black points, unfudged = red, fudged = blue
    hd->SetMarkerStyle(20); hd->SetMarkerSize(0.7);
    hd->SetMarkerColor(kBlack); hd->SetLineColor(kBlack);
    hu->SetLineColor(kRed + 1); hu->SetLineWidth(2);
    hf->SetLineColor(kBlue + 1); hf->SetLineWidth(2);

    // Upper pad
    TPad* pUp = new TPad(Form("up%d", gPadCtr), "", 0, 0.3, 1, 1);
    pUp->SetBottomMargin(0.02);
    pUp->SetLeftMargin(0.15);
    pUp->SetRightMargin(0.04);
    pUp->Draw();
    pUp->cd();

    double ymax = 1.55 * std::max({hd->GetMaximum(), hu->GetMaximum(), hf->GetMaximum()});
    hd->GetYaxis()->SetRangeUser(0, ymax);
    hd->GetYaxis()->SetTitle("Normalized");
    hd->GetYaxis()->SetTitleSize(0.06);
    hd->GetYaxis()->SetLabelSize(0.05);
    hd->GetXaxis()->SetLabelSize(0);
    hd->GetXaxis()->SetTitleSize(0);
    hd->Draw("EP");
    hu->Draw("HIST SAME");
    hf->Draw("HIST SAME");
    hd->Draw("EP SAME"); // redraw on top

    if (showLegend) {
        TLegend* leg = new TLegend(0.55, 0.65, 0.93, 0.93);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.055);
        leg->AddEntry(hd, "Data", "lep");
        leg->AddEntry(hu, "MC (unfudged)", "l");
        leg->AddEntry(hf, "MC (fudged)", "l");
        leg->Draw();
    }

    // Labels
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.050);
    tex.DrawLatex(0.17, 0.88, binLabel);
    // Chi2/ndf
    double c_u = chi2ndf(hd, hu), c_f = chi2ndf(hd, hf);
    tex.SetTextSize(0.042);
    tex.SetTextColor(kRed + 1);
    tex.DrawLatex(0.17, 0.81, Form("#chi^{2}/ndf (unf) = %.2f", c_u));
    tex.SetTextColor(kBlue + 1);
    tex.DrawLatex(0.17, 0.74, Form("#chi^{2}/ndf (fud) = %.2f", c_f));
    tex.SetTextColor(kBlack);
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.17, 0.67, Form("N_{data} = %.0f", hd_raw->GetEntries()));

    // Ratio pad
    mother->cd(padIdx);
    TPad* pLo = new TPad(Form("lo%d", gPadCtr), "", 0, 0, 1, 0.3);
    pLo->SetTopMargin(0.02);
    pLo->SetBottomMargin(0.35);
    pLo->SetLeftMargin(0.15);
    pLo->SetRightMargin(0.04);
    pLo->Draw();
    pLo->cd();

    TH1D* r_u = (TH1D*)hu->Clone(Form("ru%d", gPadCtr)); r_u->SetDirectory(0);
    TH1D* r_f = (TH1D*)hf->Clone(Form("rf%d", gPadCtr)); r_f->SetDirectory(0);
    r_u->Divide(hd);
    r_f->Divide(hd);

    r_u->GetYaxis()->SetRangeUser(0.75, 1.25);
    r_u->GetYaxis()->SetTitle("MC / Data");
    r_u->GetYaxis()->SetTitleSize(0.12);
    r_u->GetYaxis()->SetTitleOffset(0.45);
    r_u->GetYaxis()->SetLabelSize(0.10);
    r_u->GetYaxis()->SetNdivisions(505);
    r_u->GetXaxis()->SetTitle(varLabel);
    r_u->GetXaxis()->SetTitleSize(0.14);
    r_u->GetXaxis()->SetLabelSize(0.10);
    r_u->Draw("HIST");
    r_f->Draw("HIST SAME");

    double xlo = hd->GetXaxis()->GetXmin(), xhi = hd->GetXaxis()->GetXmax();
    TLine* line = new TLine(xlo, 1.0, xhi, 1.0);
    line->SetLineStyle(3);
    line->SetLineColor(kGray + 2);
    line->Draw();
}

// ===================================================================
// Fill histograms from ntuple
// ===================================================================
void fillData(const char* filename, Histos& hData, Long64_t maxEv = -1) {
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }
    TTree* t = (TTree*)f->Get(kTreeName);
    if (!t) { std::cerr << "No tree in " << filename << std::endl; return; }

    Float_t photonPt = 0, eta2 = 0;
    Double_t mll = 0, mllg = 0, dr1 = 0, dr2 = 0, drll = 0;
    Int_t cellSize = 0;
    Int_t tightID = 0;
    Bool_t isConv = false;
    Float_t reta = 0, rphi = 0, weta2 = 0;

    t->SetBranchStatus("*", 0);
    auto on = [&](const char* name, void* addr) {
        if (t->GetBranch(name)) { t->SetBranchStatus(name, 1); t->SetBranchAddress(name, addr); }
    };
    on(kPhotonPtBranch, &photonPt);
    on(kEtaBranch, &eta2);
    on(kMllBranch, &mll);
    on(kMllgBranch, &mllg);
    on(kDRLepton1Branch, &dr1);
    on(kDRLepton2Branch, &dr2);
    on(kDRllBranch, &drll);
    on(kCellSizeBranch, &cellSize);
    on(kTightIDBranch, &tightID);
    on(kIsConvBranch, &isConv);
    on(kRetaFudBranch, &reta);   // for data, reta = standard reco value
    on(kRphiFudBranch, &rphi);
    on(kWeta2FudBranch, &weta2);

    Long64_t nEntries = t->GetEntries();
    if (maxEv > 0 && maxEv < nEntries) nEntries = maxEv;
    std::cout << "  Data: " << filename << " (" << nEntries << " entries)" << std::endl;

    Long64_t nUsed = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (i > 0 && i % 2000000 == 0)
            std::cout << "    " << i << " / " << nEntries << std::endl;

        if (!passSelection(photonPt, eta2, mll, mllg, dr1, dr2, drll,
                           cellSize, true, false, tightID, isConv))
            continue;

        double abseta = std::fabs(eta2);
        int eb = tune27::findEta(abseta);
        int pb = tune27::findPt(photonPt);
        if (eb < 0 || pb < 0) continue;

        hData.fill(0, eb, pb, reta);
        hData.fill(1, eb, pb, rphi);
        hData.fill(2, eb, pb, weta2);
        nUsed++;
    }
    std::cout << "  Selected: " << nUsed << " data events" << std::endl;
    f->Close();
}

void fillMC(const char* filename, Histos& hUnf, Histos& hFud, Long64_t maxEv = -1) {
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }
    TTree* t = (TTree*)f->Get(kTreeName);
    if (!t) { std::cerr << "No tree in " << filename << std::endl; return; }

    Float_t photonPt = 0, eta2 = 0;
    Double_t mll = 0, mllg = 0, dr1 = 0, dr2 = 0, drll = 0;
    Int_t cellSize = 0;
    Int_t tightID = 0;
    Bool_t isConv = false;
    Bool_t truthMatch = false;
    Float_t mcwgt = 1, muwgt = 1;
    Float_t reta_f = 0, rphi_f = 0, weta2_f = 0;
    Float_t reta_u = 0, rphi_u = 0, weta2_u = 0;

    t->SetBranchStatus("*", 0);
    auto on = [&](const char* name, void* addr) {
        if (t->GetBranch(name)) { t->SetBranchStatus(name, 1); t->SetBranchAddress(name, addr); }
    };
    on(kPhotonPtBranch, &photonPt);
    on(kEtaBranch, &eta2);
    on(kMllBranch, &mll);
    on(kMllgBranch, &mllg);
    on(kDRLepton1Branch, &dr1);
    on(kDRLepton2Branch, &dr2);
    on(kDRllBranch, &drll);
    on(kCellSizeBranch, &cellSize);
    on(kTightIDBranch, &tightID);
    on(kIsConvBranch, &isConv);
    on(kTruthMatchBranch, &truthMatch);
    on(kMCWeightBranch, &mcwgt);
    on(kPUWeightBranch, &muwgt);
    on(kRetaFudBranch, &reta_f);
    on(kRphiFudBranch, &rphi_f);
    on(kWeta2FudBranch, &weta2_f);
    on(kRetaUnfBranch, &reta_u);
    on(kRphiUnfBranch, &rphi_u);
    on(kWeta2UnfBranch, &weta2_u);

    Long64_t nEntries = t->GetEntries();
    if (maxEv > 0 && maxEv < nEntries) nEntries = maxEv;
    std::cout << "  MC: " << filename << " (" << nEntries << " entries)" << std::endl;

    Long64_t nUsed = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (i > 0 && i % 2000000 == 0)
            std::cout << "    " << i << " / " << nEntries << std::endl;

        if (!passSelection(photonPt, eta2, mll, mllg, dr1, dr2, drll,
                           cellSize, truthMatch, true, tightID, isConv))
            continue;

        double abseta = std::fabs(eta2);
        int eb = tune27::findEta(abseta);
        int pb = tune27::findPt(photonPt);
        if (eb < 0 || pb < 0) continue;

        double w = muwgt * mcwgt;
        hUnf.fill(0, eb, pb, reta_u, w);
        hUnf.fill(1, eb, pb, rphi_u, w);
        hUnf.fill(2, eb, pb, weta2_u, w);
        hFud.fill(0, eb, pb, reta_f, w);
        hFud.fill(1, eb, pb, rphi_f, w);
        hFud.fill(2, eb, pb, weta2_f, w);
        nUsed++;
    }
    std::cout << "  Selected: " << nUsed << " MC events" << std::endl;
    f->Close();
}

// ===================================================================
// Plotting functions
// ===================================================================
void plotInclusive(const Histos& data, const Histos& unf, const Histos& fud,
                   const char* outFile, const char* chanLabel) {
    TCanvas c("c_inc", "", 1500, 600);
    c.Divide(3, 1);
    for (int v = 0; v < nVar; ++v) {
        TH1D* hd = data.sumAll(v);
        TH1D* hu = unf.sumAll(v);
        TH1D* hf = fud.sumAll(v);
        if (!hd || hd->GetEntries() < 10) continue;
        TString lbl = Form("%s  Inclusive", chanLabel);
        drawPanel(&c, v + 1, hd, hu, hf, kVars[v].label, lbl.Data(), v == 0);
    }
    c.SaveAs(outFile);
    std::cout << "  Saved: " << outFile << std::endl;
}

void plotEtaBins(const Histos& data, const Histos& unf, const Histos& fud,
                 const char* outFile, const char* chanLabel) {
    TCanvas c("c_eta", "", 1500, 600);
    c.Print(Form("%s[", outFile)); // open multi-page PDF

    for (int e = 0; e < tune27::nEta; ++e) {
        if (e == tune27::skipEta) continue;
        c.Clear();
        c.Divide(3, 1);
        bool drawn = false;
        for (int v = 0; v < nVar; ++v) {
            TH1D* hd = data.sumPt(v, e);
            TH1D* hu = unf.sumPt(v, e);
            TH1D* hf = fud.sumPt(v, e);
            if (!hd || hd->GetEntries() < 10) continue;
            TString lbl = Form("%s  %s", chanLabel, tune27::etaLabel(e).Data());
            drawPanel(&c, v + 1, hd, hu, hf, kVars[v].label, lbl.Data(), v == 0);
            drawn = true;
        }
        if (drawn) c.Print(outFile);
    }

    c.Print(Form("%s]", outFile)); // close PDF
    std::cout << "  Saved: " << outFile << std::endl;
}

void plotEtaPtBins(const Histos& data, const Histos& unf, const Histos& fud,
                   const char* outFile, const char* chanLabel) {
    TCanvas c("c_ep", "", 1500, 600);
    c.Print(Form("%s[", outFile));
    const int minEntries = 50;

    for (int e = 0; e < tune27::nEta; ++e) {
        if (e == tune27::skipEta) continue;
        for (int p = 0; p < tune27::nPt; ++p) {
            // Check if any variable has enough data
            bool ok = false;
            for (int v = 0; v < nVar; ++v)
                if (data.h[v][e][p]->GetEntries() >= minEntries) { ok = true; break; }
            if (!ok) continue;

            c.Clear();
            c.Divide(3, 1);
            for (int v = 0; v < nVar; ++v) {
                if (data.h[v][e][p]->GetEntries() < 10) continue;
                TString lbl = Form("%s  %s  %s", chanLabel,
                    tune27::etaLabel(e).Data(), tune27::ptLabel(p).Data());
                drawPanel(&c, v + 1,
                          data.h[v][e][p], unf.h[v][e][p], fud.h[v][e][p],
                          kVars[v].label, lbl.Data(), v == 0);
            }
            c.Print(outFile);
        }
    }

    c.Print(Form("%s]", outFile));
    std::cout << "  Saved: " << outFile << std::endl;
}

// ===================================================================
// Main
// ===================================================================
void plot_fudge_binned(Long64_t maxEv = -1, int isoLevel = 0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const char* isoTags[] = {"no_iso", "iso_tight", "iso_track", "iso_full"};
    const char* isoTag = (isoLevel >= 0 && isoLevel <= 3) ? isoTags[isoLevel] : "no_iso";
    std::cout << "Isolation level: " << isoLevel << " (" << isoTag << ")" << std::endl;

    const char* zeeg_data   = "/dcache/atlas/mfernand/qt_ntuples/data24/egam3.root";
    const char* zeeg_mc     = "/dcache/atlas/mfernand/qt_ntuples/data24/mc_eegamma.root";
    const char* zmumug_data = "/dcache/atlas/mfernand/qt_ntuples/data24/egam4.root";
    const char* zmumug_mc   = "/dcache/atlas/mfernand/qt_ntuples/data24/mc_mumugamma.root";

    TString outBase = Form("../../output/data_mc_reweighting_%s/", isoTag);
    gSystem->mkdir(outBase + "Zeeg", true);
    gSystem->mkdir(outBase + "Zmumug", true);
    gSystem->mkdir(outBase + "combined", true);

    // ---- Z→eeγ ----
    std::cout << "\n========== Z->eeg ==========\n" << std::endl;
    Histos zee_d, zee_u, zee_f;
    zee_d.create("zd"); zee_u.create("zu"); zee_f.create("zf");
    fillData(zeeg_data, zee_d, maxEv);
    fillMC(zeeg_mc, zee_u, zee_f, maxEv);

    plotInclusive(zee_d, zee_u, zee_f,
        outBase + "Zeeg/fudge_inclusive.pdf", "Z#rightarrowee#gamma");
    plotEtaBins(zee_d, zee_u, zee_f,
        outBase + "Zeeg/fudge_etabins.pdf", "Z#rightarrowee#gamma");
    plotEtaPtBins(zee_d, zee_u, zee_f,
        outBase + "Zeeg/fudge_eta_pt.pdf", "Z#rightarrowee#gamma");

    // ---- Z→μμγ ----
    std::cout << "\n========== Z->mumug ==========\n" << std::endl;
    Histos mmu_d, mmu_u, mmu_f;
    mmu_d.create("md"); mmu_u.create("mu_"); mmu_f.create("mf");
    fillData(zmumug_data, mmu_d, maxEv);
    fillMC(zmumug_mc, mmu_u, mmu_f, maxEv);

    plotInclusive(mmu_d, mmu_u, mmu_f,
        outBase + "Zmumug/fudge_inclusive.pdf", "Z#rightarrow#mu#mu#gamma");
    plotEtaBins(mmu_d, mmu_u, mmu_f,
        outBase + "Zmumug/fudge_etabins.pdf", "Z#rightarrow#mu#mu#gamma");
    plotEtaPtBins(mmu_d, mmu_u, mmu_f,
        outBase + "Zmumug/fudge_eta_pt.pdf", "Z#rightarrow#mu#mu#gamma");

    // ---- Combined ----
    std::cout << "\n========== Combined ==========\n" << std::endl;
    Histos comb_d, comb_u, comb_f;
    comb_d.create("cd"); comb_u.create("cu"); comb_f.create("cf");
    comb_d.add(zee_d); comb_d.add(mmu_d);
    comb_u.add(zee_u); comb_u.add(mmu_u);
    comb_f.add(zee_f); comb_f.add(mmu_f);

    plotInclusive(comb_d, comb_u, comb_f,
        outBase + "combined/fudge_inclusive.pdf", "Z#rightarrowll#gamma");
    plotEtaBins(comb_d, comb_u, comb_f,
        outBase + "combined/fudge_etabins.pdf", "Z#rightarrowll#gamma");
    plotEtaPtBins(comb_d, comb_u, comb_f,
        outBase + "combined/fudge_eta_pt.pdf", "Z#rightarrowll#gamma");

    // ---- Update old fudge_comparison_*.pdf (side-by-side Zeeg | Zmumug) ----
    std::cout << "\n========== Updating fudge_comparison PDFs ==========\n" << std::endl;
    for (int v = 0; v < nVar; ++v) {
        TCanvas c2("c2", "", 1400, 600);
        c2.Divide(2, 1);
        {
            TH1D* hd = zee_d.sumAll(v);
            TH1D* hu = zee_u.sumAll(v);
            TH1D* hf = zee_f.sumAll(v);
            drawPanel(&c2, 1, hd, hu, hf, kVars[v].label,
                      "Z#rightarrowee#gamma  Inclusive", true);
        }
        {
            TH1D* hd = mmu_d.sumAll(v);
            TH1D* hu = mmu_u.sumAll(v);
            TH1D* hf = mmu_f.sumAll(v);
            drawPanel(&c2, 2, hd, hu, hf, kVars[v].label,
                      "Z#rightarrow#mu#mu#gamma  Inclusive", false);
        }
        TString out = Form("%sfudge_comparison_%s.pdf", outBase.Data(), kVars[v].name);
        c2.SaveAs(out);
        std::cout << "  Updated: " << out << std::endl;
    }

    std::cout << "\nAll done!" << std::endl;
}
