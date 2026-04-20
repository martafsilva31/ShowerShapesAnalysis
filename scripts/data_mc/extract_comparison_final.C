///////////////////////////////////////////////////////////////////////////////
// extract_comparison_final.C
//
// Compares chi-squared results across the 4 variants:
//   eta_loose, eta_tight, eta_pt_loose, eta_pt_tight
//
// Reads the chi2_summary.tex from each variant's report/ and produces a
// combined LaTeX comparison table in report/chi2_variant_comparison.tex.
//
// Also reads histograms.root directly if chi2 tex files are not yet produced.
//
// Usage:
//   cd scripts/data_mc
//   root -l -b -q 'extract_comparison_final.C()'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace config;

// ======================================================================
// Chi-squared / ndf between two area-normalised histograms.
// ======================================================================
double chi2ndf(TH1D* h1, TH1D* h2) {
    if (!h1 || !h2) return -1;
    double chi2 = 0;
    int ndf = 0;
    for (int b = 1; b <= h1->GetNbinsX(); ++b) {
        double v1 = h1->GetBinContent(b);
        double v2 = h2->GetBinContent(b);
        if (v2 > 0) {
            chi2 += (v1 - v2) * (v1 - v2) / v2;
            ndf++;
        }
    }
    return (ndf > 1) ? chi2 / ndf : -1;
}

TH1D* normClone(TH1D* h, const char* suffix) {
    if (!h || h->GetEntries() < 50) return nullptr;
    TH1D* hc = (TH1D*)h->Clone(Form("%s_%s", h->GetName(), suffix));
    hc->SetDirectory(0);
    double s = hc->Integral();
    if (s > 0) hc->Scale(1.0 / s);
    return hc;
}

static std::string fmtChi2(double v) {
    if (v < 0) return "---";
    char buf[32];
    snprintf(buf, sizeof(buf), "%.4f", v);
    return buf;
}

// ======================================================================
// Main
// ======================================================================
void extract_comparison_final(
    const char* outBase   = "../../output/Layer_2",
    const char* reportDir = "../../report")
{
    const int nVar = 4;
    const char* variants[nVar] = {"eta_loose", "eta_tight",
                                  "eta_pt_loose", "eta_pt_tight"};
    const char* variantTeX[nVar] = {
        "$\\eta$ loose", "$\\eta$ tight",
        "$\\eta \\times p_{\\mathrm{T}}$ loose",
        "$\\eta \\times p_{\\mathrm{T}}$ tight"};

    const int nSS = 3;
    const char* ssNames[nSS]  = {"reta", "rphi", "weta2"};
    const char* ssTeX[nSS]    = {"$R_{\\eta}$", "$R_{\\phi}$", "$w_{\\eta 2}$"};

    const int nMeth = 3;  // MC, M1, M2
    const char* methTeX[nMeth] = {"MC", "M1", "M2"};

    const int nScen = 3;
    const char* scenarios[nScen] = {"unconverted", "converted", "inclusive"};
    const char* scenTeX[nScen]   = {"Unconverted", "Converted", "Inclusive"};

    const char* channel = "llgamma";

    // Result storage: avg chi2 per [variant][scenario][SS][method]
    double avg[nVar][nScen][nSS][nMeth];
    for (auto& a : avg)
        for (auto& b : a)
            for (auto& c : b)
                for (auto& d : c) d = -1;

    // ============================================================
    // Extract chi2 from histograms.root for each variant/scenario
    // ============================================================
    for (int iv = 0; iv < nVar; ++iv) {
        for (int is = 0; is < nScen; ++is) {
            TString path = Form("%s/%s/%s/%s/histograms.root",
                                outBase, variants[iv], channel, scenarios[is]);
            TFile* f = TFile::Open(path, "READ");
            if (!f || f->IsZombie()) {
                std::cerr << "WARNING: Cannot open " << path << std::endl;
                if (f) delete f;
                continue;
            }

            for (int iss = 0; iss < nSS; ++iss) {
                double sum[nMeth] = {};
                int cnt[nMeth] = {};

                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    TString suf = Form("_eta%02d", ie);
                    TH1D* hd  = (TH1D*)f->Get(
                        Form("h_%s_data%s", ssNames[iss], suf.Data()));
                    TH1D* hmc = (TH1D*)f->Get(
                        Form("h_%s_mc%s", ssNames[iss], suf.Data()));
                    TH1D* hm1 = (TH1D*)f->Get(
                        Form("h_%s_mc_M1%s", ssNames[iss], suf.Data()));
                    TH1D* hm2 = (TH1D*)f->Get(
                        Form("h_%s_mc_M2%s", ssNames[iss], suf.Data()));

                    if (!hd || hd->GetEntries() < 100) continue;

                    TH1D* nd  = normClone(hd,  "nd");
                    TH1D* nmc = normClone(hmc, "nmc");
                    TH1D* nm1 = normClone(hm1, "nm1");
                    TH1D* nm2 = normClone(hm2, "nm2");

                    double v[3];
                    v[0] = (nd && nmc) ? chi2ndf(nmc, nd) : -1;
                    v[1] = (nd && nm1) ? chi2ndf(nm1, nd) : -1;
                    v[2] = (nd && nm2) ? chi2ndf(nm2, nd) : -1;

                    for (int im = 0; im < nMeth; ++im) {
                        if (v[im] >= 0) { sum[im] += v[im]; cnt[im]++; }
                    }

                    delete nd; delete nmc; delete nm1; delete nm2;
                }

                for (int im = 0; im < nMeth; ++im)
                    if (cnt[im] > 0) avg[iv][is][iss][im] = sum[im] / cnt[im];
            }

            f->Close();
            delete f;
        }
    }

    // ============================================================
    // Print summary to stdout
    // ============================================================
    std::cout << "\n========== CROSS-VARIANT CHI2/NDF COMPARISON ==========\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int iv = 0; iv < nVar; ++iv) {
        std::cout << "\n--- " << variants[iv] << " ---\n";
        for (int is = 0; is < nScen; ++is) {
            std::cout << "  " << scenarios[is] << ":";
            for (int iss = 0; iss < nSS; ++iss) {
                std::cout << "  " << ssNames[iss] << "(";
                for (int im = 0; im < nMeth; ++im) {
                    if (im > 0) std::cout << "/";
                    std::cout << fmtChi2(avg[iv][is][iss][im]);
                }
                std::cout << ")";
            }
            std::cout << "\n";
        }
    }

    // ============================================================
    // Write LaTeX comparison table
    // ============================================================
    TString fn = Form("%s/chi2_variant_comparison.tex", reportDir);
    std::ofstream out(fn.Data());
    out << "% Auto-generated by extract_comparison_final.C\n";
    out << "\\begin{table}[htbp]\n";
    out << "\\centering\n";
    out << "\\caption{Average $\\chisq$ comparison across binning and isolation variants.\n";
    out << "  Each cell shows $\\langle\\chisq\\rangle$ averaged over all populated $|\\eta|$ bins\n";
    out << "  for the $Z\\to\\ell\\ell\\gamma$ channel.  Lower is better.}\n";
    out << "\\label{tab:chi2-variant-comparison}\n";
    out << "\\resizebox{\\textwidth}{!}{%\n";
    out << "\\begin{tabular}{ll";
    for (int iss = 0; iss < nSS; ++iss) out << "|rrr";
    out << "}\n";
    out << "\\toprule\n";
    out << " & ";
    for (int iss = 0; iss < nSS; ++iss)
        out << "& \\multicolumn{3}{"
            << (iss < nSS - 1 ? "c|" : "c")
            << "}{" << ssTeX[iss] << "} ";
    out << "\\\\\n";
    out << "Variant & Scenario ";
    for (int iss = 0; iss < nSS; ++iss)
        for (int im = 0; im < nMeth; ++im)
            out << "& " << methTeX[im] << " ";
    out << "\\\\\n";
    out << "\\midrule\n";

    for (int iv = 0; iv < nVar; ++iv) {
        for (int is = 0; is < nScen; ++is) {
            out << variantTeX[iv] << " & " << scenTeX[is];
            for (int iss = 0; iss < nSS; ++iss)
                for (int im = 0; im < nMeth; ++im)
                    out << " & " << fmtChi2(avg[iv][is][iss][im]);
            out << " \\\\\n";
        }
        if (iv < nVar - 1) out << "\\midrule\n";
    }

    out << "\\bottomrule\n";
    out << "\\end{tabular}}\n";
    out << "\\end{table}\n";
    out.close();
    std::cout << "\nWritten: " << fn << std::endl;
}
