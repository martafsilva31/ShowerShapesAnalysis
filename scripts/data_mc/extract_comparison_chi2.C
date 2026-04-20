///////////////////////////////////////////////////////////////////////////////
// extract_comparison_chi2.C
//
// Cross-source chi-squared comparison table for llgamma unconverted:
//   Francisco/run3  — 3 cols: MC | Corrected | Fudge
//   hist_sigma_wide — 4 cols: MC | M1 | M2 | Fudge
//   hist_sigma_cellrange — 4 cols: MC | M1 | M2 | Fudge
//
// Note: Francisco's histograms.root was produced by convert_histograms.py
// which maps their single cell reweighting to M2 (and clones it as M1).
// Therefore M1 is NOT shown for Francisco — it is not an independent method.
//
// Produces:
//   {reportDir}/chi2_comparison.tex         — per-eta comparison table
//   {reportDir}/chi2_comparison_summary.tex — averaged comparison table
//
// Usage:
//   cd scripts/data_mc
//   root -l -b -q 'extract_comparison_chi2.C()'
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

// ── Utilities ──────────────────────────────────────────────────────────────

double chi2ndf(TH1D* h1, TH1D* h2) {
    if (!h1 || !h2) return -1;
    double chi2 = 0; int ndf = 0;
    for (int b = 1; b <= h1->GetNbinsX(); ++b) {
        double v1 = h1->GetBinContent(b);
        double v2 = h2->GetBinContent(b);
        if (v2 > 0) { chi2 += (v1 - v2) * (v1 - v2) / v2; ndf++; }
    }
    return (ndf > 1) ? chi2 / ndf : -1;
}

TH1D* normClone(TH1D* h, const char* sfx) {
    if (!h || h->GetEntries() < 50) return nullptr;
    TH1D* hc = (TH1D*)h->Clone(Form("%s_%s", h->GetName(), sfx));
    hc->SetDirectory(0);
    double s = hc->Integral();
    if (s > 0) hc->Scale(1.0 / s);
    return hc;
}

static std::string fmt(double v) {
    if (v < 0) return "---";
    char buf[32]; snprintf(buf, sizeof(buf), "%.4f", v); return buf;
}

// ── Per-eta result struct ───────────────────────────────────────────────────

struct EtaResult {
    double mc, m1, m2, fudge;
    bool valid = false;
};

// ── Extract results from one histograms.root ────────────────────────────────

void extractFromFile(
    TFile* f, int nEtaBins, EtaResult res[][3])  // res[etaBin][var]
{
    const char* vars[3] = {"reta", "rphi", "weta2"};

    for (int ie = 0; ie < nEtaBins; ++ie) {
        TString suf = Form("_eta%02d", ie);

        for (int iv = 0; iv < 3; ++iv) {
            TH1D* hData = (TH1D*)f->Get(Form("h_%s_data%s",      vars[iv], suf.Data()));
            TH1D* hMC   = (TH1D*)f->Get(Form("h_%s_mc%s",        vars[iv], suf.Data()));
            TH1D* hM1   = (TH1D*)f->Get(Form("h_%s_mc_M1%s",     vars[iv], suf.Data()));
            TH1D* hM2   = (TH1D*)f->Get(Form("h_%s_mc_M2%s",     vars[iv], suf.Data()));
            TH1D* hFudg = (TH1D*)f->Get(Form("h_%s_mc_fudged%s", vars[iv], suf.Data()));

            TH1D* nd  = normClone(hData, "nd");
            TH1D* nmc = normClone(hMC,   "nmc");
            TH1D* nm1 = normClone(hM1,   "nm1");
            TH1D* nm2 = normClone(hM2,   "nm2");
            TH1D* nfg = normClone(hFudg, "nfg");

            if (!nd || (!nmc && !nm1 && !nm2)) continue;

            res[ie][iv].mc    = chi2ndf(nmc, nd);
            res[ie][iv].m1    = chi2ndf(nm1, nd);
            res[ie][iv].m2    = chi2ndf(nm2, nd);
            res[ie][iv].fudge = chi2ndf(nfg, nd);
            res[ie][iv].valid = true;

            delete nd; delete nmc; delete nm1; delete nm2; delete nfg;
        }
    }
}

// ── Main ───────────────────────────────────────────────────────────────────

void extract_comparison_chi2(
    const char* franciscoFile = "/project/atlas/users/mfernand/QT/ShowerShapes/photoncellbasedrw_run3/output/histograms.root",
    const char* wideDir       = "../../output/Layer_2/hist_sigma_wide",
    const char* cellrangeDir  = "../../output/Layer_2/hist_sigma_cellrange",
    const char* reportDir     = "../../output/Layer_2/comparison_report",
    const char* scenario      = "unconverted",
    const char* channel       = "llgamma")
{
    gSystem->mkdir(reportDir, true);

    std::cout << "Processing extract_comparison_chi2.C\n"
              << "  scenario = " << scenario << "\n"
              << "  channel  = " << channel  << "\n";

    // ── Open files ──────────────────────────────────────────────────────────
    TFile* fFranc = TFile::Open(franciscoFile, "READ");
    TFile* fWide  = TFile::Open(Form("%s/%s/%s/histograms.root", wideDir,       channel, scenario), "READ");
    TFile* fCell  = TFile::Open(Form("%s/%s/%s/histograms.root", cellrangeDir,  channel, scenario), "READ");

    struct Source { TFile* f; const char* label; bool isFrancisco; };
    Source sources[3] = {
        {fFranc, "Francisco / run3", true },
        {fWide,  "Wide $[-1,2]$",    false},
        {fCell,  "Cell-dep.",        false},
    };

    const int nSrc  = 3;
    const int kNE   = kNEtaBins;
    const int kCrack = 8;   // [1.37, 1.52) — skip

    // res[src][eta][var]
    EtaResult res[nSrc][14][3];
    for (int s = 0; s < nSrc; ++s) {
        if (!sources[s].f || sources[s].f->IsZombie()) {
            std::cerr << "WARNING: Cannot open file for source " << sources[s].label << "\n";
            continue;
        }
        extractFromFile(sources[s].f, kNE, res[s]);
    }

    const char* etaTeX[14] = {
        "$[0.00, 0.20)$", "$[0.20, 0.40)$", "$[0.40, 0.60)$", "$[0.60, 0.80)$",
        "$[0.80, 1.00)$", "$[1.00, 1.20)$", "$[1.20, 1.30)$", "$[1.30, 1.37)$",
        "\\textit{crack}",
        "$[1.52, 1.60)$", "$[1.60, 1.80)$", "$[1.80, 2.00)$", "$[2.00, 2.20)$",
        "$[2.20, 2.40)$"
    };

    const char* vars[3]   = {"reta", "rphi", "weta2"};
    const char* varTeX[3] = {"$R_{\\eta}$", "$R_{\\phi}$", "$w_{\\eta 2}$"};

    // Helper: print one source's columns into a row
    // For Francisco (isFrancisco=true): MC | Corr. | Fudge  (3 values)
    // For others:                       MC | M1    | M2 | Fudge  (4 values)
    auto printRowCols = [&](std::ostream& os, int s, int ie, int iv) {
        auto& r = res[s][ie][iv];
        if (sources[s].isFrancisco) {
            os << " & " << fmt(r.mc)
               << " & " << fmt(r.m2)    // their single correction
               << " & " << fmt(r.fudge);
        } else {
            os << " & " << fmt(r.mc)
               << " & " << fmt(r.m1)
               << " & " << fmt(r.m2)
               << " & " << fmt(r.fudge);
        }
    };

    // ── Write per-var, per-eta tables ────────────────────────────────────────
    // Column layout: l | rrr | rrrr | rrrr   (1+3+4+4 = 12 cols)
    TString outFile = Form("%s/chi2_comparison.tex", reportDir);
    std::ofstream out(outFile.Data());
    out << "% Auto-generated by extract_comparison_chi2.C -- do not edit by hand\n";
    out << "% Scenario: " << scenario << "  Channel: " << channel << "\n\n";

    for (int iv = 0; iv < 3; ++iv) {
        out << "\\begin{table}[htbp]\n\\centering\n";
        out << "\\caption{$\\chi^2/\\text{ndf}$ comparison for " << varTeX[iv]
            << " (" << scenario << ", " << channel << ")."
               " Francisco has a single cell-based correction (shown as Corr.);$"
               " M1/M2 are defined for Wide and Cell-dep.\\ only.}\n";
        out << "\\label{tab:chi2-cmp-" << vars[iv] << "-" << scenario << "}\n";
        out << "\\small\n";
        out << "\\begin{tabular}{l|rrr|rrrr|rrrr}\n\\toprule\n";
        out << " & \\multicolumn{3}{c|}{Francisco / run3}"
               " & \\multicolumn{4}{c|}{Wide $[-1,2]$}"
               " & \\multicolumn{4}{c}{Cell-dep.} \\\\\n";
        out << "$|\\eta|$ & MC & Corr. & Fudge"
               " & MC & M1 & M2 & Fudge"
               " & MC & M1 & M2 & Fudge \\\\\n";
        out << "\\midrule\n";

        // sums: [src][col]  col: 0=mc,1=m1(or corr),2=m2,3=fudge
        double sumV[nSrc][4] = {};
        int    cntV[nSrc]    = {};

        for (int ie = 0; ie < kNE; ++ie) {
            if (ie == kCrack) continue;
            out << etaTeX[ie];
            for (int s = 0; s < nSrc; ++s) {
                printRowCols(out, s, ie, iv);
                auto& r = res[s][ie][iv];
                if (r.valid) {
                    sumV[s][0] += r.mc;
                    sumV[s][1] += r.m1;  // for Francisco, m1==m2 but we accumulate m2 below
                    sumV[s][2] += r.m2;
                    sumV[s][3] += r.fudge;
                    cntV[s]++;
                }
            }
            out << " \\\\\n";
        }

        out << "\\midrule\n\\textbf{Average}";
        for (int s = 0; s < nSrc; ++s) {
            int c = cntV[s];
            if (sources[s].isFrancisco) {
                out << " & \\textbf{" << fmt(c ? sumV[s][0]/c : -1) << "}"
                    << " & \\textbf{" << fmt(c ? sumV[s][2]/c : -1) << "}"  // corr = m2
                    << " & \\textbf{" << fmt(c ? sumV[s][3]/c : -1) << "}";
            } else {
                out << " & \\textbf{" << fmt(c ? sumV[s][0]/c : -1) << "}"
                    << " & \\textbf{" << fmt(c ? sumV[s][1]/c : -1) << "}"
                    << " & \\textbf{" << fmt(c ? sumV[s][2]/c : -1) << "}"
                    << " & \\textbf{" << fmt(c ? sumV[s][3]/c : -1) << "}";
            }
        }
        out << " \\\\\n";
        out << "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n";
    }
    out.close();
    std::cout << "Written: " << outFile << "\n";

    // ── Write summary table ──────────────────────────────────────────────────
    // Cols per variable: Francisco (3) | Wide (4) | Cell (4) = 11 + label = 12
    TString sumFile = Form("%s/chi2_comparison_summary.tex", reportDir);
    std::ofstream sout(sumFile.Data());
    sout << "% Auto-generated by extract_comparison_chi2.C -- do not edit by hand\n";
    sout << "% Scenario: " << scenario << "  Channel: " << channel << "\n\n";
    sout << "\\begin{table}[htbp]\n\\centering\n";
    sout << "\\caption{Average $\\chi^2/\\text{ndf}$ across $|\\eta|$ bins ("
         << scenario << ", " << channel << ")."
            " Francisco has a single cell-based correction (Corr.); M1/M2 apply to Wide/Cell-dep.\\ only."
            " Lower is better.}\n";
    sout << "\\label{tab:chi2-cmp-summary}\n";
    sout << "\\small\n";
    // Per variable: Francisco has 3 cols, others have 4 cols
    // Variables: 3.  Total data cols = 3*(3+4+4) = 33.  Plus label col.
    // We'll group per variable block:
    //   Reta: [Franc: MC,Corr,Fudge | Wide: MC,M1,M2,Fudge | Cell: MC,M1,M2,Fudge]
    sout << "\\begin{tabular}{l|rrr|rrrr|rrrr|rrr|rrrr|rrrr|rrr|rrrr|rrrr}\n\\toprule\n";

    // Variable-level header
    sout << " & \\multicolumn{11}{c|}{$R_{\\eta}$}"
            " & \\multicolumn{11}{c|}{$R_{\\phi}$}"
            " & \\multicolumn{11}{c}{$w_{\\eta 2}$} \\\\\n";
    // Source-level header (within each variable block)
    std::string srcHdr3 = "\\multicolumn{3}{c|}{Franc.}";
    std::string srcHdr4w = "\\multicolumn{4}{c|}{Wide}";
    std::string srcHdr4c = "\\multicolumn{4}{c}{Cell-dep.}";
    std::string srcHdr4c_mid = "\\multicolumn{4}{c|}{Cell-dep.}";
    for (int vv = 0; vv < 3; ++vv) {
        sout << " &";
        sout << srcHdr3 << " & " << srcHdr4w << " & ";
        if (vv < 2) sout << srcHdr4c_mid;
        else        sout << srcHdr4c;
    }
    sout << " \\\\\n";
    // Column header (MC, Corr./M1/M2, Fudge)
    std::string colsFranc = "MC & Corr. & Fudge";
    std::string colsOther = "MC & M1 & M2 & Fudge";
    sout << "Source";
    for (int vv = 0; vv < 3; ++vv) {
        sout << " & " << colsFranc << " & " << colsOther << " & " << colsOther;
    }
    sout << " \\\\\n\\midrule\n";

    const char* srcLab[3] = {"Francisco / run3", "Wide $[-1,2]$", "Cell-dep."};

    // Compute averages per source, per var
    // avg[src][var][col]: col 0=mc, 1=m1, 2=m2, 3=fudge
    double avg[nSrc][3][4] = {};
    int    cnt[nSrc][3]    = {};
    for (int s = 0; s < nSrc; ++s) {
        for (int ie = 0; ie < kNE; ++ie) {
            if (ie == kCrack) continue;
            for (int iv = 0; iv < 3; ++iv) {
                auto& r = res[s][ie][iv];
                if (!r.valid) continue;
                avg[s][iv][0] += r.mc;
                avg[s][iv][1] += r.m1;
                avg[s][iv][2] += r.m2;
                avg[s][iv][3] += r.fudge;
                cnt[s][iv]++;
            }
        }
    }

    // Print each source as its own row
    for (int s = 0; s < nSrc; ++s) {
        sout << srcLab[s];
        for (int iv = 0; iv < 3; ++iv) {
            int c = cnt[s][iv];
            if (sources[s].isFrancisco) {
                sout << " & " << fmt(c ? avg[s][iv][0]/c : -1)   // MC
                     << " & " << fmt(c ? avg[s][iv][2]/c : -1)   // Corr. (=M2)
                     << " & " << fmt(c ? avg[s][iv][3]/c : -1);  // Fudge
            } else {
                sout << " & " << fmt(c ? avg[s][iv][0]/c : -1)   // MC
                     << " & " << fmt(c ? avg[s][iv][1]/c : -1)   // M1
                     << " & " << fmt(c ? avg[s][iv][2]/c : -1)   // M2
                     << " & " << fmt(c ? avg[s][iv][3]/c : -1);  // Fudge
            }
        }
        sout << " \\\\\n";
    }

    sout << "\\bottomrule\n\\end{tabular}\n\\end{table}\n";
    sout.close();
    std::cout << "Written: " << sumFile << "\n";

    if (fFranc) fFranc->Close();
    if (fWide)  fWide->Close();
    if (fCell)  fCell->Close();

    std::cout << "=== extract_comparison_chi2 complete ===\n";
}
