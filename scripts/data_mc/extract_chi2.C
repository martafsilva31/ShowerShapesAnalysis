///////////////////////////////////////////////////////////////////////////////
// extract_chi2.C
//
// Extracts chi-squared comparison tables from histograms.root output files
// for all 9 scenarios (3 channels x 3 conversion types) and writes LaTeX
// table fragments for inclusion in the physics report.
//
// Output (in report/ directory):
//   chi2_yields.tex   — Event yield table
//   chi2_summary.tex  — Average chi-squared summary table
//   chi2_tables.tex   — Per-eta chi-squared tables (all scenarios)
//
// Usage:
//   cd scripts/data_mc
//   root -l -b -q 'extract_chi2.C'
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
// h1 = test (MC variant), h2 = reference (data).
// Returns sum_b (h1_b - h2_b)^2 / h2_b  /  ndf.
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

// ======================================================================
// Area-normalise a cloned histogram.
// ======================================================================
TH1D* normClone(TH1D* h, const char* suffix) {
    if (!h || h->GetEntries() < 50) return nullptr;
    TH1D* hc = (TH1D*)h->Clone(Form("%s_%s", h->GetName(), suffix));
    hc->SetDirectory(0);
    double s = hc->Integral();
    if (s > 0) hc->Scale(1.0 / s);
    return hc;
}

// ======================================================================
// Format a chi2 value for LaTeX output.
// ======================================================================
static std::string fmtChi2(double v) {
    if (v < 0) return "---";
    char buf[32];
    snprintf(buf, sizeof(buf), "%.4f", v);
    return buf;
}

// ======================================================================
// Main
// ======================================================================
void extract_chi2(
    const char* baseDir   = "../../output/cell_energy_reweighting_Francisco_method/data24",
    const char* reportDir = "../../report")
{
    const int nChan = 3;
    const char* channels[nChan] = {"eegamma", "mumugamma", "llgamma"};
    const char* chanTeX[nChan]  = {
        "$Z \\to ee\\gamma$",
        "$Z \\to \\mu\\mu\\gamma$",
        "$Z \\to \\ell\\ell\\gamma$"};

    const int nScen = 3;
    const char* scenarios[nScen] = {"unconverted", "converted", "inclusive"};
    const char* scenTeX[nScen]   = {"Unconverted", "Converted", "Inclusive"};

    const int nVar = 3;
    const char* vars[nVar]    = {"reta", "rphi", "weta2"};
    const char* varTeX[nVar]  = {"$R_{\\eta}$", "$R_{\\phi}$", "$w_{\\eta 2}$"};

    const int nMeth = 4;  // MC, M1, M2, Fudged

    // ============================================================
    // Result storage
    // ============================================================
    struct EtaResult {
        double v[4];  // chi2/ndf for MC, M1, M2, Fudged
        bool valid;
    };
    struct ScenResult {
        EtaResult eta[15][3]; // [etaBin][var]
        double nData, nMC;
    };
    ScenResult R[3][3]; // [chan][scen]

    // Initialise
    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is) {
            R[ic][is].nData = R[ic][is].nMC = 0;
            for (int ie = 0; ie < kNEtaBins; ++ie)
                for (int iv = 0; iv < nVar; ++iv)
                    R[ic][is].eta[ie][iv] = {{-1, -1, -1, -1}, false};
        }

    // ============================================================
    // Extract chi-squared from all histograms.root files
    // ============================================================
    for (int ic = 0; ic < nChan; ++ic) {
        for (int is = 0; is < nScen; ++is) {
            TString path = Form("%s/%s/%s/histograms.root",
                                baseDir, channels[ic], scenarios[is]);
            TFile* f = TFile::Open(path, "READ");
            if (!f || f->IsZombie()) {
                std::cerr << "WARNING: Cannot open " << path << std::endl;
                if (f) delete f;
                continue;
            }
            std::cout << "Reading: " << path << std::endl;

            // Event yields: sum entries across eta bins
            for (int ie = 0; ie < kNEtaBins; ++ie) {
                TH1D* hd = (TH1D*)f->Get(Form("h_reta_data_eta%02d", ie));
                TH1D* hm = (TH1D*)f->Get(Form("h_reta_mc_eta%02d", ie));
                if (hd) R[ic][is].nData += hd->GetEntries();
                if (hm) R[ic][is].nMC   += hm->GetEntries();
            }

            // Chi-squared per eta bin per variable
            for (int iv = 0; iv < nVar; ++iv) {
                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    TString suf = Form("_eta%02d", ie);

                    // Set B: cell-computed
                    TH1D* hd  = (TH1D*)f->Get(
                        Form("h_%s_data%s", vars[iv], suf.Data()));
                    TH1D* hmc = (TH1D*)f->Get(
                        Form("h_%s_mc%s",   vars[iv], suf.Data()));
                    TH1D* hm1 = (TH1D*)f->Get(
                        Form("h_%s_mc_M1%s", vars[iv], suf.Data()));
                    TH1D* hm2 = (TH1D*)f->Get(
                        Form("h_%s_mc_M2%s", vars[iv], suf.Data()));

                    if (!hd || hd->GetEntries() < 100) continue;

                    TH1D* nd  = normClone(hd,  "nd");
                    TH1D* nmc = normClone(hmc, "nmc");
                    TH1D* nm1 = normClone(hm1, "nm1");
                    TH1D* nm2 = normClone(hm2, "nm2");

                    EtaResult& r = R[ic][is].eta[ie][iv];
                    r.v[0] = (nd && nmc) ? chi2ndf(nmc, nd) : -1;
                    r.v[1] = (nd && nm1) ? chi2ndf(nm1, nd) : -1;
                    r.v[2] = (nd && nm2) ? chi2ndf(nm2, nd) : -1;

                    // Fudge factor chi2: Data(stored) vs MC(fudged)
                    TH1D* hds = (TH1D*)f->Get(
                        Form("h_%s_data_stored%s", vars[iv], suf.Data()));
                    TH1D* hmf = (TH1D*)f->Get(
                        Form("h_%s_mc_fudged%s", vars[iv], suf.Data()));
                    TH1D* nds = normClone(hds, "nds");
                    TH1D* nmf = normClone(hmf, "nmf");
                    r.v[3] = (nds && nmf) ? chi2ndf(nmf, nds) : -1;
                    r.valid = (r.v[0] >= 0);

                    delete nd; delete nmc; delete nm1; delete nm2;
                    delete nds; delete nmf;
                }
            }
            f->Close();
            delete f;
        }
    }

    // ============================================================
    // Compute averages: avg[chan][scen][var][method]
    // ============================================================
    double avg[3][3][3][4] = {};
    int    cnt[3][3][3][4] = {};

    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is)
            for (int iv = 0; iv < nVar; ++iv)
                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    EtaResult& r = R[ic][is].eta[ie][iv];
                    if (!r.valid) continue;
                    for (int im = 0; im < nMeth; ++im)
                        if (r.v[im] >= 0) {
                            avg[ic][is][iv][im] += r.v[im];
                            cnt[ic][is][iv][im]++;
                        }
                }

    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is)
            for (int iv = 0; iv < nVar; ++iv)
                for (int im = 0; im < nMeth; ++im)
                    if (cnt[ic][is][iv][im] > 0)
                        avg[ic][is][iv][im] /= cnt[ic][is][iv][im];

    // ============================================================
    // Print summary to stdout
    // ============================================================
    std::cout << "\n========== AVERAGE CHI2/NDF SUMMARY ==========\n";
    std::cout << std::fixed << std::setprecision(3);
    for (int ic = 0; ic < nChan; ++ic) {
        for (int is = 0; is < nScen; ++is) {
            std::cout << channels[ic] << " / " << scenarios[is]
                      << "  (Data=" << (int)R[ic][is].nData
                      << ", MC=" << (int)R[ic][is].nMC << ")\n";
            for (int iv = 0; iv < nVar; ++iv)
                std::cout << "  " << vars[iv]
                          << ":  MC=" << avg[ic][is][iv][0]
                          << "  M1=" << avg[ic][is][iv][1]
                          << "  M2=" << avg[ic][is][iv][2]
                          << "  Fudged=" << avg[ic][is][iv][3] << "\n";
        }
    }

    // ============================================================
    // Write chi2_yields.tex
    // ============================================================
    {
        TString fn = Form("%s/chi2_yields.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2.C — do not edit by hand\n";
        out << "\\begin{table}[htbp]\n";
        out << "\\centering\n";
        out << "\\caption{Event yields after all selection cuts for each channel "
               "and conversion scenario.  Data counts are unweighted; MC counts "
               "are unweighted (the MC event weight $w = w_{\\text{MC}} \\times "
               "w_\\mu \\times \\sigma$ is applied to histograms but not to the "
               "event count).}\n";
        out << "\\label{tab:yields}\n";
        out << "\\begin{tabular}{llrr}\n";
        out << "\\toprule\n";
        out << "Channel & Scenario & Data events & MC events \\\\\n";
        out << "\\midrule\n";
        for (int ic = 0; ic < nChan; ++ic) {
            for (int is = 0; is < nScen; ++is) {
                out << chanTeX[ic] << " & " << scenTeX[is]
                    << " & " << (int)R[ic][is].nData
                    << " & " << (int)R[ic][is].nMC
                    << " \\\\\n";
            }
            if (ic < nChan - 1) out << "\\midrule\n";
        }
        out << "\\bottomrule\n";
        out << "\\end{tabular}\n";
        out << "\\end{table}\n";
        out.close();
        std::cout << "\nWritten: " << fn << std::endl;
    }

    // ============================================================
    // Write chi2_summary.tex — average chi2 across eta bins
    // ============================================================
    {
        TString fn = Form("%s/chi2_summary.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2.C — do not edit by hand\n";
        out << "\\begin{table}[htbp]\n";
        out << "\\centering\n";
        out << "\\caption{Average $\\chisq$ across all non-empty $|\\eta|$ bins "
               "for each correction method.  Lower is better.}\n";
        out << "\\label{tab:chi2-summary}\n";
        out << "\\resizebox{\\textwidth}{!}{%\n";
        out << "\\begin{tabular}{ll";
        for (int iv = 0; iv < nVar; ++iv) out << "|rrrr";
        out << "}\n";
        out << "\\toprule\n";
        out << " & ";
        for (int iv = 0; iv < nVar; ++iv)
            out << "& \\multicolumn{4}{"
                << (iv < nVar - 1 ? "c|" : "c")
                << "}{" << varTeX[iv] << "} ";
        out << "\\\\\n";
        out << "Channel & Scenario ";
        for (int iv = 0; iv < nVar; ++iv)
            out << "& MC & M1 & M2 & Fudged ";
        out << "\\\\\n";
        out << "\\midrule\n";
        for (int ic = 0; ic < nChan; ++ic) {
            for (int is = 0; is < nScen; ++is) {
                out << chanTeX[ic] << " & " << scenTeX[is];
                for (int iv = 0; iv < nVar; ++iv)
                    for (int im = 0; im < nMeth; ++im)
                        out << " & " << fmtChi2(avg[ic][is][iv][im]);
                out << " \\\\\n";
            }
            if (ic < nChan - 1) out << "\\midrule\n";
        }
        out << "\\bottomrule\n";
        out << "\\end{tabular}}\n";
        out << "\\end{table}\n";
        out.close();
        std::cout << "Written: " << fn << std::endl;
    }

    // ============================================================
    // Write chi2_tables.tex — per-eta tables for all scenarios
    // ============================================================
    {
        TString fn = Form("%s/chi2_tables.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2.C — do not edit by hand\n";
        out << "% Per-$|\\eta|$-bin $\\chi^2/n_{\\mathrm{df}}$ tables "
               "for all channel/scenario combinations.\n\n";

        for (int ic = 0; ic < nChan; ++ic) {
            for (int is = 0; is < nScen; ++is) {
                TString label = Form("tab:chi2-%s-%s",
                                     channels[ic], scenarios[is]);

                out << "\\begin{table}[htbp]\n";
                out << "\\centering\n";
                out << "\\caption{Per-$|\\eta|$-bin $\\chisq$ for "
                    << chanTeX[ic] << ", " << scenTeX[is]
                    << " scenario.  Columns show $\\chisq$ between "
                       "each MC variant and data (cell-computed).}\n";
                out << "\\label{" << label << "}\n";
                out << "\\resizebox{\\textwidth}{!}{%\n";
                out << "\\begin{tabular}{r|rrrr|rrrr|rrrr}\n";
                out << "\\toprule\n";
                out << "$|\\eta|$ bin ";
                for (int iv = 0; iv < nVar; ++iv)
                    out << "& \\multicolumn{4}{"
                        << (iv < nVar - 1 ? "c|" : "c")
                        << "}{" << varTeX[iv] << "} ";
                out << "\\\\\n";
                out << "(range) ";
                for (int iv = 0; iv < nVar; ++iv)
                    out << "& MC & M1 & M2 & Fud ";
                out << "\\\\\n";
                out << "\\midrule\n";

                // Per-eta rows
                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    if (!R[ic][is].eta[ie][0].valid) continue;
                    out << Form("{[}%.2f, %.2f)",
                                kEtaLimits[ie], kEtaLimits[ie + 1]);
                    for (int iv = 0; iv < nVar; ++iv)
                        for (int im = 0; im < 4; ++im)  // MC, M1, M2, Fud
                            out << " & "
                                << fmtChi2(R[ic][is].eta[ie][iv].v[im]);
                    out << " \\\\\n";
                }

                // Average row
                out << "\\midrule\n";
                out << "Average";
                for (int iv = 0; iv < nVar; ++iv)
                    for (int im = 0; im < 4; ++im)
                        out << " & " << fmtChi2(avg[ic][is][iv][im]);
                out << " \\\\\n";
                out << "\\bottomrule\n";
                out << "\\end{tabular}}\n";
                out << "\\end{table}\n\n";
            }
        }

        out.close();
        std::cout << "Written: " << fn << std::endl;
    }

    // ============================================================
    // Compute INCLUSIVE chi2 from integrated histograms
    // ============================================================
    double incl[3][3][3][4] = {};  // [chan][scen][var][method]
    bool   inclValid[3][3] = {};

    for (int ic = 0; ic < nChan; ++ic) {
        for (int is = 0; is < nScen; ++is) {
            TString path = Form("%s/%s/%s/histograms.root",
                                baseDir, channels[ic], scenarios[is]);
            TFile* f = TFile::Open(path, "READ");
            if (!f || f->IsZombie()) { if (f) delete f; continue; }

            for (int iv = 0; iv < nVar; ++iv) {
                // Cell-computed inclusive
                TH1D* hd  = (TH1D*)f->Get(Form("h_%s_data_computed", vars[iv]));
                TH1D* hmc = (TH1D*)f->Get(Form("h_%s_mc_computed",   vars[iv]));
                TH1D* hm1 = (TH1D*)f->Get(Form("h_%s_mc_M1",        vars[iv]));
                TH1D* hm2 = (TH1D*)f->Get(Form("h_%s_mc_M2",        vars[iv]));

                if (!hd || hd->GetEntries() < 100) continue;
                TH1D* nd  = normClone(hd,  "ni_d");
                TH1D* nmc = normClone(hmc, "ni_mc");
                TH1D* nm1 = normClone(hm1, "ni_m1");
                TH1D* nm2 = normClone(hm2, "ni_m2");

                incl[ic][is][iv][0] = (nd && nmc) ? chi2ndf(nmc, nd) : -1;
                incl[ic][is][iv][1] = (nd && nm1) ? chi2ndf(nm1, nd) : -1;
                incl[ic][is][iv][2] = (nd && nm2) ? chi2ndf(nm2, nd) : -1;

                // Fudge: stored data vs fudged MC
                TH1D* hds = (TH1D*)f->Get(Form("h_%s_data_stored", vars[iv]));
                TH1D* hmf = (TH1D*)f->Get(Form("h_%s_mc_fudged",  vars[iv]));
                TH1D* nds = normClone(hds, "ni_ds");
                TH1D* nmf = normClone(hmf, "ni_mf");
                incl[ic][is][iv][3] = (nds && nmf) ? chi2ndf(nmf, nds) : -1;

                inclValid[ic][is] = true;
                delete nd; delete nmc; delete nm1; delete nm2;
                delete nds; delete nmf;
            }
            f->Close();
            delete f;
        }
    }

    // ============================================================
    // Write chi2_inclusive.tex
    // ============================================================
    {
        TString fn = Form("%s/chi2_inclusive.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2.C — do not edit by hand\n";
        out << "\\begin{table}[htbp]\n";
        out << "\\centering\n";
        out << "\\caption{Inclusive (all $|\\eta|$ bins combined) $\\chisq$ "
               "for each correction method, computed on integrated histograms.  "
               "Lower is better.}\n";
        out << "\\label{tab:chi2-inclusive}\n";
        out << "\\resizebox{\\textwidth}{!}{%\n";
        out << "\\begin{tabular}{ll";
        for (int iv = 0; iv < nVar; ++iv) out << "|rrrr";
        out << "}\n";
        out << "\\toprule\n";
        out << " & ";
        for (int iv = 0; iv < nVar; ++iv)
            out << "& \\multicolumn{4}{"
                << (iv < nVar - 1 ? "c|" : "c")
                << "}{" << varTeX[iv] << "} ";
        out << "\\\\\n";
        out << "Channel & Scenario ";
        for (int iv = 0; iv < nVar; ++iv)
            out << "& MC & M1 & M2 & Fudged ";
        out << "\\\\\n";
        out << "\\midrule\n";
        for (int ic = 0; ic < nChan; ++ic) {
            for (int is = 0; is < nScen; ++is) {
                if (!inclValid[ic][is]) continue;
                out << chanTeX[ic] << " & " << scenTeX[is];
                for (int iv = 0; iv < nVar; ++iv)
                    for (int im = 0; im < nMeth; ++im)
                        out << " & " << fmtChi2(incl[ic][is][iv][im]);
                out << " \\\\\n";
            }
            if (ic < nChan - 1) out << "\\midrule\n";
        }
        out << "\\bottomrule\n";
        out << "\\end{tabular}}\n";
        out << "\\end{table}\n";
        out.close();
        std::cout << "Written: " << fn << std::endl;
    }

    std::cout << "\n=== extract_chi2 complete ===" << std::endl;
}
