///////////////////////////////////////////////////////////////////////////////
// extract_chi2_layer1.C
//
// Layer-1 version of extract_chi2.C: 5 variables (weta1, wstot, fside,
// deltae, eratio) instead of 3.
///////////////////////////////////////////////////////////////////////////////

#include "config_layer1.h"

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace config_l1;

static double chi2ndf(TH1D* h1, TH1D* h2) {
    if (!h1 || !h2) return -1;
    double chi2 = 0; int ndf = 0;
    for (int b = 1; b <= h1->GetNbinsX(); ++b) {
        double v1 = h1->GetBinContent(b);
        double v2 = h2->GetBinContent(b);
        if (v2 > 0) { chi2 += (v1 - v2) * (v1 - v2) / v2; ndf++; }
    }
    return (ndf > 1) ? chi2 / ndf : -1;
}

static TH1D* normClone(TH1D* h, const char* suffix) {
    if (!h || h->GetEntries() < 50) return nullptr;
    TH1D* hc = (TH1D*)h->Clone(Form("%s_%s", h->GetName(), suffix));
    hc->SetDirectory(0);
    double s = hc->Integral();
    if (s > 0) hc->Scale(1.0 / s);
    return hc;
}

static std::string fmtChi2(double v) {
    if (v < 0) return "---";
    char buf[32]; snprintf(buf, sizeof(buf), "%.4f", v);
    return buf;
}

void extract_chi2_layer1(
    const char* baseDir   = "../../output/Layer_1/eta_loose",
    const char* reportDir = "../../output/Layer_1/eta_loose/report")
{
    const int nChan = 3;
    const char* channels[nChan] = {"eegamma", "mumugamma", "llgamma"};
    const char* chanTeX[nChan]  = {
        "$Z \\to ee\\gamma$", "$Z \\to \\mu\\mu\\gamma$", "$Z \\to \\ell\\ell\\gamma$"};

    const int nScen = 3;
    const char* scenarios[nScen] = {"unconverted", "converted", "inclusive"};
    const char* scenTeX[nScen]   = {"Unconverted", "Converted", "Conv.~\\& Unconv."};

    const int nVar = 5;
    const char* vars[nVar]    = {"weta1", "wstot", "fside", "deltae", "eratio"};
    const char* varTeX[nVar]  = {
        "$w_{\\eta 1}$", "$w_{s,\\mathrm{tot}}$", "$f_{\\mathrm{side}}$",
        "$\\Delta E$", "$E_{\\mathrm{ratio}}$"};
    const int nMeth = 6;  // MC, M1, M2, M3, M4, Fudged

    struct EtaResult { double v[6]; bool valid; };
    struct ScenResult { EtaResult eta[15][5]; double nData, nMC; };
    ScenResult R[3][3];

    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is) {
            R[ic][is].nData = R[ic][is].nMC = 0;
            for (int ie = 0; ie < kNEtaBins; ++ie)
                for (int iv = 0; iv < nVar; ++iv)
                    R[ic][is].eta[ie][iv] = {{-1, -1, -1, -1, -1, -1}, false};
        }

    // ===== Per-eta =====
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

            for (int ie = 0; ie < kNEtaBins; ++ie) {
                TH1D* hd = (TH1D*)f->Get(Form("h_weta1_data_eta%02d", ie));
                TH1D* hm = (TH1D*)f->Get(Form("h_weta1_mc_eta%02d",   ie));
                if (hd) R[ic][is].nData += hd->GetEntries();
                if (hm) R[ic][is].nMC   += hm->GetEntries();
            }

            for (int iv = 0; iv < nVar; ++iv) {
                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    TString suf = Form("_eta%02d", ie);
                    TH1D* hd  = (TH1D*)f->Get(Form("h_%s_data%s",  vars[iv], suf.Data()));
                    TH1D* hmc = (TH1D*)f->Get(Form("h_%s_mc%s",    vars[iv], suf.Data()));
                    TH1D* hm1 = (TH1D*)f->Get(Form("h_%s_mc_M1%s", vars[iv], suf.Data()));
                    TH1D* hm2 = (TH1D*)f->Get(Form("h_%s_mc_M2%s", vars[iv], suf.Data()));
                    TH1D* hm3 = (TH1D*)f->Get(Form("h_%s_mc_M3%s", vars[iv], suf.Data()));
                    TH1D* hm4 = (TH1D*)f->Get(Form("h_%s_mc_M4%s", vars[iv], suf.Data()));
                    if (!hd || hd->GetEntries() < 100) continue;

                    TH1D* nd  = normClone(hd,  "nd");
                    TH1D* nmc = normClone(hmc, "nmc");
                    TH1D* nm1 = normClone(hm1, "nm1");
                    TH1D* nm2 = normClone(hm2, "nm2");
                    TH1D* nm3 = normClone(hm3, "nm3");
                    TH1D* nm4 = normClone(hm4, "nm4");
                    EtaResult& r = R[ic][is].eta[ie][iv];
                    r.v[0] = (nd && nmc) ? chi2ndf(nmc, nd) : -1;
                    r.v[1] = (nd && nm1) ? chi2ndf(nm1, nd) : -1;
                    r.v[2] = (nd && nm2) ? chi2ndf(nm2, nd) : -1;
                    r.v[3] = (nd && nm3) ? chi2ndf(nm3, nd) : -1;
                    r.v[4] = (nd && nm4) ? chi2ndf(nm4, nd) : -1;

                    TH1D* hds = (TH1D*)f->Get(Form("h_%s_data_stored%s", vars[iv], suf.Data()));
                    TH1D* hmf = (TH1D*)f->Get(Form("h_%s_mc_fudged%s",   vars[iv], suf.Data()));
                    TH1D* nds = normClone(hds, "nds");
                    TH1D* nmf = normClone(hmf, "nmf");
                    r.v[5] = (nds && nmf) ? chi2ndf(nmf, nds) : -1;
                    r.valid = (r.v[0] >= 0);

                    delete nd; delete nmc; delete nm1; delete nm2;
                    delete nm3; delete nm4;
                    delete nds; delete nmf;
                }
            }
            f->Close(); delete f;
        }
    }

    // Averages
    double avg[3][3][5][6] = {}; int cnt[3][3][5][6] = {};
    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is)
            for (int iv = 0; iv < nVar; ++iv)
                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    EtaResult& r = R[ic][is].eta[ie][iv];
                    if (!r.valid) continue;
                    for (int im = 0; im < nMeth; ++im)
                        if (r.v[im] >= 0) { avg[ic][is][iv][im] += r.v[im]; cnt[ic][is][iv][im]++; }
                }
    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is)
            for (int iv = 0; iv < nVar; ++iv)
                for (int im = 0; im < nMeth; ++im)
                    if (cnt[ic][is][iv][im] > 0) avg[ic][is][iv][im] /= cnt[ic][is][iv][im];

    std::cout << "\n========== AVERAGE CHI2/NDF SUMMARY (Layer 1) ==========\n"
              << std::fixed << std::setprecision(3);
    for (int ic = 0; ic < nChan; ++ic)
        for (int is = 0; is < nScen; ++is) {
            std::cout << channels[ic] << " / " << scenarios[is]
                      << "  (Data=" << (int)R[ic][is].nData
                      << ", MC=" << (int)R[ic][is].nMC << ")\n";
            for (int iv = 0; iv < nVar; ++iv)
                std::cout << "  " << vars[iv]
                          << ":  MC=" << avg[ic][is][iv][0]
                          << "  M1=" << avg[ic][is][iv][1]
                          << "  M2=" << avg[ic][is][iv][2]
                          << "  M3=" << avg[ic][is][iv][3]
                          << "  M4=" << avg[ic][is][iv][4]
                          << "  Fudged=" << avg[ic][is][iv][5] << "\n";
        }

    // Ensure report dir exists
    gSystem->mkdir(reportDir, true);

    // chi2_yields.tex
    {
        TString fn = Form("%s/chi2_yields.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2_layer1.C — do not edit by hand\n"
            << "\\begin{table}[htbp]\n\\centering\n"
            << "\\caption{Layer-1 event yields after all selection cuts.}\n"
            << "\\label{tab:yields-l1}\n"
            << "\\begin{tabular}{llrr}\n\\toprule\n"
            << "Channel & Scenario & Data events & MC events \\\\\n\\midrule\n";
        for (int ic = 0; ic < nChan; ++ic) {
            for (int is = 0; is < nScen; ++is)
                out << chanTeX[ic] << " & " << scenTeX[is]
                    << " & " << (int)R[ic][is].nData
                    << " & " << (int)R[ic][is].nMC << " \\\\\n";
            if (ic < nChan - 1) out << "\\midrule\n";
        }
        out << "\\bottomrule\n\\end{tabular}\n\\end{table}\n";
        out.close();
        std::cout << "\nWritten: " << fn << std::endl;
    }

    // chi2_summary.tex
    {
        TString fn = Form("%s/chi2_summary.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2_layer1.C — do not edit by hand\n"
            << "\\begin{table}[htbp]\n\\centering\n"
            << "\\caption{Layer-1 average $\\chisq$ across all non-empty $|\\eta|$ bins.}\n"
            << "\\label{tab:chi2-summary-l1}\n"
            << "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{ll";
        for (int iv = 0; iv < nVar; ++iv) out << "|rrrrrr";
        out << "}\n\\toprule\n & ";
        for (int iv = 0; iv < nVar; ++iv)
            out << "& \\multicolumn{6}{" << (iv < nVar - 1 ? "c|" : "c")
                << "}{" << varTeX[iv] << "} ";
        out << "\\\\\nChannel & Scenario ";
        for (int iv = 0; iv < nVar; ++iv) out << "& MC & M1 & M2 & M3 & M4 & Fud ";
        out << "\\\\\n\\midrule\n";
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
        out << "\\bottomrule\n\\end{tabular}}\n\\end{table}\n";
        out.close();
        std::cout << "Written: " << fn << std::endl;
    }

    // chi2_tables.tex
    {
        TString fn = Form("%s/chi2_tables.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2_layer1.C — do not edit by hand\n\n";
        for (int ic = 0; ic < nChan; ++ic) {
            for (int is = 0; is < nScen; ++is) {
                TString label = Form("tab:chi2-%s-%s", channels[ic], scenarios[is]);
                out << "\\begin{table}[htbp]\n\\centering\n"
                    << "\\caption{Layer-1 per-$|\\eta|$-bin $\\chisq$ for "
                    << chanTeX[ic] << ", " << scenTeX[is] << " scenario.}\n"
                    << "\\label{" << label << "}\n"
                    << "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{r";
                for (int iv = 0; iv < nVar; ++iv) out << "|rrrrrr";
                out << "}\n\\toprule\n$|\\eta|$ bin ";
                for (int iv = 0; iv < nVar; ++iv)
                    out << "& \\multicolumn{6}{" << (iv < nVar - 1 ? "c|" : "c")
                        << "}{" << varTeX[iv] << "} ";
                out << "\\\\\n(range) ";
                for (int iv = 0; iv < nVar; ++iv) out << "& MC & M1 & M2 & M3 & M4 & Fud ";
                out << "\\\\\n\\midrule\n";
                for (int ie = 0; ie < kNEtaBins; ++ie) {
                    if (!R[ic][is].eta[ie][0].valid) continue;
                    out << Form("{[}%.2f, %.2f)", kEtaLimits[ie], kEtaLimits[ie + 1]);
                    for (int iv = 0; iv < nVar; ++iv)
                        for (int im = 0; im < nMeth; ++im)
                            out << " & " << fmtChi2(R[ic][is].eta[ie][iv].v[im]);
                    out << " \\\\\n";
                }
                out << "\\midrule\nAverage";
                for (int iv = 0; iv < nVar; ++iv)
                    for (int im = 0; im < nMeth; ++im)
                        out << " & " << fmtChi2(avg[ic][is][iv][im]);
                out << " \\\\\n\\bottomrule\n\\end{tabular}}\n\\end{table}\n\n";
            }
        }
        out.close();
        std::cout << "Written: " << fn << std::endl;
    }

    // Inclusive
    double incl[3][3][5][6] = {};
    bool   inclValid[3][3]  = {};
    for (int ic = 0; ic < nChan; ++ic) {
        for (int is = 0; is < nScen; ++is) {
            TString path = Form("%s/%s/%s/histograms.root",
                                baseDir, channels[ic], scenarios[is]);
            TFile* f = TFile::Open(path, "READ");
            if (!f || f->IsZombie()) { if (f) delete f; continue; }

            for (int iv = 0; iv < nVar; ++iv) {
                TH1D* hd  = (TH1D*)f->Get(Form("h_%s_data_computed", vars[iv]));
                TH1D* hmc = (TH1D*)f->Get(Form("h_%s_mc_computed",   vars[iv]));
                TH1D* hm1 = (TH1D*)f->Get(Form("h_%s_mc_M1",         vars[iv]));
                TH1D* hm2 = (TH1D*)f->Get(Form("h_%s_mc_M2",         vars[iv]));
                TH1D* hm3 = (TH1D*)f->Get(Form("h_%s_mc_M3",         vars[iv]));
                TH1D* hm4 = (TH1D*)f->Get(Form("h_%s_mc_M4",         vars[iv]));
                if (!hd || hd->GetEntries() < 100) continue;
                TH1D* nd  = normClone(hd,  "ni_d");
                TH1D* nmc = normClone(hmc, "ni_mc");
                TH1D* nm1 = normClone(hm1, "ni_m1");
                TH1D* nm2 = normClone(hm2, "ni_m2");
                TH1D* nm3 = normClone(hm3, "ni_m3");
                TH1D* nm4 = normClone(hm4, "ni_m4");
                incl[ic][is][iv][0] = (nd && nmc) ? chi2ndf(nmc, nd) : -1;
                incl[ic][is][iv][1] = (nd && nm1) ? chi2ndf(nm1, nd) : -1;
                incl[ic][is][iv][2] = (nd && nm2) ? chi2ndf(nm2, nd) : -1;
                incl[ic][is][iv][3] = (nd && nm3) ? chi2ndf(nm3, nd) : -1;
                incl[ic][is][iv][4] = (nd && nm4) ? chi2ndf(nm4, nd) : -1;

                TH1D* hds = (TH1D*)f->Get(Form("h_%s_data_stored", vars[iv]));
                TH1D* hmf = (TH1D*)f->Get(Form("h_%s_mc_fudged",   vars[iv]));
                TH1D* nds = normClone(hds, "ni_ds");
                TH1D* nmf = normClone(hmf, "ni_mf");
                incl[ic][is][iv][5] = (nds && nmf) ? chi2ndf(nmf, nds) : -1;
                inclValid[ic][is] = true;
                delete nd; delete nmc; delete nm1; delete nm2;
                delete nm3; delete nm4;
                delete nds; delete nmf;
            }
            f->Close(); delete f;
        }
    }

    {
        TString fn = Form("%s/chi2_inclusive.tex", reportDir);
        std::ofstream out(fn.Data());
        out << "% Auto-generated by extract_chi2_layer1.C — do not edit by hand\n"
            << "\\begin{table}[htbp]\n\\centering\n"
            << "\\caption{Layer-1 inclusive $\\chisq$ on integrated histograms.}\n"
            << "\\label{tab:chi2-inclusive-l1}\n"
            << "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{ll";
        for (int iv = 0; iv < nVar; ++iv) out << "|rrrrrr";
        out << "}\n\\toprule\n & ";
        for (int iv = 0; iv < nVar; ++iv)
            out << "& \\multicolumn{6}{" << (iv < nVar - 1 ? "c|" : "c")
                << "}{" << varTeX[iv] << "} ";
        out << "\\\\\nChannel & Scenario ";
        for (int iv = 0; iv < nVar; ++iv) out << "& MC & M1 & M2 & M3 & M4 & Fud ";
        out << "\\\\\n\\midrule\n";
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
        out << "\\bottomrule\n\\end{tabular}}\n\\end{table}\n";
        out.close();
        std::cout << "Written: " << fn << std::endl;
    }

    std::cout << "\n=== extract_chi2_layer1 complete ===" << std::endl;
}

#ifndef __CLING__
#include <TROOT.h>
#include <TApplication.h>
int main(int argc, char** argv) {
    TApplication app("app", nullptr, nullptr);
    gROOT->SetBatch(true);
    const char* bd  = (argc > 1) ? argv[1] : "../../output/Layer_1/eta_loose";
    const char* rd  = (argc > 2) ? argv[2] : "../../output/Layer_1/eta_loose/report";
    extract_chi2_layer1(bd, rd);
    return 0;
}
#endif
