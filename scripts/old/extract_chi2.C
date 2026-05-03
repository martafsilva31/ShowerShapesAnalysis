#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include <cstdio>

double computeChi2(TH1D* h1, TH1D* h2) {
    if (!h1 || !h2) return -1;
    double c = 0;
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double v1 = h1->GetBinContent(i), v2 = h2->GetBinContent(i);
        if (v2 > 0) c += (v1 - v2) * (v1 - v2) / v2;
    }
    return c;
}

void extract_chi2() {
    const char* files[] = {
        "../../output/data_mc_reweighting/Zeeg/histos.root",
        "../../output/data_mc_reweighting/Zmumug/histos.root"
    };
    const char* labels[] = {"Zeeg", "Zmumug"};

    for (int f = 0; f < 2; ++f) {
        TFile* tf = TFile::Open(files[f]);
        if (!tf || tf->IsZombie()) { printf("Cannot open %s\n", files[f]); continue; }

        printf("\n=== %s: Set A (branch) ===\n", labels[f]);
        printf("EtaBin |  R_eta: unf    fud  |  R_phi: unf    fud  | w_eta2: unf    fud\n");
        printf("-------|---------------------|---------------------|--------------------\n");
        for (int n = 0; n < 14; ++n) {
            TString s = Form("_%d", n);
            TH1D* hrd = (TH1D*)tf->Get("reta_data_br" + s);
            if (!hrd || hrd->Integral() < 1e-12) continue;
            printf("%6d | %6.2f %6.2f | %6.2f %6.2f | %6.2f %6.2f\n", n,
                computeChi2((TH1D*)tf->Get("reta_mc_unf" + s), hrd),
                computeChi2((TH1D*)tf->Get("reta_mc_fud" + s), hrd),
                computeChi2((TH1D*)tf->Get("rphi_mc_unf" + s), (TH1D*)tf->Get("rphi_data_br" + s)),
                computeChi2((TH1D*)tf->Get("rphi_mc_fud" + s), (TH1D*)tf->Get("rphi_data_br" + s)),
                computeChi2((TH1D*)tf->Get("weta2_mc_unf" + s), (TH1D*)tf->Get("weta2_data_br" + s)),
                computeChi2((TH1D*)tf->Get("weta2_mc_fud" + s), (TH1D*)tf->Get("weta2_data_br" + s)));
        }

        printf("\n=== %s: Set B (cell-computed) ===\n", labels[f]);
        printf("EtaBin |  R_eta:  MC    M1    M2  |  R_phi:  MC    M1    M2  | w_eta2:  MC    M1    M2\n");
        printf("-------|-------------------------|-------------------------|------------------------\n");
        for (int n = 0; n < 14; ++n) {
            TString s = Form("_%d", n);
            TH1D* hrd = (TH1D*)tf->Get("reta_data_cell" + s);
            if (!hrd || hrd->Integral() < 1e-12) continue;
            printf("%6d | %6.2f %5.2f %5.2f | %6.2f %5.2f %5.2f | %6.2f %5.2f %5.2f\n", n,
                computeChi2((TH1D*)tf->Get("reta_mc_cell" + s), hrd),
                computeChi2((TH1D*)tf->Get("reta_m1" + s), hrd),
                computeChi2((TH1D*)tf->Get("reta_m2" + s), hrd),
                computeChi2((TH1D*)tf->Get("rphi_mc_cell" + s), (TH1D*)tf->Get("rphi_data_cell" + s)),
                computeChi2((TH1D*)tf->Get("rphi_m1" + s), (TH1D*)tf->Get("rphi_data_cell" + s)),
                computeChi2((TH1D*)tf->Get("rphi_m2" + s), (TH1D*)tf->Get("rphi_data_cell" + s)),
                computeChi2((TH1D*)tf->Get("weta2_mc_cell" + s), (TH1D*)tf->Get("weta2_data_cell" + s)),
                computeChi2((TH1D*)tf->Get("weta2_m1" + s), (TH1D*)tf->Get("weta2_data_cell" + s)),
                computeChi2((TH1D*)tf->Get("weta2_m2" + s), (TH1D*)tf->Get("weta2_data_cell" + s)));
        }
        tf->Close();
    }
}
