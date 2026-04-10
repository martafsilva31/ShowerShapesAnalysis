// check_fudge.C — Check if fudge factors are actually present in MC
#include "TFile.h"
#include "TTree.h"
#include <iostream>

void check_fudge() {
    const char* files[] = {
        "/dcache/atlas/mfernand/qt_ntuples/data24/mc_eegamma.root",
        "/dcache/atlas/mfernand/qt_ntuples/data24/mc_mumugamma.root"
    };
    const char* labels[] = {"mc_eegamma", "mc_mumugamma"};

    for (int ifile = 0; ifile < 2; ++ifile) {
        TFile* f = TFile::Open(files[ifile], "READ");
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << files[ifile] << std::endl; continue; }
        TTree* t = (TTree*)f->Get("tree");
        if (!t) { std::cerr << "No tree 'tree' in " << files[ifile] << std::endl; continue; }

        Float_t reta_fud = 0, reta_unf = 0;
        Float_t rphi_fud = 0, rphi_unf = 0;
        Float_t weta2_fud = 0, weta2_unf = 0;

        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("photon.reta", 1);
        t->SetBranchStatus("photon.rphi", 1);
        t->SetBranchStatus("photon.weta2", 1);
        t->SetBranchStatus("photon.unfudged_reta", 1);
        t->SetBranchStatus("photon.unfudged_rphi", 1);
        t->SetBranchStatus("photon.unfudged_weta2", 1);

        t->SetBranchAddress("photon.reta", &reta_fud);
        t->SetBranchAddress("photon.rphi", &rphi_fud);
        t->SetBranchAddress("photon.weta2", &weta2_fud);
        t->SetBranchAddress("photon.unfudged_reta", &reta_unf);
        t->SetBranchAddress("photon.unfudged_rphi", &rphi_unf);
        t->SetBranchAddress("photon.unfudged_weta2", &weta2_unf);

        Long64_t nEntries = t->GetEntries();
        Long64_t nCheck = std::min(nEntries, (Long64_t)100000);

        int nDiff_reta = 0, nDiff_rphi = 0, nDiff_weta2 = 0;
        double sumAbsDiff_reta = 0, sumAbsDiff_rphi = 0, sumAbsDiff_weta2 = 0;

        std::cout << "\n========================================" << std::endl;
        std::cout << labels[ifile] << "  (" << nEntries << " total, checking " << nCheck << ")" << std::endl;
        std::cout << "========================================" << std::endl;

        // Print first 10 entries
        std::cout << "\nFirst 10 entries (fudged vs unfudged):" << std::endl;
        std::cout << "  Entry   reta_fud   reta_unf   rphi_fud   rphi_unf   weta2_fud  weta2_unf" << std::endl;

        for (Long64_t i = 0; i < nCheck; ++i) {
            t->GetEntry(i);

            if (i < 10) {
                printf("  %5lld   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f   %8.6f\n",
                       i, reta_fud, reta_unf, rphi_fud, rphi_unf, weta2_fud, weta2_unf);
            }

            if (reta_fud != reta_unf) { nDiff_reta++; sumAbsDiff_reta += std::fabs(reta_fud - reta_unf); }
            if (rphi_fud != rphi_unf) { nDiff_rphi++; sumAbsDiff_rphi += std::fabs(rphi_fud - rphi_unf); }
            if (weta2_fud != weta2_unf) { nDiff_weta2++; sumAbsDiff_weta2 += std::fabs(weta2_fud - weta2_unf); }
        }

        std::cout << "\nSummary (out of " << nCheck << " checked):" << std::endl;
        std::cout << "  reta  differs in " << nDiff_reta << " events";
        if (nDiff_reta > 0) std::cout << "  (mean |diff| = " << sumAbsDiff_reta / nDiff_reta << ")";
        std::cout << std::endl;
        std::cout << "  rphi  differs in " << nDiff_rphi << " events";
        if (nDiff_rphi > 0) std::cout << "  (mean |diff| = " << sumAbsDiff_rphi / nDiff_rphi << ")";
        std::cout << std::endl;
        std::cout << "  weta2 differs in " << nDiff_weta2 << " events";
        if (nDiff_weta2 > 0) std::cout << "  (mean |diff| = " << sumAbsDiff_weta2 / nDiff_weta2 << ")";
        std::cout << std::endl;

        f->Close();
    }
}
