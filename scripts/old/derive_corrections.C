///////////////////////////////////////////////////////////////////////////////
// derive_corrections.C
//
// Derives cell-energy reweighting corrections following Francisco's code
// (ATL-COM-PHYS-2021-640) adapted for Run-3 ntuples.
//
// -----------------------------------------------------------------------
// METHOD 1 — Flat shift per cell (CalCoef.C approach, document §5.1):
//   Delta_k = <e_data>_k - <e_MC>_k    (population mean fraction difference)
//   Applied as: E_new_k = E_MC_k + Delta_k * E_total_MC
//   sum(Delta_k) = 0 by construction => total energy preserved, NO renorm.
//   Data and MC are processed independently (no matching needed).
//
// METHOD 2 — DeltaR-matched TProfile (Francisco's tree.C approach):
//   Match data-MC photon pairs by DeltaR < 0.1 in (eta, phi).
//   For each cell k, eta bin n:
//     alpha(f_MC) = <f_data - f_MC> as a function of f_MC
//   Stored as TProfile in corrections.root.
//   Application: E_new_k = E_old_k + alpha(f_MC_k) * E_total
//
// METHOD 3 — Shift+Stretch per cell (document §5.2):
//   f'_k = (sigma_data_k / sigma_MC_k) * (f_MC_k - mu_MC_k) + mu_data_k
//   Matches both mean and variance of each cell's fraction distribution.
//   Renormalization applied after to preserve total energy.
//
// -----------------------------------------------------------------------
// Selection (note §3.5, relaxed for cell reweighting):
//   FSR mass windows, pT > 7 GeV, |eta| < 2.37 (no crack),
//   DeltaR(gamma,l) > 0.4, DeltaR(l,l) > 0.2, 77 cells,
//   central cell hottest (raw), truth match (MC only).
//   No photon ID, no isolation.
//
// Reference: ATL-COM-PHYS-2021-640, ATL-COM-PHYS-2025-662
//
// Usage:
//   root -l -b -q 'derive_corrections.C("data.root","mc.root","out.root")'
//   root -l -b -q 'derive_corrections.C("data.root","mc.root","out.root",500000)'
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TProfile.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace datamc;


// ======================================================================
// M1: Accumulator for per-cell fraction statistics (population means)
// Following CalCoef.C — processes data and MC independently.
// ======================================================================
struct CellStats {
    double sumF[kClusterSize][kNEtaBins];   // sum of (weight * cell fraction)
    double sumF2[kClusterSize][kNEtaBins];  // sum of (weight * fraction^2)
    double sumW[kNEtaBins];                 // sum of weights (=count for data)

    void reset() {
        for (int k = 0; k < kClusterSize; ++k)
            for (int n = 0; n < kNEtaBins; ++n) {
                sumF[k][n] = 0;
                sumF2[k][n] = 0;
            }
        for (int n = 0; n < kNEtaBins; ++n)
            sumW[n] = 0;
    }

    double mean(int k, int n) const {
        return sumW[n] > 0 ? sumF[k][n] / sumW[n] : 0;
    }

    double rms(int k, int n) const {
        if (sumW[n] <= 0) return 0;
        double m2 = sumF2[k][n] / sumW[n];
        double m1 = sumF[k][n] / sumW[n];
        double var = m2 - m1 * m1;
        return var > 0 ? std::sqrt(var) : 0;
    }
};


// ======================================================================
// M1: Process one sample — accumulate population fraction statistics.
//
// This replicates what NTUP.C does when filling energyProfile histograms
// (lines 700-780): for each event, compute E_k/E_tot per cell and
// accumulate the sum.  CalCoef.C then takes the difference of the means.
// ======================================================================
void processM1(const char* filename, bool isMC,
               CellStats& stats, Long64_t maxEvents) {
    stats.reset();

    TChain* t = makeChain(filename);
    if (!t || t->GetEntries() == 0) {
        std::cerr << "ERROR: Cannot open or empty chain for " << filename << std::endl;
        delete t; return;
    }

    // Branch variables
    std::vector<double>* cellE = nullptr;
    Int_t cellSize = 0;
    Float_t eta2 = 0, photonPt = 0;
    Double_t mll = 0, mllg = 0;
    Double_t dR_l1_ph = 0, dR_l2_ph = 0, dR_ll = 0;
    Bool_t truthMatch = false;
    Float_t mcWgt = 1, puWgt = 1;
    Int_t tightID = 0;
    Bool_t isConv = false;

    t->SetBranchStatus("*", 0);
    auto enable = [&](const char* name, void* addr) {
        if (t->GetBranch(name)) {
            t->SetBranchStatus(name, 1);
            t->SetBranchAddress(name, addr);
        }
    };

    enable(kCellBranch, &cellE);
    enable(kCellSizeBranch, &cellSize);
    enable(kEtaBranch, &eta2);
    enable(kPhotonPtBranch, &photonPt);
    enable(kMllBranch, &mll);
    enable(kMllgBranch, &mllg);
    enable(kDRLepton1Branch, &dR_l1_ph);
    enable(kDRLepton2Branch, &dR_l2_ph);
    enable(kDRllBranch, &dR_ll);
    enable(kTightIDBranch, &tightID);
    enable(kIsConvBranch, &isConv);
    if (isMC) {
        enable(kTruthMatchBranch, &truthMatch);
        enable(kMCWeightBranch, &mcWgt);
        enable(kPUWeightBranch, &puWgt);
    }

    Long64_t nEntries = t->GetEntries();
    if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
    std::cout << "M1 " << (isMC ? "MC" : "Data") << ": " << filename
              << " (" << nEntries << " entries)" << std::endl;

    Long64_t nUsed = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (i > 0 && i % 2000000 == 0)
            std::cout << "  " << i << " / " << nEntries << std::endl;

        if (!passSelection(photonPt, eta2, mll, mllg,
                           dR_l1_ph, dR_l2_ph, dR_ll,
                           cellSize, truthMatch, isMC,
                           tightID, isConv))
            continue;

        if (!cellE || (int)cellE->size() != kClusterSize) continue;

        // Quality: central cell (k=38) must be hottest in RAW cluster
        // (matches Francisco's maxK == 38 check, NTUP.cxx line 149)
        if (!isHealthyCluster(*cellE)) continue;

        int etaBin = findEtaBin(std::fabs(eta2));
        if (etaBin < 0) continue;

        // Total energy (no clamping, same as Francisco)
        double Etot = 0;
        for (int k = 0; k < kClusterSize; ++k) Etot += cellE->at(k);
        if (Etot <= 0) continue;

        // UNWEIGHTED for M1 derivation — Francisco's NTUP.C L740 accumulates
        // raw fractions without eventweight; normalization via Scale(1/Integral).
        // MC generator weights (~1e7) would corrupt the population mean.
        double w = 1.0;

        // Accumulate weighted fraction and fraction^2 for each cell
        for (int k = 0; k < kClusterSize; ++k) {
            double frac = cellE->at(k) / Etot;
            stats.sumF[k][etaBin] += w * frac;
            stats.sumF2[k][etaBin] += w * frac * frac;
        }
        stats.sumW[etaBin] += w;
        nUsed++;
    }

    std::cout << "  Selected: " << nUsed << " events" << std::endl;
    delete t;
}


// ======================================================================
// M2: DeltaR-matched TProfile correction (Francisco's tree.C approach)
//
// 1. Read all data and MC photons passing selection
// 2. Match data-MC pairs by DeltaR < 0.1 (best match per MC photon)
// 3. Dedup: keep only the closest MC partner per data photon
// 4. Fill TProfile: x = mc_frac, y = (data_frac - mc_frac)
//    This gives alpha(f_MC) — the mean shift as a function of MC fraction.
// 5. Write TProfiles to output file
//
// Application: E_new_k = E_old_k + alpha(f_MC_k) * E_total
//
// Reference: Francisco's tree.C + NTUP.cxx (ATL-COM-PHYS-2021-640)
// ======================================================================
void deriveM2Matched(const char* dataFile, const char* mcFile,
                     Long64_t maxEvents, TFile* fout) {

    std::cout << "\n=== Method 2: DeltaR-Matched TProfile (Francisco) ===" << std::endl;

    // Photon info for matching
    struct PhotonInfo {
        double eta, phi;
        double cellE[kClusterSize];
        double Etot;
        int etaBin;
    };

    // Lambda to read photons from a file
    auto readPhotons = [&](const char* filename, bool isMC) -> std::vector<PhotonInfo> {
        std::vector<PhotonInfo> photons;
        TChain* t = makeChain(filename);
        if (!t || t->GetEntries() == 0) { delete t; return photons; }

        std::vector<double>* cellE = nullptr;
        Int_t cellSize = 0;
        Float_t eta2 = 0, photonPt = 0, photonEta = 0, photonPhi = 0;
        Double_t mll = 0, mllg = 0;
        Double_t dR_l1_ph = 0, dR_l2_ph = 0, dR_ll = 0;
        Bool_t truthMatch = false;
        Int_t tightID = 0;
        Bool_t isConv = false;
        Float_t mcWgt = 1, puWgt = 1;

        t->SetBranchStatus("*", 0);
        auto enable = [&](const char* name, void* addr) {
            if (t->GetBranch(name)) {
                t->SetBranchStatus(name, 1);
                t->SetBranchAddress(name, addr);
            }
        };
        enable(kCellBranch, &cellE);
        enable(kCellSizeBranch, &cellSize);
        enable(kEtaBranch, &eta2);
        enable(kPhotonPtBranch, &photonPt);
        enable(kPhotonEtaBranch, &photonEta);
        enable(kPhotonPhiBranch, &photonPhi);
        enable(kMllBranch, &mll);
        enable(kMllgBranch, &mllg);
        enable(kDRLepton1Branch, &dR_l1_ph);
        enable(kDRLepton2Branch, &dR_l2_ph);
        enable(kDRllBranch, &dR_ll);
        enable(kTightIDBranch, &tightID);
        enable(kIsConvBranch, &isConv);
        if (isMC) {
            enable(kTruthMatchBranch, &truthMatch);
            enable(kMCWeightBranch, &mcWgt);
            enable(kPUWeightBranch, &puWgt);
        }

        Long64_t nEntries = t->GetEntries();
        if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;

        for (Long64_t i = 0; i < nEntries; ++i) {
            t->GetEntry(i);
            if (!passSelection(photonPt, eta2, mll, mllg,
                               dR_l1_ph, dR_l2_ph, dR_ll,
                               cellSize, truthMatch, isMC,
                               tightID, isConv))
                continue;
            if (!cellE || (int)cellE->size() != kClusterSize) continue;
            if (!isHealthyCluster(*cellE)) continue;

            int etaBin = findEtaBin(std::fabs(eta2));
            if (etaBin < 0) continue;

            double Etot = 0;
            for (int k = 0; k < kClusterSize; ++k) Etot += cellE->at(k);
            if (Etot <= 0) continue;

            PhotonInfo p;
            p.eta = photonEta;
            p.phi = photonPhi;
            for (int k = 0; k < kClusterSize; ++k) p.cellE[k] = cellE->at(k);
            p.Etot = Etot;
            p.etaBin = etaBin;
            photons.push_back(p);
        }
        std::cout << "  " << (isMC ? "MC" : "Data") << ": " << photons.size()
                  << " photons from " << filename << std::endl;
        delete t;
        return photons;
    };

    auto dataPhotons = readPhotons(dataFile, false);
    auto mcPhotons   = readPhotons(mcFile,   true);

    // Match: for each MC photon, find nearest data photon by DeltaR < 0.1
    // (Francisco's NTUP.cxx matching algorithm)
    // Optimized with eta-grid spatial index: O(N_MC * K) instead of O(N_MC * N_data)
    const double kMatchDR = 0.1;
    struct Match {
        int mcIdx, dataIdx;
        double dR;
    };
    std::vector<Match> matches;

    // Build spatial index: bin data photons by eta into grid cells
    // Grid cell size = kMatchDR, so only need to search ±1 neighboring cell
    const double kEtaGridLo = -3.0, kEtaGridHi = 3.0;
    const double kGridCell = kMatchDR;
    const int kNGridBins = (int)((kEtaGridHi - kEtaGridLo) / kGridCell) + 1;
    std::vector<std::vector<int>> etaGrid(kNGridBins);

    for (int id = 0; id < (int)dataPhotons.size(); ++id) {
        int bin = (int)((dataPhotons[id].eta - kEtaGridLo) / kGridCell);
        if (bin >= 0 && bin < kNGridBins)
            etaGrid[bin].push_back(id);
    }

    for (int imc = 0; imc < (int)mcPhotons.size(); ++imc) {
        if (imc > 0 && imc % 500000 == 0)
            std::cout << "  Matching: " << imc << " / " << mcPhotons.size() << std::endl;

        int etaBin0 = (int)((mcPhotons[imc].eta - kEtaGridLo) / kGridCell);
        double bestDR = 999;
        int bestData = -1;

        // Search only neighboring eta grid cells (±1)
        for (int db = -1; db <= 1; ++db) {
            int gbin = etaBin0 + db;
            if (gbin < 0 || gbin >= kNGridBins) continue;
            for (int id : etaGrid[gbin]) {
                double deta = mcPhotons[imc].eta - dataPhotons[id].eta;
                double dphi = mcPhotons[imc].phi - dataPhotons[id].phi;
                while (dphi >  M_PI) dphi -= 2 * M_PI;
                while (dphi < -M_PI) dphi += 2 * M_PI;
                double dR = std::sqrt(deta * deta + dphi * dphi);
                if (dR < bestDR) { bestDR = dR; bestData = id; }
            }
        }
        if (bestDR < kMatchDR && bestData >= 0)
            matches.push_back({imc, bestData, bestDR});
    }

    std::cout << "  Raw matches (DeltaR < " << kMatchDR << "): "
              << matches.size() << std::endl;

    // Dedup: keep only the closest MC partner per data photon
    // Sort by dataIdx, then dR ascending; keep first per dataIdx.
    std::sort(matches.begin(), matches.end(),
              [](const Match& a, const Match& b) {
                  if (a.dataIdx != b.dataIdx) return a.dataIdx < b.dataIdx;
                  return a.dR < b.dR;
              });

    std::vector<Match> uniqueMatches;
    int lastData = -1;
    for (auto& m : matches) {
        if (m.dataIdx != lastData) {
            uniqueMatches.push_back(m);
            lastData = m.dataIdx;
        }
    }
    matches.swap(uniqueMatches);

    std::cout << "  After dedup (unique data photons): "
              << matches.size() << std::endl;

    // Create TProfile per cell per eta bin
    TProfile* prof[kClusterSize][kNEtaBins];
    for (int k = 0; k < kClusterSize; ++k) {
        for (int n = 0; n < kNEtaBins; ++n) {
            TString name = Form("m2_alpha_Cell_%d_Eta_%1.2f_%1.2f",
                                k + 1, kEtaLimits[n], kEtaLimits[n + 1]);
            prof[k][n] = new TProfile(name, name,
                                      kNFracBins, kFracBins);
        }
    }

    // Fill TProfiles: alpha = data_frac - mc_frac as function of mc_frac
    int nFilled[kNEtaBins] = {};
    for (auto& m : matches) {
        const PhotonInfo& mc   = mcPhotons[m.mcIdx];
        const PhotonInfo& data = dataPhotons[m.dataIdx];
        int etaBin = mc.etaBin;
        for (int k = 0; k < kClusterSize; ++k) {
            double mc_frac   = mc.cellE[k] / mc.Etot;
            double data_frac = data.cellE[k] / data.Etot;
            double alpha     = data_frac - mc_frac;
            prof[k][etaBin]->Fill(mc_frac, alpha);
        }
        nFilled[etaBin]++;
    }

    // Write and print summary
    fout->cd();
    std::cout << "\nEtaBin | Matched pairs | <|alpha|> cell 38" << std::endl;
    std::cout << "-------|---------------|------------------" << std::endl;
    for (int n = 0; n < kNEtaBins; ++n) {
        for (int k = 0; k < kClusterSize; ++k)
            prof[k][n]->Write();
        if (nFilled[n] > 0) {
            double sumAbsAlpha = 0;
            int nBins = prof[kCentralCell][n]->GetNbinsX();
            int nUsedBins = 0;
            for (int b = 1; b <= nBins; ++b) {
                if (prof[kCentralCell][n]->GetBinEntries(b) > 0) {
                    sumAbsAlpha += std::fabs(prof[kCentralCell][n]->GetBinContent(b));
                    nUsedBins++;
                }
            }
            std::cout << Form("%2d     |  %6d       | %.6f",
                              n, nFilled[n],
                              nUsedBins > 0 ? sumAbsAlpha / nUsedBins : 0.)
                      << std::endl;
        }
    }

    // Cleanup
    for (int k = 0; k < kClusterSize; ++k)
        for (int n = 0; n < kNEtaBins; ++n)
            delete prof[k][n];
}


// ======================================================================
// Main entry point
// ======================================================================
int derive_corrections(const char* dataFile,
                       const char* mcFile,
                       const char* outputFile,
                       Long64_t maxEvents = -1) {

    std::cout << "========================================" << std::endl;
    std::cout << "  Deriving Cell Corrections (M1-M3)" << std::endl;
    std::cout << "========================================" << std::endl;

    // ==================================================================
    // M1: Population statistics
    // ==================================================================
    std::cout << "\n=== Method 1: Population Statistics ===" << std::endl;

    CellStats dataStats, mcStats;
    processM1(dataFile, false, dataStats, maxEvents);
    processM1(mcFile,   true,  mcStats,   maxEvents);

    // ==================================================================
    // Write M1 correction histograms
    // ==================================================================
    TFile* fout = TFile::Open(outputFile, "RECREATE");

    std::cout << "\n--- M1 Correction Summary ---" << std::endl;
    std::cout << "EtaBin |    W_data       W_MC    |  Checksum(Delta)"
              << std::endl;
    std::cout << "-------|--------------------------|-------------------"
              << std::endl;

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        // cellsDelta: TH2D with flat shift per cell (phi x eta layout)
        // Same object as CalCoef.C produces.
        // Naming: "cellsDelta_Eta_LO_HI" (matches Francisco)
        TH2D* h_delta = new TH2D(
            "cellsDelta_" + suffix,
            "#Delta_{k} = <e_{data}> - <e_{MC}>;phi;eta",
            kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);

        // Eta and phi projections (same as CalCoef.C)
        TH1D* h_etaDelta = new TH1D(
            "etaProfileDelta_" + suffix,
            "Eta projection of #Delta;eta;#Delta",
            kEtaSize, 0, kEtaSize);
        TH1D* h_phiDelta = new TH1D(
            "phiProfileDelta_" + suffix,
            "Phi projection of #Delta;phi;#Delta",
            kPhiSize, 0, kPhiSize);

        double checksum = 0;

        // Fill: Delta_k = <e_data>_k - <e_MC>_k
        // Loop (phi, eta) same order as Francisco's CalCoef.C
        for (int phi = 1; phi <= kPhiSize; ++phi) {
            for (int eta = 1; eta <= kEtaSize; ++eta) {
                int k = (phi - 1) + kPhiSize * (eta - 1);
                double delta = dataStats.mean(k, n) - mcStats.mean(k, n);
                h_delta->SetBinContent(phi, eta, delta);
                checksum += delta;
            }
        }

        // Eta projection: sum over phi for each eta column
        for (int eta = 0; eta < kEtaSize; ++eta) {
            double sum = 0;
            for (int phi = 0; phi < kPhiSize; ++phi) {
                int k = phi + kPhiSize * eta;
                sum += dataStats.mean(k, n) - mcStats.mean(k, n);
            }
            h_etaDelta->SetBinContent(eta + 1, sum);
        }

        // Phi projection: sum over eta for each phi row
        for (int phi = 0; phi < kPhiSize; ++phi) {
            double sum = 0;
            for (int eta = 0; eta < kEtaSize; ++eta) {
                int k = phi + kPhiSize * eta;
                sum += dataStats.mean(k, n) - mcStats.mean(k, n);
            }
            h_phiDelta->SetBinContent(phi + 1, sum);
        }

        h_delta->Write();
        h_etaDelta->Write();
        h_phiDelta->Write();

        std::cout << Form("%2d     | %10.1f  %10.1f  |  %.6e",
                          n, dataStats.sumW[n], mcStats.sumW[n], checksum)
                  << std::endl;
    }

    // ==================================================================
    // M3: Store shift+stretch parameters (document §5.2)
    //
    // f'_k = (sigma_data / sigma_MC) * (f_MC - mu_MC) + mu_data
    //
    // We write three TH2D per eta bin:
    //   cellsMeanData  — mu_data per cell
    //   cellsMeanMC    — mu_MC per cell
    //   cellsSigmaRatio — sigma_data / sigma_MC per cell
    // ==================================================================
    std::cout << "\n--- M3 Shift+Stretch Parameters ---" << std::endl;
    std::cout << "EtaBin |  <sigma_ratio>  min    max" << std::endl;
    std::cout << "-------|----------------------------" << std::endl;

    for (int n = 0; n < kNEtaBins; ++n) {
        TString suffix = Form("Eta_%1.2f_%1.2f",
                              kEtaLimits[n], kEtaLimits[n + 1]);

        TH2D* h_meanData = new TH2D(
            "cellsMeanData_" + suffix,
            "<f_{data}>;phi;eta",
            kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);
        TH2D* h_meanMC = new TH2D(
            "cellsMeanMC_" + suffix,
            "<f_{MC}>;phi;eta",
            kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);
        TH2D* h_sigRatio = new TH2D(
            "cellsSigmaRatio_" + suffix,
            "#sigma_{data}/#sigma_{MC};phi;eta",
            kPhiSize, 0, kPhiSize, kEtaSize, 0, kEtaSize);

        double sumRatio = 0;
        double minRatio = 1e9, maxRatio = -1e9;
        int nCells = 0;

        for (int phi = 1; phi <= kPhiSize; ++phi) {
            for (int eta = 1; eta <= kEtaSize; ++eta) {
                int k = (phi - 1) + kPhiSize * (eta - 1);
                double muD  = dataStats.mean(k, n);
                double muMC = mcStats.mean(k, n);
                double sD   = dataStats.rms(k, n);
                double sMC  = mcStats.rms(k, n);

                double ratio = (sMC > 1e-12) ? sD / sMC : 1.0;

                // Cap sigma ratio to prevent extreme corrections
                // in low-statistics bins
                if (ratio < kSigmaRatioMin || ratio > kSigmaRatioMax) {
                    if (n != 8) // don't warn for empty crack bin
                        std::cout << "  WARNING: Capping sigma_ratio "
                                  << ratio << " -> ["
                                  << kSigmaRatioMin << ", "
                                  << kSigmaRatioMax << "] for cell "
                                  << k << " eta bin " << n << std::endl;
                    ratio = std::max(kSigmaRatioMin,
                                     std::min(kSigmaRatioMax, ratio));
                }

                h_meanData->SetBinContent(phi, eta, muD);
                h_meanMC->SetBinContent(phi, eta, muMC);
                h_sigRatio->SetBinContent(phi, eta, ratio);

                sumRatio += ratio;
                if (ratio < minRatio) minRatio = ratio;
                if (ratio > maxRatio) maxRatio = ratio;
                nCells++;
            }
        }

        h_meanData->Write();
        h_meanMC->Write();
        h_sigRatio->Write();

        if (nCells > 0) {
            std::cout << Form("%2d     |  %8.4f  %6.4f %6.4f",
                              n, sumRatio / nCells, minRatio, maxRatio)
                      << std::endl;
        }
    }

    // ==================================================================
    // M2: DeltaR-matched TProfile corrections
    // ==================================================================
    deriveM2Matched(dataFile, mcFile, maxEvents, fout);

    fout->Close();
    delete fout;

    std::cout << "\nCorrections written to: " << outputFile << std::endl;
    return 0;
}
