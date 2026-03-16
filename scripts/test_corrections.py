#!/usr/bin/env python3
"""Test different correction strategies for weta2."""
import ROOT

ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace qw {
    double P0A[]={0.0045,0.005375,-0.0562},P1A[]={-0.0016,-0.0215,0.114},P2A[]={-0.0866,0.0215,-0.053};
    double P0B[]={0.0039,0.005075,-0.0324},P1B[]={0.00816,-0.0203,0.0653},P2B[]={-0.145,0.0203,-0.0286};
    double P0C[]={0.0047,0.0035},P1C[]={-0.0184,-0.0139},P2C[]={0.0180,0.0137};

    double RelPos(double eta, double etacell) {
        return std::fmod(std::fabs(eta - etacell - 0.0125), 0.025) / 0.025;
    }
    double Correct(double eta, double etacell, double w) {
        double ae = std::fabs(eta), er = RelPos(eta, etacell), co = 0;
        if (ae<0.8){int j=(er<0.1)?0:(er<0.9?1:2);co=P0A[j]+P1A[j]*er+P2A[j]*er*er;}
        else if(ae<1.8){int j=(er<0.1)?0:(er<0.9?1:2);co=P0B[j]+P1B[j]*er+P2B[j]*er*er;}
        else if(ae<2.0){co=P0C[0]+P1C[0]*er+P2C[0]*er*er;}
        else if(ae<2.5){co=P0C[1]+P1C[1]*er+P2C[1]*er*er;}
        return w - co;
    }
}

struct Stats {
    double sd,sd2; long long n;
    Stats():sd(0),sd2(0),n(0){}
    void add(double d){sd+=d;sd2+=d*d;n++;}
    double mean()const{return n>0?sd/n:0;}
    double stdev()const{return n>1?std::sqrt(std::max(0.0,sd2/n-std::pow(sd/n,2))):0;}
};

void test_corrections(TTree* tree) {
    std::vector<double>*cl=nullptr,*ce=nullptr;
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("photon.7x11ClusterLr2E",1);
    tree->SetBranchStatus("photon.7x11ClusterLr2Eta",1);
    tree->SetBranchStatus("photon.unfudged_weta2",1);
    tree->SetBranchStatus("photon_cluster.eta2",1);
    tree->SetBranchAddress("photon.7x11ClusterLr2E",&cl);
    tree->SetBranchAddress("photon.7x11ClusterLr2Eta",&ce);
    auto*le=tree->GetLeaf("photon_cluster.eta2");
    auto*lu=tree->GetLeaf("photon.unfudged_weta2");
    
    Stats s_raw, s_corr_hottest, s_corr_centre, s_corr_emean, s_nocorr;
    
    Long64_t nev = tree->GetEntries();
    for (Long64_t i=0; i<nev; i++) {
        tree->GetEntry(i);
        if(!cl||cl->size()!=77||!ce||ce->size()!=77) continue;
        double eta2 = le->GetValue();
        double wu = lu->GetValue();
        
        // Find hottest
        double mx=-1e99; int si=-1;
        for(int j=0;j<77;j++) if((*cl)[j]>mx){mx=(*cl)[j];si=j;}
        if(mx<=0) continue;
        int sr=si/11,sc=si%11;
        
        // Compute raw weta2 in 3x5 around hottest
        int r0=std::max(0,sr-1),r1=std::min(6,sr+1);
        int c0=std::max(0,sc-2),c1=std::min(10,sc+2);
        double sE=0,sEe=0,sEe2=0;
        for(int r=r0;r<=r1;r++) for(int c=c0;c<=c1;c++){
            int k=r*11+c; double e=(*cl)[k],et=(*ce)[k];
            sE+=e; sEe+=e*et; sEe2+=e*et*et;
        }
        if(sE<=0) continue;
        double var = sEe2/sE - std::pow(sEe/sE,2);
        if(var<0) continue;
        double raw = std::sqrt(var);
        double eta_mean = sEe/sE;  // energy-weighted mean from our 3x5 window
        
        s_raw.add(raw - wu);
        
        // Strategy A: correction with etacell = hottest cell (current approach)
        double wA = qw::Correct(eta2, (*ce)[si], raw);
        s_corr_hottest.add(wA - wu);
        
        // Strategy B: correction with etacell = grid centre eta (= eta0)
        double wB = qw::Correct(eta2, (*ce)[38], raw);
        s_corr_centre.add(wB - wu);
        
        // Strategy C: correction with eta = computed mean, etacell = hottest
        double wC = qw::Correct(eta_mean, (*ce)[si], raw);
        s_corr_emean.add(wC - wu);
        
        // Strategy D: no correction at all
        s_nocorr.add(raw - wu);
        
        if((i+1)%1000000==0) printf("  %lld/%lld\n",i+1,nev);
    }
    
    printf("\n%-40s %12s %12s %12s\n", "Strategy", "Mean Δ", "Std Δ", "N");
    printf("%-40s %12s %12s %12s\n", "--------", "------", "-----", "-");
    printf("%-40s %+12.6f %12.6f %12lld\n", "Raw CellEta (no correction)", s_raw.mean(), s_raw.stdev(), s_raw.n);
    printf("%-40s %+12.6f %12.6f %12lld\n", "Corr: etacell=hottest (current)", s_corr_hottest.mean(), s_corr_hottest.stdev(), s_corr_hottest.n);
    printf("%-40s %+12.6f %12.6f %12lld\n", "Corr: etacell=grid_centre (eta0)", s_corr_centre.mean(), s_corr_centre.stdev(), s_corr_centre.n);
    printf("%-40s %+12.6f %12.6f %12lld\n", "Corr: eta=computed_mean, cell=hottest", s_corr_emean.mean(), s_corr_emean.stdev(), s_corr_emean.n);
}
""")

f = ROOT.TFile.Open('../ntuples/mc23e_v3/mc23e_v3_Zeeg.root')
tree = f.Get('tree')
ROOT.test_corrections(tree)
f.Close()
