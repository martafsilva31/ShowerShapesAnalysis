#!/usr/bin/env python3
"""Quick eta-binned weta2 discrepancy analysis."""
import ROOT

ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>

struct EBS {
    double sda,sda2,sdc,sdc2,sdr,sdr2;
    long long n, noff;
    EBS():sda(0),sda2(0),sdc(0),sdc2(0),sdr(0),sdr2(0),n(0),noff(0){}
};

void run_eta_analysis(TTree* tree, Long64_t maxev) {
    const int NB=8;
    double ed[NB+1]={0,0.4,0.8,1.2,1.37,1.52,1.8,2.0,2.5};
    EBS st[NB];
    
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
    
    double P0A[]={0.0045,0.005375,-0.0562},P1A[]={-0.0016,-0.0215,0.114},P2A[]={-0.0866,0.0215,-0.053};
    double P0B[]={0.0039,0.005075,-0.0324},P1B[]={0.00816,-0.0203,0.0653},P2B[]={-0.145,0.0203,-0.0286};
    double P0C[]={0.0047,0.0035},P1C[]={-0.0184,-0.0139},P2C[]={0.0180,0.0137};
    
    Long64_t nev=std::min(maxev,tree->GetEntries());
    for(Long64_t i=0;i<nev;i++){
        tree->GetEntry(i);
        if(!cl||cl->size()!=77||!ce||ce->size()!=77) continue;
        double e2=le->GetValue(),wu=lu->GetValue(),ae=std::fabs(e2);
        int b=-1;
        for(int j=0;j<NB;j++) if(ae>=ed[j]&&ae<ed[j+1]){b=j;break;}
        if(b<0) continue;
        
        double mx=-1e99; int si=-1;
        for(int j=0;j<77;j++) if((*cl)[j]>mx){mx=(*cl)[j];si=j;}
        if(mx<=0) continue;
        st[b].n++;
        if(si!=38) st[b].noff++;
        
        // Approx
        double sE=0,sEe=0,sEe2=0;
        for(int ie=0;ie<3;ie++){
            double ec=e2+(ie-1)*0.025;
            for(int ip=3;ip<8;ip++){double e=(*cl)[(ie+2)*11+ip];sE+=e;sEe+=e*ec;sEe2+=e*ec*ec;}
        }
        double wa=(sE>0&&sEe2/sE-std::pow(sEe/sE,2)>=0)?std::sqrt(sEe2/sE-std::pow(sEe/sE,2)):-999;
        
        // CellEta
        int sr=si/11,sc=si%11;
        int r0=std::max(0,sr-1),r1=std::min(6,sr+1),c0=std::max(0,sc-2),c1=std::min(10,sc+2);
        sE=sEe=sEe2=0;
        for(int r=r0;r<=r1;r++) for(int c=c0;c<=c1;c++){
            int k=r*11+c;double e=(*cl)[k],et=(*ce)[k];sE+=e;sEe+=e*et;sEe2+=e*et*et;
        }
        double wc=(sE>0&&sEe2/sE-std::pow(sEe/sE,2)>=0)?std::sqrt(sEe2/sE-std::pow(sEe/sE,2)):-999;
        
        // Correction
        double ecl=(*ce)[si],x=std::fabs(e2-ecl-0.0125),er=std::fmod(x,0.025)/0.025,co=0;
        if(ae<0.8){int j=(er<0.1)?0:(er<0.9?1:2);co=P0A[j]+P1A[j]*er+P2A[j]*er*er;}
        else if(ae<1.8){int j=(er<0.1)?0:(er<0.9?1:2);co=P0B[j]+P1B[j]*er+P2B[j]*er*er;}
        else if(ae<2.0){co=P0C[0]+P1C[0]*er+P2C[0]*er*er;}
        else if(ae<2.5){co=P0C[1]+P1C[1]*er+P2C[1]*er*er;}
        double wr=(wc>0)?wc-co:-999;
        
        auto&s=st[b];
        if(wa>0){double d=wa-wu;s.sda+=d;s.sda2+=d*d;}
        if(wc>0){double d=wc-wu;s.sdc+=d;s.sdc2+=d*d;}
        if(wr>0){double d=wr-wu;s.sdr+=d;s.sdr2+=d*d;}
        if((i+1)%1000000==0) printf("  %lld/%lld\n",i+1,nev);
    }
    printf("|eta| bin   |    N      | %%off  | Approx mean / std    | CellEta mean / std   | Corrected mean / std\n");
    printf("--------------------------------------------------------------------------------------------------------------\n");
    for(int b=0;b<NB;b++){
        auto&s=st[b]; if(s.n<100) continue;
        double p=100.0*s.noff/s.n;
        double ma=s.sda/s.n,sa=std::sqrt(std::max(0.0,s.sda2/s.n-ma*ma));
        double mc=s.sdc/s.n,sc=std::sqrt(std::max(0.0,s.sdc2/s.n-mc*mc));
        double mr=s.sdr/s.n,sr=std::sqrt(std::max(0.0,s.sdr2/s.n-mr*mr));
        printf("%.2f-%.2f | %9lld | %5.1f | %+.6f / %.6f | %+.6f / %.6f | %+.6f / %.6f\n",
               ed[b],ed[b+1],s.n,p,ma,sa,mc,sc,mr,sr);
    }
}
""")

f = ROOT.TFile.Open('../ntuples/mc23e_v3/mc23e_v3_Zeeg.root')
tree = f.Get('tree')
ROOT.run_eta_analysis(tree, tree.GetEntries())
f.Close()
