void plot_reta_rphi() {
    TFile* f = TFile::Open("ntuples_mc23e/user.femarta.700770.Zeeg_mc23e.v1_EXT0/user.femarta.48437820.EXT0._000001.root");
    TTree* t = (TTree*)f->Get("tree");

    vector<double>* cells = nullptr;
    t->SetBranchAddress("photon.7x11ClusterLr2E", &cells);

    // Create histograms
    TH2F* h_reta_2d = new TH2F("h_reta_2d", "Computed vs Unfudged R_{#eta};Unfudged R_{#eta};Computed R_{#eta}", 100, 0.5, 1.05, 100, 0.5, 1.05);
    TH2F* h_rphi_2d = new TH2F("h_rphi_2d", "Computed vs Unfudged R_{#phi};Unfudged R_{#phi};Computed R_{#phi}", 100, 0.5, 1.05, 100, 0.5, 1.05);
    TH1F* h_reta_diff = new TH1F("h_reta_diff", "R_{#eta}: Computed - Unfudged;#Delta R_{#eta};Events", 100, -0.02, 0.02);
    TH1F* h_rphi_diff = new TH1F("h_rphi_diff", "R_{#phi}: Computed - Unfudged;#Delta R_{#phi};Events", 100, -0.02, 0.02);

    int nentries = t->GetEntries();
    cout << "Processing " << nentries << " events..." << endl;

    for (int ev=0; ev<nentries; ev++) {
        t->GetEntry(ev);
        
        double reta_unfudged = t->GetLeaf("photon.unfudged_reta")->GetValue();
        double rphi_unfudged = t->GetLeaf("photon.unfudged_rphi")->GetValue();
        
        if (cells->size() != 77) continue;
        
        int neta=7, nphi=11;
        double E_7x7=0, E_3x7=0, E_3x3=0;
        for (int ie=0; ie<7; ie++) 
            for (int ip=2; ip<=8; ip++) E_7x7 += (*cells)[ie*nphi + ip];
        for (int ie=2; ie<=4; ie++) 
            for (int ip=2; ip<=8; ip++) E_3x7 += (*cells)[ie*nphi + ip];
        for (int ie=2; ie<=4; ie++) 
            for (int ip=4; ip<=6; ip++) E_3x3 += (*cells)[ie*nphi + ip];
        
        if (E_7x7 <= 0 || E_3x7 <= 0) continue;
        
        double reta_comp = E_3x7/E_7x7;
        double rphi_comp = E_3x3/E_3x7;
        
        if (reta_unfudged > 0.5 && reta_unfudged < 1.05 && reta_comp > 0.5 && reta_comp < 1.05) {
            h_reta_2d->Fill(reta_unfudged, reta_comp);
            h_reta_diff->Fill(reta_comp - reta_unfudged);
        }
        if (rphi_unfudged > 0.5 && rphi_unfudged < 1.05 && rphi_comp > 0.5 && rphi_comp < 1.05) {
            h_rphi_2d->Fill(rphi_unfudged, rphi_comp);
            h_rphi_diff->Fill(rphi_comp - rphi_unfudged);
        }
    }
    f->Close();

    cout << "Reta 2D entries: " << h_reta_2d->GetEntries() << endl;
    cout << "Rphi 2D entries: " << h_rphi_2d->GetEntries() << endl;

    // Plot
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "Reta Rphi Validation", 1200, 1000);
    c->Divide(2, 2);

    c->cd(1);
    gPad->SetLogz();
    h_reta_2d->Draw("COLZ");
    TLine* line1 = new TLine(0.5, 0.5, 1.05, 1.05);
    line1->SetLineColor(kRed);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw();

    c->cd(2);
    gPad->SetLogz();
    h_rphi_2d->Draw("COLZ");
    TLine* line2 = new TLine(0.5, 0.5, 1.05, 1.05);
    line2->SetLineColor(kRed);
    line2->SetLineStyle(2);
    line2->SetLineWidth(2);
    line2->Draw();

    c->cd(3);
    gPad->SetLogy();
    h_reta_diff->SetLineColor(kBlue);
    h_reta_diff->SetLineWidth(2);
    h_reta_diff->Draw();

    c->cd(4);
    gPad->SetLogy();
    h_rphi_diff->SetLineColor(kBlue);
    h_rphi_diff->SetLineWidth(2);
    h_rphi_diff->Draw();

    c->SaveAs("reta_rphi_validation.png");
    c->SaveAs("reta_rphi_validation.pdf");
    cout << "Saved: reta_rphi_validation.png and .pdf" << endl;
}
