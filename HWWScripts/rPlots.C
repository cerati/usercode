{

  gROOT->Reset();
  gStyle->SetOptStat(0);

  TString label = "cut";

  TFile f("outRFile_"+label+".root");

  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  int masses[] = {0,115,120,125,130,140,145,150,160,170,180,190,200,250,300};
  int nmasses = sizeof(masses)/sizeof(int);

  TCanvas c;

  for (int j=0;j<njetbins;++j) {
    int njets = jetbins[j];
    for (int jj=0;jj<nmasses;++jj) {
      int mass = masses[jj];

      TH1F* mc = f.Get(Form("mc_mh%i_%ij",mass,njets));
      mc->SetLineColor(kRed);
      mc->SetMarkerColor(kRed);
      mc->SetMarkerStyle(22);
      mc->SetLineWidth(2);
      mc->SetTitle(Form("mh=%i njets=%i",mass,njets));
      mc->GetYaxis()->SetTitle("R(out/in)");
      if (mass>0&&mass<=140) mc->GetXaxis()->SetTitle("dymva bin");
      else mc->GetXaxis()->SetTitle("min-MET bin");
      mc->GetYaxis()->SetRangeUser(0,1);
      mc->GetYaxis()->SetRangeUser(-0.5,3);
      mc->Draw("PE");

      TH1F* dd = f.Get(Form("data_mh%i_%ij",mass,njets));
      dd->SetLineColor(kBlue);
      dd->SetMarkerColor(kBlue);
      dd->SetMarkerStyle(21);
      dd->SetLineWidth(2);
      dd->Draw("PE,SAME");

      TLegend* leg = new TLegend(0.5,0.92,0.9,1.0);
      leg->SetNColumns(2);
      leg->SetFillColor(kWhite);
      leg->SetLineColor(kWhite);
      leg->AddEntry(dd,"data", "PL"); 
      leg->AddEntry(mc,"MC", "PL"); 
      leg->Draw();

      c.SaveAs(Form("routin_mh%i_%ij_%s.png",mass,njets,label.Data()));

    }
  }



}
