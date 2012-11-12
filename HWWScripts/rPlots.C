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
      if (njets<2) {
	mc->GetXaxis()->SetTitle("dymva bin");
	mc->GetXaxis()->SetBinLabel(1,"[-0.9,-0.85)");
	mc->GetXaxis()->SetBinLabel(2,"[-0.85,-0.6)");
	mc->GetXaxis()->SetBinLabel(3,"[-0.6,WP)");
	mc->GetXaxis()->SetBinLabel(4,"[WP,1)");
      } else {
	mc->GetXaxis()->SetTitle("min-MET bin");
	mc->GetXaxis()->SetBinLabel(1,"[20,25)");
	mc->GetXaxis()->SetBinLabel(2,"[25,30)");
	mc->GetXaxis()->SetBinLabel(3,"[30,45)");
	mc->GetXaxis()->SetBinLabel(4,"[45,#infty)");
      }
      mc->GetXaxis()->SetLabelSize(0.06);
      mc->GetYaxis()->SetRangeUser(0,1);
      mc->GetYaxis()->SetRangeUser(-0.5,3);
      mc->Draw("PE");

      TH1F* dd = f.Get(Form("data_mh%i_%ij",mass,njets));
      dd->SetLineColor(kBlue);
      dd->SetMarkerColor(kBlue);
      dd->SetMarkerStyle(21);
      dd->SetLineWidth(2);
      dd->SetBinContent(4,-999);
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
