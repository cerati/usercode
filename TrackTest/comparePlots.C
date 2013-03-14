void comparePlots(TString histoName, int algo){

  gStyle->SetOptStat(0);

  TFile* f_00eta09 = TFile::Open(Form("trktest_hp/histo_algo%02d_00eta09.root",algo));
  TFile* f_09eta16 = TFile::Open(Form("trktest_hp/histo_algo%02d_09eta16.root",algo));
  TFile* f_16eta25 = TFile::Open(Form("trktest_hp/histo_algo%02d_16eta25.root",algo));

  TH1F* h_00eta09 = (TH1F*) f_00eta09->Get("demo/"+histoName);
  TH1F* h_09eta16 = (TH1F*) f_09eta16->Get("demo/"+histoName);
  TH1F* h_16eta25 = (TH1F*) f_16eta25->Get("demo/"+histoName);

  h_00eta09->SetMarkerColor(kRed);
  h_09eta16->SetMarkerColor(kBlue);
  h_16eta25->SetMarkerColor(kBlack);

  h_00eta09->SetLineColor(kRed);
  h_09eta16->SetLineColor(kBlue);
  h_16eta25->SetLineColor(kBlack);

  h_00eta09->SetMarkerStyle(20);
  h_09eta16->SetMarkerStyle(21);
  h_16eta25->SetMarkerStyle(22);

  TCanvas c;

  h_00eta09->Scale(1./h_00eta09->Integral());
  h_09eta16->Scale(1./h_09eta16->Integral());
  h_16eta25->Scale(1./h_16eta25->Integral());

  h_00eta09->SetTitle("");
  h_00eta09->GetXaxis()->SetTitle(histoName);
  h_00eta09->GetYaxis()->SetTitle("event fraction");
  h_00eta09->GetYaxis()->SetRangeUser(0,1.1*TMath::Max(h_00eta09->GetMaximum(),TMath::Max(h_09eta16->GetMaximum(),h_16eta25->GetMaximum())));

  h_00eta09->Draw("PH");
  h_09eta16->Draw("samePH");
  h_16eta25->Draw("samePH");

  TLegend* leg = new TLegend(0.15,0.92,0.85,0.96);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetNColumns(3);
  leg->AddEntry(h_00eta09,"0.0<|#eta|<0.9","PL");
  leg->AddEntry(h_09eta16,"0.9<|#eta|<1.6","PL");
  leg->AddEntry(h_16eta25,"1.6<|#eta|<2.5","PL");
  leg->Draw();

  c.SaveAs(Form("trktest_hp/algo%02d_",algo)+histoName+".png");

}

/*
root -b -q comparePlots.C\(\"nlayer_strip_3d\",10\)
for plot in {algo,isHighPurity,pT,eta,phi,dxy,dz,nlayer_strip_3d,nlayer_tot_3d,nlayer_strip,nlayer_pixel}; do for algo in {9,10}; do root -b -q comparePlots.C\(\"${plot}\",${algo}\); done; done
*/
