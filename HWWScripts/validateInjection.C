{

  gROOT->Reset();
  gStyle->SetOptStat(0);

  TString dirtest = "cards_inj_stat_ww_alter";
  int nj = 0;

  TCanvas c1;

  TFile *_file_ref = TFile::Open(Form("cards_def/125/hwwof_%ij.input_8TeV.root",nj));
  THStack hs("hs","processes");
  TH1F* ZH = (TH1F*) _file_ref->Get("histo_ZH");
  TH1F* WH = (TH1F*) _file_ref->Get("histo_WH");
  TH1F* ggH = (TH1F*) _file_ref->Get("histo_ggH");
  TH1F* qqH = (TH1F*) _file_ref->Get("histo_qqH");
  TH1F* ggWW = (TH1F*) _file_ref->Get("histo_ggWW");
  TH1F* qqWW = (TH1F*) _file_ref->Get("histo_qqWW");
  TH1F* VV = (TH1F*) _file_ref->Get("histo_VV");
  TH1F* Top = (TH1F*) _file_ref->Get("histo_Top");
  TH1F* Wjets = (TH1F*) _file_ref->Get("histo_Wjets");
  TH1F* Zjets = (TH1F*) _file_ref->Get("histo_Zjets");
  TH1F* Wgamma = (TH1F*) _file_ref->Get("histo_Wgamma");
  ZH->SetLineColor(kRed);
  WH->SetLineColor(kRed);
  ggH->SetLineColor(kRed);
  qqH->SetLineColor(kRed);
  ggWW->SetLineColor(kAzure-9);
  qqWW->SetLineColor(kAzure-9);
  VV->SetLineColor(kAzure-2);
  Top->SetLineColor(kYellow);
  Wjets->SetLineColor(kGray+1);
  Zjets->SetLineColor(kGreen+2);
  Wgamma->SetLineColor(kGray+1);
  ZH->SetFillColor(kRed);
  WH->SetFillColor(kRed);
  ggH->SetFillColor(kRed);
  qqH->SetFillColor(kRed);
  ggWW->SetFillColor(kAzure-9);
  qqWW->SetFillColor(kAzure-9);
  VV->SetFillColor(kAzure-2);
  Top->SetFillColor(kYellow);
  Wjets->SetFillColor(kGray+1);
  Zjets->SetFillColor(kGreen+2);
  Wgamma->SetFillColor(kGray+1);
  ZH->SetFillStyle(1001);
  WH->SetFillStyle(1001);
  ggH->SetFillStyle(1001);
  qqH->SetFillStyle(1001);
  ggWW->SetFillStyle(1001);
  qqWW->SetFillStyle(1001);
  VV->SetFillStyle(1001);
  Top->SetFillStyle(1001);
  Wjets->SetFillStyle(1001);
  Zjets->SetFillStyle(1001);
  Wgamma->SetFillStyle(1001);
  hs.Add(ggWW);
  hs.Add(qqWW);
  hs.Add(VV);
  hs.Add(Top);
  hs.Add(Wjets);
  hs.Add(Wgamma);
  hs.Add(Zjets);
  hs.Add(ZH);
  hs.Add(WH);
  hs.Add(ggH);
  hs.Add(qqH);
  hs.Draw("hist");

  TH1F* hinj = new TH1F("hinj","hinj",80,-1.,1.);
  for (int i=0;i<1000;++i) {
    TFile *_file_inj = TFile::Open(dirtest+Form("/125/hwwof_%ij_shape_8TeV_PseudoData_sb.root",nj));
    TH1F* h = (TH1F*) _file_inj->Get(Form("j%iof_%i",nj,i));
    hinj->Add(h);
    _file_inj->Close();
  }
  hinj->Scale(1./1000.);
  hinj->SetLineWidth(2);
  hinj->SetTitle();
  hinj->GetXaxis()->SetTitle("unrolled");
  hinj->GetYaxis()->SetTitle("events");
  hinj->GetYaxis()->SetRangeUser(0.,1.1*TMath::Max(hinj->GetMaximum(),hs.GetMaximum()));
  hinj->Draw();

  hs.Draw("hist,same");
  hinj->Draw("same");

  c1.SaveAs("toys_"+dirtest+Form("_%i.png",nj));

//   TH1F* hnew = new TH1F("hnew","hnew",80,-1.,1.);
//   for (int i=0;i<100;++i) {    
//     TFile *_file_ref = TFile::Open(Form("testcards/cards_125_N%i/250/hwwof_0j.input_8TeV.root",i));
//     TH1F* h = (TH1F*) _file_ref->Get("histo_Data");
//     hnew->Add(h);
//     _file_ref->Close();
//   }
//   hnew->Scale(1./100.);
//   hnew->SetMarkerStyle(20);
//   hnew->Draw("p,same");

  /*
  TLegend* leg = new TLegend(0.15,0.91,0.85,0.99);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->AddEntry(&hs, "ggH mH=160 GeV","l");
  leg->AddEntry(hinj,"ggH mH=125 GeV","l");
  leg->Draw();

  c1.SaveAs("inj125_bdt160.png");
  */
}
