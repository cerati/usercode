void validateInjection(TString dirtest, int nj){

  gROOT->Reset();
  gStyle->SetOptStat(0);

  TCanvas c1;

  TFile *_file_ref = TFile::Open(Form("cards/125/hwwof_%ij.input_8TeV.root",nj));
  THStack hs("hs","processes");
  TH1F* ZH = (TH1F*) _file_ref->Get("histo_ZH");
  TH1F* WH = (TH1F*) _file_ref->Get("histo_WH");
  TH1F* ggH = (TH1F*) _file_ref->Get("histo_ggH");
  TH1F* qqH = (TH1F*) _file_ref->Get("histo_qqH");
  TH1F* ggWW = (TH1F*) _file_ref->Get("histo_ggWW");
  TH1F* qqWW = (TH1F*) _file_ref->Get("histo_qqWW");
  TH1F* VV = (TH1F*) _file_ref->Get("histo_VV");
  TH1F* Top = (TH1F*) _file_ref->Get("histo_Top");
  TH1F* WjetsE = (TH1F*) _file_ref->Get("histo_WjetsE");
  TH1F* WjetsM = (TH1F*) _file_ref->Get("histo_WjetsM");
  TH1F* Zjets = (TH1F*) _file_ref->Get("histo_Zjets");
  TH1F* Wgamma = (TH1F*) _file_ref->Get("histo_Wgamma");
  TH1F* Wg3l = (TH1F*) _file_ref->Get("histo_Wg3l");
  TH1F* Ztt = (TH1F*) _file_ref->Get("histo_Ztt");
  ZH->SetLineColor(kRed);
  WH->SetLineColor(kRed);
  ggH->SetLineColor(kRed);
  qqH->SetLineColor(kRed);
  ggWW->SetLineColor(kAzure-9);
  qqWW->SetLineColor(kAzure-9);
  VV->SetLineColor(kAzure-2);
  Top->SetLineColor(kYellow);
  WjetsE->SetLineColor(kGray+1);
  WjetsM->SetLineColor(kGray+1);
  Zjets->SetLineColor(kGreen+2);
  Ztt->SetLineColor(kGreen+2);
  Wgamma->SetLineColor(kGray+1);
  Wg3l->SetLineColor(kGray+1);
  ZH->SetFillColor(kRed);
  WH->SetFillColor(kRed);
  ggH->SetFillColor(kRed);
  qqH->SetFillColor(kRed);
  ggWW->SetFillColor(kAzure-9);
  qqWW->SetFillColor(kAzure-9);
  VV->SetFillColor(kAzure-2);
  Top->SetFillColor(kYellow);
  WjetsE->SetFillColor(kGray+1);
  WjetsM->SetFillColor(kGray+1);
  Zjets->SetFillColor(kGreen+2);
  Ztt->SetFillColor(kGreen+2);
  Wgamma->SetFillColor(kGray+1);
  Wg3l->SetFillColor(kGray+1);
  ZH->SetFillStyle(1001);
  WH->SetFillStyle(1001);
  ggH->SetFillStyle(1001);
  qqH->SetFillStyle(1001);
  ggWW->SetFillStyle(1001);
  qqWW->SetFillStyle(1001);
  VV->SetFillStyle(1001);
  Top->SetFillStyle(1001);
  WjetsE->SetFillStyle(1001);
  WjetsM->SetFillStyle(1001);
  Zjets->SetFillStyle(1001);
  Ztt->SetFillStyle(1001);
  Wgamma->SetFillStyle(1001);
  Wg3l->SetFillStyle(1001);
  hs.Add(ggWW);
  hs.Add(qqWW);
  hs.Add(VV);
  hs.Add(Top);
  hs.Add(WjetsE);
  hs.Add(WjetsM);
  hs.Add(Wgamma);
  hs.Add(Wg3l);
  hs.Add(Zjets);
  hs.Add(Ztt);
  hs.Add(ZH);
  hs.Add(WH);
  hs.Add(ggH);
  hs.Add(qqH);
  hs.Draw("hist");

  TFile *_file_inj = TFile::Open(dirtest+Form("/125/hwwof_%ij_shape_8TeV_PseudoData_sb.root",nj));
  TH1F* hinj = (TH1F*) _file_inj->Get(Form("j%iof_%i",nj,0));
  for (int i=1;i<1000;++i) {
    TH1F* h = (TH1F*) _file_inj->Get(Form("j%iof_%i",nj,i));
    hinj->Add(h);
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

  gSystem->Exec("mkdir -p "+dirtest+"/plots");
  c1.SaveAs(dirtest+"/plots/toys_"+dirtest+Form("_%ij.png",nj));
  _file_inj->Close();

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

/*
root -b -q validateInjection.C\(\"cards_inj_stat\",0\)
*/
