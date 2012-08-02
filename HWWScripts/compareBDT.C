{

  gStyle->SetOptStat(0);

  TCanvas c1;

  TFile *_file_160_160 = TFile::Open("/smurf/yygao/inputLimits/ana_v9_5098pb/160/hwwof_0j.input_8TeV.root");
  TH1F* ggh_160_160 = (TH1F*) _file_160_160->Get("histo_ggH");
  ggh_160_160->SetLineWidth(2);
  ggh_160_160->SetTitle();
  ggh_160_160->GetXaxis()->SetTitle("BDT output");
  ggh_160_160->Draw();


  TFile *_file_125_160 = TFile::Open("cards/160/hwwof_0j_injecthww125_5098ipb.root");
  TH1F* ggh_125_160 = (TH1F*) _file_125_160->Get("histo_ggH");
  ggh_125_160->SetLineColor(kRed);
  ggh_125_160->SetLineWidth(2);
  ggh_125_160->Draw("same");

  TLegend* leg = new TLegend(0.15,0.91,0.85,0.99);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->AddEntry(ggh_160_160,"ggH mH=160 GeV","l");
  leg->AddEntry(ggh_125_160,"ggH mH=125 GeV","l");
  leg->Draw();

  TPaveText * labelcms  = new TPaveText(0.11,0.75,0.38,0.88,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.03);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS, #sqrt{s} = 8 TeV, L_{int} = 5.1 fb^{-1}");
  labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu mH=160 GeV");
  labelcms->AddText("0-jet,OF shape based analyis");
  labelcms->SetBorderSize(0);
  labelcms->SetTextFont(62);
  labelcms->SetLineWidth(2);
  labelcms->Draw();

  c1.SaveAs("inj125_bdt160.png");

}
