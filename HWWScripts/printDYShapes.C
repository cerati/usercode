void printDYShapes(int nj,int mass){

  gROOT->Reset();
  gStyle->SetOptStat(0);
  TFile *_file = TFile::Open(Form("cards/%i/hwwsf_%ij.input.root",mass,nj));

  TH1F* nominal = (TH1F*) _file->Get("histo_Zjets");
  TH1F* up      = (TH1F*) _file->Get(Form("histo_Zjets_CMS_MVAZBounding_hwwsf_%ijUp",nj));
  TH1F* down    = (TH1F*) _file->Get(Form("histo_Zjets_CMS_MVAZBounding_hwwsf_%ijDown",nj));

  TH1F* mc      = (TH1F*) _file->Get(Form("histo_Zjets_CMS_MVAZBounding_hwwsf_%ijUpOld",nj));
  TH1F* of      = (TH1F*) _file->Get(Form("histo_Zjets_CMS_MVAZBounding_hwwsf_%ijUpOFHiMet",nj));

  TH1F* syst = nominal->Clone("syst");
  for (int bb=1;bb<=syst->GetNbinsX();++bb){
    syst->SetBinError(bb,sqrt(pow(2*fabs(nominal->GetBinContent(bb)-up->GetBinContent(bb)),2) + pow(nominal->GetBinError(bb),2)));
  }
  syst->SetFillColor(kBlack);
  syst->SetFillStyle(3005);

  TCanvas c1;

  float max1 = nominal->GetBinContent(nominal->GetMaximumBin());
  float max2 = up->GetBinContent(up->GetMaximumBin());
  float max3 = down->GetBinContent(down->GetMaximumBin());
  float max4 = mc->GetBinContent(mc->GetMaximumBin());
  float max5 = of->GetBinContent(of->GetMaximumBin());
  float max = TMath::Max(max1,TMath::Max(max2,TMath::Max(max3,TMath::Max(max4,max5))));

  mc->SetLineColor(kMagenta);
  of->SetLineColor(kRed);

//   up->SetLineColor(kRed);
//   down->SetLineColor(kBlue);

  nominal->GetYaxis()->SetRangeUser(-2,1.4*max);

  nominal->SetLineWidth(2.);
  mc->SetLineWidth(2.);
  of->SetLineWidth(2.);
  up->SetLineWidth(2.);
  down->SetLineWidth(2.);

  nominal->SetTitle(Form("DY shape mH=%i %i-jet",mass,nj));
  nominal->GetXaxis()->SetTitle("WW BDT");

  nominal->Draw("E1");
  syst->Draw("E2 same");
  nominal->Draw("E1 same");
  of->Draw("E1 same");
//   up->Draw("HIST same");
//   down->Draw("HIST same");

  TLegend* leg = new TLegend(0.63,0.70,0.85,0.88);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->AddEntry(nominal," Default+Stat Band","le");
  //leg->AddEntry(up," Syst Up","l");
  //leg->AddEntry(down," Syst Down","l");
  leg->AddEntry(syst," Syst+Stat Band","f");
  //leg->AddEntry(mc," MC High MET","le");
  leg->AddEntry(of," Data OF,VZ subtr.","le");
  leg->Draw();

  c1.SaveAs(Form("zjets_shape_mh%i_%ijet.png",mass,nj));

}
