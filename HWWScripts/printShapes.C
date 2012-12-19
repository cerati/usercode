void printShapes(int nj,int mass, TString nameN, TString nameA, TString nameF){

  gROOT->Reset();
  gStyle->SetOptStat(0);
  TFile *_file = TFile::Open(Form("cards_def/%i/hwwof_%ij.input_8TeV.root",mass,nj));

  TH1F* nominal = (TH1F*) _file->Get(nameN);
  TH1F* up      = (TH1F*) _file->Get(nameA+"Up");
  TH1F* down    = (TH1F*) _file->Get(nameA+"Down");

  TCanvas c1;

  float max1 = nominal->GetBinContent(nominal->GetMaximumBin());
  float max2 = up->GetBinContent(up->GetMaximumBin());
  float max3 = down->GetBinContent(down->GetMaximumBin());
  float max = TMath::Max(max1,TMath::Max(max2,max3));

  up->SetLineColor(kRed);
  down->SetLineColor(kBlue);

  nominal->GetYaxis()->SetRangeUser(-2,1.4*max);

  nominal->SetLineWidth(2.);
  up->SetLineWidth(2.);
  down->SetLineWidth(2.);

  nominal->SetTitle(Form("%s shape mH=%i %i-jet",nameN.ReplaceAll("histo_","").Data(),mass,nj));
  nominal->GetXaxis()->SetTitle("unrolled bins");

  nominal->Draw("E1");
  up->Draw("HIST same");
  down->Draw("HIST same");

  TLegend* leg = new TLegend(0.63,0.70,0.85,0.88);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->AddEntry(nominal," Default+Stat Band","le");
  leg->AddEntry(up," Syst Up","l");
  leg->AddEntry(down," Syst Down","l");
  leg->Draw();

  c1.SaveAs(Form("%s_shape_mh%i_%ijet.png",nameF.Data(),mass,nj));

}
