void plotOFVZsub(TString var, int nb, float min, float max, bool norm, TString mycut, int njets, TString fs) {

  gROOT->Reset();
  gStyle->SetOptStat(0);

  TString dir = "/smurf/cerati/skims/Run2012_Summer12_SmurfV9_53X/test/skim_dy/";

  TChain *ph = new TChain("tree");
  ph->Add(dir+"./dyll.root");

  TFile *_da = TFile::Open(dir+"./data.root");
  TTree* da = (TTree*) _da->Get("tree");

  TFile *_zz = TFile::Open(dir+"./zz.root");
  TTree* zz = (TTree*) _zz->Get("tree");
  TFile *_wz = TFile::Open(dir+"./wz.root");
  TTree* wz = (TTree*) _wz->Get("tree");

  TFile *_hw = TFile::Open(dir+"./hww125.root");
  TTree* hw = (TTree*) _hw->Get("tree");

  TH1F* h_zz = new TH1F("h_zz","h_zz",nb,min,max);
  TH1F* h_wz = new TH1F("h_wz","h_wz",nb,min,max);
  TH1F* h_dy = new TH1F("h_dy","h_dy",nb,min,max);
  TH1F* h_da = new TH1F("h_da","h_da",nb,min,max);
  TH1F* h_of = new TH1F("h_of","h_of",nb,min,max);
  TH1F* h_hw = new TH1F("h_hw","h_hw",nb,min,max);

  h_dy->SetFillColor(kGreen);
  h_of->SetFillColor(kRed);
  h_zz->SetFillColor(kBlue);
  h_wz->SetFillColor(kBlue);

  h_da->SetMarkerStyle(20);
  h_da->SetLineWidth(2);

  h_hw->SetLineColor(kCyan);
  h_hw->SetLineWidth(2);

  float lumi = 11.9;
  float lumicorr = 1.04;
  float dysf = 1.0;
  float dyer = 1.0;
  if (njets==0) dysf = 10.75;
  if (njets==1) dysf =  7.68;
  if (njets==0) dyer = 0.11;
  if (njets==1) dyer = 0.11;

  TString Met20 = Form("((cuts & 4719111)==4719111)&&njets==%i&&lep1.pt()>20.&&lep2.pt()>10.&&(dstype!=0 || (cuts & 1073741824)==1073741824) && met>20  && dilep.mass()>12. && min(pmet,pTrackMet)>20. && mt>80. && dilep.pt()>45.",njets);

  TString cut = Met20;

  TString minmet = " && min(pmet,pTrackMet)>45.  && (jet1.pt()<15 || dPhiDiLepJet1<165.*TMath::Pi()/180. )"; 
  TString dymva = " && ((njets==0 && dymva>0.88) || (njets==1 && dymva>0.84))";

  if (mycut.Contains("Zp")) cut+="&& abs(dilep.mass()-91)<7.5";
  if (mycut.Contains("oZ")) cut+="&& abs(dilep.mass()-91)>15";

  if (mycut.Contains("ptll45")) {
    cut+="&& dilep.pt()>45.";
    Met20+="&& dilep.pt()>45.";
  }

  if (mycut.Contains("DyMva")) cut+=dymva;

  if (mycut.Contains("MetGt45")) cut+=minmet;
  if (mycut.Contains("MetLt45")) cut+="&& min(pmet,pTrackMet)<45.";

  if (mycut.Contains("mll70")) cut+="&& dilep.mass()<70.";

  if (mycut.Contains("metsig25")) cut+="&& met/sqrt(sumet)<2.5";

  if (mycut.Contains("Mva05")) cut+="&& dymva>0.5";
  if (mycut.Contains("Mva02")) cut+="&& dymva>0.2";

  if (mycut.Contains("HiPU")) cut+="&& nvtx>=10";
  if (mycut.Contains("LoPU")) cut+="&& nvtx<10";

  if (mycut.Contains("HWW125")) {
    cut+="&& lep1.pt()>23 &&  lep2.pt()>10 && dPhi<100.*TMath::Pi()/180. && mt>80 && mt<123"+dymva;
    if (njets==0) dysf = 8.40;
    if (njets==1) dysf = 1.0;
    if (njets==0) dyer = 2.87/8.40;
    if (njets==1) dyer = 1.0;
  }

  if (mycut.Contains("HWW145")) {
    cut+="&& lep1.pt()>25 &&  lep2.pt()>15 && dPhi<90.*TMath::Pi()/180. && mt>80 && mt<130"+dymva;
    if (njets==0) dysf = 7.42;
    if (njets==1) dysf = 3.1;
    if (njets==0) dyer = 1.71/7.42;
    if (njets==1) dyer = 0.5;
  }

  if (mycut.Contains("HWW150")) {
    cut+="&& lep1.pt()>27 &&  lep2.pt()>25 && dPhi<90.*TMath::Pi()/180. && mt>80 && mt<150"+dymva;
    if (njets==0) dysf = 15.78;
    if (njets==1) dysf = 2.8;
    if (njets==0) dyer = 4.20/15.78;
    if (njets==1) dyer = 0.5;
  }

  if (mycut.Contains("HWW160")) {
    cut+="&& lep1.pt()>30 &&  lep2.pt()>25 && dPhi<60.*TMath::Pi()/180. && mt>90 && mt<160"+dymva;
    if (njets==0) dysf = 11.85;
    if (njets==1) dysf = 3.3;
    if (njets==0) dyer = 4.65/11.85;
    if (njets==1) dyer = 0.4;
  }

  if (mycut.Contains("HWW170")) {
    cut+="&& lep1.pt()>34 &&  lep2.pt()>25 && dPhi<60.*TMath::Pi()/180. && mt>110 && mt<170"+dymva;
    if (njets==0) dysf = 5.50;
    if (njets==1) dysf = 3.8;
    if (njets==0) dyer = 3.62/5.50;
    if (njets==1) dyer = 0.4;
  }

  if (mycut.Contains("HWW180")) {
    cut+="&& lep1.pt()>36 &&  lep2.pt()>25 && dPhi<70.*TMath::Pi()/180. && mt>120 && mt<180"+dymva;
    if (njets==0) dysf = 1.0;
    if (njets==1) dysf = 5.0;
    if (njets==0) dyer = 1.0;
    if (njets==1) dyer = 0.4;
  }

  if (mycut.Contains("HWW190")) {
    cut+="&& lep1.pt()>38 &&  lep2.pt()>25 && dPhi<90.*TMath::Pi()/180. && mt>120 && mt<190"+dymva;
    if (njets==0) dysf = 11.78;
    if (njets==1) dysf = 5.2;
    if (njets==0) dyer = 4.12/11.78;
    if (njets==1) dyer = 0.4;
  }

  if (mycut.Contains("HWW200")) {
    cut+="&& lep1.pt()>40 &&  lep2.pt()>25 && dPhi<100.*TMath::Pi()/180. && mt>120 && mt<200"+dymva;
    if (njets==0) dysf = 10.18;
    if (njets==1) dysf = 4.9;
    if (njets==0) dyer = 3.10/10.18;
    if (njets==1) dyer = 0.4;
  }

  TString sf = "(type==0 || type==3)";
  TString of = "(type==1 || type==2)";
  if (fs=="mm") {
    sf = "type==0";
    of = "type==1";
  } else if (fs=="ee") {
    sf = "type==3";
    of = "type==2";
  }

  TString cutsf = "("+cut+"&&"+sf+")";
  TString cutof = "("+cut+"&&"+of+")";
  TString cutmc = "scale1fb*sfWeightPU*sfWeightTrig*sfWeightEff*"+cutsf;

  cout << cutmc << endl;

  TCanvas c1;
  //c1.SetLogy();

  zz->Draw(var+">>h_zz",cutmc); 
  wz->Draw(var+">>h_wz",cutmc); 
  ph->Draw(var+">>h_dy",cutmc); 
  da->Draw(var+">>h_da",cutsf); 
  da->Draw(var+">>h_of",cutof); 
  hw->Draw(var+">>h_hw","scale1fb*sfWeightPU*sfWeightTrig*sfWeightEff*("+Met20+"&&"+sf+")"); 

  cout << h_dy->GetEntries() << endl;

  //add
  h_wz->Scale(lumi);
  h_zz->Scale(lumi);
  h_dy->Scale(lumi*dysf);
  h_of->Scale(lumicorr);
  h_hw->Scale(lumi);

  THStack hs("hs","stack");
  hs.Add(h_of);
  hs.Add(h_zz);
  hs.Add(h_wz);
  hs.Add(h_dy);

  TH1F* herr = new TH1F("herr","herr",nb,min,max);
  for (int bin=1;bin<max+1;bin++) {
    herr->SetBinContent(bin, ((TH1*)(hs.GetStack()->Last()))->GetBinContent(bin) );
    herr->SetBinError(bin,sqrt( pow(0.15*h_wz->GetBinContent(bin),2) + pow(0.15*h_wz->GetBinContent(bin),2) + 
				h_of->GetBinContent(bin) + dyer*h_dy->GetBinContent(bin) ) );
  }
  TGraphErrors* gerr = new TGraphErrors(herr);
  gerr->SetFillColor(kBlack);
  gerr->SetFillStyle(3244);

  var.ReplaceAll(".","");
  var.ReplaceAll("(","");
  var.ReplaceAll(")","");

  h_da->GetYaxis()->SetRangeUser(0.05,1.2*TMath::Max(h_da->GetBinContent(h_da->GetMaximumBin()), 
					      ((TH1*)(hs.GetStack()->Last()))->GetBinContent(((TH1*)(hs.GetStack()->Last()))->GetMaximumBin()) ));//3*
  h_da->GetXaxis()->SetTitle(var);
  h_da->SetTitle("");

  if (!norm) {
    h_da->Draw("PE");
    hs.Draw("same");
    h_da->Draw("PE,same");
    //h_hw->Draw("same");
    gerr->Draw("2");
  } else {
    h_da->Sumw2();
    h_da->DrawNormalized("PE");
    hs.DrawNormalized("same");
    //h_hw->DrawNormalized("same");
  }

  TLegend* leg = new TLegend(0.1,0.91,0.9,0.96);
  /*if (var.Contains("dPhi")) {
    delete leg;
    leg = new TLegend(0.55,0.11,0.88,0.31);
    }*/
  leg->SetFillColor(kWhite);
  leg->SetNColumns(5);
  leg->SetLineWidth(0);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->AddEntry(h_da,"SF data","p");
  leg->AddEntry(h_dy,Form("%2.1f x DY MC",dysf),"f");
  leg->AddEntry(h_of,"OF data","f");
  leg->AddEntry(h_wz,"VZ MC","f");
  //leg->AddEntry(h_hw,"HWW120","l");
  leg->Draw();

  c1.RedrawAxis();

  if (var.Contains("/")) var.ReplaceAll("/","");
  var = var+Form("_%ij",njets);

  gSystem->Exec("mkdir -p "+mycut+Form("_%ij_",njets)+fs);
  if (!norm) c1.SaveAs(mycut+Form("_%ij_",njets)+fs+"/"+var+".png");
  else c1.SaveAs(mycut+Form("_%ij_",njets)+fs+"/"+var+"_norm.png");
}
