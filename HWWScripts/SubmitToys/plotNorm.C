void plotNorm(TString proc,TString inj,int jet, int mH, TString dir=".") {

  //process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt
  //rate   0.783   2.685   0.980  77.831  1635.926  106.477  51.070  219.236   6.044  200.709  48.820   0.000

  gROOT->Reset();
  gSystem->Load("tdrStyle_C.so");
  setTDRStyle();

//   TString proc = "Wgamma";
//   TString inj  = "def";

  float input = 0;
  float maxx  = 0;
  float hcp_s = 0;
  float hcp_b = 0;
  //THIS IS FOR 0J mH=125
  if (jet==0 && mH==125) {
    if (proc=="ZH")          { input = 0.783   ; maxx=3.0; hcp_s = 0.679   ; hcp_b = 0.000   ;	 }
    else if (proc=="WH")     { input = 2.685   ; maxx=3.0; hcp_s = 2.276   ; hcp_b = 0.000   ;	 }
    else if (proc=="qqH")    { input = 0.980   ; maxx=3.0; hcp_s = 0.825   ; hcp_b = 0.000   ;	 }
    else if (proc=="ggH")    { input = 77.831  ; maxx=3.0; hcp_s = 65.418  ; hcp_b = 0.000   ;	 }
    else if (proc=="qqWW")   { input = 1635.926; maxx=0.5; hcp_s = 1654.087; hcp_b = 1685.300;	 }
    else if (proc=="ggWW")   { input = 106.477 ; maxx=1.0; hcp_s = 124.564 ; hcp_b = 127.410 ;	 }
    else if (proc=="VV")     { input = 51.070  ; maxx=2.0; hcp_s = 53.537  ; hcp_b = 54.083  ;	 }
    else if (proc=="Top")    { input = 219.236 ; maxx=1.0; hcp_s = 246.561 ; hcp_b = 234.307 ;	 }
    else if (proc=="Zjets")  { input = 6.044   ; maxx=2.0; hcp_s = 5.399   ; hcp_b = 5.015   ;	 }
    else if (proc=="Wjets")  { input = 200.709 ; maxx=2.0; hcp_s = 155.072 ; hcp_b = 182.066 ;	 }
    else if (proc=="Wgamma") { input = 48.820  ; maxx=3.0; hcp_s = 49.600  ; hcp_b = 64.858  ;   }
    else if (proc=="Ztt")    { input = 0.000   ; maxx=2.0; hcp_s = 0.0     ; hcp_b = 0.0     ;   }
    else return;
  }
  //THIS IS FOR 1J mH=125
  if (jet==1 && mH==125) {
    if (proc=="ZH")          { input = 0.871  ; maxx=3.0; ; hcp_s =   0.518; hcp_b =    0.000;}
    else if (proc=="WH")     { input = 2.975  ; maxx=3.0; ; hcp_s =   1.762; hcp_b =    0.000;}
    else if (proc=="qqH")    { input = 4.355  ; maxx=3.0; ; hcp_s =   2.504; hcp_b =    0.000;}
    else if (proc=="ggH")    { input = 37.086 ; maxx=3.0; ; hcp_s =  21.473; hcp_b =    0.000;}
    else if (proc=="qqWW")   { input = 556.323; maxx=0.5; ; hcp_s = 612.738; hcp_b =  625.489;}
    else if (proc=="ggWW")   { input = 35.480 ; maxx=1.0; ; hcp_s =  39.930; hcp_b =   41.105;}
    else if (proc=="VV")     { input = 50.720 ; maxx=2.0; ; hcp_s =  52.349; hcp_b =   52.503;}
    else if (proc=="Top")    { input = 709.708; maxx=1.0; ; hcp_s = 684.513; hcp_b =  679.790;}
    else if (proc=="Zjets")  { input = 13.123 ; maxx=2.0; ; hcp_s =  11.654; hcp_b =   12.336;}
    else if (proc=="Wjets")  { input = 162.199; maxx=2.0; ; hcp_s = 119.309; hcp_b =  133.393;}
    else if (proc=="Wgamma") { input = 19.220 ; maxx=3.0; ; hcp_s =   7.492; hcp_b =    8.220;}
    else if (proc=="Ztt")    { input = 0.000  ; maxx=2.0; ; hcp_s =     0.0; hcp_b =      0.0;}
    else return; 
  }
  //THIS IS FOR 0J mH=200
  if (jet==0 && mH==200) {
    if (proc=="ZH")          { input = 0.943   ; maxx=3.0; }
    else if (proc=="WH")     { input = 1.316   ; maxx=3.0; }
    else if (proc=="qqH")    { input = 5.467   ; maxx=3.0; }
    else if (proc=="ggH")    { input = 269.728 ; maxx=3.0; }
    else if (proc=="qqWW")   { input = 1635.926; maxx=0.5; }
    else if (proc=="ggWW")   { input = 106.477 ; maxx=1.0; }
    else if (proc=="VV")     { input = 51.070  ; maxx=2.0; }
    else if (proc=="Top")    { input = 219.236 ; maxx=1.0; }
    else if (proc=="Zjets")  { input = 6.044   ; maxx=2.0; }
    else if (proc=="Wjets")  { input = 200.709 ; maxx=2.0; }
    else if (proc=="Wgamma") { input = 48.820  ; maxx=3.0; }
    else if (proc=="Ztt")    { input = 0.000   ; maxx=2.0; }
    else return;
  }
  //THIS IS FOR 1J mH=200
  if (jet==1 && mH==200) {
    if (proc=="ZH")          { input = 0.938   ; maxx=3.0; }
    else if (proc=="WH")     { input = 3.603   ; maxx=3.0; }
    else if (proc=="qqH")    { input = 21.207  ; maxx=3.0; }
    else if (proc=="ggH")    { input = 143.661 ; maxx=3.0; }
    else if (proc=="qqWW")   { input = 556.323; maxx=0.5; }
    else if (proc=="ggWW")   { input = 35.480 ; maxx=1.0; }
    else if (proc=="VV")     { input = 50.720 ; maxx=2.0; }
    else if (proc=="Top")    { input = 709.708; maxx=1.0; }
    else if (proc=="Zjets")  { input = 13.123 ; maxx=2.0; }
    else if (proc=="Wjets")  { input = 162.199; maxx=2.0; }
    else if (proc=="Wgamma") { input = 19.220 ; maxx=3.0; }
    else if (proc=="Ztt")    { input = 0.000  ; maxx=2.0; }
    else return; 
  }
  //THIS IS FOR 0J mH=250
  if (jet==0 && mH==250) {
    if (proc=="ZH")          { input = 0.329   ; maxx=3.0; }
    else if (proc=="WH")     { input = 0.408   ; maxx=3.0; }
    else if (proc=="qqH")    { input = 3.566   ; maxx=3.0; }
    else if (proc=="ggH")    { input = 158.586 ; maxx=3.0; }
    else if (proc=="qqWW")   { input = 1635.926; maxx=0.5; }
    else if (proc=="ggWW")   { input = 106.477 ; maxx=1.0; }
    else if (proc=="VV")     { input = 51.070  ; maxx=2.0; }
    else if (proc=="Top")    { input = 219.236 ; maxx=1.0; }
    else if (proc=="Zjets")  { input = 6.044   ; maxx=2.0; }
    else if (proc=="Wjets")  { input = 200.709 ; maxx=2.0; }
    else if (proc=="Wgamma") { input = 48.820  ; maxx=3.0; }
    else if (proc=="Ztt")    { input = 0.000   ; maxx=2.0; }
    else return;
  }
  //THIS IS FOR 1J mH=250
  if (jet==1 && mH==250) {
    if (proc=="ZH")          { input = 0.385   ; maxx=3.0; }
    else if (proc=="WH")     { input = 1.282   ; maxx=3.0; }
    else if (proc=="qqH")    { input = 13.242  ; maxx=3.0; }
    else if (proc=="ggH")    { input = 92.900 ; maxx=3.0; }
    else if (proc=="qqWW")   { input = 556.323; maxx=0.5; }
    else if (proc=="ggWW")   { input = 35.480 ; maxx=1.0; }
    else if (proc=="VV")     { input = 50.720 ; maxx=2.0; }
    else if (proc=="Top")    { input = 709.708; maxx=1.0; }
    else if (proc=="Zjets")  { input = 13.123 ; maxx=2.0; }
    else if (proc=="Wjets")  { input = 162.199; maxx=2.0; }
    else if (proc=="Wgamma") { input = 19.220 ; maxx=3.0; }
    else if (proc=="Ztt")    { input = 0.000  ; maxx=2.0; }
    else return; 
  }

  gSystem->Exec(Form("grep -h %s %s/logsNorm/%i/logNorm_%s_%i_hwwof_%ij_shape_8TeV_*.log >& tmp.txt",proc.Data(),dir.Data(),mH,inj.Data(),mH,jet));
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  TCanvas c1;
  TH1F* h_proc_b = new TH1F(proc+"_bkg_fit",proc,50,-1.*maxx,maxx);
  TH1F* h_proc_s = new TH1F(proc+"_sig_fit",proc,50,-1.*maxx,maxx);

  ifstream in;
  in.open("tmp.txt");

  Float_t nbfit,nsfit;
  TString s1,s2;
  Int_t nlines = 0;
  
  cout << "input=" << input << endl;
  while (1) {
    in >> s1 >> s2 >> nsfit >> nbfit;
    if (nlines < 5) printf("s1=%s s2=%s nbfit=%.3f, nsfit=%.3f\n",s1.Data(),s2.Data(),nbfit,nsfit);
    if (!in.good()) break;
    //cout << (nsfit-input)/input << endl;
    h_proc_b->Fill( (nbfit-input)/input );
    h_proc_s->Fill( (nsfit-input)/input );
    nlines++;
  }
  printf(" found %d points\n",nlines);
  
  in.close();
  gSystem->Exec(Form("mkdir -p %s/plots",dir.Data()));
  gSystem->Exec("rm tmp.txt");

  h_proc_b->SetTitle("");
  h_proc_b->GetXaxis()->SetTitle("(N_{fit,b}-N_{in})/N_{in}");
  h_proc_b->GetYaxis()->SetTitle("toys/bin");
  h_proc_b->Draw();
  h_proc_b->Fit("gaus");
  TLine line_b((hcp_b-input)/input,0,(hcp_b-input)/input,h_proc_b->GetBinContent(h_proc_b->GetMaximumBin()));
  line_b.SetLineColor(kMagenta);
  line_b.SetLineWidth(2);
  line_b.Draw("same");
  c1.SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_bfit_%s.png",dir.Data(),inj.Data(),jet,mH,proc.Data()));

  h_proc_s->SetTitle("");
  h_proc_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_proc_s->GetYaxis()->SetTitle("toys/bin");
  h_proc_s->Draw();
  h_proc_s->Fit("gaus");
  TLine line_s((hcp_s-input)/input,0,(hcp_s-input)/input,h_proc_s->GetBinContent(h_proc_s->GetMaximumBin()));
  line_s.SetLineColor(kMagenta);
  line_s.SetLineWidth(2);
  line_s.Draw("same");
  gPad->Update();
  if (gPad->GetFrame()->GetY2()>1000) h_proc_s->GetYaxis()->SetRangeUser(0,1000);
  c1.SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_%s.png",dir.Data(),inj.Data(),jet,mH,proc.Data()));
  
  
}
 
/*
root -b -q plotNorm.C\(\"qqWW\",\"125\",0,125\)
for proc in ggH qqWW ggWW Top Wjets Wgamma; do root -b -q plotNorm.C\(\"${proc}\",\"125\",0,125\); done
for inj in def 125 200; do for proc in ggH qqWW ggWW Top Wjets Wgamma; do root -b -q plotNorm.C\(\"${proc}\",\"${inj}\",0,125\); done; done
for nj in {0,1}; do for proc in ggH qqWW ggWW Top Wjets Wgamma; do root -b -q plotNorm.C\(\"${proc}\",\"125\",${nj},125\); done; done
*/
