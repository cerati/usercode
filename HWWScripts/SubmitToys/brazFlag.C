#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TKey.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TFrame.h>
#include <TGaxis.h>
#include "tdrStyle.C"

#include <sstream>

double string_to_double( const std::string& s ) {
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return 0;
  return x;
} 

void brazFlag(TString mode,int mH) {
  
  gROOT->Reset();
  setTDRStyle();  //gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);  

  int n=21;
  //if (mode=="allcomb") n=21;
  double mass[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double inj[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double i1l[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double i1h[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double exp[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s1l[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s1h[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s2l[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s2h[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double x0e[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


  double obt[] =  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  double obcut[]  =  {6.97,3.39,2.00,1.43,1.14,0.91,0.70,0.51,0.26,0.23,0.37,0.61,0.86,1.28,1.10,0.74,0.72,0.67,0.76,0.95,1.23};
  double obcomb[] =  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  double oballb[] =  {5.25,2.76,1.68,1.15,0.85,0.69,0.56,0.36,0.21,0.21,0.26,0.36,0.40,0.44,0.46,0.30,0.27,0.30,0.33,0.62,0.81};

  double* obs = obt;
  if (mode=="allcomb") obs = oballb;
  if (mode=="combine") obs = obcomb;
  if (mode=="allcut")     obs = obcut;

  int in = 0;
  ifstream inlimits(mH>0 ? Form("toysignif_%s_inj%i.txt",mode.Data(),mH) : Form("toysignif_%s_injdef.txt",mode.Data()));
  string line;
  if (inlimits.is_open()) {
    while ( inlimits.good() && in<n ) {
      getline (inlimits,line);
      //cout << line << endl;
      TString myline(line);
      myline.ReplaceAll('[',"");
      myline.ReplaceAll(']',"");
      myline.ReplaceAll(',',' ');
      myline.ReplaceAll("-",' ');
      myline.ReplaceAll("+",' ');
      //cout << myline << endl;
      if ((*myline.Tokenize(' ')).GetEntries()<11) continue;
      float imas = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[0])->GetString().Data() );
      float iinj = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[1])->GetString().Data() );
      float ii1l = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[2])->GetString().Data() );
      float ii1h = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[3])->GetString().Data() );
      float iexp = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[4])->GetString().Data() );
      float il1s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[5])->GetString().Data() );
      float ih1s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[6])->GetString().Data() );
      float il2s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[7])->GetString().Data() );
      float ih2s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[8])->GetString().Data() );
      //cout << imas << " " << iinj << "-" << ii1l  << "+" << ii1h << " " << iexp << " [" << il1s << "," << ih1s << "] [" << il2s << "," << ih2s << "]" << endl;
      mass[in] = imas;
      inj[in]  = iinj;
      i1l[in]  = ii1l;
      i1h[in]  = ii1h;
      exp[in]  = iexp;
      s1l[in]  = iexp-il1s;
      s1h[in]  = ih1s-iexp;
      s2l[in]  = iexp-il2s;
      s2h[in]  = ih2s-iexp;
      ++in;
    }
    inlimits.close();
  }

  TCanvas c1;
  c1.SetLogx();

  //----> ALL THIS IS TO GET PROPER AXIS
  Double_t xlow   = 99.999;
  Double_t xup    = 600.01; // try 400.01 or 150.01
  Int_t    nbinsx = Int_t(xup - xlow);
  TH1F* h2 = new TH1F("h2", "h2", nbinsx, xlow, xup);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->SetTitle("");
  h2->GetXaxis()->SetTitle("m_{H} [GeV]");
  h2->GetYaxis()->SetTitle("95% C.L. Limit on #sigma/#sigma_{SM}");
  h2->GetYaxis()->SetTitleOffset(1);
  //h2->GetXaxis()->SetRangeUser(100,610);
  h2->GetYaxis()->SetRangeUser(0.03,100);
  h2->Draw();
  // http://root.cern.ch/root/html/TF1.html#TF1:TF1@1
  // http://root.cern.ch/root/html/TFormula.html#TFormula:Analyze
  TF1 *f_h2_log10_x_axis = 0;
  f_h2_log10_x_axis = new TF1("f_h2_log10_x_axis", // name
			      "log10(x)", // formula
			      h2->GetXaxis()->GetXmin(), // xmin
			      h2->GetXaxis()->GetXmax()); // xmax
  // http://root.cern.ch/root/html/TGaxis.html#TGaxis:TGaxis@3
  // http://root.cern.ch/root/html/TGaxis.html#TGaxis:PaintAxis
  TGaxis *a = new TGaxis(h2->GetXaxis()->GetXmin(), // xmin
          h2->GetYaxis()->GetXmin(), // ymin
          h2->GetXaxis()->GetXmax(), // xmax
          h2->GetYaxis()->GetXmin(), // ymax
          "f_h2_log10_x_axis", // funcname
          50206, // ndiv (try 100006 or 506, don't try 1006)
          "BS", // chopt (try "BS" or "UBS")
          0.0); // gridlength
  // a->SetTickSize(h2->GetTickLength("X")); // use "the same" size
  a->SetTickSize(1.5 * h2->GetTickLength("X")); // make it bigger
  h2->SetTickLength(0.0, "X"); // get rid of "original" ticks
  if (!(TString(a->GetOption())).Contains("U")) {
    a->SetLabelFont(h2->GetLabelFont("X")); // use "the same" font
    a->SetLabelSize(h2->GetLabelSize("X")); // use "the same" size
    a->SetLabelOffset(0.015); // use "the same" size
    h2->SetLabelSize(0.0, "X"); // get rid of "original" labels
  }
  a->Draw();
  gPad->Modified(); gPad->Update(); // make sure it's redrawn
  //<---- ALL THIS IS TO GET PROPER AXIS

  c1.SetLogy();
  
  TGraphAsymmErrors * gr_s2b = new TGraphAsymmErrors(in,mass,exp,x0e,x0e,s2l,s2h);
  gr_s2b->SetFillColor(kYellow);
  gr_s2b->SetTitle("2s band");
  gr_s2b->Draw("3");
  
  TGraphAsymmErrors * gr_s1b = new TGraphAsymmErrors(in,mass,exp,x0e,x0e,s1l,s1h);
  gr_s1b->SetTitle("1s band");
  gr_s1b->SetFillColor(kGreen);
  gr_s1b->Draw("3");

  TLine l(100,1,600.,1);
  l.SetLineColor(kBlack);
  l.SetLineWidth(2);
  l.Draw();  

  TGraph * gr_exp = new TGraph(in,mass,exp);
  gr_exp->SetTitle("exp lim");
  gr_exp->SetLineColor(kBlack);
  gr_exp->SetLineWidth(2);
  gr_exp->SetLineStyle(7);
  gr_exp->Draw("L");

  TGraphAsymmErrors * gr_inj = new TGraphAsymmErrors(in,mass,inj,x0e,x0e,i1l,i1h);
  gr_inj->SetTitle("inj band");
  gr_inj->SetFillColor(kRed);
  gr_inj->SetFillStyle(3356);
  gr_inj->SetLineColor(kRed);
  gr_inj->SetLineWidth(2);
  gr_inj->Draw("3");
  gr_inj->Draw("XL");

  TGraph * gr_obs = new TGraph(in,mass,obs);
  gr_obs->SetLineColor(kBlack);
  gr_obs->SetLineWidth(4);
  gr_obs->SetMarkerColor(kBlack);
  gr_obs->SetMarkerStyle(20);
  gr_obs->Draw("L");

  //TLegend* leg = new TLegend(0.60,0.35,0.89,0.65);
  //TLegend* leg = new TLegend(0.60,0.55,0.89,0.85);
  TLegend* leg = new TLegend(0.60,0.70,0.89,0.93);
  leg->SetTextFont(22);
  leg->SetFillColor(kWhite);
  //leg->SetNColumns(3);
  leg->SetLineWidth(0);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->AddEntry(gr_obs,"Observed","L");
  leg->AddEntry(gr_exp,"Median Expected","L");
  leg->AddEntry(gr_s1b,"Expected #pm 1#sigma","F");
  leg->AddEntry(gr_s2b,"Expected #pm 2#sigma","F");
  leg->AddEntry(gr_inj,Form("Injection m_{H}=%i GeV #pm 1#sigma",mH),"L");
  //leg->AddEntry(gr_inj,"Injection #pm 1#sigma","F");
  leg->SetTextFont(62);
  leg->Draw();
  
  //TPaveText * labelcms  = new TPaveText(0.40,0.68,0.88,0.88,"NDCBR");
  //TPaveText * labelcms  = new TPaveText(0.10,0.68,0.58,0.88,"NDCBR");
  TPaveText * labelcms  = new TPaveText(0.15,0.78,0.58,0.93,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.030);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary");
  if (mode=="allcomb" || mode=="allcut") {
    labelcms->AddText("#sqrt{s}=7 TeV, L = 4.9 fb^{-1}");
    labelcms->AddText("#sqrt{s}=8 TeV, L = 19.5 fb^{-1}");
  } else
    labelcms->AddText("#sqrt{s} = 8 TeV, L = 19.5 fb^{-1}");
  if (mode=="allcut"||mode=="cut")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1-jet cut based");
  if (mode=="shape")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1-jet e#mu 2D based");
  if (mode=="combine")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1-jet 2D+cut based");
  if (mode=="allcomb")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1-jet");
  //labelcms->AddText(Form("Signal Injection m_{H}=%i GeV",mH));
  labelcms->SetBorderSize(0);
  labelcms->SetTextFont(62);
  labelcms->SetLineWidth(2);
  labelcms->Draw();

  c1.RedrawAxis();

  c1.SaveAs(mH>0 ? Form("limit_%s_inj%i.eps",mode.Data(),mH) : Form("limit_%s_injdef.eps",mode.Data()));
  
  return; 
 
}

/*
root -b -q brazFlag.C+\(\"allcomb\",125\)
for mode in combine cut shape; do root -b -q brazFlag.C+\(\"${mode}\",125\); done
*/
