#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TKey.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
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

void signifPlotWithDataNew(TString mode,int mH) {
  
  gROOT->Reset();
  setTDRStyle();  //gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);  

  int n=22;
  if (mode=="allcomb") n=21;
  if (mode=="shape") n=13;
  double mass[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double sig[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s1l[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s1h[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s2l[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double s2h[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double x0e[]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double obs[] = {2.1,2.3,2.4,2.6,2.7,2.9,2.7,2.4,2.1,1.9,1.9,1.6,0.0};
  double exp[] = {0.6,1.2,2.1,3.2,4.5,6.9,9.9,14.8,12.9,9.7,6.7,5.2,2.9};

  int in = 0;
  ifstream inlimits(mH>0 ? Form("toysignifNew_%s_inj%i.txt",mode.Data(),mH) : Form("toysignifNew_%s_injdef.txt",mode.Data()));
  string line;
  if (inlimits.is_open()) {
    while ( inlimits.good() && in<n ) {
      getline (inlimits,line);
      //cout << line << endl;
      float imas = string_to_double(line.substr(0,3));
      float isig = string_to_double(line.substr(5,9));
      float is1l = string_to_double(line.substr(14,18));
      float is1h = string_to_double(line.substr(23,27));
      float is2l = string_to_double(line.substr(32,36));
      float is2h = string_to_double(line.substr(41,45));
      //cout << imas << " " << iobs << " " << iosd << " " << iexp << " [" << il1s << "," << ih1s << "] [" << il2s << "," << ih2s << "]" << endl;
      mass[in] = imas;
      sig[in]  = isig;
      s1l[in]  = is1l;
      s1h[in]  = is1h;
      s2l[in]  = is2l;
      s2h[in]  = is2h;
      ++in;
    }
    inlimits.close();
  }

  TCanvas c1;  
  c1.SetLogx();

  //----> ALL THIS IS TO GET PROPER AXIS
  Double_t xlow   = 99.999;
  Double_t xup    = 300.01; // try 400.01 or 150.01
  Int_t    nbinsx = Int_t(xup - xlow);
  TH1F* h2 = new TH1F("h2", "h2", nbinsx, xlow, xup);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->SetTitle("");
  h2->GetXaxis()->SetTitle("m_{H} [GeV]");
  h2->GetYaxis()->SetTitle("significance");
  //h2->GetXaxis()->SetRangeUser(100,610);
  h2->GetYaxis()->SetRangeUser(0.,20.);
  h2->GetYaxis()->SetTitleOffset(1);
  //if (mH==200) h2->GetYaxis()->SetRangeUser(0.,10.);
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

  TGraphAsymmErrors * gr_s2b = new TGraphAsymmErrors(in,mass,sig,x0e,x0e,s2l,s2h);
  gr_s2b->SetFillColor(kYellow);
  gr_s2b->SetLineWidth(2);
  gr_s2b->Draw("3");

  TGraphAsymmErrors * gr_inj = new TGraphAsymmErrors(in,mass,sig,x0e,x0e,s1l,s1h);
  gr_inj->SetFillColor(kGreen);
  gr_inj->SetMarkerColor(kBlack);
  gr_inj->SetLineWidth(2);
  gr_inj->SetMarkerStyle(20);
  gr_inj->Draw("3");
  gr_inj->Draw("XL");

  TGraph * gr_obs = new TGraph(in,mass,obs);
  gr_obs->SetLineColor(kBlue);
  gr_obs->SetMarkerColor(kBlue);
  gr_obs->SetLineWidth(2);
  gr_obs->SetMarkerStyle(20);
  gr_obs->Draw("PL");

  TGraph * gr_exp = new TGraph(in,mass,exp);
  gr_exp->SetLineColor(kBlack);
  gr_exp->SetMarkerColor(kBlack);
  gr_exp->SetLineWidth(2);
  gr_exp->SetLineStyle(2);
  gr_exp->Draw("L");

  TLegend* leg = new TLegend(0.60,0.70,0.89,0.93);
  leg->SetTextFont(22);
  leg->SetFillColor(kWhite);
  //leg->SetNColumns(3);
  leg->SetLineWidth(0);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->AddEntry(gr_exp,"Expected","L");
  leg->AddEntry(gr_obs,"Observed","PL");
  leg->AddEntry(gr_inj,Form("Injection m_{H}=%i GeV",mH),"L");
  leg->AddEntry(gr_inj,"Injection #pm 1#sigma","F");
  leg->AddEntry(gr_s2b,"Injection #pm 2#sigma","F");
  leg->SetTextFont(62);
  leg->Draw();
  
  TPaveText * labelcms  = new TPaveText(0.15,0.80,0.58,0.93,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.025);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary");
  if (mode=="allcomb")
    labelcms->AddText("L = 4.9 fb^{-1} (7 TeV) + L = 12.1 fb^{-1} (8 TeV)");
  else
    labelcms->AddText("#sqrt{s} = 8 TeV, L = 12.1 fb^{-1}");
  if (mode=="cut")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1/2-jet cut based");
  if (mode=="shape")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1-jet e#mu 2D based");
  if (mode=="combine")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1/2-jet 2D+cut based");
  if (mode=="allcomb")
    labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu 0/1/2-jet");
  //labelcms->AddText(Form("Signal Injection m_{H}=%i GeV",mH));
  labelcms->SetBorderSize(0);
  labelcms->SetTextFont(62);
  labelcms->SetLineWidth(2);
  labelcms->Draw();

  c1.RedrawAxis();
  a->Draw();
  gPad->Modified(); gPad->Update(); // make sure it's redrawn

  c1.SaveAs(mH>0 ? Form("signifNew_%s_inj%i_data.eps",mode.Data(),mH) : Form("signifNew_%s_injdef_data.eps",mode.Data()));
  
  return; 
 
}

/*
for mode in allcomb combine cut shape; do root -b -q signifPlotWithDataNew.C+\(\"${mode}\",125\); done
*/
