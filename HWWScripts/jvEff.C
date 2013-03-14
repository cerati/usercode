#include <iostream>
#include "TFile.h"
#include "TMath.h" 
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "tdrStyle.C"

#include "histTools.C"

void jvEff(float lumi=19.5) {

  gROOT->Reset();
  gStyle->SetOptStat(0);

  setTDRStyle();

  TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);

  TString cut = "((cuts & 4718595)==4718595)&&lep1.pt()>20.0&&lep2.pt()>20.0&&(dstype!=0 || (cuts & 1073741824)==1073741824)&&abs(dilep.mass()-91.2)<7.5&&type!=1&&type!=2";

  TFile* f_dy = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/dyll.root");
  TFile* f_da = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/data.root");

  TTree* dy = (TTree*) f_dy->Get("tree");
  TTree* da = (TTree*) f_da->Get("tree");

  TH1F* nj_dy = new TH1F("nj_dy","nj_dy",6,-0.5,5.5);
  TH1F* nj_da = new TH1F("nj_da","nj_da",6,-0.5,5.5);
  dy->Draw("njets>>nj_dy",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s)",lumi,cut.Data()),"g");
  da->Draw("njets>>nj_da",Form("(%s)",cut.Data()),"g");
  nj_dy->SetLineColor(kBlack);
  nj_dy->SetLineWidth(2);
  nj_dy->SetMarkerColor(kBlack);
  nj_dy->SetMarkerStyle(21);
  nj_dy->SetTitle("");
  nj_dy->GetYaxis()->SetTitle("Fraction Of Events");
  nj_dy->GetYaxis()->SetTitleOffset(1.1);
  nj_dy->GetXaxis()->SetTitle("N_{jets}");
  nj_dy->GetXaxis()->SetNdivisions(505);
  nj_da->SetLineColor(kRed);
  nj_da->SetLineWidth(2);
  nj_da->SetMarkerColor(kRed);
  nj_da->SetMarkerStyle(20);
  nj_dy->DrawNormalized("HIST");
  nj_dy->Sumw2();
  nj_da->Sumw2();
  nj_da->DrawNormalized("PEX0,same");
  TLegend* leg_nj = new TLegend(0.55,0.55,0.80,0.80);
  leg_nj->SetTextFont(42);
  leg_nj->SetFillColor(kWhite);
  leg_nj->SetNColumns(1);
  leg_nj->SetLineWidth(0);
  leg_nj->SetLineColor(kWhite);
  leg_nj->SetShadowColor(kWhite);
  leg_nj->AddEntry(nj_dy,"Z MC","L");
  leg_nj->AddEntry(nj_da,"Z Data","P");
  leg_nj->Draw();
  c1->RedrawAxis();
  c1->SaveAs("Znjets.pdf");

  TH1F* pt_dy = new TH1F("pt_dy","pt_dy",12,0,180);
  TH1F* pt_da = new TH1F("pt_da","pt_da",12,0,180);
  dy->Draw("jet1.pt()>>pt_dy",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s)",lumi,cut.Data()),"g");
  da->Draw("jet1.pt()>>pt_da",Form("(%s)",cut.Data()),"g");
  overFlowInLastBin(pt_dy);
  overFlowInLastBin(pt_da);
  pt_dy->SetLineColor(kBlack);
  pt_dy->SetLineWidth(2);
  pt_dy->SetMarkerColor(kBlack);
  pt_dy->SetMarkerStyle(21);
  pt_dy->SetTitle("");
  pt_dy->GetYaxis()->SetTitle("Events/15 GeV");
  pt_dy->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
  pt_dy->GetYaxis()->SetNdivisions(505);
  pt_dy->GetXaxis()->SetNdivisions(505);
  pt_da->SetLineColor(kRed);
  pt_da->SetLineWidth(2);
  pt_da->SetMarkerColor(kRed);
  pt_da->SetMarkerStyle(20);
  pt_dy->Draw("HIST");
  pt_dy->Sumw2();
  pt_da->Sumw2();
  pt_da->Draw("PEX0,same");
  TLegend* leg_pt = new TLegend(0.55,0.55,0.80,0.80);
  leg_pt->SetTextFont(42);
  leg_pt->SetFillColor(kWhite);
  leg_pt->SetNColumns(1);
  leg_pt->SetLineWidth(0);
  leg_pt->SetLineColor(kWhite);
  leg_pt->SetShadowColor(kWhite);
  leg_pt->AddEntry(pt_dy,"Z MC","L");
  leg_pt->AddEntry(pt_da,"Z Data","P");
  leg_pt->Draw();
  //c1->SetLogy();
  c1->RedrawAxis();
  c1->SaveAs("ZjetpT.pdf");

  TH1F* pt_e_dy = new TH1F("pt_e_dy","pt_e_dy",12,0,180);
  TH1F* pt_e_da = new TH1F("pt_e_da","pt_e_da",12,0,180);
  makeVetoEfficHisto(pt_dy,pt_e_dy);
  makeVetoEfficHisto(pt_da,pt_e_da);
  pt_e_dy->SetLineColor(kBlack);
  pt_e_dy->SetLineWidth(2);
  pt_e_dy->SetMarkerColor(kBlack);
  pt_e_dy->SetMarkerStyle(21);
  pt_e_dy->SetTitle("");
  pt_e_dy->GetYaxis()->SetTitle("Jet Veto Efficiency");
  pt_e_dy->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
  pt_e_dy->GetYaxis()->SetNdivisions(505);
  pt_e_dy->GetXaxis()->SetNdivisions(505);
  pt_e_dy->GetYaxis()->SetRangeUser(0.,1.);
  pt_e_da->SetLineColor(kRed);
  pt_e_da->SetLineWidth(2);
  pt_e_da->SetMarkerColor(kRed);
  pt_e_da->SetMarkerStyle(20);
  pt_e_dy->Draw("HIST");
  pt_e_da->Draw("PEX0,same");
  TLegend* leg_pt_e = new TLegend(0.55,0.55,0.80,0.80);
  leg_pt_e->SetTextFont(42);
  leg_pt_e->SetFillColor(kWhite);
  leg_pt_e->SetNColumns(1);
  leg_pt_e->SetLineWidth(0);
  leg_pt_e->SetLineColor(kWhite);
  leg_pt_e->SetShadowColor(kWhite);
  leg_pt_e->AddEntry(pt_e_dy,"Z MC","L");
  leg_pt_e->AddEntry(pt_e_da,"Z Data","P");
  leg_pt_e->Draw();
  c1->RedrawAxis();
  c1->SaveAs("Zjetvetoeff.pdf");

  TH1F* nv_dy = new TH1F("nv_dy","nv_dy",9,-0.5,44.5);
  TH1F* nv_da = new TH1F("nv_da","nv_da",9,-0.5,44.5);
  dy->Draw("nvtx>>nv_dy",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s)",lumi,cut.Data()),"g");
  da->Draw("nvtx>>nv_da",Form("(%s)",cut.Data()),"g");
  overFlowInLastBin(nv_dy);
  overFlowInLastBin(nv_da);
  nv_dy->SetLineColor(kBlack);
  nv_dy->SetLineWidth(2);
  nv_dy->SetMarkerColor(kBlack);
  nv_dy->SetMarkerStyle(21);
  nv_dy->SetTitle("");
  nv_dy->GetYaxis()->SetTitle("Number Of Events");
  nv_dy->GetXaxis()->SetTitle("N_{vtx}");
  nv_dy->GetYaxis()->SetNdivisions(505);
  nv_dy->GetXaxis()->SetNdivisions(505);
  nv_da->SetLineColor(kRed);
  nv_da->SetLineWidth(2);
  nv_da->SetMarkerColor(kRed);
  nv_da->SetMarkerStyle(20);
  nv_dy->Draw("HIST");
  nv_dy->Sumw2();
  nv_da->Sumw2();
  nv_da->Draw("PEX0,same");
  TLegend* leg_nv = new TLegend(0.55,0.55,0.80,0.80);
  leg_nv->SetTextFont(42);
  leg_nv->SetFillColor(kWhite);
  leg_nv->SetNColumns(1);
  leg_nv->SetLineWidth(0);
  leg_nv->SetLineColor(kWhite);
  leg_nv->SetShadowColor(kWhite);
  leg_nv->AddEntry(nv_dy,"Z MC","L");
  leg_nv->AddEntry(nv_da,"Z Data","P");
  leg_nv->Draw();
  c1->RedrawAxis();
  c1->SaveAs("Zvtx.pdf");
 
  TH1F* nv0_dy = new TH1F("nv0_dy","nv0_dy",9,-0.5,44.5);
  TH1F* nv0_da = new TH1F("nv0_da","nv0_da",9,-0.5,44.5);
  dy->Draw("nvtx>>nv0_dy",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && njets==0)",lumi,cut.Data()),"g");
  da->Draw("nvtx>>nv0_da",Form("(%s && njets==0)",cut.Data()),"g");
  overFlowInLastBin(nv0_dy);
  overFlowInLastBin(nv0_da);
  divideHisto(nv0_dy,nv_dy);
  divideHisto(nv0_da,nv_da);
  nv0_dy->SetLineColor(kBlack);
  nv0_dy->SetLineWidth(2);
  nv0_dy->SetMarkerColor(kBlack);
  nv0_dy->SetMarkerStyle(21);
  nv0_dy->SetTitle("");
  nv0_dy->GetYaxis()->SetTitle("Jet Veto Efficiency");
  nv0_dy->GetXaxis()->SetTitle("N_{vtx}");
  nv0_dy->GetYaxis()->SetNdivisions(505);
  nv0_dy->GetXaxis()->SetNdivisions(505);
  nv0_dy->GetYaxis()->SetRangeUser(0.,1.);
  nv0_da->SetLineColor(kRed);
  nv0_da->SetLineWidth(2);
  nv0_da->SetMarkerColor(kRed);
  nv0_da->SetMarkerStyle(20);
  nv0_dy->Draw("HIST");
  nv0_da->Draw("PEX0,same");
  TLegend* leg_nv0 = new TLegend(0.55,0.45,0.80,0.70);
  leg_nv0->SetTextFont(42);
  leg_nv0->SetFillColor(kWhite);
  leg_nv0->SetNColumns(1);
  leg_nv0->SetLineWidth(0);
  leg_nv0->SetLineColor(kWhite);
  leg_nv0->SetShadowColor(kWhite);
  leg_nv0->AddEntry(nv0_dy,"Z MC","L");
  leg_nv0->AddEntry(nv0_da,"Z Data","P");
  leg_nv0->Draw();
  c1->RedrawAxis();
  c1->SaveAs("Zjetvetoeff_vs_nvtx.pdf");

  TH1F* nv1_dy = new TH1F("nv1_dy","nv1_dy",9,-0.5,44.5);
  TH1F* nv1_da = new TH1F("nv1_da","nv1_da",9,-0.5,44.5);
  dy->Draw("nvtx>>nv1_dy",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && njets==1)",lumi,cut.Data()),"g");
  da->Draw("nvtx>>nv1_da",Form("(%s && njets==1)",cut.Data()),"g");
  overFlowInLastBin(nv1_dy);
  overFlowInLastBin(nv1_da);
  divideHisto(nv1_dy,nv_dy);
  divideHisto(nv1_da,nv_da);
  nv1_dy->SetLineColor(kBlack);
  nv1_dy->SetLineWidth(2);
  nv1_dy->SetMarkerColor(kBlack);
  nv1_dy->SetMarkerStyle(21);
  nv1_dy->SetTitle("");
  nv1_dy->GetYaxis()->SetTitle("Fraction Of 1-Jet Events");
  nv1_dy->GetYaxis()->SetRangeUser(0.,1.);
  nv1_dy->GetXaxis()->SetTitle("N_{vtx}");
  nv1_dy->GetYaxis()->SetNdivisions(505);
  nv1_dy->GetXaxis()->SetNdivisions(505);
  nv1_da->SetLineColor(kRed);
  nv1_da->SetLineWidth(2);
  nv1_da->SetMarkerColor(kRed);
  nv1_da->SetMarkerStyle(20);
  nv1_dy->Draw("");
  nv1_da->Draw("PEX0,same");
  TLegend* leg_nv1 = new TLegend(0.55,0.55,0.80,0.80);
  leg_nv1->SetTextFont(42);
  leg_nv1->SetFillColor(kWhite);
  leg_nv1->SetNColumns(1);
  leg_nv1->SetLineWidth(0);
  leg_nv1->SetLineColor(kWhite);
  leg_nv1->SetShadowColor(kWhite);
  leg_nv1->AddEntry(nv1_dy,"Z MC","L");
  leg_nv1->AddEntry(nv1_da,"Z Data","P");
  leg_nv1->Draw();
  c1->RedrawAxis();
  c1->SaveAs("Zonejeteff_vs_nvtx.pdf");

}
