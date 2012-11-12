#include <TH1F.h>
#include <TCut.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TPad.h>
#include <TMath.h>
#include <TLegend.h>
#include <TTree.h>
#include <TChain.h>
#include <THStack.h>
#include <TFile.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TFrame.h>

#include "tdrStyle.C"
#include "common.C"

#include <iostream>
#include <vector>

void fillOverflowInLastBin(TH1F* h) {
  h->SetBinContent(h->GetNbinsX(), h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1) );
}

int getIndexForProcess(TString* mcs, int nMC, TString proc) {
  for (int j=0;j<nMC;++j) if (mcs[j]==proc) return j;
  return -1;
}

void addInQuadrature(TH1* hnew, TH1* hold, float scale=1) {
  for (int bin=1;bin<=hold->GetNbinsX();bin++) {
    hnew->SetBinContent(bin, hnew->GetBinContent(bin)+pow(hold->GetBinContent(bin)*scale,2) );
  }
}

void plotBaby(float lumi=3.553, int njets=0, int mass=0, TString fs="", bool dodata=1, bool useSF=true, bool logy=0){
  //lumi is in /fb

  gROOT->Reset();
  setTDRStyle();
  gStyle->SetOptStat(0);

  bool wjetsFromData=1;

  bool compareData = 0;

  bool doSignal = 1;
  if (mass==0) doSignal = 0;

  bool doRatio = 1;

  bool dyFromLowMet = 1;
  //if (mass<=1) dyFromLowMet = 0;

  bool withZpeak = 0;

  bool blindData = 0;

  bool doSyst = 1;

  bool doWWxsec = 0;

  bool topTagged = 0;
  bool sameSign = 0;

  TString extensions[] = {".png"};
  //TString extensions[] = {".png",".eps",".root",".C"};

  TString plot[]    = {
    "dilep.mass()"
//     "dPhi*180./TMath::Pi()",
//     "dilep.pt()",
//     "lep1.pt()",
//     "lep2.pt()",
//     "type",
//     "pmet",
//     "pTrackMet",
//     "mt",
//     "jet1.pt()",
//     "jet2.pt()",
//     "nvtx",
//     "dymva",
//     "jet1.eta()",
//     "jet2.eta()",
//     "jet1Btag"
//     "jet1ProbBtag",
//     "njets"
//     "abs(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px())))*180./TMath::Pi()",
//     "acos(cos(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px()))-dilep.phi()))*180./TMath::Pi()",
//     "dPhiDiLepJet1*180./TMath::Pi()",
//     "abs(jet1.eta()-jet2.eta())"
  };
  int nPlots = sizeof(plot)/sizeof(TString);

  TString mcs[] = {"qqww","ggww","dyll","ttbar_powheg","tw","wz","zz","www","wjets","wgamma","zgamma","wglll","hww125"};//
  int nMC = sizeof(mcs)/sizeof(TString);

  vector<int> colors;
  for (int j=0;j<nMC;++j) colors.push_back(kBlack);
  if (getIndexForProcess(mcs,nMC,"qqww")>=0)   colors[getIndexForProcess(mcs,nMC,"qqww")]   = kAzure-9;
  if (getIndexForProcess(mcs,nMC,"ggww")>=0)   colors[getIndexForProcess(mcs,nMC,"ggww")]   = kAzure-9;
  if (getIndexForProcess(mcs,nMC,"dyll")>=0)   colors[getIndexForProcess(mcs,nMC,"dyll")]   = kGreen+2;
  if (getIndexForProcess(mcs,nMC,"ttbar_powheg")>=0)  colors[getIndexForProcess(mcs,nMC,"ttbar_powheg")]  = kYellow;
  if (getIndexForProcess(mcs,nMC,"tw")>=0)     colors[getIndexForProcess(mcs,nMC,"tw")]     = kYellow;
  if (getIndexForProcess(mcs,nMC,"wz")>=0)     colors[getIndexForProcess(mcs,nMC,"wz")]     = kAzure-2;
  if (getIndexForProcess(mcs,nMC,"zz")>=0)     colors[getIndexForProcess(mcs,nMC,"zz")]     = kAzure-2;
  if (getIndexForProcess(mcs,nMC,"www")>=0)    colors[getIndexForProcess(mcs,nMC,"www")]     = kAzure-2;
  if (getIndexForProcess(mcs,nMC,"wjets")>=0)  colors[getIndexForProcess(mcs,nMC,"wjets")]  = kGray+1;
  if (getIndexForProcess(mcs,nMC,"wgamma")>=0) colors[getIndexForProcess(mcs,nMC,"wgamma")] = kGray+1;
  if (getIndexForProcess(mcs,nMC,"zgamma")>=0) colors[getIndexForProcess(mcs,nMC,"zgamma")] = kGray+1;
  if (getIndexForProcess(mcs,nMC,"wglll")>=0)  colors[getIndexForProcess(mcs,nMC,"wglll")]  = kGray+1;
  if (getIndexForProcess(mcs,nMC,"hww125")>=0) colors[getIndexForProcess(mcs,nMC,"hww125")] = kRed;

  TCut runrange("run>0");//Full dataset

  TString dir = "/smurf/cerati/skims/Run2012_Summer12_SmurfV9_53X/skim_topww/";
  if (withZpeak) dir+="../skim_dy/";
  if (sameSign) dir+="../skim_wj/";

  vector<float> sfs;
  for (int j=0;j<nMC;++j) sfs.push_back(1.0);
  vector<float> ks;
  for (int j=0;j<nMC;++j) ks.push_back(1.0);

  if (useSF) {

    if (getIndexForProcess(mcs,nMC,"qqww")>=0)  sfs[getIndexForProcess(mcs,nMC,"qqww")]  = WWBkgScaleFactorCutBased(max(115,mass),min(njets,1));
    if (getIndexForProcess(mcs,nMC,"ggww")>=0)  sfs[getIndexForProcess(mcs,nMC,"ggww")]  = WWBkgScaleFactorCutBased(max(115,mass),min(njets,1));
    if (getIndexForProcess(mcs,nMC,"dyll")>=0)  sfs[getIndexForProcess(mcs,nMC,"dyll")]  = DYBkgScaleFactor(mass,njets);
    if (getIndexForProcess(mcs,nMC,"ttbar_powheg")>=0) sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = TopBkgScaleFactor(njets);
    if (getIndexForProcess(mcs,nMC,"tw")>=0)    sfs[getIndexForProcess(mcs,nMC,"tw")]    = TopBkgScaleFactor(njets);
    if (getIndexForProcess(mcs,nMC,"wglll")>=0) sfs[getIndexForProcess(mcs,nMC,"wglll")] = WGstarScaleFactor();
    if (wjetsFromData==0 && getIndexForProcess(mcs,nMC,"wjets")>=0) sfs[getIndexForProcess(mcs,nMC,"wjets")] = WJetsMCScaleFactor();
    
    if (getIndexForProcess(mcs,nMC,"qqww")>=0)  ks[getIndexForProcess(mcs,nMC,"qqww")]  = WWBkgScaleFactorKappaCutBased(max(115,mass),min(njets,1));
    if (getIndexForProcess(mcs,nMC,"ggww")>=0)  ks[getIndexForProcess(mcs,nMC,"ggww")]  = WWBkgScaleFactorKappaCutBased(max(115,mass),min(njets,1));
    if (getIndexForProcess(mcs,nMC,"dyll")>=0)  ks[getIndexForProcess(mcs,nMC,"dyll")]  = DYBkgScaleFactorKappa(mass,njets);
    if (getIndexForProcess(mcs,nMC,"ttbar_powheg")>=0) ks[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = TopBkgScaleFactorKappa(njets);
    if (getIndexForProcess(mcs,nMC,"tw")>=0)    ks[getIndexForProcess(mcs,nMC,"tw")]    = TopBkgScaleFactorKappa(njets);
    if (getIndexForProcess(mcs,nMC,"wglll")>=0) ks[getIndexForProcess(mcs,nMC,"wglll")] = 1.+WGstarScaleFactorSyst();
    if (getIndexForProcess(mcs,nMC,"wjets")>=0) ks[getIndexForProcess(mcs,nMC,"wjets")] = 1.360;
    
    if (dyFromLowMet) {
      //put the yields by hand,
      if (doWWxsec) {
	if (mass==0&&njets==0) sfs[getIndexForProcess(mcs,nMC,"dyll")] = 73.6;
	if (mass==0&&njets==1) sfs[getIndexForProcess(mcs,nMC,"dyll")] = 30.3;
	if (mass==0&&njets==2) sfs[getIndexForProcess(mcs,nMC,"dyll")] = 359.2;
      } else {
	if (mass==0&&njets==0) sfs[getIndexForProcess(mcs,nMC,"dyll")] = 119.74;
	if (mass==0&&njets==1) sfs[getIndexForProcess(mcs,nMC,"dyll")] = 46.99;
	if (mass==0&&njets==2) sfs[getIndexForProcess(mcs,nMC,"dyll")] = 528.88;
      }
    }

    /*
    //need to set mass dependent scale factors by hand //fixme
    if (njets==2) {
      //top 2-jet 
      if (mass==115) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 7.1; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==120) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 4.9; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==125) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 5.2; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==130) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 4.7; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==140) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 4.1; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==150) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 4.1; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==160) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 4.4; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==300) {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 1.4; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
    }
    */
    if (doWWxsec) {
      if (getIndexForProcess(mcs,nMC,"qqww")>=0)  sfs[getIndexForProcess(mcs,nMC,"qqww")]  = WWBkgScaleFactorMVA(max(115,mass),min(njets,1));
      if (getIndexForProcess(mcs,nMC,"ggww")>=0)  sfs[getIndexForProcess(mcs,nMC,"ggww")]  = WWBkgScaleFactorMVA(max(115,mass),min(njets,1));
      if (mass==0&&njets==0)   {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 1.02; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==0&&njets==1)   {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 1.07; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
      if (mass==0&&njets==2)   {sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")] = 1.17; sfs[getIndexForProcess(mcs,nMC,"tw")]=sfs[getIndexForProcess(mcs,nMC,"ttbar_powheg")];}
    }

  }
    
  TString mh = Form("%i",mass); 
  TString nj = Form("%i",njets); 

  float lep1pt_val=0.,lep2pt_val=0.,dPhi_val=0.,mll_val=0.,mtL_val=0.,mtH_val=0.,himass_val=0.;
  getCutValues(mass,lep1pt_val,lep2pt_val,dPhi_val,mll_val,mtL_val,mtH_val,himass_val);

  TCut lep1pt = Form("lep1.pt()>%.1f",lep1pt_val);
  TCut lep2pt = Form("lep2.pt()>%.1f",lep2pt_val);
  TCut dPhi = Form("dPhi<TMath::Pi()*%.1f/180.",dPhi_val);
  TCut mll = Form("dilep.mass()<%.1f",mll_val);
  TCut mt = Form("mt>%.1f&&mt<%.1f",mtL_val,mtH_val);
  TCut himass = Form("dilep.mass()>%.1f",himass_val);

  TCut ptreg   = lep1pt && lep2pt;
  TCut massreg = lep1pt && lep2pt && mll;
  TCut mtreg   = lep1pt && lep2pt && mll && mt;
  TCut sigreg  = lep1pt && lep2pt && mll && mt && dPhi;

  TCut sideband = lep1pt && lep2pt && himass;
  TCut mtside   = lep1pt && lep2pt && himass && mt;
  TCut dphiside = lep1pt && lep2pt && himass && mt && dPhi;

  if (withZpeak) wwSelNoMetNoLep = wwSelNoMetNoLep&~ZVeto;
  if (topTagged) {
    wwSelNoMetNoLep = wwSelNoMetNoLep&~TopVeto;
    wwSelNoMetNoLep = wwSelNoMetNoLep|TopTag;
  }
  if (sameSign) {
    wwSelNoMetNoLep = wwSelNoMetNoLep&~ChargeMatch;
  }

  TCut base(Form("(cuts & %i)==%i",wwSelNoMetNoLep,wwSelNoMetNoLep));
  TCut leps(Form("(cuts & %i)==%i && (cuts & %i)==%i",Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection));
  TCut lplusf(Form("(((cuts & %i)==%i && (cuts & %i)!=%i)||((cuts & %i)!=%i && (cuts & %i)==%i))",
		   Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection,Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection));
  TCut notTagNotInJets(Form("(cuts & %i)!=%i",TopTagNotInJets,TopTagNotInJets));
  TCut trig(Form("dstype!=0 || (cuts & %i)==%i",Trigger,Trigger));
  TCut newcuts = "";//"type==1 || type==2 || ( min(pmet,pTrackMet)>45. && (jet1.pt()<15 || dPhiDiLepJet1*180./TMath::Pi()<165.) )";
  if (njets==2) newcuts = "type==1 || type==2 || (met>45. && acos(cos( atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px())) - dilep.phi()))*180./TMath::Pi()<165.)";
  if (njets==0) newcuts = "type==1 || type==2 || dymva>0.88";
  if (njets==1) newcuts = "type==1 || type==2 || dymva>0.84";
  TCut njcut(Form("njets==%i",njets));
  if (njets==2) njcut = "(njets==2 || (njets==3 && !((jet1.eta()-jet3.eta() > 0 && jet2.eta()-jet3.eta() < 0) || (jet2.eta()-jet3.eta() > 0 && jet1.eta()-jet3.eta() < 0)) ))&&!(TMath::Abs(jet1.eta())>= 4.7||TMath::Abs(jet2.eta()) >= 4.7)";
  if (njets==-1) njcut = "njets==0 || njets==1 || (njets==2 || (njets==3 && !((jet1.eta()-jet3.eta() > 0 && jet2.eta()-jet3.eta() < 0) || (jet2.eta()-jet3.eta() > 0 && jet1.eta()-jet3.eta() < 0)) ))&&!(TMath::Abs(jet1.eta())>= 4.5||TMath::Abs(jet2.eta()) >= 4.5)";
  TCut kincuts = "dilep.pt()>45.";
  if (njets==2 && mass>0) kincuts = "dilep.pt()>45. && TMath::Abs(jet1.eta()-jet2.eta())>3.5 && (((jet1.eta()-lep1.eta() > 0 && jet2.eta()-lep1.eta() < 0) || (jet2.eta()-lep1.eta() > 0 && jet1.eta()-lep1.eta() < 0)) && ((jet1.eta()-lep2.eta() > 0 && jet2.eta()-lep2.eta() < 0) || (jet2.eta()-lep2.eta() > 0 && jet1.eta()-lep2.eta() < 0))) && sqrt(2*jet1.pt()*jet2.pt()*(TMath::CosH(jet1.eta()-jet2.eta())-TMath::Cos(jet1.phi()-jet2.phi())))>500.";
  TCut blindDataCut = "(dstype!=0 || dilep.mass()>70)";
  if (blindData) kincuts = kincuts&&blindDataCut;
  TCut wwXSecCut = "lep2.pt()>20";
  if (doWWxsec) kincuts = kincuts&&wwXSecCut;
  TCut ssCut = "lq1*lq2>0";
  if (sameSign) kincuts = kincuts&&ssCut;

  TCut flav = "";
  TString flavstr = "";
  TCut sf = "type!=1 && type!=2";
  TCut of = "type!=0 && type!=3";
  TCut mm = "type==0";
  TCut me = "type==1";
  TCut em = "type==2";
  TCut ee = "type==3";
  if (fs=="of") {
    flav = of;
    flavstr = "_of";
  } else if (fs=="sf") {
    flav = sf;
    flavstr = "_sf";
  } else if (fs=="mm") {
    flav = mm;
    flavstr = "_mm";
  } else if (fs=="me") {
    flav = me;
    flavstr = "_me";
  } else if (fs=="em") {
    flav = em;
    flavstr = "_em";
  } else if (fs=="ee") {
    flav = ee;
    flavstr = "_ee";
  } else if (fs!=""){
    cout << "final state not supported: " << fs << endl;
    return;
  }

  TCut cut = base&&njcut&&newcuts&&sigreg&&trig&&kincuts&&flav;
  //TCut cut = base&&njcut&&newcuts&&dphiside&&trig&&kincuts&&flav;
  cut.SetName("mh"+mh+"_nj"+nj+flavstr);

  cout << "cut: " << cut.GetTitle() << endl;

  vector<TString> binning;
  for (int j=0;j<nPlots;++j) binning.push_back("");
  if (getIndexForProcess(plot,nPlots,"dilep.mass()")>=0) binning[getIndexForProcess(plot,nPlots,"dilep.mass()")] = "25,0.,250.";
  if (getIndexForProcess(plot,nPlots,"dPhi*180./TMath::Pi()")>=0) binning[getIndexForProcess(plot,nPlots,"dPhi*180./TMath::Pi()")] = "18,0.,180.";
  if (getIndexForProcess(plot,nPlots,"dPhiDiLepJet1*180./TMath::Pi()")>=0) binning[getIndexForProcess(plot,nPlots,"dPhiDiLepJet1*180./TMath::Pi()")] = "9,0.,180.";
  if (getIndexForProcess(plot,nPlots,"abs(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px())))*180./TMath::Pi()")>=0) 
    binning[getIndexForProcess(plot,nPlots,"abs(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px())))*180./TMath::Pi()")] = "9,0.,180.";
  if (getIndexForProcess(plot,nPlots,"acos(cos(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px()))-dilep.phi()))*180./TMath::Pi()")>=0) 
    binning[getIndexForProcess(plot,nPlots,"acos(cos(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px()))-dilep.phi()))*180./TMath::Pi()")] = "9,0.,180.";
  if (getIndexForProcess(plot,nPlots,"dilep.pt()")>=0) binning[getIndexForProcess(plot,nPlots,"dilep.pt()")] = "15,0.,150.";
  if (getIndexForProcess(plot,nPlots,"lep1.pt()")>=0) binning[getIndexForProcess(plot,nPlots,"lep1.pt()")] = "30,0.,150.";
  if (getIndexForProcess(plot,nPlots,"lep2.pt()")>=0) binning[getIndexForProcess(plot,nPlots,"lep2.pt()")] = "30,0.,150.";
  if (getIndexForProcess(plot,nPlots,"type")>=0) binning[getIndexForProcess(plot,nPlots,"type")] = "4,0,4";
  if (getIndexForProcess(plot,nPlots,"pmet")>=0) binning[getIndexForProcess(plot,nPlots,"pmet")] = "15,0.,150.";
  if (getIndexForProcess(plot,nPlots,"pTrackMet")>=0) binning[getIndexForProcess(plot,nPlots,"pTrackMet")] = "15,0.,150.";
  if (getIndexForProcess(plot,nPlots,"mt")>=0) binning[getIndexForProcess(plot,nPlots,"mt")] = "30,0.,300.";
  if (getIndexForProcess(plot,nPlots,"jet1.pt()")>=0) binning[getIndexForProcess(plot,nPlots,"jet1.pt()")] = njets!=0 ? "20,0.,200." : "6,0.,30.";
  if (getIndexForProcess(plot,nPlots,"jet2.pt()")>=0) binning[getIndexForProcess(plot,nPlots,"jet2.pt()")] = njets==2 ? "20,0.,200." : "6,0.,30.";
  if (getIndexForProcess(plot,nPlots,"jet1.eta()")>=0) binning[getIndexForProcess(plot,nPlots,"jet1.eta()")] = "10,-5.,5.";
  if (getIndexForProcess(plot,nPlots,"jet2.eta()")>=0) binning[getIndexForProcess(plot,nPlots,"jet2.eta()")] = "10,-5.,5.";
  if (getIndexForProcess(plot,nPlots,"jet1Btag")>=0) binning[getIndexForProcess(plot,nPlots,"jet1Btag")] = "24,-4,20";
  if (getIndexForProcess(plot,nPlots,"jet1ProbBtag")>=0) binning[getIndexForProcess(plot,nPlots,"jet1ProbBtag")] = "20,0,2";
  if (getIndexForProcess(plot,nPlots,"nvtx")>=0) binning[getIndexForProcess(plot,nPlots,"nvtx")] = "25,0,50";
  if (getIndexForProcess(plot,nPlots,"dymva")>=0) binning[getIndexForProcess(plot,nPlots,"dymva")] = "21,-1,1.1";
  if (getIndexForProcess(plot,nPlots,"abs(jet1.eta()-jet2.eta())")>=0) binning[getIndexForProcess(plot,nPlots,"abs(jet1.eta()-jet2.eta())")] = "10,0,5.0";
  if (getIndexForProcess(plot,nPlots,"njets")>=0) binning[getIndexForProcess(plot,nPlots,"njets")] = "11,-0.5,10.5";

  vector<TString> xtitle;
  for (int j=0;j<nPlots;++j) xtitle.push_back("");
  if (getIndexForProcess(plot,nPlots,"dilep.mass()")>=0) xtitle[getIndexForProcess(plot,nPlots,"dilep.mass()")] = "m_{l,l} [GeV]";
  if (getIndexForProcess(plot,nPlots,"dPhi*180./TMath::Pi()")>=0) xtitle[getIndexForProcess(plot,nPlots,"dPhi*180./TMath::Pi()")] = "#Delta#phi_{l,l} [#circ]";
  if (getIndexForProcess(plot,nPlots,"dPhiDiLepJet1*180./TMath::Pi()")>=0) xtitle[getIndexForProcess(plot,nPlots,"dPhiDiLepJet1*180./TMath::Pi()")] = "#Delta#phi_{ll,j1} [#circ]";
  if (getIndexForProcess(plot,nPlots,"abs(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px())))*180./TMath::Pi()")>=0) 
    xtitle[getIndexForProcess(plot,nPlots,"abs(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px())))*180./TMath::Pi()")] = "#Delta#phi_{j1,j2} [#circ]";
  if (getIndexForProcess(plot,nPlots,"acos(cos(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px()))-dilep.phi()))*180./TMath::Pi()")>=0) 
    xtitle[getIndexForProcess(plot,nPlots,"acos(cos(atan2((jet1.py()+jet2.py()),(jet1.px()+jet2.px()))-dilep.phi()))*180./TMath::Pi()")] = "#Delta#phi_{ll,jj} [#circ]";
  if (getIndexForProcess(plot,nPlots,"dilep.pt()")>=0) xtitle[getIndexForProcess(plot,nPlots,"dilep.pt()")] = "pT_{l,l} [GeV]";
  if (getIndexForProcess(plot,nPlots,"lep1.pt()")>=0) xtitle[getIndexForProcess(plot,nPlots,"lep1.pt()")] = "pT_{l1} [GeV]";
  if (getIndexForProcess(plot,nPlots,"lep2.pt()")>=0) xtitle[getIndexForProcess(plot,nPlots,"lep2.pt()")] = "pT_{l2} [GeV]";
  if (getIndexForProcess(plot,nPlots,"type")>=0) xtitle[getIndexForProcess(plot,nPlots,"type")] = "type";
  if (getIndexForProcess(plot,nPlots,"pmet")>=0) xtitle[getIndexForProcess(plot,nPlots,"pmet")] = "projected #slash{E}_{T} [GeV]";
  if (getIndexForProcess(plot,nPlots,"pTrackMet")>=0) xtitle[getIndexForProcess(plot,nPlots,"pTrackMet")] = "projected track #slash{E}_{T} [GeV]";
  if (getIndexForProcess(plot,nPlots,"mt")>=0) xtitle[getIndexForProcess(plot,nPlots,"mt")] = "m_{T} [GeV]";
  if (getIndexForProcess(plot,nPlots,"jet1.pt()")>=0) xtitle[getIndexForProcess(plot,nPlots,"jet1.pt()")] = "pT_{j1} [GeV]";
  if (getIndexForProcess(plot,nPlots,"jet2.pt()")>=0) xtitle[getIndexForProcess(plot,nPlots,"jet2.pt()")] = "pT_{j2} [GeV]";
  if (getIndexForProcess(plot,nPlots,"jet1.eta()")>=0) xtitle[getIndexForProcess(plot,nPlots,"jet1.eta()")] = "#eta_{j1}";
  if (getIndexForProcess(plot,nPlots,"jet2.eta()")>=0) xtitle[getIndexForProcess(plot,nPlots,"jet2.eta()")] = "#eta_{j2}";
  if (getIndexForProcess(plot,nPlots,"jet1Btag")>=0) xtitle[getIndexForProcess(plot,nPlots,"jet1Btag")] = "TCHE discriminator";
  if (getIndexForProcess(plot,nPlots,"jet1ProbBtag")>=0) xtitle[getIndexForProcess(plot,nPlots,"jet1ProbBtag")] = "JetProb. discriminator";
  if (getIndexForProcess(plot,nPlots,"nvtx")>=0) xtitle[getIndexForProcess(plot,nPlots,"nvtx")] = "N_{vtx}";
  if (getIndexForProcess(plot,nPlots,"dymva")>=0) xtitle[getIndexForProcess(plot,nPlots,"dymva")] = "dymva output";
  if (getIndexForProcess(plot,nPlots,"abs(jet1.eta()-jet2.eta())")>=0) xtitle[getIndexForProcess(plot,nPlots,"abs(jet1.eta()-jet2.eta())")] = "#Delta#eta(j1.j2)";
  if (getIndexForProcess(plot,nPlots,"njets")>=0) xtitle[getIndexForProcess(plot,nPlots,"njets")] = "N_{jets}";

  //if we want to plot signal
  TFile *_signal = 0;
  TTree * signal = 0;
  if (doSignal) {
    _signal = TFile::Open(dir+"/hww"+mh+".root");
    signal = (TTree*) _signal->Get("tree");
  }

  TFile *_data = TFile::Open(dir+"/data.root");
  TTree * data = (TTree*) _data->Get("tree");

  TFile *_spill = TFile::Open(dir+"../skim_wj/data_spill.root");
  TTree * spill = (TTree*) _spill->Get("tree");
  TFile *_fakes = TFile::Open(dir+"../skim_wj/data.root");
  TTree * fakes = (TTree*) _fakes->Get("tree");

  bool allleg = false;

  for (int pl=0;pl<nPlots;++pl) {
    //cout << plot[pl] << endl;
    TString bin = binning[pl];
    int firstc = bin.First(',');
    int lastc  = bin.Last(',');
    TString nbins = bin;
    nbins.Remove(firstc,bin.Sizeof()-firstc);
    TString maxbin  = bin;
    maxbin.Remove(0,lastc+1);
    TString minbin  = bin;
    minbin.Remove(0,firstc+1);
    minbin.Remove(lastc-firstc-1,bin.Sizeof()-lastc);

    //do N-1 plots
    float x1=-10000,x2=+10000;
    if (mass>0) {
      cut = base&&njcut&&newcuts&&sigreg&&trig&&kincuts&&flav;
      if (plot[pl]=="dilep.mass()") cut = TCut(TString(cut.GetTitle()).ReplaceAll(mll,"dilep.mass()<999"));
      if (plot[pl]=="dPhi*180./TMath::Pi()") cut = TCut(TString(cut.GetTitle()).ReplaceAll(dPhi,"dPhi*180./TMath::Pi()<180."));
      if (plot[pl]=="mt") cut = TCut(TString(cut.GetTitle()).ReplaceAll(mt,"mt>80&&mt<999"));
      if (plot[pl]=="lep1.pt()") cut = TCut(TString(cut.GetTitle()).ReplaceAll(lep1pt,"lep1.pt()>20."));
      if (plot[pl]=="lep2.pt()") cut = TCut(TString(cut.GetTitle()).ReplaceAll(lep2pt,"lep2.pt()>10."));
      cut.SetName("mh"+mh+"_nj"+nj+flavstr);

      if (plot[pl]=="dilep.mass()") {
	x1 = 12.;
	x2 = TString(mll.GetTitle()).ReplaceAll("dilep.mass()<","").Atof();
      }
      if (plot[pl]=="dPhi*180./TMath::Pi()") {
	x1 = -10000.;
	x2 = TString(dPhi.GetTitle()).ReplaceAll("dPhi<TMath::Pi()*","").ReplaceAll("/180.","").Atof();
      }
      if (plot[pl]=="lep1.pt()") {
	x1 = TString(lep1pt.GetTitle()).ReplaceAll("lep1.pt()>","").Atof();
	x2 = +10000.;
      }
      if (plot[pl]=="lep2.pt()") {
	x1 = TString(lep2pt.GetTitle()).ReplaceAll("lep2.pt()>","").Atof();
	x2 = +10000.;
      }
      if (plot[pl]=="mt") {
	TString tit = TString(mt.GetTitle());
	TString s1 = TString(mt.GetTitle()).Remove(tit.First('&'),tit.Length()-1);
	TString s2 = TString(mt.GetTitle()).Remove(0,tit.Last('&')+1);
	x1 = s1.ReplaceAll("mt>","").ReplaceAll("&","").Atof();
	x2 = s2.ReplaceAll("mt<","").ReplaceAll("&","").Atof();
      }
    }

    TCanvas c;
    THStack* hs = new THStack("hs",plot[pl]);
    TLegend* leg = new TLegend(0.15,0.74,0.50,0.94);
    leg->SetTextFont(42);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(2);
    leg->SetLineWidth(0);
    leg->SetLineColor(kWhite);
    leg->SetShadowColor(kWhite);

    TH1F* plotmc = 0;
    TH1F* plotSignal = 0;
    TH1F* plotData = 0;
    TH1F* h_band = 0;
    TH1F* plotSyst = new TH1F("plotSyst","plotSyst",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
    TH1F h2("h2","h2",2,0,2);
    TH1F h3("h3","h3",2,0,2);
    TH1F h4("h4","h4",2,0,2);
    TH1F h5("h5","h5",2,0,2);
    TH1F h6("h6","h6",2,0,2);
    TH1F h7("h7","h7",2,0,2);
    TLine l1,l2;
    if (!compareData) {
      //standard plots
      for (int i=0;i<nMC;++i){    
	//cout << mcs[i] << " " << sfs[i] << endl;
	TFile *_mc = 0;
	if (dyFromLowMet && mcs[i]=="dyll")  _mc = TFile::Open(dir+"/../skim_dy/"+mcs[i]+".root");
	else _mc = TFile::Open(dir+"/"+mcs[i]+".root");
	TTree * mc = (TTree*) _mc->Get("tree");
	float correction = sfs[i];
	if (lumi>0.) {
	  if (mcs[i]=="dyll" && (mass>1 || fs=="of")) { 
	    mc->SetWeight(lumi);
	  } else mc->SetWeight(lumi*correction);
	}
	plotmc = new TH1F("plotmc_"+mcs[i],"plotmc",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	if (mcs[i].Contains("wjets") && wjetsFromData) {
	  TH1F* hspill = new TH1F("plotmc_spill","plotmc_spill",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	  spill->SetWeight(lumi);
	  spill->Draw(plot[pl]+">>plotmc_spill",Form("scale1fb*sfWeightTrig*sfWeightPU*sfWeightFR*( %s && %s && sfWeightFR>-99. && ( !((abs(lep1McId)==11||abs(lep1McId)==13)&&(abs(lep2McId)==11||abs(lep2McId)==13) )||dstype==51 ) )",cut.GetTitle(),lplusf.GetTitle()),"g");
	  fakes->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("sfWeightFR*(%s && %s && %s && sfWeightFR>-99.)",cut.GetTitle(),lplusf.GetTitle(),runrange.GetTitle()),"g");
	  plotmc->Add(hspill);
	} else if (mcs[i].Contains("dyll") && dyFromLowMet && plot[pl]!="mt" && plot[pl]!="pmet" && plot[pl]!="pTrackMet" && plot[pl]!="dymva") {
	  TString lm_cut = cut.GetTitle();
	  if (njets<2) {
	    lm_cut = lm_cut.ReplaceAll("dymva>","dymva>-0.9&&dymva<");
	  } else {
	    lm_cut = lm_cut.ReplaceAll("met>45.","pmet>20.&&met<45.");
	  }
	  mc->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && %s)",lm_cut.Data(),leps.GetTitle()),"g");
	} else {
	  mc->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && %s)",cut.GetTitle(),leps.GetTitle()),"g");	
	}
	float integral = plotmc->Integral(0,plotmc->GetNbinsX()+1);//integral includes under/over flow bins
	if ( mcs[i]=="dyll" && dyFromLowMet && fs!="of" ) {
	  float integral_cut = plotmc->Integral(plotmc->FindBin(x1),plotmc->FindBin(x2));//integral in signal region
	  if (integral_cut>0) plotmc->Scale(sfs[i]/integral_cut);
	  integral = plotmc->Integral(0,plotmc->GetNbinsX()+1);
	}
	//cout << mcs[i] << " " << integral << endl;
	plotmc->SetLineColor(colors[i]);
	plotmc->SetFillColor(colors[i]);
	fillOverflowInLastBin(plotmc);
	hs->Add(plotmc);
	float systerr = ks[i]-1;
	addInQuadrature(plotSyst,plotmc,systerr);
	if (allleg) leg->AddEntry(plotmc,TString(mcs[i]+": "+Form("%.2f",integral)),"f");
      }
      if (logy) {
	c.SetLogy();
	hs->SetMinimum(0.1);
      }
      hs->SetTitle();
      hs->Draw();
      TString unit = "";
      if (xtitle[pl].Contains("[")) {
	unit = xtitle[pl];
	unit.Remove(0,unit.First("[")+1);
	unit.Remove(unit.First("]"),unit.First("]")+1);
      }
      hs->GetYaxis()->SetTitle(Form("events/%.1f %s",hs->GetXaxis()->GetBinWidth(1),unit.Data()));
      hs->GetXaxis()->SetTitle(xtitle[pl]);
      hs->GetXaxis()->SetNdivisions(505);
      hs->SetMaximum(hs->GetMaximum()*1.4);
      if (logy) {
	hs->SetMaximum(hs->GetMaximum()*1.4);
      }

      hs->Draw();

      if (doSignal) {
	//if we want to plot signal
	plotSignal = new TH1F("plotSignal","plotSignal",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	signal->Draw(plot[pl]+">>plotSignal",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && %s)",lumi,cut.GetTitle(),leps.GetTitle()),"same");
	//cout << plotSignal->Integral(0,plotSignal->GetNbinsX()+1) << endl;
	plotSignal->SetLineWidth(2);
	plotSignal->SetLineStyle(1);
	plotSignal->SetLineColor(kRed);
	fillOverflowInLastBin(plotSignal);
	leg->AddEntry(plotSignal," m_{H}="+mh+" GeV","l");
      }

      if (dodata) {
	plotData = new TH1F("plotData","plotData",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	data->Draw(plot[pl]+">>plotData",cut&&leps&&runrange,"EP,same");
	c.Update();  
	//cout << "data" << endl;
	//float integral = plotData->Integral(0,plotmc->GetNbinsX()+1);//integral includes under/over flow bins
	//cout << integral << endl;
	plotData->SetMarkerStyle(20);
	fillOverflowInLastBin(plotData);
	int maxb=plotData->GetMaximumBin();
	double max=plotData->GetBinContent(maxb)+plotData->GetBinError(maxb);
	if (doSignal) {
	  max = TMath::Max(max,plotSignal->GetBinContent(plotSignal->GetMaximumBin()));
	}
	if (max>(hs->GetMaximum())) {
	  hs->SetMaximum(max*1.4);
	}

	if (doRatio) {
	  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
	  pad2->SetBottomMargin(0.4);
	  pad2->SetTopMargin(0.);
	  pad2->Draw();
	  pad2->cd();
	
	  TH1F* ratio = new TH1F("ratio","",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	  for (unsigned int ibin=1;ibin<=nbins.Atoi();ibin++) {
	    float den = ((TH1*)(hs->GetStack()->Last()))->GetBinContent(ibin);
	    float num = plotData->GetBinContent(ibin);
	    if (den>0 && num>0) {
	      if (plotData->GetBinError(ibin)/den<0.75*num/den || fabs(num/den-1.)<1.) {
		ratio->SetBinContent(ibin,num/den);
		ratio->SetBinError(ibin,plotData->GetBinError(ibin)/den);
	      }
	    } 
	  }
	  ratio->GetXaxis()->SetTitle(xtitle[pl]);
	  ratio->GetXaxis()->SetTitleSize(0.15);
	  ratio->GetXaxis()->SetLabelSize(0.15);
	  ratio->GetYaxis()->SetTitle("Data/MC");
	  ratio->GetYaxis()->SetTitleOffset(0.5);
	  ratio->GetYaxis()->SetTitleSize(0.12);
	  ratio->GetYaxis()->SetLabelSize(0.1);
	  ratio->GetYaxis()->SetNdivisions(505);
	  ratio->GetYaxis()->SetRangeUser(0.,2.);
	  pad2->SetGridy();
	  ratio->SetMarkerStyle(21);
	  ratio->Draw("EP");
	  /*
	  TString plotcopy = plot[pl];
	  plotcopy.ReplaceAll(".","");
	  plotcopy.ReplaceAll("(","");
	  plotcopy.ReplaceAll(")","");
	  ratio->SaveAs("dirplots/mh"+mh+"_nj"+nj+"/"+plotcopy+"_"+cut.GetName()+".C");
	  */
	  if (doSyst) {
	    h_band = new TH1F("hr_band","hr_band",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	    for (int mbin=1;mbin<=h_band->GetNbinsX();++mbin) {  
	      float den = ((TH1*)(hs->GetStack()->Last()))->GetBinContent(mbin);
	      if (den>0) {
		h_band->SetBinContent(mbin,1);
		h_band->SetBinError(mbin,sqrt(plotSyst->GetBinContent(mbin))/den);
	      }
	    }
	    h_band->SetMarkerSize(0);
	    h_band->SetFillColor(kBlack);
	    h_band->SetLineColor(kBlack);
	    h_band->SetFillStyle(3005);//
	    h_band->Draw("E2 same");
	  }
	  TLine line(minbin.Atof(),1,maxbin.Atof(),1);
	  line.SetLineColor(kRed);
	  line.Draw("same");
	  c.cd();
	
	  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
	  pad1->SetBottomMargin(0.);
	  pad1->Draw();
	  pad1->cd();
	  hs->GetXaxis()->SetTitle("");
	  hs->GetXaxis()->SetLabelSize(0.);
	} else {
	  c.SetBottomMargin(0.15);
	  c.SetTopMargin(0.05);
	}
	hs->Draw("HIST");
	
	c.Update();  
	if (doSyst) {
	  h_band = new TH1F("h_band","h_band",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	  for (int mbin=1;mbin<=h_band->GetNbinsX();++mbin) {  
	    h_band->SetBinContent(mbin,((TH1*)(hs->GetStack()->Last()))->GetBinContent(mbin));
	    h_band->SetBinError(mbin,sqrt(plotSyst->GetBinContent(mbin)));
	  }
	  h_band->SetMarkerSize(0);
	  h_band->SetFillColor(kBlack);
	  h_band->SetLineColor(kBlack);
	  h_band->SetFillStyle(3005);//
	  h_band->Draw("E2 same");
	}
	if (dodata) plotData->Draw("EP,same");
	if (doSignal) plotSignal->Draw("same");

	l1 = TLine(x1,0,x1,max*1.1);
	l2 = TLine(x2,0,x2,max*1.1);
	l1.SetLineWidth(2.);
	l2.SetLineWidth(2.);
	l1.SetLineStyle(2);
	l2.SetLineStyle(2);
	//l1.SetLineColor(kMagenta);
	//l2.SetLineColor(kMagenta);
	l1.Draw();//"same"
	l2.Draw();//"same"

	c.cd();
	
	leg->AddEntry(plotData," data","pl");
	
      }
      if (!allleg) {

	if (getIndexForProcess(mcs,nMC,"qqww")>=0 || getIndexForProcess(mcs,nMC,"ggww")>=0) {
	  h2.SetFillColor(colors[max(getIndexForProcess(mcs,nMC,"qqww"),getIndexForProcess(mcs,nMC,"ggww"))]); 
	  leg->AddEntry(&h2," WW","f");
	}
	if (getIndexForProcess(mcs,nMC,"dyll")>=0)  {
	  h3.SetFillColor(colors[getIndexForProcess(mcs,nMC,"dyll")]); 
	  leg->AddEntry(&h3," Z/#gamma*","f");
	}
	if (getIndexForProcess(mcs,nMC,"ttbar_powheg")>=0 || getIndexForProcess(mcs,nMC,"tw")>=0) {
	  h4.SetFillColor(colors[max(getIndexForProcess(mcs,nMC,"ttbar_powheg"),getIndexForProcess(mcs,nMC,"tw"))]);  
	  leg->AddEntry(&h4," Top","f");
	}
	if (getIndexForProcess(mcs,nMC,"wz")>=0 || getIndexForProcess(mcs,nMC,"zz") || getIndexForProcess(mcs,nMC,"www")>=0) {
	  h5.SetFillColor(colors[max(getIndexForProcess(mcs,nMC,"wz"),max(getIndexForProcess(mcs,nMC,"zz"),getIndexForProcess(mcs,nMC,"www")))]); 
	  leg->AddEntry(&h5," VZ","f");
	}
	if (getIndexForProcess(mcs,nMC,"wjets")>=0 || getIndexForProcess(mcs,nMC,"wgamma")>=0 || getIndexForProcess(mcs,nMC,"zgamma")>=0 || getIndexForProcess(mcs,nMC,"wglll")>=0) {
	  h6.SetFillColor(colors[max(getIndexForProcess(mcs,nMC,"wjets"),max(getIndexForProcess(mcs,nMC,"wgamma"),max(getIndexForProcess(mcs,nMC,"zgamma"),getIndexForProcess(mcs,nMC,"wglll"))))]);  
	  leg->AddEntry(&h6," W+jets","f");
	} 
	if (getIndexForProcess(mcs,nMC,"hww125")>=0) {
	  h7.SetFillColor(colors[getIndexForProcess(mcs,nMC,"hww125")]);  
	  leg->AddEntry(&h7," HWW125","f");
	}
	if (doSyst) {
	  leg->AddEntry(h_band," syst.","f");
	}
      }
    } else {
      //if we want to compare two different data samples
      TH1F* plotData2011A = new TH1F("plotData2011A","",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
      data->Draw(plot[pl]+">>plotData2011A",leps&&cut+"run<=173692","EP");
      plotData2011A->SetLineWidth(2);
      TH1F* plotData2011B = new TH1F("plotData2011B","",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
      data->Draw(plot[pl]+">>plotData2011B",leps&&cut+"run>173692","EP");
      plotData2011B->SetMarkerStyle(21);
      plotData2011B->SetMarkerColor(kRed);
      plotData2011A->Sumw2();
      plotData2011B->Sumw2();
      plotData2011A->Scale(1.);
      plotData2011B->Scale(1.*1.2);
      int maxb2011A=plotData2011A->GetMaximumBin();
      double max2011A=plotData2011A->GetBinContent(maxb2011A);
      int maxb2011B=plotData2011B->GetMaximumBin();
      double max2011B=plotData2011B->GetBinContent(maxb2011B);
      //cout << max2011A << " " << max2011B << endl;
      plotData2011A->GetYaxis()->SetRangeUser(0,TMath::Max(max2011A,max2011B)*1.3);
      plotData2011A->GetYaxis()->SetTitle("events scaled to 2011A");
      plotData2011A->GetXaxis()->SetTitle(xtitle[pl]);
      for (int b = 1;b<=plotData2011A->GetXaxis()->GetNbins();++b) {
	plotData2011B->SetBinError(b,sqrt( pow(plotData2011B->GetBinError(b),2) + pow(plotData2011A->GetBinError(b),2) ));
	//cout << plotData2011A->GetBinContent(b) << " " << plotData2011B->GetBinContent(b) << " " << plotData2011B->GetBinError(b) << endl;
	plotData2011A->SetBinError(b,0.);
      }
      plotData2011A->Draw("");
      plotData2011B->Draw("same,EP");
      leg->AddEntry(plotData2011A,"2011A data","l");
      leg->AddEntry(plotData2011B,"2011B data","p");
    }
    
    leg->Draw();
    if (!compareData) {
      TPaveText* labelcms  = new TPaveText(0.60,0.80,0.94,0.94,"NDCBR");
      labelcms->SetTextAlign(12);
      labelcms->SetTextSize(0.035);
      labelcms->SetFillColor(kWhite);
      labelcms->AddText("CMS Preliminary");
      labelcms->AddText("#sqrt{s} = 8 TeV");
      labelcms->AddText(Form("L = %.2f fb^{-1}",lumi));
      labelcms->SetBorderSize(0);
      labelcms->SetTextFont(42);
      labelcms->SetLineWidth(2);
      labelcms->Draw();
    }

    //add label for FS
    TPaveText* labelfs  = new TPaveText(0.60,0.75,0.94,0.80,"NDCBR");
    labelfs->SetTextAlign(12);
    labelfs->SetTextSize(0.035);
    labelfs->SetFillColor(kWhite);
    labelfs->SetBorderSize(0);
    labelfs->SetTextFont(42);
    labelfs->SetLineWidth(2);
    if (njets>=0) {
      if (fs=="of") {
	labelfs->AddText(Form("%i-j, #mue+e#mu final state",njets));
      } else if (fs=="sf") {
	labelfs->AddText(Form("%i-j, #mu#mu+ee final state",njets));
      } else if (fs=="mm") {
	labelfs->AddText(Form("%i-j, #mu#mu final state",njets));
      } else if (fs=="me") {
	labelfs->AddText(Form("%i-j, #mue final state",njets));
      } else if (fs=="em") {
	labelfs->AddText(Form("%i-j, e#mu final state",njets));
      } else if (fs=="ee") {
	labelfs->AddText(Form("%i-j, ee final state",njets));
      } else {
	labelfs->AddText(Form("%i-j, all final states",njets));
      } 
    } else {
      if (fs=="of") {
	labelfs->AddText("#mue+e#mu final state");
      } else if (fs=="sf") {
	labelfs->AddText("#mu#mu+ee final state");
      } else if (fs=="mm") {
	labelfs->AddText("#mu#mu final state");
      } else if (fs=="me") {
	labelfs->AddText("#mue final state");
      } else if (fs=="em") {
	labelfs->AddText("e#mu final state");
      } else if (fs=="ee") {
	labelfs->AddText("ee final state");
      } else {
	labelfs->AddText("all final states");
      } 
    }
    labelfs->Draw();

    c.Update();  
    if (!doRatio) {
      c.GetFrame()->DrawClone();
      c.RedrawAxis();
    }

    TString cutName = cut.GetName();
    plot[pl].ReplaceAll(".","");
    plot[pl].ReplaceAll("(","");
    plot[pl].ReplaceAll(")","");
    plot[pl].ReplaceAll("+","");
    plot[pl].ReplaceAll("-","");
    plot[pl].ReplaceAll(",","");
    plot[pl].ReplaceAll("*180/TMath::Pi","");
    gSystem->Exec("mkdir -p dirplots/mh"+mh+"_nj"+nj);
    int nExt = sizeof(extensions)/sizeof(TString);
    for (int ext=0;ext<nExt;++ext) {
      TString extension = extensions[ext];
      if (logy) extension="_log"+extension;
      if (topTagged) extension="_tt"+extension;
      if (sameSign) extension="_ss"+extension;
      if (cutName.Length()==0) c.SaveAs("dirplots/mh"+mh+"_nj"+nj+"/"+plot[pl]+extension);
      else c.SaveAs("dirplots/mh"+mh+"_nj"+nj+"/"+plot[pl]+"_"+cutName+extension);
    }
    if (!compareData) {
      delete plotmc;
      if (plotData) delete plotData; 
      delete hs;
    } else {
      //delete plotData2011A;//fixme
      //delete plotData2011B;
    }
    delete leg;
  }

}

