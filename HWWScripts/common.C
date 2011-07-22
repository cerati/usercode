#ifndef COMMON_C
#define COMMON_C

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h> 
#include <string.h>
#include <vector>
#include <utility>
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
#include "TRandom3.h"
#include "Smurf/Core/SmurfTree.h"
#include <vector>

//copy here to avoid SmurfTree::
enum Selection {
  BaseLine          = 1UL<<0,  // pt(reco)>20/10, acceptance,!STA muon, mll>12
  ChargeMatch       = 1UL<<1,  // q1*q2<0
  Lep1FullSelection = 1UL<<2,  // full id, isolation, d0, dz etc
  Lep1LooseEleV1    = 1UL<<3,  // electron fakeable object selection is passed V1
  Lep1LooseEleV2    = 1UL<<4,  // electron fakeable object selection is passed V2
  Lep1LooseEleV3    = 1UL<<5,  // electron fakeable object selection is passed V3
  Lep1LooseEleV4    = 1UL<<6,  // electron fakeable object selection is passed V4
  Lep1LooseMuV1     = 1UL<<7,  // muon fakeable object selection (relIso<1.0)
  Lep1LooseMuV2     = 1UL<<8,  // muon fakeable object selection (relIso<0.4)
  Lep2FullSelection = 1UL<<9,  // full id, isolation, d0, dz etc
  Lep2LooseEleV1    = 1UL<<10, // electron fakeable object selection is passed V1
  Lep2LooseEleV2    = 1UL<<11, // electron fakeable object selection is passed V2
  Lep2LooseEleV3    = 1UL<<12, // electron fakeable object selection is passed V3
  Lep2LooseEleV4    = 1UL<<13, // electron fakeable object selection is passed V4
  Lep2LooseMuV1     = 1UL<<14, // muon fakeable object selection (relIso<1.0)
  Lep2LooseMuV2     = 1UL<<15, // muon fakeable object selection (relIso<0.4)
  FullMET           = 1UL<<16, // full met selection
  ZVeto             = 1UL<<17, // event is not in the Z-mass peak for ee/mm final states
  TopTag            = 1UL<<18, // soft muon and b-jet tagging for the whole event regardless of n-jets (non-zero means tagged)
  TopVeto           = 1UL<<19, // soft muon and b-jet tagging for the whole event regardless of n-jets (zero means tagged)
  OneBJet           = 1UL<<20, // 1-jet events, where the jet is b-tagged (top control sample with one b-quark missing)
  TopTagNotInJets   = 1UL<<21, // soft muon and b-jet tagging for areas outside primary jets (non-zero means tagged)
  ExtraLeptonVeto   = 1UL<<22  // extra lepton veto, DR(muon-electron)>=0.3
};
unsigned int wwSelection     = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelNoELV      = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ZVeto|TopVeto;
unsigned int wwSelectionELV  = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ExtraLeptonVeto;
unsigned int wwSelectionNoZV = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ExtraLeptonVeto|TopVeto;
unsigned int wwSelNoZVNoMet  = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ExtraLeptonVeto|TopVeto;
unsigned int wwSelectionNoTV = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ZVeto|ExtraLeptonVeto;
unsigned int wwSelNoMet      = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelectionFOe1 = BaseLine|ChargeMatch|Lep1LooseEleV4|Lep2FullSelection|FullMET|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelectionFOe2 = BaseLine|ChargeMatch|Lep1FullSelection|Lep2LooseEleV4|FullMET|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelectionFOm1 = BaseLine|ChargeMatch|Lep1LooseMuV2|Lep2FullSelection|FullMET|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelectionFOm2 = BaseLine|ChargeMatch|Lep1FullSelection|Lep2LooseMuV2|FullMET|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelectionNoLep= BaseLine|ChargeMatch|FullMET|ZVeto|ExtraLeptonVeto|TopVeto;
unsigned int wwSelNoLepNoTV  = BaseLine|ChargeMatch|FullMET|ZVeto|ExtraLeptonVeto;
unsigned int wwSelLepOnly    = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ExtraLeptonVeto;
unsigned int noVeto          = 1UL<<30;
unsigned int noCut           = 1UL<<0;

TString dir_mc         = "/smurf/data/EPS/tas/";
//TString dir_mc         = "/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/";
TString dir_mc_mit     = "/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/";
TString data_file      = "/smurf/data/EPS/tas/data-met20-1092ipb";
TString fr_file_mit    = "/smurf/data/EPS/auxiliar/FakeRates_SmurfV6.root";
TString fr_file_el_tas = "/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/ww_el_fr.root";
TString fr_file_mu_tas = "/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/ww_mu_fr.root";
TString eff_file       = "/smurf/data/EPS/auxiliar/efficiency_results_v6.root";

bool passJson(int run, int lumi) {
  ifstream ifs( "goodruns.txt" );
  string temp;
  istringstream instream;
  TString cutstr = "";
  while( getline( ifs, temp ) ) {
    instream.clear();
    instream.str(temp);
    int runJson;
    int lumiJsonMin;
    int lumiJsonMax;
    instream>>runJson;
    instream>>lumiJsonMin;
    instream>>lumiJsonMax;
    if (runJson==run && lumi>=lumiJsonMin && lumi<=lumiJsonMax) return true;
  }
  return false;
}

float getScale1fb(TString sample) {
  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();
  dataEvent->tree_->GetEntry(0);
  float result = dataEvent->scale1fb_;
  delete dataEvent;
  return result;
}

void getCutValues(int mass, float& lep1pt,float& lep2pt,float& dPhi,float& mll,float& mtL,float& mtH,float& himass, bool doMVA=false){

  if (mass==0) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 180.;
    mll    = 9999.;
    mtL    = 0.;
    mtH    = 9999.;
    himass = 0.;
  } else if (mass==115) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 115.;
    mll    = 40.;
    mtL    = 70.;
    mtH    = 110.;
    himass = 100.;
  } else if (mass==120) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 115.;
    mll    = 40.;
    mtL    = 70.;
    mtH    = 120.;
    himass = 100.;
  } else if (mass==130) {
    lep1pt = 25.;
    lep2pt = 10.;
    dPhi   = 90.;
    mll    = 45.;
    mtL    = 75.;
    mtH    = 125.;
    himass = 100.;
  } else if (mass==140) {
    lep1pt = 25.;
    lep2pt = 15.;
    dPhi   = 90.;
    mll    = 45.;
    mtL    = 80.;
    mtH    = 130.;
    himass = 100.;
  } else if (mass==150) {
    lep1pt = 27.;
    lep2pt = 25.;
    dPhi   = 90.;
    mll    = 50.;
    mtL    = 80.;
    mtH    = 150.;
    himass = 100.;
  } else if (mass==160) {
    lep1pt = 30.;
    lep2pt = 25.;
    dPhi   = 60.;
    mll    = 50.;
    mtL    = 90.;
    mtH    = 160.;
    himass = 100.;
  } else if (mass==170) {
    lep1pt = 34.;
    lep2pt = 25.;
    dPhi   = 60.;
    mll    = 50.;
    mtL    = 110.;
    mtH    = 170.;
    himass = 100.;
  } else if (mass==180) {
    lep1pt = 36.;
    lep2pt = 25.;
    dPhi   = 70.;
    mll    = 60.;
    mtL    = 120.;
    mtH    = 180.;
    himass = 100.;
  } else if (mass==190) {
    lep1pt = 38.;
    lep2pt = 25.;
    dPhi   = 90.;
    mll    = 80.;
    mtL    = 120.;
    mtH    = 190.;
    himass = 100.;
  } else if (mass==200) {
    lep1pt = 40.;
    lep2pt = 25.;
    dPhi   = 100.;
    mll    = 90.;
    mtL    = 120.;
    mtH    = 200.;
    himass = 100.;
  } else if (mass==250) {
    lep1pt = 55.;
    lep2pt = 25.;
    dPhi   = 140.;
    mll    = 150.;
    mtL    = 120.;
    mtH    = 250.;
    himass = 100.;
  } else if (mass==300) {
    lep1pt = 70.;
    lep2pt = 25.;
    dPhi   = 175.;
    mll    = 200.;
    mtL    = 120.;
    mtH    = 300.;
    himass = 100.;
  } else {
    cout << "MASS POINT NOT SUPPORTED!!!!!! mH=" << mass << endl;
  }

  if (doMVA){
    if (mass==115) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 70.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==120) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 70.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==130) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 80.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==140) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 90.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==150) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 100.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==160) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 100.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==170) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 100.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==180) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 110.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==190) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 120.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else if (mass==200) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mll    = 130.;
      mtL    = 0.;
      mtH    = 999.;
      himass = 100.;
    } else {
      cout << "MASS POINT NOT SUPPORTED!!!!!! mH=" << mass << endl;
    }
  }

}

void getCutMasks(unsigned int njets, unsigned int& baseline_toptag, unsigned int& control_top, unsigned int& control_toptag, unsigned int& veto, unsigned int& nj_top){
  //0-jet bin
  unsigned int j1bt   = OneBJet;
  unsigned int ttag0j = TopTagNotInJets;
  unsigned int veto0j = noVeto;
  //1-jet bin
  //unsigned int ttaghj = OneBJet;//not working in dimas ntuple (requires  jets==1)
  unsigned int veto1j = ttag0j;
  //cook it
  unsigned int top_sel = 0;
  unsigned int ttag = 0;
  nj_top = 0;
  veto = 0;
  if (njets==0) {
    nj_top = 1;
    top_sel = j1bt;
    ttag = ttag0j;
    veto = veto0j;
  }
  else if (njets==1) {
    nj_top = 2;
    top_sel = noCut;
    ttag = noCut;//ttaghj
    veto = veto1j;
  }  
  baseline_toptag = wwSelectionNoTV|ttag;
  control_top     = wwSelectionNoTV|top_sel;
  control_toptag  = control_top|ttag;
}

pair<float, float> getYield(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString region, float lumi, bool useJson=0, 
			    bool applyEff=true, bool doFake=false) {

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);

  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");

  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;

  //efficiency correction stuff
  TFile *effs=0;
  TH2F *effElSel=0,*effMuSel=0,*effElDbl=0,*effElSgl=0,*effMuDbl=0,*effMuSgl=0;
  if (isMC&&applyEff) {
    effs = TFile::Open(eff_file);
    effElSel = (TH2F*) effs->Get("h2_results_electron_selection");
    effMuSel = (TH2F*) effs->Get("h2_results_muon_selection");
    effElDbl = (TH2F*) effs->Get("h2_results_electron_double");
    effElSgl = (TH2F*) effs->Get("h2_results_electron_single");
    effMuDbl = (TH2F*) effs->Get("h2_results_muon_double");
    effMuSgl = (TH2F*) effs->Get("h2_results_muon_single");
  }

  //fake stuff
  bool useMit = true;
  TFile *fE=0, *fM=0;
  TH2F *eFR=0,*mFR=0;
  TH2F *elPFY=0,*muPFY=0;
  //FIXME: should be gotten from files
  float maxPt=34.9, maxEta=2.4;
  if (doFake) {
    if (useMit) {
      fE = TFile::Open(fr_file_mit);
      eFR = (TH2F*) fE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta");
      fM = TFile::Open(fr_file_mit);
      mFR = (TH2F*) fM->Get("MuonFakeRate_M2_ptThreshold15_PtEta");
    } else {
      fE = TFile::Open(fr_file_el_tas);
      eFR = (TH2F*) fE->Get("el_fr_v4_35");
      fM = TFile::Open(fr_file_mu_tas);
      mFR = (TH2F*) fM->Get("mu_fr_m2_15");
    }
    elPFY = new TH2F("elPFY","elPFY",eFR->GetXaxis()->GetNbins(),eFR->GetXaxis()->GetXmin(),eFR->GetXaxis()->GetXmax(),
		                     eFR->GetYaxis()->GetNbins(),eFR->GetYaxis()->GetXmin(),eFR->GetYaxis()->GetXmax());
    muPFY = new TH2F("muPFY","muPFY",mFR->GetXaxis()->GetNbins(),mFR->GetXaxis()->GetXmin(),mFR->GetXaxis()->GetXmax(),
			             mFR->GetYaxis()->GetNbins(),mFR->GetYaxis()->GetXmin(),mFR->GetYaxis()->GetXmax());
  }

  float weight = 1.;
  float yield = 0.;
  float error = 0.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    if (useJson && !passJson(dataEvent->run_,dataEvent->lumi_)) continue;    
    if ( dataEvent->njets_!=njets) continue;
    if ( (dataEvent->cuts_ & cut) != cut ) continue;
    if ( (dataEvent->cuts_ & veto) == veto ) continue;
    //WARNING: do not define region names that are subset of others!!!!!!!
    //final states
    if ( region.Contains("mmfs") && dataEvent->type_!=0 ) continue;
    if ( region.Contains("mefs") && dataEvent->type_!=1 ) continue;
    if ( region.Contains("emfs") && dataEvent->type_!=2 ) continue;
    if ( region.Contains("eefs") && dataEvent->type_!=3 ) continue;
    //mass regions
    if ( (region.Contains("leppts")||region.Contains("sideband")||region.Contains("mtside")||region.Contains("dphiside")||region.Contains("massreg")||
	  region.Contains("mtreg")||region.Contains("dphireg")) && (dataEvent->lep1_.pt()<lep1pt || dataEvent->lep2_.pt()<lep2pt) ) continue;
    if ( (region.Contains("sideband")||region.Contains("mtside")||region.Contains("dphiside")) && (dataEvent->dilep_.mass()<himass) ) continue;
    if ( (region.Contains("massreg")||region.Contains("mtreg")||region.Contains("dphireg")) && (dataEvent->dilep_.mass()>mll) ) continue;
    if ( (region.Contains("mtside")||region.Contains("dphiside")||region.Contains("mtreg")||region.Contains("mtcut")||region.Contains("dphireg")) 
	 && (dataEvent->mt_<mtL || dataEvent->mt_>mtH) ) continue;
    if ( (region.Contains("dphiside")||region.Contains("dphireg")||region.Contains("dphicut")) && (dataEvent->dPhi_>dPhi*TMath::Pi()/180.) ) continue;
    if ( (region.Contains("zregion")) && fabs(dataEvent->dilep_.mass()-91.1876)>15. ) continue;
    if ( (region.Contains("zvetoall")) && fabs(dataEvent->dilep_.mass()-91.1876)<15. ) continue;
    //additional higgs cuts
    if ( region.Contains("dphijet")   && ( dataEvent->type_!=1 && dataEvent->type_!=2 && (dataEvent->jet1_.pt()>15&&dataEvent->dPhiDiLepJet1_>165.*TMath::Pi()/180.)) ) continue;
    if ( region.Contains("minmet40")  && ( dataEvent->type_!=1 && dataEvent->type_!=2 && min(dataEvent->pmet_,dataEvent->pTrackMet_)<40.0) ) continue;
    if ( region.Contains("dpjallfs")  && dataEvent->jet1_.pt()>15 && dataEvent->dPhiDiLepJet1_>165.*TMath::Pi()/180. ) continue;
    if ( region.Contains("mm40allfs") && min(dataEvent->pmet_,dataEvent->pTrackMet_)<40.0 ) continue;
    //cuts for Rout/in
    if ( region.Contains("met2022") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<20.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>22.0) ) continue;
    if ( region.Contains("met2226") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<22.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>26.0) ) continue;
    if ( region.Contains("met2631") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<26.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>31.0) ) continue;
    if ( region.Contains("met3140") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<31.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>40.0) ) continue;
    if ( region.Contains("met40up") &&  min(dataEvent->pmet_,dataEvent->pTrackMet_)<40.0 ) continue;
    //cuts on btag
    if ( region.Contains("btag1or2")  && dataEvent->jet1Btag_<2.1 && dataEvent->jet2Btag_<2.1 ) continue;
    if ( region.Contains("btag1and2") && (dataEvent->jet1Btag_<2.1 || dataEvent->jet2Btag_<2.1) ) continue;
    if ( region.Contains("btagJet1")  && dataEvent->jet1Btag_<2.1 ) continue;
    if ( region.Contains("btagJet2")  && dataEvent->jet2Btag_<2.1 ) continue;
    if ( region.Contains("nobJet1")   && dataEvent->jet1Btag_>2.1 ) continue;
    if ( region.Contains("nobJet2")   && dataEvent->jet2Btag_>2.1 ) continue;

    if (!doFake) {
      float effSF=1.;
      if (isMC&&applyEff) {
	float selLep1=1.,selLep2=1.,dblLep1=1.,dblLep2=1.,sglLep1=1.,sglLep2=1.;
	float pt1 = min(dataEvent->lep1_.pt(),49.);
	float eta1 = min(fabs(dataEvent->lep1_.eta()),2.4);
	float pt2 = min(dataEvent->lep2_.pt(),49.);
	float eta2 = min(fabs(dataEvent->lep2_.eta()),2.4);
	if (abs(dataEvent->lid1_)==11) {
	  selLep1 = effElSel->GetBinContent(effElSel->FindBin(pt1,eta1));
	  dblLep1 = effElDbl->GetBinContent(effElDbl->FindBin(pt1,eta1));
	  sglLep1 = effElSgl->GetBinContent(effElSgl->FindBin(pt1,eta1));
	} else {
	  selLep1 = effMuSel->GetBinContent(effMuSel->FindBin(pt1,eta1));
	  dblLep1 = effMuDbl->GetBinContent(effMuDbl->FindBin(pt1,eta1));
	  sglLep1 = effMuSgl->GetBinContent(effMuSgl->FindBin(pt1,eta1));      
	}
	if (abs(dataEvent->lid2_)==11) {
	  selLep2 = effElSel->GetBinContent(effElSel->FindBin(pt2,eta2));
	  dblLep2 = effElDbl->GetBinContent(effElDbl->FindBin(pt2,eta2));
	  sglLep2 = effElSgl->GetBinContent(effElSgl->FindBin(pt2,eta2));
	} else {
	  selLep2 = effMuSel->GetBinContent(effMuSel->FindBin(pt2,eta2));
	  dblLep2 = effMuDbl->GetBinContent(effMuDbl->FindBin(pt2,eta2));
	  sglLep2 = effMuSgl->GetBinContent(effMuSgl->FindBin(pt2,eta2));      
	}
	float effSelSF = selLep1*selLep2;
	float effTrgSF = dblLep1*dblLep2+sglLep1*(1-dblLep1)+sglLep2*(1-dblLep2);
	effSF = effSelSF*effTrgSF;
      }
      yield = yield + weight*effSF;
      //fixme: should consider the error of effSF
      error = error + pow(weight*effSF,2);
    } 
    else {//DO FAKES!!!
      if ( (dataEvent->cuts_ & Lep1LooseEleV4)    == Lep1LooseEleV4    &&
	   (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection && 
	   (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) {
	float pt = dataEvent->lep1_.pt();
	float eta = fabs(dataEvent->lep1_.eta());
	if (pt>maxPt) pt=maxPt;
	if (eta>maxEta) eta=maxEta;
	useMit ? elPFY->Fill(pt,eta,weight) : elPFY->Fill(eta,pt,weight);
      }
      if ( (dataEvent->cuts_ & Lep2LooseEleV4)    == Lep2LooseEleV4    &&
	   (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection &&
	   (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) {
	float pt = dataEvent->lep2_.pt();
	float eta = fabs(dataEvent->lep2_.eta());
	if (pt>maxPt) pt=maxPt;
	if (eta>maxEta) eta=maxEta;
	useMit ? elPFY->Fill(pt,eta,weight) : elPFY->Fill(eta,pt,weight);
      }
      if ( (dataEvent->cuts_ & Lep1LooseMuV2)     == Lep1LooseMuV2     &&
	   (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection &&
	   (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) {
	float pt = dataEvent->lep1_.pt();
	float eta = fabs(dataEvent->lep1_.eta());
	if (pt>maxPt) pt=maxPt;
	if (eta>maxEta) eta=maxEta;
	useMit ? muPFY->Fill(pt,eta,weight) : muPFY->Fill(eta,pt,weight);
      }
      if ( (dataEvent->cuts_ & Lep2LooseMuV2)     == Lep2LooseMuV2     &&
	   (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection &&
	   (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) {
	float pt = dataEvent->lep2_.pt();
	float eta = fabs(dataEvent->lep2_.eta());
	if (pt>maxPt) pt=maxPt;
	if (eta>maxEta) eta=maxEta;
	useMit ? muPFY->Fill(pt,eta,weight) : muPFY->Fill(eta,pt,weight);
      }
    }
  }
  if (doFake) {
    float elPFYield = 0;
    float elPFError = 0;
    for (int ix=1;ix<elPFY->GetNbinsX()+1;++ix) {
      for (int iy=1;iy<elPFY->GetNbinsY()+1;++iy) {
	float fr = eFR->GetBinContent(ix,iy);
	float fe = eFR->GetBinError(ix,iy);
	float nb = elPFY->GetBinContent(ix,iy);
	elPFYield += nb*fr/(1-fr);
	elPFError += pow(fr/(1-fr),2)*nb+pow(fe,2)*pow(nb,2)/pow(1-fr,4);
      }
    }
    float muPFYield = 0;
    float muPFError = 0;
    for (int ix=1;ix<muPFY->GetNbinsX()+1;++ix) {
      for (int iy=1;iy<muPFY->GetNbinsY()+1;++iy) {
	float fr = mFR->GetBinContent(ix,iy);
	float fe = mFR->GetBinError(ix,iy);
	float nb = muPFY->GetBinContent(ix,iy);
	muPFYield += nb*fr/(1-fr);
	muPFError += pow(fr/(1-fr),2)*nb+pow(fe,2)*pow(nb,2)/pow(1-fr,4);
      }
    }
    yield = elPFYield+muPFYield;
    error = sqrt(elPFError+muPFError);
  } else error = sqrt(error);
  if (isMC&&applyEff) effs->Close();
  delete dataEvent;
  //cout << yield << " " << error << endl;
  return make_pair<float, float>(yield,error);
}


#endif
