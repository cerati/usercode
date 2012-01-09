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
#include "Smurf/Core/LeptonScaleLookup.cc"
#include "CMS2/NtupleMacros/Tools/goodrun.cc"
#include "Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors.h"  
#include "Smurf/Analysis/HWWlvlv/DYBkgScaleFactors.h"  
#include "Smurf/Analysis/HWWlvlv/TopBkgScaleFactors.h"	
#include "Smurf/Analysis/HWWlvlv/WWBkgScaleFactors.h"	
#include "BDTG.class.C"//add soft link to compute bdtg values

TString main_dir    = "/smurf/cerati/skims/Run2011_Summer11_SmurfV7_42X/4700ipbWeights/";
TString topww_dir   = "wwSelNoLepNoTV/";
TString dy_dir      = "wwSelNoMetNoZVminMET20/";

TString fr_file_mit    = "/smurf/data/Winter11_4700ipb/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";

TString eff_file       = "/smurf/data/Winter11_4700ipb/auxiliar/efficiency_results_v7_42x_Full2011_4700ipb.root";
TString puw_file       = "/smurf/data/Winter11_4700ipb/auxiliar/PileupReweighting.Summer11DYmm_To_Full2011.root";
TString jsonFile       = "";//"hww.Full2011.json";

//+++ deprecated
TString fr_file_el_tas = "/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/ww_el_fr.root";
TString fr_file_mu_tas = "/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/ww_mu_fr.root";
//+++

//fixme, Fall11 test
// TString main_dir = "/smurf/cerati/skims/Run2011_Fall11_SmurfV7_42X/";
// TString puw_file = "/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root";

bool redoWeights  = 0;
bool checkWeights = 0;
bool doMVA = 0;
bool doVBF = 0;

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
  ExtraLeptonVeto   = 1UL<<22, // extra lepton veto, DR(muon-electron)>=0.3
  Lep3FullSelection = 1UL<<23,  // full id, isolation, d0, dz etc
  Lep3LooseEleV1    = 1UL<<24, // electron fakeable object selection is passed V1
  Lep3LooseEleV2    = 1UL<<25, // electron fakeable object selection is passed V2
  Lep3LooseEleV3    = 1UL<<26, // electron fakeable object selection is passed V3
  Lep3LooseEleV4    = 1UL<<27, // electron fakeable object selection is passed V4
  Lep3LooseMuV1     = 1UL<<28, // muon fakeable object selection (relIso<1.0)
  Lep3LooseMuV2     = 1UL<<29, // muon fakeable object selection (relIso<0.4)
  Trigger           = 1UL<<30  // passed a set of triggers
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
unsigned int noVeto          = 1UL<<31;
unsigned int noCut           = 1UL<<0;

bool passJson(int run, int lumi) {
  //this is ok but slow: use goodrun.cc from tas
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
  //warning this method does not work if the sample is composite 
  //(i.e. more than one scale1fb value, e.g. dy=[10-20]+[20-inf])
  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();
  dataEvent->tree_->GetEntry(0);
  float result = dataEvent->scale1fb_;
  dataEvent->tree_->Delete();
  delete dataEvent;
  return result;
}

void getCutValues(int mass, float& lep1pt,float& lep2pt,float& dPhi,float& mll,float& mtL,float& mtH,float& himass/*, bool doMVA=false*/){

  if (mass==0) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 180.;
    mll    = 9999.;
    mtL    = 0.;
    mtH    = 9999.;
    himass = 100.;
  } else if (mass==110) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 115.;
    mll    = 40.;
    mtL    = 80.;
    mtH    = 110.;
    himass = 100.;
  } else if (mass==115) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 115.;
    mll    = 40.;
    mtL    = 80.;
    mtH    = 110.;
    himass = 100.;
  } else if (mass==118) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 115.;
    mll    = 40.;
    mtL    = 80.;
    mtH    = 115.;
    himass = 100.;
  } else if (mass==120) {
    lep1pt = 20.;
    lep2pt = 10.;
    dPhi   = 115.;
    mll    = 40.;
    mtL    = 80.;
    mtH    = 120.;
    himass = 100.;
  } else if (mass==122) {
    lep1pt = 21.;
    lep2pt = 10.;
    dPhi   = 110.;
    mll    = 41.;
    mtL    = 80.;
    mtH    = 121.;
    himass = 100.;
  } else if (mass==124) {
    lep1pt = 22.;
    lep2pt = 10.;
    dPhi   = 105.;
    mll    = 42.;
    mtL    = 80.;
    mtH    = 122.;
    himass = 100.;
  } else if (mass==126) {
    lep1pt = 23.;
    lep2pt = 10.;
    dPhi   = 100.;
    mll    = 43.;
    mtL    = 80.;
    mtH    = 123.;
    himass = 100.;
  } else if (mass==128) {
    lep1pt = 24.;
    lep2pt = 10.;
    dPhi   = 95.;
    mll    = 44.;
    mtL    = 80.;
    mtH    = 124.;
    himass = 100.;
  } else if (mass==130) {
    lep1pt = 25.;
    lep2pt = 10.;
    dPhi   = 90.;
    mll    = 45.;
    mtL    = 80.;
    mtH    = 125.;
    himass = 100.;
  } else if (mass==135) {
    lep1pt = 25.;
    lep2pt = 12.;
    dPhi   = 90.;
    mll    = 45.;
    mtL    = 80.;
    mtH    = 128.;
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

  if (doMVA || doVBF) {
      lep1pt = 20.;
      lep2pt = 10.;
      dPhi   = 180.;
      mtL    = 80.;
      himass = 100.;
    if (mass==0) {
      mll    = 9999.;
      mtH    = 9999.;
    } else if (mass==110||mass==115||mass==118||mass==120||mass==122||mass==124 ) {
      mll    = 70.;
      mtH    = ((float) mass);
    } else if (mass==126||mass==128||mass==130) {
      mll    = 80.;
      mtH    = ((float) mass);
    } else if (mass==135||mass==140) {
      mll    = 90.;
      mtH    = ((float) mass);
    } else if (mass==150) {
      mll    = 100.;
      mtH    = ((float) mass);
    } else if (mass==160) {
      mll    = 100.;
      mtH    = ((float) mass);
    } else if (mass==170) {
      mll    = 100.;
      mtH    = ((float) mass);
    } else if (mass==180) {
      mll    = 110.;
      mtH    = ((float) mass);
    } else if (mass==190) {
      mll    = 120.;
      mtH    = ((float) mass);
    } else if (mass==200) {
      mll    = 130.;
      mtH    = ((float) mass);
    } else if (mass==250) {
      mll    = ((float) mass);
      mtH    = ((float) mass);
    } else if (mass==300) {
      mll    = ((float) mass);
      mtH    = ((float) mass);
    } else {
      cout << "MASS POINT NOT SUPPORTED!!!!!! mH=" << mass << endl;
    }
    if (doVBF) {
      mtL    = 0.;
      mtH    = 9999.;
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
  } else if (njets==2) {
    veto = veto0j;
  }
  baseline_toptag = wwSelectionNoTV|ttag;
  control_top     = wwSelectionNoTV|top_sel;
  control_toptag  = control_top|ttag;
}

//from /UserCode/Smurf/Analysis/HZZllvv/PileupReweighting.h
float getPileupReweightFactor(int nvtx, TH1F* puweights=0) {
  if (puweights) {
    return puweights->GetBinContent(puweights->FindBin(nvtx));
  } else {
    //this is deprecated
    float weight = 1.0;
    if     (nvtx == 0) weight = 0.00000;
    else if(nvtx == 1) weight = 0.23111;
    else if(nvtx == 2) weight = 0.82486;
    else if(nvtx == 3) weight = 1.47045;
    else if(nvtx == 4) weight = 1.83450;
    else if(nvtx == 5) weight = 1.79663;
    else if(nvtx == 6) weight = 1.46496;
    else if(nvtx == 7) weight = 1.05682;
    else if(nvtx == 8) weight = 0.70823;
    else if(nvtx == 9) weight = 0.47386;
    else if(nvtx ==10) weight = 0.32382;
    else if(nvtx ==11) weight = 0.22383;
    else if(nvtx ==12) weight = 0.17413;
    else if(nvtx ==13) weight = 0.10930;
    else if(nvtx ==14) weight = 0.09563;
    else if(nvtx ==15) weight = 0.08367;
    else if(nvtx ==16) weight = 0.05418;
    else if(nvtx ==17) weight = 0.04891;
    else if(nvtx ==18) weight = 0.03515;
    else if(nvtx >=19) weight = 0.01000;  
    return weight;
  }
}

bool passEvent(SmurfTree *dataEvent, unsigned int njets, unsigned int cut, unsigned int veto, TString region, 
	       float lep1pt,float lep2pt,float dPhi,float mll,float mtL,float mtH,float himass,
	       bool isMC, bool useJson){

    if (!isMC && (dataEvent->cuts_ & Trigger) != Trigger ) return 0;
    //if (!isMC && dataEvent->run_>172802) return 0;//FIXME LP TEST
    if (!isMC && useJson && !goodrun(dataEvent->run_,dataEvent->lumi_)) return 0;    
    if (doVBF && njets==2) {
      if (dataEvent->njets_!=2 && dataEvent->njets_!=3) return 0;
      if ( !(((dataEvent->jet1_.eta()-dataEvent->lep1_.eta() > 0 && dataEvent->jet2_.eta()-dataEvent->lep1_.eta() < 0) ||
	      (dataEvent->jet2_.eta()-dataEvent->lep1_.eta() > 0 && dataEvent->jet1_.eta()-dataEvent->lep1_.eta() < 0)) &&
	     ((dataEvent->jet1_.eta()-dataEvent->lep2_.eta() > 0 && dataEvent->jet2_.eta()-dataEvent->lep2_.eta() < 0) ||
	      (dataEvent->jet2_.eta()-dataEvent->lep2_.eta() > 0 && dataEvent->jet1_.eta()-dataEvent->lep2_.eta() < 0))) ) return 0;
      if (dataEvent->njets_==3 && dataEvent->jet3_.pt()>30 && ((dataEvent->jet1_.eta()-dataEvent->jet3_.eta() > 0 && dataEvent->jet2_.eta()-dataEvent->jet3_.eta() < 0) ||
							       (dataEvent->jet2_.eta()-dataEvent->jet3_.eta() > 0 && dataEvent->jet1_.eta()-dataEvent->jet3_.eta() < 0)) ) return 0;
      if ( TMath::Abs(dataEvent->jet1_.eta()) >= 4.5 || TMath::Abs(dataEvent->jet2_.eta()) >= 4.5) return 0;
      if ( TMath::Abs(dataEvent->jet1_.eta()-dataEvent->jet2_.eta())<=3.5) return 0;
      if ((dataEvent->jet1_+dataEvent->jet2_).mass()<=450) return 0;
    } else if ( dataEvent->njets_!=njets) return 0;
    if ( (dataEvent->cuts_ & cut) != cut ) return 0;
    if ( veto!=noVeto && (dataEvent->cuts_ & veto) == veto ) return 0;
    //WARNING: do not define region names that are subset of others!!!!!!!
    //final states
    if ( region.Contains("mmfs") && dataEvent->type_!=0 ) return 0;
    if ( region.Contains("mefs") && dataEvent->type_!=1 ) return 0;
    if ( region.Contains("emfs") && dataEvent->type_!=2 ) return 0;
    if ( region.Contains("eefs") && dataEvent->type_!=3 ) return 0;
    if ( region.Contains("sffs") && (dataEvent->type_==1 || dataEvent->type_==2) ) return 0;
    if ( region.Contains("offs") && (dataEvent->type_==0 || dataEvent->type_==3) ) return 0;
    //mass regions
    if ( (region.Contains("leppts")||region.Contains("sideband")||region.Contains("mtside")||region.Contains("dphiside")||region.Contains("massreg")||
	  region.Contains("mtreg")||region.Contains("dphireg")) && (dataEvent->lep1_.pt()<lep1pt || dataEvent->lep2_.pt()<lep2pt) ) return 0;
    if ( (region.Contains("sideband")||region.Contains("mtside")||region.Contains("dphiside")) && (dataEvent->dilep_.mass()<himass) ) return 0;
    if ( (region.Contains("massreg")||region.Contains("masscut")||region.Contains("mtreg")||region.Contains("dphireg")) && (dataEvent->dilep_.mass()>mll) ) return 0;
    if ( (region.Contains("mtside")||region.Contains("dphiside")||region.Contains("mtreg")||region.Contains("mtcut")||region.Contains("dphireg")) 
	 && (dataEvent->mt_<mtL || dataEvent->mt_>mtH) ) return 0;
    if ( (region.Contains("dphiside")||region.Contains("dphireg")||region.Contains("dphicut")) && (dataEvent->dPhi_>dPhi*TMath::Pi()/180.) ) return 0;
    if ( (region.Contains("zregion")) && fabs(dataEvent->dilep_.mass()-91.1876)>7.5 ) return 0;
    if ( (region.Contains("zreg15")) && fabs(dataEvent->dilep_.mass()-91.1876)>15. ) return 0;
    if ( (region.Contains("zvetoall")) && fabs(dataEvent->dilep_.mass()-91.1876)<15. ) return 0;
    //additional higgs cuts
    if ( region.Contains("dphijet")   && ( dataEvent->type_!=1 && dataEvent->type_!=2 &&
					   ( (dataEvent->njets_<2 && dataEvent->jet1_.pt()>15&&dataEvent->dPhiDiLepJet1_>165.*TMath::Pi()/180.) || 
					     (dataEvent->njets_>=2 && fabs(ROOT::Math::VectorUtil::DeltaPhi((dataEvent->jet1_+dataEvent->jet2_),dataEvent->dilep_))*180.0/TMath::Pi() > 165.) )
					   ) ) return 0;
    if ( region.Contains("minmet40")  && ( dataEvent->type_!=1 && dataEvent->type_!=2 && min(dataEvent->pmet_,dataEvent->pTrackMet_)<40.0)) return 0;
    if ( region.Contains("dpjallfs")  && ( (dataEvent->njets_<2 && dataEvent->jet1_.pt()>15 && dataEvent->dPhiDiLepJet1_>165.*TMath::Pi()/180. ) ||
					   (dataEvent->njets_>=2 && fabs(ROOT::Math::VectorUtil::DeltaPhi((dataEvent->jet1_+dataEvent->jet2_),dataEvent->dilep_))*180.0/TMath::Pi() > 165.) ) ) return 0;
    if ( region.Contains("mm40allfs") && min(dataEvent->pmet_,dataEvent->pTrackMet_)<40.0 ) return 0;
    if ( region.Contains("minmetvtx") && ( dataEvent->type_!=1 && dataEvent->type_!=2 && min(dataEvent->pmet_,dataEvent->pTrackMet_)<(37.+dataEvent->nvtx_/2.)) ) return 0;
    if ( region.Contains("mmvtxallfs") && ( min(dataEvent->pmet_,dataEvent->pTrackMet_)<(37.+dataEvent->nvtx_/2.)) ) return 0;
    if ( region.Contains("met2040") && dataEvent->type_!=1 && dataEvent->type_!=2 && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<20.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>40.0) ) return 0;
    if ( region.Contains("lep2pt15") && ( dataEvent->type_!=1 && dataEvent->type_!=2 && dataEvent->lep2_.pt()<15.) ) return 0;
    if ( region.Contains("lep2pt15allfs") && (dataEvent->lep2_.pt()<15.) ) return 0;
    if ( region.Contains("ptll45") && ( dataEvent->dilep_.pt()<45.) ) return 0;
    //cuts for Rout/in
    if ( region.Contains("met2025") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<20.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>25.0) ) return 0;
    if ( region.Contains("met2530") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<25.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>30.0) ) return 0;
    if ( region.Contains("met3037") && (min(dataEvent->pmet_,dataEvent->pTrackMet_)<30.0||min(dataEvent->pmet_,dataEvent->pTrackMet_)>37.0) ) return 0;
    if ( region.Contains("met37up") &&  min(dataEvent->pmet_,dataEvent->pTrackMet_)<37.0 ) return 0;
    //cuts on btag
    if ( region.Contains("btag1or2")  && dataEvent->jet1Btag_<2.1 && dataEvent->jet2Btag_<2.1 ) return 0;
    if ( region.Contains("btag1and2") && (dataEvent->jet1Btag_<2.1 || dataEvent->jet2Btag_<2.1) ) return 0;
    if ( region.Contains("btagJet1")  && dataEvent->jet1Btag_<2.1 ) return 0;
    if ( region.Contains("btagJet2")  && dataEvent->jet2Btag_<2.1 ) return 0;
    if ( region.Contains("nobJet1")   && dataEvent->jet1Btag_>2.1 ) return 0;
    if ( region.Contains("nobJet2")   && dataEvent->jet2Btag_>2.1 ) return 0;
    //check peaking at MC level
    if ( region.Contains("fromZ") && (dataEvent->lep1MotherMcId_!=23 || dataEvent->lep2MotherMcId_!=23) ) return 0;
    if ( region.Contains("notZ") && !(dataEvent->lep1MotherMcId_!=23 || dataEvent->lep2MotherMcId_!=23) ) return 0;
    //spillage
    if ( region.Contains("spill") && ( !( abs(dataEvent->lep1McId_)==11 || abs(dataEvent->lep1McId_)==13 ) || !( abs(dataEvent->lep2McId_)==11 || abs(dataEvent->lep2McId_)==13 ) ) ) return 0;

    //mll>20 cut
    if ( region.Contains("mll20") && dataEvent->type_!=1 && dataEvent->type_!=2 && dataEvent->dilep_.mass()<20.) return 0;
    //higgs processId
    if ( region.Contains("ggH") && ( dataEvent->processId_!=10010 ) ) return 0;
    if ( region.Contains("qqH") && ( dataEvent->processId_!=10001 ) ) return 0;
    if ( region.Contains("ZH") && ( dataEvent->processId_!=24 ) ) return 0;
    if ( region.Contains("WH") && ( dataEvent->processId_!=26 ) ) return 0;

    return 1;  
}


float getFR(SmurfTree *dataEvent, TH2F* eFR, TH2F* mFR) {
  //FIXME: should be gotten from files
  float maxPt=34.9, maxEta=2.4;
  if ( (dataEvent->cuts_ & Lep1LooseEleV4)    == Lep1LooseEleV4    &&
       (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection && 
       (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) {
    float pt = dataEvent->lep1_.pt();
    float eta = fabs(dataEvent->lep1_.eta());
    if (eta==2) eta=1.9999;//test to sync with G's weights
    if (pt>maxPt) pt=maxPt;
    if (eta>maxEta) eta=maxEta;
    return eFR->GetBinContent(eFR->FindBin(pt,eta));
  } else if ( (dataEvent->cuts_ & Lep2LooseEleV4)    == Lep2LooseEleV4    &&
	      (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection &&
	      (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) {
    float pt = dataEvent->lep2_.pt();
    float eta = fabs(dataEvent->lep2_.eta());
    if (eta==2) eta=1.9999;//test to sync with G's weights
    if (pt>maxPt) pt=maxPt;
    if (eta>maxEta) eta=maxEta;
    return eFR->GetBinContent(eFR->FindBin(pt,eta));
  } else if ( (dataEvent->cuts_ & Lep1LooseMuV2)     == Lep1LooseMuV2     &&
	      (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection &&
	      (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) {
    float pt = dataEvent->lep1_.pt();
    float eta = fabs(dataEvent->lep1_.eta());
    if (eta==2) eta=1.9999;//test to sync with G's weights
    if (pt>maxPt) pt=maxPt;
    if (eta>maxEta) eta=maxEta;
    return mFR->GetBinContent(eFR->FindBin(pt,eta));
  } else if ( (dataEvent->cuts_ & Lep2LooseMuV2)     == Lep2LooseMuV2     &&
	      (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection &&
	      (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) {
    float pt = dataEvent->lep2_.pt();
    float eta = fabs(dataEvent->lep2_.eta());
    if (eta==2) eta=1.9999;//test to sync with G's weights
    if (pt>maxPt) pt=maxPt;
    if (eta>maxEta) eta=maxEta;
    return mFR->GetBinContent(eFR->FindBin(pt,eta));
  }
  return 0;
}

pair<float, float> getYield(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString region, float lumi, 
			    bool useJson=0, bool applyEff=true, bool doFake=false, bool doPUw=false) {

//   cout << sample << " " << cut << " " << veto << " " << mass << " " << njets << " " << region << " " << lumi << " " 
//        << useJson << " " << applyEff << " " << doFake << " " << doPUw << endl;

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);

  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");

  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;

  //PU reweighting
  TFile* puwf=0;
  TH1F* puweights=0;
  if (isMC&&doPUw&&redoWeights) {
    puwf = TFile::Open(puw_file);
    puweights = (TH1F*) puwf->Get("puWeights");
  }

  //fake stuff
  bool useMit = true;
  TFile *fE=0, *fM=0;
  TH2F *eFR=0,*mFR=0;
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
  }

  if (!isMC && useJson && jsonFile!=""){
    if (jsonFile.Contains(".txt")) set_goodrun_file(jsonFile);
    else set_goodrun_file_json(jsonFile);
  } else if (jsonFile=="") useJson=false;

  float weight = 1.;
  float yield = 0.;
  float error = 0.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    TString dataSetName(dataEvent->name(dataEvent->dstype_).c_str());
    if (dataSetName.Contains("hww")) weight*=dataEvent->sfWeightHPt_;
    if (region.Contains("embed")) {
      weight=lumi*ZttScaleFactor(dataEvent->nvtx_,2,dataEvent->scale1fb_);//fixme period
    }
 
    if (!passEvent(dataEvent, njets, cut, veto, region,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass,isMC, useJson)) continue;

//     cout << "run, lumi, evt: " << dataEvent->run_ << " " << dataEvent->lumi_ << " " << dataEvent->event_;
//     cout << " - type, nvtx, njets: " << dataEvent->type_ << " " << dataEvent->nvtx_ << " " << dataEvent->njets_;
//     cout << " - dphi: " << " " << ROOT::Math::VectorUtil::DeltaPhi((dataEvent->jet1_+dataEvent->jet2_),dataEvent->dilep_)*180.0/TMath::Pi();
//     cout << " - pt1, pt2, pmet: " << dataEvent->lep1_.pt() << " " << dataEvent->lep2_.pt() << " " << dataEvent->pmet_;
//     cout << " - mll, ptll, mt: " << dataEvent->dilep_.mass() << " " << dataEvent->dilep_.pt() << " " << dataEvent->mt_ << endl;
//     cout << endl;

    float puw = 1.;
    if (isMC&&doPUw) {
      ////puw = getPileupReweightFactor(dataEvent->nvtx_,puweights);
      if (redoWeights) puw = getPileupReweightFactor(dataEvent->npu_,puweights);
      else puw = dataEvent->sfWeightPU_;
      if (checkWeights) assert(fabs((getPileupReweightFactor(dataEvent->npu_,puweights)-dataEvent->sfWeightPU_)/dataEvent->sfWeightPU_)<0.0001);
      //cout << puw << " " << dataEvent->sfWeightPU_ << endl;
    }

    if (!doFake) {
      float effSF=1.;
      if (isMC&&applyEff) {
	if (redoWeights) {
	  float pt1 = min(dataEvent->lep1_.pt(),49.);
	  float eta1 = min(fabs(dataEvent->lep1_.eta()),2.4);
	  float pt2 = min(dataEvent->lep2_.pt(),49.);
	  float eta2 = min(fabs(dataEvent->lep2_.eta()),2.4);
	  LeptonScaleLookup lsl(&*eff_file);
	  float effSelSF = lsl.GetExpectedLeptonSF(eta1,pt1,dataEvent->lid1_)*lsl.GetExpectedLeptonSF(eta2,pt2,dataEvent->lid2_);
	  float effTrgSF = lsl.GetExpectedTriggerEfficiency(eta1,pt1,eta2,pt2,dataEvent->lid1_,dataEvent->lid2_);
	  effSF = effSelSF*effTrgSF;
	  if (checkWeights) {
	    if (dataEvent->sfWeightEff_>0) assert(fabs((effSelSF-dataEvent->sfWeightEff_)/dataEvent->sfWeightEff_)<0.0001);
	    assert(fabs((effTrgSF-dataEvent->sfWeightTrig_)/dataEvent->sfWeightTrig_)<0.0001);
	  }
	} else {
	  effSF = dataEvent->sfWeightEff_ * dataEvent->sfWeightTrig_;
	}
      }
      yield = yield + weight*effSF*puw;
      //fixme: should consider the error of effSF
      error = error + pow(weight*effSF*puw,2);
    } 
    else {//DO FAKES!!!      
      if ( ( (dataEvent->cuts_ & Lep1LooseEleV4)    == Lep1LooseEleV4    &&
	     (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection && 
	     (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) ||
	   ( (dataEvent->cuts_ & Lep2LooseEleV4)    == Lep2LooseEleV4    &&
	     (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection && 
	     (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) ||
	   ( (dataEvent->cuts_ & Lep1LooseMuV2)     == Lep1LooseMuV2     &&
	     (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection &&
	     (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) ||
	   ( (dataEvent->cuts_ & Lep2LooseMuV2)     == Lep2LooseMuV2     &&
	     (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection &&
	     (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) ) {
	float frW = 1;
	if (redoWeights || region.Contains("alternativeFR")) {
	  float fr = getFR(dataEvent,eFR,mFR);
	  frW = fr/(1.-fr);
	  if (checkWeights) assert((fabs(fr/(1.-fr))-fabs(dataEvent->sfWeightFR_))/fabs(dataEvent->sfWeightFR_)<0.0001);
	} else frW = fabs(dataEvent->sfWeightFR_);//warning, using fabs, does not work with mixed samples!!!!!
	yield += weight*puw*frW;
	error += pow(weight*puw*frW,2);//warning! not including error on FR here...
      } 
    }
  }
  error = sqrt(error);
  //if (isMC&&applyEff&&redoWeights) effs->Close();
  if (isMC&&doPUw&&redoWeights) puwf->Close();
  if (doFake) {
    delete eFR;
    delete mFR;
    fE->Close();
    fM->Close();
  }
  dataEvent->tree_->Delete();
  delete dataEvent;
  //cout << yield << " " << error << endl;
  return make_pair<float, float>(yield,error);
}

void fillPlot(ReadBDTG* rbdtg, TH1F* h, TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString region, float lumi, 
	      bool useJson=0, bool applyEff=true, bool doFake=false, bool doPUw=false, TString syst="") {

//   cout << sample << " " << cut << " " << veto << " " << mass << " " << njets << " " << region << " " << lumi << " " 
//        << useJson << " " << applyEff << " " << doFake << " " << doPUw << endl;

  h->Sumw2();

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);

  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");

  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;

  //PU reweighting
  TFile* puwf=0;
  TH1F* puweights=0;
  if (isMC&&doPUw&&redoWeights) {
    puwf = TFile::Open(puw_file);
    puweights = (TH1F*) puwf->Get("puWeights");
  }

  //fake stuff
  bool useMit = true;
  TFile *fE=0, *fM=0;
  TH2F *eFR=0,*mFR=0;
  TH2F *elPFY=0,*muPFY=0;
  if (doFake) {
    if (useMit) {
      fE = TFile::Open(fr_file_mit);
      if (region.Contains("alternativeFR")==0) eFR = (TH2F*) fE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta");//default
      else eFR = (TH2F*) fE->Get("ElectronFakeRate_V4_ptThreshold50_PtEta");
      fM = TFile::Open(fr_file_mit);
      if (region.Contains("alternativeFR")==0) mFR = (TH2F*) fM->Get("MuonFakeRate_M2_ptThreshold15_PtEta");//default
      else  mFR = (TH2F*) fM->Get("MuonFakeRate_M2_ptThreshold30_PtEta");
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

  if (!isMC && useJson && jsonFile!=""){
    if (jsonFile.Contains(".txt")) set_goodrun_file(jsonFile);
    else set_goodrun_file_json(jsonFile);
  } else if (jsonFile=="") useJson=false;

  float weight = 1.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    TString dataSetName(dataEvent->name(dataEvent->dstype_).c_str());
    if (dataSetName.Contains("hww")) weight*=dataEvent->sfWeightHPt_;
    if (region.Contains("embed")) {
      weight=lumi*ZttScaleFactor(dataEvent->nvtx_,2,dataEvent->scale1fb_);//fixme period
    }

    if (!passEvent(dataEvent, njets, cut, veto, region,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass,isMC, useJson)) continue;

    //change dataEvent for syst studies
    //fixme should go before passEvent: it affects the cuts also, but this is to reproduce Guillelmo's result
    if (syst=="metSmear") {//met variation
      double metx=0.0;double mety=0.0;double trkmetx=0.0;double trkmety=0.0;
      if	(dataEvent->njets_ == 0){
	metx    = dataEvent->met_*cos(dataEvent->metPhi_)+gRandom->Gaus(0.0,3.2);
	mety    = dataEvent->met_*sin(dataEvent->metPhi_)+gRandom->Gaus(0.0,3.2);
	trkmetx = dataEvent->trackMet_*cos(dataEvent->trackMetPhi_)+gRandom->Gaus(0.0,2.1);
	trkmety = dataEvent->trackMet_*sin(dataEvent->trackMetPhi_)+gRandom->Gaus(0.0,2.1);
      }
      else if(dataEvent->njets_ == 1){
	metx    = dataEvent->met_*cos(dataEvent->metPhi_)+gRandom->Gaus(0.0,3.6);
	mety    = dataEvent->met_*sin(dataEvent->metPhi_)+gRandom->Gaus(0.0,3.6);
	trkmetx = dataEvent->trackMet_*cos(dataEvent->trackMetPhi_)+gRandom->Gaus(0.0,7.6);
	trkmety = dataEvent->trackMet_*sin(dataEvent->trackMetPhi_)+gRandom->Gaus(0.0,7.6);
      }
      else if(dataEvent->njets_ >= 2){
	metx    = dataEvent->met_*cos(dataEvent->metPhi_)+gRandom->Gaus(0.0,4.3);
	mety    = dataEvent->met_*sin(dataEvent->metPhi_)+gRandom->Gaus(0.0,4.3);
	trkmetx = dataEvent->trackMet_*cos(dataEvent->trackMetPhi_)+gRandom->Gaus(0.0,12.4);
	trkmety = dataEvent->trackMet_*sin(dataEvent->trackMetPhi_)+gRandom->Gaus(0.0,12.4);
      }
      double newMet      = sqrt(metx*metx+mety*mety);
      double newTrackMet = sqrt(trkmetx*trkmetx+trkmety*trkmety);
      double deltaPhiA[3] = {TMath::Abs(dataEvent->lep1_.Phi()-TMath::ATan2(mety,metx)),TMath::Abs(dataEvent->lep2_.Phi()-TMath::ATan2(mety,metx)),0.0};
      while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
      while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
      deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
      double pmetA = newMet;
      if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);      
      double deltaPhiB[3] = {TMath::Abs(dataEvent->lep1_.Phi()-TMath::ATan2(trkmety,trkmetx)),TMath::Abs(dataEvent->lep2_.Phi()-TMath::ATan2(trkmety,trkmetx)),0.0};
      while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
      while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
      deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
      double pmetB = newTrackMet;
      if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);
      double oldMet = dataEvent->met_;      
      dataEvent->pmet_ = pmetA; 
      dataEvent->met_  = newMet; 
      dataEvent->pTrackMet_ = pmetB; 
      dataEvent->trackMet_  = newTrackMet; 
      dataEvent->mt_  = dataEvent->mt_*sqrt(newMet/oldMet); 
      dataEvent->mt1_ = dataEvent->mt1_*sqrt(newMet/oldMet);
      dataEvent->mt2_ = dataEvent->mt2_*sqrt(newMet/oldMet);      
      dataEvent->dPhiLep1MET_ = TMath::Abs(dataEvent->lep1_.phi()-TMath::ATan2(mety,metx));
      while(dataEvent->dPhiLep1MET_>TMath::Pi()) dataEvent->dPhiLep1MET_ = TMath::Abs(dataEvent->dPhiLep1MET_ - 2*TMath::Pi());
      dataEvent->dPhiLep2MET_ = TMath::Abs(dataEvent->lep2_.phi()-TMath::ATan2(mety,metx));
      while(dataEvent->dPhiLep2MET_>TMath::Pi()) dataEvent->dPhiLep2MET_ = TMath::Abs(dataEvent->dPhiLep2MET_ - 2*TMath::Pi());
      dataEvent->dPhiDiLepMET_ = TMath::Abs(dataEvent->dilep_.phi()-TMath::ATan2(mety,metx));
      while(dataEvent->dPhiDiLepMET_>TMath::Pi()) dataEvent->dPhiDiLepMET_ = TMath::Abs(dataEvent->dPhiDiLepMET_ - 2*TMath::Pi());

    } else if (syst.Contains("momScale")) {//momentum scale and resolution up/down
      double corr[2] = {1.0, 1.0};
      if (syst.Contains("Up")) {
	if (TMath::Abs(dataEvent->lid1_) == 13){
	  corr[0] = 1.01 + gRandom->Gaus(0.00,0.01);
	} else if (TMath::Abs(dataEvent->lid1_) == 11 && TMath::Abs(dataEvent->lep1_.eta()) <  1.479){
	  corr[0] = 1.01 + gRandom->Gaus(0.00,0.02);
	} else if (TMath::Abs(dataEvent->lid1_) == 11 && TMath::Abs(dataEvent->lep1_.eta()) >= 1.479){
	  corr[0] = 1.06 + gRandom->Gaus(0.00,0.06);
	}
	if (TMath::Abs(dataEvent->lid2_) == 13){
	  corr[1] = 1.01 + gRandom->Gaus(0.00,0.01);
	} else if (TMath::Abs(dataEvent->lid2_) == 11 && TMath::Abs(dataEvent->lep2_.eta()) <  1.479){
	  corr[1] = 1.01 + gRandom->Gaus(0.00,0.02);
	} else if (TMath::Abs(dataEvent->lid2_) == 11 && TMath::Abs(dataEvent->lep2_.eta()) >= 1.479){
	  corr[1] = 1.06 + gRandom->Gaus(0.00,0.06);
	}
      } else if (syst.Contains("Down")) {
	if (TMath::Abs(dataEvent->lid1_) == 13){
	  corr[0] = 0.99 - gRandom->Gaus(0.00,0.01);
	} else if (TMath::Abs(dataEvent->lid1_) == 11 && TMath::Abs(dataEvent->lep1_.eta()) <  1.479){
	  corr[0] = 0.99 - gRandom->Gaus(0.00,0.02);
	} else if (TMath::Abs(dataEvent->lid1_) == 11 && TMath::Abs(dataEvent->lep1_.eta()) >= 1.479){
	  corr[0] = 0.94 - gRandom->Gaus(0.00,0.06);
	}
	if(TMath::Abs(dataEvent->lid2_) == 13){
	  corr[1] = 0.99 - gRandom->Gaus(0.00,0.01);
	} else if (TMath::Abs(dataEvent->lid2_) == 11 && TMath::Abs(dataEvent->lep2_.eta()) <  1.479){
	  corr[1] = 0.99 - gRandom->Gaus(0.00,0.02);
	} else if (TMath::Abs(dataEvent->lid2_) == 11 && TMath::Abs(dataEvent->lep2_.eta()) >= 1.479){
	  corr[1] = 0.94 - gRandom->Gaus(0.00,0.06);
	}
      } else {
	cout << "Wrong name for momScale systematic! Should contain Up or Down" << endl;
	return;
      }
//       cout << "new evt" << endl;
//       cout << dataEvent->lep1_.pt() << " " << dataEvent->lep2_.pt() << " " << dataEvent->dilep_.mass() << " " << dataEvent->dPhi_ << " "
// 	   << dataEvent->mt_ << " " << dataEvent->mt1_ << " " << dataEvent->mt2_ << " " << dataEvent->dPhiDiLepMET_ << " " << dataEvent->dPhiDiLepJet1_
// 	   << endl;
      SmurfTree::LorentzVector lep1Old = dataEvent->lep1_;
      SmurfTree::LorentzVector lep2Old = dataEvent->lep2_;
      SmurfTree::LorentzVector dilepOld = dataEvent->dilep_;
      SmurfTree::LorentzVector lep1New = dataEvent->lep1_*corr[0];
      SmurfTree::LorentzVector lep2New = dataEvent->lep2_*corr[1];
      SmurfTree::LorentzVector dilepNew = lep1New+lep2New;

      dataEvent->lep1_ = lep1New;
      dataEvent->lep2_ = lep2New;
      dataEvent->dilep_ = dilepNew;
      dataEvent->dPhi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(lep1New,lep2New));
      dataEvent->dR_ = ROOT::Math::VectorUtil::DeltaR(lep1New,lep2New);

      dataEvent->mt_  = dataEvent->mt_*sqrt(dilepNew.pt()/dilepOld.pt()); 
      dataEvent->mt1_ = dataEvent->mt1_*sqrt(corr[0]); 
      dataEvent->mt2_ = dataEvent->mt2_*sqrt(corr[1]); 

      dataEvent->dPhiDiLepMET_  = fabs(ROOT::Math::VectorUtil::DeltaPhi(dilepNew,SmurfTree::LorentzVector(dataEvent->met_*cos(dataEvent->metPhi_),dataEvent->met_*sin(dataEvent->metPhi_),0,dataEvent->met_)));
      if (dataEvent->dPhiDiLepJet1_<998) dataEvent->dPhiDiLepJet1_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(dilepNew,dataEvent->jet1_));
//       cout << dataEvent->lep1_.pt() << " " << dataEvent->lep2_.pt() << " " << dataEvent->dilep_.mass() << " " << dataEvent->dPhi_ << " "
// 	   << dataEvent->mt_ << " " << dataEvent->mt1_ << " " << dataEvent->mt2_ << " " << dataEvent->dPhiDiLepMET_ << " " << dataEvent->dPhiDiLepJet1_
// 	   << endl;
    }

//     cout << "run, lumi, evt: " << dataEvent->run_ << " " << dataEvent->lumi_ << " " << dataEvent->event_;
//     cout << " - type, nvtx, njets: " << dataEvent->type_ << " " << dataEvent->nvtx_ << " " << dataEvent->njets_;
//     cout << " - dphi: " << " " << ROOT::Math::VectorUtil::DeltaPhi((dataEvent->jet1_+dataEvent->jet2_),dataEvent->dilep_)*180.0/TMath::Pi();
//     cout << " - pt1, pt2, pmet: " << dataEvent->lep1_.pt() << " " << dataEvent->lep2_.pt() << " " << dataEvent->pmet_;
//     cout << " - mll, ptll, mt: " << dataEvent->dilep_.mass() << " " << dataEvent->dilep_.pt() << " " << dataEvent->mt_ << endl;
//     cout << endl;

    float puw = 1.;
    if (isMC&&doPUw) {
      ////puw = getPileupReweightFactor(dataEvent->nvtx_,puweights);
      if (redoWeights) puw = getPileupReweightFactor(dataEvent->npu_,puweights);
      else puw = dataEvent->sfWeightPU_;
      if (checkWeights) assert(fabs((getPileupReweightFactor(dataEvent->npu_,puweights)-dataEvent->sfWeightPU_)/dataEvent->sfWeightPU_)<0.0001);
      //cout << puw << " " << dataEvent->sfWeightPU_ << endl;
    }

    if (!doFake) {
      float effSF=1.;
      if (isMC&&applyEff) {
	if (redoWeights) {
	  float pt1 = min(dataEvent->lep1_.pt(),49.);
	  float eta1 = min(fabs(dataEvent->lep1_.eta()),2.4);
	  float pt2 = min(dataEvent->lep2_.pt(),49.);
	  float eta2 = min(fabs(dataEvent->lep2_.eta()),2.4);
	  LeptonScaleLookup lsl(&*eff_file);
	  float effSelSF = lsl.GetExpectedLeptonSF(eta1,pt1,dataEvent->lid1_)*lsl.GetExpectedLeptonSF(eta2,pt2,dataEvent->lid2_);
	  float effTrgSF = lsl.GetExpectedTriggerEfficiency(eta1,pt1,eta2,pt2,dataEvent->lid1_,dataEvent->lid2_);
	  effSF = effSelSF*effTrgSF;
	  if (checkWeights) {
	    if (dataEvent->sfWeightEff_>0) assert(fabs((effSelSF-dataEvent->sfWeightEff_)/dataEvent->sfWeightEff_)<0.0001);
	    assert(fabs((effTrgSF-dataEvent->sfWeightTrig_)/dataEvent->sfWeightTrig_)<0.0001);
	  }
	} else {
	  effSF = dataEvent->sfWeightEff_ * dataEvent->sfWeightTrig_;
	}
      }
      //h->Fill(dataEvent->dilep_.mass(),weight*effSF*puw);
      std::vector<double> theInputVals;
      if (njets==0) {
	const double inputVals[] = { dataEvent->lep1_.pt(), dataEvent->lep2_.pt(), dataEvent->dPhi_, dataEvent->dR_, dataEvent->dilep_.mass(), dataEvent->type_, dataEvent->mt_};
	for (int i=0;i<7;++i) theInputVals.push_back(inputVals[i]);
      } else if (njets==1) {
	const double inputVals[] = { dataEvent->lep1_.pt(), dataEvent->lep2_.pt(), dataEvent->dPhi_, dataEvent->dR_, dataEvent->dilep_.mass(), dataEvent->type_, dataEvent->mt_, dataEvent->dPhiDiLepMET_, dataEvent->dPhiDiLepJet1_};
	for (int i=0;i<9;++i) theInputVals.push_back(inputVals[i]);	
      } else {
	cout << "njets not supported: " << njets << endl;
	return;
      }
      h->Fill(rbdtg->GetMvaValue(theInputVals),weight*effSF*puw);
    } else {//DO FAKES!!!      
      if ( ( (dataEvent->cuts_ & Lep1LooseEleV4)    == Lep1LooseEleV4    &&
	     (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection && 
	     (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) ||
	   ( (dataEvent->cuts_ & Lep2LooseEleV4)    == Lep2LooseEleV4    &&
	     (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection && 
	     (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) ||
	   ( (dataEvent->cuts_ & Lep1LooseMuV2)     == Lep1LooseMuV2     &&
	     (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection &&
	     (dataEvent->cuts_ & Lep2FullSelection) == Lep2FullSelection ) ||
	   ( (dataEvent->cuts_ & Lep2LooseMuV2)     == Lep2LooseMuV2     &&
	     (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection &&
	     (dataEvent->cuts_ & Lep1FullSelection) == Lep1FullSelection ) ) {
	const double inputVals[] = { dataEvent->lep1_.pt(), dataEvent->lep2_.pt(), dataEvent->dPhi_, dataEvent->dR_, dataEvent->dilep_.mass(), dataEvent->type_, dataEvent->mt_};
	std::vector<double> theInputVals;
	for (int i=0;i<7;++i) theInputVals.push_back(inputVals[i]);
	float frW = 1;
	if (redoWeights || region.Contains("alternativeFR")) {
	  float fr = getFR(dataEvent,eFR,mFR);
	  frW = fr/(1.-fr);
	} else frW = fabs(dataEvent->sfWeightFR_);//warning, using fabs, does not work with mixed samples!!!!!
	h->Fill(rbdtg->GetMvaValue(theInputVals),weight*puw*frW);
      } 
    }
  }
  if (isMC&&doPUw&&redoWeights) puwf->Close();
  if (doFake) {
    delete eFR;
    delete mFR;
    delete muPFY;
    delete elPFY;
    fE->Close();
    fM->Close();
  }
  dataEvent->tree_->Delete();
  delete dataEvent;
}

#endif
