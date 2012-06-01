/*
from Yanyan Gao
*/

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRint.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCut.h"
#include <TSystem.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

// steal from ../../Core/SmurfTree.h
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

using namespace std;

//###################
//# main function
//###################
void skim(TString smurfFDir, TString fileName, TString outputDir, TString cutstring) {

  TFile* fin = new TFile(smurfFDir+fileName);
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 
  
  TString outputFileName = outputDir + fileName;
  if (cutstring == "PassFail")   outputFileName.ReplaceAll(".root","_PassFail.root");
  if (cutstring == "LooseMET")   outputFileName.ReplaceAll(".root","_LooseMET.root");
  
  TFile *newfile= new TFile(outputFileName,"recreate");
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // get event based branches..
  unsigned int njets_ = 0;
  unsigned int cuts_ = 0;
  LorentzVector*  lep1_ = 0;
  LorentzVector*  lep2_ = 0;
  LorentzVector*  dilep_ = 0;
  int type_ = 0; // 0/1/2/3 for mm/me/em/ee
  float pmet_ = 0.0;
  float pTrackMet_ = 0.0;
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  float mt_ = 0.;
  float dPhi_ = 0.;
  LorentzVector*  jet1_ = 0;
  float dPhiDiLepJet1_ = 0.;
  float met_ = 0.0;
  float trackMet_ = 0.0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  float jet1Btag_ = 0;
  float jet2Btag_ = 0;
  float jet3Btag_ = 0;
  unsigned int nSoftMuons_ = 0;
  
  ch->SetBranchAddress( "njets"     , &njets_     );     
  ch->SetBranchAddress( "cuts"      , &cuts_     );     
  ch->SetBranchAddress( "lep1"      , &lep1_      );   
  ch->SetBranchAddress( "lep2"      , &lep2_      );   
  ch->SetBranchAddress( "dilep"      , &dilep_      );   
  ch->SetBranchAddress( "type"      , &type_     );     
  ch->SetBranchAddress( "pmet"      , &pmet_     );     
  ch->SetBranchAddress( "pTrackMet"      , &pTrackMet_     );   
  ch->SetBranchAddress( "run"     , &run_     );     
  ch->SetBranchAddress( "lumi"     , &lumi_     );     
  ch->SetBranchAddress( "dPhi"      , &dPhi_     );     
  ch->SetBranchAddress( "mt"      , &mt_     );     
  ch->SetBranchAddress( "jet1"      , &jet1_      );   
  ch->SetBranchAddress( "jet2"      , &jet2_      );   
  ch->SetBranchAddress( "jet3"      , &jet3_      );   
  ch->SetBranchAddress( "dPhiDiLepJet1"      , &dPhiDiLepJet1_     );     
  ch->SetBranchAddress( "met"      , &met_     );     
  ch->SetBranchAddress( "trackMet"      , &trackMet_     );   
  ch->SetBranchAddress( "jet1Btag"      , &jet1Btag_      );   
  ch->SetBranchAddress( "jet2Btag"      , &jet2Btag_      );   
  ch->SetBranchAddress( "jet3Btag"      , &jet3Btag_      );   
  ch->SetBranchAddress( "nSoftMuons"      , &nSoftMuons_      );   

  float scale1fb = 0.0;
  ch->SetBranchAddress( "scale1fb"      , &scale1fb     );   

  unsigned int nvtx_ = 0;
  ch->SetBranchAddress( "nvtx"     , &nvtx_     );     
  unsigned int npu_ = 0;
  ch->SetBranchAddress( "npu"     , &npu_     );     

  float sfWeightPU_ = 1.;
  ch->SetBranchAddress("sfWeightPU",    &sfWeightPU_ );  

  
  //==========================================
  // Loop All Events
  //==========================================
  
  cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; \n";

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt); 

    if ( int(njets_) > 3 ) continue;

    if ( dilep_->mass() < 12.0) continue;
    
    //this is for dy
    if (cutstring=="dy") {
      if ((cuts_ & wwSelNoZVNoMet) != wwSelNoZVNoMet) continue;
      if (min(pmet_,pTrackMet_)<20.0) continue;
    }

    if (cutstring=="mm20") {
      unsigned int selBLCMELV    = BaseLine|ChargeMatch|ExtraLeptonVeto;
      if ((cuts_ & selBLCMELV) != selBLCMELV) continue;
      if (min(pmet_,pTrackMet_)<20.0) continue;
    }

    //this is for top/ww
    if (cutstring=="topww") {
      if ((cuts_ & wwSelNoLepNoTV) != wwSelNoLepNoTV) continue;
    }
    
    evt_tree->Fill();
  }   //nevent
  
  cout << outputDir + fileName << " has " << evt_tree->GetEntries() << " entries; \n";
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
}  

void skim(TString smurfFDir = "/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/", TString outputDir = "/smurf/cerati/skims/", TString cut = "topww") {
  if (cut!="dy" && cut!="topww" && cut!="mm20") {
    cout << "cut not supported. please use dy or topww or mm20" << endl;
    return; 
  }
  if (cut=="topww" && !outputDir.Contains("wwSelNoLepNoTV")) outputDir+="/wwSelNoLepNoTV/";
  else if (cut=="dy" && !outputDir.Contains("wwSelNoMetNoZVminMET20")) outputDir+="/wwSelNoMetNoZVminMET20/";
  gSystem->Exec("mkdir -p "+outputDir);
  skim(smurfFDir,"data.root",outputDir,cut);
  skim(smurfFDir,"dyll.root",outputDir,cut);
  skim(smurfFDir,"ggww.root",outputDir,cut);
  skim(smurfFDir,"hww110.root",outputDir,cut);
  skim(smurfFDir,"hww115.root",outputDir,cut);
  skim(smurfFDir,"hww120.root",outputDir,cut);
  skim(smurfFDir,"hww125.root",outputDir,cut);
  skim(smurfFDir,"hww130.root",outputDir,cut);
  skim(smurfFDir,"hww135.root",outputDir,cut);
  skim(smurfFDir,"hww140.root",outputDir,cut);
  skim(smurfFDir,"hww145.root",outputDir,cut);
  skim(smurfFDir,"hww150.root",outputDir,cut);
  skim(smurfFDir,"hww160.root",outputDir,cut);
  skim(smurfFDir,"hww170.root",outputDir,cut);
  skim(smurfFDir,"hww180.root",outputDir,cut);
  skim(smurfFDir,"hww190.root",outputDir,cut);
  skim(smurfFDir,"hww200.root",outputDir,cut);
  skim(smurfFDir,"hww250.root",outputDir,cut);
  skim(smurfFDir,"hww300.root",outputDir,cut);
  skim(smurfFDir,"hww350.root",outputDir,cut);
  skim(smurfFDir,"hww400.root",outputDir,cut);
  skim(smurfFDir,"hww450.root",outputDir,cut);
  skim(smurfFDir,"hww500.root",outputDir,cut);
  skim(smurfFDir,"hww550.root",outputDir,cut);
  skim(smurfFDir,"hww600.root",outputDir,cut);
  skim(smurfFDir,"qqww_py.root",outputDir,cut);
  skim(smurfFDir,"qqww.root",outputDir,cut);
  skim(smurfFDir,"ttbar.root",outputDir,cut);
  skim(smurfFDir,"tw.root",outputDir,cut);
  skim(smurfFDir,"wjets.root",outputDir,cut);
  skim(smurfFDir,"wz.root",outputDir,cut);
  skim(smurfFDir,"zz.root",outputDir,cut);
  return;
}

/*
backgroundA.root        data_skim2.root  hww110.root  hww125.root  hww140.root  hww160.root  hww190.root  hww300.root  hww450.root  hww600.root          qqww_py.root  tw.root     zz.root
backgroundA_skim2.root  dyll.root        hww115.root  hww130.root  hww145.root  hww170.root  hww200.root  hww350.root  hww500.root  hww_syst.root        qqww.root     wjets.root
data.root               ggww.root        hww120.root  hww135.root  hww150.root  hww180.root  hww250.root  hww400.root  hww550.root  hww_syst_skim2.root  ttbar.root    wz.root
*/
