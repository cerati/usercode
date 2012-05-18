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

#include "/smurf/cerati/weightsDYMVA/TMVA_BDTG_0j_mll12.class.C"
#include "/smurf/cerati/weightsDYMVA/TMVA_BDTG_1j_mll12.class.C"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

using namespace std;

float recoilvar(float met, float metPhi, LorentzVector* dilep)
{
  float px = met*cos(metPhi) + dilep->px();       
  float py = met*sin(metPhi) + dilep->py();
  return sqrt(px*px+py*py);
}

//###################
//# main function
//###################
void addDYMVA(TString smurfFDir, TString fileName, TString outputDir) {

  TFile* fin = new TFile(smurfFDir+fileName);
  TTree* ch=(TTree*)fin->Get("tree"); 
  if (ch==0x0) return; 

  TString outputFileName = outputDir + fileName;  
  TFile *newfile= new TFile(outputFileName,"recreate");
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

  std::vector<std::string> theInputVars;
  const char* inputVars[] = { "pmet", "pTrackMet","nvtx", "dilpt", "jet1pt", "metSig", "dPhiDiLepJet1", "dPhiDiLepMET", "dPhiJet1MET", "recoil", "mt" };
  for (int i=0;i<11;++i) theInputVars.push_back(inputVars[i]);
  dymva_0j_Zveto::ReadBDTG* rbdtgDy_0j = new dymva_0j_Zveto::ReadBDTG(theInputVars);
  dymva_1j_Zveto::ReadBDTG* rbdtgDy_1j = new dymva_1j_Zveto::ReadBDTG(theInputVars);
  
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
  float dPhiDiLepMET_ = 0.;
  float met_ = 0.0;
  float trackMet_ = 0.0;
  LorentzVector*  jet2_ = 0;
  LorentzVector*  jet3_ = 0;
  float jet1Btag_ = 0;
  float jet2Btag_ = 0;
  float jet3Btag_ = 0;
  unsigned int nSoftMuons_ = 0;
  float scale1fb = 0.0;
  unsigned int nvtx_ = 0;
  unsigned int npu_ = 0;
  float sfWeightPU_ = 1.;
  float sfWeightFR_ = 1.;
  float sfWeightTrig_ = 1.;
  float sfWeightEff_ = 1.;
  int   lid1_;
  int   lid2_;
  float sumet_ = 0.0;
  //float metSig_ = 0.0;
  float metPhi_ = 0.0;
  int dstype_ = 0;
  int processId_ = 0;
  float higgsPt_ = 0.0;
  float sfWeightHPt_ = 1.;
  
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
  ch->SetBranchAddress( "dPhiDiLepMET"      , &dPhiDiLepMET_     );     
  ch->SetBranchAddress( "met"      , &met_     );     
  ch->SetBranchAddress( "trackMet"      , &trackMet_     );   
  ch->SetBranchAddress( "jet1Btag"      , &jet1Btag_      );   
  ch->SetBranchAddress( "jet2Btag"      , &jet2Btag_      );   
  ch->SetBranchAddress( "jet3Btag"      , &jet3Btag_      );   
  ch->SetBranchAddress( "nSoftMuons"      , &nSoftMuons_      );   
  ch->SetBranchAddress( "scale1fb"      , &scale1fb     );   
  ch->SetBranchAddress( "nvtx"     , &nvtx_     );     
  ch->SetBranchAddress( "npu"     , &npu_     );     
  ch->SetBranchAddress("sfWeightPU",    &sfWeightPU_ );  
  ch->SetBranchAddress("sfWeightFR",    &sfWeightFR_ );  
  ch->SetBranchAddress("sfWeightTrig",    &sfWeightTrig_ );  
  ch->SetBranchAddress("sfWeightEff",    &sfWeightEff_ );    
  ch->SetBranchAddress( "lid1"      , &lid1_      );   
  ch->SetBranchAddress( "lid2"      , &lid2_      );   
  ch->SetBranchAddress( "sumet"      , &sumet_     );     
  ch->SetBranchAddress( "metPhi"      , &metPhi_     );     
  ch->SetBranchAddress( "dstype"      , &dstype_     );     
  ch->SetBranchAddress( "processId"      , &processId_     );     
  ch->SetBranchAddress( "higgsPt"      , &higgsPt_     );     
  ch->SetBranchAddress( "sfWeightHPt"      , &sfWeightHPt_     );     

  //==========================================
  // Loop All Events
  //==========================================  
  cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; \n";

  float dymva_= -999.;
  evt_tree->Branch("dymva", &dymva_, "dymva/F");
  float recoil_= -999.;
  evt_tree->Branch("recoil", &recoil_, "recoil/F");
  float dPhiJet1MET_= -999.;
  evt_tree->Branch("dPhiJet1MET", &dPhiJet1MET_, "dPhiJet1MET/F");

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt); 

    LorentzVector metlv( met_*cos(metPhi_), met_*sin(metPhi_), 0, met_ );
    assert((metlv.phi()-metPhi_)<0.001);
    assert((metlv.pt()-met_)<0.001);

    recoil_ = recoilvar(met_, metPhi_, dilep_);
    dPhiJet1MET_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(*jet1_,metlv));

    //const char* inputVars[] = { "pmet","pTrackMet","nvtx","dilpt","jet1pt","metSig","dPhiDiLepJet1","dPhiDiLepMET","dPhiJet1MET","recoil","mt" };
    std::vector<double> theInputVals;
    const double inputVals[] = { pmet_,pTrackMet_,nvtx_,dilep_->pt(),
				 std::max(15.,jet1_->pt()),
				 met_/sqrt(sumet_),
				 (jet1_->pt()<15. ? -0.1 : dPhiDiLepJet1_ ),
				 dPhiDiLepMET_,
				 (jet1_->pt()<15. ? -0.1 : dPhiJet1MET_ ),
				 recoil_,
				 mt_};
    for (int i=0;i<11;++i) theInputVals.push_back(inputVals[i]);
    if (njets_==0) dymva_ = rbdtgDy_0j->GetMvaValue(theInputVals);
    else if (njets_==1) dymva_ = rbdtgDy_1j->GetMvaValue(theInputVals);
    else dymva_= -999.;

    evt_tree->Fill();
  }  
  cout << outputDir + fileName << " has " << evt_tree->GetEntries() << " entries; \n";
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
}  

void addDYMVA(TString smurfFDir, TString outputDir) {
  gSystem->Exec("mkdir -p "+outputDir);
  addDYMVA(smurfFDir,"dymm.root",outputDir);
  addDYMVA(smurfFDir,"dyee.root",outputDir);
  addDYMVA(smurfFDir,"dymm_mg.root",outputDir);
  addDYMVA(smurfFDir,"dyee_mg.root",outputDir);
  addDYMVA(smurfFDir,"dytt.root",outputDir);
  addDYMVA(smurfFDir,"qqww.root",outputDir);
  addDYMVA(smurfFDir,"ggww.root",outputDir);
  addDYMVA(smurfFDir,"ttbar.root",outputDir);
  addDYMVA(smurfFDir,"tw.root",outputDir);
  addDYMVA(smurfFDir,"wjets.root",outputDir);
  addDYMVA(smurfFDir,"wgamma.root",outputDir);
  addDYMVA(smurfFDir,"wz.root",outputDir);
  addDYMVA(smurfFDir,"zz_py.root",outputDir);
  addDYMVA(smurfFDir,"data.root",outputDir);
  addDYMVA(smurfFDir,"hww115.root",outputDir);
  addDYMVA(smurfFDir,"hww120.root",outputDir);
  addDYMVA(smurfFDir,"hww130.root",outputDir);
  addDYMVA(smurfFDir,"hww140.root",outputDir);
  addDYMVA(smurfFDir,"hww150.root",outputDir);
  addDYMVA(smurfFDir,"hww160.root",outputDir);
  addDYMVA(smurfFDir,"hww170.root",outputDir);
  addDYMVA(smurfFDir,"hww180.root",outputDir);
  addDYMVA(smurfFDir,"hww190.root",outputDir);
  addDYMVA(smurfFDir,"hww200.root",outputDir);
  addDYMVA(smurfFDir,"hww250.root",outputDir);
  addDYMVA(smurfFDir,"hww300.root",outputDir);
  addDYMVA(smurfFDir,"ttbar_mg.root",outputDir);
  addDYMVA(smurfFDir,"tw_ds.root",outputDir);
  addDYMVA(smurfFDir,"ww_mcnlo.root",outputDir);
  addDYMVA(smurfFDir,"ww_mcnlo_up.root",outputDir);
  addDYMVA(smurfFDir,"ww_mcnlo_down.root",outputDir);
  addDYMVA(smurfFDir,"wg3l.root",outputDir);
  addDYMVA(smurfFDir,"data-emb-tau123.root",outputDir);
  addDYMVA(smurfFDir,"data-emb-tau123.root",outputDir);
  addDYMVA(smurfFDir,"mc_lfake.root",outputDir);

  //addDYMVA(smurfFDir,"wz_py.root",outputDir);
  //addDYMVA(smurfFDir,"qqww_py.root",outputDir);

  //addDYMVA(smurfFDir,"data-emb-tau122.root",outputDir);
  //addDYMVA(smurfFDir,"data-emb-tau121.root",outputDir);
  return;
}
