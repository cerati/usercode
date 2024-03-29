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

#include "/smurf/cerati/weightsDYMVA/TMVA_0j_metshift_BDTG.class.C"
#include "/smurf/cerati/weightsDYMVA/TMVA_1j_metshift_BDTG.class.C"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

using namespace std;

float recoilvar(float met, float metPhi, LorentzVector* dilep)
{
  float px = met*cos(metPhi) + dilep->px();       
  float py = met*sin(metPhi) + dilep->py();
  return sqrt(px*px+py*py);
}

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

float deltaPhi( float phi1 , float phi2 ) {
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

double projectedMet(double phiL1, double phiL2, double met, double metPhi)
{
  double tightDPhi = fabs(deltaPhi(phiL1,metPhi));
  double looseDPhi = fabs(deltaPhi(phiL2,metPhi));
  double DeltaPhi = TMath::Min(tightDPhi, looseDPhi);
  if (DeltaPhi < TMath::Pi()/2) return met*TMath::Sin(DeltaPhi);
  return met;
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
  dymva_0j_ms::ReadBDTG* rbdtgDy_0j = new dymva_0j_ms::ReadBDTG(theInputVars);
  dymva_1j_ms::ReadBDTG* rbdtgDy_1j = new dymva_1j_ms::ReadBDTG(theInputVars);
  
  // get event based branches..
  unsigned int njets_ = 0;
  unsigned int cuts_ = 0;
  LorentzVector*  lep1_ = 0;
  LorentzVector*  lep2_ = 0;
  LorentzVector*  dilep_ = 0;
  int type_ = 0; // 0/1/2/3 for mm/me/em/ee
  float pTrackMet_ = 0.0;
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  float dPhi_ = 0.;
  LorentzVector*  jet1_ = 0;
  float dPhiDiLepJet1_ = 0.;
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
  int dstype_ = 0;
  int processId_ = 0;
  float higgsPt_ = 0.0;
  float sfWeightHPt_ = 1.;
  float pmet_ = 0.0;
  float met_ = 0.0;
  float metPhi_ = 0.0;
  float mt_ = 0.;
  float dPhiDiLepMET_ = 0.;
  float dPhiLep1MET_;
  float dPhiLep2MET_;
  
  ch->SetBranchAddress( "njets"     , &njets_     );     
  ch->SetBranchAddress( "cuts"      , &cuts_     );     
  ch->SetBranchAddress( "lep1"      , &lep1_      );   
  ch->SetBranchAddress( "lep2"      , &lep2_      );   
  ch->SetBranchAddress( "dilep"      , &dilep_      );   
  ch->SetBranchAddress( "type"      , &type_     );     
  ch->SetBranchAddress( "pTrackMet"      , &pTrackMet_     );   
  ch->SetBranchAddress( "run"     , &run_     );     
  ch->SetBranchAddress( "lumi"     , &lumi_     );     
  ch->SetBranchAddress( "dPhi"      , &dPhi_     );     
  ch->SetBranchAddress( "jet1"      , &jet1_      );   
  ch->SetBranchAddress( "jet2"      , &jet2_      );   
  ch->SetBranchAddress( "jet3"      , &jet3_      );   
  ch->SetBranchAddress( "dPhiDiLepJet1"      , &dPhiDiLepJet1_     );     
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
  ch->SetBranchAddress( "dstype"      , &dstype_     );     
  ch->SetBranchAddress( "processId"      , &processId_     );     
  ch->SetBranchAddress( "higgsPt"      , &higgsPt_     );     
  ch->SetBranchAddress( "sfWeightHPt"      , &sfWeightHPt_     );     

  ch->SetBranchAddress( "pmet"      , &pmet_     );     
  ch->SetBranchAddress( "mt"      , &mt_     );     
  ch->SetBranchAddress( "met"      , &met_     );     
  ch->SetBranchAddress( "metPhi"      , &metPhi_     );     
  ch->SetBranchAddress( "dPhiDiLepMET"      , &dPhiDiLepMET_     );     
  ch->SetBranchAddress( "dPhiLep1MET"      , &dPhiLep1MET_     );     
  ch->SetBranchAddress( "dPhiLep2MET"      , &dPhiLep2MET_     );     

  //==========================================
  // Loop All Events
  //==========================================  
  cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; \n";

  float dymva_= -999.;
  float recoil_= -999.;
  float dPhiJet1MET_= -999.;
  //if branches have to be added:
  //evt_tree->Branch("dymva", &dymva_, "dymva/F");
  //evt_tree->Branch("recoil", &recoil_, "recoil/F");
  //evt_tree->Branch("dPhiJet1MET", &dPhiJet1MET_, "dPhiJet1MET/F");
  //if branches already present:
  ch->SetBranchAddress( "dymva"      , &dymva_     );
  ch->SetBranchAddress( "recoil"      , &recoil_     );
  ch->SetBranchAddress( "dPhiJet1MET"      , &dPhiJet1MET_     );

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++){
    ch->GetEntry(ievt); 

    //from JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py
    float metx = met_ * cos(metPhi_);
    float mety = met_ * sin(metPhi_);
    if (dstype_==0) {
      //2012runAvsNvtx_data
      //px = cms.string("+3.54233e-01 + 2.65299e-01*Nvtx"),
      //py = cms.string("+1.88923e-01 - 1.66425e-01*Nvtx")
      metx -= (+3.54233e-01 + 2.65299e-01*nvtx_);
      mety -= (+1.88923e-01 - 1.66425e-01*nvtx_);
    } else {
      //2012runAvsNvtx_mc
      //px = cms.string("-2.99576e-02 - 6.61932e-02*Nvtx"),
      //py = cms.string("+3.70819e-01 - 1.48617e-01*Nvtx")
      metx -= (-2.99576e-02 - 6.61932e-02*nvtx_);
      mety -= (+3.70819e-01 - 1.48617e-01*nvtx_);
    }
    met_ = sqrt( metx*metx + mety*mety );
    metPhi_ = atan2(mety,metx);
    pmet_ = projectedMet(lep1_->phi(), lep2_->phi(), met_, metPhi_);
    mt_ = mt(dilep_->pt(),met_,deltaPhi(dilep_->phi(),metPhi_));
    LorentzVector metlv( met_*cos(metPhi_), met_*sin(metPhi_), 0, met_ );
    assert((metlv.phi()-metPhi_)<0.001);
    assert((metlv.pt()-met_)<0.001);
    recoil_ = recoilvar(met_, metPhi_, dilep_);
    dPhiJet1MET_  = fabs(ROOT::Math::VectorUtil::DeltaPhi(*jet1_,metlv));
    dPhiDiLepMET_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(*dilep_,metlv));
    dPhiLep1MET_  = fabs(ROOT::Math::VectorUtil::DeltaPhi(*lep1_,metlv));
    dPhiLep2MET_  = fabs(ROOT::Math::VectorUtil::DeltaPhi(*lep2_,metlv));
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
  addDYMVA(smurfFDir,"data.root",outputDir);
  addDYMVA(smurfFDir,"dyll.root",outputDir);
  addDYMVA(smurfFDir,"ggww.root",outputDir);
  addDYMVA(smurfFDir,"hww110.root",outputDir);
  addDYMVA(smurfFDir,"hww115.root",outputDir);
  addDYMVA(smurfFDir,"hww120.root",outputDir);
  addDYMVA(smurfFDir,"hww125.root",outputDir);
  addDYMVA(smurfFDir,"hww130.root",outputDir);
  addDYMVA(smurfFDir,"hww135.root",outputDir);
  addDYMVA(smurfFDir,"hww140.root",outputDir);
  addDYMVA(smurfFDir,"hww145.root",outputDir);
  addDYMVA(smurfFDir,"hww150.root",outputDir);
  addDYMVA(smurfFDir,"hww160.root",outputDir);
  addDYMVA(smurfFDir,"hww170.root",outputDir);
  addDYMVA(smurfFDir,"hww180.root",outputDir);
  addDYMVA(smurfFDir,"hww190.root",outputDir);
  addDYMVA(smurfFDir,"hww200.root",outputDir);
  addDYMVA(smurfFDir,"hww250.root",outputDir);
  addDYMVA(smurfFDir,"hww300.root",outputDir);
  addDYMVA(smurfFDir,"hww350.root",outputDir);
  addDYMVA(smurfFDir,"hww400.root",outputDir);
  addDYMVA(smurfFDir,"hww450.root",outputDir);
  addDYMVA(smurfFDir,"hww500.root",outputDir);
  addDYMVA(smurfFDir,"hww550.root",outputDir);
  addDYMVA(smurfFDir,"hww600.root",outputDir);
  addDYMVA(smurfFDir,"qqww_py.root",outputDir);
  addDYMVA(smurfFDir,"qqww.root",outputDir);
  addDYMVA(smurfFDir,"ttbar.root",outputDir);
  addDYMVA(smurfFDir,"tw.root",outputDir);
  addDYMVA(smurfFDir,"wjets.root",outputDir);
  addDYMVA(smurfFDir,"wz.root",outputDir);
  addDYMVA(smurfFDir,"zz.root",outputDir);
  addDYMVA(smurfFDir,"wgamma.root",outputDir);
  addDYMVA(smurfFDir,"wglll.root",outputDir);
  addDYMVA(smurfFDir,"ttbar_powheg.root",outputDir);
  addDYMVA(smurfFDir,"wwmcnlodown.root",outputDir);
  addDYMVA(smurfFDir,"wwmcnlo.root",outputDir);
  addDYMVA(smurfFDir,"wwmcnloup.root",outputDir);
  return;
}
