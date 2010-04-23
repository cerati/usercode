#include <TVector.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
using namespace std;
void test2 () {
  //gROOT->Reset();
  //gROOT->SetStyle("Plain");
  //gStyle->SetOptFit(1);
  
  
  TFile * f1 = TFile::Open("ntuple_Zanalisys7TeV_data7TeV.root");
  TTree *t1 = (TTree*)f1->Get("ntuple");

  int run, event, lumi;
  t1->SetBranchAddress("run",&run);
  t1->SetBranchAddress("event",&event);
  t1->SetBranchAddress("lumi",&lumi);

  string listhlt[140] = {"HLT_Activity_L1A","HLT_Activity_PixelClusters","HLT_Activity_DT","HLT_Activity_DT_Tuned","HLT_Activity_Ecal","HLT_Activity_EcalREM","HLT_SelectEcalSpikes_L1R","HLT_SelectEcalSpikesHighEt_L1R","HLT_L1Jet6U","HLT_L1Jet6U_NoBPTX","HLT_L1Jet10U","HLT_L1Jet10U_NoBPTX","HLT_Jet15U","HLT_Jet30U","HLT_Jet50U","HLT_L1SingleForJet","HLT_L1SingleForJet_NoBPTX","HLT_L1SingleCenJet","HLT_L1SingleCenJet_NoBPTX","HLT_L1SingleTauJet","HLT_L1SingleTauJet_NoBPTX","HLT_FwdJet20U","HLT_DiJetAve15U_8E29","HLT_DiJetAve30U_8E29","HLT_DoubleJet15U_ForwardBackward","HLT_QuadJet15U","HLT_L1MET20","HLT_MET45","HLT_MET100","HLT_HT100U","HLT_L1MuOpen","HLT_L1MuOpen_NoBPTX","HLT_L1Mu","HLT_L1Mu20","HLT_L2Mu9","HLT_L2Mu11","HLT_IsoMu3","HLT_Mu3","HLT_Mu5","HLT_Mu9","HLT_L1DoubleMuOpen","HLT_DoubleMu0","HLT_DoubleMu3","HLT_Mu0_L1MuOpen","HLT_Mu3_L1MuOpen","HLT_Mu5_L1MuOpen","HLT_Mu0_Track0_Jpsi","HLT_Mu3_Track0_Jpsi","HLT_Mu5_Track0_Jpsi","HLT_L1SingleEG1","HLT_L1SingleEG1_NoBPTX","HLT_L1SingleEG2","HLT_L1SingleEG2_NoBPTX","HLT_L1SingleEG5","HLT_L1SingleEG5_NoBPTX","HLT_L1SingleEG8","HLT_L1SingleEG20_NoBPTX","HLT_L1DoubleEG5","HLT_EgammaSuperClusterOnly_L1R","HLT_Ele10_LW_L1R","HLT_Ele10_LW_EleId_L1R","HLT_Ele15_LW_L1R","HLT_Ele15_SC10_LW_L1R","HLT_Ele15_SiStrip_L1R","HLT_Ele20_LW_L1R","HLT_DoubleEle5_SW_L1R","HLT_Photon10_L1R","HLT_Photon15_L1R","HLT_Photon15_TrackIso_L1R","HLT_Photon15_LooseEcalIso_L1R","HLT_Photon20_L1R","HLT_Photon30_L1R_8E29","HLT_DoublePhoton5_eeRes_L1R","HLT_DoublePhoton5_Jpsi_L1R","HLT_DoublePhoton5_Upsilon_L1R","HLT_DoublePhoton10_L1R","HLT_SingleLooseIsoTau20","HLT_DoubleLooseIsoTau15","HLT_BTagIP_Jet50U","HLT_BTagMu_Jet10U","HLT_StoppedHSCP_8E29","HLT_L1Mu14_L1SingleEG10","HLT_L1Mu14_L1SingleJet6U","HLT_L1Mu14_L1ETM30","HLT_ZeroBias","HLT_MinBias","HLT_MinBiasBSC","HLT_MinBiasBSC_NoBPTX","HLT_MinBiasBSC_OR","HLT_MinBiasHcal","HLT_MinBiasEcal","HLT_ZeroBiasPixel_SingleTrack","HLT_MinBiasPixel_SingleTrack","HLT_MinBiasPixel_DoubleTrack","HLT_MinBiasPixel_DoubleIsoTrack5","HLT_CSCBeamHalo","HLT_CSCBeamHaloOverlapRing1","HLT_CSCBeamHaloOverlapRing2","HLT_CSCBeamHaloRing2or3","HLT_BackwardBSC","HLT_ForwardBSC","HLT_HighMultiplicityBSC","HLT_SplashBSC","HLT_L1_BSC","HLT_L1_BscMinBiasOR_BptxPlusORMinus","HLT_L1_BscMinBiasOR_BptxPlusORMinus_NoBPTX","HLT_RPCBarrelCosmics","HLT_TrackerCosmics","HLT_L1Tech_RPC_TTU_RBst1_collisions","HLT_IsoTrackHE_8E29","HLT_IsoTrackHB_8E29","HLT_HcalPhiSym","HLT_HcalNZS_8E29","AlCa_EcalPhiSym","AlCa_EcalPi0_8E29","AlCa_EcalEta_8E29","AlCa_RPCMuonNoHits","AlCa_RPCMuonNoTriggers","AlCa_RPCMuonNormalisation","HLT_DTErrors","HLT_HighMult40","HLT_Calibration","HLT_EcalCalibration","HLT_HcalCalibration","HLT_Random","HLT_L1_HFtech","HLT_L1Tech_HCAL_HF_coincidence_PM","HLT_HFThreshold3","HLT_HFThreshold10","HLT_GlobalRunHPDNoise","HLT_TechTrigHCALNoise","HLT_L1_BPTX","HLT_L1_BPTX_MinusOnly","HLT_L1_BPTX_PlusOnly","HLT_L2Mu0_NoVertex","HLT_TkMu3_NoVertex","HLT_PhysicsDeclared","HLT_LogMonitor","DQM_FEDIntegrity","HLTriggerFinalPath"};
  std::vector<string> listHLT(listhlt,listhlt+140);
  std::vector<bool> *trigBits = 0;
  t1->SetBranchAddress("trigBits",&trigBits);

  std::vector<bool> *vtxFake=0;
  std::vector<float> *vtxX=0,*vtxY=0,*vtxZ=0,*vtxXE=0,*vtxYE=0,*vtxZE=0;
  std::vector<int> *vtxNdof=0,*vtxNtks=0;
  t1->SetBranchAddress("vtxFake",&vtxFake);
  t1->SetBranchAddress("vtxX",&vtxX);
  t1->SetBranchAddress("vtxY",&vtxY);
  t1->SetBranchAddress("vtxZ",&vtxZ);
  t1->SetBranchAddress("vtxXE",&vtxXE);
  t1->SetBranchAddress("vtxYE",&vtxYE);
  t1->SetBranchAddress("vtxZE",&vtxZE);
  t1->SetBranchAddress("vtxNdof",&vtxNdof);
  t1->SetBranchAddress("vtxNtks",&vtxNtks);

  std::vector<bool> *muGLB=0,*muTRK=0;
  std::vector<float> *muPt=0, *muPtE=0, *muEta=0, *muPhi=0, *muDxy=0, *muDxyE=0, *muDz=0, *muDzE=0, 
    *muIsoTRK03=0, *muIsoECAL03=0, *muIsoHCAL03=0, *muIsoHO03=0, *muIsoTRK05=0, *muIsoECAL05=0, *muIsoHCAL05=0, *muIsoHO05=0;
  std::vector<int> *muQ=0;
  t1->SetBranchAddress("muGLB",&muGLB);
  t1->SetBranchAddress("muTRK",&muTRK);
  t1->SetBranchAddress("muPt",&muPt);
  t1->SetBranchAddress("muPtE",&muPtE);
  t1->SetBranchAddress("muEta",&muEta);
  t1->SetBranchAddress("muPhi",&muPhi);
  t1->SetBranchAddress("muDxy",&muDxy);
  t1->SetBranchAddress("muDxyE",&muDxyE);
  t1->SetBranchAddress("muDz",&muDz);
  t1->SetBranchAddress("muDzE",&muDzE);
  t1->SetBranchAddress("muIsoTRK03",&muIsoTRK03);
  t1->SetBranchAddress("muIsoECAL03",&muIsoECAL03);
  t1->SetBranchAddress("muIsoHCAL03",&muIsoHCAL03);
  t1->SetBranchAddress("muIsoHO03",&muIsoHO03);
  t1->SetBranchAddress("muIsoTRK05",&muIsoTRK05);
  t1->SetBranchAddress("muIsoECAL05",&muIsoECAL05);
  t1->SetBranchAddress("muIsoHCAL05",&muIsoHCAL05);
  t1->SetBranchAddress("muIsoHO05",&muIsoHO05);
  t1->SetBranchAddress("muQ",&muQ);

  std::vector<bool> *elIdLoose=0,*elIdTight=0, *elIdRobLo=0, *elIdRobTi=0, *elIdRobHE=0;
  std::vector<float> *elESC=0, *elEoP=0, *elFbrem=0, *elPt=0, *elPtE=0, *elEta=0, *elPhi=0, *elDxy=0, *elDxyE=0, *elDz=0, *elDzE=0, 
    *elIsoTRK03=0, *elIsoECAL03=0, *elIsoHCAL03=0, *elIsoTRK04=0, *elIsoECAL04=0, *elIsoHCAL04=0;
  std::vector<int> *elQ=0, *elClass=0;
  t1->SetBranchAddress("elIdLoose",&elIdLoose);
  t1->SetBranchAddress("elIdTight",&elIdTight);
  t1->SetBranchAddress("elIdRobLo",&elIdRobLo);
  t1->SetBranchAddress("elIdRobTi",&elIdRobTi);
  t1->SetBranchAddress("elIdRobHE",&elIdRobHE);
  t1->SetBranchAddress("elESC",&elESC);
  t1->SetBranchAddress("elEoP",&elEoP);
  t1->SetBranchAddress("elFbrem",&elFbrem);
  t1->SetBranchAddress("elPt",&elPt);
  t1->SetBranchAddress("elPtE",&elPtE);
  t1->SetBranchAddress("elEta",&elEta);
  t1->SetBranchAddress("elPhi",&elPhi);
  t1->SetBranchAddress("elDxy",&elDxy);
  t1->SetBranchAddress("elDxyE",&elDxyE);
  t1->SetBranchAddress("elDz",&elDz);
  t1->SetBranchAddress("elDzE",&elDzE);
  t1->SetBranchAddress("elIsoTRK03",&elIsoTRK03);
  t1->SetBranchAddress("elIsoECAL03",&elIsoECAL03);
  t1->SetBranchAddress("elIsoHCAL03",&elIsoHCAL03);
  t1->SetBranchAddress("elIsoTRK04",&elIsoTRK04);
  t1->SetBranchAddress("elIsoECAL04",&elIsoECAL04);
  t1->SetBranchAddress("elIsoHCAL04",&elIsoHCAL04);
  t1->SetBranchAddress("elQ",&elQ);
  t1->SetBranchAddress("elClass",&elClass);

  float metEt, metPhi, metEtSig, metSig, metUncorr, metGenPt, metGenPhi;
  t1->SetBranchAddress("metEt",&metEt);
  t1->SetBranchAddress("metPhi",&metPhi);
  t1->SetBranchAddress("metEtSig",&metEtSig);
  t1->SetBranchAddress("metSig",&metSig);
  t1->SetBranchAddress("metUncorr",&metUncorr);
  t1->SetBranchAddress("metGenPt",&metGenPt);
  t1->SetBranchAddress("metGenPhi",&metGenPhi);

  float pfMetEt, pfMetPhi, pfMetEtSig, pfMetSig, pfMetUncorr, pfMetGenPt, pfMetGenPhi;
  t1->SetBranchAddress("pfMetEt",&pfMetEt);
  t1->SetBranchAddress("pfMetPhi",&pfMetPhi);
  t1->SetBranchAddress("pfMetEtSig",&pfMetEtSig);
  t1->SetBranchAddress("pfMetSig",&pfMetSig);
  t1->SetBranchAddress("pfMetUncorr",&pfMetUncorr);
  t1->SetBranchAddress("pfMetGenPt",&pfMetGenPt);
  t1->SetBranchAddress("pfMetGenPhi",&pfMetGenPhi);

  std::vector<float> *jetPt=0, *jetEta=0, *jetPhi=0, *jetCorrFact=0, 
    *jetDiscTCHE=0,*jetDiscTCHP=0,*jetDiscSFTM=0,*jetDiscSMTMPT=0,*jetDiscSFTMIP3=0,
    *jetDiscSFTE=0,*jetDiscJETP=0,*jetDiscJETBP=0,*jetDiscSSV=0,*jetDiscCSV=0,*jetDiscCSVMVA=0,*jetQ=0;
  std::vector<int> *jetGenId=0, *jetGenPt=0, *jetTks=0;
  t1->SetBranchAddress("jetPt",&jetPt);
  t1->SetBranchAddress("jetEta",&jetEta);
  t1->SetBranchAddress("jetPhi",&jetPhi);
  t1->SetBranchAddress("jetCorrFact",&jetCorrFact);
  t1->SetBranchAddress("jetDiscTCHE",&jetDiscTCHE);
  t1->SetBranchAddress("jetDiscTCHP",&jetDiscTCHP);
  t1->SetBranchAddress("jetDiscSFTM",&jetDiscSFTM);
  t1->SetBranchAddress("jetDiscSMTMPT",&jetDiscSMTMPT);
  t1->SetBranchAddress("jetDiscSFTMIP3",&jetDiscSFTMIP3);
  t1->SetBranchAddress("jetDiscSFTE",&jetDiscSFTE);
  t1->SetBranchAddress("jetDiscJETP",&jetDiscJETP);
  t1->SetBranchAddress("jetDiscJETBP",&jetDiscJETBP);
  t1->SetBranchAddress("jetDiscSSV",&jetDiscSSV);
  t1->SetBranchAddress("jetDiscCSV",&jetDiscCSV);
  t1->SetBranchAddress("jetDiscCSVMVA",&jetDiscCSVMVA);
  t1->SetBranchAddress("jetGenId",&jetGenId);
  t1->SetBranchAddress("jetGenPt",&jetGenPt);
  t1->SetBranchAddress("jetQ",&jetQ);
  t1->SetBranchAddress("jetTks",&jetTks);

  int n = (int) t1->GetEntries();
  cout << " Tree Entries " << n << endl; 
	
  for (int evt = 0; evt < n; ++evt) { //event loop
    //for (int evt = 25000; evt < 30000; ++evt) { //event loop
    t1->GetEntry(evt); //fill the variables with the i-th event

    cout << "run=" << run << " evt=" << event << " lumi=" << lumi << endl;
    //cout << trigBits->size() << endl;
    cout << "triggers passed= ";
    for (unsigned int i=0;i!=trigBits->size();++i) if (trigBits->at(i)) cout << listHLT[i] << ", ";
    cout << endl;

    cout << "N vertices=" << vtxFake->size() << endl;
    for (unsigned int i=0;i!=vtxFake->size();++i) {
      cout << "vtx fake=" << vtxFake->at(i) 
	   << " x=" << vtxX->at(i) << " y=" << vtxY->at(i) << " z=" << vtxZ->at(i)
	   << " xE=" << vtxXE->at(i) << " yE=" << vtxYE->at(i) << " zE=" << vtxZE->at(i)
	   << " ndof=" << vtxNdof->at(i) << " ntks=" << vtxNtks->at(i)
	   << endl;
    }

    cout << "N muons=" << muGLB->size() << endl;
    for (unsigned int i=0;i!=muGLB->size();++i) {
      cout << "mu global=" << muGLB->at(i) << " tracker=" << muTRK->at(i)
	   << " pt=" << muPt->at(i)
	   << " ptE=" << muPtE->at(i)
	   << " eta=" << muEta->at(i)
	   << " phi=" << muPhi->at(i)
	   << " dxy=" << muDxy->at(i)
	   << " dxyE=" << muDxyE->at(i)
	   << " dz=" << muDz->at(i)
	   << " dzE=" << muDzE->at(i)
	   << " isoTrk03=" << muIsoTRK03->at(i)
	   << " isoEcal03=" << muIsoECAL03->at(i)
	   << " isoHcal03=" << muIsoHCAL03->at(i)
	   << " isoHo03=" << muIsoHO03->at(i)
	   << " isoTrk05=" << muIsoTRK05->at(i)
	   << " isoEcal05=" << muIsoECAL05->at(i)
	   << " isoHcal05=" << muIsoHCAL05->at(i)
	   << " isoHo05=" << muIsoHO05->at(i)
	   << " q=" << muQ->at(i)
	   << endl;
    }

    cout << "N electrons=" << elIdLoose->size() << endl;
    for (unsigned int i=0;i!=elIdLoose->size();++i) {
      cout << "el idLoose =" << elIdLoose->at(i)
	   << " idTight=" << elIdTight->at(i)
	   << " idRobLo=" << elIdRobLo->at(i)
	   << " idRobTi=" << elIdRobTi->at(i)
	   << " idRobHE=" << elIdRobHE->at(i)
	   << " eSC=" << elESC->at(i)
	   << " e/p=" << elEoP->at(i)
	   << " fbrem=" << elFbrem->at(i)
	   << " pt=" << elPt->at(i)
	   << " ptE=" << elPtE->at(i)
	   << " eta=" << elEta->at(i)
	   << " phi=" << elPhi->at(i)
	   << " dxy=" << elDxy->at(i)
	   << " dxyE=" << elDxyE->at(i)
	   << " dz=" << elDz->at(i)
	   << " dzE=" << elDzE->at(i)
	   << " isoTrk03=" << elIsoTRK03->at(i)
	   << " isoEcal03=" << elIsoECAL03->at(i)
	   << " isoHcal03=" << elIsoHCAL03->at(i)
	   << " isoTrk04=" << elIsoTRK04->at(i)
	   << " isoEcal04=" << elIsoECAL04->at(i)
	   << " isoHcal04=" << elIsoHCAL04->at(i)
	   << " q=" << elQ->at(i)
	   << " class=" << elClass->at(i)
	   << endl;
    }

    cout << "met et=" << metEt << " phi=" << metPhi << " etSig=" << metEtSig << " sig=" << metSig 
	 << " uncorr=" << metUncorr << " genPt=" << metGenPt << " genPhi=" << metGenPhi << endl;;

    cout << "pfMet et=" << pfMetEt << " phi=" << pfMetPhi << " etSig=" << pfMetEtSig << " sig=" << pfMetSig 
	 << " uncorr=" << pfMetUncorr << " genPt=" << pfMetGenPt << " genPhi=" << pfMetGenPhi << endl;;

    cout << "N jets=" << jetPt->size() << endl;
    for (unsigned int i=0;i!=jetPt->size();++i) {
      cout << "jet pt=" << jetPt->at(i)
	   << " eta=" << jetEta->at(i)
	   << " phi=" << jetPhi->at(i)
	   << " corrFact=" << jetCorrFact->at(i)
	   << " disc_tcHE=" << jetDiscTCHE->at(i)
	   << " disc_tcHP=" << jetDiscTCHP->at(i)
	   << " disc_sftM=" << jetDiscSFTM->at(i)
	   << " disc_sftMPt=" << jetDiscSMTMPT->at(i)
	   << " disc_sftMIP3=" << jetDiscSFTMIP3->at(i)
	   << " disc_sftE=" << jetDiscSFTE->at(i)
	   << " disc_jetP=" << jetDiscJETP->at(i)
	   << " disc_jetBP=" << jetDiscJETBP->at(i)
	   << " disc_ssv=" << jetDiscSSV->at(i)
	   << " disc_csv==" << jetDiscCSV->at(i)
	   << " disc_cvsMVA=" << jetDiscCSVMVA->at(i)
	   << " genId=" << jetGenId->at(i)
	   << " genPt=" << jetGenPt->at(i)
	   << " charge=" << jetQ->at(i)
	   << " nTks=" << jetTks->at(i)
	   <<endl;
    }

  } // end event loop


  f1->Close();
}
