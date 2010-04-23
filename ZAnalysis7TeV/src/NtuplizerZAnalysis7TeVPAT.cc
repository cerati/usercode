// -*- C++ -*-
//
// Package:    ZAnalysis7TeV
// Class:      NtuplizerZAnalysis7TeVPAT
// 
/**\class NtuplizerZAnalysis7TeVPAT NtuplizerZAnalysis7TeVPAT.cc Tests/ZAnalysis7TeV/src/NtuplizerZAnalysis7TeVPAT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giuseppe Cerati,28 S-012,+41227678302,
//         Created:  Tue Apr 20 12:53:58 CEST 2010
// $Id$
//
//

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>

typedef struct { 
  int run, event, lumi;

  std::vector<bool> trigBits;

  std::vector<bool> vtxFake;
  std::vector<float> vtxX,vtxY,vtxZ,vtxXE,vtxYE,vtxZE;
  std::vector<int> vtxNdof,vtxNtks;

  std::vector<bool> muGLB,muTRK;
  std::vector<float> muPt, muPtE, muEta, muPhi, muDxy, muDxyE, muDz, muDzE, 
    muIsoTRK03, muIsoECAL03, muIsoHCAL03, muIsoHO03, muIsoTRK05, muIsoECAL05, muIsoHCAL05, muIsoHO05;
  std::vector<int> muQ;

  std::vector<bool> elIdLoose,elIdTight, elIdRobLo, elIdRobTi, elIdRobHE;
  std::vector<float> elESC, elEoP, elFbrem, elPt, elPtE, elEta, elPhi, elDxy, elDxyE, elDz, elDzE, 
    elIsoTRK03, elIsoECAL03, elIsoHCAL03, elIsoTRK04, elIsoECAL04, elIsoHCAL04;
  std::vector<int> elQ, elClass;

  float metEt, metPhi, metEtSig, metSig, metUncorr, metGenPt, metGenPhi;
  float pfMetEt, pfMetPhi, pfMetEtSig, pfMetSig, pfMetUncorr, pfMetGenPt, pfMetGenPhi;

  std::vector<float> jetPt, jetEta, jetPhi, jetCorrFact, 
    jetDiscTCHE,jetDiscTCHP,jetDiscSFTM,jetDiscSMTMPT,jetDiscSFTMIP3,
    jetDiscSFTE,jetDiscJETP,jetDiscJETBP,jetDiscSSV,jetDiscCSV,jetDiscCSVMVA, jetQ;
  std::vector<int> jetGenId, jetGenPt, jetTks;

} BRANCH;

class NtuplizerZAnalysis7TeVPAT : public edm::EDAnalyzer {
public:
  explicit NtuplizerZAnalysis7TeVPAT(const edm::ParameterSet&);
  ~NtuplizerZAnalysis7TeVPAT();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  TFile * file;
  std::string outfile;
  BRANCH branch;
  TTree *ntuple;
};

NtuplizerZAnalysis7TeVPAT::NtuplizerZAnalysis7TeVPAT(const edm::ParameterSet& iConfig):
  outfile(iConfig.getParameter<std::string>("outfile")),
  ntuple(0)
{}

NtuplizerZAnalysis7TeVPAT::~NtuplizerZAnalysis7TeVPAT(){}


void NtuplizerZAnalysis7TeVPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace pat;
  using namespace std;

  branch.run=iEvent.run();
  branch.event=iEvent.id().event();
  branch.lumi=iEvent.luminosityBlock();
  cout << "event=" << iEvent.id().event() << " run=" << iEvent.run() << " lumi=" << iEvent.luminosityBlock() << endl;

  std::vector<bool> trigBits;
  edm::Handle<edm::TriggerResults> hltresults;
  InputTag tag("TriggerResults","","HLT");//InputTag tag("TriggerResults");
  iEvent.getByLabel(tag,hltresults);
  edm::TriggerNames triggerNames_ = iEvent.triggerNames(*hltresults);
  int ntrigs=hltresults->size();
  vector<string> triggernames = triggerNames_.triggerNames();
  //int trigjolly = 0;
  for (int itrig = 0; itrig != ntrigs; ++itrig){
     string trigName=triggerNames_.triggerName(itrig);
     bool accept = hltresults->accept(itrig);
     cout << trigName << " " << accept << endl;
     trigBits.push_back(accept);
     //      if (accept) {
     //        if (trigName=="HLT1MuonIso") trigjolly+=1;
     //        else if (trigName=="HLT1MuonNonIso") trigjolly+=10;
     //        else if (trigName=="HLT1Electron") trigjolly+=100;
     //        else if (trigName=="HLT1ElectronRelaxed") trigjolly+=1000;
     //     }
  }
  branch.trigBits=trigBits;


  std::vector<bool> vtxFake;
  std::vector<float> vtxX,vtxY,vtxZ,vtxXE,vtxYE,vtxZE;
  std::vector<int> vtxNdof,vtxNtks;
  Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel("offlinePrimaryVertices", vertexHandle);
  cout << "nVertices=" << vertexHandle->size() << endl;
  for (reco::VertexCollection::const_iterator vtx=vertexHandle->begin();vtx!=vertexHandle->end();++vtx){
    cout << "vtx isFake=" << vtx->isFake() <<" pos=" << vtx->position() << " ndof=" << vtx->ndof() << " nTks=" << vtx->tracksSize() << endl; 
    vtxFake.push_back(vtx->isFake());
    vtxX.push_back(vtx->position().x());
    vtxY.push_back(vtx->position().y());
    vtxZ.push_back(vtx->position().z());
    vtxXE.push_back(vtx->xError());
    vtxYE.push_back(vtx->yError());
    vtxZE.push_back(vtx->zError());
    vtxNdof.push_back(vtx->ndof());
    vtxNtks.push_back(vtx->tracksSize());
  }
  branch.vtxFake=vtxFake;
  branch.vtxX=vtxX;
  branch.vtxY=vtxY;
  branch.vtxZ=vtxZ;
  branch.vtxXE=vtxXE;
  branch.vtxYE=vtxYE;
  branch.vtxZE=vtxZE;
  branch.vtxNdof=vtxNdof;
  branch.vtxNtks=vtxNtks;

  std::vector<bool> muGLB,muTRK;
  std::vector<float> muPt, muPtE, muEta, muPhi, muDxy, muDxyE, muDz, muDzE, 
    muIsoTRK03, muIsoECAL03, muIsoHCAL03, muIsoHO03, muIsoTRK05, muIsoECAL05, muIsoHCAL05, muIsoHO05;
  std::vector<int> muQ;
  Handle<MuonCollection> muonHandle;
  iEvent.getByLabel("cleanPatMuons", muonHandle);
  cout << "nMuons=" << muonHandle->size() << endl;
  for (MuonCollection::const_iterator mu=muonHandle->begin();mu!=muonHandle->end();++mu){
    cout << "mu global=" << mu->isGlobalMuon() << " tracker=" << mu->isTrackerMuon() 
	 << " pt=" << mu->pt() << " eta=" << mu->eta() << " phi=" << mu->phi() << " q=" << mu->charge()
	 << " iso03.sumPt=" << mu->isolationR03().sumPt << " iso03.emEt=" << mu->isolationR03().emEt 
	 << " iso03.hadEt=" << mu->isolationR03().hadEt << " iso03.hoEt=" << mu->isolationR03().hoEt 
	 << " iso05.sumPt=" << mu->isolationR05().sumPt << " iso05.emEt=" << mu->isolationR05().emEt 
	 << " iso05.hadEt=" << mu->isolationR05().hadEt << " iso05.hoEt=" << mu->isolationR05().hoEt 
	 << endl;
    muGLB.push_back(mu->isGlobalMuon());
    muTRK.push_back(mu->isTrackerMuon());
    muPt.push_back(mu->pt());
    muEta.push_back(mu->eta());
    muPhi.push_back(mu->phi());
    muIsoTRK03.push_back(mu->isolationR03().sumPt);
    muIsoECAL03.push_back(mu->isolationR03().emEt);
    muIsoHCAL03.push_back(mu->isolationR03().hadEt);
    muIsoHO03.push_back(mu->isolationR03().hoEt);
    muIsoTRK05.push_back(mu->isolationR05().sumPt);
    muIsoECAL05.push_back(mu->isolationR05().emEt);
    muIsoHCAL05.push_back(mu->isolationR05().hadEt);
    muIsoHO05.push_back(mu->isolationR05().hoEt);
    muQ.push_back(mu->charge());
    if (mu->track().isNonnull()){
      muPtE.push_back(mu->track()->ptError());
      muDxy.push_back(mu->track()->dxy(vertexHandle->begin()->position()));
      muDxyE.push_back(mu->track()->dxyError());
      muDz.push_back(mu->track()->dz(vertexHandle->begin()->position()));
      muDzE.push_back(mu->track()->dzError());
    } else {
      muPtE.push_back(mu->standAloneMuon()->ptError());
      muDxy.push_back(mu->standAloneMuon()->dxy(vertexHandle->begin()->position()));
      muDxyE.push_back(mu->standAloneMuon()->dxyError());
      muDz.push_back(mu->standAloneMuon()->dz(vertexHandle->begin()->position()));
      muDzE.push_back(mu->standAloneMuon()->dzError());
    }
  }
  branch.muGLB=muGLB;
  branch.muTRK=muTRK;
  branch.muPt=muPt;
  branch.muPtE=muPtE;
  branch.muEta=muEta;
  branch.muPhi=muPhi;
  branch.muDxy=muDxy;
  branch.muDxyE=muDxyE;
  branch.muDz=muDz;
  branch.muDzE=muDzE;
  branch.muIsoTRK03=muIsoTRK03;
  branch.muIsoECAL03=muIsoECAL03;
  branch.muIsoHCAL03=muIsoHCAL03;
  branch.muIsoHO03=muIsoHO03;
  branch.muIsoTRK05=muIsoTRK05;
  branch.muIsoECAL05=muIsoECAL05;
  branch.muIsoHCAL05=muIsoHCAL05;
  branch.muIsoHO05=muIsoHO05;
  branch.muQ=muQ;

  std::vector<bool> elIdLoose,elIdTight, elIdRobLo, elIdRobTi, elIdRobHE;
  std::vector<float> elESC, elEoP, elFbrem, elPt, elPtE, elEta, elPhi, elDxy, elDxyE, elDz, elDzE, 
    elIsoTRK03, elIsoECAL03, elIsoHCAL03, elIsoTRK04, elIsoECAL04, elIsoHCAL04;
  std::vector<int> elQ, elClass;
  Handle<ElectronCollection> electronHandle;
  iEvent.getByLabel("cleanPatElectrons", electronHandle);
  cout << "nElectrons=" << electronHandle->size() << endl;
  for (ElectronCollection::const_iterator el=electronHandle->begin();el!=electronHandle->end();++el) {
    cout << "el e=" << el->superCluster()->energy() << " e/p=" << el->eSuperClusterOverP() 
	 << " fbrem=" << el->fbrem() << " class=" << el->classification()
	 << " pt=" << el->pt() << " eta=" << el->eta() << " phi=" << el->phi() << " q=" << el->charge()
	 << " iso03_tk=" << el->dr03TkSumPt() 
	 << " iso03_ecal=" << el->dr03EcalRecHitSumEt() 
	 << " iso03_hcal=" << el->dr03HcalTowerSumEt()
	 << " iso04_tk=" << el->dr04TkSumPt() 
	 << " iso04_ecal=" << el->dr04EcalRecHitSumEt() 
	 << " iso04_hcal=" << el->dr04HcalTowerSumEt()
	 << " idLoose=" << el->electronID("eidLoose")
	 << " idTight=" << el->electronID("eidTight")
	 << " idRobLo=" << el->electronID("eidRobustLoose")
	 << " idRobTi=" << el->electronID("eidRobustTight")
	 << " idRobHE=" << el->electronID("eidRobustHighEnergy")
	 << endl;
    elIdLoose.push_back(el->electronID("eidLoose"));
    elIdTight.push_back(el->electronID("eidTight"));
    elIdRobLo.push_back(el->electronID("eidRobustLoose"));
    elIdRobTi.push_back(el->electronID("eidRobustTight"));
    elIdRobHE.push_back(el->electronID("eidRobustHighEnergy"));
    elESC.push_back(el->superCluster()->energy());
    elEoP.push_back(el->eSuperClusterOverP());
    elFbrem.push_back(el->fbrem());
    elPt.push_back(el->pt());
    elEta.push_back(el->eta());
    elPhi.push_back(el->phi());
    elIsoTRK03.push_back(el->dr03TkSumPt());
    elIsoECAL03.push_back(el->dr03EcalRecHitSumEt());
    elIsoHCAL03.push_back(el->dr03HcalTowerSumEt());
    elIsoTRK04.push_back(el->dr04TkSumPt());
    elIsoECAL04.push_back(el->dr04EcalRecHitSumEt());
    elIsoHCAL04.push_back(el->dr04HcalTowerSumEt());
    elQ.push_back(el->charge());
    elClass.push_back(el->classification());
    if (el->gsfTrack().isNonnull()){
      elPtE.push_back(el->gsfTrack()->ptError());
      elDxy.push_back(el->gsfTrack()->dxy(vertexHandle->begin()->position()));
      elDxyE.push_back(el->gsfTrack()->dxyError());
      elDz.push_back(el->gsfTrack()->dz(vertexHandle->begin()->position()));
      elDzE.push_back(el->gsfTrack()->dzError());
    } else {
      elPtE.push_back(0);
      elDxy.push_back(0);
      elDxyE.push_back(0);
      elDz.push_back(0);
      elDzE.push_back(0);
    }
  }
  branch.elIdLoose=elIdLoose;
  branch.elIdTight=elIdTight;
  branch.elIdRobLo=elIdRobLo;
  branch.elIdRobTi=elIdRobTi;
  branch.elIdRobHE=elIdRobHE;
  branch.elESC=elESC;
  branch.elEoP=elEoP;
  branch.elFbrem=elFbrem;
  branch.elPt=elPt;
  branch.elPtE=elPtE;
  branch.elEta=elEta;
  branch.elPhi=elPhi;
  branch.elDxy=elDxy;
  branch.elDxyE=elDxyE;
  branch.elDz=elDz;
  branch.elDzE=elDzE;
  branch.elIsoTRK03=elIsoTRK03;
  branch.elIsoECAL03=elIsoECAL03;
  branch.elIsoHCAL03=elIsoHCAL03;
  branch.elIsoTRK04=elIsoTRK04;
  branch.elIsoECAL04=elIsoECAL04;
  branch.elIsoHCAL04=elIsoHCAL04;
  branch.elQ=elQ;
  branch.elClass=elClass;

  Handle<METCollection> metHandle;
  iEvent.getByLabel("patMETs", metHandle);
  cout << "nMETs=" << metHandle->size() << endl;
  float genMetPhi = 0, genMetPt = 0;
  if (metHandle->begin()->genMET()) {
    genMetPhi = metHandle->begin()->genMET()->phi();
    genMetPt = metHandle->begin()->genMET()->pt();
  }
  cout << "met=" << metHandle->begin()->sumEt() << " phi=" << metHandle->begin()->phi() 
       << " EtSig=" << metHandle->begin()->mEtSig() << " sig=" << metHandle->begin()->significance()
       << " uncorr=" << metHandle->begin()->uncorrectedPt()
       << " genPt=" << genMetPt << " genPhi=" << genMetPhi
       << endl;
  branch.metEt=metHandle->begin()->sumEt();
  branch.metPhi=metHandle->begin()->phi();
  branch.metEtSig=metHandle->begin()->mEtSig();
  branch.metSig=metHandle->begin()->significance();
  branch.metUncorr=metHandle->begin()->uncorrectedPt();
  branch.metGenPt=genMetPt;
  branch.metGenPhi=genMetPhi;

  Handle<METCollection> pfmetHandle;
  iEvent.getByLabel("patMETsPF", pfmetHandle);
  cout << "nPFMETs=" << pfmetHandle->size() << endl;
  cout << "met=" << pfmetHandle->begin()->sumEt() << " phi=" << pfmetHandle->begin()->phi() 
       << " EtSig=" << pfmetHandle->begin()->mEtSig() << " sig=" << pfmetHandle->begin()->significance() 
       << " uncorr=" << pfmetHandle->begin()->uncorrectedPt()
       << " genPt=" << genMetPt << " genPhi=" << genMetPhi
       << endl;
  branch.pfMetEt=pfmetHandle->begin()->sumEt();
  branch.pfMetPhi=pfmetHandle->begin()->phi();
  branch.pfMetEtSig=pfmetHandle->begin()->mEtSig();
  branch.pfMetSig=pfmetHandle->begin()->significance();
  branch.pfMetUncorr=pfmetHandle->begin()->uncorrectedPt();
  branch.pfMetGenPt=genMetPt;
  branch.pfMetGenPhi=genMetPhi;

  std::vector<float> jetPt, jetEta, jetPhi, jetCorrFact, 
    jetDiscTCHE,jetDiscTCHP,jetDiscSFTM,jetDiscSMTMPT,jetDiscSFTMIP3,
    jetDiscSFTE,jetDiscJETP,jetDiscJETBP,jetDiscSSV,jetDiscCSV,jetDiscCSVMVA, jetQ;
  std::vector<int> jetGenId, jetGenPt, jetTks;
  Handle<JetCollection> jetHandle;
  iEvent.getByLabel("selectedPatJets", jetHandle);
  cout << "nJets=" << jetHandle->size() << endl;
  for (JetCollection::const_iterator jet=jetHandle->begin();jet!=jetHandle->end();++jet) {
    int id=0;
    if (jet->genParton()) id = jet->genParton()->pdgId();
    float genPt=0;
    if (jet->genJet()) genPt = jet->genJet()->pt();
    cout << "jet pt=" << jet->pt() << " eta=" << jet->eta() << " phi=" << jet->phi()
	 << " corrFact=" << jet->corrFactor("raw")
	 << " disc_tcHE=" << jet->bDiscriminator("trackCountingHighEffBJetTags")
	 << " disc_tcHP=" << jet->bDiscriminator("trackCountingHighPurBJetTags")
	 << " disc_sftM=" << jet->bDiscriminator("softMuonBJetTags")
	 << " disc_sftMPt=" << jet->bDiscriminator("softMuonByPtBJetTags")
	 << " disc_sftMIP3=" << jet->bDiscriminator("softMuonByIP3dBJetTags")
	 << " disc_sftE=" << jet->bDiscriminator("softElectronBJetTags")
	 << " disc_jetP=" << jet->bDiscriminator("jetProbabilityBJetTags")
	 << " disc_jetBP=" << jet->bDiscriminator("jetBProbabilityBJetTags")
	 << " disc_ssv=" << jet->bDiscriminator("simpleSecondaryVertexBJetTags")
	 << " disc_csv=" << jet->bDiscriminator("combinedSecondaryVertexBJetTags")
	 << " disc_cvsMVA=" << jet->bDiscriminator("combinedSecondaryVertexMVABJetTags")
	 << " genId=" << id << " genPt=" << genPt
	 << endl;
    jetPt.push_back(jet->pt());
    jetEta.push_back(jet->eta());
    jetPhi.push_back(jet->phi());
    jetCorrFact.push_back(jet->corrFactor("raw"));
    jetDiscTCHE.push_back(jet->bDiscriminator("trackCountingHighEffBJetTags"));
    jetDiscTCHP.push_back(jet->bDiscriminator("trackCountingHighPurBJetTags"));
    jetDiscSFTM.push_back(jet->bDiscriminator("softMuonBJetTags"));
    jetDiscSMTMPT.push_back(jet->bDiscriminator("softMuonByPtBJetTags"));
    jetDiscSFTMIP3.push_back(jet->bDiscriminator("softMuonByIP3dBJetTags"));
    jetDiscSFTE.push_back(jet->bDiscriminator("softElectronBJetTags"));
    jetDiscJETP.push_back(jet->bDiscriminator("jetProbabilityBJetTags"));
    jetDiscJETBP.push_back(jet->bDiscriminator("jetBProbabilityBJetTags"));
    jetDiscSSV.push_back(jet->bDiscriminator("simpleSecondaryVertexBJetTags"));
    jetDiscCSV.push_back(jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    jetDiscCSVMVA.push_back(jet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
    jetGenId.push_back(id);
    jetGenPt.push_back(genPt);
    jetQ.push_back(jet->jetCharge());
    jetTks.push_back(jet->associatedTracks().size());
  }
  branch.jetPt=jetPt;
  branch.jetEta=jetEta;
  branch.jetPhi=jetPhi;
  branch.jetCorrFact=jetCorrFact;
  branch.jetDiscTCHE=jetDiscTCHE;
  branch.jetDiscTCHP=jetDiscTCHP;
  branch.jetDiscSFTM=jetDiscSFTM;
  branch.jetDiscSMTMPT=jetDiscSMTMPT;
  branch.jetDiscSFTMIP3=jetDiscSFTMIP3;
  branch.jetDiscSFTE=jetDiscSFTE;
  branch.jetDiscJETP=jetDiscJETP;
  branch.jetDiscJETBP=jetDiscJETBP;
  branch.jetDiscSSV=jetDiscSSV;
  branch.jetDiscCSV=jetDiscCSV;
  branch.jetDiscCSVMVA=jetDiscCSVMVA;
  branch.jetGenId=jetGenId;
  branch.jetGenPt=jetGenPt;
  branch.jetQ=jetQ;
  branch.jetTks=jetTks;

  ntuple->Fill();
}


void NtuplizerZAnalysis7TeVPAT::beginJob(){
  file = new TFile(outfile.c_str(),"recreate");
  
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  ntuple = new TTree("ntuple","sim2reco");
  ntuple->Branch("run",&(branch.run),"run/I");
  ntuple->Branch("event",&(branch.event),"event/I");
  ntuple->Branch("lumi",&(branch.lumi),"lumi/I");

  ntuple->Branch("trigBits",&(branch.trigBits));

  ntuple->Branch("vtxFake",&(branch.vtxFake));
  ntuple->Branch("vtxX",&(branch.vtxX));
  ntuple->Branch("vtxY",&(branch.vtxY));
  ntuple->Branch("vtxZ",&(branch.vtxZ));
  ntuple->Branch("vtxXE",&(branch.vtxXE));
  ntuple->Branch("vtxYE",&(branch.vtxYE));
  ntuple->Branch("vtxZE",&(branch.vtxZE));
  ntuple->Branch("vtxNdof",&(branch.vtxNdof));
  ntuple->Branch("vtxNtks",&(branch.vtxNtks));
  
  ntuple->Branch("muGLB",&(branch.muGLB));
  ntuple->Branch("muTRK",&(branch.muTRK));
  ntuple->Branch("muPt",&(branch.muPt));
  ntuple->Branch("muPtE",&(branch.muPtE));
  ntuple->Branch("muEta",&(branch.muEta));
  ntuple->Branch("muPhi",&(branch.muPhi));
  ntuple->Branch("muDxy",&(branch.muDxy));
  ntuple->Branch("muDxyE",&(branch.muDxyE));
  ntuple->Branch("muDz",&(branch.muDz));
  ntuple->Branch("muDzE",&(branch.muDzE));
  ntuple->Branch("muIsoTRK03",&(branch.muIsoTRK03));
  ntuple->Branch("muIsoECAL03",&(branch.muIsoECAL03));
  ntuple->Branch("muIsoHCAL03",&(branch.muIsoHCAL03));
  ntuple->Branch("muIsoHO03",&(branch.muIsoHO03));
  ntuple->Branch("muIsoTRK05",&(branch.muIsoTRK05));
  ntuple->Branch("muIsoECAL05",&(branch.muIsoECAL05));
  ntuple->Branch("muIsoHCAL05",&(branch.muIsoHCAL05));
  ntuple->Branch("muIsoHO05",&(branch.muIsoHO05));
  ntuple->Branch("muQ",&(branch.muQ));

  ntuple->Branch("elIdLoose",&(branch.elIdLoose));
  ntuple->Branch("elIdTight",&(branch.elIdTight));
  ntuple->Branch("elIdRobLo",&(branch.elIdRobLo));
  ntuple->Branch("elIdRobTi",&(branch.elIdRobTi));
  ntuple->Branch("elIdRobHE",&(branch.elIdRobHE));
  ntuple->Branch("elESC",&(branch.elESC));
  ntuple->Branch("elEoP",&(branch.elEoP));
  ntuple->Branch("elFbrem",&(branch.elFbrem));
  ntuple->Branch("elPt",&(branch.elPt));
  ntuple->Branch("elPtE",&(branch.elPtE));
  ntuple->Branch("elEta",&(branch.elEta));
  ntuple->Branch("elPhi",&(branch.elPhi));
  ntuple->Branch("elDxy",&(branch.elDxy));
  ntuple->Branch("elDxyE",&(branch.elDxyE));
  ntuple->Branch("elDz",&(branch.elDz));
  ntuple->Branch("elDzE",&(branch.elDzE));
  ntuple->Branch("elIsoTRK03",&(branch.elIsoTRK03));
  ntuple->Branch("elIsoECAL03",&(branch.elIsoECAL03));
  ntuple->Branch("elIsoHCAL03",&(branch.elIsoHCAL03));
  ntuple->Branch("elIsoTRK04",&(branch.elIsoTRK04));
  ntuple->Branch("elIsoECAL04",&(branch.elIsoECAL04));
  ntuple->Branch("elIsoHCAL04",&(branch.elIsoHCAL04));
  ntuple->Branch("elQ",&(branch.elQ));
  ntuple->Branch("elClass",&(branch.elClass));

  ntuple->Branch("metEt",&(branch.metEt));
  ntuple->Branch("metPhi",&(branch.metPhi));
  ntuple->Branch("metEtSig",&(branch.metEtSig));
  ntuple->Branch("metSig",&(branch.metSig));
  ntuple->Branch("metUncorr",&(branch.metUncorr));
  ntuple->Branch("metGenPt",&(branch.metGenPt));
  ntuple->Branch("metGenPhi",&(branch.metGenPhi));

  ntuple->Branch("pfMetEt",&(branch.pfMetEt));
  ntuple->Branch("pfMetPhi",&(branch.pfMetPhi));
  ntuple->Branch("pfMetEtSig",&(branch.pfMetEtSig));
  ntuple->Branch("pfMetSig",&(branch.pfMetSig));
  ntuple->Branch("pfMetUncorr",&(branch.pfMetUncorr));
  ntuple->Branch("pfMetGenPt",&(branch.pfMetGenPt));
  ntuple->Branch("pfMetGenPhi",&(branch.pfMetGenPhi));

  ntuple->Branch("jetPt",&(branch.jetPt));
  ntuple->Branch("jetEta",&(branch.jetEta));
  ntuple->Branch("jetPhi",&(branch.jetPhi));
  ntuple->Branch("jetCorrFact",&(branch.jetCorrFact));
  ntuple->Branch("jetDiscTCHE",&(branch.jetDiscTCHE));
  ntuple->Branch("jetDiscTCHP",&(branch.jetDiscTCHP));
  ntuple->Branch("jetDiscSFTM",&(branch.jetDiscSFTM));
  ntuple->Branch("jetDiscSMTMPT",&(branch.jetDiscSMTMPT));
  ntuple->Branch("jetDiscSFTMIP3",&(branch.jetDiscSFTMIP3));
  ntuple->Branch("jetDiscSFTE",&(branch.jetDiscSFTE));
  ntuple->Branch("jetDiscJETP",&(branch.jetDiscJETP));
  ntuple->Branch("jetDiscJETBP",&(branch.jetDiscJETBP));
  ntuple->Branch("jetDiscSSV",&(branch.jetDiscSSV));
  ntuple->Branch("jetDiscCSV",&(branch.jetDiscCSV));
  ntuple->Branch("jetDiscCSVMVA",&(branch.jetDiscCSVMVA));
  ntuple->Branch("jetGenId",&(branch.jetGenId));
  ntuple->Branch("jetGenPt",&(branch.jetGenPt));
  ntuple->Branch("jetQ",&(branch.jetQ));
  ntuple->Branch("jetTks",&(branch.jetTks));

  TH1::AddDirectory(oldAddDir);
}

void NtuplizerZAnalysis7TeVPAT::endJob() {
  file->Write();
  file->Close();
}

DEFINE_FWK_MODULE(NtuplizerZAnalysis7TeVPAT);
