#include "common.C"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

//these are applied only for plots, yields are produced with veto for WW and region for for ZZ
bool applyZVeto = false;
bool applyZReg  = true;
float minPtLW = 20.;
float kzzwz = 0.30/(0.38*0.60);

void arbitrate(LorentzVector l1, LorentzVector l2, LorentzVector l3,  LorentzVector met, 
	       int lid1, int lid2, int lid3, LorentzVector& lw, LorentzVector& l1z, LorentzVector& l2z) {
  //it is eee or mmm... we have to arbitrate :-(
  float mtboxhigh = 85.;
  float mtboxlow = 60.;
  //arbitrate between the SS
  if (lid1*lid3>0) {
    // 1 and 3 are SS
    //check mT
    float mt1 = 2*sqrt(l1.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l1,met)/2));
    float mt3 = 2*sqrt(l3.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l3,met)/2));
    if (TMath::Max(mt1,mt3)>mtboxhigh) {
      //if one has too large mT it is not from W
      if (mt1>mt3) {
	l1z = l1;
	l2z = l2;
	lw = l3;
      } else {
	l1z = l2;
	l2z = l3;
	lw = l1;
      }
    } else if (TMath::Min(mt1,mt3)<mtboxlow) {
      //if one has low mT, take it as fromZ
      if (mt1<mt3) {
	l1z = l1;
	l2z = l2;
	lw = l3;
      } else {
	l1z = l2;
	l2z = l3;
	lw = l1;
      }
    } else {
      //now... if both mTs are very close we do not know
      //both are in the box
      //for now do the same as previous case
      if (mt1<mt3) {
	l1z = l1;
	l2z = l2;
	lw = l3;
      } else {
	l1z = l2;
	l2z = l3;
	lw = l1;
      }
    }
  } else if (lid1*lid2>0) { 
    // 1 and 2 are SS
    //check mT
    float mt1 = 2*sqrt(l1.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l1,met)/2));
    float mt2 = 2*sqrt(l2.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l2,met)/2));
    if (TMath::Max(mt1,mt2)>mtboxhigh) {
      //if one has too large mT it is not from W
      if (mt1>mt2) {
	l1z = l1;
	l2z = l3;
	lw = l2;
      } else {
	l1z = l3;
	l2z = l2;
	lw = l1;
      }
    } else if (TMath::Min(mt1,mt2)<mtboxlow) {
      //if both have reasonable mT, take the largest mT as from W
      if (mt1<mt2) {
	l1z = l1;
	l2z = l3;
	lw = l2;
      } else {
	l1z = l3;
	l2z = l2;
	lw = l1;
      }
    } else {
      //now... if both mTs are very close we do not know
      //both are in the box
      //for now do the same as previous case
      if (mt1<mt2) {
	l1z = l1;
	l2z = l3;
	lw = l2;
      } else {
	l1z = l3;
	l2z = l2;
	lw = l1;
      }
    }
  } else if (lid2*lid3>0) {
    // 1 and 2 are SS
    //check mT
    float mt3 = 2*sqrt(l3.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l3,met)/2));
    float mt2 = 2*sqrt(l2.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l2,met)/2));
    if (TMath::Max(mt3,mt2)>mtboxhigh) {
      //if one has too large mT it is not from W
      if (mt3>mt2) {
	l1z = l3;
	l2z = l1;
	lw = l2;
      } else {
	l1z = l1;
	l2z = l2;
	lw = l3;
      }
    } else if (TMath::Min(mt3,mt2)<mtboxlow) {
      //if both have reasonable mT, take the largest mT as from W
      if (mt3<mt2) {
	l1z = l3;
	l2z = l1;
	lw = l2;
      } else {
	l1z = l1;
	l2z = l2;
	lw = l3;
      } 
    } else {
      //now... if both mTs are very close we do not know
      //both are in the box
      //for now do the same as previous case
      if (mt3<mt2) {
	l1z = l3;
	l2z = l1;
	lw = l2;
      } else {
	l1z = l1;
	l2z = l2;
	lw = l3;
      } 
    }
  }
} 

void plotMetForZZ(TString sample="/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/zz.root", 
		  int mass=0, unsigned int njets=0, float lumi=1.) {

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);

  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");
  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;

  TFile *fZZKFactorFile = TFile::Open("/smurf/sixie/KFactors/ZZ_KFactor.root");
  TH1D* ZZKFactor_ = (TH1D*)(fZZKFactorFile->Get("KFactorZZ_DileptonPt"));
  if (ZZKFactor_) {
    ZZKFactor_->SetDirectory(0);
  }
  assert(ZZKFactor_);
  fZZKFactorFile->Close();
  delete fZZKFactorFile;


  TFile* outf = TFile::Open("plotsFromZZ.root","RECREATE");

  TH1F* znnpt = new TH1F("znnpt","znnpt",20,0,200);
  TH1F* zllmt = new TH1F("zllmt","zllmt",20,0,200);
  TH1F* zllmt_zz = new TH1F("zllmt_zz","zllmt_zz",20,180,340);
  //TH1F* zllmt_zz = new TH1F("zllmt_zz","zllmt_zz",20,0,400);
  TH1F* zllm = new TH1F("zllm","zllm",20,0,200);
  TH1F* zlldphi = new TH1F("zlldphi","zlldphi",20,0,200);
  TH1F* zllminmet = new TH1F("zllminmet","zllminmet",20,0,200);
  TH1F* zllminpmet = new TH1F("zllminpmet","zllminpmet",20,0,200);
  znnpt->Sumw2();
  zllmt->Sumw2();
  zllmt_zz->Sumw2();
  zllm->Sumw2();
  zlldphi->Sumw2();
  zllminmet->Sumw2();
  zllminpmet->Sumw2();

  float weight = 1.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    //if (useJson && !passJson(dataEvent->run_,dataEvent->lumi_)) continue;    
    //if ( dataEvent->njets_!=njets) continue;

    //consider only di-lepton events
    if (abs(dataEvent->lq3_)==1) continue;

    if ( (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection ) continue;
    if (max(dataEvent->lep1_.pt(),dataEvent->lep2_.pt())<20.) continue;
    if (min(dataEvent->lep1_.pt(),dataEvent->lep2_.pt())<20.) continue;

    //require to pass WW selection but MET
    if ( (dataEvent->cuts_ & wwSelNoZVNoMet) != wwSelNoZVNoMet ) continue;
    if ( applyZVeto && ((dataEvent->cuts_ & ZVeto) != ZVeto) ) continue;
    if ( applyZReg  && ((dataEvent->cuts_ & ZVeto) == ZVeto) ) continue;

    //consider SF OS events only
    int id1 = dataEvent->lid1_;
    int id2 = dataEvent->lid2_;
    if ( (id1+id2)!=0 ) continue;

    //ok, I cannot be sure it is ZZ->2l2nu so I try to get rid of the other ZZ in this way
    //if (dataEvent->genmet_<10.) continue; 

    //double PtKFactor = ZZKFactor_->GetBinContent( ZZKFactor_->GetXaxis()->FindFixBin(dataEvent->dilep_.Pt()));
    //weight = weight * PtKFactor;

    float met    = dataEvent->met_;
    float metPhi = dataEvent->metPhi_;
    LorentzVector met_(met*cos(metPhi),met*sin(metPhi),0,met);

    if (met<68.5) continue;//new
    if (dataEvent->dilep_.pt()<25) continue;//new

    znnpt->Fill(met_.pt(),weight);
    zlldphi->Fill(dataEvent->dPhi_*180./TMath::Pi(),weight);
    zllmt->Fill(dataEvent->mt_,weight);
    zllm->Fill(dataEvent->dilep_.mass(),weight);
    float termA = sqrt( dataEvent->dilep_.Pt()*dataEvent->dilep_.Pt() + dataEvent->dilep_.M()*dataEvent->dilep_.M() );
    float termB = sqrt( dataEvent->met_*dataEvent->met_ + dataEvent->dilep_.M()*dataEvent->dilep_.M() );
    float newX = dataEvent->dilep_.Px() + dataEvent->met_ * cos(dataEvent->metPhi_);
    float newY = dataEvent->dilep_.Py() + dataEvent->met_ * sin(dataEvent->metPhi_);
    float termC = newX*newX + newY*newY;
    float mtZZ = sqrt( pow(termA + termB, 2) - termC );
    zllmt_zz->Fill(mtZZ,weight);

    zllminmet->Fill(min(dataEvent->met_,dataEvent->trackMet_),weight);
    zllminpmet->Fill(min(dataEvent->pmet_,dataEvent->pTrackMet_),weight);

  }
  outf->Write();
  outf->Close();
}

bool getWZcandidate(LorentzVector l1, LorentzVector l2, LorentzVector l3,  LorentzVector met, 
		    int id1, int id2, int id3, LorentzVector& lw, LorentzVector& l1z, LorentzVector& l2z) {
  bool isArbitrated = false;
  if (abs(id1)!=abs(id2) || abs(id1)!=abs(id3)) {;
    //case with eem/mme
    if ( (id1+id2)==0 ) {
      l1z = l1;
      l2z = l2;
      lw = l3;
    } else if ( (id1+id3)==0 ){
      l1z = l1;
      l2z = l3;
      lw = l2;
    } else if ( (id2+id3)==0 ){
      l1z = l2;
      l2z = l3;
      lw = l1;
    } else {
      //cout << "WTF! check charge" << endl;
      l1z = LorentzVector(0,0,0,0);
      l2z = LorentzVector(0,0,0,0);
      lw  = LorentzVector(0,0,0,0);
      return 0;
    }
    return isArbitrated;
  } else {
    //consider case eee or mmm
    isArbitrated = true;
    //run arbitration
    arbitrate(l1, l2, l3, met, id1, id2, id3, lw, l1z, l2z);
    return isArbitrated;    
  }//end of eee/mmm
}

bool checkArbitration(LorentzVector l1, LorentzVector l2, LorentzVector l3,
		      int l1MotherMcId, int l2MotherMcId, int l3MotherMcId, 
		      LorentzVector z) {
    //check if arbitration is success
    if (fabs(z.mass()-(l1+l2).mass())<0.0001) {
      //ok the choice is z=1+2
      if (l1MotherMcId==l2MotherMcId) {
	//it was the good choice
	return 1;
      } else {
	//it was NOT the good choice
	return 0;
      }
    } else if (fabs(z.mass()-(l1+l3).mass())<0.001) {
      //ok the choice is z=1+3
      if (l1MotherMcId==l3MotherMcId) {
	//it was the good choice
	return 1;
      } else {
	//it was NOT the good choice
	return 0;
      }
    } else if (fabs(z.mass()-(l2+l3).mass())<0.001) {
      //ok the choice is z=2+3
      if (l2MotherMcId==l3MotherMcId) {
	//it was the good choice
	return 1;
      } else {
	//it was NOT the good choice
	return 0;
       }
    } 
    cout << "WTF! should not be here" << endl;
    return 0;
}

float getMinProjMet(LorentzVector l1z, LorentzVector l2z, LorentzVector lw, float met_, float metPhi_, float trackMet_, float trackMetPhi_) {

  LorentzVector met(met_*cos(metPhi_),met_*sin(metPhi_),0,met_);
  LorentzVector trackMet(trackMet_*cos(trackMetPhi_),trackMet_*sin(trackMetPhi_),0,trackMet_);

  LorentzVector metw      = met+lw;
  LorentzVector trackMetw = trackMet+lw;

  float dPhiLMetw      = min( fabs(ROOT::Math::VectorUtil::DeltaPhi(l1z,metw)),      fabs(ROOT::Math::VectorUtil::DeltaPhi(l2z,metw)) );
  float dPhiLTrackMetw = min( fabs(ROOT::Math::VectorUtil::DeltaPhi(l1z,trackMetw)), fabs(ROOT::Math::VectorUtil::DeltaPhi(l2z,trackMetw)) );

  float pmetw      = dPhiLMetw<TMath::Pi()/2      ? metw.pt()*sin(dPhiLMetw)           : metw.pt();  
  float ptrackMetw = dPhiLTrackMetw<TMath::Pi()/2 ? trackMetw.pt()*sin(dPhiLTrackMetw) : trackMetw.pt();  

  return min(pmetw,ptrackMetw);

}

float getMinMet(LorentzVector lw, float met_, float metPhi_, float trackMet_, float trackMetPhi_) {

  LorentzVector met(met_*cos(metPhi_),met_*sin(metPhi_),0,met_);
  LorentzVector trackMet(trackMet_*cos(trackMetPhi_),trackMet_*sin(trackMetPhi_),0,trackMet_);

  LorentzVector metw      = met+lw;
  LorentzVector trackMetw = trackMet+lw;

  return min(metw.pt(),trackMetw.pt());

}

void plotMetForWZ(TString sample="/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/wz.root", 
		  int mass=0, unsigned int njets=0, float lumi=1.) {

  //plotMetForZZ();

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);


  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");
  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;

  TFile* outf;
  if (sample.Contains("wz")) outf = TFile::Open("plotsFromWZ_wz.root","RECREATE");
  else if (sample.Contains("dymm")) outf = TFile::Open("plotsFromWZ_dymm.root","RECREATE");
  else if (sample.Contains("dyee")) outf = TFile::Open("plotsFromWZ_dyee.root","RECREATE");
  else if (sample.Contains("dytt")) outf = TFile::Open("plotsFromWZ_dytt.root","RECREATE");
  else if (sample.Contains("ttbar_powheg"))outf = TFile::Open("plotsFromWZ_ttbar.root","RECREATE");
  else if (sample.Contains("qqww")) outf = TFile::Open("plotsFromWZ_qqww.root","RECREATE");
  else if (sample.Contains("wjets"))outf = TFile::Open("plotsFromWZ_wjets.root","RECREATE");
  else outf = TFile::Open("plotsFromWZ_data.root","RECREATE");

  TH1F* wpt = new TH1F("wpt","wpt",20,0,200);
  TH1F* zm = new TH1F("zm","zm",20,0,200);
  TH1F* zmt = new TH1F("zmt","zmt",20,0,200);
  TH1F* zmt_zz = new TH1F("zmt_zz","zmt_zz",20,180,340);
  //TH1F* zmt_zz = new TH1F("zmt_zz","zmt_zz",20,0,400);
  TH1F* zdphi = new TH1F("zdphi","zdphi",20,0,200);
  TH1F* zdeta = new TH1F("zdeta","zdeta",20,0,5);
  TH1F* zminmet = new TH1F("zminmet","zminmet",20,0,200);
  TH1F* zminpmet = new TH1F("zminpmet","zminpmet",20,0,200);
  wpt->Sumw2();
  zm->Sumw2();
  zmt->Sumw2();
  zmt_zz->Sumw2();
  zdphi->Sumw2();
  zdeta->Sumw2();
  zminmet->Sumw2();
  zminpmet->Sumw2();
  TH1F* arbeff = new TH1F("arbeff","arbeff",2,0,2);

  float weight = 1.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    else weight = 1.;
    weight = weight * kzzwz;

    //if (useJson && !passJson(dataEvent->run_,dataEvent->lumi_)) continue;    
    //if ( dataEvent->njets_!=njets) continue;
    //consider only tri-lepton events
    if (abs(dataEvent->lq3_)!=1) continue;
    if ( (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep3FullSelection) != Lep3FullSelection ) continue;

    //total charge has to be 1
    if (abs(dataEvent->lq1_+dataEvent->lq2_+dataEvent->lq3_)!=1) continue;
    //pt 20/10/10
    if (max(dataEvent->lep1_.pt(),max(dataEvent->lep2_.pt(),dataEvent->lep3_.pt()))<20.) continue;
    if (min(dataEvent->lep1_.pt(),min(dataEvent->lep2_.pt(),dataEvent->lep3_.pt()))<20.) continue;

    int id1 = dataEvent->lid1_;
    int id2 = dataEvent->lid2_;
    int id3 = dataEvent->lid3_;

    LorentzVector l1z;
    LorentzVector l2z;
    LorentzVector lw;
    float met    = dataEvent->met_;
    float metPhi = dataEvent->metPhi_;
    LorentzVector met_(met*cos(metPhi),met*sin(metPhi),0,met);

    bool arbitrated = getWZcandidate(dataEvent->lep1_, dataEvent->lep2_, dataEvent->lep3_,  met_, id1, id2, id3, lw, l1z, l2z);
    LorentzVector z = l1z+l2z;
    if (z.mass()<1E-5) continue;//e+e+m- or e-e-m+ or m+m+e- or m-m-e+
    float dphill = fabs(ROOT::Math::VectorUtil::DeltaPhi(l1z,l2z))*180./TMath::Pi(); 
    float detall = fabs(l2z.Eta() - l1z.Eta());
    LorentzVector w = lw+met_;
    if (arbitrated) arbeff->Fill(checkArbitration(dataEvent->lep1_, dataEvent->lep2_, dataEvent->lep3_,
						  dataEvent->lep1MotherMcId_, dataEvent->lep2MotherMcId_, dataEvent->lep3MotherMcId_, 
						  z));

    float zmt_ = 2*sqrt(z.pt()*w.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(z,w)/2));

    //require the lepton from W to be high pT to reject fakes (from DY)
    if (lw.pt()<minPtLW) continue;

    //apply other cuts to mimic WW selection
    if ( (dataEvent->cuts_ & TopVeto) != TopVeto ) continue;
    if ( applyZVeto && fabs(z.mass()-91.1876)<15.) continue;//z veto
    if ( applyZReg  && fabs(z.mass()-91.1876)>15.) continue;//z region

    if (w.pt()<68.5) continue;//new
    if (z.pt()<25) continue;//new

    wpt->Fill(w.pt(),weight);
    zdphi->Fill(dphill,weight);
    zdeta->Fill(detall,weight);
    zmt->Fill( zmt_,weight );
    zm->Fill(z.mass(),weight);
    //mt cut a la ZZ
    float termA = sqrt( z.pt()*z.pt() + z.mass()*z.mass() );
    float termB = sqrt( w.pt()*w.pt() + z.mass()*z.mass() );
    float newX = z.px() + w.pt() * cos(w.phi());
    float newY = z.py() + w.pt() * sin(w.phi());
    float termC = newX*newX + newY*newY;
    float mtZZ = sqrt( pow(termA + termB, 2) - termC );
    zmt_zz->Fill( mtZZ,weight );
    float minmet = getMinMet(lw, dataEvent->met_, dataEvent->metPhi_, dataEvent->trackMet_, dataEvent->trackMetPhi_);
    float minpmet = getMinProjMet(l1z, l2z, lw, dataEvent->met_, dataEvent->metPhi_, dataEvent->trackMet_, dataEvent->trackMetPhi_);
    zminmet->Fill( minmet,weight);
    zminpmet->Fill( minpmet,weight);

  }
  outf->Write();
  outf->Close();
}

void getYieldZZatWWlevelNoMet(TString sample="/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/zz.root", 
			      int mass=0, int njets=0, float lumi=1.) {

  pair<float, float> zz = getYield(sample, wwSelNoMet, noVeto, mass, njets, "", lumi, 0, 0, 0);
  pair<float, float> zzGen10 = getYield(sample, wwSelNoMet, noVeto, mass, njets, "genMetGt10", lumi, 0, 0, 0);
  pair<float, float> zzGen10mm40 = getYield(sample, wwSelNoMet, noVeto, mass, njets, "genMetGt10,minmet40", lumi, 0, 0, 0);
  pair<float, float> zzGen10MH140_1 = getYield(sample, wwSelNoMet, noVeto, 140, njets, "genMetGt10,minmet40,leppts", lumi, 0, 0, 0);
  pair<float, float> zzGen10MH140_2 = getYield(sample, wwSelNoMet, noVeto, 140, njets, "genMetGt10,minmet40,masscut,leppts", lumi, 0, 0, 0);
  pair<float, float> zzGen10MH140_3 = getYield(sample, wwSelNoMet, noVeto, 140, njets, "genMetGt10,minmet40,dphicut,masscut,leppts", lumi, 0, 0, 0);
  pair<float, float> zzGen10MH140_4 = getYield(sample, wwSelNoMet, noVeto, 140, njets, "genMetGt10,minmet40,dphicut,mtcut,masscut,leppts", lumi, 0, 0, 0);


  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC"                    ,
	       zz.first, zz.second, 
	       100*zz.first/zz.first, 0.) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (genMet>10)",
	       zzGen10.first, zzGen10.second, 
	       100*zzGen10.first/zzGen10.first, 100*sqrt( zzGen10.first/zzGen10.first*(1-zzGen10.first/zzGen10.first)/zzGen10.first )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + minpmet>40)",
	       zzGen10mm40.first, zzGen10mm40.second, 
	       100*zzGen10mm40.first/zzGen10.first, 100*sqrt( zzGen10mm40.first/zzGen10.first*(1-zzGen10mm40.first/zzGen10.first)/zzGen10.first )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + MH140 leppt)",
	       zzGen10MH140_1.first, zzGen10MH140_1.second, 
	       100*zzGen10MH140_1.first/zzGen10mm40.first, 100*sqrt( zzGen10MH140_1.first/zzGen10mm40.first*(1-zzGen10MH140_1.first/zzGen10mm40.first)/zzGen10mm40.first )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + MH140 mll)",
	       zzGen10MH140_2.first, zzGen10MH140_2.second, 
	       100*zzGen10MH140_2.first/zzGen10MH140_1.first, 100*sqrt( zzGen10MH140_2.first/zzGen10MH140_1.first*(1-zzGen10MH140_2.first/zzGen10MH140_1.first)/zzGen10MH140_1.first )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + MH140 dphill)",
	       zzGen10MH140_3.first, zzGen10MH140_3.second, 
	       100*zzGen10MH140_3.first/zzGen10MH140_2.first, 100*sqrt( zzGen10MH140_3.first/zzGen10MH140_2.first*(1-zzGen10MH140_3.first/zzGen10MH140_2.first)/zzGen10MH140_2.first )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + MH140 mt)",
	       zzGen10MH140_4.first, zzGen10MH140_4.second, 
	       100*zzGen10MH140_4.first/zzGen10MH140_3.first, 100*sqrt( zzGen10MH140_4.first/zzGen10MH140_3.first*(1-zzGen10MH140_4.first/zzGen10MH140_3.first)/zzGen10MH140_3.first )) << endl;

}

void getYieldZZatWWlevelNoMetFromWZ(TString sample="/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/wz.root", 
				    int mass=0, unsigned int njets=0, float lumi=1.) {

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);

  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");
  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;


  float weight = 1.;
  float yield = 0.;
  float error = 0.;
  float yieldGt40 = 0.;
  float errorGt40 = 0.;
  float yieldMH140_1 = 0.;
  float errorMH140_1 = 0.;
  float yieldMH140_2 = 0.;
  float errorMH140_2 = 0.;
  float yieldMH140_3 = 0.;
  float errorMH140_3 = 0.;
  float yieldMH140_4 = 0.;
  float errorMH140_4 = 0.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    else weight = 1.;
    weight = weight * kzzwz;

    //if (useJson && !passJson(dataEvent->run_,dataEvent->lumi_)) continue;    
    if ( dataEvent->njets_!=njets) continue;
    //consider only tri-lepton events
    if (abs(dataEvent->lq3_)!=1) continue;
    if ( (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep3FullSelection) != Lep3FullSelection ) continue;

    //total charge has to be 1
    if (abs(dataEvent->lq1_+dataEvent->lq2_+dataEvent->lq3_)!=1) continue;
    //pt 20/10/10
    if (max(dataEvent->lep1_.pt(),max(dataEvent->lep2_.pt(),dataEvent->lep3_.pt()))<20.) continue;
    if (min(dataEvent->lep1_.pt(),min(dataEvent->lep2_.pt(),dataEvent->lep3_.pt()))<10.) continue;

    int id1 = dataEvent->lid1_;
    int id2 = dataEvent->lid2_;
    int id3 = dataEvent->lid3_;

    LorentzVector l1z;
    LorentzVector l2z;
    LorentzVector lw;
    float met    = dataEvent->met_;
    float metPhi = dataEvent->metPhi_;
    LorentzVector met_(met*cos(metPhi),met*sin(metPhi),0,met);

    bool arbitrated = false;
    arbitrated = getWZcandidate(dataEvent->lep1_, dataEvent->lep2_, dataEvent->lep3_,  met_, id1, id2, id3, lw, l1z, l2z);
    LorentzVector z = l1z+l2z;
    if (z.mass()<1E-5) continue;//e+e+m- or e-e-m+ or m+m+e- or m-m-e+
    float dphill = fabs(ROOT::Math::VectorUtil::DeltaPhi(l1z,l2z))*180./TMath::Pi();
    LorentzVector w = lw+met_;

    //require the lepton from W to be high pT to reject fakes (from DY)
    if (lw.pt()<minPtLW) continue;

    //apply other cuts to mimic WW selection
    if ( (dataEvent->cuts_ & TopVeto) != TopVeto ) continue;
    if ( fabs(z.mass()-91.1876)<15.) continue;//z veto
    //dphijet
    if (dataEvent->jet1_.pt()>15 && z.phi()>165.*TMath::Pi()/180.) continue;

    yield = yield+weight;
    error = error+pow(weight,2);

    float mtllmet = 2*sqrt(z.pt()*w.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(z,w)/2));
    float minpmet = getMinProjMet(l1z, l2z, lw, dataEvent->met_, dataEvent->metPhi_, dataEvent->trackMet_, dataEvent->trackMetPhi_);

    if (minpmet>40.) {
      yieldGt40 = yieldGt40+weight;
      errorGt40 = errorGt40+pow(weight,2);
      if ( min(l1z.pt(),l2z.pt())>15. && max(l1z.pt(),l2z.pt())>25. ) {
	yieldMH140_1 = yieldMH140_1+weight;
	errorMH140_1 = errorMH140_1+pow(weight,2);
	if (z.mass()<45.) {
	  yieldMH140_2 = yieldMH140_2+weight;
	  errorMH140_2 = errorMH140_2+pow(weight,2);
	  if (dphill<90.) {
	    yieldMH140_3 = yieldMH140_3+weight;
	    errorMH140_3 = errorMH140_3+pow(weight,2);
	    if (mtllmet>80. && mtllmet<130.) {
	      yieldMH140_4 = yieldMH140_4+weight;
	      errorMH140_4 = errorMH140_4+pow(weight,2);
	    }
	  }
	}
      }
    }
    
  }

  error = sqrt(error);
  errorGt40 = sqrt(errorGt40);
  errorMH140_1 = sqrt(errorMH140_1);
  errorMH140_2 = sqrt(errorMH140_2);
  errorMH140_3 = sqrt(errorMH140_3);
  errorMH140_4 = sqrt(errorMH140_4);

  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (R-corrected)"               ,
	       yield       ,error       ,100*yield/yield              ,0.) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + met>40)"      ,
	       yieldGt40   ,errorGt40   ,100*yieldGt40/yield          ,100*sqrt( yieldGt40/yield*(1-yieldGt40/yield)/yield )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + MH140 leppt)" ,
	       yieldMH140_1,errorMH140_1,100*yieldMH140_1/yieldGt40   ,100*sqrt( yieldMH140_1/yieldGt40*(1-yieldMH140_1/yieldGt40)/yieldGt40 )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + MH140 mll)"   ,
	       yieldMH140_2,errorMH140_2,100*yieldMH140_2/yieldMH140_1,100*sqrt( yieldMH140_2/yieldMH140_1*(1-yieldMH140_2/yieldMH140_1)/yieldMH140_1 )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + MH140 dphill)",
	       yieldMH140_3,errorMH140_3,100*yieldMH140_3/yieldMH140_2,100*sqrt( yieldMH140_3/yieldMH140_2*(1-yieldMH140_3/yieldMH140_2)/yieldMH140_2 )) << endl;
  cout << Form("%-40s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + MH140 mt)"    ,
	       yieldMH140_4,errorMH140_4,100*yieldMH140_4/yieldMH140_3,100*sqrt( yieldMH140_4/yieldMH140_3*(1-yieldMH140_4/yieldMH140_3)/yieldMH140_3 )) << endl;

}


void getYieldZZforZZ(TString sample="/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/zz.root", 
		     int mass=0, int njets=0, float lumi=1.) {

  unsigned int baseZZSel = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|ExtraLeptonVeto|TopVeto; 

  pair<float, float> zz   = getYield(sample, baseZZSel, ZVeto, mass, njets, "dphijet,dileppt40,l12pt2020", lumi, 0, 0, 0);
  pair<float, float> zz_1 = getYield(sample, baseZZSel, ZVeto, mass, njets, "dphijet,dileppt40,minmetGt50,l12pt2020", lumi, 0, 0, 0);
  pair<float, float> zz_2 = getYield(sample, baseZZSel, ZVeto, mass, njets, "dphijet,dileppt40,minmetGt50,l12pt2020,mtZZmH300", lumi, 0, 0, 0);
  pair<float, float> zz_3 = getYield(sample, baseZZSel, ZVeto, mass, njets, "dphijet,dileppt40,minmetGt70,l12pt2020,mtZZmH300", lumi, 0, 0, 0);

  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC"                    ,
	       zz.first       ,zz.second       ,100*zz.first/zz.first              ,0.) << endl;
  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + met>50)",
	       zz_1.first       ,zz_1.second       ,100*zz_1.first/zz.first              ,100*sqrt( zz_1.first/zz.first*(1-zz_1.first/zz.first)/zz.first )) << endl;
  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + mt cut)",
	       zz_2.first       ,zz_2.second       ,100*zz_2.first/zz_1.first              ,100*sqrt( zz_2.first/zz_1.first*(1-zz_2.first/zz_1.first)/zz_1.first )) << endl;
  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from MC (as above + met>70)",
	       zz_3.first       ,zz_3.second       ,100*zz_3.first/zz_2.first              ,100*sqrt( zz_3.first/zz_2.first*(1-zz_3.first/zz_2.first)/zz_2.first )) << endl;

}

void getYieldZZforZZFromWZ(TString sample="/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/wz.root", 
			   int mass=0, unsigned int njets=0, float lumi=1.) {

  float lep1pt=0.,lep2pt=0.,dPhi=0.,mll=0.,mtL=0.,mtH=0.,himass=0.;
  getCutValues(mass,lep1pt,lep2pt,dPhi,mll,mtL,mtH,himass);
  
  if (sample.Contains(".root")) sample=sample.ReplaceAll(".root","");
  SmurfTree *dataEvent = new SmurfTree();
  dataEvent->LoadTree(sample+".root");
  dataEvent->InitTree();

  bool isMC = lumi>1E-5;

  float weight = 1.;
  float yield = 0.;
  float error = 0.;
  float yieldGt50 = 0.;
  float errorGt50 = 0.;
  float yieldMH300_1 = 0.;
  float errorMH300_1 = 0.;
  float yieldMH300_2 = 0.;
  float errorMH300_2 = 0.;
  for(UInt_t n=0; n < dataEvent->tree_->GetEntries() ; ++n) {
    
    dataEvent->tree_->GetEntry(n);
    if (isMC) weight = lumi*dataEvent->scale1fb_;
    else weight = 1.;
    weight = weight * kzzwz;

    //if (useJson && !passJson(dataEvent->run_,dataEvent->lumi_)) continue;    
    if ( dataEvent->njets_!=njets) continue;
    //consider only tri-lepton events
    if (abs(dataEvent->lq3_)!=1) continue;
    if ( (dataEvent->cuts_ & Lep1FullSelection) != Lep1FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep2FullSelection) != Lep2FullSelection ) continue;
    if ( (dataEvent->cuts_ & Lep3FullSelection) != Lep3FullSelection ) continue;

    //total charge has to be 1
    if (abs(dataEvent->lq1_+dataEvent->lq2_+dataEvent->lq3_)!=1) continue;
    //pt 20/10/10
    if (max(dataEvent->lep1_.pt(),max(dataEvent->lep2_.pt(),dataEvent->lep3_.pt()))<20.) continue;
    if (min(dataEvent->lep1_.pt(),min(dataEvent->lep2_.pt(),dataEvent->lep3_.pt()))<10.) continue;

    int id1 = dataEvent->lid1_;
    int id2 = dataEvent->lid2_;
    int id3 = dataEvent->lid3_;

    LorentzVector l1z;
    LorentzVector l2z;
    LorentzVector lw;
    float met    = dataEvent->met_;
    float metPhi = dataEvent->metPhi_;
    LorentzVector met_(met*cos(metPhi),met*sin(metPhi),0,met);

    bool arbitrated = false;
    arbitrated = getWZcandidate(dataEvent->lep1_, dataEvent->lep2_, dataEvent->lep3_,  met_, id1, id2, id3, lw, l1z, l2z);
    LorentzVector z = l1z+l2z;
    if (z.mass()<1E-5) continue;//e+e+m- or e-e-m+ or m+m+e- or m-m-e+
    LorentzVector w = lw+met_;

    //require the lepton from W to be high pT to reject fakes (from DY)
    if (lw.pt()<minPtLW) continue;

    //require 20/20 for Z
    if (l1z.pt()<20.) continue;
    if (l2z.pt()<20.) continue;

    //require zpt>40
    if (z.pt()<40.) continue;

    //dphijet
    if (dataEvent->jet1_.pt()>15 && z.phi()>165.*TMath::Pi()/180.) continue;

    //apply other cuts to mimic WW selection
    if ( (dataEvent->cuts_ & TopVeto) != TopVeto ) continue;
    if ( fabs(z.mass()-91.1876)>15.) continue;//z requirement

    yield = yield+weight;
    error = error+pow(weight,2);

    float minmet = getMinMet(lw, dataEvent->met_, dataEvent->metPhi_, dataEvent->trackMet_, dataEvent->trackMetPhi_);

    if (minmet>50.) {
      yieldGt50 = yieldGt50+weight;
      errorGt50 = errorGt50+pow(weight,2);

      //mt cut a la ZZ
      float termA = sqrt( z.pt()*z.pt() + z.mass()*z.mass() );
      float termB = sqrt( w.pt()*w.pt() + z.mass()*z.mass() );
      float newX = z.px() + w.pt() * cos(w.phi());
      float newY = z.py() + w.pt() * sin(w.phi());
      float termC = newX*newX + newY*newY;
      float mt = sqrt( pow(termA + termB, 2) - termC );
      if ( mt>260. && mt<320. ) {
	yieldMH300_1 = yieldMH300_1+weight;
	errorMH300_1 = errorMH300_1+pow(weight,2);
	if (minmet>70.) {
	  yieldMH300_2 = yieldMH300_2+weight;
	  errorMH300_2 = errorMH300_2+pow(weight,2);
	}
      }

    }
    
  }

  error = sqrt(error);
  errorGt50 = sqrt(errorGt50);
  errorMH300_1 = sqrt(errorMH300_1);
  errorMH300_2 = sqrt(errorMH300_2);

  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (R-corrected)"         ,
	       yield       ,error       ,100*yield/yield              ,0.) << endl;
  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + met>50)",
	       yieldGt50   ,errorGt50   ,100*yieldGt50/yield          ,100*sqrt( yieldGt50/yield*(1-yieldGt50/yield)/yield )) << endl;
  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + mt cut)",
	       yieldMH300_1,errorMH300_1,100*yieldMH300_1/yieldGt50   ,100*sqrt( yieldMH300_1/yieldGt50*(1-yieldMH300_1/yieldGt50)/yieldGt50 )) << endl;
  cout << Form("%-30s: %6.2f +/- %-6.2f - eff: %6.2f +/- %-6.2f","zz from WZ (as above + met>70)",
	       yieldMH300_2,errorMH300_2,100*yieldMH300_2/yieldMH300_1,100*sqrt( yieldMH300_2/yieldMH300_1*(1-yieldMH300_2/yieldMH300_1)/yieldMH300_1 )) << endl;
  
}



/*
void arbitrateOld(LorentzVector l1, LorentzVector l2, LorentzVector l3,  LorentzVector met, 
		  int lid1, int lid2, int lid3, LorentzVector& w, LorentzVector& z, float& dphill) {
  //it is eee or mmm... we have to arbitrate :-(
  float zwindow = 0;
  float mtthreshold = 85.;
  //float mtboxlow = 60.;
  //arbitrate between the SS
  if (lid1*lid3>0) {
    // 1 and 3 are SS
    //check mass
    float m12 = (l1+l2).mass();
    float m23 = (l2+l3).mass();
    if ( TMath::Min( fabs(m12-91.), fabs(m23-91.) )<zwindow ) {
      //ok, one has mass close to z, take it as z
      if ( fabs(m12-91.) < fabs(m23-91.) ) {
	z = l1+l2;
	w = l3+met;
	dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l2);
      } else {
	z = l2+l3;
	w = l1+met;
	dphill = ROOT::Math::VectorUtil::DeltaPhi(l2,l3);
      }
    } else {
      //none is close to z mass, pick it by mT
      float mt1 = 2*sqrt(l1.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l1,met)/2));
      float mt3 = 2*sqrt(l3.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l3,met)/2));
      if (TMath::Max(mt1,mt3)>mtthreshold) {
	//if one has too large mT it is not from W
	if (mt1>mt3) {
	  z = l1+l2;
	  w = l3+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l2);
	} else {
	  z = l2+l3;
	  w = l1+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l2,l3);
	}
      } else {
	//if both have reasonable mT, take the largest mT as from W
	if (mt1<mt3) {
	  z = l1+l2;
	  w = l3+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l2);
	} else {
	  z = l2+l3;
	  w = l1+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l2,l3);
	}
      }
    }
  } else if (lid1*lid2>0) { 
    // 1 and 2 are SS
    //check mass
    float m13 = (l1+l3).mass();
    float m32 = (l3+l2).mass();
    if ( TMath::Min( fabs(m13-91.), fabs(m32-91.) )<zwindow ) {
      //ok, one has mass close to z, take it as z
      if ( fabs(m13-91.) < fabs(m32-91.) ) {
	z = l1+l3;
	w = l2+met;
	dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l3);
      } else {
	z = l3+l2;
	w = l1+met;
	dphill = ROOT::Math::VectorUtil::DeltaPhi(l2,l3);
      }
    } else {
      //none is close to z mass, pick it by mT
      float mt1 = 2*sqrt(l1.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l1,met)/2));
      float mt2 = 2*sqrt(l2.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l2,met)/2));
      if (TMath::Max(mt1,mt2)>mtthreshold) {
	//if one has too large mT it is not from W
	if (mt1>mt2) {
	  z = l1+l3;
	  w = l2+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l3);
	} else {
	  z = l3+l2;
	  w = l1+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l2,l3);
	}
      } else {
	//if both have reasonable mT, take the largest mT as from W
	if (mt1<mt2) {
	  z = l1+l3;
	  w = l2+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l3);
	} else {
	  z = l3+l2;
	  w = l1+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l2,l3);
	}
      }
    }
  } else if (lid2*lid3>0) {
    // 1 and 2 are SS
    //check mass
    float m31 = (l3+l1).mass();
    float m12 = (l1+l2).mass();
    if ( TMath::Min( fabs(m31-91.), fabs(m12-91.) )<zwindow ) {
      //ok, one has mass close to z, take it as z
      if ( fabs(m31-91.) < fabs(m12-91.) ) {
	z = l3+l1;
	w = l2+met;
	dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l3);
      } else {
	z = l1+l2;
	w = l3+met;
	dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l2);
      }
    } else {
      //none is close to z mass, pick it by mT
      float mt3 = 2*sqrt(l3.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l3,met)/2));
      float mt2 = 2*sqrt(l2.pt()*met.pt())*fabs(sin(ROOT::Math::VectorUtil::DeltaPhi(l2,met)/2));
      if (TMath::Max(mt3,mt2)>mtthreshold) {
	//if one has too large mT it is not from W
	if (mt3>mt2) {
	  z = l3+l1;
	  w = l2+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l3);
	} else {
	  z = l1+l2;
	  w = l3+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l2);
	}
      } else {
	//if both have reasonable mT, take the largest mT as from W
	if (mt3<mt2) {
	  z = l3+l1;
	  w = l2+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l3);
	} else {
	  z = l1+l2;
	  w = l3+met;
	  dphill = ROOT::Math::VectorUtil::DeltaPhi(l1,l2);
	}
      }
    }
  }
} 
*/
