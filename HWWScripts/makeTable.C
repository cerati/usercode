void makeTable(float lumi=2.121, int njets=0, int mass=0, bool useSF=false, bool dodata=false){
  //lumi is in /fb

  gROOT->Reset();

  TCut lep1pt,lep2pt,dPhi,mll,mt,himass;
  if (mass==0) {
    lep1pt = "lep1.pt()>20.";
    lep2pt = "lep2.pt()>10.";
    dPhi = "dPhi<TMath::Pi()*180./180.";
    mll = "dilep.mass()<999";
    mt = "mt>0&&mt<999";
    himass = "dilep.mass()>100.";
  } else if (mass==120) {
    lep1pt = "lep1.pt()>20.";
    lep2pt = "lep2.pt()>10.";
    dPhi = "dPhi<2.0";
    mll = "dilep.mass()<40";
    mt = "mt>70&&mt<120";
    himass = "dilep.mass()>100.";
  } else if (mass==130) {
    lep1pt = "lep1.pt()>25.";
    lep2pt = "lep2.pt()>10.";
    dPhi = "dPhi<1.5";
    mll = "dilep.mass()<45";
    mt = "mt>75&&mt<125";
    himass = "dilep.mass()>100.";
  } else if (mass==140) {
    lep1pt = "lep1.pt()>25.";
    lep2pt = "lep2.pt()>15.";
    dPhi = "dPhi<1.57";
    mll = "dilep.mass()<45";
    mt = "mt>80&&mt<130";
    himass = "dilep.mass()>100.";
  } else if (mass==150) {
    lep1pt = "lep1.pt()>27.";
    lep2pt = "lep2.pt()>25.";
    dPhi = "dPhi<1.57";
    mll = "dilep.mass()<50";
    mt = "mt>80&&mt<150";
    himass = "dilep.mass()>100.";
  } else if (mass==160) {
    lep1pt = "lep1.pt()>30.";
    lep2pt = "lep2.pt()>25.";
    dPhi = "dPhi<TMath::Pi()*60./180.";
    mll = "dilep.mass()<50";
    mt = "mt>90&&mt<160";
    himass = "dilep.mass()>100.";
  } else if (mass==200) {
    lep1pt = "lep1.pt()>40.";
    lep2pt = "lep2.pt()>25.";
    dPhi = "dPhi<TMath::Pi()*100./180.";
    mll = "dilep.mass()<90";
    mt = "mt>120&&mt<200";
    himass = "dilep.mass()>100.";
  } else if (mass==250) {
    lep1pt = "lep1.pt()>55.";
    lep2pt = "lep2.pt()>25.";
    dPhi = "dPhi<2.44";
    mll = "dilep.mass()<150";
    mt = "mt>120&&mt<250";
    himass = "dilep.mass()>100.";
  } else if (mass==300) {
    lep1pt = "lep1.pt()>70.";
    lep2pt = "lep2.pt()>25.";
    dPhi = "dPhi<3.05";
    mll = "dilep.mass()<200";
    mt = "mt>120&&mt<300";
    himass = "dilep.mass()>100.";
  }

  TCut ptreg   = lep1pt && lep2pt;
  TCut massreg = lep1pt && lep2pt && mll;
  TCut mtreg   = lep1pt && lep2pt && mll && mt;
  TCut sigreg  = lep1pt && lep2pt && mll && mt && dPhi;

  TCut sideband = lep1pt && lep2pt && himass;
  TCut mtside   = lep1pt && lep2pt && himass && mt;
  TCut dphiside = lep1pt && lep2pt && himass && mt && dPhi;

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
  unsigned int wwSelectionNoTV = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ZVeto|ExtraLeptonVeto;
  unsigned int wwSelection     = wwSelectionNoTV|TopVeto;

  TCut jets(Form("njets==%i",njets));
  TCut base(Form("(cuts & %i)==%i",wwSelection,wwSelection));
  TCut baseNTV(Form("(cuts & %i)==%i",wwSelectionNoTV,wwSelectionNoTV));
  TCut notTagNotInJets(Form("(cuts & %i)!=%i",TopTagNotInJets,TopTagNotInJets));
  TCut trig(Form("dstype!=0 || (cuts & %i)==%i",Trigger,Trigger));
  TCut newcuts = "type==1 || type==2 || ( lep2.pt()>15. && min(pmet,pTrackMet)>(37.+nvtx/2.) && (jet1.pt()<15 || dPhiDiLepJet1*180./TMath::Pi()<165.) )";
  TCut kincuts = "dilep.pt()>45.";

  TCut cut = base&&jets&&newcuts&&sigreg&&trig&&kincuts;
  //cout << "cut: " << cut.GetTitle() << endl;

  TString mcs[] = {"qqww","ggww","dyee","dymm","dytt","ttbar","tw","wz","zz","wjets","wgamma"};
  int  colors[] = {kCyan,kCyan,kGreen,kGreen,kGreen,kYellow,kYellow,kBlue,kBlue,kGray,kGray};
  float sfs0j[] = { 1.0,   1.0,   3.0,   3.0,   1.0,   1.5,    1.5, 1.0, 1.0, 1.1,    1.0};
  float sfs1j[] = { 1.0,   1.0,   2.8,   2.8,   1.0,   1.2,    1.2, 1.0, 1.0, 2.5,    1.0};
  float nosfs[] = { 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0, 1.0, 1.0, 1.0,    1.0};
  float* sfs;
  if (useSF&&njets==0) sfs = sfs0j;
  else if (useSF&&njets==1) sfs = sfs1j;
  else sfs = nosfs;

  TString dir = "/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets/";
  //TString dir = "/smurf/cerati/skims/Run2011_Spring11_SmurfV7_42X/2121ipb/wwSelNoLepNoTV/";
  TString data_file = "/smurf/data/LP2011/tas/data.root";
  int nMC = sizeof(mcs)/sizeof(TString);

  TCanvas c;

  float nSigMuMu=0,nSigElMu=0,nSigMuEl=0,nSigElEl=0,nSigAll=0;
  float nBckMuMu=0,nBckElMu=0,nBckMuEl=0,nBckElEl=0,nBckAll=0;
  TString header = TString(Form("%25s | %10s | %10s | %10s | %10s | %10s |","sample","MuMu","ElMu","MuEl","ElEl","All"));
  cout << header << endl;
  TString line = TString("--------------------------------------------------------------------------------------------");
  cout << line << endl;
  for (int i=0;i<nMC;++i){    
    TFile *_mc = TFile::Open(dir+mcs[i]+".root");
    TTree * mc = (TTree*) _mc->Get("tree");
    float correction = sfs[i];
    if (lumi>0.) mc->SetWeight(lumi*correction);
    //cout << scale1fb_ << endl;
    mc->Draw("type>>plotmc(4,0,4)",Form("scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s)",cut.GetTitle()),"g");
    float nMuMu = plotmc->GetBinContent(1);
    float nElMu = plotmc->GetBinContent(2);
    float nMuEl = plotmc->GetBinContent(3);
    float nElEl = plotmc->GetBinContent(4);
    float nAll  = nMuMu+nElMu+nMuEl+nElEl;
    string test = mcs[i];
    TString print = "";
    print = Form("%25s | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f |",test.c_str(),nMuMu,nElMu,nMuEl,nElEl,nAll);
    cout << print << endl;
    if (mcs[i].Contains("hww")) {
      nSigMuMu+=nMuMu;
      nSigElMu+=nElMu;
      nSigMuEl+=nMuEl;
      nSigElEl+=nElEl;
      nSigAll +=nAll;
    } else {
      nBckMuMu+=nMuMu;
      nBckElMu+=nElMu;
      nBckMuEl+=nMuEl;
      nBckElEl+=nElEl;
      nBckAll +=nAll;
    }
  }
  cout << line << endl;
  TString bck = TString(Form("%25s | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f |","tot. bck.",nBckMuMu,nBckElMu,nBckMuEl,nBckElEl,nBckAll));
  cout << bck << endl;
  cout << line << endl;

  if (dodata) {
    TFile *_data = TFile::Open(data_file);
    TTree * data = (TTree*) _data->Get("tree");
    data->Draw("type>>plotdata(4,0,4)",cut,"g");
    int ndataMuMu = plotdata->GetBinContent(1);
    int ndataMuEl = plotdata->GetBinContent(2);
    int ndataElMu = plotdata->GetBinContent(3);
    int ndataElEl = plotdata->GetBinContent(4);
    int ndataAll  = ndataMuMu+ndataMuEl+ndataElMu+ndataElEl;
    TString printdata = "";
    printdata = Form("%25s | %10i | %10i | %10i | %10i | %10i |","data",ndataMuMu,ndataMuEl,ndataElMu,ndataElEl,ndataAll);    
    cout << printdata << endl;
    cout << line << endl;
  }

  TString sig = TString(Form("%25s | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f |","tot. sig.",nSigMuMu,nSigElMu,nSigMuEl,nSigElEl,nSigAll));
  cout << sig << endl;
  cout << line << endl;

  if (!dodata) {
    TString sOb = TString(Form("%25s | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f |","s/b",nSigMuMu/nBckMuMu,nSigElMu/nBckElMu,nSigMuEl/nBckMuEl,nSigElEl/nBckElEl,nSigAll/nBckAll));
    cout << sOb << endl;
    cout << line << endl;
    TString fom = TString(Form("%25s | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f |","s/sqrt(s+b)",
			       nSigMuMu/sqrt(nSigMuMu+nBckMuMu),nSigElMu/sqrt(nSigElMu+nBckElMu),
			       nSigMuEl/sqrt(nSigMuEl+nBckMuEl),nSigElEl/sqrt(nSigElEl+nBckElEl),nSigAll/sqrt(nSigAll+nBckAll)));
    cout << fom << endl;
    cout << line << endl;
    TString fom = TString(Form("%25s | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f |","s/sqrt(s+b+(0.35*b)^2)",
			       nSigMuMu/sqrt(nSigMuMu+nBckMuMu+0.35*0.35*nBckMuMu*nBckMuMu),
			       nSigElMu/sqrt(nSigElMu+nBckElMu+0.35*0.35*nBckElMu*nBckElMu),
			       nSigMuEl/sqrt(nSigMuEl+nBckMuEl+0.35*0.35*nBckMuEl*nBckMuEl),
			       nSigElEl/sqrt(nSigElEl+nBckElEl+0.35*0.35*nBckElEl*nBckElEl),
			       nSigAll/sqrt(nSigAll+nBckAll+0.35*0.35*nBckAll*nBckAll)));
    cout << fom << endl;
    cout << line << endl;
  }

}

