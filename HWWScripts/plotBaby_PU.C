void plotBaby_PU(float lumi=19.5, int njets=-1, int mass=0, TString fs="sf", bool dodata=1, bool useSF=false, bool logy=0){
  //lumi is in /fb

  gROOT->Reset();
  gROOT->LoadMacro("tdrStyle.C");
  setTDRStyle();
  gStyle->SetOptStat(0);

  bool wjetsFromData=1;

  bool compareData = 0;

  bool doSignal = 0;

  bool doRatio = false;

  TString extension = ".root";

  TString mcs[] = {"qqww","ggww","dyll","ttbar_powheg","tw","wz","zz","wjets"};
  int  colors[] = {kAzure-9,kAzure-9,kGreen+2,kYellow,kYellow,kAzure-2,kAzure-2,kGray+1};

  TCut runrange("run>0");//Full2011
  TString dir = "/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/";
  float sfs0j[] = { 1.37,  1.37,   4.8,   0.9,   0.9,   1.0,    1.0, 1.0};
  float sfs1j[] = { 1.26,  1.26,   4.3,   0.9,   0.9,   1.0,    1.0, 1.0};
  float sfs2j[] = { 1.0,   1.0,    1.9,   1.5,   1.5,   1.0,    1.0, 1.0};

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
  unsigned int wwSelNoLep = BaseLine|ChargeMatch|TopVeto|ExtraLeptonVeto;//FullMET|ZVeto|
  
  TString mh = Form("%i",mass); 
  TString nj = Form("%i",njets); 

  TCut lep1pt,lep2pt,dPhi,mll,mt,himass;
  if (mass==0) {
    lep1pt = "lep1.pt()>20.";
    lep2pt = "lep2.pt()>20.";
    dPhi = "dPhi<TMath::Pi()*180./180.";
    mll = "dilep.mass()<999";
    mt = "mt>0&&mt<99999.";
    himass = "dilep.mass()>100.";
  } else if (mass==115) {
    lep1pt = "lep1.pt()>20.";
    lep2pt = "lep2.pt()>10.";
    dPhi = "dPhi<TMath::Pi()*115./180.";
    mll = "dilep.mass()<40";
    mt = "mt>80&&mt<110";
    himass = "dilep.mass()>100.";
  } else if (mass==120) {
    lep1pt = "lep1.pt()>20.";
    lep2pt = "lep2.pt()>10.";
    dPhi = "dPhi<TMath::Pi()*115./180.";
    mll = "dilep.mass()<40";
    mt = "mt>80&&mt<120";
    himass = "dilep.mass()>100.";
  } else if (mass==130) {
    lep1pt = "lep1.pt()>25.";
    lep2pt = "lep2.pt()>10.";
    dPhi = "dPhi<TMath::Pi()*90./180.";
    mll = "dilep.mass()<45";
    mt = "mt>80&&mt<125";
    himass = "dilep.mass()>100.";
  } else if (mass==140) {
    lep1pt = "lep1.pt()>25.";
    lep2pt = "lep2.pt()>15.";
    dPhi = "dPhi<TMath::Pi()*90./180.";
    mll = "dilep.mass()<45";
    mt = "mt>80&&mt<130";
    himass = "dilep.mass()>100.";
  } else if (mass==150) {
    lep1pt = "lep1.pt()>27.";
    lep2pt = "lep2.pt()>25.";
    dPhi = "dPhi<TMath::Pi()*90./180.";
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
  } else if (mass==180) {
    lep1pt = "lep1.pt()>36.";
    lep2pt = "lep2.pt()>25.";
    dPhi   = "dPhi<TMath::Pi()*70./180.";
    mll    = "dilep.mass()<60.";
    mt     = "mt>120&&mt<180";
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

  TCut base(Form("(cuts & %i)==%i",wwSelNoLep,wwSelNoLep));
  TCut leps(Form("(cuts & %i)==%i && (cuts & %i)==%i",Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection));
  TCut lplusf(Form("(((cuts & %i)==%i && (cuts & %i)!=%i)||((cuts & %i)!=%i && (cuts & %i)==%i))",
		   Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection,Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection));
  TCut notTagNotInJets(Form("(cuts & %i)!=%i",TopTagNotInJets,TopTagNotInJets));
  TCut trig(Form("dstype!=0 || (cuts & %i)==%i",Trigger,Trigger));

  TCut newcuts = "";
  TCut kincuts = "abs(dilep.mass()-91.2)<15.";

  TCut njcut(Form("njets==%i",njets));
  if (njets==-1) njcut = "";

  TCut flav = "";
  TString flavstr = "";
  TCut sf = "type!=1 && type!=2";
  //TCut sf = "type==0";
  TCut of = "type!=0 && type!=3";
  if (fs=="of") {
    flav = of;
    flavstr = "_of";
  } else if (fs=="sf") {
    flav = sf;
    flavstr = "_sf";
  } else if (fs!=""){
    cout << "final state not supported: " << fs << endl;
    return;
  }

  TCut cut = base&&njcut&&newcuts&&sigreg&&trig&&kincuts&&flav;
  cut.SetName("mh"+mh+"_nj"+nj+flavstr);

  cout << cut.GetTitle() << endl;

  TString plot[]    = {
//     "dilep.mass()",
//     "dPhi",
//     "dilep.pt()",
//     "lep1.pt()","lep2.pt()",
//     "type",
//     "pmet","pTrackMet",
//     "mt",
//     "jet1.pt()",
//     "jet1Btag",
//     "jet1ProbBtag",
    "nvtx"//,
//     "dymva"
  };
  TString binning[] = {
//     "30,0.,300.",
//     "20,0.,3.2",
//     "40,0.,200.",
//     "40,0.,200.","40,0.,200.",
//     "4,0,4",
//     "40,0.,200.","40,0.,200.",
//     "30,0.,300.",
//     "40,0.,200.",
//     "20,0.,4",
//     "20,0,2",
    "40,0,40",
    "21,-1,1.1"
  };
  TString xtitle[]  = {
//     "m_{l,l} [GeV/c^{2}]",
//     "#Delta#phi_{l,l} [rad]",
//     "pT_{l,l} [GeV/c]",
//     "pT_{l1} [GeV/c]","pT_{l2} [GeV/c]",
//     "type",
//     "projected #slash{E}_{T} [GeV/c^{2}]","projected track #slash{E}_{T} [GeV/c^{2}]",
//     "m_{T} [GeV/c^{2}]",
//     "pT_{j1} [GeV/c]",
//     "TCHE discriminator",
//     "JetProb. discriminator",
    "N_{vtx}",
    "dymva output"
  };

  //TString plot[]    = {"type"};
  //TString binning[] = {"4,0,4"};
  //TString xtitle[]  = {"type"};

  //if we want to plot signal
  TFile *_signal = 0;
  TTree * signal = 0;
  if (doSignal) {
    _signal = TFile::Open(dir+"/hww"+mh+".root");
    signal = (TTree*) _signal->Get("tree");
  }

  float nosfs[] = { 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0, 1.0, 1.0, 1.0,    1.0};
  float* sfs;
  if (useSF&&njets==0) sfs = sfs0j;
  else if (useSF&&njets==1) sfs = sfs1j;
  else if (useSF&&njets==2) sfs = sfs2j;
  else sfs = nosfs;

  TFile *_data = TFile::Open(dir+"/data.root");
  TTree * data = (TTree*) _data->Get("tree");

  bool allleg = false;

  int nMC = sizeof(mcs)/sizeof(TString);
  int nPlots = sizeof(plot)/sizeof(TString);
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

    TCanvas c;
    gPad->SetLeftMargin(0.125);
    THStack* hs = new THStack("hs",plot[pl]);
    TLegend* leg = new TLegend(0.15,0.84,0.48,0.94);
    leg->SetTextFont(42);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(3);
    leg->SetLineWidth(0);
    leg->SetLineColor(kWhite);
    leg->SetShadowColor(kWhite);

    if (!compareData) {
      //standard plots
      for (int i=0;i<nMC;++i){    
	TFile *_mc = TFile::Open(dir+"/"+mcs[i]+".root");
	TTree * mc = (TTree*) _mc->Get("tree");
	float correction = sfs[i];
	if (lumi>0.) mc->SetWeight(lumi*correction);
	//cout << mcs[i] << endl;
	TH1F* plotmc = new TH1F("plotmc_"+mcs[i],"plotmc",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	if (mcs[i].Contains("wjets")==0 || wjetsFromData==0) {
	  mc->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && %s)",cut.GetTitle(),leps.GetTitle()),"g");	
	  //cout << plotmc->Integral(0,plotmc->GetNbinsX()+1) << endl;
	} else {
	  //cout << plot[pl]+">>plotmc_"+mcs[i] << " " << Form("sfWeightFR*(%s && %s && %s)",cut.GetTitle(),lplusf.GetTitle(),runrange.GetTitle()) << endl;
	  data->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("sfWeightFR*(%s && %s && %s && sfWeightFR>-99.)",cut.GetTitle(),lplusf.GetTitle(),runrange.GetTitle()),"EP,same");
	  //cout << plotmc->Integral(0,plotmc->GetNbinsX()+1) << endl;
	}
	float integral = plotmc->Integral(0,plotmc->GetNbinsX()+1);//integral includes under/over flow bins
	//cout << integral << endl;
	plotmc->SetLineColor(colors[i]);
	plotmc->SetFillColor(colors[i]);
	hs->Add(plotmc);
	//if (allleg) leg->AddEntry(plotmc,TString(mcs[i]+": "+Form("%.2f",integral)),"f");
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
      hs->SetMaximum(hs->GetMaximum()*1.3);
      if (logy) {
	hs->SetMaximum(hs->GetMaximum()*1.3);
      }
      hs->Draw();

      TH1F* plotSignal = 0;
      if (doSignal) {
	//if we want to plot signal
	plotSignal = new TH1F("plotSignal","plotSignal",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	signal->Draw(plot[pl]+">>plotSignal",Form("%f*scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*sfWeightHPt*(%s && %s)",lumi,cut.GetTitle(),leps.GetTitle()),"same");
	//cout << plotSignal->Integral(0,plotSignal->GetNbinsX()+1) << endl;
	plotSignal->SetLineWidth(2);
	plotSignal->SetLineStyle(1);
	plotSignal->SetLineColor(kRed);
	leg->AddEntry(plotSignal,"HWW"+mh,"l");
      }

      TH1F* plotData = 0;
      if (dodata) {
	plotData = new TH1F("plotData","plotData",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	data->Draw(plot[pl]+">>plotData",cut&&leps&&runrange,"EP,same");
	c.Update();  
	//cout << "data" << endl;
	//float integral = plotData->Integral(0,plotmc->GetNbinsX()+1);//integral includes under/over flow bins
	//cout << integral << endl;
	plotData->SetMarkerStyle(20);
	int maxb=plotData->GetMaximumBin();
	double max=plotData->GetBinContent(maxb);
	if (max>(hs->GetMaximum())) {
	  hs->SetMaximum(max*1.3);
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
	  pad2->SetGridy();
	  ratio->SetMarkerStyle(21);
	  ratio->Draw("EP");
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
	hs->Draw();
	if (dodata) plotData->Draw("EP,same");
	if (doSignal) plotSignal->Draw("same");
	c.cd();
	
	leg->AddEntry(plotData,"data","p");
	
      }
      if (!allleg) {
	//TH1F h1("h1","h1",2,0,2);h1.SetFillColor(kBlack); h1.SetLineColor(kBlack); h1.SetLineWidth(2); h1.SetLineStyle(1); leg->AddEntry(&h1,"H("+mh+")#rightarrow WW x 10","f");
	TH1F h2("h2","h2",2,0,2);h2.SetFillColor(kAzure-9); leg->AddEntry(&h2,"WW","f");
	TH1F h3("h3","h3",2,0,2);h3.SetFillColor(kGreen+2); leg->AddEntry(&h3,"Z/#gamma*","f");
	TH1F h4("h4","h4",2,0,2);h4.SetFillColor(kYellow); leg->AddEntry(&h4,"top","f");
	TH1F h5("h5","h5",2,0,2);h5.SetFillColor(kAzure-2); leg->AddEntry(&h5,"VZ","f");
	TH1F h6("h6","h6",2,0,2);h6.SetFillColor(kGray+1); leg->AddEntry(&h6,"W+jets","f");
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
      labelcms  = new TPaveText(0.48,0.91,0.48,0.91,"NDCBR");
      labelcms->SetTextAlign(12);
      labelcms->SetTextSize(0.035);
      labelcms->SetFillColor(kWhite);
      labelcms->AddText(Form("CMS, #sqrt{s} = 8 TeV, L_{int} = %.2f fb^{-1}",lumi));
      labelcms->SetBorderSize(0);
      labelcms->SetTextFont(42);
      labelcms->SetLineWidth(2);
      labelcms->Draw();
    }

    c.Update();  

    TString cutName = cut.GetName();
    plot[pl].ReplaceAll(".","");
    plot[pl].ReplaceAll("(","");
    plot[pl].ReplaceAll(")","");
    gSystem->Exec("mkdir -p dirplots/mh"+mh+"_nj"+nj);
    if (logy) extension="_log"+extension;
    if (cutName.Length()==0) c.SaveAs("dirplots/mh"+mh+"_nj"+nj+"/"+plot[pl]+extension);
    else c.SaveAs("dirplots/mh"+mh+"_nj"+nj+"/"+plot[pl]+"_"+cutName+extension);
    if (!compareData) {
      delete plotmc;
      if (plotData) delete plotData; 
      delete hs;
    } else {
      delete plotData2011A;
      delete plotData2011B;
    }
    delete leg;
  }

}

