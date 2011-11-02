void plotBaby(float lumi=2.121, int njets=0, int mass=0, bool dodata=1, bool useSF=false, bool logy=0){
  //lumi is in /fb

  gROOT->Reset();
  gStyle->SetOptStat(0);

  bool wjetsFromData=1;

  bool compareData = 0;

  TCut runrange("run>0");//Full2011
  TString dir = "/smurf/cerati/skims/Run2011_Spring11_SmurfV7_42X/4ipbFull2011/wwSelNoLepNoTV/";
  float sfs0j[] = { 1.0,   1.0,   3.4,   3.4,   1.0,   1.4,    1.4, 1.0, 1.0, 2.8,    1.0};
  float sfs1j[] = { 1.0,   1.0,   4.2,   4.2,   1.0,   1.2,    1.2, 1.0, 1.0, 2.8,    1.0};
  //TString dir = "/smurf/cerati/skims/Run2011_Spring11_SmurfV7_42X/4ipbFull2011/wwSelNoMetNoZVminMET20/";

  //TCut runrange("run<=172802");//LP
  //float sfs0j[] = { 1.0,   1.0,   3.2,   3.2,   1.0,   1.0,    1.0, 1.0, 1.0, 2.7,    1.0};
  //float sfs1j[] = { 1.0,   1.0,   4.0,   4.0,   1.0,   1.1,    1.1, 1.0, 1.0, 3.2,    1.0};
  //TCut runrange("run<=173692");//2011A
  //float sfs0j[] = { 1.0,   1.0,   2.9,   2.9,   1.0,   0.9,    0.9, 1.0, 1.0, 2.7,    1.0};
  //float sfs1j[] = { 1.0,   1.0,   3.8,   3.8,   1.0,   1.1,    1.1, 1.0, 1.0, 3.3,    1.0};
  //TString dir = "/smurf/cerati/skims/Run2011_Spring11_SmurfV7_42X/4ipbRun2011A/wwSelNoLepNoTV/";

  //TCut runrange("run>173692");//2011B
  //float sfs0j[] = { 1.0,   1.0,   3.9,   3.9,   1.0,   2.1,    2.1, 1.0, 1.0, 3.2,    1.0};
  //float sfs1j[] = { 1.0,   1.0,   4.6,   4.6,   1.0,   1.3,    1.3, 1.0, 1.0, 2.4,    1.0};
  //TString dir = "/smurf/cerati/skims/Run2011_Spring11_SmurfV7_42X/4ipbRun2011B/wwSelNoLepNoTV/";

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
  unsigned int wwSelNoLep = BaseLine|ChargeMatch|FullMET|ZVeto|TopVeto|ExtraLeptonVeto;
  
  TString mh = Form("%i",mass); 
  TString nj = Form("%i",njets); 

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

  TCut base(Form("(cuts & %i)==%i",wwSelNoLep,wwSelNoLep));
  TCut leps(Form("(cuts & %i)==%i && (cuts & %i)==%i",Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection));
  TCut lplusf(Form("(((cuts & %i)==%i && (cuts & %i)!=%i)||((cuts & %i)!=%i && (cuts & %i)==%i))",
		   Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection,Lep1FullSelection,Lep1FullSelection,Lep2FullSelection,Lep2FullSelection));
  TCut notTagNotInJets(Form("(cuts & %i)!=%i",TopTagNotInJets,TopTagNotInJets));
  TCut trig(Form("dstype!=0 || (cuts & %i)==%i",Trigger,Trigger));
  TCut newcuts = "type==1 || type==2 || ( lep2.pt()>15. && min(pmet,pTrackMet)>(37.+nvtx/2.) && (jet1.pt()<15 || dPhiDiLepJet1*180./TMath::Pi()<165.) )";
  TCut kincuts = "dilep.pt()>45.";
  TCut njcut(Form("njets==%i",njets));
  if (njets==-1) njcut = "";

  TCut cut = base&&njcut&&newcuts&&sigreg&&trig&&kincuts;
  cut.SetName("mh"+mh+"_nj"+nj);

  TString plot[]    = {
    "dPhi",
    "dilep.mass()","dilep.pt()",
    "lep1.pt()","lep2.pt()",
    "type",
    "pmet","pTrackMet",
    "mt",
    "jet1.pt()"
  };
  TString binning[] = {
    "10,0.,3.2",
    "30,0.,300.","20,0.,200.",
    "20,0.,200.","20,0.,200.",
    "4,0,4",
    "20,0.,200.","20,0.,200.",
    "30,0.,300.",
    "75,0.,150."
  };
  TString xtitle[]  = {
    "#Delta#phi_{l,l}",
    "m_{l,l} [GeV]","pT_{l,l} [GeV]",
    "pT_{l1} [GeV]","pT_{l2} [GeV]",
    "type",
    "projected #slash{E}_{T} [GeV]","projected track #slash{E}_{T} [GeV]",
    "m_{T} [GeV]",
    "pT_{j1} [GeV]"
  };

  //TString plot[]    = {"type"};
  //TString binning[] = {"4,0,4"};
  //TString xtitle[]  = {"type"};
  //TString mcs[] = {"wz_py"};
  //int  colors[] = {kBlue};

  TString mcs[] = {"qqww","ggww","dyee","dymm","dytt","ttbar","tw","wz","zz","wjets","wgamma"};
  int  colors[] = {kCyan,kCyan,kGreen,kGreen,kGreen,kYellow,kYellow,kBlue,kBlue,kGray,kGray};
  float nosfs[] = { 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,    1.0, 1.0, 1.0, 1.0,    1.0};
  float* sfs;
  if (useSF&&njets==0) sfs = sfs0j;
  else if (useSF&&njets==1) sfs = sfs1j;
  else sfs = nosfs;

  ////if we want to plot signal
  //TFile *_signal = TFile::Open("/smurf/data/Run2011_Spring11_SmurfV3/tas-alljets/hww"+mh+".root");
  //TTree * signal = (TTree*) _signal->Get("tree");
  //float sig_correction = 10.;
  //if (lumi>0.) signal->SetWeight(lumi*sig_correction);

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
    THStack* hs = new THStack("hs",plot[pl]);
    TLegend* leg = new TLegend(0.1,0.9,0.9,1.);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(3);

    if (!compareData) {
      //standard plots
      for (int i=0;i<nMC;++i){    
	TFile *_mc = TFile::Open(dir+"/"+mcs[i]+".root");
	TTree * mc = (TTree*) _mc->Get("tree");
	float correction = sfs[i];
	if (lumi>0.) mc->SetWeight(lumi*correction);
	//cout << mcs[i] << endl;
	TH1F* plotmc = new TH1F("plotmc_"+mcs[i],"plotmc",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
	if (mcs[i]!="wjets" || wjetsFromData==0) {
	  mc->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s && %s)",cut.GetTitle(),leps.GetTitle()),"g");	
	} else {
	  data->Draw(plot[pl]+">>plotmc_"+mcs[i],Form("sfWeightFR*(%s && %s && %s)",cut.GetTitle(),lplusf.GetTitle(),runrange.GetTitle()),"EP,same");
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
      hs->GetYaxis()->SetTitle(Form("events/%.1f",hs->GetXaxis()->GetBinWidth(1)));
      hs->GetXaxis()->SetTitle(xtitle[pl]);
      hs->SetMaximum(hs->GetMaximum()*1.5.);
      if (logy) {
	hs->SetMaximum(hs->GetMaximum()*1.5);
      }
      hs->Draw();
      ////if we want to plot signal
      //signal->Draw(plot[pl]+">>plotSignal("+binning[pl]+")",Form("scale1fb*sfWeightTrig*sfWeightEff*sfWeightPU*(%s)",cut.GetTitle()),"same");
      //plotSignal->SetLineWidth(2);
      //plotSignal->SetLineStyle(1);
      TH1F* plotData = 0;
      if (dodata) {
	plotData = new TH1F("plotData","plotData",nbins.Atoi(),minbin.Atof(),maxbin.Atof());
      data->Draw(plot[pl]+">>plotData",cut&&leps&&runrange,"EP,same");//Full2011
      c.Update();  
      plotData->SetMarkerStyle(20);
      int maxb=plotData->GetMaximumBin();
      double max=plotData->GetBinContent(maxb);
      if (max>(hs->GetMaximum()*1.5)) {
	hs->SetMaximum(max*1.5);
	hs->Draw();
	plotData->Draw("EP,same");
      }
      leg->AddEntry(plotData,"data","p");
    }
    if (!allleg) {
      //TH1F h1("h1","h1",2,0,2);h1.SetFillColor(kBlack); h1.SetLineColor(kBlack); h1.SetLineWidth(2); h1.SetLineStyle(1); leg->AddEntry(&h1,"H("+mh+")#rightarrow WW x 10","f");
      TH1F h2("h2","h2",2,0,2);h2.SetFillColor(kCyan); leg->AddEntry(&h2,"WW","f");
      TH1F h3("h3","h3",2,0,2);h3.SetFillColor(kGreen); leg->AddEntry(&h3,"Drell-Yan","f");
      TH1F h4("h4","h4",2,0,2);h4.SetFillColor(kYellow); leg->AddEntry(&h4,"t#bar{t}, tW","f");
      TH1F h5("h5","h5",2,0,2);h5.SetFillColor(kBlue); leg->AddEntry(&h5,"di-boson","f");
      TH1F h6("h6","h6",2,0,2);h6.SetFillColor(kGray); leg->AddEntry(&h6,"W+jets","f");
    }
    labelcms  = new TPaveText(0.28,0.85,0.28,0.85,"NDCBR");
    labelcms->SetTextAlign(12);
    labelcms->SetTextSize(0.05);
    labelcms->SetFillColor(kWhite);
    labelcms->AddText(Form("CMS, #sqrt{s} = 7 TeV, L_{ int} = %.2f fb^{-1}",lumi));
    labelcms->SetBorderSize(0);
    labelcms->SetTextFont(132);
    labelcms->SetLineWidth(2);
    labelcms->Draw();
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
    //hs->SetMaximum(31);  

    c.Update();  

    TString cutName = cut.GetName();
    plot[pl].ReplaceAll(".","");
    plot[pl].ReplaceAll("(","");
    plot[pl].ReplaceAll(")","");
    gSystem->Exec("mkdir -p dirplots/mh"+mh+"_nj"+nj);
    TString extension = ".png";
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

