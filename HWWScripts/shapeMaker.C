#include "common.C"

void printBins(TH1F* h) {
  cout << h->GetName() << endl;
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    cout << bin << " " << h->GetBinContent(bin) << endl;
  }
}

void avoidNegativeBins(TH1F* h) {
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    if (h->GetBinContent(bin)<0) h->SetBinContent(bin,0);
  }
}

void fillDownMirrorUp(TH1F* central,TH1F* up,TH1F* down) {
  down->Add(up);
  down->Scale(-1);
  down->Add(central);
  down->Add(central);
  //need to avoid negative values...
  avoidNegativeBins(down);
}

void divideHisto(TH1F* num,TH1F* den){
  for (int bin=1;bin<=num->GetNbinsX();++bin) {
    if (fabs(den->GetBinContent(bin)>0) ) {
      //fixme: avoid large flututaions
      if (num->GetBinContent(bin)/den->GetBinContent(bin)<5.)
	num->SetBinContent(bin,num->GetBinContent(bin)/den->GetBinContent(bin));
      else {
	cout << "divideHisto - WARNING: setting bin to zero because dividing by a very small number - bin: " << bin <<  " - num " << num->GetName() << " " << num->GetBinContent(bin) 
	     << " den " << den->GetName() << " " << den->GetBinContent(bin) << endl;
	num->SetBinContent(bin,0);
      }
    } else num->SetBinContent(bin,0);
  }
}

void multiplyHisto(TH1F* num,TH1F* den){
  for (int bin=1;bin<=num->GetNbinsX();++bin) {
    num->SetBinContent(bin,num->GetBinContent(bin)*den->GetBinContent(bin));
  }
}

void scaleIntegral(TH1F* central,TH1F* other) {
  if (other->Integral()>0) other->Scale(central->Integral()/other->Integral());
}

void writeStatUpDown(TH1F* central,bool down, int njets, TString fs) {
  TString proc = TString(central->GetName());
  proc.ReplaceAll("histo_","");
  TString updown = "Up";
  if (down) updown = "Down";
  TH1F* statUpDown = new TH1F(TString(central->GetName())+"_CMS_MVA"+proc+Form("StatBounding_hww%s_%ij",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets)+updown,
			  TString(central->GetTitle())+"_CMS_MVA"+proc+Form("StatBounding_hww%s_%ij",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets)+updown,
			  central->GetNbinsX(),central->GetXaxis()->GetXmin(),central->GetXaxis()->GetXmax());
  for (int bin=1;bin<=statUpDown->GetNbinsX();++bin) {
    float val = down ? (central->GetBinContent(bin)-central->GetBinError(bin)) : (central->GetBinContent(bin)+central->GetBinError(bin));
    if (val>0) statUpDown->SetBinContent(bin,val);
    else statUpDown->SetBinContent(bin,1E-6);
  }
  statUpDown->Write();
}

void writeStatUpDown(TH1F* central, int njets, TString fs) {
  writeStatUpDown(central,true, njets, fs);
  writeStatUpDown(central,false, njets, fs);
}

void shapeMaker(float lumi=4.7, int njets=0, int mass=130, TString fs="sffs") {

  TString dir = main_dir+topww_dir;
  TString dirdy = main_dir+dy_dir;

  bool useJson = 0;
  bool applyEff=true;
  bool doFake=false; 
  bool doPUw=true;

  int nbins = 20;
  float minx = -1.;
  float maxx = 1.;

  TString sigreg = "=dphireg=dphijet=dymvacut=ptll45=";
  TString sigreg_lowmet = sigreg+"=zvetoall=";
  sigreg_lowmet.ReplaceAll("=dymvacut=","=loosedymva=");//fixme
  sigreg_lowmet.ReplaceAll("=dphijet=","=dpjallfs=");//fixme
  TString sigreg_himet = sigreg+"=zvetoall=";
  sigreg_himet.ReplaceAll("=dymvacut=","=dymvaallfs=");//fixme
  sigreg_himet.ReplaceAll("=dphijet=","=dpjallfs=");//fixme

  //cout << sigreg << endl;
  //cout << sigreg_lowmet << endl;

  std::vector<std::string> theInputVars;
  if (njets==0) {
    const char* inputVars[] = { "lep1pt", "lep2pt", "dPhi", "dR", "dilmass", "type", "mt" };
    for (int i=0;i<7;++i) theInputVars.push_back(inputVars[i]);
  } else if(njets==1) {
    const char* inputVars[] = { "lep1pt", "lep2pt", "dPhi", "dR", "dilmass", "type", "mt", "dPhiDiLepMET", "dPhiDiLepJet1" };
    for (int i=0;i<9;++i) theInputVars.push_back(inputVars[i]);
  } else {
    cout << "njets not supported: " << njets << endl;
    return;
  }
  rbdtg = new ReadBDTG(theInputVars);

  doMVA = 1;

  TString suffix = "" ;

  fs=fs+"=";

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  //Data
  TH1F* data_h = new TH1F("histo_Data","histo_Data",nbins,minx,maxx);
  fillPlot("bdtg",data_h, dir+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, 0, useJson, false, false, false);

  //qqWW
  TH1F* qqww_h = new TH1F("histo_qqWW","histo_qqWW",nbins,minx,maxx);
  fillPlot("bdtg",qqww_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  qqww_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  //shape variation: 1- mg vs mc@nlo (down mirror);
  TH1F* qqww_mcnlo_h = new TH1F("histo_qqww_mcnlo","histo_qqww_mcnlo",nbins,minx,maxx);
  fillPlot("bdtg",qqww_mcnlo_h, dir+"wwmcnlo"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  qqww_mcnlo_h->Scale(1.23);// fixme 7 TeV xsec
  avoidNegativeBins(qqww_mcnlo_h);
  scaleIntegral(qqww_h,qqww_mcnlo_h);
  TH1F* qqww_h_up = new TH1F("histo_qqWW_CMS_MVAWWBounding_hwwUp","histo_qqWW_CMS_MVAWWBounding_hwwUp",nbins,minx,maxx);
  qqww_h_up->Add(qqww_mcnlo_h);
  TH1F* qqww_h_down = new TH1F("histo_qqWW_CMS_MVAWWBounding_hwwDown","histo_qqWW_CMS_MVAWWBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(qqww_h,qqww_h_up,qqww_h_down);  
  //shape variation: 2- ratio from mc@nlo w.r.t. QCD up and down
  TH1F* qqww_mcnlo_up_h = new TH1F("histo_qqww_mcnlo_up","histo_qqww_mcnlo_up",nbins,minx,maxx);
  fillPlot("bdtg",qqww_mcnlo_up_h, dir+"wwmcnloup"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  qqww_mcnlo_up_h->Scale(1.23);// fixme 7 TeV xsec
  avoidNegativeBins(qqww_mcnlo_up_h);
  scaleIntegral(qqww_h,qqww_mcnlo_up_h);
  divideHisto(qqww_mcnlo_up_h,qqww_mcnlo_h);
  TH1F* qqww_mcnlo_down_h = new TH1F("histo_qqww_mcnlo_down","histo_qqww_mcnlo_down",nbins,minx,maxx);
  fillPlot("bdtg",qqww_mcnlo_down_h, dir+"wwmcnlodown"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  qqww_mcnlo_down_h->Scale(1.23);// fixme 7 TeV xsec
  avoidNegativeBins(qqww_mcnlo_down_h);
  scaleIntegral(qqww_h,qqww_mcnlo_down_h);
  divideHisto(qqww_mcnlo_down_h,qqww_mcnlo_h);
  TH1F* qqww_h_nlo_up = new TH1F("histo_qqWW_CMS_MVAWWNLOBounding_hwwUp","histo_qqWW_CMS_MVAWWNLOBounding_hwwUp",nbins,minx,maxx);
  qqww_h_nlo_up->Add(qqww_h);
  multiplyHisto(qqww_h_nlo_up,qqww_mcnlo_up_h);
  TH1F* qqww_h_nlo_down = new TH1F("histo_qqWW_CMS_MVAWWNLOBounding_hwwDown","histo_qqWW_CMS_MVAWWNLOBounding_hwwDown",nbins,minx,maxx);
  qqww_h_nlo_down->Add(qqww_h);
  multiplyHisto(qqww_h_nlo_down,qqww_mcnlo_down_h);

  //ggWW
  TH1F* ggww_h = new TH1F("histo_ggWW","histo_ggWW",nbins,minx,maxx);
  fillPlot("bdtg",ggww_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ggww_h->Scale(WWBkgScaleFactorMVA(mass,njets));

  //VV
  TH1F* wz_h = new TH1F("histo_wz","histo_wz",nbins,minx,maxx);
  fillPlot("bdtg",wz_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* zz_h = new TH1F("histo_zz","histo_zz",nbins,minx,maxx);
  fillPlot("bdtg",zz_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* vv_h = new TH1F("histo_VV","histo_VV",nbins,minx,maxx);
  vv_h->Add(wz_h);
  vv_h->Add(zz_h);

  //Top
  TH1F* ttbar_h = new TH1F("histo_ttbar","histo_ttbar",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ttbar_h->Scale(TopBkgScaleFactor(njets));
  TH1F* tw_h = new TH1F("histo_tw","histo_tw",nbins,minx,maxx);
  fillPlot("bdtg",tw_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  tw_h->Scale(TopBkgScaleFactor(njets));
  TH1F* top_h = new TH1F("histo_Top","histo_Top",nbins,minx,maxx);
  top_h->Add(ttbar_h);
  top_h->Add(tw_h);
  //shape variation: use madgraph ttbar and ds tw
  TH1F* ttbar_var_h = new TH1F("histo_ttbar_var","histo_ttbar_var",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_var_h, dir+"ttbar_mg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ttbar_var_h->Scale(1.43);//fixme scale 7 TeV xsec
  //scaleIntegral(ttbar_h,ttbar_var_h);
  TH1F* tw_ds_h = new TH1F("histo_tw_ds","histo_tw_ds",nbins,minx,maxx);
  fillPlot("bdtg",tw_ds_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);//fixme tw_ds
  //scaleIntegral(tw_h,tw_ds_h);
  TH1F* top_h_up = new TH1F("histo_Top_CMS_MVATopBounding_hwwUp","histo_Top_CMS_MVATopBounding_hwwUp",nbins,minx,maxx);
  top_h_up->Add(ttbar_var_h);
  top_h_up->Add(tw_ds_h);
  scaleIntegral(top_h,top_h_up);
  TH1F* top_h_down = new TH1F("histo_Top_CMS_MVATopBounding_hwwDown","histo_Top_CMS_MVATopBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(top_h,top_h_up,top_h_down);  

  //Wgamma
  TH1F* wg_h = new TH1F("histo_wg","histo_wg",nbins,minx,maxx);
  fillPlot("bdtg",wg_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* wg3l_h = new TH1F("histo_wg3l","histo_wg3l",nbins,minx,maxx);
  fillPlot("bdtg",wg3l_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  wg3l_h->Scale(WGstarScaleFactor());
  TH1F* wgamma_h = new TH1F("histo_Wgamma","histo_Wgamma",nbins,minx,maxx);
  wgamma_h->Add(wg_h);
  wgamma_h->Add(wg3l_h);

  //Zjets
  float dysf = 1.;
  TH1F* dyll_lowmet_h = new TH1F("histo_dyll_lowmet","histo_dyll_lowmet",nbins,minx,maxx);
  fillPlot("bdtg",dyll_lowmet_h, dirdy+"dyll"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+fs, lumi, useJson, applyEff, doFake, doPUw);
  float dyY = 0;
  if (fs.Contains("sffs")) {
    dyY = DYBkgScaleFactorBDT(mass,njets);
    dyll_lowmet_h->Scale(dyY/dyll_lowmet_h->Integral());
  } else {
    TH1F* dyll_vtx_h = new TH1F("histo_dyll_vtx","histo_dyll_vtx",nbins,minx,maxx);
    fillPlot("bdtg",dyll_vtx_h, dir+"dyll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
    dyll_vtx_h->Scale(dysf);
    scaleIntegral(dyll_vtx_h,dyll_lowmet_h);
  }
  TH1F* pwz_h = new TH1F("histo_pwz","histo_pwz",nbins,minx,maxx);
  fillPlot("bdtg",pwz_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=fromZ="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* pzz_h = new TH1F("histo_pzz","histo_pzz",nbins,minx,maxx);
  fillPlot("bdtg",pzz_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=fromZ="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* zjets_h = new TH1F("histo_Zjets","histo_Zjets",nbins,minx,maxx);
  zjets_h->Add(dyll_lowmet_h);
  zjets_h->Add(pwz_h);
  zjets_h->Add(pzz_h);
  //new shape variation: loose MET cuts in data
  TH1F* zjets_h_up = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUp",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUp",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    float kee = 0.8;//getK(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, 0., useJson, false, doFake, false);
    float lumicorr = (1.+kee*kee)/(2.*kee);
    TH1F* sffs_lowmet_new_h = new TH1F("histo_sffs_lowmet_new","histo_sffs_lowmet_new",nbins,minx,maxx);
    fillPlot("bdtg",sffs_lowmet_new_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", 0, useJson, false, false, false,"");//fixme add zeta method
    TH1F* offs_lowmet_new_h = new TH1F("histo_offs_lowmet_new","histo_offs_lowmet_new",nbins,minx,maxx);
    fillPlot("bdtg",offs_lowmet_new_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=offs=", 0, useJson, false, false, false,"");
    TH1F* wz_lowmet_new_h = new TH1F("histo_wz_lowmet_new","histo_wz_lowmet_new",nbins,minx,maxx);
    fillPlot("bdtg",wz_lowmet_new_h, dirdy+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    TH1F* zz_lowmet_new_h = new TH1F("histo_zz_lowmet_new","histo_zz_lowmet_new",nbins,minx,maxx);
    fillPlot("bdtg",zz_lowmet_new_h, dirdy+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    zjets_h_up->Add(sffs_lowmet_new_h);
    zjets_h_up->Add(offs_lowmet_new_h,-1.*lumicorr);
    zjets_h_up->Add(wz_lowmet_new_h,-1);
    zjets_h_up->Add(zz_lowmet_new_h,-1);
  }
  avoidNegativeBins(zjets_h_up);
  zjets_h_up->Scale(dyY/zjets_h_up->Integral());
  zjets_h_up->Add(pwz_h);
  zjets_h_up->Add(pzz_h);
  TH1F* zjets_h_down = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijDown",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijDown",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  fillDownMirrorUp(zjets_h,zjets_h_up,zjets_h_down);

  //old shape variation: use full MET cuts
  TH1F* zjets_h_old_up = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUpOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
                              Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUpOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  TH1F* zjets_h_old_down = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijDownOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
                                Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijDownOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    TH1F* dyll_vtx_h = new TH1F("histo_dyll_vtx","histo_dyll_vtx",nbins,minx,maxx);
    fillPlot("bdtg",dyll_vtx_h, dir+"dyll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
    zjets_h_old_up->Add(dyll_vtx_h);
    zjets_h_old_up->Add(pwz_h);
    zjets_h_old_up->Add(pzz_h);
    scaleIntegral(zjets_h,zjets_h_old_up);
    fillDownMirrorUp(zjets_h,zjets_h_old_up,zjets_h_old_down);
  }

  // OF,VZ subtraction
  TH1F* zjets_h_up_himet = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUpOFHiMet",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUpOFHiMet",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    float kee = 0.8;//getK(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, 0., useJson, false, doFake, false);
    float lumicorr = (1.+kee*kee)/(2.*kee);
    TH1F* sffs_himet_h = new TH1F("histo_sffs_himet","histo_sffs_himet",nbins,minx,maxx);
    fillPlot("bdtg",sffs_himet_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=sffs=", 0, useJson, false, false, false,"");
    TH1F* offs_himet_h = new TH1F("histo_offs_himet","histo_offs_himet",nbins,minx,maxx);
    fillPlot("bdtg",offs_himet_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=offs=", 0, useJson, false, false, false,"");
    TH1F* wz_himet_h = new TH1F("histo_wz_himet","histo_wz_himet",nbins,minx,maxx);
    fillPlot("bdtg",wz_himet_h, dirdy+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    TH1F* zz_himet_h = new TH1F("histo_zz_himet","histo_zz_himet",nbins,minx,maxx);
    fillPlot("bdtg",zz_himet_h, dirdy+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    zjets_h_up_himet->Add(sffs_himet_h);
    zjets_h_up_himet->Add(offs_himet_h,-1.*lumicorr);
    zjets_h_up_himet->Add(wz_himet_h,-1);
    zjets_h_up_himet->Add(zz_himet_h,-1);
  }
  avoidNegativeBins(zjets_h_up_himet);
  zjets_h_up_himet->Scale(dyY/zjets_h_up_himet->Integral());
  zjets_h_up_himet->Add(pwz_h);
  zjets_h_up_himet->Add(pzz_h);

  // histo_Zjets_CMS_MVAZBounding_hwwsf_0jUpZeta
  // zeta method
  TH1F* zjets_h_up_zeta = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUpZeta",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUpZeta",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    float kee = 0.8;//getK(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, 0., useJson, false, doFake, false);
    float lumicorr = (1.+kee*kee)/(2.*kee);
    TH1F* sffs_lowmet_zeta_h = new TH1F("histo_sffs_lowmet_zeta","histo_sffs_lowmet_zeta",nbins,minx,maxx);
    fillPlot("bdtg",sffs_lowmet_zeta_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", 0, useJson, false, false, false,"zeta");
    TH1F* offs_lowmet_zeta_h = new TH1F("histo_offs_lowmet_zeta","histo_offs_lowmet_zeta",nbins,minx,maxx);
    fillPlot("bdtg",offs_lowmet_zeta_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=offs=", 0, useJson, false, false, false,"zeta");
    TH1F* wz_lowmet_zeta_h = new TH1F("histo_wz_lowmet_zeta","histo_wz_lowmet_zeta",nbins,minx,maxx);
    fillPlot("bdtg",wz_lowmet_zeta_h, dirdy+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    TH1F* zz_lowmet_zeta_h = new TH1F("histo_zz_lowmet_zeta","histo_zz_lowmet_zeta",nbins,minx,maxx);
    fillPlot("bdtg",zz_lowmet_zeta_h, dirdy+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    zjets_h_up_zeta->Add(sffs_lowmet_zeta_h);
    zjets_h_up_zeta->Add(offs_lowmet_zeta_h,-1.*lumicorr);
    zjets_h_up_zeta->Add(wz_lowmet_zeta_h,-1);
    zjets_h_up_zeta->Add(zz_lowmet_zeta_h,-1);
  }
  avoidNegativeBins(zjets_h_up_zeta);
  zjets_h_up_zeta->Scale(dyY/zjets_h_up_zeta->Integral());
  zjets_h_up_zeta->Add(pwz_h);
  zjets_h_up_zeta->Add(pzz_h);

  //Ztt
  TH1F* dytt_1_h = new TH1F("histo_dytt_1","histo_dytt_1",nbins,minx,maxx);
  //fillPlot("bdtg",dytt_1_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false);
  TH1F* dytt_2_h = new TH1F("histo_dytt_2","histo_dytt_2",nbins,minx,maxx);
  //fillPlot("bdtg",dytt_2_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false);
  TH1F* dytt_3_h = new TH1F("histo_dytt_3","histo_dytt_3",nbins,minx,maxx);
  //fillPlot("bdtg",dytt_3_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false);
  TH1F* ztt_h = new TH1F("histo_Ztt","histo_Ztt",nbins,minx,maxx);
  ztt_h->Add(dytt_1_h);
  ztt_h->Add(dytt_2_h);
  ztt_h->Add(dytt_3_h);

  //Wjets
  TH1F* datafake_h = new TH1F("datafake","datafake",nbins,minx,maxx);
  fillPlot("bdtg",datafake_h,dir+"data"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, 0, useJson, false, true, false);
  TH1F* qqwwfake_h = new TH1F("qqwwfake","qqwwfake",nbins,minx,maxx);
  fillPlot("bdtg",qqwwfake_h,dir+"qqww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_h = new TH1F("ggwwfake","ggwwfake",nbins,minx,maxx);
  fillPlot("bdtg",ggwwfake_h,dir+"ggww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_h = new TH1F("ttbarfake","ttbarfake",nbins,minx,maxx);
  fillPlot("bdtg",ttbarfake_h,dir+"ttbar"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_h = new TH1F("twfake","twfake",nbins,minx,maxx);
  fillPlot("bdtg",twfake_h,dir+"tw"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_h = new TH1F("wzfake","wzfake",nbins,minx,maxx);
  fillPlot("bdtg",wzfake_h,dir+"wz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_h = new TH1F("zzfake","zzfake",nbins,minx,maxx);
  fillPlot("bdtg",zzfake_h,dir+"zz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wjets_h = new TH1F("histo_Wjets","histo_Wjets",nbins,minx,maxx);
  wjets_h->Add(datafake_h);
  wjets_h->Add(qqwwfake_h,-1.);
  wjets_h->Add(ggwwfake_h,-1.);
  wjets_h->Add(ttbarfake_h,-1.);
  wjets_h->Add(twfake_h,-1.);
  wjets_h->Add(wzfake_h,-1.);
  wjets_h->Add(zzfake_h,-1.);
  avoidNegativeBins(wjets_h);
  //syst 1: MC closure test
  TH1F* wjets_mc_up_h = new TH1F("histo_Wjets_CMS_MVAWMCBounding_hwwUp","histo_Wjets_CMS_MVAWMCBounding_hwwUp",nbins,minx,maxx);
  fillPlot("bdtg",wjets_mc_up_h,dir+"wjets"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  scaleIntegral(wjets_h,wjets_mc_up_h);
  TH1F* wjets_mc_down_h = new TH1F("histo_Wjets_CMS_MVAWMCBounding_hwwDown","histo_Wjets_CMS_MVAWMCBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(wjets_h,wjets_mc_up_h,wjets_mc_down_h);
  //syst 2: alternative fakebale object definition
  TH1F* datafake_fr_up_h = new TH1F("datafake_fr_up","datafake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",datafake_fr_up_h,dir+"data"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, 0, useJson, false, true, false);
  TH1F* qqwwfake_fr_up_h = new TH1F("qqwwfake_fr_up","qqwwfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",qqwwfake_fr_up_h,dir+"qqww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_fr_up_h = new TH1F("ggwwfake_fr_up","ggwwfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",ggwwfake_fr_up_h,dir+"ggww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_fr_up_h = new TH1F("ttbarfake_fr_up","ttbarfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",ttbarfake_fr_up_h,dir+"ttbar"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_fr_up_h = new TH1F("twfake_fr_up","twfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",twfake_fr_up_h,dir+"tw"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_fr_up_h = new TH1F("wzfake_fr_up","wzfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",wzfake_fr_up_h,dir+"wz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_fr_up_h = new TH1F("zzfake_fr_up","zzfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",zzfake_fr_up_h,dir+"zz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+"=alternativeFR="+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wjets_fr_up_h = new TH1F("histo_Wjets_CMS_MVAWBounding_hwwUp","histo_Wjets_CMS_MVAWBounding_hwwUp",nbins,minx,maxx);
  wjets_fr_up_h->Add(datafake_fr_up_h);
  wjets_fr_up_h->Add(qqwwfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(ggwwfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(ttbarfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(twfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(wzfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(zzfake_fr_up_h,-1.);
  avoidNegativeBins(wjets_fr_up_h);
  TH1F* wjets_fr_down_h = new TH1F("histo_Wjets_CMS_MVAWBounding_hwwDown","histo_Wjets_CMS_MVAWBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(wjets_h,wjets_fr_up_h,wjets_fr_down_h);

  //Higgs
  TH1F* ggH_h = new TH1F("histo_ggH","histo_ggH",nbins,minx,maxx);
  fillPlot("bdtg",ggH_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* qqH_h = new TH1F("histo_qqH","histo_qqH",nbins,minx,maxx);
  fillPlot("bdtg",qqH_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* WH_h = new TH1F("histo_WH","histo_WH",nbins,minx,maxx);
  fillPlot("bdtg",WH_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* ZH_h = new TH1F("histo_ZH","histo_ZH",nbins,minx,maxx);
  fillPlot("bdtg",ZH_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw);

  //ggH k-factor syst  
  TH1F* ggH_up_h = new TH1F("histo_ggH_CMS_MVAggHBoundingUp","histo_ggH_CMS_MVAggHBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ggH_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "ggH_k_syst_up");
  TH1F* ggH_down_h = new TH1F("histo_ggH_CMS_MVAggHBoundingDown","histo_ggH_CMS_MVAggHBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",ggH_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "ggH_k_syst_down");

  //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
  TH1F *ggH_metres_up_h=0, *ggH_metres_down_h=0, *qqH_metres_up_h=0, *qqH_metres_down_h=0, *WH_metres_up_h=0, *WH_metres_down_h=0, *ZH_metres_up_h=0, *ZH_metres_down_h=0, *qqww_metres_up_h=0, *qqww_metres_down_h=0, 
    *ggww_metres_up_h=0, *ggww_metres_down_h=0, *vv_metres_up_h=0, *vv_metres_down_h=0, *top_metres_up_h=0, *top_metres_down_h=0, *wgamma_metres_up_h=0, *wgamma_metres_down_h=0, *ztt_metres_up_h=0, *  ztt_metres_down_h=0;
  //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
  TH1F *ggH_lepres_up_h=0, *ggH_lepres_down_h=0, *qqH_lepres_up_h=0, *qqH_lepres_down_h=0, *WH_lepres_up_h=0, *WH_lepres_down_h=0, *ZH_lepres_up_h=0, *ZH_lepres_down_h=0, *qqww_lepres_up_h=0, *qqww_lepres_down_h=0, 
    *ggww_lepres_up_h=0, *ggww_lepres_down_h=0, *vv_lepres_up_h=0, *vv_lepres_down_h=0, *top_lepres_up_h=0, *top_lepres_down_h=0, *wgamma_lepres_up_h=0, *wgamma_lepres_down_h=0, *ztt_lepres_up_h=0, *ztt_lepres_down_h=0;
  //JES systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
  TH1F *ggH_jes_up_h=0, *ggH_jes_down_h=0, *qqH_jes_up_h=0, *qqH_jes_down_h=0, *WH_jes_up_h=0, *WH_jes_down_h=0, *ZH_jes_up_h=0, *ZH_jes_down_h=0, *qqww_jes_up_h=0, *qqww_jes_down_h=0, *ggww_jes_up_h=0, 
    *ggww_jes_down_h=0, *vv_jes_up_h=0, *vv_jes_down_h=0, *top_jes_up_h=0, *top_jes_down_h=0, *wgamma_jes_up_h=0, *wgamma_jes_down_h=0, *ztt_jes_up_h=0, *ztt_jes_down_h=0;
  //Lepton efficiency systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Wgamma and Ztt components (why no Top???)
  TH1F *ggH_lepeff_up_h=0, *ggH_lepeff_down_h=0, *qqH_lepeff_up_h=0, *qqH_lepeff_down_h=0, *WH_lepeff_up_h=0, *WH_lepeff_down_h=0, *ZH_lepeff_up_h=0, *ZH_lepeff_down_h=0, *qqww_lepeff_up_h=0, *qqww_lepeff_down_h=0, *ggww_lepeff_up_h=0, *ggww_lepeff_down_h=0, *vv_lepeff_up_h=0, *vv_lepeff_down_h=0, *wgamma_lepeff_up_h=0, *wgamma_lepeff_down_h=0, *ztt_lepeff_up_h=0, *ztt_lepeff_down_h=0;

  if (doResEffSyst) {
    //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
    //Higgs
    ggH_metres_up_h = new TH1F("histo_ggH_CMS_MVAMETResBoundingUp","histo_ggH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    scaleIntegral(ggH_h,ggH_metres_up_h);//fixme
    ggH_metres_down_h = new TH1F("histo_ggH_CMS_MVAMETResBoundingDown","histo_ggH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ggH_h,ggH_metres_up_h,ggH_metres_down_h);
    qqH_metres_up_h = new TH1F("histo_qqH_CMS_MVAMETResBoundingUp","histo_qqH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    qqH_metres_down_h = new TH1F("histo_qqH_CMS_MVAMETResBoundingDown","histo_qqH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(qqH_h,qqH_metres_up_h,qqH_metres_down_h);
    WH_metres_up_h = new TH1F("histo_WH_CMS_MVAMETResBoundingUp","histo_WH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",WH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    WH_metres_down_h = new TH1F("histo_WH_CMS_MVAMETResBoundingDown","histo_WH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(WH_h,WH_metres_up_h,WH_metres_down_h);
    ZH_metres_up_h = new TH1F("histo_ZH_CMS_MVAMETResBoundingUp","histo_ZH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ZH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    ZH_metres_down_h = new TH1F("histo_ZH_CMS_MVAMETResBoundingDown","histo_ZH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ZH_h,ZH_metres_up_h,ZH_metres_down_h);
    //WW
    qqww_metres_up_h = new TH1F("histo_qqWW_CMS_MVAMETResBoundingUp","histo_qqWW_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqww_metres_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    qqww_metres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    qqww_metres_down_h = new TH1F("histo_qqWW_CMS_MVAMETResBoundingDown","histo_qqWW_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(qqww_h,qqww_metres_up_h,qqww_metres_down_h);
    ggww_metres_up_h = new TH1F("histo_ggWW_CMS_MVAMETResBoundingUp","histo_ggWW_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggww_metres_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    ggww_metres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_metres_down_h = new TH1F("histo_ggWW_CMS_MVAMETResBoundingDown","histo_ggWW_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ggww_h,ggww_metres_up_h,ggww_metres_down_h);
    //VV
    TH1F* wz_metres_up_h = new TH1F("histo_wz_metres_up","histo_wz_metres_up",nbins,minx,maxx);
    fillPlot("bdtg",wz_metres_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    TH1F* zz_metres_up_h = new TH1F("histo_zz_metres_up","histo_zz_metres_up",nbins,minx,maxx);
    fillPlot("bdtg",zz_metres_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    vv_metres_up_h = new TH1F("histo_VV_CMS_MVAMETResBoundingUp","histo_VV_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    vv_metres_up_h->Add(wz_metres_up_h);
    vv_metres_up_h->Add(zz_metres_up_h);
    vv_metres_down_h = new TH1F("histo_VV_CMS_MVAMETResBoundingDown","histo_VV_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(vv_h,vv_metres_up_h,vv_metres_down_h);
    //Top
    TH1F* ttbar_metres_up_h = new TH1F("histo_ttbar_metres_up","histo_ttbar_metres_up",nbins,minx,maxx);
    fillPlot("bdtg",ttbar_metres_up_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    ttbar_metres_up_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_metres_up_h = new TH1F("histo_tw_metres_up","histo_tw_metres_up",nbins,minx,maxx);
    fillPlot("bdtg",tw_metres_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    tw_metres_up_h->Scale(TopBkgScaleFactor(njets));
    top_metres_up_h = new TH1F("histo_Top_CMS_MVAMETResBoundingUp","histo_Top_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    top_metres_up_h->Add(ttbar_metres_up_h);
    top_metres_up_h->Add(tw_metres_up_h);
    top_metres_down_h = new TH1F("histo_Top_CMS_MVAMETResBoundingDown","histo_Top_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(top_h,top_metres_up_h,top_metres_down_h);
    //Wgamma
    TH1F* wg_metres_up_h = new TH1F("histo_wg_metres_up","histo_wg_metres_up",nbins,minx,maxx);
    fillPlot("bdtg",wg_metres_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    TH1F* wg3l_metres_up_h = new TH1F("histo_wg3l_metres_up","histo_wg3l_metres_up",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_metres_up_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    wg3l_metres_up_h->Scale(WGstarScaleFactor());
    wgamma_metres_up_h = new TH1F("histo_Wgamma_CMS_MVAMETResBoundingUp","histo_Wgamma_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    wgamma_metres_up_h->Add(wg_metres_up_h);
    wgamma_metres_up_h->Add(wg3l_metres_up_h);
    wgamma_metres_down_h = new TH1F("histo_Wgamma_CMS_MVAMETResBoundingDown","histo_Wgamma_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(wgamma_h,wgamma_metres_up_h,wgamma_metres_down_h);
    //Ztt
    TH1F* dytt_1_metres_up_h = new TH1F("histo_dytt_metres_up_1","histo_dytt_metres_up_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_metres_up_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "metSmear");
    TH1F* dytt_2_metres_up_h = new TH1F("histo_dytt_metres_up_2","histo_dytt_metres_up_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_metres_up_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "metSmear");
    TH1F* dytt_3_metres_up_h = new TH1F("histo_dytt_metres_up_3","histo_dytt_metres_up_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_metres_up_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "metSmear");
    ztt_metres_up_h = new TH1F("histo_Ztt_CMS_MVAMETResBoundingUp","histo_Ztt_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
    ztt_metres_up_h->Add(dytt_1_metres_up_h);
    ztt_metres_up_h->Add(dytt_2_metres_up_h);
    ztt_metres_up_h->Add(dytt_3_metres_up_h);
    ztt_metres_down_h = new TH1F("histo_Ztt_CMS_MVAMETResBoundingDown","histo_Ztt_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ztt_h,ztt_metres_up_h,ztt_metres_down_h);
    //when WW and Top are from data need to normalize histo!
    if (mass<200){
      scaleIntegral(ggww_h,ggww_metres_down_h);
      scaleIntegral(ggww_h,ggww_metres_up_h);
      scaleIntegral(qqww_h,qqww_metres_down_h);
      scaleIntegral(qqww_h,qqww_metres_up_h);
    }
    scaleIntegral(top_h,top_metres_down_h);
    scaleIntegral(top_h,top_metres_up_h);

    //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
    ggH_lepres_up_h = new TH1F("histo_ggH_CMS_MVALepResBoundingUp","histo_ggH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ggH_lepres_down_h = new TH1F("histo_ggH_CMS_MVALepResBoundingDown","histo_ggH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ggH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    qqH_lepres_up_h = new TH1F("histo_qqH_CMS_MVALepResBoundingUp","histo_qqH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    qqH_lepres_down_h = new TH1F("histo_qqH_CMS_MVALepResBoundingDown","histo_qqH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",qqH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    WH_lepres_up_h = new TH1F("histo_WH_CMS_MVALepResBoundingUp","histo_WH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",WH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    WH_lepres_down_h = new TH1F("histo_WH_CMS_MVALepResBoundingDown","histo_WH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",WH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    ZH_lepres_up_h = new TH1F("histo_ZH_CMS_MVALepResBoundingUp","histo_ZH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ZH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ZH_lepres_down_h = new TH1F("histo_ZH_CMS_MVALepResBoundingDown","histo_ZH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ZH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    //WW
    qqww_lepres_up_h = new TH1F("histo_qqWW_CMS_MVALepResBoundingUp","histo_qqWW_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqww_lepres_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    qqww_lepres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_lepres_up_h = new TH1F("histo_ggWW_CMS_MVALepResBoundingUp","histo_ggWW_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggww_lepres_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ggww_lepres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    qqww_lepres_down_h = new TH1F("histo_qqWW_CMS_MVALepResBoundingDown","histo_qqWW_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",qqww_lepres_down_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    qqww_lepres_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_lepres_down_h = new TH1F("histo_ggWW_CMS_MVALepResBoundingDown","histo_ggWW_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ggww_lepres_down_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    ggww_lepres_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    //VV
    TH1F* wz_lepres_up_h = new TH1F("histo_wz_lepres_up","histo_wz_lepres_up",nbins,minx,maxx);
    fillPlot("bdtg",wz_lepres_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    TH1F* zz_lepres_up_h = new TH1F("histo_zz_lepres_up","histo_zz_lepres_up",nbins,minx,maxx);
    fillPlot("bdtg",zz_lepres_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    vv_lepres_up_h = new TH1F("histo_VV_CMS_MVALepResBoundingUp","histo_VV_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    vv_lepres_up_h->Add(wz_lepres_up_h);
    vv_lepres_up_h->Add(zz_lepres_up_h);
    TH1F* wz_lepres_down_h = new TH1F("histo_wz_lepres_down","histo_wz_lepres_down",nbins,minx,maxx);
    fillPlot("bdtg",wz_lepres_down_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    TH1F* zz_lepres_down_h = new TH1F("histo_zz_lepres_down","histo_zz_lepres_down",nbins,minx,maxx);
    fillPlot("bdtg",zz_lepres_down_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    vv_lepres_down_h = new TH1F("histo_VV_CMS_MVALepResBoundingDown","histo_VV_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    vv_lepres_down_h->Add(wz_lepres_down_h);
    vv_lepres_down_h->Add(zz_lepres_down_h);
    //Top
    TH1F* ttbar_lepres_up_h = new TH1F("histo_ttbar_lepres_up","histo_ttbar_lepres_up",nbins,minx,maxx);
    fillPlot("bdtg",ttbar_lepres_up_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ttbar_lepres_up_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_lepres_up_h = new TH1F("histo_tw_lepres_up","histo_tw_lepres_up",nbins,minx,maxx);
    fillPlot("bdtg",tw_lepres_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    tw_lepres_up_h->Scale(TopBkgScaleFactor(njets));
    top_lepres_up_h = new TH1F("histo_Top_CMS_MVALepResBoundingUp","histo_Top_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    top_lepres_up_h->Add(ttbar_lepres_up_h);
    top_lepres_up_h->Add(tw_lepres_up_h);
    TH1F* ttbar_lepres_down_h = new TH1F("histo_ttbar_lepres_down","histo_ttbar_lepres_down",nbins,minx,maxx);
    fillPlot("bdtg",ttbar_lepres_down_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    ttbar_lepres_down_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_lepres_down_h = new TH1F("histo_tw_lepres_down","histo_tw_lepres_down",nbins,minx,maxx);
    fillPlot("bdtg",tw_lepres_down_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    tw_lepres_down_h->Scale(TopBkgScaleFactor(njets));
    top_lepres_down_h = new TH1F("histo_Top_CMS_MVALepResBoundingDown","histo_Top_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    top_lepres_down_h->Add(ttbar_lepres_down_h);
    top_lepres_down_h->Add(tw_lepres_down_h);
    //Wgamma
    TH1F* wg_lepres_up_h = new TH1F("histo_wg_lepres_up","histo_wg_lepres_up",nbins,minx,maxx);
    fillPlot("bdtg",wg_lepres_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    TH1F* wg3l_lepres_up_h = new TH1F("histo_wg3l_lepres_up","histo_wg3l_lepres_up",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_lepres_up_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    wg3l_lepres_up_h->Scale(WGstarScaleFactor());
    wgamma_lepres_up_h = new TH1F("histo_Wgamma_CMS_MVALepResBoundingUp","histo_Wgamma_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    wgamma_lepres_up_h->Add(wg_lepres_up_h);
    wgamma_lepres_up_h->Add(wg3l_lepres_up_h);
    TH1F* wg_lepres_down_h = new TH1F("histo_wg_lepres_down","histo_wg_lepres_down",nbins,minx,maxx);
    fillPlot("bdtg",wg_lepres_down_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    TH1F* wg3l_lepres_down_h = new TH1F("histo_wg3l_lepres_down","histo_wg3l_lepres_down",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_lepres_down_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    wg3l_lepres_down_h->Scale(WGstarScaleFactor());
    wgamma_lepres_down_h = new TH1F("histo_Wgamma_CMS_MVALepResBoundingDown","histo_Wgamma_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    wgamma_lepres_down_h->Add(wg_lepres_down_h);
    wgamma_lepres_down_h->Add(wg3l_lepres_down_h);
    //Ztt
    TH1F* dytt_1_lepres_up_h = new TH1F("histo_dytt_lepres_up_1","histo_dytt_lepres_up_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_lepres_up_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "momScaleUp");
    TH1F* dytt_2_lepres_up_h = new TH1F("histo_dytt_lepres_up_2","histo_dytt_lepres_up_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_lepres_up_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "momScaleUp");
    TH1F* dytt_3_lepres_up_h = new TH1F("histo_dytt_lepres_up_3","histo_dytt_lepres_up_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_lepres_up_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "momScaleUp");
    ztt_lepres_up_h = new TH1F("histo_Ztt_CMS_MVALepResBoundingUp","histo_Ztt_CMS_MVALepResBoundingUp",nbins,minx,maxx);
    ztt_lepres_up_h->Add(dytt_1_lepres_up_h);
    ztt_lepres_up_h->Add(dytt_2_lepres_up_h);
    ztt_lepres_up_h->Add(dytt_3_lepres_up_h);
    TH1F* dytt_1_lepres_down_h = new TH1F("histo_dytt_lepres_down_1","histo_dytt_lepres_down_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_lepres_down_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "momScaleDown");
    TH1F* dytt_2_lepres_down_h = new TH1F("histo_dytt_lepres_down_2","histo_dytt_lepres_down_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_lepres_down_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "momScaleDown");
    TH1F* dytt_3_lepres_down_h = new TH1F("histo_dytt_lepres_down_3","histo_dytt_lepres_down_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_lepres_down_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "momScaleDown");
    ztt_lepres_down_h = new TH1F("histo_Ztt_CMS_MVALepResBoundingDown","histo_Ztt_CMS_MVALepResBoundingDown",nbins,minx,maxx);
    ztt_lepres_down_h->Add(dytt_1_lepres_down_h);
    ztt_lepres_down_h->Add(dytt_2_lepres_down_h);
    ztt_lepres_down_h->Add(dytt_3_lepres_down_h);
    //when WW and Top are from data need to normalize histo!
    if (mass<200){
      scaleIntegral(ggww_h,ggww_lepres_down_h);
      scaleIntegral(ggww_h,ggww_lepres_up_h);
      scaleIntegral(qqww_h,qqww_lepres_down_h);
      scaleIntegral(qqww_h,qqww_lepres_up_h);
    }
    scaleIntegral(top_h,top_lepres_down_h);
    scaleIntegral(top_h,top_lepres_up_h);

    //JES systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
    ggH_jes_up_h = new TH1F("histo_ggH_CMS_MVAJESBoundingUp","histo_ggH_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggH_jes_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ggH_jes_down_h = new TH1F("histo_ggH_CMS_MVAJESBoundingDown","histo_ggH_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ggH_jes_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    qqH_jes_up_h = new TH1F("histo_qqH_CMS_MVAJESBoundingUp","histo_qqH_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqH_jes_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    qqH_jes_down_h = new TH1F("histo_qqH_CMS_MVAJESBoundingDown","histo_qqH_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",qqH_jes_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    WH_jes_up_h = new TH1F("histo_WH_CMS_MVAJESBoundingUp","histo_WH_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",WH_jes_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    WH_jes_down_h = new TH1F("histo_WH_CMS_MVAJESBoundingDown","histo_WH_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",WH_jes_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    ZH_jes_up_h = new TH1F("histo_ZH_CMS_MVAJESBoundingUp","histo_ZH_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ZH_jes_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ZH_jes_down_h = new TH1F("histo_ZH_CMS_MVAJESBoundingDown","histo_ZH_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ZH_jes_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    //WW
    qqww_jes_up_h = new TH1F("histo_qqWW_CMS_MVAJESBoundingUp","histo_qqWW_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqww_jes_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    qqww_jes_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_jes_up_h = new TH1F("histo_ggWW_CMS_MVAJESBoundingUp","histo_ggWW_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggww_jes_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ggww_jes_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    qqww_jes_down_h = new TH1F("histo_qqWW_CMS_MVAJESBoundingDown","histo_qqWW_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",qqww_jes_down_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    qqww_jes_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_jes_down_h = new TH1F("histo_ggWW_CMS_MVAJESBoundingDown","histo_ggWW_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ggww_jes_down_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    ggww_jes_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    //VV
    TH1F* wz_jes_up_h = new TH1F("histo_wz_jes_up","histo_wz_jes_up",nbins,minx,maxx);
    fillPlot("bdtg",wz_jes_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    TH1F* zz_jes_up_h = new TH1F("histo_zz_jes_up","histo_zz_jes_up",nbins,minx,maxx);
    fillPlot("bdtg",zz_jes_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    vv_jes_up_h = new TH1F("histo_VV_CMS_MVAJESBoundingUp","histo_VV_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    vv_jes_up_h->Add(wz_jes_up_h);
    vv_jes_up_h->Add(zz_jes_up_h);
    TH1F* wz_jes_down_h = new TH1F("histo_wz_jes_down","histo_wz_jes_down",nbins,minx,maxx);
    fillPlot("bdtg",wz_jes_down_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ"+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    TH1F* zz_jes_down_h = new TH1F("histo_zz_jes_down","histo_zz_jes_down",nbins,minx,maxx);
    fillPlot("bdtg",zz_jes_down_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    vv_jes_down_h = new TH1F("histo_VV_CMS_MVAJESBoundingDown","histo_VV_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    vv_jes_down_h->Add(wz_jes_down_h);
    vv_jes_down_h->Add(zz_jes_down_h);
    //Top
    TH1F* ttbar_jes_up_h = new TH1F("histo_ttbar_jes_up","histo_ttbar_jes_up",nbins,minx,maxx);
    fillPlot("bdtg",ttbar_jes_up_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ttbar_jes_up_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_jes_up_h = new TH1F("histo_tw_jes_up","histo_tw_jes_up",nbins,minx,maxx);
    fillPlot("bdtg",tw_jes_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    tw_jes_up_h->Scale(TopBkgScaleFactor(njets));
    top_jes_up_h = new TH1F("histo_Top_CMS_MVAJESBoundingUp","histo_Top_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    top_jes_up_h->Add(ttbar_jes_up_h);
    top_jes_up_h->Add(tw_jes_up_h);
    TH1F* ttbar_jes_down_h = new TH1F("histo_ttbar_jes_down","histo_ttbar_jes_down",nbins,minx,maxx);
    fillPlot("bdtg",ttbar_jes_down_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    ttbar_jes_down_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_jes_down_h = new TH1F("histo_tw_jes_down","histo_tw_jes_down",nbins,minx,maxx);
    fillPlot("bdtg",tw_jes_down_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    tw_jes_down_h->Scale(TopBkgScaleFactor(njets));
    top_jes_down_h = new TH1F("histo_Top_CMS_MVAJESBoundingDown","histo_Top_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    top_jes_down_h->Add(ttbar_jes_down_h);
    top_jes_down_h->Add(tw_jes_down_h);
    //Wgamma
    TH1F* wg_jes_up_h = new TH1F("histo_wg_jes_up","histo_wg_jes_up",nbins,minx,maxx);
    fillPlot("bdtg",wg_jes_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    TH1F* wg3l_jes_up_h = new TH1F("histo_wg3l_jes_up","histo_wg3l_jes_up",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_jes_up_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    wg3l_jes_up_h->Scale(WGstarScaleFactor());
    wgamma_jes_up_h = new TH1F("histo_Wgamma_CMS_MVAJESBoundingUp","histo_Wgamma_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    wgamma_jes_up_h->Add(wg_jes_up_h);
    wgamma_jes_up_h->Add(wg3l_jes_up_h);
    TH1F* wg_jes_down_h = new TH1F("histo_wg_jes_down","histo_wg_jes_down",nbins,minx,maxx);
    fillPlot("bdtg",wg_jes_down_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    TH1F* wg3l_jes_down_h = new TH1F("histo_wglll_jes_down","histo_wg3l_jes_down",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_jes_down_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    wg3l_jes_down_h->Scale(WGstarScaleFactor());
    wgamma_jes_down_h = new TH1F("histo_Wgamma_CMS_MVAJESBoundingDown","histo_Wgamma_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    wgamma_jes_down_h->Add(wg_jes_down_h);
    wgamma_jes_down_h->Add(wg3l_jes_down_h);
    //Ztt
    TH1F* dytt_1_jes_up_h = new TH1F("histo_dytt_jes_up_1","histo_dytt_jes_up_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_jes_up_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "jesUp");
    TH1F* dytt_2_jes_up_h = new TH1F("histo_dytt_jes_up_2","histo_dytt_jes_up_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_jes_up_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "jesUp");
    TH1F* dytt_3_jes_up_h = new TH1F("histo_dytt_jes_up_3","histo_dytt_jes_up_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_jes_up_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "jesUp");
    ztt_jes_up_h = new TH1F("histo_Ztt_CMS_MVAJESBoundingUp","histo_Ztt_CMS_MVAJESBoundingUp",nbins,minx,maxx);
    ztt_jes_up_h->Add(dytt_1_jes_up_h);
    ztt_jes_up_h->Add(dytt_2_jes_up_h);
    ztt_jes_up_h->Add(dytt_3_jes_up_h);
    TH1F* dytt_1_jes_down_h = new TH1F("histo_dytt_jes_down_1","histo_dytt_jes_down_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_jes_down_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "jesDown");
    TH1F* dytt_2_jes_down_h = new TH1F("histo_dytt_jes_down_2","histo_dytt_jes_down_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_jes_down_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "jesDown");
    TH1F* dytt_3_jes_down_h = new TH1F("histo_dytt_jes_down_3","histo_dytt_jes_down_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_jes_down_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "jesDown");
    ztt_jes_down_h = new TH1F("histo_Ztt_CMS_MVAJESBoundingDown","histo_Ztt_CMS_MVAJESBoundingDown",nbins,minx,maxx);
    ztt_jes_down_h->Add(dytt_1_jes_down_h);
    ztt_jes_down_h->Add(dytt_2_jes_down_h);
    ztt_jes_down_h->Add(dytt_3_jes_down_h);
    //when WW and Top are from data need to normalize histo!
    if (mass<200){
      scaleIntegral(ggww_h,ggww_jes_down_h);
      scaleIntegral(ggww_h,ggww_jes_up_h);
      scaleIntegral(qqww_h,qqww_jes_down_h);
      scaleIntegral(qqww_h,qqww_jes_up_h);
    }
    scaleIntegral(top_h,top_jes_down_h);
    scaleIntegral(top_h,top_jes_up_h);

    //Lepton efficiency systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Wgamma and Ztt components (no Top???)
    ggH_lepeff_up_h = new TH1F("histo_ggH_CMS_MVALepEffBoundingUp","histo_ggH_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggH_lepeff_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    ggH_lepeff_down_h = new TH1F("histo_ggH_CMS_MVALepEffBoundingDown","histo_ggH_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ggH_lepeff_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    qqH_lepeff_up_h = new TH1F("histo_qqH_CMS_MVALepEffBoundingUp","histo_qqH_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqH_lepeff_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    qqH_lepeff_down_h = new TH1F("histo_qqH_CMS_MVALepEffBoundingDown","histo_qqH_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",qqH_lepeff_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    WH_lepeff_up_h = new TH1F("histo_WH_CMS_MVALepEffBoundingUp","histo_WH_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",WH_lepeff_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    WH_lepeff_down_h = new TH1F("histo_WH_CMS_MVALepEffBoundingDown","histo_WH_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",WH_lepeff_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    ZH_lepeff_up_h = new TH1F("histo_ZH_CMS_MVALepEffBoundingUp","histo_ZH_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ZH_lepeff_up_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    ZH_lepeff_down_h = new TH1F("histo_ZH_CMS_MVALepEffBoundingDown","histo_ZH_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ZH_lepeff_down_h, dir+Form("hww%i",mass)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //WW
    qqww_lepeff_up_h = new TH1F("histo_qqWW_CMS_MVALepEffBoundingUp","histo_qqWW_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",qqww_lepeff_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    qqww_lepeff_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_lepeff_up_h = new TH1F("histo_ggWW_CMS_MVALepEffBoundingUp","histo_ggWW_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot("bdtg",ggww_lepeff_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    ggww_lepeff_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    qqww_lepeff_down_h = new TH1F("histo_qqWW_CMS_MVALepEffBoundingDown","histo_qqWW_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",qqww_lepeff_down_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    qqww_lepeff_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    ggww_lepeff_down_h = new TH1F("histo_ggWW_CMS_MVALepEffBoundingDown","histo_ggWW_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot("bdtg",ggww_lepeff_down_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    ggww_lepeff_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
    //VV
    TH1F* wz_lepeff_up_h = new TH1F("histo_wz_lepeff_up","histo_wz_lepeff_up",nbins,minx,maxx);
    fillPlot("bdtg",wz_lepeff_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    TH1F* zz_lepeff_up_h = new TH1F("histo_zz_lepeff_up","histo_zz_lepeff_up",nbins,minx,maxx);
    fillPlot("bdtg",zz_lepeff_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    vv_lepeff_up_h = new TH1F("histo_VV_CMS_MVALepEffBoundingUp","histo_VV_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    vv_lepeff_up_h->Add(wz_lepeff_up_h);
    vv_lepeff_up_h->Add(zz_lepeff_up_h);
    TH1F* wz_lepeff_down_h = new TH1F("histo_wz_lepeff_down","histo_wz_lepeff_down",nbins,minx,maxx);
    fillPlot("bdtg",wz_lepeff_down_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    TH1F* zz_lepeff_down_h = new TH1F("histo_zz_lepeff_down","histo_zz_lepeff_down",nbins,minx,maxx);
    fillPlot("bdtg",zz_lepeff_down_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=notZ="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    vv_lepeff_down_h = new TH1F("histo_VV_CMS_MVALepEffBoundingDown","histo_VV_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    vv_lepeff_down_h->Add(wz_lepeff_down_h);
    vv_lepeff_down_h->Add(zz_lepeff_down_h);
    //Top
    //TH1F* ttbar_lepeff_up_h = new TH1F("histo_ttbar_lepeff_up","histo_ttbar_lepeff_up",nbins,minx,maxx);
    //fillPlot("bdtg",ttbar_lepeff_up_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    //ttbar_lepeff_up_h->Scale(TopBkgScaleFactor(njets));
    //TH1F* tw_lepeff_up_h = new TH1F("histo_tw_lepeff_up","histo_tw_lepeff_up",nbins,minx,maxx);
    //fillPlot("bdtg",tw_lepeff_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    //tw_lepeff_up_h->Scale(TopBkgScaleFactor(njets));
    //top_lepeff_up_h = new TH1F("histo_Top_CMS_MVALepEffBoundingUp","histo_Top_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    //top_lepeff_up_h->Add(ttbar_lepeff_up_h);
    //top_lepeff_up_h->Add(tw_lepeff_up_h);
    //TH1F* ttbar_lepeff_down_h = new TH1F("histo_ttbar_lepeff_down","histo_ttbar_lepeff_down",nbins,minx,maxx);
    //fillPlot("bdtg",ttbar_lepeff_down_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //ttbar_lepeff_down_h->Scale(TopBkgScaleFactor(njets));
    //TH1F* tw_lepeff_down_h = new TH1F("histo_tw_lepeff_down","histo_tw_lepeff_down",nbins,minx,maxx);
    //fillPlot("bdtg",tw_lepeff_down_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //tw_lepeff_down_h->Scale(TopBkgScaleFactor(njets));
    //top_lepeff_down_h = new TH1F("histo_Top_CMS_MVALepEffBoundingDown","histo_Top_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    //top_lepeff_down_h->Add(ttbar_lepeff_down_h);
    //top_lepeff_down_h->Add(tw_lepeff_down_h);
    //Wgamma
    TH1F* wg_lepeff_up_h = new TH1F("histo_wg_lepeff_up","histo_wg_lepeff_up",nbins,minx,maxx);
    fillPlot("bdtg",wg_lepeff_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    TH1F* wg3l_lepeff_up_h = new TH1F("histo_wg3l_lepeff_up","histo_wg3l_lepeff_up",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_lepeff_up_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    wg3l_lepeff_up_h->Scale(WGstarScaleFactor());
    wgamma_lepeff_up_h = new TH1F("histo_Wgamma_CMS_MVALepEffBoundingUp","histo_Wgamma_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    wgamma_lepeff_up_h->Add(wg_lepeff_up_h);
    wgamma_lepeff_up_h->Add(wg3l_lepeff_up_h);
    TH1F* wg_lepeff_down_h = new TH1F("histo_wg_lepeff_down","histo_wg_lepeff_down",nbins,minx,maxx);
    fillPlot("bdtg",wg_lepeff_down_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    TH1F* wg3l_lepeff_down_h = new TH1F("histo_wg3l_lepeff_down","histo_wg3l_lepeff_down",nbins,minx,maxx);
    fillPlot("bdtg",wg3l_lepeff_down_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    wg3l_lepeff_down_h->Scale(WGstarScaleFactor());
    wgamma_lepeff_down_h = new TH1F("histo_Wgamma_CMS_MVALepEffBoundingDown","histo_Wgamma_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    wgamma_lepeff_down_h->Add(wg_lepeff_down_h);
    wgamma_lepeff_down_h->Add(wg3l_lepeff_down_h);
    //Ztt
    TH1F* dytt_1_lepeff_up_h = new TH1F("histo_dytt_lepeff_up_1","histo_dytt_lepeff_up_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_lepeff_up_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "lepeffUp");
    TH1F* dytt_2_lepeff_up_h = new TH1F("histo_dytt_lepeff_up_2","histo_dytt_lepeff_up_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_lepeff_up_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "lepeffUp");
    TH1F* dytt_3_lepeff_up_h = new TH1F("histo_dytt_lepeff_up_3","histo_dytt_lepeff_up_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_lepeff_up_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "lepeffUp");
    ztt_lepeff_up_h = new TH1F("histo_Ztt_CMS_MVALepEffBoundingUp","histo_Ztt_CMS_MVALepEffBoundingUp",nbins,minx,maxx);
    ztt_lepeff_up_h->Add(dytt_1_lepeff_up_h);
    ztt_lepeff_up_h->Add(dytt_2_lepeff_up_h);
    ztt_lepeff_up_h->Add(dytt_3_lepeff_up_h);
    TH1F* dytt_1_lepeff_down_h = new TH1F("histo_dytt_lepeff_down_1","histo_dytt_lepeff_down_1",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_1_lepeff_down_h, dir+"data-emb-tau121"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "lepeffDown");
    TH1F* dytt_2_lepeff_down_h = new TH1F("histo_dytt_lepeff_down_2","histo_dytt_lepeff_down_2",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_2_lepeff_down_h, dir+"data-emb-tau122"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "lepeffDown");
    TH1F* dytt_3_lepeff_down_h = new TH1F("histo_dytt_lepeff_down_3","histo_dytt_lepeff_down_3",nbins,minx,maxx);
    //fillPlot("bdtg",dytt_3_lepeff_down_h, dir+"data-emb-tau123"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false, "lepeffDown");
    ztt_lepeff_down_h = new TH1F("histo_Ztt_CMS_MVALepEffBoundingDown","histo_Ztt_CMS_MVALepEffBoundingDown",nbins,minx,maxx);
    ztt_lepeff_down_h->Add(dytt_1_lepeff_down_h);
    ztt_lepeff_down_h->Add(dytt_2_lepeff_down_h);
    ztt_lepeff_down_h->Add(dytt_3_lepeff_down_h);
    //when WW and Top are from data need to normalize histo!
    if (mass<200){
      scaleIntegral(ggww_h,ggww_lepeff_down_h);
      scaleIntegral(ggww_h,ggww_lepeff_up_h);
      scaleIntegral(qqww_h,qqww_lepeff_down_h);
      scaleIntegral(qqww_h,qqww_lepeff_up_h);
    }
    //scaleIntegral(top_h,top_lepeff_down_h);
    //scaleIntegral(top_h,top_lepeff_up_h);

  }

//   printBins(data_h);
//   printBins(ggww_h);
//   printBins(qqww_h);
//   printBins(top_h);
//   printBins(zjets_h);
//   printBins(wjets_h);
//   printBins(vv_h);
//   printBins(wgamma_h);
//   printBins(ztt_h);
//   printBins(qqH_h);
//   printBins(ggH_h);
//   printBins(WH_h);
//   printBins(ZH_h);

//   printBins(ggH_up_h);
//   printBins(ggH_down_h);

//   printBins(ggH_metres_up_h);
//   printBins(ggH_metres_down_h);
//   printBins(qqH_metres_up_h);
//   printBins(qqH_metres_down_h);
//   printBins(WH_metres_up_h);
//   printBins(WH_metres_down_h);
//   printBins(ZH_metres_up_h);
//   printBins(ZH_metres_down_h);
//   printBins(qqww_metres_up_h);
//   printBins(qqww_metres_down_h);
//   printBins(ggww_metres_up_h);
//   printBins(ggww_metres_down_h);
//   printBins(vv_metres_up_h);
//   printBins(vv_metres_down_h);
//   printBins(top_metres_up_h);
//   printBins(top_metres_down_h);
//   printBins(wgamma_metres_up_h);
//   printBins(wgamma_metres_down_h);
//   printBins(ztt_metres_up_h);
//   printBins(ztt_metres_down_h);

//   printBins(ggH_lepres_up_h);
//   printBins(ggH_lepres_down_h);
//   printBins(qqH_lepres_up_h);
//   printBins(qqH_lepres_down_h);
//   printBins(WH_lepres_up_h);
//   printBins(WH_lepres_down_h);
//   printBins(ZH_lepres_up_h);
//   printBins(ZH_lepres_down_h);
//   printBins(qqww_lepres_up_h);
//   printBins(qqww_lepres_down_h);
//   printBins(ggww_lepres_up_h);
//   printBins(ggww_lepres_down_h);
//   printBins(vv_lepres_up_h);
//   printBins(vv_lepres_down_h);
//   printBins(top_lepres_up_h);
//   printBins(top_lepres_down_h);
//   printBins(wgamma_lepres_up_h);
//   printBins(wgamma_lepres_down_h);
//   printBins(ztt_lepres_up_h);
//   printBins(ztt_lepres_down_h);

//   printBins(ggH_jes_up_h);
//   printBins(ggH_jes_down_h);
//   printBins(qqH_jes_up_h);
//   printBins(qqH_jes_down_h);
//   printBins(WH_jes_up_h);
//   printBins(WH_jes_down_h);
//   printBins(ZH_jes_up_h);
//   printBins(ZH_jes_down_h);
//   printBins(qqww_jes_up_h);
//   printBins(qqww_jes_down_h);
//   printBins(ggww_jes_up_h);
//   printBins(ggww_jes_down_h);
//   printBins(vv_jes_up_h);
//   printBins(vv_jes_down_h);
//   printBins(top_jes_up_h);
//   printBins(top_jes_down_h);
//   printBins(wgamma_jes_up_h);
//   printBins(wgamma_jes_down_h);
//   printBins(ztt_jes_up_h);
//   printBins(ztt_jes_down_h);

//   printBins(ggH_lepeff_up_h);
//   printBins(ggH_lepeff_down_h);
//   printBins(qqH_lepeff_up_h);
//   printBins(qqH_lepeff_down_h);
//   printBins(WH_lepeff_up_h);
//   printBins(WH_lepeff_down_h);
//   printBins(ZH_lepeff_up_h);
//   printBins(ZH_lepeff_down_h);
//   printBins(qqww_lepeff_up_h);
//   printBins(qqww_lepeff_down_h);
//   printBins(ggww_lepeff_up_h);
//   printBins(ggww_lepeff_down_h);
//   printBins(vv_lepeff_up_h);
//   printBins(vv_lepeff_down_h);
//   printBins(wgamma_lepeff_up_h);
//   printBins(wgamma_lepeff_down_h);
//   printBins(ztt_lepeff_up_h);
//   printBins(ztt_lepeff_down_h);

//   printBins(qqww_h_up);
//   printBins(qqww_h_down);
//   printBins(qqww_h_nlo_up);
//   printBins(qqww_h_nlo_down);

//   printBins(top_h_up);
//   printBins(top_h_down);

//   printBins(wjets_fr_up_h);
//   printBins(wjets_fr_down_h);
//   printBins(wjets_mc_up_h);
//   printBins(wjets_mc_down_h);

 //Write histos to file
  TString outfname = Form("hww%s_%ij.input.root",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets);
  TFile* outfile = TFile::Open(outfname,"RECREATE");

  //Nominal Histos
  data_h->Write();
  ggww_h->Write();
  qqww_h->Write();
  top_h->Write();
  zjets_h->Write();
  wjets_h->Write();
  vv_h->Write();
  wgamma_h->Write();
  ztt_h->Write();
  qqH_h->Write();
  ggH_h->Write();
  WH_h->Write();
  ZH_h->Write();

  //ggH k-factor syst  
  ggH_up_h->Write();
  ggH_down_h->Write();

  //Stat uncertainty
  writeStatUpDown(ggww_h,njets,fs);
  writeStatUpDown(qqww_h,njets,fs);
  writeStatUpDown(top_h,njets,fs);
  writeStatUpDown(zjets_h,njets,fs);
  writeStatUpDown(wjets_h,njets,fs);
  writeStatUpDown(vv_h,njets,fs);
  writeStatUpDown(wgamma_h,njets,fs);
  writeStatUpDown(ztt_h,njets,fs);
  writeStatUpDown(qqH_h,njets,fs);
  writeStatUpDown(ggH_h,njets,fs);
  writeStatUpDown(WH_h,njets,fs);
  writeStatUpDown(ZH_h,njets,fs);

  if (doResEffSyst) {
    //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
    ggH_metres_up_h->Write();
    ggH_metres_down_h->Write();
    qqH_metres_up_h->Write();
    qqH_metres_down_h->Write();
    WH_metres_up_h->Write();
    WH_metres_down_h->Write();
    ZH_metres_up_h->Write();
    ZH_metres_down_h->Write();
    qqww_metres_up_h->Write();
    qqww_metres_down_h->Write();
    ggww_metres_up_h->Write();
    ggww_metres_down_h->Write();
    vv_metres_up_h->Write();
    vv_metres_down_h->Write();
    top_metres_up_h->Write();
    top_metres_down_h->Write();
    wgamma_metres_up_h->Write();
    wgamma_metres_down_h->Write();
    ztt_metres_up_h->Write();
    ztt_metres_down_h->Write();

    //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
    ggH_lepres_up_h->Write();
    ggH_lepres_down_h->Write();
    qqH_lepres_up_h->Write();
    qqH_lepres_down_h->Write();
    WH_lepres_up_h->Write();
    WH_lepres_down_h->Write();
    ZH_lepres_up_h->Write();
    ZH_lepres_down_h->Write();
    qqww_lepres_up_h->Write();
    qqww_lepres_down_h->Write();
    ggww_lepres_up_h->Write();
    ggww_lepres_down_h->Write();
    vv_lepres_up_h->Write();
    vv_lepres_down_h->Write();
    top_lepres_up_h->Write();
    top_lepres_down_h->Write();
    wgamma_lepres_up_h->Write();
    wgamma_lepres_down_h->Write();
    ztt_lepres_up_h->Write();
    ztt_lepres_down_h->Write();

    //JES systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
    ggH_jes_up_h->Write();
    ggH_jes_down_h->Write();
    qqH_jes_up_h->Write();
    qqH_jes_down_h->Write();
    WH_jes_up_h->Write();
    WH_jes_down_h->Write();
    ZH_jes_up_h->Write();
    ZH_jes_down_h->Write();
    qqww_jes_up_h->Write();
    qqww_jes_down_h->Write();
    ggww_jes_up_h->Write();
    ggww_jes_down_h->Write();
    vv_jes_up_h->Write();
    vv_jes_down_h->Write();
    top_jes_up_h->Write();
    top_jes_down_h->Write();
    wgamma_jes_up_h->Write();
    wgamma_jes_down_h->Write();
    ztt_jes_up_h->Write();
    ztt_jes_down_h->Write();

    //Lepton efficiency systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Wgamma and Ztt components (why no Top???)
    ggH_lepeff_up_h->Write();
    ggH_lepeff_down_h->Write();
    qqH_lepeff_up_h->Write();
    qqH_lepeff_down_h->Write();
    WH_lepeff_up_h->Write();
    WH_lepeff_down_h->Write();
    ZH_lepeff_up_h->Write();
    ZH_lepeff_down_h->Write();
    qqww_lepeff_up_h->Write();
    qqww_lepeff_down_h->Write();
    ggww_lepeff_up_h->Write();
    ggww_lepeff_down_h->Write();
    vv_lepeff_up_h->Write();
    vv_lepeff_down_h->Write();
    //top_lepeff_up_h->Write();
    //top_lepeff_down_h->Write();
    wgamma_lepeff_up_h->Write();
    wgamma_lepeff_down_h->Write();
    ztt_lepeff_up_h->Write();
    ztt_lepeff_down_h->Write();
  }

  //Other systematics
  qqww_h_up->Write();
  qqww_h_down->Write();
  qqww_h_nlo_up->Write();
  qqww_h_nlo_down->Write();

  top_h_up->Write();
  top_h_down->Write();

  if (fs.Contains("sffs")) {
    zjets_h_up->Write();
    zjets_h_down->Write();

    zjets_h_old_up->Write();
    zjets_h_old_down->Write();
    zjets_h_up_zeta->Write();
    zjets_h_up_himet->Write();
  }

  wjets_fr_up_h->Write();
  wjets_fr_down_h->Write();
  wjets_mc_up_h->Write();
  wjets_mc_down_h->Write();
  
  outfile->Close();
  delete rbdtg;
}
