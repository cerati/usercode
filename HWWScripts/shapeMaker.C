#include "common.C"

void writeStatUpDown(TH1F* central,bool down, int njets, TString fs) {
  TString proc = TString(central->GetName());
  proc.ReplaceAll("histo_","");
  TString updown = "Up";
  if (down) updown = "Down";
  TH1F* statUpDown = new TH1F(Form("%s_CMS_hww%s_%ij_MVA%sStatBounding_8TeV%s",central->GetName(),TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets,proc.Data(),updown.Data()),
			      Form("%s_CMS_hww%s_%ij_MVA%sStatBounding_8TeV%s",central->GetName(),TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets,proc.Data(),updown.Data()),
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

void shapeMaker(float lumi=4.7, int njets=0, int mass=130, TString fs="sffs", TString plotvar="bdtg") {

  bool inj125 = 0;

  TString dir = main_dir+topww_dir;
  TString dirdy = main_dir+dy_dir;
  TString dirwj = main_dir+wj_dir;

  bool useJson = 0;
  bool applyEff=true;
  bool doFake=false; 
  bool doPUw=true;

  int nbins = 20;
  float minx = -1.;
  float maxx = 1.;
  if (plotvar=="mtmll2D") {
    nbins = 126;
    double  mTBinning[15]    = {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};    
    double  mllBinning[10]   = {12,30,45,60,75,100,125,150,175,200}; 
    mtmll2d_lom = new TH2F("mtmll2d_lom","mtmll2d_lom",14,mTBinning,9,mllBinning);
    mtmll2d_him = new TH2F("mtmll2d_him","mtmll2d_him",10,80.,380.,8,0.,450.);
  }

  TString sigreg = "=dphireg=dphijet=dymvacut=ptll3045=";
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

  int mH = mass;
  //different in case of new signal injection test
  if (inj125) mH=125;

  TFile *weightPDFShapeFILE=0;   
  if(mass < 300.) weightPDFShapeFILE = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/PDFUncertainty_LowMass.root"); 
  else weightPDFShapeFILE = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/PDFUncertainty_HighMass.root"); 
  TH2F *weightPDFShapeUp   = (TH2F*)(weightPDFShapeFILE->Get( Form("qqWW_DF_%ij_alternateUp", njets) ));
  TH2F *weightPDFShapeDown = (TH2F*)(weightPDFShapeFILE->Get( Form("qqWW_DF_%ij_alternateDown", njets) ));

  //Data
  TH1F* data_h = new TH1F("histo_Data","histo_Data",nbins,minx,maxx);
  fillPlot(plotvar,data_h, dir+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, 0, useJson, false, false, false);

  //qqWW
  TH1F* qqww_h = new TH1F("histo_qqWW","histo_qqWW",nbins,minx,maxx);
  fillPlot(plotvar,qqww_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  qqww_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
  //shape variation: 1- mg vs mc@nlo (down mirror);
  TH1F* qqww_mcnlo_h = new TH1F("histo_qqww_mcnlo","histo_qqww_mcnlo",nbins,minx,maxx);
  fillPlot(plotvar,qqww_mcnlo_h, dir+"wwmcnlo"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  avoidNegativeBins(qqww_mcnlo_h);
  scaleIntegral(qqww_h,qqww_mcnlo_h);
  TH1F* qqww_h_up = new TH1F("histo_qqWW_CMS_hww_MVAWWBoundingUp","histo_qqWW_CMS_hww_MVAWWBoundingUp",nbins,minx,maxx);
  qqww_h_up->Add(qqww_mcnlo_h);
  TH1F* qqww_h_down = new TH1F("histo_qqWW_CMS_hww_MVAWWBoundingDown","histo_qqWW_CMS_hww_MVAWWBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(qqww_h,qqww_h_up,qqww_h_down);  
  //shape variation: 2- ratio from mc@nlo w.r.t. QCD up and down
  TH1F* qqww_mcnlo_up_h = new TH1F("histo_qqww_mcnlo_up","histo_qqww_mcnlo_up",nbins,minx,maxx);
  fillPlot(plotvar,qqww_mcnlo_up_h, dir+"wwmcnloup"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  avoidNegativeBins(qqww_mcnlo_up_h);
  scaleIntegral(qqww_h,qqww_mcnlo_up_h);
  divideHistoProtected(qqww_mcnlo_up_h,qqww_mcnlo_h);
  TH1F* qqww_mcnlo_down_h = new TH1F("histo_qqww_mcnlo_down","histo_qqww_mcnlo_down",nbins,minx,maxx);
  fillPlot(plotvar,qqww_mcnlo_down_h, dir+"wwmcnlodown"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  avoidNegativeBins(qqww_mcnlo_down_h);
  scaleIntegral(qqww_h,qqww_mcnlo_down_h);
  divideHistoProtected(qqww_mcnlo_down_h,qqww_mcnlo_h);
  TH1F* qqww_h_nlo_up = new TH1F("histo_qqWW_CMS_hww_MVAWWNLOBoundingUp","histo_qqWW_CMS_hww_MVAWWNLOBoundingUp",nbins,minx,maxx);
  qqww_h_nlo_up->Add(qqww_h);
  multiplyHisto(qqww_h_nlo_up,qqww_mcnlo_up_h);
  TH1F* qqww_h_nlo_down = new TH1F("histo_qqWW_CMS_hww_MVAWWNLOBoundingDown","histo_qqWW_CMS_hww_MVAWWNLOBoundingDown",nbins,minx,maxx);
  qqww_h_nlo_down->Add(qqww_h);
  multiplyHisto(qqww_h_nlo_down,qqww_mcnlo_down_h);
  //shape variation: 3- pdf uncertainty
  TH1F* qqww_pdf_up_h	= unrollHisto2DTo1D(weightPDFShapeUp,"histo_qqWW_CMS_hww_PDFqqWWUp");
  multiplyHisto(qqww_pdf_up_h,qqww_h);
  TH1F* qqww_pdf_down_h	= unrollHisto2DTo1D(weightPDFShapeDown,"histo_qqWW_CMS_hww_PDFqqWWDown");
  multiplyHisto(qqww_pdf_down_h,qqww_h);

  //ggWW
  TH1F* ggww_h = new TH1F("histo_ggWW","histo_ggWW",nbins,minx,maxx);
  fillPlot(plotvar,ggww_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ggww_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
  //shape variation: pdf uncertainty
  TH1F* ggww_pdf_up_h	= unrollHisto2DTo1D(weightPDFShapeUp,"histo_ggWW_CMS_hww_PDFggWWUp");
  multiplyHisto(ggww_pdf_up_h,ggww_h);
  TH1F* ggww_pdf_down_h	= unrollHisto2DTo1D(weightPDFShapeDown,"histo_ggWW_CMS_hww_PDFggWWDown");
  multiplyHisto(ggww_pdf_down_h,ggww_h);

  //VV
  TH1F* wz_h = new TH1F("histo_wz","histo_wz",nbins,minx,maxx);
  fillPlot(plotvar,wz_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* zz_h = new TH1F("histo_zz","histo_zz",nbins,minx,maxx);
  fillPlot(plotvar,zz_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* www_h = new TH1F("histo_www","histo_www",nbins,minx,maxx);
  fillPlot(plotvar,www_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* vv_h = new TH1F("histo_VV","histo_VV",nbins,minx,maxx);
  vv_h->Add(wz_h);
  vv_h->Add(zz_h);
  vv_h->Add(www_h);

  //Top
  TH1F* ttbar_h = new TH1F("histo_ttbar","histo_ttbar",nbins,minx,maxx);
  fillPlot(plotvar,ttbar_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ttbar_h->Scale(TopBkgScaleFactor(njets));
  TH1F* tw_h = new TH1F("histo_tw","histo_tw",nbins,minx,maxx);
  fillPlot(plotvar,tw_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  tw_h->Scale(TopBkgScaleFactor(njets));
  TH1F* top_h = new TH1F("histo_Top","histo_Top",nbins,minx,maxx);
  top_h->Add(ttbar_h);
  top_h->Add(tw_h);
  //shape variation: use madgraph ttbar and ds tw
  TH1F* ttbar_var_h = new TH1F("histo_ttbar_var","histo_ttbar_var",nbins,minx,maxx);
  fillPlot(plotvar,ttbar_var_h, dir+"ttbar"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ttbar_var_h->Scale(1.43);//fixme scale 7 TeV xsec
  //scaleIntegral(ttbar_h,ttbar_var_h);
  TH1F* tw_ds_h = new TH1F("histo_tw_ds","histo_tw_ds",nbins,minx,maxx);
  fillPlot(plotvar,tw_ds_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);//fixme tw_ds
  //scaleIntegral(tw_h,tw_ds_h);
  TH1F* top_h_up = new TH1F("histo_Top_CMS_hww_MVATopBoundingUp","histo_Top_CMS_hww_MVATopBoundingUp",nbins,minx,maxx);
  top_h_up->Add(ttbar_var_h);
  top_h_up->Add(tw_ds_h);
  scaleIntegral(top_h,top_h_up);
  TH1F* top_h_down = new TH1F("histo_Top_CMS_hww_MVATopBoundingDown","histo_Top_CMS_hww_MVATopBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(top_h,top_h_up,top_h_down);  

  //Wgamma
  TH1F* wg_h = new TH1F("histo_wg","histo_wg",nbins,minx,maxx);
  fillPlot(plotvar,wg_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* wgfo_h = new TH1F("histo_wgfo","histo_wgfo",nbins,minx,maxx);
  fillPlot(plotvar,wgfo_h, main_dir+wj_dir+"wgammafo"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw,"wgammafo");
  scaleIntegral(wg_h,wgfo_h);
  TH1F* zg_h = new TH1F("histo_zg","histo_zg",nbins,minx,maxx);
  fillPlot(plotvar,zg_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* wgamma_h = new TH1F("histo_Wgamma","histo_Wgamma",nbins,minx,maxx);
  wgamma_h->Add(wgfo_h);
  wgamma_h->Add(zg_h);

  //Wg3l
  TH1F* wg3l_h = new TH1F("histo_Wg3l","histo_Wg3l",nbins,minx,maxx);
  fillPlot(plotvar,wg3l_h, dir+"wglll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  //Zjets
  float dysf = 1.;
  TH1F* zjets_h = new TH1F("histo_Zjets","histo_Zjets",nbins,minx,maxx);
  TH1F* dyll_lowmet_h = new TH1F("histo_dyll_lowmet","histo_dyll_lowmet",nbins,minx,maxx);
  float dyY = 0;
  if (fs.Contains("sffs")) {
    fillPlot(plotvar,dyll_lowmet_h, dirdy+"dyll"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+fs, lumi, useJson, applyEff, doFake, doPUw);
    dyY = DYBkgScaleFactorBDT(mass,njets);
    dyll_lowmet_h->Scale(dyY/dyll_lowmet_h->Integral());
    zjets_h->Add(dyll_lowmet_h);
  } 
  //new shape variation: loose MET cuts in data
  TH1F* zjets_h_up = new TH1F(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUp",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
			      Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUp",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    float kee = 0.8;//getK(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, 0., useJson, false, doFake, false);
    float lumicorr = (1.+kee*kee)/(2.*kee);
    TH1F* sffs_lowmet_new_h = new TH1F("histo_sffs_lowmet_new","histo_sffs_lowmet_new",nbins,minx,maxx);
    fillPlot(plotvar,sffs_lowmet_new_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", 0, useJson, false, false, false,"");//fixme add zeta method
    TH1F* offs_lowmet_new_h = new TH1F("histo_offs_lowmet_new","histo_offs_lowmet_new",nbins,minx,maxx);
    fillPlot(plotvar,offs_lowmet_new_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=offs=", 0, useJson, false, false, false,"");
    TH1F* wz_lowmet_new_h = new TH1F("histo_wz_lowmet_new","histo_wz_lowmet_new",nbins,minx,maxx);
    fillPlot(plotvar,wz_lowmet_new_h, dirdy+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    TH1F* zz_lowmet_new_h = new TH1F("histo_zz_lowmet_new","histo_zz_lowmet_new",nbins,minx,maxx);
    fillPlot(plotvar,zz_lowmet_new_h, dirdy+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    zjets_h_up->Add(sffs_lowmet_new_h);
    zjets_h_up->Add(offs_lowmet_new_h,-1.*lumicorr);
    zjets_h_up->Add(wz_lowmet_new_h,-1);
    zjets_h_up->Add(zz_lowmet_new_h,-1);
    avoidNegativeBins(zjets_h_up);
    zjets_h_up->Scale(dyY/zjets_h_up->Integral());
  }
  TH1F* zjets_h_down = new TH1F(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingDown",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
				Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingDown",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  fillDownMirrorUp(zjets_h,zjets_h_up,zjets_h_down);

  //old shape variation: use full MET cuts
  TH1F* zjets_h_old_up = new TH1F(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUpOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
				  Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUpOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  TH1F* zjets_h_old_down = new TH1F(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingDownOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
				    Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingDownOld",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    TH1F* dyll_vtx_h = new TH1F("histo_dyll_vtx","histo_dyll_vtx",nbins,minx,maxx);
    fillPlot(plotvar,dyll_vtx_h, dir+"dyll"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
    zjets_h_old_up->Add(dyll_vtx_h);
    scaleIntegral(zjets_h,zjets_h_old_up);
    fillDownMirrorUp(zjets_h,zjets_h_old_up,zjets_h_old_down);
  }

  // OF,VZ subtraction
  TH1F* zjets_h_up_himet = new TH1F(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUpOFHiMet",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
				    Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUpOFHiMet",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    float kee = 0.8;//getK(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, 0., useJson, false, doFake, false);
    float lumicorr = (1.+kee*kee)/(2.*kee);
    TH1F* sffs_himet_h = new TH1F("histo_sffs_himet","histo_sffs_himet",nbins,minx,maxx);
    fillPlot(plotvar,sffs_himet_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=sffs=", 0, useJson, false, false, false,"");
    TH1F* offs_himet_h = new TH1F("histo_offs_himet","histo_offs_himet",nbins,minx,maxx);
    fillPlot(plotvar,offs_himet_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=offs=", 0, useJson, false, false, false,"");
    TH1F* wz_himet_h = new TH1F("histo_wz_himet","histo_wz_himet",nbins,minx,maxx);
    fillPlot(plotvar,wz_himet_h, dirdy+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    TH1F* zz_himet_h = new TH1F("histo_zz_himet","histo_zz_himet",nbins,minx,maxx);
    fillPlot(plotvar,zz_himet_h, dirdy+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_himet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    zjets_h_up_himet->Add(sffs_himet_h);
    zjets_h_up_himet->Add(offs_himet_h,-1.*lumicorr);
    zjets_h_up_himet->Add(wz_himet_h,-1);
    zjets_h_up_himet->Add(zz_himet_h,-1);
    avoidNegativeBins(zjets_h_up_himet);
    zjets_h_up_himet->Scale(dyY/zjets_h_up_himet->Integral());
  }

  // histo_Zjets_CMS_hww_MVAZBounding_hwwsf_0jUpZeta
  // zeta method
  TH1F* zjets_h_up_zeta = new TH1F(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUpZeta",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),
				   Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUpZeta",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets),nbins,minx,maxx);
  if (fs.Contains("sffs")) {
    float kee = 0.8;//getK(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, 0., useJson, false, doFake, false);
    float lumicorr = (1.+kee*kee)/(2.*kee);
    TH1F* sffs_lowmet_zeta_h = new TH1F("histo_sffs_lowmet_zeta","histo_sffs_lowmet_zeta",nbins,minx,maxx);
    fillPlot(plotvar,sffs_lowmet_zeta_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", 0, useJson, false, false, false,"zeta");
    TH1F* offs_lowmet_zeta_h = new TH1F("histo_offs_lowmet_zeta","histo_offs_lowmet_zeta",nbins,minx,maxx);
    fillPlot(plotvar,offs_lowmet_zeta_h, dirdy+"data"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=offs=", 0, useJson, false, false, false,"zeta");
    TH1F* wz_lowmet_zeta_h = new TH1F("histo_wz_lowmet_zeta","histo_wz_lowmet_zeta",nbins,minx,maxx);
    fillPlot(plotvar,wz_lowmet_zeta_h, dirdy+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    TH1F* zz_lowmet_zeta_h = new TH1F("histo_zz_lowmet_zeta","histo_zz_lowmet_zeta",nbins,minx,maxx);
    fillPlot(plotvar,zz_lowmet_zeta_h, dirdy+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg_lowmet+"=sffs=", lumi, useJson, applyEff, doFake, doPUw,"");
    zjets_h_up_zeta->Add(sffs_lowmet_zeta_h);
    zjets_h_up_zeta->Add(offs_lowmet_zeta_h,-1.*lumicorr);
    zjets_h_up_zeta->Add(wz_lowmet_zeta_h,-1);
    zjets_h_up_zeta->Add(zz_lowmet_zeta_h,-1);
    avoidNegativeBins(zjets_h_up_zeta);
    zjets_h_up_zeta->Scale(dyY/zjets_h_up_zeta->Integral());
  }

  //Ztt
  TH1F* dytt_h = new TH1F("histo_dytt","histo_dytt",nbins,minx,maxx);
  fillPlot(plotvar,dytt_h, dir+"data_ztt"+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false);
  TH1F* ztt_h = new TH1F("histo_Ztt","histo_Ztt",nbins,minx,maxx);
  ztt_h->Add(dytt_h);

  //WjetsE
  TH1F* datafake_e_h = new TH1F("datafake_e","datafake_e",nbins,minx,maxx);
  fillPlot(plotvar,datafake_e_h,dirwj+"data"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=", 0, useJson, false, true, false);
  TH1F* qqwwfake_e_h = new TH1F("qqwwfake_e","qqwwfake_e",nbins,minx,maxx);
  fillPlot(plotvar,qqwwfake_e_h,dirwj+"qqww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_e_h = new TH1F("ggwwfake_e","ggwwfake_e",nbins,minx,maxx);
  fillPlot(plotvar,ggwwfake_e_h,dirwj+"ggww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_e_h = new TH1F("ttbarfake_e","ttbarfake_e",nbins,minx,maxx);
  fillPlot(plotvar,ttbarfake_e_h,dirwj+"ttbar_powheg"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_e_h = new TH1F("twfake_e","twfake_e",nbins,minx,maxx);
  fillPlot(plotvar,twfake_e_h,dirwj+"tw"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_e_h = new TH1F("wzfake_e","wzfake_e",nbins,minx,maxx);
  fillPlot(plotvar,wzfake_e_h,dirwj+"wz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_e_h = new TH1F("zzfake_e","zzfake_e",nbins,minx,maxx);
  fillPlot(plotvar,zzfake_e_h,dirwj+"zz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wgfake_e_h = new TH1F("wgfake_e","wgfake_e",nbins,minx,maxx);
  fillPlot(plotvar,wgfake_e_h,dirwj+"wgamma"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wg3lfake_e_h = new TH1F("wg3lfake_e","wg3lfake_e",nbins,minx,maxx);
  fillPlot(plotvar,wg3lfake_e_h,dirwj+"wglll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* dyllfake_e_h = new TH1F("dyllfake_e","dyllfake_e",nbins,minx,maxx);
  fillPlot(plotvar,dyllfake_e_h,dirwj+"dyll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wjetsE_h = new TH1F("histo_WjetsE","histo_WjetsE",nbins,minx,maxx);
  wjetsE_h->Add(datafake_e_h);
  wjetsE_h->Add(qqwwfake_e_h,-1.);
  wjetsE_h->Add(ggwwfake_e_h,-1.);
  wjetsE_h->Add(ttbarfake_e_h,-1.);
  wjetsE_h->Add(twfake_e_h,-1.);
  wjetsE_h->Add(wzfake_e_h,-1.);
  wjetsE_h->Add(zzfake_e_h,-1.);
  wjetsE_h->Add(wgfake_e_h,-1.);
  wjetsE_h->Add(wg3lfake_e_h,-1.);
  wjetsE_h->Add(dyllfake_e_h,-1.);
  float intgr_wj_e = wjetsE_h->Integral();
  avoidNegativeBins(wjetsE_h);
  wjetsE_h->Scale(intgr_wj_e/wjetsE_h->Integral());
  //syst 1: MC closure test
  /*
  TH1F* wjetsE_mc_up_h = new TH1F("histo_WjetsE_CMS_hww_MVAWMCBoundingUp","histo_WjetsE_CMS_hww_MVAWMCBoundingUp",nbins,minx,maxx);
  fillPlot(plotvar,wjetsE_mc_up_h,dirwj+"wjetsE"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);//fixme is this correct??? shouldn't be doFake=0?
  scaleIntegral(wjetsE_h,wjetsE_mc_up_h);
  TH1F* wjetsE_mc_down_h = new TH1F("histo_WjetsE_CMS_hww_MVAWMCBoundingDown","histo_WjetsE_CMS_hww_MVAWMCBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(wjetsE_h,wjetsE_mc_up_h,wjetsE_mc_down_h);
  */
  //syst 2: alternative fakebale object definition
  TH1F* datafake_e_fr_up_h = new TH1F("datafake_e_fr_up","datafake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,datafake_e_fr_up_h,dirwj+"data"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", 0, useJson, false, true, false);
  TH1F* qqwwfake_e_fr_up_h = new TH1F("qqwwfake_e_fr_up","qqwwfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,qqwwfake_e_fr_up_h,dirwj+"qqww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_e_fr_up_h = new TH1F("ggwwfake_e_fr_up","ggwwfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,ggwwfake_e_fr_up_h,dirwj+"ggww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_e_fr_up_h = new TH1F("ttbarfake_e_fr_up","ttbarfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,ttbarfake_e_fr_up_h,dirwj+"ttbar_powheg"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_e_fr_up_h = new TH1F("twfake_e_fr_up","twfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,twfake_e_fr_up_h,dirwj+"tw"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_e_fr_up_h = new TH1F("wzfake_e_fr_up","wzfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,wzfake_e_fr_up_h,dirwj+"wz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_e_fr_up_h = new TH1F("zzfake_e_fr_up","zzfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,zzfake_e_fr_up_h,dirwj+"zz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wgfake_e_fr_up_h = new TH1F("wgfake_e_fr_up","wgfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,wgfake_e_fr_up_h,dirwj+"wgamma"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wg3lfake_e_fr_up_h = new TH1F("wg3lfake_e_fr_up","wg3lfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,wg3lfake_e_fr_up_h,dirwj+"wglll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* dyllfake_e_fr_up_h = new TH1F("dyllfake_e_fr_up","dyllfake_e_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,dyllfake_e_fr_up_h,dirwj+"dyll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=elfake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wjetsE_fr_up_h = new TH1F("histo_WjetsE_CMS_hww_MVAWEBoundingUp","histo_WjetsE_CMS_hww_MVAWEBoundingUp",nbins,minx,maxx);
  wjetsE_fr_up_h->Add(datafake_e_fr_up_h);
  wjetsE_fr_up_h->Add(qqwwfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(ggwwfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(ttbarfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(twfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(wzfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(zzfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(wgfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(wg3lfake_e_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(dyllfake_e_fr_up_h,-1.);
  avoidNegativeBins(wjetsE_fr_up_h);
  scaleIntegral(wjetsE_h,wjetsE_fr_up_h);
  TH1F* wjetsE_fr_down_h = new TH1F("histo_WjetsE_CMS_hww_MVAWEBoundingDown","histo_WjetsE_CMS_hww_MVAWEBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(wjetsE_h,wjetsE_fr_up_h,wjetsE_fr_down_h);

  //WjetsM
  TH1F* datafake_m_h = new TH1F("datafake_m","datafake_m",nbins,minx,maxx);
  fillPlot(plotvar,datafake_m_h,dirwj+"data"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=", 0, useJson, false, true, false);
  TH1F* qqwwfake_m_h = new TH1F("qqwwfake_m","qqwwfake_m",nbins,minx,maxx);
  fillPlot(plotvar,qqwwfake_m_h,dirwj+"qqww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_m_h = new TH1F("ggwwfake_m","ggwwfake_m",nbins,minx,maxx);
  fillPlot(plotvar,ggwwfake_m_h,dirwj+"ggww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_m_h = new TH1F("ttbarfake_m","ttbarfake_m",nbins,minx,maxx);
  fillPlot(plotvar,ttbarfake_m_h,dirwj+"ttbar_powheg"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_m_h = new TH1F("twfake_m","twfake_m",nbins,minx,maxx);
  fillPlot(plotvar,twfake_m_h,dirwj+"tw"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_m_h = new TH1F("wzfake_m","wzfake_m",nbins,minx,maxx);
  fillPlot(plotvar,wzfake_m_h,dirwj+"wz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_m_h = new TH1F("zzfake_m","zzfake_m",nbins,minx,maxx);
  fillPlot(plotvar,zzfake_m_h,dirwj+"zz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wgfake_m_h = new TH1F("wgfake_m","wgfake_m",nbins,minx,maxx);
  fillPlot(plotvar,wgfake_m_h,dirwj+"wgamma"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wg3lfake_m_h = new TH1F("wg3lfake_m","wg3lfake_m",nbins,minx,maxx);
  fillPlot(plotvar,wg3lfake_m_h,dirwj+"wglll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* dyllfake_m_h = new TH1F("dyllfake_m","dyllfake_m",nbins,minx,maxx);
  fillPlot(plotvar,dyllfake_m_h,dirwj+"dyll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=spill=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wjetsM_h = new TH1F("histo_WjetsM","histo_WjetsM",nbins,minx,maxx);
  wjetsM_h->Add(datafake_m_h);
  wjetsM_h->Add(qqwwfake_m_h,-1.);
  wjetsM_h->Add(ggwwfake_m_h,-1.);
  wjetsM_h->Add(ttbarfake_m_h,-1.);
  wjetsM_h->Add(twfake_m_h,-1.);
  wjetsM_h->Add(wzfake_m_h,-1.);
  wjetsM_h->Add(zzfake_m_h,-1.);
  wjetsM_h->Add(wgfake_m_h,-1.);
  wjetsM_h->Add(wg3lfake_m_h,-1.);
  wjetsM_h->Add(dyllfake_m_h,-1.);
  float intgr_wj_m = wjetsM_h->Integral();
  avoidNegativeBins(wjetsM_h);
  wjetsM_h->Scale(intgr_wj_m/wjetsM_h->Integral());
  //syst 1: MC closure test
  /*
  TH1F* wjetsM_mc_up_h = new TH1F("histo_WjetsM_CMS_hww_MVAWMCBoundingUp","histo_WjetsM_CMS_hww_MVAWMCBoundingUp",nbins,minx,maxx);
  fillPlot(plotvar,wjetsM_mc_up_h,dirwj+"wjetsM"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);//fixme is this correct??? shouldn't be doFake=0?
  scaleIntegral(wjetsM_h,wjetsM_mc_up_h);
  TH1F* wjetsM_mc_down_h = new TH1F("histo_WjetsM_CMS_hww_MVAWMCBoundingDown","histo_WjetsM_CMS_hww_MVAWMCBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(wjetsM_h,wjetsM_mc_up_h,wjetsM_mc_down_h);
  */
  //syst 2: alternative fakebale object definition
  TH1F* datafake_m_fr_up_h = new TH1F("datafake_m_fr_up","datafake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,datafake_m_fr_up_h,dirwj+"data"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", 0, useJson, false, true, false);
  TH1F* qqwwfake_m_fr_up_h = new TH1F("qqwwfake_m_fr_up","qqwwfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,qqwwfake_m_fr_up_h,dirwj+"qqww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_m_fr_up_h = new TH1F("ggwwfake_m_fr_up","ggwwfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,ggwwfake_m_fr_up_h,dirwj+"ggww"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_m_fr_up_h = new TH1F("ttbarfake_m_fr_up","ttbarfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,ttbarfake_m_fr_up_h,dirwj+"ttbar_powheg"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_m_fr_up_h = new TH1F("twfake_m_fr_up","twfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,twfake_m_fr_up_h,dirwj+"tw"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_m_fr_up_h = new TH1F("wzfake_m_fr_up","wzfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,wzfake_m_fr_up_h,dirwj+"wz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_m_fr_up_h = new TH1F("zzfake_m_fr_up","zzfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,zzfake_m_fr_up_h,dirwj+"zz"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wgfake_m_fr_up_h = new TH1F("wgfake_m_fr_up","wgfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,wgfake_m_fr_up_h,dirwj+"wgamma"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wg3lfake_m_fr_up_h = new TH1F("wg3lfake_m_fr_up","wg3lfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,wg3lfake_m_fr_up_h,dirwj+"wglll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* dyllfake_m_fr_up_h = new TH1F("dyllfake_m_fr_up","dyllfake_m_fr_up",nbins,minx,maxx);
  fillPlot(plotvar,dyllfake_m_fr_up_h,dirwj+"dyll"+suffix, wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=alternativeFR=mufake=", lumi, useJson, applyEff, true, doPUw);
  TH1F* wjetsM_fr_up_h = new TH1F("histo_WjetsM_CMS_hww_MVAWMBoundingUp","histo_WjetsM_CMS_hww_MVAWMBoundingUp",nbins,minx,maxx);
  wjetsM_fr_up_h->Add(datafake_m_fr_up_h);
  wjetsM_fr_up_h->Add(qqwwfake_m_fr_up_h,-1.);
  wjetsM_fr_up_h->Add(ggwwfake_m_fr_up_h,-1.);
  wjetsM_fr_up_h->Add(ttbarfake_m_fr_up_h,-1.);
  wjetsM_fr_up_h->Add(twfake_m_fr_up_h,-1.);
  wjetsM_fr_up_h->Add(wzfake_m_fr_up_h,-1.);
  wjetsM_fr_up_h->Add(zzfake_m_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(wgfake_m_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(wg3lfake_m_fr_up_h,-1.);
  wjetsE_fr_up_h->Add(dyllfake_m_fr_up_h,-1.);
  avoidNegativeBins(wjetsM_fr_up_h);
  scaleIntegral(wjetsM_h,wjetsM_fr_up_h);
  TH1F* wjetsM_fr_down_h = new TH1F("histo_WjetsM_CMS_hww_MVAWMBoundingDown","histo_WjetsM_CMS_hww_MVAWMBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(wjetsM_h,wjetsM_fr_up_h,wjetsM_fr_down_h);

  //Higgs
  TH1F* ggH_h = new TH1F("histo_ggH","histo_ggH",nbins,minx,maxx);
  fillPlot(plotvar,ggH_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw);
  //ggH_h->Scale(InterfgHHSystematics(mass));
  TH1F* qqH_h = new TH1F("histo_qqH","histo_qqH",nbins,minx,maxx);
  fillPlot(plotvar,qqH_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* WH_h = new TH1F("histo_WH","histo_WH",nbins,minx,maxx);
  fillPlot(plotvar,WH_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* ZH_h = new TH1F("histo_ZH","histo_ZH",nbins,minx,maxx);
  fillPlot(plotvar,ZH_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw);

  /*
  //ggH k-factor syst  
  TH1F* ggH_up_h = new TH1F("histo_ggH_CMS_hww_MVAggHBoundingUp","histo_ggH_CMS_hww_MVAggHBoundingUp",nbins,minx,maxx);
  fillPlot(plotvar,ggH_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "ggH_k_syst_up");
  ggH_up_h->Scale(InterfgHHSystematics(mass));
  TH1F* ggH_down_h = new TH1F("histo_ggH_CMS_hww_MVAggHBoundingDown","histo_ggH_CMS_hww_MVAggHBoundingDown",nbins,minx,maxx);
  fillPlot(plotvar,ggH_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "ggH_k_syst_down");
  ggH_down_h->Scale(InterfgHHSystematics(mass));
  */

  //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
  TH1F *ggH_metres_up_h=0, *ggH_metres_down_h=0, *qqH_metres_up_h=0, *qqH_metres_down_h=0, *WH_metres_up_h=0, *WH_metres_down_h=0, *ZH_metres_up_h=0, *ZH_metres_down_h=0, *qqww_metres_up_h=0, *qqww_metres_down_h=0, 
    *ggww_metres_up_h=0, *ggww_metres_down_h=0, *vv_metres_up_h=0, *vv_metres_down_h=0, *top_metres_up_h=0, *top_metres_down_h=0, *wgamma_metres_up_h=0, *wgamma_metres_down_h=0;
  //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
  TH1F *ggH_lepres_up_h=0, *ggH_lepres_down_h=0, *qqH_lepres_up_h=0, *qqH_lepres_down_h=0, *WH_lepres_up_h=0, *WH_lepres_down_h=0, *ZH_lepres_up_h=0, *ZH_lepres_down_h=0, *qqww_lepres_up_h=0, *qqww_lepres_down_h=0, 
    *ggww_lepres_up_h=0, *ggww_lepres_down_h=0, *vv_lepres_up_h=0, *vv_lepres_down_h=0, *top_lepres_up_h=0, *top_lepres_down_h=0, *wgamma_lepres_up_h=0, *wgamma_lepres_down_h=0;
  //JES systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
  TH1F *ggH_jes_up_h=0, *ggH_jes_down_h=0, *qqH_jes_up_h=0, *qqH_jes_down_h=0, *WH_jes_up_h=0, *WH_jes_down_h=0, *ZH_jes_up_h=0, *ZH_jes_down_h=0, *qqww_jes_up_h=0, *qqww_jes_down_h=0, *ggww_jes_up_h=0, 
    *ggww_jes_down_h=0, *vv_jes_up_h=0, *vv_jes_down_h=0, *top_jes_up_h=0, *top_jes_down_h=0, *wgamma_jes_up_h=0, *wgamma_jes_down_h=0;
  //Lepton efficiency systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Wgamma (why no Top???)
  TH1F *ggH_lepeff_up_h=0, *ggH_lepeff_down_h=0, *qqH_lepeff_up_h=0, *qqH_lepeff_down_h=0, *WH_lepeff_up_h=0, *WH_lepeff_down_h=0, *ZH_lepeff_up_h=0, *ZH_lepeff_down_h=0, *qqww_lepeff_up_h=0, *qqww_lepeff_down_h=0, *ggww_lepeff_up_h=0, *ggww_lepeff_down_h=0, *vv_lepeff_up_h=0, *vv_lepeff_down_h=0, *wgamma_lepeff_up_h=0, *wgamma_lepeff_down_h=0;

  if (doResEffSyst) {
    //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
    //Higgs
    ggH_metres_up_h = new TH1F("histo_ggH_CMS_hww_MVAMETResBoundingUp","histo_ggH_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggH_metres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    scaleIntegral(ggH_h,ggH_metres_up_h);//fixme
    ggH_metres_down_h = new TH1F("histo_ggH_CMS_hww_MVAMETResBoundingDown","histo_ggH_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ggH_h,ggH_metres_up_h,ggH_metres_down_h);
    qqH_metres_up_h = new TH1F("histo_qqH_CMS_hww_MVAMETResBoundingUp","histo_qqH_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqH_metres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    qqH_metres_down_h = new TH1F("histo_qqH_CMS_hww_MVAMETResBoundingDown","histo_qqH_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(qqH_h,qqH_metres_up_h,qqH_metres_down_h);
    WH_metres_up_h = new TH1F("histo_WH_CMS_hww_MVAMETResBoundingUp","histo_WH_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,WH_metres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    WH_metres_down_h = new TH1F("histo_WH_CMS_hww_MVAMETResBoundingDown","histo_WH_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(WH_h,WH_metres_up_h,WH_metres_down_h);
    ZH_metres_up_h = new TH1F("histo_ZH_CMS_hww_MVAMETResBoundingUp","histo_ZH_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ZH_metres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    ZH_metres_down_h = new TH1F("histo_ZH_CMS_hww_MVAMETResBoundingDown","histo_ZH_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ZH_h,ZH_metres_up_h,ZH_metres_down_h);
    //WW
    qqww_metres_up_h = new TH1F("histo_qqWW_CMS_hww_MVAMETResBoundingUp","histo_qqWW_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqww_metres_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    qqww_metres_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    qqww_metres_down_h = new TH1F("histo_qqWW_CMS_hww_MVAMETResBoundingDown","histo_qqWW_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(qqww_h,qqww_metres_up_h,qqww_metres_down_h);
    ggww_metres_up_h = new TH1F("histo_ggWW_CMS_hww_MVAMETResBoundingUp","histo_ggWW_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggww_metres_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    ggww_metres_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_metres_down_h = new TH1F("histo_ggWW_CMS_hww_MVAMETResBoundingDown","histo_ggWW_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(ggww_h,ggww_metres_up_h,ggww_metres_down_h);
    //VV
    TH1F* wz_metres_up_h = new TH1F("histo_wz_metres_up","histo_wz_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,wz_metres_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    TH1F* zz_metres_up_h = new TH1F("histo_zz_metres_up","histo_zz_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,zz_metres_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    TH1F* www_metres_up_h = new TH1F("histo_www_metres_up","histo_www_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,www_metres_up_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    vv_metres_up_h = new TH1F("histo_VV_CMS_hww_MVAMETResBoundingUp","histo_VV_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    vv_metres_up_h->Add(wz_metres_up_h);
    vv_metres_up_h->Add(zz_metres_up_h);
    vv_metres_up_h->Add(www_metres_up_h);
    vv_metres_down_h = new TH1F("histo_VV_CMS_hww_MVAMETResBoundingDown","histo_VV_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(vv_h,vv_metres_up_h,vv_metres_down_h);
    //Top
    TH1F* ttbar_metres_up_h = new TH1F("histo_ttbar_metres_up","histo_ttbar_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,ttbar_metres_up_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    ttbar_metres_up_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_metres_up_h = new TH1F("histo_tw_metres_up","histo_tw_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,tw_metres_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    tw_metres_up_h->Scale(TopBkgScaleFactor(njets));
    top_metres_up_h = new TH1F("histo_Top_CMS_hww_MVAMETResBoundingUp","histo_Top_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    top_metres_up_h->Add(ttbar_metres_up_h);
    top_metres_up_h->Add(tw_metres_up_h);
    top_metres_down_h = new TH1F("histo_Top_CMS_hww_MVAMETResBoundingDown","histo_Top_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(top_h,top_metres_up_h,top_metres_down_h);
    //Wgamma
    /*
    TH1F* wg_metres_up_h = new TH1F("histo_wg_metres_up","histo_wg_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,wg_metres_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    TH1F* zg_metres_up_h = new TH1F("histo_zg_metres_up","histo_zg_metres_up",nbins,minx,maxx);
    fillPlot(plotvar,zg_metres_up_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
    wgamma_metres_up_h = new TH1F("histo_Wgamma_CMS_hww_MVAMETResBoundingUp","histo_Wgamma_CMS_hww_MVAMETResBoundingUp",nbins,minx,maxx);
    wgamma_metres_up_h->Add(wg_metres_up_h);
    wgamma_metres_up_h->Add(zg_metres_up_h);
    wgamma_metres_down_h = new TH1F("histo_Wgamma_CMS_hww_MVAMETResBoundingDown","histo_Wgamma_CMS_hww_MVAMETResBoundingDown",nbins,minx,maxx);
    fillDownMirrorUp(wgamma_h,wgamma_metres_up_h,wgamma_metres_down_h);
    */
    //when WW and Top are from data need to normalize histo!
    if (mass<=200){
      scaleIntegral(ggww_h,ggww_metres_down_h);
      scaleIntegral(ggww_h,ggww_metres_up_h);
      scaleIntegral(qqww_h,qqww_metres_down_h);
      scaleIntegral(qqww_h,qqww_metres_up_h);
    }
    scaleIntegral(top_h,top_metres_down_h);
    scaleIntegral(top_h,top_metres_up_h);

    //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
    ggH_lepres_up_h = new TH1F("histo_ggH_CMS_hww_MVALepResBoundingUp","histo_ggH_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggH_lepres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    //ggH_lepres_up_h->Scale(InterfgHHSystematics(mass));
    ggH_lepres_down_h = new TH1F("histo_ggH_CMS_hww_MVALepResBoundingDown","histo_ggH_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ggH_lepres_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    //ggH_lepres_down_h->Scale(InterfgHHSystematics(mass));
    qqH_lepres_up_h = new TH1F("histo_qqH_CMS_hww_MVALepResBoundingUp","histo_qqH_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqH_lepres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    qqH_lepres_down_h = new TH1F("histo_qqH_CMS_hww_MVALepResBoundingDown","histo_qqH_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,qqH_lepres_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    WH_lepres_up_h = new TH1F("histo_WH_CMS_hww_MVALepResBoundingUp","histo_WH_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,WH_lepres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    WH_lepres_down_h = new TH1F("histo_WH_CMS_hww_MVALepResBoundingDown","histo_WH_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,WH_lepres_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    ZH_lepres_up_h = new TH1F("histo_ZH_CMS_hww_MVALepResBoundingUp","histo_ZH_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ZH_lepres_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ZH_lepres_down_h = new TH1F("histo_ZH_CMS_hww_MVALepResBoundingDown","histo_ZH_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ZH_lepres_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    //WW
    qqww_lepres_up_h = new TH1F("histo_qqWW_CMS_hww_MVALepResBoundingUp","histo_qqWW_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqww_lepres_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    qqww_lepres_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_lepres_up_h = new TH1F("histo_ggWW_CMS_hww_MVALepResBoundingUp","histo_ggWW_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggww_lepres_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ggww_lepres_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    qqww_lepres_down_h = new TH1F("histo_qqWW_CMS_hww_MVALepResBoundingDown","histo_qqWW_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,qqww_lepres_down_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    qqww_lepres_down_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_lepres_down_h = new TH1F("histo_ggWW_CMS_hww_MVALepResBoundingDown","histo_ggWW_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ggww_lepres_down_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    ggww_lepres_down_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    //VV
    TH1F* wz_lepres_up_h = new TH1F("histo_wz_lepres_up","histo_wz_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,wz_lepres_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    TH1F* zz_lepres_up_h = new TH1F("histo_zz_lepres_up","histo_zz_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,zz_lepres_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    TH1F* www_lepres_up_h = new TH1F("histo_www_lepres_up","histo_www_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,www_lepres_up_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    vv_lepres_up_h = new TH1F("histo_VV_CMS_hww_MVALepResBoundingUp","histo_VV_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    vv_lepres_up_h->Add(wz_lepres_up_h);
    vv_lepres_up_h->Add(zz_lepres_up_h);
    vv_lepres_up_h->Add(www_lepres_up_h);
    TH1F* wz_lepres_down_h = new TH1F("histo_wz_lepres_down","histo_wz_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,wz_lepres_down_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    TH1F* zz_lepres_down_h = new TH1F("histo_zz_lepres_down","histo_zz_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,zz_lepres_down_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    TH1F* www_lepres_down_h = new TH1F("histo_www_lepres_down","histo_www_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,www_lepres_down_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    vv_lepres_down_h = new TH1F("histo_VV_CMS_hww_MVALepResBoundingDown","histo_VV_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    vv_lepres_down_h->Add(wz_lepres_down_h);
    vv_lepres_down_h->Add(zz_lepres_down_h);
    vv_lepres_down_h->Add(www_lepres_down_h);
    //Top
    TH1F* ttbar_lepres_up_h = new TH1F("histo_ttbar_lepres_up","histo_ttbar_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,ttbar_lepres_up_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    ttbar_lepres_up_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_lepres_up_h = new TH1F("histo_tw_lepres_up","histo_tw_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,tw_lepres_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    tw_lepres_up_h->Scale(TopBkgScaleFactor(njets));
    top_lepres_up_h = new TH1F("histo_Top_CMS_hww_MVALepResBoundingUp","histo_Top_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    top_lepres_up_h->Add(ttbar_lepres_up_h);
    top_lepres_up_h->Add(tw_lepres_up_h);
    TH1F* ttbar_lepres_down_h = new TH1F("histo_ttbar_lepres_down","histo_ttbar_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,ttbar_lepres_down_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    ttbar_lepres_down_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_lepres_down_h = new TH1F("histo_tw_lepres_down","histo_tw_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,tw_lepres_down_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    tw_lepres_down_h->Scale(TopBkgScaleFactor(njets));
    top_lepres_down_h = new TH1F("histo_Top_CMS_hww_MVALepResBoundingDown","histo_Top_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    top_lepres_down_h->Add(ttbar_lepres_down_h);
    top_lepres_down_h->Add(tw_lepres_down_h);
    //Wgamma
    /*
    TH1F* wg_lepres_up_h = new TH1F("histo_wg_lepres_up","histo_wg_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,wg_lepres_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    TH1F* zg_lepres_up_h = new TH1F("histo_zg_lepres_up","histo_zg_lepres_up",nbins,minx,maxx);
    fillPlot(plotvar,zg_lepres_up_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
    wgamma_lepres_up_h = new TH1F("histo_Wgamma_CMS_hww_MVALepResBoundingUp","histo_Wgamma_CMS_hww_MVALepResBoundingUp",nbins,minx,maxx);
    wgamma_lepres_up_h->Add(wg_lepres_up_h);
    wgamma_lepres_up_h->Add(zg_lepres_up_h);
    TH1F* wg_lepres_down_h = new TH1F("histo_wg_lepres_down","histo_wg_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,wg_lepres_down_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    TH1F* zg_lepres_down_h = new TH1F("histo_zg_lepres_down","histo_zg_lepres_down",nbins,minx,maxx);
    fillPlot(plotvar,zg_lepres_down_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
    wgamma_lepres_down_h = new TH1F("histo_Wgamma_CMS_hww_MVALepResBoundingDown","histo_Wgamma_CMS_hww_MVALepResBoundingDown",nbins,minx,maxx);
    wgamma_lepres_down_h->Add(wg_lepres_down_h);
    wgamma_lepres_down_h->Add(zg_lepres_down_h);
    */
    //when WW and Top are from data need to normalize histo!
    if (mass<=200){
      scaleIntegral(ggww_h,ggww_lepres_down_h);
      scaleIntegral(ggww_h,ggww_lepres_up_h);
      scaleIntegral(qqww_h,qqww_lepres_down_h);
      scaleIntegral(qqww_h,qqww_lepres_up_h);
    }
    scaleIntegral(top_h,top_lepres_down_h);
    scaleIntegral(top_h,top_lepres_up_h);

    //JES systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
    ggH_jes_up_h = new TH1F("histo_ggH_CMS_hww_MVAJESBoundingUp","histo_ggH_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggH_jes_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    //ggH_jes_up_h->Scale(InterfgHHSystematics(mass));
    ggH_jes_down_h = new TH1F("histo_ggH_CMS_hww_MVAJESBoundingDown","histo_ggH_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ggH_jes_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    //ggH_jes_down_h->Scale(InterfgHHSystematics(mass));
    qqH_jes_up_h = new TH1F("histo_qqH_CMS_hww_MVAJESBoundingUp","histo_qqH_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqH_jes_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    qqH_jes_down_h = new TH1F("histo_qqH_CMS_hww_MVAJESBoundingDown","histo_qqH_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,qqH_jes_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    WH_jes_up_h = new TH1F("histo_WH_CMS_hww_MVAJESBoundingUp","histo_WH_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,WH_jes_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    WH_jes_down_h = new TH1F("histo_WH_CMS_hww_MVAJESBoundingDown","histo_WH_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,WH_jes_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    ZH_jes_up_h = new TH1F("histo_ZH_CMS_hww_MVAJESBoundingUp","histo_ZH_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ZH_jes_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ZH_jes_down_h = new TH1F("histo_ZH_CMS_hww_MVAJESBoundingDown","histo_ZH_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ZH_jes_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    //WW
    qqww_jes_up_h = new TH1F("histo_qqWW_CMS_hww_MVAJESBoundingUp","histo_qqWW_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqww_jes_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    qqww_jes_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_jes_up_h = new TH1F("histo_ggWW_CMS_hww_MVAJESBoundingUp","histo_ggWW_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggww_jes_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ggww_jes_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    qqww_jes_down_h = new TH1F("histo_qqWW_CMS_hww_MVAJESBoundingDown","histo_qqWW_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,qqww_jes_down_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    qqww_jes_down_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_jes_down_h = new TH1F("histo_ggWW_CMS_hww_MVAJESBoundingDown","histo_ggWW_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ggww_jes_down_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    ggww_jes_down_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    //VV
    TH1F* wz_jes_up_h = new TH1F("histo_wz_jes_up","histo_wz_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,wz_jes_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    TH1F* zz_jes_up_h = new TH1F("histo_zz_jes_up","histo_zz_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,zz_jes_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    TH1F* www_jes_up_h = new TH1F("histo_www_jes_up","histo_www_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,www_jes_up_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    vv_jes_up_h = new TH1F("histo_VV_CMS_hww_MVAJESBoundingUp","histo_VV_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    vv_jes_up_h->Add(wz_jes_up_h);
    vv_jes_up_h->Add(zz_jes_up_h);
    vv_jes_up_h->Add(www_jes_up_h);
    TH1F* wz_jes_down_h = new TH1F("histo_wz_jes_down","histo_wz_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,wz_jes_down_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    TH1F* zz_jes_down_h = new TH1F("histo_zz_jes_down","histo_zz_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,zz_jes_down_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    TH1F* www_jes_down_h = new TH1F("histo_www_jes_down","histo_www_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,www_jes_down_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    vv_jes_down_h = new TH1F("histo_VV_CMS_hww_MVAJESBoundingDown","histo_VV_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    vv_jes_down_h->Add(wz_jes_down_h);
    vv_jes_down_h->Add(zz_jes_down_h);
    vv_jes_down_h->Add(www_jes_down_h);
    //Top
    TH1F* ttbar_jes_up_h = new TH1F("histo_ttbar_jes_up","histo_ttbar_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,ttbar_jes_up_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    ttbar_jes_up_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_jes_up_h = new TH1F("histo_tw_jes_up","histo_tw_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,tw_jes_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    tw_jes_up_h->Scale(TopBkgScaleFactor(njets));
    top_jes_up_h = new TH1F("histo_Top_CMS_hww_MVAJESBoundingUp","histo_Top_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    top_jes_up_h->Add(ttbar_jes_up_h);
    top_jes_up_h->Add(tw_jes_up_h);
    TH1F* ttbar_jes_down_h = new TH1F("histo_ttbar_jes_down","histo_ttbar_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,ttbar_jes_down_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    ttbar_jes_down_h->Scale(TopBkgScaleFactor(njets));
    TH1F* tw_jes_down_h = new TH1F("histo_tw_jes_down","histo_tw_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,tw_jes_down_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    tw_jes_down_h->Scale(TopBkgScaleFactor(njets));
    top_jes_down_h = new TH1F("histo_Top_CMS_hww_MVAJESBoundingDown","histo_Top_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    top_jes_down_h->Add(ttbar_jes_down_h);
    top_jes_down_h->Add(tw_jes_down_h);
    //Wgamma
    /*
    TH1F* wg_jes_up_h = new TH1F("histo_wg_jes_up","histo_wg_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,wg_jes_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    TH1F* zg_jes_up_h = new TH1F("histo_zg_jes_up","histo_zg_jes_up",nbins,minx,maxx);
    fillPlot(plotvar,zg_jes_up_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesUp");
    wgamma_jes_up_h = new TH1F("histo_Wgamma_CMS_hww_MVAJESBoundingUp","histo_Wgamma_CMS_hww_MVAJESBoundingUp",nbins,minx,maxx);
    wgamma_jes_up_h->Add(wg_jes_up_h);
    wgamma_jes_up_h->Add(zg_jes_up_h);
    TH1F* wg_jes_down_h = new TH1F("histo_wg_jes_down","histo_wg_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,wg_jes_down_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    TH1F* zg_jes_down_h = new TH1F("histo_zg_jes_down","histo_zg_jes_down",nbins,minx,maxx);
    fillPlot(plotvar,zg_jes_down_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "jesDown");
    wgamma_jes_down_h = new TH1F("histo_Wgamma_CMS_hww_MVAJESBoundingDown","histo_Wgamma_CMS_hww_MVAJESBoundingDown",nbins,minx,maxx);
    wgamma_jes_down_h->Add(wg_jes_down_h);
    wgamma_jes_down_h->Add(zg_jes_down_h);
    */
    //when WW and Top are from data need to normalize histo!
    if (mass<=200){
      scaleIntegral(ggww_h,ggww_jes_down_h);
      scaleIntegral(ggww_h,ggww_jes_up_h);
      scaleIntegral(qqww_h,qqww_jes_down_h);
      scaleIntegral(qqww_h,qqww_jes_up_h);
    }
    scaleIntegral(top_h,top_jes_down_h);
    scaleIntegral(top_h,top_jes_up_h);
    //Lepton efficiency systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Wgamma (no Top???)
    ggH_lepeff_up_h = new TH1F("histo_ggH_CMS_hww_MVALepEffBoundingUp","histo_ggH_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggH_lepeff_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    //ggH_lepeff_up_h->Scale(InterfgHHSystematics(mass));
    ggH_lepeff_down_h = new TH1F("histo_ggH_CMS_hww_MVALepEffBoundingDown","histo_ggH_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ggH_lepeff_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //ggH_lepeff_down_h->Scale(InterfgHHSystematics(mass));
    qqH_lepeff_up_h = new TH1F("histo_qqH_CMS_hww_MVALepEffBoundingUp","histo_qqH_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqH_lepeff_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    qqH_lepeff_down_h = new TH1F("histo_qqH_CMS_hww_MVALepEffBoundingDown","histo_qqH_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,qqH_lepeff_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    WH_lepeff_up_h = new TH1F("histo_WH_CMS_hww_MVALepEffBoundingUp","histo_WH_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,WH_lepeff_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    WH_lepeff_down_h = new TH1F("histo_WH_CMS_hww_MVALepEffBoundingDown","histo_WH_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,WH_lepeff_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    ZH_lepeff_up_h = new TH1F("histo_ZH_CMS_hww_MVALepEffBoundingUp","histo_ZH_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ZH_lepeff_up_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    ZH_lepeff_down_h = new TH1F("histo_ZH_CMS_hww_MVALepEffBoundingDown","histo_ZH_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ZH_lepeff_down_h, dir+Form("hww%i",mH)+suffix, wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //WW
    qqww_lepeff_up_h = new TH1F("histo_qqWW_CMS_hww_MVALepEffBoundingUp","histo_qqWW_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,qqww_lepeff_up_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    qqww_lepeff_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_lepeff_up_h = new TH1F("histo_ggWW_CMS_hww_MVALepEffBoundingUp","histo_ggWW_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    fillPlot(plotvar,ggww_lepeff_up_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    ggww_lepeff_up_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    qqww_lepeff_down_h = new TH1F("histo_qqWW_CMS_hww_MVALepEffBoundingDown","histo_qqWW_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,qqww_lepeff_down_h, dir+"qqww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    qqww_lepeff_down_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    ggww_lepeff_down_h = new TH1F("histo_ggWW_CMS_hww_MVALepEffBoundingDown","histo_ggWW_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    fillPlot(plotvar,ggww_lepeff_down_h, dir+"ggww"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    ggww_lepeff_down_h->Scale(WWBkgScaleFactorMVA(TMath::Min(TMath::Max((int)mass,115),200),njets));
    //VV
    TH1F* wz_lepeff_up_h = new TH1F("histo_wz_lepeff_up","histo_wz_lepeff_up",nbins,minx,maxx);
    fillPlot(plotvar,wz_lepeff_up_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    TH1F* zz_lepeff_up_h = new TH1F("histo_zz_lepeff_up","histo_zz_lepeff_up",nbins,minx,maxx);
    fillPlot(plotvar,zz_lepeff_up_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    TH1F* www_lepeff_up_h = new TH1F("histo_www_lepeff_up","histo_www_lepeff_up",nbins,minx,maxx);
    fillPlot(plotvar,www_lepeff_up_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    vv_lepeff_up_h = new TH1F("histo_VV_CMS_hww_MVALepEffBoundingUp","histo_VV_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    vv_lepeff_up_h->Add(wz_lepeff_up_h);
    vv_lepeff_up_h->Add(zz_lepeff_up_h);
    vv_lepeff_up_h->Add(www_lepeff_up_h);
    TH1F* wz_lepeff_down_h = new TH1F("histo_wz_lepeff_down","histo_wz_lepeff_down",nbins,minx,maxx);
    fillPlot(plotvar,wz_lepeff_down_h, dir+"wz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    TH1F* zz_lepeff_down_h = new TH1F("histo_zz_lepeff_down","histo_zz_lepeff_down",nbins,minx,maxx);
    fillPlot(plotvar,zz_lepeff_down_h, dir+"zz"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    TH1F* www_lepeff_down_h = new TH1F("histo_www_lepeff_down","histo_www_lepeff_down",nbins,minx,maxx);
    fillPlot(plotvar,www_lepeff_down_h, dir+"www"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    vv_lepeff_down_h = new TH1F("histo_VV_CMS_hww_MVALepEffBoundingDown","histo_VV_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    vv_lepeff_down_h->Add(wz_lepeff_down_h);
    vv_lepeff_down_h->Add(zz_lepeff_down_h);
    vv_lepeff_down_h->Add(www_lepeff_down_h);
    //Top
    //TH1F* ttbar_lepeff_up_h = new TH1F("histo_ttbar_lepeff_up","histo_ttbar_lepeff_up",nbins,minx,maxx);
    //fillPlot(plotvar,ttbar_lepeff_up_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    //ttbar_lepeff_up_h->Scale(TopBkgScaleFactor(njets));
    //TH1F* tw_lepeff_up_h = new TH1F("histo_tw_lepeff_up","histo_tw_lepeff_up",nbins,minx,maxx);
    //fillPlot(plotvar,tw_lepeff_up_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    //tw_lepeff_up_h->Scale(TopBkgScaleFactor(njets));
    //top_lepeff_up_h = new TH1F("histo_Top_CMS_hww_MVALepEffBoundingUp","histo_Top_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    //top_lepeff_up_h->Add(ttbar_lepeff_up_h);
    //top_lepeff_up_h->Add(tw_lepeff_up_h);
    //TH1F* ttbar_lepeff_down_h = new TH1F("histo_ttbar_lepeff_down","histo_ttbar_lepeff_down",nbins,minx,maxx);
    //fillPlot(plotvar,ttbar_lepeff_down_h, dir+"ttbar_powheg"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //ttbar_lepeff_down_h->Scale(TopBkgScaleFactor(njets));
    //TH1F* tw_lepeff_down_h = new TH1F("histo_tw_lepeff_down","histo_tw_lepeff_down",nbins,minx,maxx);
    //fillPlot(plotvar,tw_lepeff_down_h, dir+"tw"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    //tw_lepeff_down_h->Scale(TopBkgScaleFactor(njets));
    //top_lepeff_down_h = new TH1F("histo_Top_CMS_hww_MVALepEffBoundingDown","histo_Top_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    //top_lepeff_down_h->Add(ttbar_lepeff_down_h);
    //top_lepeff_down_h->Add(tw_lepeff_down_h);
    //Wgamma
    /*
    TH1F* wg_lepeff_up_h = new TH1F("histo_wg_lepeff_up","histo_wg_lepeff_up",nbins,minx,maxx);
    fillPlot(plotvar,wg_lepeff_up_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    TH1F* zg_lepeff_up_h = new TH1F("histo_zg_lepeff_up","histo_zg_lepeff_up",nbins,minx,maxx);
    fillPlot(plotvar,zg_lepeff_up_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffUp");
    wgamma_lepeff_up_h = new TH1F("histo_Wgamma_CMS_hww_MVALepEffBoundingUp","histo_Wgamma_CMS_hww_MVALepEffBoundingUp",nbins,minx,maxx);
    wgamma_lepeff_up_h->Add(wg_lepeff_up_h);
    wgamma_lepeff_up_h->Add(zg_lepeff_up_h);
    TH1F* wg_lepeff_down_h = new TH1F("histo_wg_lepeff_down","histo_wg_lepeff_down",nbins,minx,maxx);
    fillPlot(plotvar,wg_lepeff_down_h, dir+"wgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    TH1F* zg_lepeff_down_h = new TH1F("histo_zg_lepeff_down","histo_zg_lepeff_down",nbins,minx,maxx);
    fillPlot(plotvar,zg_lepeff_down_h, dir+"zgamma"+suffix, wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "lepeffDown");
    wgamma_lepeff_down_h = new TH1F("histo_Wgamma_CMS_hww_MVALepEffBoundingDown","histo_Wgamma_CMS_hww_MVALepEffBoundingDown",nbins,minx,maxx);
    wgamma_lepeff_down_h->Add(wg_lepeff_down_h);
    wgamma_lepeff_down_h->Add(zg_lepeff_down_h);
    */
    //when WW and Top are from data need to normalize histo!
    if (mass<=200){
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

  if (1) {
    //reduce the number of channels by merging WH and ZH
    WH_h->Add(ZH_h);
    ZH_h->Reset();
    WH_metres_up_h->Add(ZH_metres_up_h);
    ZH_metres_up_h->Reset();
    WH_metres_down_h->Add(ZH_metres_down_h);
    ZH_metres_down_h->Reset();
    WH_lepres_up_h->Add(ZH_lepres_up_h);
    ZH_lepres_up_h->Reset();
    WH_lepres_down_h->Add(ZH_lepres_down_h);
    ZH_lepres_down_h->Reset();
    WH_jes_up_h->Add(ZH_jes_up_h);
    ZH_jes_up_h->Reset();
    WH_jes_down_h->Add(ZH_jes_down_h);
    ZH_jes_down_h->Reset();
    WH_lepeff_up_h->Add(ZH_lepeff_up_h);
    ZH_lepeff_up_h->Reset();
    WH_lepeff_down_h->Add(ZH_lepeff_down_h);
    ZH_lepeff_down_h->Reset();
  }


  //Write histos to file
  TString outfname = Form("hww%s_%ij.input_8TeV.root",TString(fs).ReplaceAll("fs","").ReplaceAll("=","").Data(),njets);
  TFile* outfile = TFile::Open(outfname,"RECREATE");

  //Nominal Histos
  data_h->Write();
  ggww_h->Write();
  qqww_h->Write();
  top_h->Write();
  zjets_h->Write();
  wjetsE_h->Write();
  wjetsM_h->Write();
  vv_h->Write();
  wgamma_h->Write();
  wg3l_h->Write();
  ztt_h->Write();
  qqH_h->Write();
  ggH_h->Write();
  WH_h->Write();
  ZH_h->Write();

  //ggH k-factor syst  
  /*
  ggH_up_h->Write();
  ggH_down_h->Write();
  */

  //Stat uncertainty
  writeStatUpDown(ggww_h,njets,fs);
  writeStatUpDown(qqww_h,njets,fs);
  writeStatUpDown(top_h,njets,fs);
  writeStatUpDown(zjets_h,njets,fs);
  writeStatUpDown(wjetsE_h,njets,fs);
  writeStatUpDown(wjetsM_h,njets,fs);
  writeStatUpDown(vv_h,njets,fs);
  writeStatUpDown(wgamma_h,njets,fs);
  writeStatUpDown(wg3l_h,njets,fs);
  writeStatUpDown(ztt_h,njets,fs);
  writeStatUpDown(qqH_h,njets,fs);
  writeStatUpDown(ggH_h,njets,fs);
  writeStatUpDown(WH_h,njets,fs);
  writeStatUpDown(ZH_h,njets,fs);

  if (doResEffSyst) {
    //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
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
    //wgamma_metres_up_h->Write();
    //wgamma_metres_down_h->Write();

    //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
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
    //wgamma_lepres_up_h->Write();
    //wgamma_lepres_down_h->Write();

    //JES systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma
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
    //wgamma_jes_up_h->Write();
    //wgamma_jes_down_h->Write();

    //Lepton efficiency systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Wgamma (why no Top???)
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
    //wgamma_lepeff_up_h->Write();
    //wgamma_lepeff_down_h->Write();
  }

  //Other systematics
  qqww_h_up->Write();
  qqww_h_down->Write();
  qqww_h_nlo_up->Write();
  qqww_h_nlo_down->Write();

  qqww_pdf_up_h->Write();
  qqww_pdf_down_h->Write();
  ggww_pdf_up_h->Write();
  ggww_pdf_down_h->Write();

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

  wjetsE_fr_up_h->Write();
  wjetsE_fr_down_h->Write();
  wjetsM_fr_up_h->Write();
  wjetsM_fr_down_h->Write();
  /*
  wjets_mc_up_h->Write();
  wjets_mc_down_h->Write();
  */
  
  outfile->Close();
  delete rbdtg;
  if (plotvar=="mtmll2D") {
    delete mtmll2d_lom;
    delete mtmll2d_him;
  }
  weightPDFShapeFILE->Close();   
}
