#include "common.C"

void fillDownMirrorUp(TH1F* central,TH1F* up,TH1F* down) {
  down->Add(up);
  down->Scale(-1);
  down->Add(central);
  down->Add(central);
  //need to avoid negative values...
  for (int bin=1;bin<=down->GetNbinsX();++bin) {
    if (down->GetBinContent(bin)<0) down->SetBinContent(bin,0);
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
  TH1F* statUpDown = new TH1F(TString(central->GetName())+"_CMS_MVA"+proc+Form("StatBounding_hww%s_%ij",TString(fs).ReplaceAll("fs","").Data(),njets)+updown,
			  TString(central->GetTitle())+"_CMS_MVA"+proc+Form("StatBounding_hww%s_%ij",TString(fs).ReplaceAll("fs","").Data(),njets)+updown,
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

  TString sigreg = "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20,";
  TString sigreg_met2040 = sigreg;
  sigreg_met2040.ReplaceAll("minmetvtx","met2040");

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

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  //Data
  TH1F* data_h = new TH1F("histo_Data","histo_Data",nbins,minx,maxx);
  fillPlot("bdtg",data_h, dir+"data"+suffix, wwSelection, veto, mass, njets, sigreg+fs, 0, useJson, false, false, false);

  //qqWW
  TH1F* qqww_h = new TH1F("histo_qqWW","histo_qqWW",nbins,minx,maxx);
  fillPlot("bdtg",qqww_h, dir+"qqww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  qqww_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  //shape variation: 1- mg vs mc@nlo (down mirror);
  TH1F* qqww_mcnlo_h = new TH1F("histo_qqww_mcnlo","histo_qqww_mcnlo",nbins,minx,maxx);
  fillPlot("bdtg",qqww_mcnlo_h, dir+"ww_mcnlo"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(qqww_h,qqww_mcnlo_h);
  TH1F* qqww_h_up = new TH1F("histo_qqWW_CMS_MVAWWBounding_hwwUp","histo_qqWW_CMS_MVAWWBounding_hwwUp",nbins,minx,maxx);
  qqww_h_up->Add(qqww_mcnlo_h);
  TH1F* qqww_h_down = new TH1F("histo_qqWW_CMS_MVAWWBounding_hwwDown","histo_qqWW_CMS_MVAWWBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(qqww_h,qqww_h_up,qqww_h_down);  
  //shape variation: 2- ratio from mc@nlo w.r.t. QCD up and down
  TH1F* qqww_mcnlo_up_h = new TH1F("histo_qqww_mcnlo_up","histo_qqww_mcnlo_up",nbins,minx,maxx);
  fillPlot("bdtg",qqww_mcnlo_up_h, dir+"ww_mcnlo_up"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(qqww_h,qqww_mcnlo_up_h);
  qqww_mcnlo_up_h->Divide(qqww_mcnlo_h);
  TH1F* qqww_mcnlo_down_h = new TH1F("histo_qqww_mcnlo_down","histo_qqww_mcnlo_down",nbins,minx,maxx);
  fillPlot("bdtg",qqww_mcnlo_down_h, dir+"ww_mcnlo_down"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(qqww_h,qqww_mcnlo_down_h);
  qqww_mcnlo_down_h->Divide(qqww_mcnlo_h);
  TH1F* qqww_h_nlo_up = new TH1F("histo_qqWW_CMS_MVAWWNLOBounding_hwwUp","histo_qqWW_CMS_MVAWWNLOBounding_hwwUp",nbins,minx,maxx);
  qqww_h_nlo_up->Add(qqww_h);
  qqww_h_nlo_up->Multiply(qqww_mcnlo_up_h);
  TH1F* qqww_h_nlo_down = new TH1F("histo_qqWW_CMS_MVAWWNLOBounding_hwwDown","histo_qqWW_CMS_MVAWWNLOBounding_hwwDown",nbins,minx,maxx);
  qqww_h_nlo_down->Add(qqww_h);
  qqww_h_nlo_down->Multiply(qqww_mcnlo_down_h);

  //ggWW
  TH1F* ggww_h = new TH1F("histo_ggWW","histo_ggWW",nbins,minx,maxx);
  fillPlot("bdtg",ggww_h, dir+"ggww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ggww_h->Scale(WWBkgScaleFactorMVA(mass,njets));

  //VV
  TH1F* wz_h = new TH1F("histo_wz","histo_wz",nbins,minx,maxx);
  fillPlot("bdtg",wz_h, dir+"wz"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* zz_h = new TH1F("histo_zz","histo_zz",nbins,minx,maxx);
  fillPlot("bdtg",zz_h, dir+"zz_py"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* vv_h = new TH1F("histo_VV","histo_VV",nbins,minx,maxx);
  vv_h->Add(wz_h);
  vv_h->Add(zz_h);

  //Top
  TH1F* ttbar_h = new TH1F("histo_ttbar","histo_ttbar",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_h, dir+"ttbar"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  ttbar_h->Scale(TopBkgScaleFactor(njets));
  TH1F* tw_h = new TH1F("histo_tw","histo_tw",nbins,minx,maxx);
  fillPlot("bdtg",tw_h, dir+"tw"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  tw_h->Scale(TopBkgScaleFactor(njets));
  TH1F* top_h = new TH1F("histo_Top","histo_Top",nbins,minx,maxx);
  top_h->Add(ttbar_h);
  top_h->Add(tw_h);
  //shape variation: use madgraph ttbar and ds tw
  TH1F* ttbar_mg_h = new TH1F("histo_ttbar_mg","histo_ttbar_mg",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_mg_h, dir+"ttbar_mg"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(ttbar_h,ttbar_mg_h);
  TH1F* tw_ds_h = new TH1F("histo_tw_ds","histo_tw_ds",nbins,minx,maxx);
  fillPlot("bdtg",tw_ds_h, dir+"tw_ds"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(tw_h,tw_ds_h);
  TH1F* top_h_up = new TH1F("histo_Top_CMS_MVATopBounding_hwwUp","histo_Top_CMS_MVATopBounding_hwwUp",nbins,minx,maxx);
  top_h_up->Add(ttbar_mg_h);
  top_h_up->Add(tw_ds_h);
  TH1F* top_h_down = new TH1F("histo_Top_CMS_MVATopBounding_hwwDown","histo_Top_CMS_MVATopBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(top_h,top_h_up,top_h_down);  

  //Wgamma
  TH1F* wg_h = new TH1F("histo_wg","histo_wg",nbins,minx,maxx);
  fillPlot("bdtg",wg_h, dir+"wgamma"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* wg3l_h = new TH1F("histo_wg3l","histo_wg3l",nbins,minx,maxx);
  fillPlot("bdtg",wg3l_h, dir+"wg3l"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  wg3l_h->Scale(WGstarScaleFactor());
  TH1F* wgamma_h = new TH1F("histo_Wgamma","histo_Wgamma",nbins,minx,maxx);
  wgamma_h->Add(wg_h);
  wgamma_h->Add(wg3l_h);

  //Zjets
  float dysf = 1.;
  if (fs=="sffs") dysf = DYBkgScaleFactor(0,njets);
  TH1F* dyee_vtx_h = new TH1F("histo_dyee_vtx","histo_dyee_vtx",nbins,minx,maxx);
  fillPlot("bdtg",dyee_vtx_h, dir+"dyee"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  dyee_vtx_h->Scale(dysf);
  TH1F* dymm_vtx_h = new TH1F("histo_dymm_vtx","histo_dymm_vtx",nbins,minx,maxx);
  fillPlot("bdtg",dymm_vtx_h, dir+"dymm"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  dymm_vtx_h->Scale(dysf);
  TH1F* dyee_2040_h = new TH1F("histo_dyee_2040","histo_dyee_2040",nbins,minx,maxx);
  fillPlot("bdtg",dyee_2040_h, dirdy+"dyee"+suffix, wwSelNoMet, veto, mass, njets, sigreg_met2040+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(dyee_vtx_h,dyee_2040_h);
  TH1F* dymm_2040_h = new TH1F("histo_dymm_2040","histo_dymm_2040",nbins,minx,maxx);
  fillPlot("bdtg",dymm_2040_h, dirdy+"dymm"+suffix, wwSelNoMet, veto, mass, njets, sigreg_met2040+fs, lumi, useJson, applyEff, doFake, doPUw);
  scaleIntegral(dymm_vtx_h,dymm_2040_h);
  TH1F* pwz_h = new TH1F("histo_pwz","histo_pwz",nbins,minx,maxx);
  fillPlot("bdtg",pwz_h, dir+"wz"+suffix, wwSelection, veto, mass, njets, sigreg+"fromZ,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* pzz_h = new TH1F("histo_pzz","histo_pzz",nbins,minx,maxx);
  fillPlot("bdtg",pzz_h, dir+"zz_py"+suffix, wwSelection, veto, mass, njets, sigreg+"fromZ,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* zjets_h = new TH1F("histo_Zjets","histo_Zjets",nbins,minx,maxx);
  zjets_h->Add(dyee_2040_h);
  zjets_h->Add(dymm_2040_h);
  zjets_h->Add(pwz_h);
  zjets_h->Add(pzz_h);
  //shape variation: use full MET cuts
  TH1F* zjets_h_up = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUp",TString(fs).ReplaceAll("fs","").Data(),njets),Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijUp",TString(fs).ReplaceAll("fs","").Data(),njets),nbins,minx,maxx);
  TH1F* zjets_h_down = new TH1F(Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijDown",TString(fs).ReplaceAll("fs","").Data(),njets),Form("histo_Zjets_CMS_MVAZBounding_hww%s_%ijDown",TString(fs).ReplaceAll("fs","").Data(),njets),nbins,minx,maxx);
  if (fs=="sffs") {
    zjets_h_up->Add(dyee_vtx_h);
    zjets_h_up->Add(dymm_vtx_h);
    zjets_h_up->Add(pwz_h);
    zjets_h_up->Add(pzz_h);
    fillDownMirrorUp(zjets_h,zjets_h_up,zjets_h_down);
  }

  //Ztt
  TH1F* dytt_1_h = new TH1F("histo_dytt_1","histo_dytt_1",nbins,minx,maxx);
  fillPlot("bdtg",dytt_1_h, dir+"data-emb-tau121"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  TH1F* dytt_2_h = new TH1F("histo_dytt_2","histo_dytt_2",nbins,minx,maxx);
  fillPlot("bdtg",dytt_2_h, dir+"data-emb-tau122"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  TH1F* dytt_3_h = new TH1F("histo_dytt_3","histo_dytt_3",nbins,minx,maxx);
  fillPlot("bdtg",dytt_3_h, dir+"data-emb-tau123"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  TH1F* ztt_h = new TH1F("histo_Ztt","histo_Ztt",nbins,minx,maxx);
  ztt_h->Add(dytt_1_h);
  ztt_h->Add(dytt_2_h);
  ztt_h->Add(dytt_3_h);

  //Wjets
  TH1F* datafake_h = new TH1F("datafake","datafake",nbins,minx,maxx);
  fillPlot("bdtg",datafake_h,dir+"data"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, 0, useJson, false, true, false);
  TH1F* qqwwfake_h = new TH1F("qqwwfake","qqwwfake",nbins,minx,maxx);
  fillPlot("bdtg",qqwwfake_h,dir+"qqww"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_h = new TH1F("ggwwfake","ggwwfake",nbins,minx,maxx);
  fillPlot("bdtg",ggwwfake_h,dir+"ggww"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_h = new TH1F("ttbarfake","ttbarfake",nbins,minx,maxx);
  fillPlot("bdtg",ttbarfake_h,dir+"ttbar"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_h = new TH1F("twfake","twfake",nbins,minx,maxx);
  fillPlot("bdtg",twfake_h,dir+"tw"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_h = new TH1F("wzfake","wzfake",nbins,minx,maxx);
  fillPlot("bdtg",wzfake_h,dir+"wz"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_h = new TH1F("zzfake","zzfake",nbins,minx,maxx);
  fillPlot("bdtg",zzfake_h,dir+"zz_py"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wjets_h = new TH1F("histo_Wjets","histo_Wjets",nbins,minx,maxx);
  wjets_h->Add(datafake_h);
  wjets_h->Add(qqwwfake_h,-1.);
  wjets_h->Add(ggwwfake_h,-1.);
  wjets_h->Add(ttbarfake_h,-1.);
  wjets_h->Add(twfake_h,-1.);
  wjets_h->Add(wzfake_h,-1.);
  wjets_h->Add(zzfake_h,-1.);
  //syst 1: MC closure test
  TH1F* wjets_mc_up_h = new TH1F("histo_Wjets_CMS_MVAWMCBounding_hwwUp","histo_Wjets_CMS_MVAWMCBounding_hwwUp",nbins,minx,maxx);
  fillPlot("bdtg",wjets_mc_up_h,dir+"wjets"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, true, doPUw);
  scaleIntegral(wjets_h,wjets_mc_up_h);
  TH1F* wjets_mc_down_h = new TH1F("histo_Wjets_CMS_MVAWMCBounding_hwwDown","histo_Wjets_CMS_MVAWMCBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(wjets_h,wjets_mc_up_h,wjets_mc_down_h);
  //syst 2: alternative fakebale object definition
  TH1F* datafake_fr_up_h = new TH1F("datafake_fr_up","datafake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",datafake_fr_up_h,dir+"data"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, 0, useJson, false, true, false);
  TH1F* qqwwfake_fr_up_h = new TH1F("qqwwfake_fr_up","qqwwfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",qqwwfake_fr_up_h,dir+"qqww"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ggwwfake_fr_up_h = new TH1F("ggwwfake_fr_up","ggwwfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",ggwwfake_fr_up_h,dir+"ggww"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* ttbarfake_fr_up_h = new TH1F("ttbarfake_fr_up","ttbarfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",ttbarfake_fr_up_h,dir+"ttbar"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* twfake_fr_up_h = new TH1F("twfake_fr_up","twfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",twfake_fr_up_h,dir+"tw"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wzfake_fr_up_h = new TH1F("wzfake_fr_up","wzfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",wzfake_fr_up_h,dir+"wz"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* zzfake_fr_up_h = new TH1F("zzfake_fr_up","zzfake_fr_up",nbins,minx,maxx);
  fillPlot("bdtg",zzfake_fr_up_h,dir+"zz_py"+suffix, wwSelectionNoLep, veto, mass, njets, sigreg+"alternativeFR,"+fs, lumi, useJson, applyEff, true, doPUw);
  TH1F* wjets_fr_up_h = new TH1F("histo_Wjets_CMS_MVAWBounding_hwwUp","histo_Wjets_CMS_MVAWBounding_hwwUp",nbins,minx,maxx);
  wjets_fr_up_h->Add(datafake_fr_up_h);
  wjets_fr_up_h->Add(qqwwfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(ggwwfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(ttbarfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(twfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(wzfake_fr_up_h,-1.);
  wjets_fr_up_h->Add(zzfake_fr_up_h,-1.);
  TH1F* wjets_fr_down_h = new TH1F("histo_Wjets_CMS_MVAWBounding_hwwDown","histo_Wjets_CMS_MVAWBounding_hwwDown",nbins,minx,maxx);
  fillDownMirrorUp(wjets_h,wjets_fr_up_h,wjets_fr_down_h);

  //Higgs
  TH1F* ggH_h = new TH1F("histo_ggH","histo_ggH",nbins,minx,maxx);
  fillPlot("bdtg",ggH_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ggH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* qqH_h = new TH1F("histo_qqH","histo_qqH",nbins,minx,maxx);
  fillPlot("bdtg",qqH_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"qqH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* WH_h = new TH1F("histo_WH","histo_WH",nbins,minx,maxx);
  fillPlot("bdtg",WH_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"WH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  TH1F* ZH_h = new TH1F("histo_ZH","histo_ZH",nbins,minx,maxx);
  fillPlot("bdtg",ZH_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ZH,"+fs, lumi, useJson, applyEff, doFake, doPUw);

  //MET RESOLUTION SYSTEMATICS: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
  //Higgs
  TH1F* ggH_metres_up_h = new TH1F("histo_ggH_CMS_MVAMETResBoundingUp","histo_ggH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ggH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ggH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  scaleIntegral(ggH_h,ggH_metres_up_h);//fixme
  TH1F* ggH_metres_down_h = new TH1F("histo_ggH_CMS_MVAMETResBoundingDown","histo_ggH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(ggH_h,ggH_metres_up_h,ggH_metres_down_h);
  TH1F* qqH_metres_up_h = new TH1F("histo_qqH_CMS_MVAMETResBoundingUp","histo_qqH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",qqH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"qqH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  TH1F* qqH_metres_down_h = new TH1F("histo_qqH_CMS_MVAMETResBoundingDown","histo_qqH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(qqH_h,qqH_metres_up_h,qqH_metres_down_h);
  TH1F* WH_metres_up_h = new TH1F("histo_WH_CMS_MVAMETResBoundingUp","histo_WH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",WH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"WH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  TH1F* WH_metres_down_h = new TH1F("histo_WH_CMS_MVAMETResBoundingDown","histo_WH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(WH_h,WH_metres_up_h,WH_metres_down_h);
  TH1F* ZH_metres_up_h = new TH1F("histo_ZH_CMS_MVAMETResBoundingUp","histo_ZH_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ZH_metres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ZH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  TH1F* ZH_metres_down_h = new TH1F("histo_ZH_CMS_MVAMETResBoundingDown","histo_ZH_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(ZH_h,ZH_metres_up_h,ZH_metres_down_h);
  //WW
  TH1F* qqww_metres_up_h = new TH1F("histo_qqWW_CMS_MVAMETResBoundingUp","histo_qqWW_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",qqww_metres_up_h, dir+"qqww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  qqww_metres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  TH1F* qqww_metres_down_h = new TH1F("histo_qqWW_CMS_MVAMETResBoundingDown","histo_qqWW_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(qqww_h,qqww_metres_up_h,qqww_metres_down_h);
  TH1F* ggww_metres_up_h = new TH1F("histo_ggWW_CMS_MVAMETResBoundingUp","histo_ggWW_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ggww_metres_up_h, dir+"ggww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  ggww_metres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  TH1F* ggww_metres_down_h = new TH1F("histo_ggWW_CMS_MVAMETResBoundingDown","histo_ggWW_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(ggww_h,ggww_metres_up_h,ggww_metres_down_h);
  //VV
  TH1F* wz_metres_up_h = new TH1F("histo_wz_metres_up","histo_wz_metres_up",nbins,minx,maxx);
  fillPlot("bdtg",wz_metres_up_h, dir+"wz"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  TH1F* zz_metres_up_h = new TH1F("histo_zz_metres_up","histo_zz_metres_up",nbins,minx,maxx);
  fillPlot("bdtg",zz_metres_up_h, dir+"zz_py"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  TH1F* vv_metres_up_h = new TH1F("histo_VV_CMS_MVAMETResBoundingUp","histo_VV_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  vv_metres_up_h->Add(wz_metres_up_h);
  vv_metres_up_h->Add(zz_metres_up_h);
  TH1F* vv_metres_down_h = new TH1F("histo_VV_CMS_MVAMETResBoundingDown","histo_VV_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(vv_h,vv_metres_up_h,vv_metres_down_h);
  //Top
  TH1F* ttbar_metres_up_h = new TH1F("histo_ttbar_metres_up","histo_ttbar_metres_up",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_metres_up_h, dir+"ttbar"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  ttbar_metres_up_h->Scale(TopBkgScaleFactor(njets));
  TH1F* tw_metres_up_h = new TH1F("histo_tw_metres_up","histo_tw_metres_up",nbins,minx,maxx);
  fillPlot("bdtg",tw_metres_up_h, dir+"tw"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  tw_metres_up_h->Scale(TopBkgScaleFactor(njets));
  TH1F* top_metres_up_h = new TH1F("histo_Top_CMS_MVAMETResBoundingUp","histo_Top_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  top_metres_up_h->Add(ttbar_metres_up_h);
  top_metres_up_h->Add(tw_metres_up_h);
  TH1F* top_metres_down_h = new TH1F("histo_Top_CMS_MVAMETResBoundingDown","histo_Top_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(top_h,top_metres_up_h,top_metres_down_h);
  //Wgamma
  TH1F* wg_metres_up_h = new TH1F("histo_wg_metres_up","histo_wg_metres_up",nbins,minx,maxx);
  fillPlot("bdtg",wg_metres_up_h, dir+"wgamma"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  TH1F* wg3l_metres_up_h = new TH1F("histo_wg3l_metres_up","histo_wg3l_metres_up",nbins,minx,maxx);
  fillPlot("bdtg",wg3l_metres_up_h, dir+"wg3l"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "metSmear");
  wg3l_metres_up_h->Scale(WGstarScaleFactor());
  TH1F* wgamma_metres_up_h = new TH1F("histo_Wgamma_CMS_MVAMETResBoundingUp","histo_Wgamma_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  wgamma_metres_up_h->Add(wg_metres_up_h);
  wgamma_metres_up_h->Add(wg3l_metres_up_h);
  TH1F* wgamma_metres_down_h = new TH1F("histo_Wgamma_CMS_MVAMETResBoundingDown","histo_Wgamma_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(wgamma_h,wgamma_metres_up_h,wgamma_metres_down_h);
  //Ztt
  TH1F* dytt_1_metres_up_h = new TH1F("histo_dytt_metres_up_1","histo_dytt_metres_up_1",nbins,minx,maxx);
  fillPlot("bdtg",dytt_1_metres_up_h, dir+"data-emb-tau121"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "metSmear");
  TH1F* dytt_2_metres_up_h = new TH1F("histo_dytt_metres_up_2","histo_dytt_metres_up_2",nbins,minx,maxx);
  fillPlot("bdtg",dytt_2_metres_up_h, dir+"data-emb-tau122"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "metSmear");
  TH1F* dytt_3_metres_up_h = new TH1F("histo_dytt_metres_up_3","histo_dytt_metres_up_3",nbins,minx,maxx);
  fillPlot("bdtg",dytt_3_metres_up_h, dir+"data-emb-tau123"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "metSmear");
  TH1F* ztt_metres_up_h = new TH1F("histo_Ztt_CMS_MVAMETResBoundingUp","histo_Ztt_CMS_MVAMETResBoundingUp",nbins,minx,maxx);
  ztt_metres_up_h->Add(dytt_1_metres_up_h);
  ztt_metres_up_h->Add(dytt_2_metres_up_h);
  ztt_metres_up_h->Add(dytt_3_metres_up_h);
  TH1F* ztt_metres_down_h = new TH1F("histo_Ztt_CMS_MVAMETResBoundingDown","histo_Ztt_CMS_MVAMETResBoundingDown",nbins,minx,maxx);
  fillDownMirrorUp(ztt_h,ztt_metres_up_h,ztt_metres_down_h);

  //Lepton energy resolution and scale systematics: Affecting ZH, WH, qqH, ggH, qqWW, ggWW, VV, Top, Wgamma and Ztt components
  TH1F* ggH_lepres_up_h = new TH1F("histo_ggH_CMS_MVALepResBoundingUp","histo_ggH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ggH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ggH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* ggH_lepres_down_h = new TH1F("histo_ggH_CMS_MVALepResBoundingDown","histo_ggH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",ggH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ggH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  TH1F* qqH_lepres_up_h = new TH1F("histo_qqH_CMS_MVALepResBoundingUp","histo_qqH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",qqH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"qqH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* qqH_lepres_down_h = new TH1F("histo_qqH_CMS_MVALepResBoundingDown","histo_qqH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",qqH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"qqH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  TH1F* WH_lepres_up_h = new TH1F("histo_WH_CMS_MVALepResBoundingUp","histo_WH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",WH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"WH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* WH_lepres_down_h = new TH1F("histo_WH_CMS_MVALepResBoundingDown","histo_WH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",WH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"WH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  TH1F* ZH_lepres_up_h = new TH1F("histo_ZH_CMS_MVALepResBoundingUp","histo_ZH_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ZH_lepres_up_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ZH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* ZH_lepres_down_h = new TH1F("histo_ZH_CMS_MVALepResBoundingDown","histo_ZH_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",ZH_lepres_down_h, dir+Form("hww%i",mass)+suffix, wwSelection, veto, mass, njets, sigreg+"ZH,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  //WW
  TH1F* qqww_lepres_up_h = new TH1F("histo_qqWW_CMS_MVALepResBoundingUp","histo_qqWW_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",qqww_lepres_up_h, dir+"qqww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  qqww_lepres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  TH1F* ggww_lepres_up_h = new TH1F("histo_ggWW_CMS_MVALepResBoundingUp","histo_ggWW_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  fillPlot("bdtg",ggww_lepres_up_h, dir+"ggww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  ggww_lepres_up_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  TH1F* qqww_lepres_down_h = new TH1F("histo_qqWW_CMS_MVALepResBoundingDown","histo_qqWW_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",qqww_lepres_down_h, dir+"qqww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  qqww_lepres_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  TH1F* ggww_lepres_down_h = new TH1F("histo_ggWW_CMS_MVALepResBoundingDown","histo_ggWW_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  fillPlot("bdtg",ggww_lepres_down_h, dir+"ggww"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  ggww_lepres_down_h->Scale(WWBkgScaleFactorMVA(mass,njets));
  //VV
  TH1F* wz_lepres_up_h = new TH1F("histo_wz_lepres_up","histo_wz_lepres_up",nbins,minx,maxx);
  fillPlot("bdtg",wz_lepres_up_h, dir+"wz"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* zz_lepres_up_h = new TH1F("histo_zz_lepres_up","histo_zz_lepres_up",nbins,minx,maxx);
  fillPlot("bdtg",zz_lepres_up_h, dir+"zz_py"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* vv_lepres_up_h = new TH1F("histo_VV_CMS_MVALepResBoundingUp","histo_VV_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  vv_lepres_up_h->Add(wz_lepres_up_h);
  vv_lepres_up_h->Add(zz_lepres_up_h);
  TH1F* wz_lepres_down_h = new TH1F("histo_wz_lepres_down","histo_wz_lepres_down",nbins,minx,maxx);
  fillPlot("bdtg",wz_lepres_down_h, dir+"wz"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  TH1F* zz_lepres_down_h = new TH1F("histo_zz_lepres_down","histo_zz_lepres_down",nbins,minx,maxx);
  fillPlot("bdtg",zz_lepres_down_h, dir+"zz_py"+suffix, wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  TH1F* vv_lepres_down_h = new TH1F("histo_VV_CMS_MVALepResBoundingDown","histo_VV_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  vv_lepres_down_h->Add(wz_lepres_down_h);
  vv_lepres_down_h->Add(zz_lepres_down_h);
  //Top
  TH1F* ttbar_lepres_up_h = new TH1F("histo_ttbar_lepres_up","histo_ttbar_lepres_up",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_lepres_up_h, dir+"ttbar"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  ttbar_lepres_up_h->Scale(TopBkgScaleFactor(njets));
  TH1F* tw_lepres_up_h = new TH1F("histo_tw_lepres_up","histo_tw_lepres_up",nbins,minx,maxx);
  fillPlot("bdtg",tw_lepres_up_h, dir+"tw"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  tw_lepres_up_h->Scale(TopBkgScaleFactor(njets));
  TH1F* top_lepres_up_h = new TH1F("histo_Top_CMS_MVALepResBoundingUp","histo_Top_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  top_lepres_up_h->Add(ttbar_lepres_up_h);
  top_lepres_up_h->Add(tw_lepres_up_h);
  TH1F* ttbar_lepres_down_h = new TH1F("histo_ttbar_lepres_down","histo_ttbar_lepres_down",nbins,minx,maxx);
  fillPlot("bdtg",ttbar_lepres_down_h, dir+"ttbar"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  ttbar_lepres_down_h->Scale(TopBkgScaleFactor(njets));
  TH1F* tw_lepres_down_h = new TH1F("histo_tw_lepres_down","histo_tw_lepres_down",nbins,minx,maxx);
  fillPlot("bdtg",tw_lepres_down_h, dir+"tw"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  tw_lepres_down_h->Scale(TopBkgScaleFactor(njets));
  TH1F* top_lepres_down_h = new TH1F("histo_Top_CMS_MVALepResBoundingDown","histo_Top_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  top_lepres_down_h->Add(ttbar_lepres_down_h);
  top_lepres_down_h->Add(tw_lepres_down_h);
  //Wgamma
  TH1F* wg_lepres_up_h = new TH1F("histo_wg_lepres_up","histo_wg_lepres_up",nbins,minx,maxx);
  fillPlot("bdtg",wg_lepres_up_h, dir+"wgamma"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  TH1F* wg3l_lepres_up_h = new TH1F("histo_wg3l_lepres_up","histo_wg3l_lepres_up",nbins,minx,maxx);
  fillPlot("bdtg",wg3l_lepres_up_h, dir+"wg3l"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleUp");
  wg3l_lepres_up_h->Scale(WGstarScaleFactor());
  TH1F* wgamma_lepres_up_h = new TH1F("histo_Wgamma_CMS_MVALepResBoundingUp","histo_Wgamma_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  wgamma_lepres_up_h->Add(wg_lepres_up_h);
  wgamma_lepres_up_h->Add(wg3l_lepres_up_h);
  TH1F* wg_lepres_down_h = new TH1F("histo_wg_lepres_down","histo_wg_lepres_down",nbins,minx,maxx);
  fillPlot("bdtg",wg_lepres_down_h, dir+"wgamma"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  TH1F* wg3l_lepres_down_h = new TH1F("histo_wg3l_lepres_down","histo_wg3l_lepres_down",nbins,minx,maxx);
  fillPlot("bdtg",wg3l_lepres_down_h, dir+"wg3l"+suffix, wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw, "momScaleDown");
  wg3l_lepres_down_h->Scale(WGstarScaleFactor());
  TH1F* wgamma_lepres_down_h = new TH1F("histo_Wgamma_CMS_MVALepResBoundingDown","histo_Wgamma_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  wgamma_lepres_down_h->Add(wg_lepres_down_h);
  wgamma_lepres_down_h->Add(wg3l_lepres_down_h);
  //Ztt
  TH1F* dytt_1_lepres_up_h = new TH1F("histo_dytt_lepres_up_1","histo_dytt_lepres_up_1",nbins,minx,maxx);
  fillPlot("bdtg",dytt_1_lepres_up_h, dir+"data-emb-tau121"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "momScaleUp");
  TH1F* dytt_2_lepres_up_h = new TH1F("histo_dytt_lepres_up_2","histo_dytt_lepres_up_2",nbins,minx,maxx);
  fillPlot("bdtg",dytt_2_lepres_up_h, dir+"data-emb-tau122"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "momScaleUp");
  TH1F* dytt_3_lepres_up_h = new TH1F("histo_dytt_lepres_up_3","histo_dytt_lepres_up_3",nbins,minx,maxx);
  fillPlot("bdtg",dytt_3_lepres_up_h, dir+"data-emb-tau123"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "momScaleUp");
  TH1F* ztt_lepres_up_h = new TH1F("histo_Ztt_CMS_MVALepResBoundingUp","histo_Ztt_CMS_MVALepResBoundingUp",nbins,minx,maxx);
  ztt_lepres_up_h->Add(dytt_1_lepres_up_h);
  ztt_lepres_up_h->Add(dytt_2_lepres_up_h);
  ztt_lepres_up_h->Add(dytt_3_lepres_up_h);
  TH1F* dytt_1_lepres_down_h = new TH1F("histo_dytt_lepres_down_1","histo_dytt_lepres_down_1",nbins,minx,maxx);
  fillPlot("bdtg",dytt_1_lepres_down_h, dir+"data-emb-tau121"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "momScaleDown");
  TH1F* dytt_2_lepres_down_h = new TH1F("histo_dytt_lepres_down_2","histo_dytt_lepres_down_2",nbins,minx,maxx);
  fillPlot("bdtg",dytt_2_lepres_down_h, dir+"data-emb-tau122"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "momScaleDown");
  TH1F* dytt_3_lepres_down_h = new TH1F("histo_dytt_lepres_down_3","histo_dytt_lepres_down_3",nbins,minx,maxx);
  fillPlot("bdtg",dytt_3_lepres_down_h, dir+"data-emb-tau123"+suffix, wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false, "momScaleDown");
  TH1F* ztt_lepres_down_h = new TH1F("histo_Ztt_CMS_MVALepResBoundingDown","histo_Ztt_CMS_MVALepResBoundingDown",nbins,minx,maxx);
  ztt_lepres_down_h->Add(dytt_1_lepres_down_h);
  ztt_lepres_down_h->Add(dytt_2_lepres_down_h);
  ztt_lepres_down_h->Add(dytt_3_lepres_down_h);

  //Write histos to file
  TString outfname = Form("hww%s_%ij.input.root",TString(fs).ReplaceAll("fs","").Data(),njets);
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

  //Other systematics
  qqww_h_up->Write();
  qqww_h_down->Write();
  qqww_h_nlo_up->Write();
  qqww_h_nlo_down->Write();

  top_h_up->Write();
  top_h_down->Write();

  if (fs=="sffs") {
    zjets_h_up->Write();
    zjets_h_down->Write();
  }

  wjets_fr_up_h->Write();
  wjets_fr_down_h->Write();
  wjets_mc_up_h->Write();
  wjets_mc_down_h->Write();
  
  outfile->Close();
  delete rbdtg;
}
