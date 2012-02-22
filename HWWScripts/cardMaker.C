#include "common.C"
#include "fakeBg.C"
#include "Smurf/Analysis/HWWlvlv/HiggsQCDScaleSystematics.h"
#include "Smurf/Analysis/HWWlvlv/PDFgHHSystematics.h"
#include "Smurf/Analysis/HWWlvlv/PSUESystematics.h"

#include "TSystem.h"

void cardMaker(float lumi, int mass, unsigned int njets, TString fs, TString mode, bool saveFile=false) {

  if (fs!="sffs" && fs!="offs" && fs!="" ) {
    cout << "final state not supported" << endl;
    return;
  }

  if (mode!="cut" && mode!="shape") {
    cout << "mode not supported" << endl;
    return;
  }

  if (mode=="shape") doMVA = 1;
  else doMVA = 0;

  if (njets==2) doVBF = 1;
  else doVBF = 0;

  ofstream myfile;
  if (saveFile) {
    TString fname = "hww";
    if (fs=="sffs") fname=fname+"sf_";
    else if (fs=="offs") fname=fname+"of_";
    else fname=fname+"_";
    fname = Form("%s%ij_%s.txt",fname.Data(),njets,mode.Data());
    TString dname = Form("cards/%i/",mass);
    gSystem->Exec("mkdir -p "+dname);
    myfile.open(dname+fname);
  }
  ostream &out = saveFile ? myfile : cout;

  TString dir = main_dir+topww_dir;

  bool useJson = 1;
  bool applyEff=true;
  bool doFake=false; 
  bool doPUw=true;

  TString sigreg = "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20,";

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  pair<float, float> wwSF = make_pair<float, float>(1.0,0.0);
  //for 2j use 1j sf
  unsigned int njetsForWW = njets;
  if (njetsForWW==2) njetsForWW=1;
  if (mass==0) {
    if (mode=="cut") wwSF = make_pair<float, float>(WWBkgScaleFactorCutBased(115,njetsForWW), WWBkgScaleFactorCutBased(115,njetsForWW)*(WWBkgScaleFactorKappaCutBased(115,njetsForWW)-1.));
    else if (mode=="shape") wwSF = make_pair<float, float>(WWBkgScaleFactorMVA(115,njetsForWW), WWBkgScaleFactorMVA(115,njetsForWW)*(WWBkgScaleFactorKappaMVA(115,njetsForWW)-1.));
  } else {
    if (mode=="cut") wwSF = make_pair<float, float>(WWBkgScaleFactorCutBased(mass,njetsForWW), WWBkgScaleFactorCutBased(mass,njetsForWW)*(WWBkgScaleFactorKappaCutBased(mass,njetsForWW)-1.));
    else if (mode=="shape") wwSF = make_pair<float, float>(WWBkgScaleFactorMVA(mass,njetsForWW), WWBkgScaleFactorMVA(mass,njetsForWW)*(WWBkgScaleFactorKappaMVA(mass,njetsForWW)-1.));
  }

  //uncertainty on population of jet bins
  float v_QCDscale_WW    = 1.000;
  float v_QCDscale_WW1in = 1.000;
  float v_QCDscale_WW2in = 1.000;
  if (mass>=200) {
    if (njets==0) {
      v_QCDscale_WW    = 1.042;
      v_QCDscale_WW1in = 0.978;
      v_QCDscale_WW2in = 1.000;
    } else if (njets==1) {
      v_QCDscale_WW    = 1.000;
      v_QCDscale_WW1in = 1.076;
      v_QCDscale_WW2in = 0.914;
    } else if (njets==2) {
      v_QCDscale_WW    = 1.000;
      v_QCDscale_WW1in = 1.000;
      v_QCDscale_WW2in = 1.420;      
    }
  }

  //unceratinty on Higgs cross section at high mass
  float theoryUncXS_HighMH = 1.000;
  if(mass>=200) {
    theoryUncXS_HighMH = 1.0+1.5*(mass/1000.0)*(mass/1000.0)*(mass/1000.0);
  }

  pair<float, float> topSF = make_pair<float, float>(TopBkgScaleFactor(njets), TopBkgScaleFactor(njets)*(TopBkgScaleFactorKappa(njets)-1.));
  //this uncertainty should work only for the dy component
  pair<float, float> dySF = make_pair<float, float>(1., 0);
  if (fs=="sffs") {
    if (mode=="cut") dySF = make_pair<float, float>(DYBkgScaleFactor(mass,njets), DYBkgScaleFactor(mass,njets)*(DYBkgScaleFactorKappa(mass,njets)-1.));
    else if (mode=="shape") dySF = make_pair<float, float>(DYBkgScaleFactor(0,njets), DYBkgScaleFactor(0,njets)*(DYBkgScaleFactorKappa(0,njets)-1.));
  } else dySF = make_pair<float, float>(0., 0);
  if (njets==2) dySF = make_pair<float, float>(DYBkgScaleFactor(0,njets), DYBkgScaleFactor(0,njets)*(DYBkgScaleFactorKappa(0,njets)-1.));

  pair<float, float> data = getYield(dir+"data", wwSelection, veto, mass, njets, sigreg+fs, 0.,   useJson, false, false, false);

  pair<float, float> qqww = getYield(dir+"qqww", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  //pair<float, float> qqww = getYield(dir+"ww_mcnlo", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> ggww = getYield(dir+"ggww", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> ttbar = getYield(dir+"ttbar", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> tw = getYield(dir+"tw", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> zz = getYield(dir+"zz_py", wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wz = getYield(dir+"wz", wwSelection, veto, mass, njets, sigreg+"notZ,"+fs, lumi, useJson, applyEff, doFake, doPUw);

  //add peaking component of VV to DY
  pair<float, float> dyee = getYield(dir+"dyee", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);//fixme should test fromZ
  pair<float, float> dymm = getYield(dir+"dymm", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> pzz = getYield(dir+"zz_py", wwSelection, veto, mass, njets, sigreg+fs+"fromZ", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> pwz = getYield(dir+"wz", wwSelection, veto, mass, njets, sigreg+fs+"fromZ", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> dytt_1 = make_pair<float, float>(0,0);//getYield(dir+"data-emb-tau121", wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  pair<float, float> dytt_2 = make_pair<float, float>(0,0);//getYield(dir+"data-emb-tau122", wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  pair<float, float> dytt_3 = getYield(dir+"data-emb-tau123", wwSelection, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  pair<float, float> dytt = make_pair<float, float>(dytt_1.first+dytt_2.first+dytt_3.first,sqrt(pow(dytt_1.second,2)+pow(dytt_2.second,2)+pow(dytt_3.second,2)));
  
  pair<float, float> wgamma = getYield(dir+"wgamma", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wg3l = getYield(dir+"wg3l", wwSelection, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> wjets = fakeBgEstimation(dir,wwSelectionNoLep, veto, mass, njets, sigreg+fs, lumi,   useJson, applyEff, doPUw);

  pair<float, float> gghww = make_pair<float, float>(0.,0.); 
  pair<float, float> qqhww = make_pair<float, float>(0.,0.); 
  pair<float, float> zhww = make_pair<float, float>(0.,0.); 
  pair<float, float> whww = make_pair<float, float>(0.,0.); 
  if (mass>0) {
    gghww = getYield(dir+Form("hww%i",mass), wwSelection, veto, mass, njets, sigreg+"ggH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
    qqhww = getYield(dir+Form("hww%i",mass), wwSelection, veto, mass, njets, sigreg+"qqH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
    zhww = getYield(dir+Form("hww%i",mass), wwSelection, veto, mass, njets, sigreg+"ZH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
    whww = getYield(dir+Form("hww%i",mass), wwSelection, veto, mass, njets, sigreg+"WH,"+fs, lumi, useJson, applyEff, doFake, doPUw);
  }

  bool makeTable = false;
  if (makeTable) {
    float tot = wwSF.first*qqww.first+ wwSF.first*ggww.first+ zz.first+wz.first+ 
      topSF.first*(ttbar.first+tw.first)+dySF.first*(dymm.first+dyee.first)+pzz.first+pwz.first+ 
      wjets.first+ wgamma.first+wg3l.first*WGstarScaleFactor()+dytt.first; 
    float err = sqrt(pow(wwSF.first*qqww.second,2)+pow(wwSF.first*ggww.second,2)+pow(sqrt(pow(zz.second,2)+pow(wz.second,2)),2)+pow(topSF.first*sqrt(pow(ttbar.second,2)+pow(tw.second,2)),2)+
		     pow(sqrt(pow(dySF.first*dymm.second,2)+pow(dySF.first*dyee.second,2)+pow(pzz.second,2)+pow(pwz.second,2)),2)+pow(wjets.second,2)+
		     pow(sqrt(pow(wgamma.second,2)+pow(wg3l.second,2)),2)+pow(dytt.second,2));
    cout << Form(" %6.0f $\\pm$ %6.0f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f \\\\",
		 data.first,data.second,
		 tot,err,
		 wwSF.first*qqww.first, wwSF.first*qqww.second,
		 wwSF.first*ggww.first, wwSF.first*ggww.second,
		 topSF.first*(ttbar.first+tw.first), topSF.first*sqrt(pow(ttbar.second,2)+pow(tw.second,2)),
		 wjets.first, wjets.second,
		 zz.first+wz.first, sqrt(pow(zz.second,2)+pow(wz.second,2)),
		 dySF.first*(dymm.first+dyee.first)+pzz.first+pwz.first, sqrt(pow(dySF.first*dymm.second,2)+pow(dySF.first*dyee.second,2)+pow(pzz.second,2)+pow(pwz.second,2)),
		 wgamma.first+wg3l.first*WGstarScaleFactor(), sqrt(pow(wgamma.second,2)+pow(wg3l.second,2)),
		 dytt.first, dytt.second) << endl;
  } else {

    //remove "fs"
    fs.ReplaceAll("fs","");
    if (fs=="") fs="ll";

    TString bn = Form("j%i%s",njets,fs.Data());

    out << Form("%s \n","imax 1 number of channels");
    out << Form("%s \n","jmax * number of background");
    out << Form("%s \n","kmax * number of nuisance parameters");
    out << Form("%s %i\n","Observation",((int)data.first));
    if (mode=="shape") {
      out << Form("shapes *   *   hww%s_%ij.input.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC \n",fs.Data(),njets);
      out << Form("shapes data_obs * hww%s_%ij.input.root  histo_Data \n",fs.Data(),njets);
    }
    out << Form("%-35s %5s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s \n","bin","",bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data());
    out << Form("%-35s %5s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s \n","process","","ZH","WH","qqH","ggH","qqWW","ggWW","VV","Top","Zjets","Wjets","Wgamma","Ztt");
    out << Form("%-35s %5s %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i \n","process","",-3,-2,-1,0,1,2,3,4,5,6,7,8);
    out << Form("%-35s %5s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n","rate","",
	   zhww.first, whww.first, qqhww.first, gghww.first, wwSF.first*qqww.first, wwSF.first*ggww.first, zz.first+wz.first, topSF.first*(ttbar.first+tw.first), 
	   dySF.first*(dymm.first+dyee.first)+pzz.first+pwz.first, wjets.first, wgamma.first+wg3l.first*WGstarScaleFactor(), dytt.first);
    out << Form("%-35s %5s   1.045   1.045   1.045   1.045     -       -     1.045     -       -       -     1.045   1.045\n","lumi","lnN");
    if (mode=="shape" && doResEffSyst) {
      //we do not have WH and ZH for mH=115
      out << Form("%-35s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000     -       -       -     %5s   %5s \n","CMS_MVALepEffBounding","shape",zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",(wgamma.first+wg3l.first)>0?"1.000":"  -  ",dytt.first>0?"1.000":"  -  ");//why no top???
      out << Form("%-35s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000   1.000     -       -     %5s   %5s \n","CMS_MVALepResBounding","shape",zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",(wgamma.first+wg3l.first)>0?"1.000":"  -  ",dytt.first>0?"1.000":"  -  ");
      out << Form("%-35s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000   1.000     -       -     %5s   %5s \n","CMS_MVAMETResBounding","shape",zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",(wgamma.first+wg3l.first)>0?"1.000":"  -  ",dytt.first>0?"1.000":"  -  ");
      out << Form("%-35s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000   1.000     -       -     %5s   %5s \n","CMS_MVAJESBounding","shape",zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",(wgamma.first+wg3l.first)>0?"1.000":"  -  ",dytt.first>0?"1.000":"  -  ");
    } else {
      out << Form("%-35s %5s   1.030   1.030   1.030   1.030   1.030   1.030   1.030     -       -       -     1.030   1.030\n","CMS_eff_m","lnN");
      out << Form("%-35s %5s   1.040   1.040   1.040   1.040   1.040   1.040   1.040     -       -       -     1.040   1.040\n","CMS_eff_e","lnN");
      out << Form("%-35s %5s   1.015   1.015   1.015   1.015   1.015   1.015   1.015     -       -       -     1.015   1.015\n","CMS_scale_m","lnN");
      out << Form("%-35s %5s   1.020   1.020   1.020   1.020   1.020   1.020   1.020     -       -       -     1.020   1.020\n","CMS_scale_e","lnN");
      out << Form("%-35s %5s   1.020   1.020   1.020   1.020   1.020   1.020   1.020     -       -       -     1.020   1.020\n","CMS_hww_met_resolution","lnN");
      if (njets==0) out << Form("%-35s %5s   1.020   1.020   1.020   1.020   1.020   1.020   1.020     -       -       -     1.020   1.020\n","CMS_scale_j","lnN");
      if (njets==1) out << Form("%-35s %5s   1.050   1.050   1.050   1.050   1.050   1.050   1.050     -       -       -     1.050   1.050\n","CMS_scale_j","lnN");
      if (njets==2) out << Form("%-35s %5s   1.100   1.100   1.100   1.100   1.100   1.100   1.100     -       -       -     1.100   1.100\n","CMS_scale_j","lnN");
    }
    out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -     1.360     -       -  \n","FakeRate","lnN");
    if (mode=="shape") {
      out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -     1.000     -       -  \n","CMS_MVAWBounding_hww","shape");
      out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -     1.000     -       -  \n","CMS_MVAWMCBounding_hww","shape");      
      out << Form("%-35s %5s     -       -       -     1.000     -       -       -       -       -       -       -       -  \n","CMS_MVAggHBounding","shape");      
    }
    out << Form("%-35s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -  \n","UEPS","lnN", 
		mass>0 ? HiggsSignalPSUESystematics(mass, njets) : 0.);
    out << Form("%-35s %5s   %5.3f   %5.3f   %5.3f   %5.3f     -       -       -       -       -       -       -       -  \n","theoryUncXS_HighMH","lnN",theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH);
    out << Form("%-35s %5s     -       -       -     %5.3f     -     1.040     -       -       -       -       -       -  \n","pdf_gg","lnN",PDFgHHSystematics(mass));
    out << Form("%-35s %5s   1.050   1.050   1.050     -     1.040     -     1.040     -       -       -     1.040   1.040\n","pdf_qqbar","lnN");
    out << Form("%-35s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -  \n","QCDscale_ggH","lnN", 
		mass>0 ? HiggsSignalQCDScaleKappa("QCDscale_ggH",mass, njets) : 0.);
    out << Form("%-35s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -  \n","QCDscale_ggH1in","lnN", 
		mass>0 ? HiggsSignalQCDScaleKappa("QCDscale_ggH1in",mass, njets) : 0.);
    out << Form("%-35s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -  \n","QCDscale_ggH2in","lnN", 
		mass>0 ? HiggsSignalQCDScaleKappa("QCDscale_ggH2in",mass, njets) : 0.);
    out << Form("%-35s %5s     -       -     1.010     -       -       -       -       -       -       -       -       -  \n","QCDscale_qqH","lnN");
    out << Form("%-35s %5s   1.020   1.020     -       -       -       -       -       -       -       -       -       -  \n","QCDscale_VH","lnN");
    out << Form("%-35s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -  \n","QCDscale_WW","lnN", 
	        v_QCDscale_WW );
    out << Form("%-35s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -  \n","QCDscale_WW1in","lnN", 
		v_QCDscale_WW1in );
    out << Form("%-35s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -  \n","QCDscale_WW2in","lnN", 
		v_QCDscale_WW2in );
    out << Form("%-35s %5s     -       -       -       -       -       -     1.040     -       -       -       -       -  \n","QCDscale_VV","lnN");
    out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -       -     1.300     -  \n","QCDscale_V","lnN");
    out << Form("%-35s %5s     -       -       -       -       -     1.300     -       -       -       -       -       -  \n","QCDscale_ggVV","lnN");
    out << Form("%-35s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -  \n","QCDscale_WW_EXTRAP","lnN",
		1.060);//this is the unceratinty for extrapolation from sideband to signal region
    out << Form("%-35s %5s     -       -       -     1.020     -       -       -       -       -       -       -       -  \n","QCDscale_ggH_ACEPT","lnN");
    out << Form("%-35s %5s     -       -     1.020     -       -       -       -       -       -       -       -       -  \n","QCDscale_qqH_ACEPT","lnN");
    out << Form("%-35s %5s   1.020   1.020     -       -       -       -       -       -       -       -       -       -  \n","QCDscale_VH_ACEPT","lnN");
    out << Form("%-35s %5s     -       -       -       -       -       -       -     %5.3f     -       -       -       -  \n",Form("CMS_hww_%ij_ttbar",njets),"lnN",
		1.+topSF.second/topSF.first);
    out << Form("%-35s %5s     -       -       -       -       -       -       -       -     %5.3f     -       -       -  \n",Form("CMS_hww%s_%ij_Z",fs.Data(),njets),"lnN",
		(1.+dySF.second*(dymm.first+dyee.first)/(dySF.first*(dymm.first+dyee.first)+pzz.first+pwz.first)));
    out << Form("%-35s %5s     -       -       -       -     %5.3f   %5.3f     -       -       -       -       -       -  \n",Form("CMS_hww_%ij_WW",njets),"lnN",
		1.+wwSF.second/wwSF.first,1.+wwSF.second/wwSF.first);
    out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -       -       -     %5.3f\n","CMS_hww_Ztt","lnN",
		ZttScaleFactorKappa());

    if (mode=="cut") {
      if (zhww.first>0.) out << Form("%-35s %5s   %5.3f     -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ZH",fs.Data(),njets),"lnN",
				     zhww.first>0 ? 1.+zhww.second/zhww.first : 1.);
      if (whww.first>0.) out << Form("%-35s %5s     -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_WH",fs.Data(),njets),"lnN",
				     whww.first>0 ? 1.+whww.second/whww.first : 1.);
      out << Form("%-35s %5s     -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_qqH",fs.Data(),njets),"lnN",
		  qqhww.first>0 ? 1.+qqhww.second/qqhww.first : 1.);
      out << Form("%-35s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ggH",fs.Data(),njets),"lnN",
		  gghww.first>0 ? 1.+gghww.second/gghww.first : 1.);
      out << Form("%-35s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_WW",fs.Data(),njets),"lnN",
		  1.+qqww.second/qqww.first);
      out << Form("%-35s %5s     -       -       -       -       -     %5.3f     -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ggWW",fs.Data(),njets),"lnN",
		  1.+ggww.second/ggww.first);
      out << Form("%-35s %5s     -       -       -       -       -       -     %5.3f     -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_VV",fs.Data(),njets),"lnN",
		  1.+sqrt(pow(zz.second,2)+pow(wz.second,2))/(zz.first+wz.first));
      out << Form("%-35s %5s     -       -       -       -       -       -       -     %5.3f     -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ttbar",fs.Data(),njets),"lnN",
		  1.+sqrt(pow(ttbar.second,2)+pow(tw.second,2))/(ttbar.first+tw.first));
      out << Form("%-35s %5s     -       -       -       -       -       -       -       -     %5.3f     -       -       -  \n",Form("CMS_hww%s_stat_%ij_Z",fs.Data(),njets),"lnN",
		  1.+sqrt(pow(dySF.first*dymm.second,2)+pow(dySF.first*dyee.second,2)+pow(pzz.second,2)+pow(pwz.second,2))/(dySF.first*dymm.first+dySF.first*dyee.first+pzz.first+pwz.first));
      out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -     %5.3f     -       -  \n",Form("CMS_hww%s_stat_%ij_Wjets",fs.Data(),njets),"lnN",
		  wjets.first>0 ? 1.+wjets.second/wjets.first : 1.);
      out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -       -     %5.3f     -  \n",Form("CMS_hww%s_stat_%ij_Wgamma",fs.Data(),njets),"lnN",
		  (wgamma.first+wg3l.first)>0 ? 1.+sqrt(pow(wgamma.second,2)+pow(wg3l.second,2))/(wgamma.first+wg3l.first) : 1.);
      out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -       -       -     %5.3f\n",Form("CMS_hww%s_stat_%ij_Ztt",fs.Data(),njets),"lnN",
		  dytt.first>0 ? 1.+dytt.second/dytt.first : 1.);
    } else if (mode=="shape") {
      if (fs=="sf") out << Form("%-35s %5s     -       -       -       -       -       -       -       -     1.000     -       -       -  \n",Form("CMS_MVAZBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      out << Form("%-35s %5s     -       -       -       -       -       -       -     1.000     -       -       -       -  \n","CMS_MVATopBounding_hww","shape");      
      out << Form("%-35s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -  \n","CMS_MVAWWBounding_hww","shape");      
      out << Form("%-35s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -  \n","CMS_MVAWWNLOBounding_hww","shape");      
      if (zhww.first>0.)
	out << Form("%-35s %5s   1.000     -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_MVAZHStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (whww.first>0.)
	out << Form("%-35s %5s     -     1.000     -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_MVAWHStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (qqhww.first>0.)
	out << Form("%-35s %5s     -       -     1.000     -       -       -       -       -       -       -       -       -  \n",Form("CMS_MVAqqHStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (gghww.first>0.)
	out << Form("%-35s %5s     -       -       -     1.000     -       -       -       -       -       -       -       -  \n",Form("CMS_MVAggHStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (qqww.first>0.)
	out << Form("%-35s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -  \n",Form("CMS_MVAqqWWStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (ggww.first>0.)
	out << Form("%-35s %5s     -       -       -       -       -     1.000     -       -       -       -       -       -  \n",Form("CMS_MVAggWWStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if ((zz.first+wz.first)>0.)
	out << Form("%-35s %5s     -       -       -       -       -       -     1.000     -       -       -       -       -  \n",Form("CMS_MVAVVStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if ((ttbar.first+tw.first)>0.)
	out << Form("%-35s %5s     -       -       -       -       -       -       -     1.000     -       -       -       -  \n",Form("CMS_MVATopStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if ((dySF.first*dymm.first+dySF.first*dyee.first+pzz.first+pwz.first)>0.)
	out << Form("%-35s %5s     -       -       -       -       -       -       -       -     1.000     -       -       -  \n",Form("CMS_MVAZjetsStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (wjets.first>0.)
	out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -     1.000     -       -  \n",Form("CMS_MVAWjetsStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if ((wgamma.first+wg3l.first)>0) 
	out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -       -     1.000     -  \n",Form("CMS_MVAWgammaStatBounding_hww%s_%ij",fs.Data(),njets),"shape");      
      if (dytt.first>0) 
	out << Form("%-35s %5s     -       -       -       -       -       -       -       -       -       -       -     1.000\n",Form("CMS_MVAZttStatBounding_hww%s_%ij",fs.Data(),njets),"shape");
    }
  }

  if (saveFile) {
    myfile.close();
  }

  doMVA = 0;
  doVBF = 0;

}

void cardMaker(float lumi, TString mode) {

  //int jets[] = {1};
  int jets[] = {0,1,2};

  //int masses[] = {120};
  int masses[] = {115,120,130,140,150,160,170,180,190,200,250,300};
  int nmasses = sizeof(masses)/sizeof(int);
  int njets = sizeof(jets)/sizeof(int);
  for (int j=0;j<nmasses;++j) {
    for (int jj=0;jj<njets;++jj) {
      int mass = masses[j];
      int jetbin = jets[jj];
      if (jetbin!=2) {
	cardMaker(lumi,mass,jetbin,"sffs",mode,true);
	cardMaker(lumi,mass,jetbin,"offs",mode,true);
      } else {
	cardMaker(lumi,mass,jetbin,"","cut",true);
      }
    }
  }

}


