#include "common.C"
#include "./fakeBg.C"
#include "../Smurf/Analysis/HWWlvlv/HiggsQCDScaleSystematics_7TeV.h"
#include "../Smurf/Analysis/HWWlvlv/PDFgHHSystematics_7TeV.h"
#include "../Smurf/Analysis/HWWlvlv/PSUESystematics_7TeV.h"

#include "TSystem.h"

void cardMaker(float lumi, int mass, unsigned int njets, TString fs, TString mode, bool saveFile=false) {

  bool inj125 = 0;

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
    fname = Form("%s%ij_%s_7TeV.txt",fname.Data(),njets,mode.Data());
    TString dname = Form("../cards/%i/",mass);
    gSystem->Exec("mkdir -p "+dname);
    myfile.open(dname+fname);
  }
  ostream &out = saveFile ? myfile : cout;

  TString dir   = main_dir+topww_dir;
  TString dirwj = main_dir+wj_dir;

  bool useJson = 1;
  bool applyEff=true;
  bool doFake=false; 
  bool doPUw=true;

  TString sigreg = "=dphireg=dphijet=minmetvtx=lep2pt15=ptll3045=mll20=";

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  pair<float, float> wwSF = make_pair<float, float>(1.0,0.0);
  int massForWW = TMath::Min(TMath::Max((int)mass,115),200);
  if (massForWW==125) massForWW=126;

  if (njets==2) {
    wwSF = make_pair<float, float>(WWVBFScaleFactor()*WWBkgScaleFactorCutBased(massForWW,1),1.0);
  } else if (mass==0) {
    if (mode=="cut") wwSF = make_pair<float, float>(WWBkgScaleFactorCutBased(115,njets), WWBkgScaleFactorCutBased(115,njets)*(WWBkgScaleFactorKappaCutBased(115,njets)-1.));
    else if (mode=="shape") wwSF = make_pair<float, float>(WWBkgScaleFactorMVA(115,njets), WWBkgScaleFactorMVA(115,njets)*(WWBkgScaleFactorKappaMVA(115,njets)-1.));
  } else {
    if (mode=="cut") wwSF = make_pair<float, float>(WWBkgScaleFactorCutBased(massForWW,njets), WWBkgScaleFactorCutBased(massForWW,njets)*(WWBkgScaleFactorKappaCutBased(massForWW,njets)-1.));
    else if (mode=="shape") wwSF = make_pair<float, float>(WWBkgScaleFactorMVA(massForWW,njets), WWBkgScaleFactorMVA(massForWW,njets)*(WWBkgScaleFactorKappaMVA(massForWW,njets)-1.));
  }
  //try this to sync with guillelmo
  if (njets==2) wwSF.second = 0.5*wwSF.first;
  else if (mass>200 && doMVA==0) {
    wwSF.second = 0.0;
  }

  int mH = mass;
  if (mH==125) mH=126;
  //different in case of new signal injection test
  if (inj125) mH=126;

  //uncertainty on population of jet bins
  float v_QCDscale_WW    = 1.000;
  float v_QCDscale_WW1in = 1.000;
  float v_QCDscale_WW2in = 1.000;
  if (mass>200) {
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
  if (mH>200) theoryUncXS_HighMH = 1.0+1.5*(mass/1000.0)*(mass/1000.0)*(mass/1000.0);//needed for 7TeV

  pair<float, float> topSF = make_pair<float, float>(TopBkgScaleFactor(njets), TopBkgScaleFactor(njets)*(TopBkgScaleFactorKappa(njets)-1.));
  //this uncertainty should work only for the dy component
  pair<float, float> dySF = make_pair<float, float>(1., 0);
  if (fs=="sffs") {
    if (mode=="cut") dySF = make_pair<float, float>(DYBkgScaleFactor(max(115,mass),njets), DYBkgScaleFactor(max(115,mass),njets)*(DYBkgScaleFactorKappa(max(115,mass),njets)-1.));
    else if (mode=="shape") dySF = make_pair<float, float>(DYBkgScaleFactor(0,njets), DYBkgScaleFactor(0,njets)*(DYBkgScaleFactorKappa(0,njets)-1.));
  } else dySF = make_pair<float, float>(1.0, 0);
  //fixme error from mass dependent do be closer to Guillelmo...
  if (njets==2) dySF = make_pair<float, float>(DYBkgScaleFactor(0,njets), DYBkgScaleFactor(0,njets)*(DYBkgScaleFactorKappa(max(115,mass),njets)-1.));
  fs=fs+"=";

  pair<float, float> data = getYield(dir+"data", wwSelNoMet, veto, mass, njets, sigreg+fs, 0.,   useJson, false, false, false);

  pair<float, float> qqww = getYield(dir+"qqww", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> ggww = getYield(dir+"ggww", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> ttbar = getYield(dir+"ttbar", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> tw = getYield(dir+"tw", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> zz = getYield(dir+"zz_py", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wz = getYield(dir+"wz", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> www = getYield(dir+"www", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  //now compute DY yield and K
  float dyY = 0.0;
  float dyMCE = 0.0;
  float dyDDK = 0.0;//for of case
  pair<float, float> dymm = getYield(dir+"dymm", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);//+"=fromZ="
  pair<float, float> dyee = getYield(dir+"dyee", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);//+"=fromZ="
  pair<float, float> dyll = make_pair<float, float>(dymm.first+dyee.first,sqrt(pow(dymm.second,2)+pow(dyee.second,2)));
  //ok, in this case dySF means the scale factor
  dyY = dySF.first*dyll.first;
  dyMCE = dySF.first*dyll.second;
  if (fs.Contains("offs")) {
    dyDDK = 1.00;
  } else {
    if (dyll.first>0.) dyDDK = 1.+dySF.second/dySF.first;
    else dyDDK = 1.0;
  }

  //pair<float, float> dytt_1;//= getYield(dir+"data-emb-tau121", wwSelNoMet, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  //pair<float, float> dytt_2;//= getYield(dir+"data-emb-tau122", wwSelNoMet, veto, mass, njets, sigreg+"embed,"+fs, lumi, false, false, false, false);
  //pair<float, float> dytt_3 = getYield(dir+"data-emb-tau123", wwSelNoMet, veto, mass, njets, sigreg+"=embed="+fs, lumi, false, false, false, false);
  //pair<float, float> dytt = make_pair<float, float>(dytt_1.first+dytt_2.first+dytt_3.first,sqrt(pow(dytt_1.second,2)+pow(dytt_2.second,2)+pow(dytt_3.second,2)));
  pair<float, float> dytt = getYield(dir+"dytt", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  
  pair<float, float> wgamma = getYield(dir+"wgamma", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wg3l = getYield(dir+"wg3l", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zgamma;//= getYield(dir+"zgamma", wwSelNoMet, veto, mass, njets, sigreg+fs, lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> wjets_e = fakeBgEstimation(dirwj,wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=elfake=", lumi,   useJson, applyEff, doPUw);
  pair<float, float> wjets_m = fakeBgEstimation(dirwj,wwSelNoMetNoLep, veto, mass, njets, sigreg+fs+"=mufake=", lumi,   useJson, applyEff, doPUw);

  pair<float, float> gghww = make_pair<float, float>(0.,0.); 
  pair<float, float> qqhww = make_pair<float, float>(0.,0.); 
  pair<float, float> zhww = make_pair<float, float>(0.,0.); 
  pair<float, float> whww = make_pair<float, float>(0.,0.); 
  if (mass>0) {
      gghww = getYield(dir+Form("hww%i",mH), wwSelNoMet, veto, mass, njets, sigreg+"=ggH="+fs, lumi, useJson, applyEff, doFake, doPUw);
      qqhww = getYield(dir+Form("hww%i",mH), wwSelNoMet, veto, mass, njets, sigreg+"=qqH="+fs, lumi, useJson, applyEff, doFake, doPUw);
      zhww = getYield(dir+Form("hww%i",mH), wwSelNoMet, veto, mass, njets, sigreg+"=ZH="+fs, lumi, useJson, applyEff, doFake, doPUw);
      whww = getYield(dir+Form("hww%i",mH), wwSelNoMet, veto, mass, njets, sigreg+"=WH="+fs, lumi, useJson, applyEff, doFake, doPUw);
      //gghww = make_pair<float, float>(InterfgHHSystematics(mass)*gghww.first,InterfgHHSystematics(mass)*gghww.second);
    if (mass==125 || inj125) {
      qqhww = make_pair<float, float>(qqhww.first*(1.211*2.16)/(1.199*2.33)  ,qqhww.second*(1.211*2.16)/(1.199*2.33)  ); 
      gghww = make_pair<float, float>(gghww.first*(15.31*2.16)/(15.06*2.33)  ,gghww.second*(15.31*2.16)/(15.06*2.33)  ); 
      whww = make_pair<float, float> (whww.first *(0.5729*2.16)/(0.5576*2.33),whww.second *(0.5729*2.16)/(0.5576*2.33)); 
      zhww = make_pair<float, float> (zhww.first *(0.3158*2.16)/(0.3077*2.33),zhww.second *(0.3158*2.16)/(0.3077*2.33)); 
    }
  }

  if (1) {
    //reduce the number of channels by merging WH and ZH
    whww = make_pair<float, float>(whww.first+zhww.first,sqrt(whww.first*whww.first+zhww.first*zhww.first));
    zhww = make_pair<float, float>(0.,0.);
  }

  bool makeTable = false;
  if (makeTable) {
    /*
    float tot = wwSF.first*qqww.first+ wwSF.first*ggww.first+ zz.first+wz.first+www.first+
      topSF.first*(ttbar.first+tw.first)+dyY+pzz.first+pwz.first+ 
      wjets.first+ wgamma.first+zgamma.first+wg3l.first+dytt.first; 
    float err = sqrt(pow(wwSF.first*qqww.second,2)+pow(wwSF.first*ggww.second,2)+pow(sqrt(pow(zz.second,2)+pow(wz.second,2)+pow(www.second,2)),2)+pow(topSF.first*sqrt(pow(ttbar.second,2)+pow(tw.second,2)),2)+
		     pow(sqrt(pow(dyY*(dyDDK-1.),2)+pow(pzz.second,2)+pow(pwz.second,2)),2)+pow(wjets.second,2)+
		     pow(sqrt(pow(wgamma.second,2)+pow(zgamma.second,2)+pow(wg3l.second,2)),2)+pow(dytt.second,2));
    cout << Form(" %6.0f $\\pm$ %6.0f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f & %6.1f $\\pm$ %6.1f \\\\",
		 data.first,data.second,
		 tot,err,
		 wwSF.first*qqww.first, wwSF.first*qqww.second,
		 wwSF.first*ggww.first, wwSF.first*ggww.second,
		 topSF.first*(ttbar.first+tw.first), topSF.first*sqrt(pow(ttbar.second,2)+pow(tw.second,2)),
		 wjets.first, wjets.second,
		 zz.first+wz.first+www.first, sqrt(pow(zz.second,2)+pow(wz.second,2)+pow(www.second,2)),
		 dyY+pzz.first+pwz.first, sqrt(pow(dyY*(dyDDK-1.),2)+pow(pzz.second,2)+pow(pwz.second,2)),
		 wgamma.first+zgamma.first+wg3l.first, sqrt(pow(wgamma.second,2)+pow(zgamma.second,2)+pow(wg3l.second,2)),
		 dytt.first, dytt.second) << endl;
    */
  } else {

    //remove "fs" and "+"
    fs.ReplaceAll("=","");
    fs.ReplaceAll("fs","");
    if (fs=="") fs="ll";

    TString bn = Form("j%i%s",njets,fs.Data());

    out << Form("%s \n","imax 1 number of channels");
    out << Form("%s \n","jmax * number of background");
    out << Form("%s \n","kmax * number of nuisance parameters");
    out << Form("%s %i\n","Observation",((int)data.first));
    if (mode=="shape") {
      out << Form("shapes *   *   hww%s_%ij.input_7TeV.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC \n",fs.Data(),njets);
      out << Form("shapes data_obs * hww%s_%ij.input_7TeV.root  histo_Data \n",fs.Data(),njets);
    }
    out << Form("%-40s %5s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s \n","bin","",bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data(),bn.Data());
    out << Form("%-40s %5s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s \n","process","","ZH","WH","qqH","ggH","qqWW","ggWW","VV","Top","Zjets","WjetsE","Wgamma","Wg3l","Ztt","WjetsM");
    out << Form("%-40s %5s %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i \n","process","",-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10);
    out << Form("%-40s %5s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n","rate","",
	   zhww.first, whww.first, qqhww.first, gghww.first, wwSF.first*qqww.first, wwSF.first*ggww.first, zz.first+wz.first+www.first, topSF.first*(ttbar.first+tw.first), 
		dyY, wjets_e.first, wgamma.first+zgamma.first,wg3l.first, dytt.first, wjets_m.first);
    if (mass>200) {
      out << Form("%-40s %5s   1.022   1.022   1.022   1.022   1.022   1.022   1.022     -       -       -     1.022   1.022   1.022     -  \n","lumi_7TeV","lnN");
    } else {
      out << Form("%-40s %5s   1.022   1.022   1.022   1.022     -       -     1.022     -       -       -     1.022   1.022   1.022     -  \n","lumi_7TeV","lnN");
    }
    if (mode=="shape" && doResEffSyst) {
      //do not consider wgamma here, too few MC events just random results (30% syst covers)
      out << Form("%-40s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000     -       -       -     %5s   %5s     -       -  \n","CMS_hww_MVALepEffBounding","shape",
		  zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",0*(wgamma.first+zgamma.first)>0?"1.000":"  -  ",0*(wg3l.first)>0?"1.000":"  -  ");//why no top???
      out << Form("%-40s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000   1.000     -       -     %5s   %5s     -       -  \n","CMS_hww_MVALepResBounding","shape",
		  zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",0*(wgamma.first+zgamma.first)>0?"1.000":"  -  ",0*(wg3l.first)>0?"1.000":"  -  ");
      out << Form("%-40s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000   1.000     -       -     %5s   %5s     -       -  \n","CMS_hww_MVAMETResBounding","shape",
		  zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",0*(wgamma.first+zgamma.first)>0?"1.000":"  -  ",0*(wg3l.first)>0?"1.000":"  -  ");
      out << Form("%-40s %5s   %5s   %5s   1.000   1.000   1.000   1.000   1.000   1.000     -       -     %5s   %5s     -       -  \n","CMS_hww_MVAJESBounding","shape",
		  zhww.first>0?"1.000":"  -  ",whww.first>0?"1.000":"  -  ",0*(wgamma.first+zgamma.first)>0?"1.000":"  -  ",0*(wg3l.first)>0?"1.000":"  -  ");
    } else {
      if (mass>200) {
	out << Form("%-40s %5s   1.030   1.030   1.030   1.030   1.030   1.030   1.030     -       -       -     1.030   1.030   1.030     -  \n","CMS_eff_m","lnN");
	out << Form("%-40s %5s   1.040   1.040   1.040   1.040   1.040   1.040   1.040     -       -       -     1.040   1.040   1.040     -  \n","CMS_eff_e","lnN");
	out << Form("%-40s %5s   1.015   1.015   1.015   1.015   1.015   1.015   1.015     -       -       -     1.015   1.015   1.015     -  \n","CMS_scale_m","lnN");
	out << Form("%-40s %5s   1.020   1.020   1.020   1.020   1.020   1.020   1.020     -       -       -     1.020   1.020   1.020     -  \n","CMS_scale_e","lnN");
	out << Form("%-40s %5s   1.020   1.020   1.020   1.020   1.020   1.020   1.020     -       -       -     1.020   1.020   1.020     -  \n","CMS_hww_met_resolution","lnN");
      } else {
       out << Form("%-40s %5s   1.030   1.030   1.030   1.030     -       -     1.030     -       -       -     1.030   1.030   1.030     -  \n","CMS_eff_m","lnN");
       out << Form("%-40s %5s   1.040   1.040   1.040   1.040     -       -     1.040     -       -       -     1.040   1.040   1.040     -  \n","CMS_eff_e","lnN");
       out << Form("%-40s %5s   1.015   1.015   1.015   1.015     -       -     1.015     -       -       -     1.015   1.015   1.015     -  \n","CMS_scale_m","lnN");
       out << Form("%-40s %5s   1.020   1.020   1.020   1.020     -       -     1.020     -       -       -     1.020   1.020   1.020     -  \n","CMS_scale_e","lnN");
       out << Form("%-40s %5s   1.020   1.020   1.020   1.020     -       -     1.020     -       -       -     1.020   1.020   1.020     -  \n","CMS_hww_met_resolution","lnN");
      }
      if (njets==0) out << Form("%-40s %5s   1.020   1.020   1.020   1.020   1.020   1.020   1.020     -       -       -     1.020   1.020   1.020     -  \n","CMS_scale_j","lnN");
      if (njets==1) out << Form("%-40s %5s   1.050   1.050   1.050   1.050   1.050   1.050   1.050     -       -       -     1.050   1.050   1.050     -  \n","CMS_scale_j","lnN");
      if (njets==2) out << Form("%-40s %5s   1.100   1.100   1.100   1.100   1.100   1.100   1.100     -       -       -     1.100   1.100   1.100     -  \n","CMS_scale_j","lnN");
    }
    if (mode=="shape") {
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -     1.360     -       -       -       -  \n","FakeRate_e","lnU");
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -       -     1.360\n","FakeRate_m","lnU");
    } else {
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -     1.360     -       -       -       -  \n","FakeRate_cut_e","lnN");
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -       -     1.360\n","FakeRate_cut_m","lnN");
    }
    if (mode=="shape") {
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -     1.000     -       -       -       -  \n","CMS_hww_MVAWEBounding","shape");
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -       -     1.000\n","CMS_hww_MVAWMBounding","shape");
    }
    out << Form("%-40s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n","UEPS","lnN",
		mass>0 ? HiggsSignalPSUESystematics(mH, njets) : 0.);
    out << Form("%-40s %5s   %5.3f   %5.3f   %5.3f   %5.3f     -       -       -       -       -       -       -       -       -       -  \n","theoryUncXS_HighMH","lnN",theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH,theoryUncXS_HighMH);
    out << Form("%-40s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n","pdf_gg","lnN",PDFgHHSystematics(mH));
    out << Form("%-40s %5s   1.050   1.050   1.050     -       -       -     1.040     -       -       -     1.040   1.040     -       -  \n","pdf_qqbar","lnN");
    out << Form("%-40s %5s     -       -       -       -       -     1.000     -       -       -       -       -       -       -       -  \n","CMS_hww_PDFggWW","shape");
    out << Form("%-40s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -       -       -  \n","CMS_hww_PDFqqWW","shape");
    out << Form("%-40s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n","QCDscale_ggH","lnN", 
                mass>0 ? HiggsSignalQCDScaleKappa("QCDscale_ggH",mH, njets) : 0.);
    out << Form("%-40s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n","QCDscale_ggH1in","lnN", 
                mass>0 ? HiggsSignalQCDScaleKappa("QCDscale_ggH1in",mH, njets) : 0.);
    out << Form("%-40s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n","QCDscale_ggH2in","lnN", 
                mass>0 ? HiggsSignalQCDScaleKappa("QCDscale_ggH2in",mH, njets) : 0.);
    out << Form("%-40s %5s     -       -     1.010     -       -       -       -       -       -       -       -       -       -       -  \n","QCDscale_qqH","lnN");
    out << Form("%-40s %5s   1.020   1.020     -       -       -       -       -       -       -       -       -       -       -       -  \n","QCDscale_VH","lnN");
    if (mode!="shape") {
      out << Form("%-40s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n","QCDscale_WW","lnN", 
                  v_QCDscale_WW );
      out << Form("%-40s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n","QCDscale_WW1in","lnN", 
                  v_QCDscale_WW1in );
      out << Form("%-40s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n","QCDscale_WW2in","lnN", 
                  v_QCDscale_WW2in );
    }
    out << Form("%-40s %5s     -       -       -       -       -       -     1.040     -       -       -       -       -       -       -  \n","QCDscale_VV","lnN");
    out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -     1.300     -       -       -  \n","QCDscale_Vgamma","lnN");
    out << Form("%-40s %5s     -       -       -       -       -     1.300     -       -       -       -       -       -       -       -  \n","QCDscale_ggVV","lnN");
    if (mode!="shape") {
      out << Form("%-40s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n","QCDscale_WW_EXTRAP","lnN",
                  1.060);//this is the unceratinty for extrapolation from sideband to signal region
    }
    out << Form("%-40s %5s     -       -       -     1.020     -       -       -       -       -       -       -       -       -       -  \n","QCDscale_ggH_ACCEPT","lnN");
    out << Form("%-40s %5s     -       -     1.020     -       -       -       -       -       -       -       -       -       -       -  \n","QCDscale_qqH_ACCEPT","lnN");
    out << Form("%-40s %5s   1.020   1.020     -       -       -       -       -       -       -       -       -       -       -       -  \n","QCDscale_VH_ACCEPT","lnN");
    out << Form("%-40s %5s     -       -       -       -       -       -       -     %5.3f     -       -       -       -       -       -  \n",Form("CMS_hww_%ij_ttbar_7TeV",njets),"lnN",
		1.+topSF.second/topSF.first);
    if (fs!="of") {
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -     %5.3f     -       -       -       -       -  \n",Form("CMS_hww%s_%ij_Z_7TeV",fs.Data(),njets),"lnN",dyDDK);
    }
    if (mode=="shape") {
      out << Form("%-40s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n", Form("CMS_hww_%ij_WW_7TeV_SHAPE",njets),"lnU",2.0);
    } else {
      //for high mass WW is from MC so it correlated!
      out << Form("%-40s %5s     -       -       -       -     %5.3f   %5.3f     -       -       -       -       -       -       -       -  \n", (mass<=200||mode=="shape") ? Form("CMS_hww_%ij_WW_7TeV",njets) : "CMS_hww_WW","lnN",
		  1.+wwSF.second/wwSF.first,1.+wwSF.second/wwSF.first);
    }
    out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -     %5.3f     -       -  \n","CMS_hww_Wg3l","lnN",1.+WGstarScaleFactorSyst());
    out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -     %5.3f     -  \n","CMS_hww_Ztt","lnN",ZttScaleFactorKappa());
    if (mode=="cut") {
      if (zhww.first>0.) out << Form("%-40s %5s   %5.3f     -       -       -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ZH_7TeV",fs.Data(),njets),"lnN",
				     zhww.first>0 ? 1.+zhww.second/zhww.first : 1.);
      if (whww.first>0.) out << Form("%-40s %5s     -     %5.3f     -       -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_WH_7TeV",fs.Data(),njets),"lnN",
				     whww.first>0 ? 1.+whww.second/whww.first : 1.);
      out << Form("%-40s %5s     -       -     %5.3f     -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_qqH_7TeV",fs.Data(),njets),"lnN",
		  qqhww.first>0 ? 1.+qqhww.second/qqhww.first : 1.);
      out << Form("%-40s %5s     -       -       -     %5.3f     -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ggH_7TeV",fs.Data(),njets),"lnN",
		  gghww.first>0 ? 1.+gghww.second/gghww.first : 1.);
      out << Form("%-40s %5s     -       -       -       -     %5.3f     -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_WW_7TeV",fs.Data(),njets),"lnN",
		  1.+qqww.second/qqww.first);
      out << Form("%-40s %5s     -       -       -       -       -     %5.3f     -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ggWW_7TeV",fs.Data(),njets),"lnN",
		  ggww.first>0 ? 1.+ggww.second/ggww.first : 1.);
      out << Form("%-40s %5s     -       -       -       -       -       -     %5.3f     -       -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_VV_7TeV",fs.Data(),njets),"lnN",
		  1.+sqrt(pow(zz.second,2)+pow(wz.second,2)+pow(www.second,2))/(zz.first+wz.first+www.first));
      out << Form("%-40s %5s     -       -       -       -       -       -       -     %5.3f     -       -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_ttbar_7TeV",fs.Data(),njets),"lnN",
		  1.+sqrt(pow(ttbar.second,2)+pow(tw.second,2))/(ttbar.first+tw.first));
      if (fs!="of") {
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -     %5.3f     -       -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_Z_7TeV",fs.Data(),njets),"lnN",
		  dyY>0 ? 1.+dyMCE/dyY : 1.0);
      }
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -     %5.3f     -       -       -       -  \n",Form("CMS_hww%s_stat_%ij_WjetsE_7TeV",fs.Data(),njets),"lnN",
		  wjets_e.first>0 ? 1.+wjets_e.second/wjets_e.first : 1.);
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -     %5.3f     -       -       -  \n",Form("CMS_hww%s_stat_%ij_Wgamma_7TeV",fs.Data(),njets),"lnN",
		  (wgamma.first+zgamma.first)>0 ? 1.+sqrt(pow(wgamma.second,2)+pow(zgamma.second,2))/(wgamma.first+zgamma.first) : 1.);
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -     %5.3f     -       -  \n",Form("CMS_hww%s_stat_%ij_Wg3l_7TeV",fs.Data(),njets),"lnN",
		  wg3l.first>0 ? 1.+wg3l.second/wg3l.first : 1.);
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -     %5.3f     -  \n",Form("CMS_hww%s_stat_%ij_Ztt_7TeV",fs.Data(),njets),"lnN",
		  dytt.first>0 ? 1.+dytt.second/dytt.first : 1.);
      out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -       -     %5.3f\n",Form("CMS_hww%s_stat_%ij_WjetsM_7TeV",fs.Data(),njets),"lnN",
		  wjets_m.first>0 ? 1.+wjets_m.second/wjets_m.first : 1.);
    } else if (mode=="shape") {
      if (fs!="of") out << Form("%-40s %5s     -       -       -       -       -       -       -       -     2.000     -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAZBounding_7TeV",fs.Data(),njets),"shape");      
      //these are from MC so are in common for 7TeV and 8TeV
      out << Form("%-40s %5s     -       -       -       -       -       -       -     1.000     -       -       -       -       -       -  \n","CMS_hww_MVATopBounding","shape");      
      out << Form("%-40s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -       -       -  \n","CMS_hww_MVAWWBounding","shape");      
      out << Form("%-40s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -       -       -  \n","CMS_hww_MVAWWNLOBounding","shape");      
      if (zhww.first>0.)
	out << Form("%-40s %5s   1.000     -       -       -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAZHStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (whww.first>0.)
	out << Form("%-40s %5s     -     1.000     -       -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAWHStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (qqhww.first>0.)
	out << Form("%-40s %5s     -       -     1.000     -       -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAqqHStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (gghww.first>0.)
	out << Form("%-40s %5s     -       -       -     1.000     -       -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAggHStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (qqww.first>0.)
	out << Form("%-40s %5s     -       -       -       -     1.000     -       -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAqqWWStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (ggww.first>0.)
	out << Form("%-40s %5s     -       -       -       -       -     1.000     -       -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAggWWStatBounding_7TeV",fs.Data(),njets),"shape");      
      if ((zz.first+wz.first+www.first)>0.)
	out << Form("%-40s %5s     -       -       -       -       -       -     1.000     -       -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAVVStatBounding_7TeV",fs.Data(),njets),"shape");      
      if ((ttbar.first+tw.first)>0.)
	out << Form("%-40s %5s     -       -       -       -       -       -       -     1.000     -       -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVATopStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (dyY>0.)
	out << Form("%-40s %5s     -       -       -       -       -       -       -       -     1.000     -       -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAZjetsStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (wjets_e.first>0.)
	out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -     1.000     -       -       -       -  \n",Form("CMS_hww%s_%ij_MVAWjetsEStatBounding_7TeV",fs.Data(),njets),"shape");      
      if ((wgamma.first+zgamma.first)>0) 
	out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -     1.000     -       -       -  \n",Form("CMS_hww%s_%ij_MVAWgammaStatBounding_7TeV",fs.Data(),njets),"shape");      
      if ((wg3l.first)>0) 
	out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -     1.000     -       -  \n",Form("CMS_hww%s_%ij_MVAWg3lStatBounding_7TeV",fs.Data(),njets),"shape");      
      if (dytt.first>0) 
	out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -     1.000     -  \n",Form("CMS_hww%s_%ij_MVAZttStatBounding_7TeV",fs.Data(),njets),"shape");
      if (wjets_m.first>0.)
	out << Form("%-40s %5s     -       -       -       -       -       -       -       -       -       -       -       -       -     1.000\n",Form("CMS_hww%s_%ij_MVAWjetsMStatBounding_7TeV",fs.Data(),njets),"shape");      
    }
  }

  if (saveFile) {
    myfile.close();
  }

  doMVA = 0;
  doVBF = 0;

}

void cardMaker(float lumi, TString mode) {

  //int jets[] = {0};
  int jets[] = {0,1};

  //int masses[] = {125};
  //int masses[] = {120,200,350};
  //int masses[] = {110,125,160,250,300,600};
  int masses[] = {110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  int nmasses = sizeof(masses)/sizeof(int);
  int njets = sizeof(jets)/sizeof(int);
  for (int j=0;j<nmasses;++j) {
    int mass = masses[j];
    for (int jj=0;jj<njets;++jj) {
      int jetbin = jets[jj];
      if (mode=="shape" && jetbin==2) continue;
      if (jetbin!=2) {
	if (mode=="cut") cardMaker(lumi,mass,jetbin,"sffs",mode,true);
	cardMaker(lumi,mass,jetbin,"offs",mode,true);
      } else cardMaker(lumi,mass,2,"",mode,true);
    }
    cout << "mH=" << mass << " completed" << endl;
  }

}


