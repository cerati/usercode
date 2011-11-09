#include "common.C"
void cardMaker(float lumi, int mass, unsigned int njets) {

  TString dir = "/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/";

  bool useJson = 1;
  bool applyEff=true;
  bool doFake=false; 
  bool doPUw=true;

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  pair<float, float> wwSF = make_pair<float, float>(1.1,0.11);
  pair<float, float> topSF = make_pair<float, float>(1.0,0.01);
  pair<float, float> dySF = make_pair<float, float>(1.0,0.01);

  pair<float, float> data = getYield(dir+"data", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", 0.,   useJson, false, false, false);

  pair<float, float> qqww = getYield(dir+"qqww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> ggww = getYield(dir+"ggww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> ttbar = getYield(dir+"ttbar", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> tw = getYield(dir+"tw", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> zz = getYield(dir+"zz", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wz = getYield(dir+"wz", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> dyee = getYield(dir+"dyee", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> dymm = getYield(dir+"dymm", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> dytt = getYield(dir+"dytt", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> wgamma = getYield(dir+"wgamma", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> wjets = getYield(dir+"data", wwSelectionNoLep, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs", 0.,   useJson, false, true, false);

  pair<float, float> gghww = getYield(dir+Form("hww%i",mass), wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs,ggH", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> qqhww = getYield(dir+Form("hww%i",mass), wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,sffs,qqH", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zhww = make_pair<float, float>(0,0);
  pair<float, float> whww = make_pair<float, float>(0,0);

  printf("%s \n","imax 1 number of channels");
  printf("%s \n","jmax * number of background");
  printf("%s \n","kmax * number of nuisance parameters");
  printf("%s %i\n","Observation",((int)data.first));
  printf("%s \n","bin 1 1 1 1 1 1 1 1 1 1 1 1");
  printf("%s \n","process ZH WH qqH ggH qqWW ggWW VV Top Zjets Wjets Wgamma Ztt");
  printf("%s \n","process -3 -2 -1 0 1 2 3 4 5 6 7 8");
  printf("%s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n","rate", 
	 zhww.first, whww.first, qqhww.first, gghww.first, wwSF.first*qqww.first, wwSF.first*ggww.first, zz.first+wz.first, topSF.first*(ttbar.first+tw.first), 
	 dySF.first*(dymm.first+dyee.first), wjets.first, wgamma.first, dytt.first);
  printf("%s \n","lumi                       lnN 1.045 1.045 1.045 1.045   -     -   1.045   -     -     -   1.045 1.045");
  printf("%s \n","CMS_eff_m                  lnN 1.030 1.030 1.030 1.030 1.030 1.030 1.030   -     -     -   1.030 1.030");
  printf("%s \n","CMS_eff_e                  lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040   -     -     -   1.040 1.040");
  printf("%s \n","CMS_scale_m                lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.015   -     -     -   1.015 1.015");
  printf("%s \n","CMS_scale_e                lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020");
  printf("%s \n","CMS_hww_met_resolution     lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020");
  printf("%s \n","CMS_scale_j                lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020");
  printf("%s \n","CMS_hww_fake_em            lnN   -     -     -     -     -     -     -     -     -   1.360   -     -  ");
  printf("%s \n","UEPS                       lnN   -     -     -   0.943   -     -     -     -     -     -     -     -  ");
  printf("%s \n","pdf_gg                     lnN   -     -     -   1.100   -   1.040   -     -     -     -     -     -  ");
  printf("%s \n","pdf_qqbar                  lnN 1.050 1.050 1.050   -   1.040   -   1.040   -     -     -   1.040 1.040");
  printf("%s \n","QCDscale_ggH               lnN   -     -     -   1.150   -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_ggH1in            lnN   -     -     -   0.900   -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_ggH2in            lnN   -     -     -   1.000   -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_qqH               lnN   -     -   1.010   -     -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_VH                lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_VV                lnN   -     -     -     -     -     -   1.040   -     -     -     -     -  ");
  printf("%s \n","QCDscale_V                 lnN   -     -     -     -     -     -     -     -     -     -     -   1.100");
  printf("%s \n","QCDscale_ggVV              lnN   -     -     -     -     -   1.500   -     -     -     -     -     -  ");
  printf("%s \n","CMS_QCDscale_WW_EXTRAP     lnN   -     -     -     -   0.954   -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_ggH_ACEPT         lnN   -     -     -   1.020   -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_qqH_ACEPT         lnN   -     -   1.020   -     -     -     -     -     -     -     -     -  ");
  printf("%s \n","QCDscale_VH_ACEPT          lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -  ");
  printf("CMS_hww_0j_ttbar           lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",1.+topSF.second/topSF.first);
  printf("CMS_hwwsf_0j_Z             lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",1.+dySF.second/dySF.first);
  printf("CMS_hww_0j_WW              lnN   -     -     -     -   %5.3f %5.3f   -     -     -     -     -     -  \n",1.+wwSF.second/wwSF.first,1.+wwSF.second/wwSF.first);
  printf("CMS_hwwsf_stat_0j_ZH       lnN %5.3f   -     -     -     -     -     -     -     -     -     -     -  \n",zhww.first>0 ? 1.+zhww.second/zhww.first : 0.);
  printf("CMS_hwwsf_stat_0j_WH       lnN   -   %5.3f   -     -     -     -     -     -     -     -     -     -  \n",whww.first>0 ? 1.+whww.second/whww.first : 0.);
  printf("CMS_hwwsf_stat_0j_qqH      lnN   -     -   %5.3f   -     -     -     -     -     -     -     -     -  \n",qqhww.first>0 ? 1.+qqhww.second/qqhww.first : 0.);
  printf("CMS_hwwsf_stat_0j_ggH      lnN   -     -     -   %5.3f   -     -     -     -     -     -     -     -  \n",gghww.first>0 ? 1.+gghww.second/gghww.first : 0.);
  printf("CMS_hwwsf_stat_0j_WW       lnN   -     -     -     -   %5.3f   -     -     -     -     -     -     -  \n",1.+qqww.second/qqww.first);
  printf("CMS_hwwsf_stat_0j_ggWW     lnN   -     -     -     -     -   %5.3f   -     -     -     -     -     -  \n",1.+ggww.second/ggww.first);
  printf("CMS_hwwsf_stat_0j_VV       lnN   -     -     -     -     -     -   %5.3f   -     -     -     -     -  \n",1.+sqrt(pow(zz.second,2)+pow(wz.second,2))/(zz.first+wz.first));
  printf("CMS_hwwsf_stat_0j_ttbar    lnN   -     -     -     -     -     -     -   %5.3f   -     -     -     -  \n",1.+sqrt(pow(ttbar.second,2)+pow(tw.second,2))/(ttbar.first+tw.first));
  printf("CMS_hwwsf_stat_0j_Z        lnN   -     -     -     -     -     -     -     -   %5.3f   -     -     -  \n",1.+sqrt(pow(dymm.second,2)+pow(dyee.second,2))/(dymm.first+dyee.first));
  printf("CMS_hwwsf_stat_0j_Wjets    lnN   -     -     -     -     -     -     -     -     -   %5.3f   -     -  \n",wjets.first>0 ? 1.+wjets.second/wjets.first : 0.);
  printf("CMS_hwwsf_stat_0j_Wgamma   lnN   -     -     -     -     -     -     -     -     -     -   %5.3f   -  \n",wgamma.first>0 ? 1.+wgamma.second/wgamma.first : 0.);
  printf("CMS_hwwsf_stat_0j_Ztt      lnN   -     -     -     -     -     -     -     -     -     -     -   %5.3f\n",dytt.first>0 ? 1.+dytt.second/dytt.first : 0.);

}
