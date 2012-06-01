#include "common.C"
#include "fakeBg.C"

void mvaSyst(float lumi) {

  bool useJson    = false;
  bool applyTnPSF = true;
  bool doPUw = true;

  TString dir = main_dir;

  int mass = 0;
  TString regionDen = "=dphiside=mm20allfs=ptll45=offs=";
  TString regionNum = "=dphiside=dymvaallfs=ptll45=offs=";

  //int jetbins[] = {1};
  int jetbins[] = {0,1};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    cout << endl;
    cout << "njets==" << njets << endl;

    unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
    getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

    float wwSF   = WWBkgScaleFactorCutBased(115,njets);
    float wwSFe  = 0.*WWBkgScaleFactorCutBased(115,njets)*(WWBkgScaleFactorKappaCutBased(115,njets)-1.);
    float topSF  = TopBkgScaleFactor(njets);
    float topSFe = 0.*TopBkgScaleFactor(njets)*(TopBkgScaleFactorKappa(njets)-1.);
    wwSF=1.;
    topSF=1.;

    pair<float, float> qqwwDen = getYield(dir+"qqww",  wwSelNoMet, veto, mass, njets, regionDen, lumi, useJson, applyTnPSF, false, doPUw);
    pair<float, float> qqwwNum = getYield(dir+"qqww",  wwSelNoMet, veto, mass, njets, regionNum, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "qqwwDen: " << wwSF*qqwwDen.first << " +/- " << sqrt(pow(wwSF*qqwwDen.second,2)+pow(wwSFe*qqwwDen.first,2))
	 << " qqwwNum: " << wwSF*qqwwNum.first << " +/- " << sqrt(pow(wwSF*qqwwNum.second,2)+pow(wwSFe*qqwwNum.first,2))
	 << " qqwwEff: " << (wwSF*qqwwNum.first)/(wwSF*qqwwDen.first) << " +/- " << efficiencyErr((wwSF*qqwwNum.first)/(wwSF*qqwwDen.first), wwSF*qqwwDen.first) << endl;

    pair<float, float> ggwwDen = getYield(dir+"ggww",  wwSelNoMet, veto, mass, njets, regionDen, lumi, useJson, applyTnPSF, false, doPUw);
    pair<float, float> ggwwNum = getYield(dir+"ggww",  wwSelNoMet, veto, mass, njets, regionNum, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "ggwwDen: " << wwSF*ggwwDen.first << " +/- " <<  sqrt(pow(wwSF*ggwwDen.second,2)+pow(wwSFe*ggwwDen.first,2))
	 << " ggwwNum: " << wwSF*ggwwNum.first << " +/- " << sqrt(pow(wwSF*ggwwNum.second,2)+pow(wwSFe*ggwwNum.first,2))
	 << " ggwwEff: " << (wwSF*ggwwNum.first)/(wwSF*ggwwDen.first) << " +/- " << efficiencyErr((wwSF*ggwwNum.first)/(wwSF*ggwwDen.first), wwSF*ggwwDen.first) << endl;

    pair<float, float> ttbarDen = getYield(dir+"ttbar",  wwSelNoMet, veto, mass, njets, regionDen, lumi, useJson, applyTnPSF, false, doPUw);
    pair<float, float> ttbarNum = getYield(dir+"ttbar",  wwSelNoMet, veto, mass, njets, regionNum, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "ttbarDen: " << topSF*ttbarDen.first << " +/- " << sqrt(pow(topSF*ttbarDen.second,2)+ pow(topSFe*ttbarDen.first,2))
	 << " ttbarNum: " << topSF*ttbarNum.first << " +/- " << sqrt(pow(topSF*ttbarNum.second,2)+ pow(topSFe*ttbarNum.first,2))
	 << " ttbarEff: " << (topSF*ttbarNum.first)/(topSF*ttbarDen.first) << " +/- " << efficiencyErr((topSF*ttbarNum.first)/(topSF*ttbarDen.first), topSF*ttbarDen.first) << endl;

    pair<float, float> twDen = getYield(dir+"tw",  wwSelNoMet, veto, mass, njets, regionDen, lumi, useJson, applyTnPSF, false, doPUw);
    pair<float, float> twNum = getYield(dir+"tw",  wwSelNoMet, veto, mass, njets, regionNum, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "twDen: " << topSF*twDen.first << " +/- " << sqrt(pow(topSF*twDen.second,2)+ pow(topSFe*twDen.first,2))
	 << " twNum: " << topSF*twNum.first << " +/- " << sqrt(pow(topSF*twNum.second,2)+ pow(topSFe*twNum.first,2))
	 << " twEff: " << (topSF*twNum.first)/(topSF*twDen.first) << " +/- " << efficiencyErr((topSF*twNum.first)/(topSF*twDen.first), topSF*twDen.first) << endl;

    pair<float, float> wzDen = getYield(dir+"wz",  wwSelNoMet, veto, mass, njets, regionDen, lumi, useJson, applyTnPSF, false, doPUw);
    pair<float, float> wzNum = getYield(dir+"wz",  wwSelNoMet, veto, mass, njets, regionNum, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "wzDen: " << wzDen.first << " +/- " << wzDen.second 
	 << " wzNum: " << wzNum.first << " +/- " << wzNum.second 
	 << " wzEff: " << wzNum.first/wzDen.first << " +/- " << efficiencyErr(wzNum.first/wzDen.first, wzDen.first) << endl;

    pair<float, float> zzDen = getYield(dir+"zz",  wwSelNoMet, veto, mass, njets, regionDen, lumi, useJson, applyTnPSF, false, doPUw);
    pair<float, float> zzNum = getYield(dir+"zz",  wwSelNoMet, veto, mass, njets, regionNum, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "zzDen: " << zzDen.first << " +/- " << zzDen.second 
	 << " zzNum: " << zzNum.first << " +/- " << zzNum.second 
	 << " zzEff: " << zzNum.first/zzDen.first << " +/- " << efficiencyErr(zzNum.first/zzDen.first, zzDen.first) << endl;


    pair<float, float> wjetsDen = fakeBgEstimation(dir, wwSelNoMet, veto, mass, njets, regionDen, lumi,   useJson, applyTnPSF, doPUw);
    pair<float, float> wjetsNum = fakeBgEstimation(dir, wwSelNoMet, veto, mass, njets, regionNum, lumi,   useJson, applyTnPSF, doPUw);
    cout << "wjetsDen: " << wjetsDen.first << " +/- " << wjetsDen.second 
	 << " wjetsNum: " << wjetsNum.first << " +/- " << wjetsNum.second 
	 << " wjetsEff: " << wjetsNum.first/wjetsDen.first << " +/- " << efficiencyErr(wjetsNum.first/wjetsDen.first, wjetsDen.first) << endl;


    ///////////////////////////////

    float mcDen = wwSF*qqwwDen.first+wwSF*ggwwDen.first+topSF*ttbarDen.first+topSF*twDen.first+wzDen.first+zzDen.first+wjetsDen.first;
    float mcNum = wwSF*qqwwNum.first+wwSF*ggwwNum.first+topSF*ttbarNum.first+topSF*twNum.first+wzNum.first+zzNum.first+wjetsNum.first;
    float mcDenErr = sqrt(pow(wwSF*qqwwDen.second,2)+pow(wwSFe*qqwwDen.first,2)+pow(wwSF*ggwwDen.second,2)+pow(wwSFe*ggwwDen.first,2)+
			  pow(topSF*ttbarDen.second,2)+pow(topSFe*ttbarDen.first,2)+pow(topSF*twDen.second,2)+pow(topSFe*twDen.first,2)+
			  pow(wzDen.second,2)+pow(zzDen.second,2)+pow(wjetsDen.second,2));
    float mcNumErr = sqrt(pow(wwSF*qqwwNum.second,2)+pow(wwSFe*qqwwNum.first,2)+pow(wwSF*ggwwNum.second,2)+pow(wwSFe*ggwwNum.first,2)+
			  pow(topSF*ttbarNum.second,2)+pow(topSFe*ttbarNum.first,2)+pow(topSF*twNum.second,2)+pow(topSFe*twNum.first,2)+
			  pow(wzNum.second,2)+pow(zzNum.second,2)+pow(wjetsNum.second,2));

    float qqwwscale = getScale1fb(dir+"qqww");

    cout << "mcDen: " << mcDen << " +/- " << mcDenErr
	 << " mcNum: " << mcNum << " +/- " << mcNumErr 
	 << " mcEff: " << mcNum/mcDen << " +/- " << efficiencyErr(mcNum/mcDen, mcDen/qqwwscale) << endl;

    pair<float, float> dataDen = getYield(dir+"data",  wwSelNoMet, veto, mass, njets, regionDen, 0., useJson, false, false, false);
    pair<float, float> dataNum = getYield(dir+"data",  wwSelNoMet, veto, mass, njets, regionNum, 0., useJson, false, false, false);
    cout << "dataDen: " << dataDen.first << " +/- " << dataDen.second 
	 << " dataNum: " << dataNum.first << " +/- " << dataNum.second 
	 << " dataEff: " << dataNum.first/dataDen.first << " +/- " << efficiencyErr(dataNum.first/dataDen.first, dataDen.first) << endl;

  }



}
