#include "common.C"

//please set useMit = false for FR

pair<float, float> getSpillage(TString dir, unsigned int cut_nolep, unsigned int veto, int mass, int njets, TString region, float lumi, bool applyEff, bool doPUw) {

  bool debug = 0;
  region = region+="=spill=";
  //get scale factors for DY
  float dySF = DYBkgScaleFactor(0,njets);
  dySF=1.;

  //correct for spillage...
  pair<float, float> qqwwFake  = getYield(dir+"qqww",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> ggwwFake  = getYield(dir+"ggww",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> ttbarFake = getYield(dir+"ttbar_powheg", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> twFake    = getYield(dir+"tw",    cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
//   pair<float, float> dymmFake  = getYield(dir+"dymm",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
//   pair<float, float> dyeeFake  = getYield(dir+"dyee",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> dyllFake  = getYield(dir+"dyll",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> wzFake    = getYield(dir+"wz",    cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> zzFake    = getYield(dir+"zz",    cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> wgFake    = getYield(dir+"wgamma",cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> wg3lFake  ;//= getYield(dir+"wglll",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  float spillYield = qqwwFake.first+ggwwFake.first+ttbarFake.first+twFake.first+wzFake.first+zzFake.first+
                     //dySF*dymmFake.first+dySF*dyeeFake.first+
                     dySF*dyllFake.first+
                     wgFake.first+wg3lFake.first;
  float spillError = sqrt(pow(qqwwFake.second,2)+pow(ggwwFake.second,2)+pow(ttbarFake.second,2)+pow(twFake.second,2)+
			  pow(wzFake.second,2)+pow(zzFake.second,2)+
			  //pow(dySF*dymmFake.second,2)+pow(dySF*dyeeFake.second,2)+
			  pow(dySF*dyllFake.second,2)+
			  pow(wgFake.second,2)+pow(wg3lFake.second,2));
  if (debug) {
    cout << "qqww: " << qqwwFake.first << " ggww: " << ggwwFake.first << " ttbar: " << ttbarFake.first << " tw: " << twFake.first 
         //<< " dymm: " << dySF*dymmFake.first << " dyee: " << dySF*dyeeFake.first 
	 << " dyll: " << dySF*dyllFake.first
	 << " wz: " << wzFake.first << " zz: " << zzFake.first 
	 << " wgamma: " << wgFake.first << " wg3l: " << wg3lFake.first << endl;
  }
  return make_pair<float, float>(spillYield,spillError);
}

pair<float, float> fakeBgEstimation(TString dir, unsigned int cut, unsigned int veto, int mass, unsigned int njets=0, TString region="dphijet,minmet40", 
				    float lumi=0., bool useJson=false, bool applyEff=true, bool doPUw=false)  {

  bool debug = 0;

  unsigned int cut_nolep = cut&~Lep1FullSelection;
  cut_nolep = cut_nolep&~Lep2FullSelection;

  pair<float, float> dataFake  = getYield(dir+"data.root", cut_nolep, veto, mass, njets, region, 0, useJson, false, true, false);
  //correct for spillage...
  pair<float, float> spill=getSpillage(dir,  cut_nolep, veto, mass, njets, region, lumi, applyEff, doPUw);
  float fakeYield = max(float(0.),dataFake.first-spill.first);
  float fakeError = sqrt(pow(dataFake.second,2)+pow(spill.second,2));
  if (debug) {
    cout << "total estimate, bare data, spill: " << fakeYield  << "+/-" << fakeError 
	 << " " << dataFake.first << "+/-" << dataFake.second << " " << spill.first << "+/-" << spill.second<< endl;
  }
  return make_pair<float, float>(fakeYield,fakeError);
}

pair<float, float> fakeBgEstimationWithSyst(TString dir, unsigned int cut, unsigned int veto, int mass, unsigned int njets=0, TString region="dphijet,minmet40", 
					    float lumi=0., bool useJson=false, bool applyEff=true, bool doPUw=false)  {

  pair<float, float> fake=fakeBgEstimation(dir, cut,veto,mass,njets,region,lumi,useJson,applyEff,doPUw);
  //assume 35% syst uncertainty
  return make_pair<float, float>(fake.first,sqrt(pow(fake.first*0.35,2)+pow(fake.second,2)));
}

void makeFakeTable(float lumi) {

  //bool printAll = 0;

  bool useJson    = true;
  bool applyTnPSF = true;
  bool doPUw = true;

  //bool doSpillage = 1;

  int mass = 0;
  TString region = "=dphireg=dphijet=dymvacut=ptll45=lep2pt20allfs=";

  int jetbins[] = {0};
  //int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    pair<float, float> wjMC   = getYield(main_dir+topww_dir+"/wjets",  wwSelNoMet, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_mm   = getYield(main_dir+wj_dir+"/wjets",  wwSelNoMet, noVeto, mass, njets, region+"mmfs", lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_em   = getYield(main_dir+wj_dir+"/wjets",  wwSelNoMet, noVeto, mass, njets, region+"emfs", lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_me   = getYield(main_dir+wj_dir+"/wjets",  wwSelNoMet, noVeto, mass, njets, region+"mefs", lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_ee   = getYield(main_dir+wj_dir+"/wjets",  wwSelNoMet, noVeto, mass, njets, region+"eefs", lumi, useJson, applyTnPSF, false, doPUw);

    bool doLatex = false;

    pair<float, float> wjdatamm = fakeBgEstimation(main_dir+wj_dir,wwSelNoMet, noVeto, mass, njets, region+"=mmfs=", lumi, useJson, applyTnPSF,doPUw);
    pair<float, float> wjdatame = fakeBgEstimation(main_dir+wj_dir,wwSelNoMet, noVeto, mass, njets, region+"=mefs=", lumi, useJson, applyTnPSF,doPUw);
    pair<float, float> wjdataem = fakeBgEstimation(main_dir+wj_dir,wwSelNoMet, noVeto, mass, njets, region+"=emfs=", lumi, useJson, applyTnPSF,doPUw);
    pair<float, float> wjdataee = fakeBgEstimation(main_dir+wj_dir,wwSelNoMet, noVeto, mass, njets, region+"=eefs=", lumi, useJson, applyTnPSF,doPUw);
    if (!doLatex) {
      cout << njets << "-jet bin" << endl;
      cout << Form("mm: %5.1f +/- %5.1f - me: %5.1f +/- %5.1f - em: %5.1f +/- %5.1f - ee: %5.1f +/- %5.1f",
		   wjdatamm.first,wjdatamm.second,
		   wjdatame.first,wjdatame.second,
		   wjdataem.first,wjdataem.second,
		   wjdataee.first,wjdataee.second) << endl;
      cout << Form("tot Data: %5.1f +/- %5.1f",wjdatamm.first+wjdatame.first+wjdataem.first+wjdataee.first,
		   sqrt(pow(wjdatamm.second,2)+pow(wjdatame.second,2)+pow(wjdataem.second,2)+pow(wjdataee.second,2))) << endl;
      cout << Form("MC: %5.1f +/- %5.1f", wjMC.first, wjMC.second) << endl;
    }else {
      cout << Form(" %i & %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f \\\\",
		   njets,
		   wjdatamm.first,wjdatamm.second,
		   wjdatame.first,wjdatame.second,
		   wjdataem.first,wjdataem.second,
		   wjdataee.first,wjdataee.second,
		   wjdatamm.first+wjdatame.first+wjdataem.first+wjdataee.first,sqrt(pow(wjdatamm.second,2)+pow(wjdatame.second,2)+pow(wjdataem.second,2)+pow(wjdataee.second,2))) 
	   << endl;
    }
  }
}

void fakeBg(float lumi) {
  makeFakeTable(lumi);
}


void makeSSTable(float lumi) {

  bool useJson    = true;
  bool applyTnPSF = true;
  bool doPUw = true;

  int mass = 0;
  TString region = "=dphireg=dphijet=dymvacut=ptll45=";

  //int jetbins[] = {0};
  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    unsigned int wwSel_noq = wwSelNoMet&~ChargeMatch;
    pair<float, float> ssData = getYield(main_dir+wj_dir+"data",  wwSel_noq, ChargeMatch, mass, njets, region, 0., useJson, false, false, false);
    cout << "SS data: " << ssData.first << " +/- " << ssData.second << endl;

    pair<float, float> ssFake = fakeBgEstimation(main_dir+wj_dir,  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, true, true);
    cout << "SS fake: " << ssFake.first << " +/- " << ssFake.second << endl;

    pair<float, float> ssWw = getYield(main_dir+wj_dir+"qqww",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    ssWw = make_pair<float, float>(WWBkgScaleFactorCutBased(115,njets)*ssWw.first,WWBkgScaleFactorCutBased(115,njets)*ssWw.second);
    cout << "SS ww: " << ssWw.first << " +/- " << ssWw.second << endl;

    pair<float, float> ssTtbar = getYield(main_dir+wj_dir+"ttbar_powheg",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    ssTtbar = make_pair<float, float>(TopBkgScaleFactor(njets)*ssTtbar.first,TopBkgScaleFactor(njets)*ssTtbar.second);
    cout << "SS ttbar: " << ssTtbar.first << " +/- " << ssTtbar.second << endl;

    pair<float, float> ssWgamma = getYield(main_dir+wj_dir+"wgamma",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS wgamma: " << ssWgamma.first << " +/- " << ssWgamma.second << endl;

    pair<float, float> ssWglll = getYield(main_dir+wj_dir+"wglll",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    ssWglll = make_pair<float, float>(WGstarScaleFactor()*ssWglll.first,WGstarScaleFactor()*ssWglll.second);
    cout << "SS wglll: " << ssWglll.first << " +/- " << ssWglll.second << endl;

    pair<float, float> ssWz = getYield(main_dir+wj_dir+"wz",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS wz: " << ssWz.first << " +/- " << ssWz.second << endl;
 
    pair<float, float> ssZz = getYield(main_dir+wj_dir+"zz",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS zz: " << ssZz.first << " +/- " << ssZz.second << endl;

    cout << "total expected: " << ssFake.first/*-spill.first*/+ssWw.first+ssTtbar.first+ssWgamma.first+ssWglll.first+ssWz.first+ssZz.first << " +/- " 
	 << sqrt(pow(ssFake.second,2)/*+pow(spill.second,2)*/+pow(ssWw.second,2)+pow(ssTtbar.second,2)+pow(ssWgamma.second,2)+pow(ssWglll.second,2)+pow(ssWz.second,2)+pow(ssZz.second,2)) << endl;


  }

}
