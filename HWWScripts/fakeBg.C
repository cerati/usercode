#include "common.C"

//please set useMit = false for FR

pair<float, float> getSpillage(TString dir, unsigned int cut_nolep, unsigned int veto, int mass, int njets, TString region, float lumi, bool applyEff, bool doPUw) {

  bool debug = 0;
  region = region+=",spill";
  //get scale factors for DY
  float dySF = DYBkgScaleFactor(0,njets);

  //correct for spillage...
  pair<float, float> qqwwFake    = getYield(dir+"qqww",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> ggwwFake    = getYield(dir+"ggww",  cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> ttbarFake = getYield(dir+"ttbar", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> twFake    = getYield(dir+"tw", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> dymmFake  = getYield(dir+"dymm", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> dyeeFake  = getYield(dir+"dyee", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> wzFake    = getYield(dir+"wz", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> zzFake    = getYield(dir+"zz_py", cut_nolep, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  float spillYield = qqwwFake.first+ggwwFake.first+ttbarFake.first+twFake.first+wzFake.first+zzFake.first+dySF*dymmFake.first+dySF*dyeeFake.first;
  float spillError = sqrt(pow(qqwwFake.second,2)+pow(ggwwFake.second,2)+pow(ttbarFake.second,2)+pow(twFake.second,2)+
			 pow(wzFake.second,2)+pow(zzFake.second,2)+pow(dySF*dymmFake.second,2)+pow(dySF*dyeeFake.second,2));
  if (debug) {
    cout << "qqww: " << qqwwFake.first << " ggww: " << ggwwFake.first << " ttbar: " << ttbarFake.first << " tw: " << twFake.first 
	 << " dymm: " << dySF*dymmFake.first << " dyee: " << dySF*dyeeFake.first << " wz: " << wzFake.first << " zz: " << zzFake.first << endl;
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
  float fakeYield = dataFake.first-spill.first;
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
  TString region = "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20";

  //int jetbins[] = {0};
  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    pair<float, float> wjMC   = getYield(main_dir+topww_dir+"/wjets",  wwSelection, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_mm   = getYield(main_dir+topww_dir+"/wjets",  wwSelection, noVeto, mass, njets, region+"mmfs", lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_em   = getYield(main_dir+topww_dir+"/wjets",  wwSelection, noVeto, mass, njets, region+"emfs", lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_me   = getYield(main_dir+topww_dir+"/wjets",  wwSelection, noVeto, mass, njets, region+"mefs", lumi, useJson, applyTnPSF, false, doPUw);
    //pair<float, float> wjMC_ee   = getYield(main_dir+topww_dir+"/wjets",  wwSelection, noVeto, mass, njets, region+"eefs", lumi, useJson, applyTnPSF, false, doPUw);

    bool doLatex = false;

    pair<float, float> wjdatamm = fakeBgEstimation(main_dir+topww_dir,wwSelection, noVeto, mass, njets, region+"mmfs", lumi, useJson, applyTnPSF,doPUw);
    pair<float, float> wjdatame = fakeBgEstimation(main_dir+topww_dir,wwSelection, noVeto, mass, njets, region+"mefs", lumi, useJson, applyTnPSF,doPUw);
    pair<float, float> wjdataem = fakeBgEstimation(main_dir+topww_dir,wwSelection, noVeto, mass, njets, region+"emfs", lumi, useJson, applyTnPSF,doPUw);
    pair<float, float> wjdataee = fakeBgEstimation(main_dir+topww_dir,wwSelection, noVeto, mass, njets, region+"eefs", lumi, useJson, applyTnPSF,doPUw);
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

  TString dir = "/smurf/data/Run2011_Spring11_SmurfV7_42X/mitf-alljets_Full2011/";

  bool doSpillage = 1;

  int mass = 0;
  TString region = "dphireg,dphijet,minmetvtx,lep2pt15,ptll45";

  int jetbins[] = {0};
  //int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    unsigned int wwSel_noq = wwSelection&~ChargeMatch;
    pair<float, float> ssData = getYield(dir+"data",  wwSel_noq, ChargeMatch, mass, njets, region, 0., useJson, false, false, false);
    cout << "SS data: " << ssData.first << " +/- " << ssData.second << endl;

    unsigned int wwSelNoLep_noq = wwSelectionNoLep&~ChargeMatch;
    pair<float, float> ssFake = getYield(dir+"data",  wwSelNoLep_noq, ChargeMatch, mass, njets, region, 0., useJson, false, true, false);
    cout << "SS fake: " << ssFake.first << " +/- " << ssFake.second << endl;

    pair<float, float> spill = make_pair<float, float>(0,0);
    if (doSpillage) {
      spill = getSpillage(dir,wwSelectionNoLep, noVeto, mass, njets, region+",mefs", lumi, applyTnPSF, doPUw);
    }
    cout << "spill: " << spill.first << " " << spill.second << endl;

    pair<float, float> ssWw = getYield(dir+"qqww",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS ww: " << ssWw.first << " +/- " << ssWw.second << endl;

    pair<float, float> ssTtbar = getYield(dir+"ttbar",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS ttbar: " << ssTtbar.first << " +/- " << ssTtbar.second << endl;

    pair<float, float> ssWgamma = getYield(dir+"wgamma",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS wgamma: " << ssWgamma.first << " +/- " << ssWgamma.second << endl;

    pair<float, float> ssWz = getYield(dir+"wz",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS wz: " << ssWz.first << " +/- " << ssWz.second << endl;

    pair<float, float> ssZz = getYield(dir+"zz_py",  wwSel_noq, ChargeMatch, mass, njets, region, lumi, useJson, applyTnPSF, false, doPUw);
    cout << "SS zz: " << ssZz.first << " +/- " << ssZz.second << endl;

    cout << "total expected: " << ssFake.first-spill.first+ssWw.first+ssTtbar.first+ssWgamma.first+ssWz.first+ssZz.first << " +/- " 
	 << sqrt(pow(ssFake.second,2)+pow(spill.second,2)+pow(ssWw.second,2)+pow(ssTtbar.second,2)+pow(ssWgamma.second,2)+pow(ssWz.second,2)+pow(ssZz.second,2)) << endl;


  }

}
