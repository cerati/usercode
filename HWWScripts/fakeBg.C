#include "common.C"

//please set useMit = false for FR

pair<float, float> getSpillage(TString dir, unsigned int cut, unsigned int veto, int mass, int njets, TString region, float lumi, bool applyEff, bool doPUw) {

  pair<float, float> qqwwFake_1  = getYield(dir+"qqww", cut|Lep2FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> qqwwFake_2  = getYield(dir+"qqww", cut|Lep1FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);

  pair<float, float> ggwwFake_1  = getYield(dir+"ggww", cut|Lep2FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> ggwwFake_2  = getYield(dir+"ggww", cut|Lep1FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);

  pair<float, float> ttbarFake_1 = getYield(dir+"ttbar",cut|Lep2FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> ttbarFake_2 = getYield(dir+"ttbar",cut|Lep1FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);

  pair<float, float> twFake_1 = getYield(dir+"tw",cut|Lep2FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> twFake_2 = getYield(dir+"tw",cut|Lep1FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);

  pair<float, float> wzFake_1 = getYield(dir+"wz",cut|Lep2FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> wzFake_2 = getYield(dir+"wz",cut|Lep1FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);

  pair<float, float> zzFake_1 = getYield(dir+"zz_py",cut|Lep2FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);
  pair<float, float> zzFake_2 = getYield(dir+"zz_py",cut|Lep1FullSelection, veto, mass, njets, region, lumi, false, applyEff, true, doPUw);

  //cout << "qqww: " << qqwwFake_1.first+qqwwFake_2.first << " ggww: " << ggwwFake_1.first+ggwwFake_2.first << " ttbar: " << ttbarFake_1.first+ttbarFake_2.first 
  //     << " tw: " << twFake_1.first+twFake_2.first << " wz: " << wzFake_1.first+wzFake_2.first << " zz: " << zzFake_1.first+zzFake_2.first << endl;

  float spillYield = qqwwFake_1.first+ggwwFake_1.first+ttbarFake_1.first+twFake_1.first+wzFake_1.first+zzFake_1.first+
    qqwwFake_2.first+ggwwFake_2.first+ttbarFake_2.first+twFake_2.first+wzFake_2.first+zzFake_2.first;
  float spillError = sqrt(pow(qqwwFake_1.second,2)+pow(ggwwFake_1.second,2)+pow(ttbarFake_1.second,2)+pow(twFake_1.second,2)+pow(wzFake_1.second,2)+pow(zzFake_1.second,2)
			  +pow(qqwwFake_2.second,2)+pow(ggwwFake_2.second,2)+pow(ttbarFake_2.second,2)+pow(twFake_2.second,2)+pow(wzFake_2.second,2)+pow(zzFake_2.second,2));

  return make_pair<float, float>(spillYield,spillError);
}

void makeFakeTable(float lumi) {

  bool printAll = 0;

  bool useJson    = true;
  bool applyTnPSF = true;
  bool doPUw = true;

  bool doSpillage = 1;

  int mass = 0;
  TString region = "dphireg,dphijet,minmetvtx,lep2pt15,ptll45";

  //int jetbins[] = {2};
  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    pair<float, float> wjMC   = getYield(main_dir+topww_dir+"wjets",  wwSelection, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, false, true);
    
    //electron+fake electron
    pair<float, float> wjData_ee1 = getYield(main_dir+topww_dir+"data.root", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", 0, useJson, false, true, false);
    pair<float, float> wjData_ee2 = getYield(main_dir+topww_dir+"data.root", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", 0, useJson, false, true, false);
    pair<float, float> spill_ee = make_pair<float, float>(0,0);
    if (doSpillage) {
      spill_ee = getSpillage(main_dir+topww_dir,wwSelectionNoLep, noVeto, mass, njets, region+",eefs", lumi, applyTnPSF, doPUw);
    }
    float ee = wjData_ee1.first+wjData_ee2.first-spill_ee.first;
    float ee_err = sqrt(pow(wjData_ee1.second,2)+pow(wjData_ee2.second,2)+pow(spill_ee.second,2));
    
    //muon+fake electron
    pair<float, float> wjData_me1 = getYield(main_dir+topww_dir+"data.root", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",emfs", 0, useJson, false, true, false);
    pair<float, float> wjData_me2 = getYield(main_dir+topww_dir+"data.root",  wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", 0, useJson, false, true, false);
    pair<float, float> spill_me = make_pair<float, float>(0,0);
    if (doSpillage) {
      spill_me = getSpillage(main_dir+topww_dir,wwSelectionNoLep, noVeto, mass, njets, region+",mefs", lumi, applyTnPSF, doPUw);
    }
    float me = wjData_me1.first+wjData_me2.first-spill_me.first;
    float me_err = sqrt(pow(wjData_me1.second,2)+pow(wjData_me2.second,2)+pow(spill_me.second,2));
    
    //electron+fake muon
    pair<float, float> wjData_em1 = getYield(main_dir+topww_dir+"data.root", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", 0, useJson, false, true, false);
    pair<float, float> wjData_em2 = getYield(main_dir+topww_dir+"data.root",  wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", 0, useJson, false, true, false);
    pair<float, float> spill_em = make_pair<float, float>(0,0);
    if (doSpillage) {
      spill_em = getSpillage(main_dir+topww_dir,wwSelectionNoLep, noVeto, mass, njets, region+",emfs", lumi, applyTnPSF, doPUw);
    }
    float em = wjData_em1.first+wjData_em2.first-spill_em.first;
    float em_err = sqrt(pow(wjData_em1.second,2)+pow(wjData_em2.second,2)+pow(spill_em.second,2));
    
    //muon+fake muon
    pair<float, float> wjData_mm1 = getYield(main_dir+topww_dir+"data.root", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", 0, useJson, false, true, false);
    pair<float, float> wjData_mm2 = getYield(main_dir+topww_dir+"data.root", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", 0, useJson, false, true, false);
    pair<float, float> spill_mm = make_pair<float, float>(0,0);
    if (doSpillage) {
      spill_mm = getSpillage(main_dir+topww_dir,wwSelectionNoLep, noVeto, mass, njets, region+",mmfs", lumi, applyTnPSF, doPUw);
    }
    float mm = wjData_mm1.first+wjData_mm2.first-spill_mm.first;
    float mm_err = sqrt(pow(wjData_mm1.second,2)+pow(wjData_mm2.second,2)+pow(spill_mm.second,2));
    
    if (printAll) {      
      cout << "bare data e + fake e: " << wjData_ee1.first+wjData_ee2.first << "+/-" << sqrt(pow(wjData_ee1.second,2)+pow(wjData_ee2.second,2)) << endl;
      cout << "spill     e + fake e: " << spill_ee.first << "+/-" << spill_ee.second << endl;
      cout << "tot data  e + fake e: " << ee << "+/-" << ee_err << endl;      
      cout << "bare data m + fake e: " << wjData_me1.first+wjData_me2.first << "+/-" << sqrt(pow(wjData_me1.second,2)+pow(wjData_me2.second,2)) << endl;
      cout << "spill     m + fake e: " << spill_me.first << "+/-" << spill_me.second << endl;
      cout << "tot data  m + fake e: " << me << "+/-" << me_err << endl;      
      cout << "bare data e + fake m: " << wjData_em1.first+wjData_em2.first << "+/-" << sqrt(pow(wjData_em1.second,2)+pow(wjData_em2.second,2)) << endl;
      cout << "spill     e + fake m: " << spill_em.first << "+/-" << spill_em.second << endl;
      cout << "tot data  e + fake m: " << em << "+/-" << em_err << endl;      
      cout << "bare data m + fake m: " << wjData_mm1.first+wjData_mm2.first << "+/-" << sqrt(pow(wjData_mm1.second,2)+pow(wjData_mm2.second,2)) << endl;
      cout << "spill     m + fake m: " << spill_mm.first << "+/-" << spill_mm.second << endl;
      cout << "tot data  m + fake m: " << mm << "+/-" << mm_err << endl;      
    }
    
    cout << "------------------- " << njets << "-jet bin -------------------------" << endl;
    cout << Form("| %-24s | %-24s |",
		 "Electron + Fake Electron",
		 "Muon + Fake Electron") 
	 << endl;
    cout << Form("|     %6.2f +/- %-6.2f    |     %6.2f +/- %-6.2f    |",
		 round(100.*ee)/100.,round(100.*ee_err)/100.,
		 round(100.*me)/100.,round(100.*me_err)/100.) 
	 << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << Form("| %-24s | %-24s |",
		 "Electron + Fake Muon",
		 "Muon + Fake Muon")
	 << endl;
    cout << Form("|     %6.2f +/- %-6.2f    |     %6.2f +/- %-6.2f    |",
		 round(100.*em)/100.,round(100.*em_err)/100.,
		 round(100.*mm)/100.,round(100.*mm_err)/100.) 
	 << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "tot Data: " << ee+me+em+mm << "+/-" << sqrt(pow(ee_err,2)+pow(me_err,2)+pow(em_err,2)+pow(mm_err,2)) << endl;
    cout << "MC: " << wjMC.first << "+/-" << wjMC.second << endl;
    cout << "-------------------------------------------------------" << endl;
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
