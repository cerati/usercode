#include "common.C"

//please set useMit = false for FR

void makeFakeTable(float lumi) {

  bool printAll = 0;

  bool useJson    = true;
  bool applyTnPSF = false;

  bool doSpillage = 1;

  int mass = 0;
  TString region = "dphireg,dphijet,minmetvtx,lep2pt15,ptll45";

  int jetbins[] = {0,1};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    pair<float, float> wjMC   = getYield(dir_mc+"wjets",  wwSelection, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, false);
    
    //electron+fake electron
    pair<float, float> wjData_ee1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", 0, useJson, false, true, false);
    pair<float, float> wjData_ee2 = getYield(data_file, wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", 0, useJson, false, true, false);
    //pair<float, float> wjData_ee1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection|Lep1LooseEleV4, Lep1FullSelection, mass, njets, region+",eefs", 0, useJson, false, true, false);
    //pair<float, float> wjData_ee2 = getYield(data_file, wwSelectionNoLep|Lep1FullSelection|Lep2LooseEleV4, Lep2FullSelection, mass, njets, region+",eefs", 0, useJson, false, true, false);
    pair<float, float> wwFake_ee1 = make_pair<float, float>(0,0);
    pair<float, float> wwFake_ee2 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_ee1 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_ee2 = make_pair<float, float>(0,0);
    if (doSpillage) {
      wwFake_ee1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
      wwFake_ee2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
      ttFake_ee1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
      ttFake_ee2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
    }
    float ee = wjData_ee1.first+wjData_ee2.first-wwFake_ee1.first-wwFake_ee2.first-ttFake_ee1.first-ttFake_ee2.first;
    float ee_err = sqrt(pow(wjData_ee1.second,2)+pow(wjData_ee2.second,2)+pow(wwFake_ee1.second,2)+pow(wwFake_ee2.second,2)+pow(ttFake_ee1.second,2)+pow(ttFake_ee2.second,2));
    
    //muon+fake electron
    pair<float, float> wjData_me1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",emfs", 0, useJson, false, true, false);
    pair<float, float> wjData_me2 = getYield(data_file,  wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", 0, useJson, false, true, false);
    pair<float, float> wwFake_me1 = make_pair<float, float>(0,0);
    pair<float, float> wwFake_me2 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_me1 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_me2 = make_pair<float, float>(0,0);
    if (doSpillage) {
      wwFake_me1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",emfs", lumi, useJson, applyTnPSF, true);
      wwFake_me2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
      ttFake_me1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
      ttFake_me2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
    }
    float me = wjData_me1.first+wjData_me2.first-wwFake_me1.first-wwFake_me2.first-ttFake_me1.first-ttFake_me2.first;
    float me_err = sqrt(pow(wjData_me1.second,2)+pow(wjData_me2.second,2)+pow(wwFake_me1.second,2)+pow(wwFake_me2.second,2)+pow(ttFake_me1.second,2)+pow(ttFake_me2.second,2));
    
    //electron+fake muon
    pair<float, float> wjData_em1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", 0, useJson, false, true, false);
    pair<float, float> wjData_em2 = getYield(data_file,  wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", 0, useJson, false, true, false);
    pair<float, float> wwFake_em1 = make_pair<float, float>(0,0);
    pair<float, float> wwFake_em2 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_em1 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_em2 = make_pair<float, float>(0,0);
    if (doSpillage) {
      wwFake_em1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
      wwFake_em2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", lumi, useJson, applyTnPSF, true);
      ttFake_em1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
      ttFake_em2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", lumi, useJson, applyTnPSF, true);
    }
    float em = wjData_em1.first+wjData_em2.first-wwFake_em1.first-wwFake_em2.first-ttFake_em1.first-ttFake_em2.first;
    float em_err = sqrt(pow(wjData_em1.second,2)+pow(wjData_em2.second,2)+pow(wwFake_em1.second,2)+pow(wwFake_em2.second,2)+pow(ttFake_em1.second,2)+pow(ttFake_em2.second,2));
    
    //muon+fake muon
    pair<float, float> wjData_mm1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", 0, useJson, false, true, false);
    pair<float, float> wjData_mm2 = getYield(data_file, wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", 0, useJson, false, true, false);
    pair<float, float> wwFake_mm1 = make_pair<float, float>(0,0);
    pair<float, float> wwFake_mm2 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_mm1 = make_pair<float, float>(0,0);
    pair<float, float> ttFake_mm2 = make_pair<float, float>(0,0);
    if (doSpillage) {
      wwFake_mm1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
      wwFake_mm2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
      ttFake_mm1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
      ttFake_mm2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
    }
    float mm = wjData_mm1.first+wjData_mm2.first-wwFake_mm1.first-wwFake_mm2.first-ttFake_mm1.first-ttFake_mm2.first;
    float mm_err = sqrt(pow(wjData_mm1.second,2)+pow(wjData_mm2.second,2)+pow(wwFake_mm1.second,2)+pow(wwFake_mm2.second,2)+pow(ttFake_mm1.second,2)+pow(ttFake_mm2.second,2));
    
    if (printAll) {
      cout << "MC: " << wjMC.first << "+/-" << wjMC.second << endl;
      
      cout << "bare data e + fake e: " << wjData_ee1.first+wjData_ee2.first << "+/-" << sqrt(pow(wjData_ee1.second,2)+pow(wjData_ee2.second,2)) << endl;
      cout << "ww spill  e + fake e: " << wwFake_ee1.first+wwFake_ee2.first << "+/-" << sqrt(pow(wwFake_ee1.second,2)+pow(wwFake_ee2.second,2)) << endl;
      cout << "tt spill  e + fake e: " << ttFake_ee1.first+ttFake_ee2.first << "+/-" << sqrt(pow(ttFake_ee1.second,2)+pow(ttFake_ee2.second,2)) << endl;
      cout << "tot data  e + fake e: " << ee << "+/-" << ee_err << endl;
      
      cout << "bare data m + fake e: " << wjData_me1.first+wjData_me2.first << "+/-" << sqrt(pow(wjData_me1.second,2)+pow(wjData_me2.second,2)) << endl;
      cout << "ww spill  m + fake e: " << wwFake_me1.first+wwFake_me2.first << "+/-" << sqrt(pow(wwFake_me1.second,2)+pow(wwFake_me2.second,2)) << endl;
      cout << "tt spill  m + fake e: " << ttFake_me1.first+ttFake_me2.first << "+/-" << sqrt(pow(ttFake_me1.second,2)+pow(ttFake_me2.second,2)) << endl;
      cout << "tot data  m + fake e: " << me << "+/-" << me_err << endl;
      
      cout << "bare data e + fake m: " << wjData_em1.first+wjData_em2.first << "+/-" << sqrt(pow(wjData_em1.second,2)+pow(wjData_em2.second,2)) << endl;
      cout << "ww spill  e + fake m: " << wwFake_em1.first+wwFake_em2.first << "+/-" << sqrt(pow(wwFake_em1.second,2)+pow(wwFake_em2.second,2)) << endl;
      cout << "tt spill  e + fake m: " << ttFake_em1.first+ttFake_em2.first << "+/-" << sqrt(pow(ttFake_em1.second,2)+pow(ttFake_em2.second,2)) << endl;
      cout << "tot data  e + fake m: " << em << "+/-" << em_err << endl;
      
      cout << "bare data m + fake m: " << wjData_mm1.first+wjData_mm2.first << "+/-" << sqrt(pow(wjData_mm1.second,2)+pow(wjData_mm2.second,2)) << endl;
      cout << "ww spill  m + fake m: " << wwFake_mm1.first+wwFake_mm2.first << "+/-" << sqrt(pow(wwFake_mm1.second,2)+pow(wwFake_mm2.second,2)) << endl;
      cout << "tt spill  m + fake m: " << ttFake_mm1.first+ttFake_mm2.first << "+/-" << sqrt(pow(ttFake_mm1.second,2)+pow(ttFake_mm2.second,2)) << endl;
      cout << "tot data  m + fake m: " << mm << "+/-" << mm_err << endl;
      
      cout << "tot Data: " << ee+me+em+mm << "+/-" << sqrt(pow(ee_err,2)+pow(me_err,2)+pow(em_err,2)+pow(mm_err,2)) << endl;
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
  }

}

void fakeBg(float lumi) {
  makeFakeTable(lumi);
}
