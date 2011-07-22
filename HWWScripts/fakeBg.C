#include "common.C"

void makeFakeTable(float lumi) {

  bool printAll = 0;

  bool useJson    = false;
  bool applyTnPSF = false;

  int mass = 0;
  TString region = "minmet40,dphijet,dphireg";

  int jetbins[] = {0,1};
  int njetbins = sizeof(jetbins)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];

    pair<float, float> wjMC   = getYield(dir_mc+"wjets",  wwSelection, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, false);
    
    //electron+fake electron
    pair<float, float> wjData_ee1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", 0, useJson, applyTnPSF, true);
    pair<float, float> wjData_ee2 = getYield(data_file, wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", 0, useJson, applyTnPSF, true);
    pair<float, float> wwFake_ee1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> wwFake_ee2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_ee1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_ee2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",eefs", lumi, useJson, applyTnPSF, true);
    float ee = wjData_ee1.first+wjData_ee2.first-wwFake_ee1.first-wwFake_ee2.first-ttFake_ee1.first-ttFake_ee2.first;
    float ee_err = sqrt(pow(wjData_ee1.second,2)+pow(wjData_ee2.second,2)+pow(wwFake_ee1.second,2)+pow(wwFake_ee2.second,2)+pow(ttFake_ee1.second,2)+pow(ttFake_ee2.second,2));
    
    //muon+fake electron
    pair<float, float> wjData_me1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",emfs", 0, useJson, applyTnPSF, true);
    pair<float, float> wjData_me2 = getYield(data_file,  wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", 0, useJson, applyTnPSF, true);
    pair<float, float> wwFake_me1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",emfs", lumi, useJson, applyTnPSF, true);
    pair<float, float> wwFake_me2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_me1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_me2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
    float me = wjData_me1.first+wjData_me2.first-wwFake_me1.first-wwFake_me2.first-ttFake_me1.first-ttFake_me2.first;
    float me_err = sqrt(pow(wjData_me1.second,2)+pow(wjData_me2.second,2)+pow(wwFake_me1.second,2)+pow(wwFake_me2.second,2)+pow(ttFake_me1.second,2)+pow(ttFake_me2.second,2));
    
    //electron+fake muon
    pair<float, float> wjData_em1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", 0, useJson, applyTnPSF, true);
    pair<float, float> wjData_em2 = getYield(data_file,  wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", 0, useJson, applyTnPSF, true);
    pair<float, float> wwFake_em1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> wwFake_em2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_em1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mefs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_em2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",emfs", lumi, useJson, applyTnPSF, true);
    float em = wjData_em1.first+wjData_em2.first-wwFake_em1.first-wwFake_em2.first-ttFake_em1.first-ttFake_em2.first;
    float em_err = sqrt(pow(wjData_em1.second,2)+pow(wjData_em2.second,2)+pow(wwFake_em1.second,2)+pow(wwFake_em2.second,2)+pow(ttFake_em1.second,2)+pow(ttFake_em2.second,2));
    
    //muon+fake muon
    pair<float, float> wjData_mm1 = getYield(data_file, wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", 0, useJson, applyTnPSF, true);
    pair<float, float> wjData_mm2 = getYield(data_file, wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", 0, useJson, applyTnPSF, true);
    pair<float, float> wwFake_mm1 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
    pair<float, float> wwFake_mm2 = getYield(dir_mc+"qqww", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_mm1 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep2FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
    pair<float, float> ttFake_mm2 = getYield(dir_mc+"ttbar", wwSelectionNoLep|Lep1FullSelection, noVeto, mass, njets, region+",mmfs", lumi, useJson, applyTnPSF, true);
    float mm = wjData_mm1.first+wjData_mm2.first-wwFake_mm1.first-wwFake_mm2.first-ttFake_mm1.first-ttFake_mm2.first;
    float mm_err = sqrt(pow(wjData_mm1.second,2)+pow(wjData_mm2.second,2)+pow(wwFake_mm1.second,2)+pow(wwFake_mm2.second,2)+pow(ttFake_mm1.second,2)+pow(ttFake_mm2.second,2));
    
    if (printAll) {
      cout << "MC: " << wjMC.first << "+/-" << wjMC.second << endl;
      
      cout << "bare data e+fe: " << wjData_ee1.first+wjData_ee2.first << "+/-" << sqrt(pow(wjData_ee1.second,2)+pow(wjData_ee2.second,2)) << endl;
      cout << "ww spill  e+fe: " << wwFake_ee1.first+wwFake_ee2.first << "+/-" << sqrt(pow(wwFake_ee1.second,2)+pow(wwFake_ee2.second,2)) << endl;
      cout << "tt spill  e+fe: " << ttFake_ee1.first+ttFake_ee2.first << "+/-" << sqrt(pow(ttFake_ee1.second,2)+pow(ttFake_ee2.second,2)) << endl;
      cout << "tot data  e+fe: " << ee << "+/-" << ee_err << endl;
      
      cout << "bare data m+fe: " << wjData_me1.first+wjData_me2.first << "+/-" << sqrt(pow(wjData_me1.second,2)+pow(wjData_me2.second,2)) << endl;
      cout << "ww spill  m+fe: " << wwFake_me1.first+wwFake_me2.first << "+/-" << sqrt(pow(wwFake_me1.second,2)+pow(wwFake_me2.second,2)) << endl;
      cout << "tt spill  m+fe: " << ttFake_me1.first+ttFake_me2.first << "+/-" << sqrt(pow(ttFake_me1.second,2)+pow(ttFake_me2.second,2)) << endl;
      cout << "tot data  m+fe: " << me << "+/-" << me_err << endl;
      
      cout << "bare data e+fm: " << wjData_em1.first+wjData_em2.first << "+/-" << sqrt(pow(wjData_em1.second,2)+pow(wjData_em2.second,2)) << endl;
      cout << "ww spill  e+fm: " << wwFake_em1.first+wwFake_em2.first << "+/-" << sqrt(pow(wwFake_em1.second,2)+pow(wwFake_em2.second,2)) << endl;
      cout << "tt spill  e+fm: " << ttFake_em1.first+ttFake_em2.first << "+/-" << sqrt(pow(ttFake_em1.second,2)+pow(ttFake_em2.second,2)) << endl;
      cout << "tot data  e+fm: " << em << "+/-" << em_err << endl;
      
      cout << "bare data m+fm: " << wjData_mm1.first+wjData_mm2.first << "+/-" << sqrt(pow(wjData_mm1.second,2)+pow(wjData_mm2.second,2)) << endl;
      cout << "ww spill  m+fm: " << wwFake_mm1.first+wwFake_mm2.first << "+/-" << sqrt(pow(wwFake_mm1.second,2)+pow(wwFake_mm2.second,2)) << endl;
      cout << "tt spill  m+fm: " << ttFake_mm1.first+ttFake_mm2.first << "+/-" << sqrt(pow(ttFake_mm1.second,2)+pow(ttFake_mm2.second,2)) << endl;
      cout << "tot data  m+fm: " << mm << "+/-" << mm_err << endl;
      
      cout << "tot Data: " << ee+me+em+mm << "+/-" << sqrt(pow(ee_err,2)+pow(me_err,2)+pow(em_err,2)+pow(mm_err,2)) << endl;
    }
    
    
    cout << "------------------- " << njets << "-jet bin -------------------------" << endl;
    cout << Form("| %-24s | %-24s |",
		 "Electron + Fake Electron",
		 "Muon + Fake Electron") 
	 << endl;
    cout << Form("|     %-6.2f +/- %-6.2f    |     %-6.2f +/- %-6.2f    |",
		 round(100.*ee)/100.,round(100.*ee_err)/100.,
		 round(100.*me)/100.,round(100.*me_err)/100.) 
	 << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << Form("| %-24s | %-24s |",
		 "Electron + Fake Muon",
		 "Muon + Fake Muon")
	 << endl;
    cout << Form("|     %-6.2f +/- %-6.2f    |     %-6.2f +/- %-6.2f    |",
		 round(100.*em)/100.,round(100.*em_err)/100.,
		 round(100.*mm)/100.,round(100.*mm_err)/100.) 
	 << endl;
    cout << "-------------------------------------------------------" << endl;
  }

}
