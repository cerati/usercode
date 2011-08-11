#include "common.C"

pair<float, float> fakeBgEstimation(unsigned int cut, unsigned int veto, int mass, unsigned int njets=0, TString region="dphijet,minmet40", 
				    float lumi=0., bool useJson=false, bool applyEff=true, bool doPUw=false)  {

  bool debug = 0;

  unsigned int cut_nolep = cut&~Lep1FullSelection;
  cut_nolep = cut_nolep&~Lep2FullSelection;

  pair<float, float> dataFake  = getYield(data_file, cut_nolep, veto, mass, njets, region, 0, useJson, false, true, false);
  //correct for spillage...
  pair<float, float> wwFake    = getYield(dir_mc+"qqww",  cut_nolep, veto, mass, njets, region+"spill", lumi, false, applyEff, true, doPUw);
  pair<float, float> ttbarFake = getYield(dir_mc+"ttbar", cut_nolep, veto, mass, njets, region+"spill", lumi, false, applyEff, true, doPUw);
  float fakeYield = dataFake.first-wwFake.first-ttbarFake.first;
  //assume 35% syst uncertainty
  float fakeError = sqrt(pow(dataFake.second,2)+pow(wwFake.second,2)+pow(ttbarFake.second,2)+pow(0.35*fakeYield,2));
  if (debug) {
    cout << "fakes data, ww, ttbar: " << dataFake.first << "+/-" << dataFake.second << " " 
	                              << wwFake.first << "+/-" << wwFake.second << " " 
	                              << ttbarFake.first << "+/-" << ttbarFake.second << endl;
    cout << "total estimate: " << fakeYield  << "+/-" << fakeError << endl;
  }
  return make_pair<float, float>(fakeYield,fakeError);
}

float expectedPurity(TString dir, unsigned int cut, unsigned int veto, int mass, int njets, TString myRegion, float lumi, bool useJson=0, 
		     bool applyEff=true, bool doFake=false, bool doPUw=false){

  //assume these scale factors for DY and W+jets
  float dySF = 1.;
  float wjSF = 1.;
  if (njets==0){
    dySF=4.19;
    wjSF=1.52;
  } else if (njets==1) {
    dySF=3.17;
    wjSF=2.71;
  }

  float sb_mc_top_tag = getYield(dir+"ttbar", cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first+
                        getYield(dir+"tw",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first/*+
                        getYield(dir_mc_mit+"stop",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first+
                        getYield(dir_mc_mit+"ttop",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first*/;
  float sb_mc_other_tag = getYield(dir+"qqww",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                          getYield(dir+"ggww",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first +
                          wjSF*getYield(dir+"wjets", cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                          dySF*getYield(dir+"dymm",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first +
                          dySF*getYield(dir+"dyee",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                          getYield(dir+"dytt",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first +
                          getYield(dir+"zz",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                          getYield(dir+"wz",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;

  //float nwj = fakeBgEstimation(cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doPUw).first;
  //cout << nwj << endl;

  float sb_tag_purity = sb_mc_top_tag/(sb_mc_top_tag+sb_mc_other_tag);
  return sb_tag_purity;
}

float evaluateBackground(TString dir, unsigned int cut, unsigned int veto, int mass, int njets, TString myRegion, float lumi, bool useJson=0, 
			 bool applyEff=true, bool doFake=false, bool doPUw=false){

  //assume these scale factors for DY 
  float dySF = 1.;
  if (njets==0){
    dySF=4.19;
  } else if (njets==1) {
    dySF=3.17;
  }else if (njets==2) {
    dySF=4.76;
  }

  float mc_other = getYield(dir+"qqww",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                   getYield(dir+"ggww",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first +
                   dySF*getYield(dir+"dymm",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first +
                   dySF*getYield(dir+"dyee",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                   getYield(dir+"dytt",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first +
                   getYield(dir+"zz",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first + 
                   getYield(dir+"wz",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;

  float nwj = fakeBgEstimation(cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doPUw).first;
  //cout << nwj << endl;

  float background = mc_other+nwj;
  return background;
}

pair<float, float> topVetoEffEstimation(int mass=160, unsigned int njets=0, float lumi=0., TString region = "dphijet,minmet40", 
					bool useJson=0, bool applyEff=true, bool doPUw=false) {
  //mass is not used
  bool debug = 0;
  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  TString region_top    = region+"btagJet1";
  TString region_toptag = region+"btagJet1";
  if (njets==1) {
//     //*********** fixme: ok we do not apply dphijet in the 2 jet bin
    if (region.Contains("dphijet"))  region.ReplaceAll("dphijet","");
    if (region.Contains("dpjallfs")) region.ReplaceAll("dpjallfs","");
//     //*********** 
    region_top    = region+"btagJet2";
    region_toptag = region+"btag1and2";
  }

  //get data yields
  float sideband_data_nj_top     = getYield(data_file, control_top,    veto, mass, nj_top, region_top,    0., useJson, false, false, false).first;
  float sideband_data_nj_top_tag = getYield(data_file, control_toptag, veto, mass, nj_top, region_toptag, 0., useJson, false, false, false).first;

//   //compute purity fraction
//   float sb_nj_top_frac     = expectedPurity(dir_mc, control_top,    veto, mass, nj_top, region_top,    lumi, false, applyEff, false, doPUw);
//   float sb_nj_top_tag_frac = expectedPurity(dir_mc, control_toptag, veto, mass, nj_top, region_toptag, lumi, false, applyEff, false, doPUw);
//   //FIXME: add uncertainty on purity
//   //correct for expected purity
//   sideband_data_nj_top*=sb_nj_top_frac;
//   sideband_data_nj_top_tag*=sb_nj_top_tag_frac;

  //subtract contamination
  float sb_nj_top_bkg     = evaluateBackground(dir_mc, control_top,    veto, mass, nj_top, region_top,    lumi, false, applyEff, false, doPUw);
  float sb_nj_top_tag_bkg = evaluateBackground(dir_mc, control_toptag, veto, mass, nj_top, region_toptag, lumi, false, applyEff, false, doPUw);
  sideband_data_nj_top-=sb_nj_top_bkg;
  sideband_data_nj_top_tag-=sb_nj_top_tag_bkg;

  //get tag and veto efficiencies
  float eff_tag_data = sideband_data_nj_top_tag/sideband_data_nj_top;
  float eff_err_tag_data = sqrt(eff_tag_data*(1.-eff_tag_data)/sideband_data_nj_top);//binomial
  //float eff_err_tag_data = eff_tag_data*sqrt(1./sideband_data_nj_top_tag + 1./sideband_data_nj_top);//poissonian
  float eff_veto_data = 0;
  float eff_err_veto_data = 0;
  if (njets==0) {
    float sideband_ttbar  = getYield(dir_mc+"ttbar", control_top, veto, mass, nj_top, region_top, lumi, false, applyEff, false, doPUw).first;
    float sideband_tw     = getYield(dir_mc+"tw",    control_top, veto, mass, nj_top, region_top, lumi, false, applyEff, false, doPUw).first;
    float fttbar = sideband_ttbar/(sideband_ttbar+sideband_tw);
    float fttbar_err = fttbar*0.17;//from CS uncertainty
    eff_veto_data = fttbar*(1 - (1-eff_tag_data)*(1-eff_tag_data)) + (1.-fttbar)*eff_tag_data;
    eff_err_veto_data = pow(eff_tag_data-pow(eff_tag_data,2),2)*pow(fttbar_err,2) + pow(fttbar-2*fttbar*eff_tag_data+1,2)*pow(eff_err_tag_data,2);
    eff_err_veto_data = sqrt(eff_err_veto_data);
  } else if (njets==1) {
    eff_veto_data = eff_tag_data;
    eff_err_veto_data = eff_err_tag_data;
  }
  if (debug) {
    cout << "nj cr: " << nj_top << endl;
    cout << "uncorr data cr, tag: " << sideband_data_nj_top+sb_nj_top_bkg << " " << sideband_data_nj_top_tag+sb_nj_top_tag_bkg << endl;
    cout << "backgound cr, tag: " << sb_nj_top_bkg << " " << sb_nj_top_tag_bkg << endl;
    cout << "corr. data cr, tag: " << sideband_data_nj_top << " " << sideband_data_nj_top_tag << endl;
    cout << "data eff tag,veto: " << eff_tag_data << "+/-" << eff_err_tag_data << " " << eff_veto_data << "+/-" << eff_err_veto_data << endl;
  }
  return make_pair(eff_veto_data, eff_err_veto_data);
}

pair<float, float> topBgEstimation(int mass=160, unsigned int njets=0, float lumi=0., TString region = "dphijet,minmet40", 
				   float eff_veto_data=0, float eff_err_veto_data=0, bool useJson=0, 
				   bool applyEff=true, bool doFake=false, bool doPUw=false)  {
  bool debug = 0;
  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  TString toptagreg = "";
  if (njets==1) toptagreg = ",btagJet1";

  float sideband_data_tag  = getYield(data_file, baseline_toptag, veto, mass, njets,  region+toptagreg, 0, useJson, 0, 0, 0).first;

  ////compute expected purity
  //float sb_tag_purity = expectedPurity(dir_mc, baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw);
  ////apply purity correction
  //sideband_data_tag*=sb_tag_purity;
  //float sideband_data_tag_err  = sqrt( pow(sb_tag_purity,2)*sideband_data_tag );//fixme: should account for purity error

  //compute background contamination
  float sb_tag_bkg = evaluateBackground(dir_mc, baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw);
  //apply purity correction
  sideband_data_tag-=sb_tag_bkg;
  float sideband_data_tag_err  = sqrt( (sideband_data_tag+sb_tag_bkg) + pow(sb_tag_bkg,2) );//fixme: should account for background error, now assume 100%

  //get the efficiency (if needed)
  if (eff_veto_data<1E-5) {
    pair<float, float> tagEff = topVetoEffEstimation(mass,njets,lumi,region,useJson,applyEff,doPUw);
    eff_veto_data = tagEff.first;
    eff_err_veto_data = tagEff.second;
  }
  //compute the number of top events
  float num_top_data = sideband_data_tag*(1.-eff_veto_data)/eff_veto_data;
  float num_top_err_data = sqrt( pow((1.-eff_veto_data)/eff_veto_data,2)*pow(sideband_data_tag_err,2) + pow(sideband_data_tag,2)*pow(eff_err_veto_data,2)/pow(eff_veto_data,4) );
  if (debug) {
    float sb_mc_top_tag = getYield(dir_mc+"ttbar", baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first+
                          getYield(dir_mc+"tw",    baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first;
    float sb_mc_other_tag = getYield(dir_mc+"qqww",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
                            getYield(dir_mc+"ggww",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first +
                            getYield(dir_mc+"wjets", baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
                            getYield(dir_mc+"dymm",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first +
                            getYield(dir_mc+"dyee",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
                            getYield(dir_mc+"dytt",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first +
                            getYield(dir_mc+"zz",    baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
                            getYield(dir_mc+"wz",    baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first;
    cout << "sb tag mc top, mc other: " << sb_mc_top_tag << " " << sb_mc_other_tag << endl;
    cout << "sb tag data unc, bkg, data corr: " << sideband_data_tag+sb_tag_bkg << " " << sb_tag_bkg << " " << sideband_data_tag << endl;
    //cout << "sb_tag_purity: " << sb_tag_purity << endl;
    //cout << "sb tag data unc, exp bkg, data corr: " << sideband_data_tag/sb_tag_purity << " " << (sideband_data_tag/sb_tag_purity)*(1.-sb_tag_purity) << " " << sideband_data_tag << endl;
    float sigreg_ttbar  = getYield(dir_mc+"ttbar", wwSelection, veto, mass, njets,  region, lumi).first;
    float sigreg_tw     = getYield(dir_mc+"tw",    wwSelection, veto, mass, njets,  region, lumi).first;
    cout << "data(notag)=" << num_top_data << " mc(notag)=" << sigreg_ttbar+sigreg_tw << " mc(tag)=" << sb_mc_top_tag << endl;
  }
  return make_pair<float, float>(num_top_data, num_top_err_data);
}

void makeTopTable(float lumi) {

  bool useJson  = false;
  bool applyEff = true;
  bool doFake   = false;
  bool doPUw    = true;

  TString anaRegion = "dphijet,minmet40";

  int mass = 0;

  pair<float, float> sigreg_ttbar_0j = getYield(dir_mc+"ttbar",    wwSelection, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_0j    = getYield(dir_mc+"tw",       wwSelection, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
//   pair<float, float> sigreg_stop_0j  = getYield(dir_mc_mit+"stop", wwSelection, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
//   pair<float, float> sigreg_ttop_0j  = getYield(dir_mc_mit+"ttop", wwSelection, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  //   ********* this does not work because of soft muons:
  //   pair<float, float> sigreg_ttbar_0j  = getYield(dir_mc+"ttbar", wwSelection, noVeto, mass, 0, "dphijet,minmet40", lumi, useJson, applyEff, doFake);
  //   pair<float, float> sigreg_tw_0j     = getYield(dir_mc+"tw",    wwSelection, noVeto, mass, 0, "dphijet,minmet40", lumi, useJson, applyEff, doFake);
  //   ********* so we need to do this to compare with the data prediction:
  pair<float, float> sigreg_ttbar_1j = getYield(dir_mc+"ttbar",    wwSelectionNoTV, TopTagNotInJets, mass, 1, anaRegion+"nobJet1", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_1j    = getYield(dir_mc+"tw",       wwSelectionNoTV, TopTagNotInJets, mass, 1, anaRegion+"nobJet1", lumi, useJson, applyEff, doFake, doPUw);
//   pair<float, float> sigreg_stop_1j  = getYield(dir_mc_mit+"stop", wwSelectionNoTV, TopTagNotInJets, mass, 1, anaRegion+"nobJet1", lumi, useJson, applyEff, doFake, doPUw);
//   pair<float, float> sigreg_ttop_1j  = getYield(dir_mc_mit+"ttop", wwSelectionNoTV, TopTagNotInJets, mass, 1, anaRegion+"nobJet1", lumi, useJson, applyEff, doFake, doPUw);

  pair<float, float> sigreg_top_0j = make_pair<float, float>(sigreg_ttbar_0j.first+sigreg_tw_0j.first/*+sigreg_stop_0j.first+sigreg_ttop_0j.first*/,
							     sqrt(pow(sigreg_ttbar_0j.second,2)+pow(sigreg_tw_0j.second,2)/*+pow(sigreg_stop_0j.second,2)+pow(sigreg_ttop_0j.second,2)*/));
  pair<float, float> sigreg_top_1j = make_pair<float, float>(sigreg_ttbar_1j.first+sigreg_tw_1j.first/*+sigreg_stop_1j.first+sigreg_ttop_1j.first*/,
							     sqrt(pow(sigreg_ttbar_1j.second,2)+pow(sigreg_tw_1j.second,2)/*+pow(sigreg_stop_1j.second,2)+pow(sigreg_ttop_1j.second,2)*/));

  pair<float, float> vetoEff0j = topVetoEffEstimation(mass, 0, lumi, anaRegion, useJson, applyEff, doPUw);
  pair<float, float> vetoEff1j = topVetoEffEstimation(mass, 1, lumi, anaRegion, useJson, applyEff, doPUw);

  pair<float, float> top_tag_data_0j = getYield(data_file, wwSelectionNoTV|TopTagNotInJets, noVeto,  mass, 0, anaRegion, 0, useJson, 0, 0, 0);
  pair<float, float> top_tag_data_1j = getYield(data_file, wwSelectionNoTV|OneBJet, TopTagNotInJets, mass, 1, anaRegion, 0, useJson, 0, 0, 0);

//   float cr_purity_0j = expectedPurity(dir_mc, wwSelectionNoTV|TopTagNotInJets, noVeto,  mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
//   float cr_purity_1j = expectedPurity(dir_mc, wwSelectionNoTV|OneBJet, TopTagNotInJets, mass, 1, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
//   pair<float, float> bkg_exp_cr_0j = make_pair<float, float>(top_tag_data_0j.first*(1.-cr_purity_0j),top_tag_data_0j.second*(1.-cr_purity_0j));
//   pair<float, float> bkg_exp_cr_1j = make_pair<float, float>(top_tag_data_1j.first*(1.-cr_purity_1j),top_tag_data_1j.second*(1.-cr_purity_1j));

  float cr_bkg_0j = evaluateBackground(dir_mc, wwSelectionNoTV|TopTagNotInJets, noVeto,  mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  float cr_bkg_1j = evaluateBackground(dir_mc, wwSelectionNoTV|OneBJet, TopTagNotInJets, mass, 1, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> bkg_exp_cr_0j = make_pair<float, float>(cr_bkg_0j,cr_bkg_0j);
  pair<float, float> bkg_exp_cr_1j = make_pair<float, float>(cr_bkg_1j,cr_bkg_1j);

  //this can be replaced with a simple count
  pair<float, float> topData0j = topBgEstimation(mass, 0, lumi, anaRegion, vetoEff0j.first, vetoEff0j.second, useJson, applyEff, doFake, doPUw);
  pair<float, float> topData1j = topBgEstimation(mass, 1, lumi, anaRegion, vetoEff1j.first, vetoEff1j.second, useJson, applyEff, doFake, doPUw);

  cout << "--------------------------------------------------------------------------------" << endl;
  cout << Form("| %40s | %-15s | %-15s |","Sample","0-jet","1-jet") << endl;
  cout << "--------------------------------------------------------------------------------" << endl;
  cout << Form("| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |",
	       "Estimated top events in simulation",
	       round(10.*sigreg_top_0j.first)/10.,round(10.*sigreg_top_0j.second)/10.,
	       round(10.*sigreg_top_1j.first)/10.,round(10.*sigreg_top_1j.second)/10.) 
       << endl;
  cout << Form("| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |",
	       "tagging efficiency (%)",
	       round(1000.*vetoEff0j.first)/10.,round(1000.*vetoEff0j.second)/10.,
	       round(1000.*vetoEff1j.first)/10.,round(1000.*vetoEff1j.second)/10.) 
       << endl;
  cout << Form("| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |",
	       "top-tagged events in data",
	       round(10.*top_tag_data_0j.first)/10.,round(10.*top_tag_data_0j.second)/10.,
	       round(10.*top_tag_data_1j.first)/10.,round(10.*top_tag_data_1j.second)/10.) 
       << endl;
  cout << Form("| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |",
	       "background events in control region",
	       round(10.*bkg_exp_cr_0j.first)/10.,round(10.*bkg_exp_cr_0j.second)/10.,
	       round(10.*bkg_exp_cr_1j.first)/10.,round(10.*bkg_exp_cr_1j.second)/10.) 
       << endl;
  cout << Form("| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |",
	       "Data-driven top background estimate",
	       round(10.*topData0j.first)/10.,round(10.*topData0j.second)/10.,
	       round(10.*topData1j.first)/10.,round(10.*topData1j.second)/10.) 
       << endl;
  cout << "--------------------------------------------------------------------------------" << endl;

}
