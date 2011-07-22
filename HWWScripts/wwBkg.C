#include "common.C"
#include "topBg.C"

pair<float,float> wwEstimationMC(int mass=160, unsigned int njets=0, float lumi = 1./*fb-1*/){

  bool printAll = 0;

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  float scale_qq = getScale1fb(dir_mc+"qqww")*lumi;
  float scale_gg = getScale1fb(dir_mc+"ggww")*lumi;
  float sigreg_qqww_l = getYield(dir_mc+"qqww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmet40", lumi).first;
  float sigreg_ggww_l = getYield(dir_mc+"ggww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmet40", lumi).first;
  float sr_ww_true = sigreg_qqww_l+sigreg_ggww_l;
  float sr_ww_true_err = sqrt(scale_qq*sigreg_qqww_l+scale_gg*sigreg_ggww_l);

  if (printAll) cout << "ww MC: " << sr_ww_true << "+/-" << sr_ww_true_err << endl;

  return make_pair<float,float>(sr_ww_true,sr_ww_true_err);

}

pair<float, float> fakeBgEstimation(int mass=160, unsigned int njets=0, TString region="dphijet,minmet40", float lumi=0., bool useJson=false, bool applyTnPSF=false)  {
  bool debug = false;
  pair<float, float> dataFake  = getYield(data_file, wwSelectionNoLep, noVeto, mass, njets, region, 0, useJson, applyTnPSF, true);
  //correct for spillage...
  pair<float, float> wwFake    = getYield(dir_mc+"qqww",  wwSelectionNoLep, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, true);
  pair<float, float> ttbarFake = getYield(dir_mc+"ttbar", wwSelectionNoLep, noVeto, mass, njets, region, lumi, useJson, applyTnPSF, true);
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

pair<float,float> wwEstimationData(int mass=160, unsigned int njets=0, float lumi = 1./*fb-1*/, float eff_veto_data=0, float eff_err_veto_data=0){

  bool printAll = 0;

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  float sideband_qqww = getYield(dir_mc+"qqww", wwSelection, veto, mass, njets, "sideband,dphijet,minmet40", 0.).first;
  float sideband_ggww = getYield(dir_mc+"ggww", wwSelection, veto, mass, njets, "sideband,dphijet,minmet40", 0.).first;

  float massreg_qqww = getYield(dir_mc+"qqww", wwSelection, veto, mass, njets, "massreg,dphijet,minmet40", 0.).first;
  float massreg_ggww = getYield(dir_mc+"ggww", wwSelection, veto, mass, njets, "massreg,dphijet,minmet40", 0.).first;

  float mtreg_qqww = getYield(dir_mc+"qqww", wwSelection, veto, mass, njets, "mtreg,dphijet,minmet40", 0.).first;
  float mtreg_ggww = getYield(dir_mc+"ggww", wwSelection, veto, mass, njets, "mtreg,dphijet,minmet40", 0.).first;

  float sigreg_qqww = getYield(dir_mc+"qqww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmet40", 0.).first;
  float sigreg_ggww = getYield(dir_mc+"ggww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmet40", 0.).first;

  float scale_qq = getScale1fb(dir_mc+"qqww")*lumi;
  float scale_gg = getScale1fb(dir_mc+"ggww")*lumi;
  float weight_gg = scale_gg/scale_qq;

  float rio_mg = (massreg_qqww+massreg_ggww*weight_gg)/(sideband_qqww+sideband_ggww*weight_gg);
  float rio_err_mg = rio_mg*sqrt(1./(massreg_qqww+massreg_ggww*weight_gg)+1./(sideband_qqww+sideband_ggww*weight_gg));
  float eff_mt_mr_mg = (mtreg_qqww+mtreg_ggww*weight_gg)/(massreg_qqww+massreg_ggww*weight_gg);
  float eff_err_mt_mr_mg = eff_mt_mr_mg*sqrt(1./(mtreg_qqww+mtreg_ggww*weight_gg)+1./(massreg_qqww+massreg_ggww*weight_gg));
  float eff_dp_mr_mg = (sigreg_qqww+sigreg_ggww*weight_gg)/(mtreg_qqww+mtreg_ggww*weight_gg);
  float eff_err_dp_mr_mg = eff_dp_mr_mg*sqrt(1./(sigreg_qqww+sigreg_ggww*weight_gg)+1./(mtreg_qqww+mtreg_ggww*weight_gg));

  float sideband_data = getYield(data_file, wwSelection, veto, mass, njets,  "sideband,dphijet,minmet40", 0).first;

  //********* get top **********
  pair<float, float> ttbar_tmp = getYield(dir_mc+"ttbar", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi);
  pair<float, float> tw_tmp    = getYield(dir_mc+"tw", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi);
  pair<float, float> top = topBgEstimation(mass, njets, lumi, "sideband,dphijet,minmet40", eff_veto_data, eff_err_veto_data);
  float num_top_data = top.first;
  float num_top_err_data = top.second;
  //********* get top **********

  //********* get w+jets **********
  pair<float, float> wj = fakeBgEstimation(mass, njets, "sideband,dphijet,minmet40", lumi);
  float num_wjets_data = wj.first;
  float num_wjets_err_data = wj.second;
  //********* get w+jets **********

  //********* get DY ********** => take it from MC
  //pair<float, float> dy = dyBkgEstimation(data_file, wwSelection, veto, mass, njets,  "sideband,dphijet,minmet40", 0, 0);
  //float num_dy_data = dy.first;
  //float num_dy_err_data = dy.second;
  pair<float, float> sb_dymm   = getYield(dir_mc+"dymm",  wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi);
  pair<float, float> sb_dyee   = getYield(dir_mc+"dyee",  wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi);
  float sideband_dymm          = sb_dymm.first;
  float sideband_dyee          = sb_dyee.first;
  float sideband_dymm_err      = sb_dymm.second;
  float sideband_dyee_err      = sb_dyee.second;
  float dySF = 2.28;
  float num_dy_data = (sideband_dyee+sideband_dymm)*dySF;
  float num_dy_err_data = num_dy_data;
  //********* get DY **********

  //********* get others **********
  float sideband_dytt          = getYield(dir_mc+"dytt",  wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
  float sideband_zz            = getYield(dir_mc+"zz",    wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
  float sideband_wz            = getYield(dir_mc+"wz",    wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
  float num_other_mc = sideband_zz+sideband_wz+sideband_dytt;
  float num_other_err_mc = num_other_mc;//100% for now
  //********* get others **********

  float sb_ww_meas_data = sideband_data-num_top_data-num_wjets_data-num_dy_data-num_other_mc;
  float sb_ww_meas_err_data = sqrt(sideband_data + num_top_err_data*num_top_err_data + num_wjets_err_data*num_wjets_err_data + num_dy_err_data*num_dy_err_data + num_other_err_mc*num_other_err_mc);

  float sr_ww_meas_data = sb_ww_meas_data*rio_mg*eff_mt_mr_mg*eff_dp_mr_mg;
  float e1 = rio_mg*eff_mt_mr_mg*eff_dp_mr_mg*rio_mg*eff_mt_mr_mg*eff_dp_mr_mg*sb_ww_meas_err_data*sb_ww_meas_err_data;
  float e2  = sb_ww_meas_data*rio_mg*eff_mt_mr_mg*sb_ww_meas_data*rio_mg*eff_mt_mr_mg*eff_err_dp_mr_mg*eff_err_dp_mr_mg;
  float e3 = sb_ww_meas_data*eff_mt_mr_mg*eff_dp_mr_mg*sb_ww_meas_data*eff_mt_mr_mg*eff_dp_mr_mg*rio_err_mg*rio_err_mg;
  float e4 = sb_ww_meas_data*rio_mg*eff_dp_mr_mg*sb_ww_meas_data*rio_mg*eff_dp_mr_mg*eff_err_mt_mr_mg*eff_err_mt_mr_mg;
  float sr_ww_meas_err_data = sqrt( e1+e2+e3+e4 );

  if (printAll) {

    float scale_ttbar = getScale1fb(dir_mc+"ttbar")*lumi;
    float scale_tw    = getScale1fb(dir_mc+"tw")*lumi;
    float scale_wjets = getScale1fb(dir_mc+"wjets")*lumi;
    float scale_dytt  = getScale1fb(dir_mc+"dytt")*lumi;
    float scale_wz    = getScale1fb(dir_mc+"wz")*lumi;
    float scale_zz    = getScale1fb(dir_mc+"zz")*lumi;

    float sideband_ttbar = getYield(dir_mc+"ttbar", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
    float sideband_tw    = getYield(dir_mc+"tw", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
    cout << "sb top data, mc: " << num_top_data << "+/-" << num_top_err_data << " " << sideband_ttbar+sideband_tw << "+/-" << sqrt(scale_ttbar*sideband_ttbar+scale_tw*sideband_tw) << endl;
    float sideband_wjets = getYield(dir_mc+"wjets", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
    cout << "sb w+jets data, mc: " << num_wjets_data << "+/-" << num_wjets_err_data << " " << sideband_wjets << "+/-" << sqrt(scale_wjets*sideband_wjets) << endl;
    cout << "sb mc dymm, dyee, dytt, wz, zz: " << sideband_dymm << "+/-" << sideband_dymm_err << " " 
	                                       << sideband_dyee << "+/-" << sideband_dyee_err << " " 
	                                       << sideband_dytt << "+/-" << sqrt(scale_dytt*sideband_dytt) << " " 
	                                       << sideband_wz << "+/-" << sqrt(scale_wz*sideband_wz) << " " 
	                                       << sideband_zz << "+/-" << sqrt(scale_zz*sideband_zz) << endl;    
    float sideband_qqww_l        = getYield(dir_mc+"qqww", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
    float sideband_ggww_l        = getYield(dir_mc+"ggww", wwSelection,     veto, mass, njets,  "sideband,dphijet,minmet40", lumi).first;
    cout << "sb data, all mc: " << sideband_data << " " 
	 << sideband_qqww_l+sideband_ggww_l+sideband_ttbar+sideband_tw+sideband_wjets+sideband_dyee+sideband_dymm+sideband_zz+sideband_wz+sideband_dytt 
	 << endl;
    cout << "sb ww data, mc: " << sb_ww_meas_data << "+/-" << sb_ww_meas_err_data << " " 
	 << sideband_qqww_l+sideband_ggww_l << "+/-" << sqrt(scale_qq*sideband_qqww_l+scale_gg*sideband_ggww_l) << endl;
    pair<float,float> wwmc = wwEstimationMC(mass, njets, lumi);
    cout << "hr ww data, mc: " << sr_ww_meas_data << "+/-" << sr_ww_meas_err_data << " " << wwmc.first << "+/-" << wwmc.second << endl;
  }

  return make_pair<float,float>(sr_ww_meas_data,sr_ww_meas_err_data);
}

void makeWWTable(float lumi=1./*fb-1*/, bool doLatex=false) {

  pair<float, float> tagEff0j = topVetoEffEstimation(120,0,lumi);
  pair<float, float> tagEff1j = topVetoEffEstimation(120,1,lumi);

  //   cout << "topVeto eff 0j: " << tagEff0j.first << "+/-" << tagEff0j.second << endl;
  //   cout << "topVeto eff 1j: " << tagEff1j.first << "+/-" << tagEff1j.second << endl;

  int masses[] = {115,120,130,140,150,160,170,180,190,200};
  //int masses[] = {115,120,130,140,150};
  //int masses[] = {160,170,180,190,200};
  //int masses[] = {160};
  int nmasses = sizeof(masses)/sizeof(int);

  if (!doLatex) cout << "| mH  |   0-j meas    |   0-j exp     |   1-j meas   |   1-j exp    |" << endl;
  for (int j=0;j<nmasses;++j) {

    int mass = masses[j];
    pair<float,float> j0dd = wwEstimationData(mass,0,lumi,tagEff0j.first,tagEff0j.second);
    pair<float,float> j0mc = wwEstimationMC  (mass,0,lumi);
    pair<float,float> j1dd = wwEstimationData(mass,1,lumi,tagEff1j.first,tagEff1j.second);
    pair<float,float> j1mc = wwEstimationMC  (mass,1,lumi);

    if (doLatex) {
      cout << Form("%i & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %.1f \\\\",mass,
		   round(j0dd.first*10.)/10.,round(j0dd.second*10.)/10.,
		   round(j0mc.first*10.)/10.,round(j0mc.second*10.)/10.,
		   round(j1dd.first*10.)/10.,round(j1dd.second*10.)/10.,
		   round(j1mc.first*10.)/10.,round(j1mc.second*10.)/10.) 
	   << endl;
    } else {
      cout << Form("| %i | %4.1f +/- %4.1f | %4.1f +/- %4.1f | %4.1f +/- %4.1f | %4.1f +/- %4.1f |",mass,
		   round(j0dd.first*10.)/10.,round(j0dd.second*10.)/10.,
		   round(j0mc.first*10.)/10.,round(j0mc.second*10.)/10.,
		   round(j1dd.first*10.)/10.,round(j1dd.second*10.)/10.,
		   round(j1mc.first*10.)/10.,round(j1mc.second*10.)/10.) 
	   << endl;
    }
    
  }

}
