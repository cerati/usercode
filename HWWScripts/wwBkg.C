#include "common.C"
#include "topBg.C"

pair<float,float> wwEstimationMC(int mass=160, unsigned int njets=0, float lumi = 1./*fb-1*/, 
				 bool applyEff=true, bool doPUw=true){

  bool printAll = 0;

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  float scale_qq = getScale1fb(main_dir+topww_dir+"qqww")*lumi;
  float scale_gg = getScale1fb(main_dir+topww_dir+"ggww")*lumi;
  float sigreg_qqww_l = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float sigreg_ggww_l = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float sr_ww_true = sigreg_qqww_l+sigreg_ggww_l;
  float sr_ww_true_err = sqrt(scale_qq*sigreg_qqww_l+scale_gg*sigreg_ggww_l);

  if (printAll) cout << "ww MC: " << sr_ww_true << "+/-" << sr_ww_true_err << endl;

  return make_pair<float,float>(sr_ww_true,sr_ww_true_err);

}

float ratioPoissErr(float nval, float nerr, float dval, float derr) {
  return sqrt( pow(nerr/dval,2) + pow(nval*derr/pow(dval,2),2) );
}

float efficiencyErr(float eff, float den) {
  return sqrt( eff*(1-eff)/den );
}


pair<float,float> getAllBkg(int mass=160, unsigned int njets=0, string region="", float lumi = 1./*fb-1*/, float eff_veto_data=0, float eff_err_veto_data=0, 
			    bool useJson=false, bool applyEff=true, bool doPUw=true) {

  bool printAll = 0;

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  //********* get top **********
  pair<float, float> top = topBgEstimation(mass, njets, lumi, region, eff_veto_data, eff_err_veto_data, useJson, applyEff, false, doPUw);
  bool useTopSF = false;
  if (useTopSF) {
    //check error...
    float topSF = TopBkgScaleFactor(njets);
    pair<float, float> tt = getYield(main_dir+topww_dir+"ttbar", wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45", lumi, false, applyEff, false, doPUw);
    pair<float, float> tw = getYield(main_dir+topww_dir+"tw",    wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45", lumi, false, applyEff, false, doPUw);
    top = make_pair<float, float>(topSF*(tt.first+tw.first),sqrt( pow(topSF*tt.second,2) + pow(topSF*tw.second,2) ));
  }
  float num_top_data = top.first;
  float num_top_err_data = top.second;
  //********* get top **********

  //********* get w+jets **********
  pair<float, float> wj = fakeBgEstimationWithSyst(main_dir+topww_dir,wwSelection, veto, mass, njets, region, lumi, useJson, applyEff, doPUw);
  float num_wjets_data = wj.first;
  float num_wjets_err_data = wj.second;
  //********* get w+jets **********

  //********* get DY ********** => take it from MC
  //pair<float, float> dy = dyBkgEstimation(main_dir+topww_dir+"data.root", wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", 0);
  //float num_dy_data = dy.first;
  //float num_dy_err_data = dy.second;
  pair<float, float> sb_dymm   = getYield(main_dir+topww_dir+"dymm",  wwSelection,     veto, mass, njets,  region, lumi, false, applyEff, false, doPUw);
  pair<float, float> sb_dyee   = getYield(main_dir+topww_dir+"dyee",  wwSelection,     veto, mass, njets,  region, lumi, false, applyEff, false, doPUw);
  float sideband_dymm          = sb_dymm.first;
  float sideband_dyee          = sb_dyee.first;
  //float sideband_dymm_err      = sb_dymm.second;
  //float sideband_dyee_err      = sb_dyee.second;
  //get scale factors for DY
  float dySF = DYBkgScaleFactor(0,njets);
  float num_dy_data = (sideband_dyee+sideband_dymm)*dySF;
  float num_dy_err_data = num_dy_data;
  //********* get DY **********

  //********* get others **********
  pair<float, float> sideband_dytt = getYield(main_dir+topww_dir+"dytt",  wwSelection,     veto, mass, njets,  region, lumi, false, applyEff, false, doPUw);
  pair<float, float> sideband_zz   = getYield(main_dir+topww_dir+"zz_py",    wwSelection,     veto, mass, njets,  region, lumi, false, applyEff, false, doPUw);
  pair<float, float> sideband_wz   = getYield(main_dir+topww_dir+"wz",    wwSelection,     veto, mass, njets,  region, lumi, false, applyEff, false, doPUw);
  float num_other_mc  = sideband_zz.first+sideband_wz.first+sideband_dytt.first;
  float num_other_err_mc = num_other_mc*0.2;//20% for now
  //********* get others **********

  if (printAll) {
    cout << Form("data driven backgrounds all: %5.1f +/- %5.1f ; top: %5.1f +/- %5.1f ; fake: %5.1f +/- %5.1f ; dy: %5.1f +/- %5.1f ; oth: %5.1f +/- %5.1f . ",
		 num_top_data+num_wjets_data+num_dy_data+num_other_mc,sqrt(pow(num_top_err_data,2)+pow(num_wjets_err_data,2)+pow(num_dy_err_data,2)+pow(num_other_err_mc,2)),
		 num_top_data,num_top_err_data,
		 num_wjets_data,num_wjets_err_data,
		 num_dy_data,num_dy_err_data,
		 num_other_mc,num_other_err_mc) << endl;

    pair<float,float> sb_ttbar  = getYield(main_dir+topww_dir+"ttbar", wwSelection, veto, mass, njets, region, lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_tw     = getYield(main_dir+topww_dir+"tw",    wwSelection, veto, mass, njets, region, lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_wjets  = getYield(main_dir+topww_dir+"wjets", wwSelection, veto, mass, njets, region, lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_dytt   = getYield(main_dir+topww_dir+"dytt",  wwSelection, veto, mass, njets, region, lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_zz     = getYield(main_dir+topww_dir+"zz_py",    wwSelection, veto, mass, njets, region, lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_wz     = getYield(main_dir+topww_dir+"wz",    wwSelection, veto, mass, njets, region, lumi, false, applyEff, false, doPUw);

    cout << Form("MC backgrounds (no sf) all: %5.1f +/- %5.1f ; top: %5.1f +/- %5.1f ; fake: %5.1f +/- %5.1f ; dy: %5.1f +/- %5.1f ; oth: %5.1f +/- %5.1f . ",
		 sb_ttbar.first+sb_tw.first+sb_wjets.first+sb_dymm.first+sb_dyee.first+sb_dytt.first+sb_zz.first+sb_wz.first, 
		 sqrt(pow(sb_ttbar.second,2)+pow(sb_tw.second,2)+pow(sb_wjets.second,2)+pow(sb_dymm.second,2)+pow(sb_dyee.second,2)+pow(sb_dytt.second,2)+pow(sb_zz.second,2)+pow(sb_wz.second,2)), 
		 sb_ttbar.first+sb_tw.first,sqrt(pow(sb_ttbar.second,2)+pow(sb_tw.second,2)),
		 sb_wjets.first,sb_wjets.second,
		 sb_dymm.first+sb_dyee.first,sqrt(pow(sb_dymm.second,2)+pow(sb_dyee.second,2)),
		 sb_dytt.first+sb_zz.first+sb_wz.first,sqrt(pow(sb_dytt.second,2)+pow(sb_zz.second,2)+pow(sb_wz.second,2))) << endl;

  }

  return make_pair<float,float>(num_top_data+num_wjets_data+num_dy_data+num_other_mc,sqrt(pow(num_top_err_data,2)+pow(num_wjets_err_data,2)+pow(num_dy_err_data,2)+pow(num_other_err_mc,2)));

}

pair<float,float> wwEstimationData(int mass=160, unsigned int njets=0, float lumi = 1./*fb-1*/, float eff_veto_data=0, float eff_err_veto_data=0, 
				   bool useJson=false, bool applyEff=true, bool doPUw=true) {

  bool printAll = 0;

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  pair<float, float> sideband_qqww = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
  pair<float, float> sideband_ggww = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);

  pair<float, float> massreg_qqww = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "massreg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
  pair<float, float> massreg_ggww = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "massreg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);

  pair<float, float> mtreg_qqww = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "mtreg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
  pair<float, float> mtreg_ggww = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "mtreg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);

  pair<float, float> sigreg_qqww = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
  pair<float, float> sigreg_ggww = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "dphireg,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);

  float scale_qq  = getScale1fb(main_dir+topww_dir+"qqww")*lumi;
  //float scale_gg  = getScale1fb(main_dir+topww_dir+"ggww")*lumi;
  //float weight_gg = scale_gg/scale_qq;

  float rio_mg = (massreg_qqww.first+massreg_ggww.first)/(sideband_qqww.first+sideband_ggww.first);
  //ok, this is poissonian (not an eff)
  float rio_err_mg = ratioPoissErr( massreg_qqww.first+massreg_ggww.first, sqrt(pow(massreg_qqww.second,2)+pow(massreg_ggww.second,2)),
				    sideband_qqww.first+sideband_ggww.first, sqrt(pow(sideband_qqww.second,2)+pow(sideband_ggww.second,2)) );
  float eff_mt_mr_mg = (mtreg_qqww.first+mtreg_ggww.first)/(massreg_qqww.first+massreg_ggww.first);
  //use binomial, it's an efficiency
  float eff_err_mt_mr_mg = efficiencyErr( eff_mt_mr_mg, (massreg_qqww.first+massreg_ggww.first)/scale_qq  );//assume sample size is qqww
  //float eff_err_mt_mr_mg = ratioPoissErr( mtreg_qqww.first+mtreg_ggww.first, sqrt(pow(mtreg_qqww.second,2)+pow(mtreg_ggww.second,2)),
  //                                        massreg_qqww.first+massreg_ggww.first, sqrt(pow(massreg_qqww.second,2)+pow(massreg_ggww.second,2)) );
  float eff_dp_mr_mg = (sigreg_qqww.first+sigreg_ggww.first)/(mtreg_qqww.first+mtreg_ggww.first);
  //use binomial, it's an efficiency
  float eff_err_dp_mr_mg = efficiencyErr( eff_dp_mr_mg, (mtreg_qqww.first+mtreg_ggww.first)/scale_qq  );//assume sample size is qqww
  //float eff_err_dp_mr_mg = ratioPoissErr( sigreg_qqww.first+sigreg_ggww.first, sqrt(pow(sigreg_qqww.second,2)+pow(sigreg_ggww.second,2)),
  //                                        mtreg_qqww.first+mtreg_ggww.first, sqrt(pow(mtreg_qqww.second,2)+pow(mtreg_ggww.second,2)) );

  float sideband_data = getYield(main_dir+topww_dir+"data.root", wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", 0., useJson, false, false, false).first;

  pair<float,float> sb_bkg = getAllBkg(mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, eff_veto_data, eff_err_veto_data, useJson, applyEff, doPUw);

  float sb_ww_meas_data = sideband_data-sb_bkg.first;
  float sb_ww_meas_err_data = sqrt(sideband_data + pow(sb_bkg.second,2));

  float sr_ww_meas_data = sb_ww_meas_data*rio_mg*eff_mt_mr_mg*eff_dp_mr_mg;
  float e1 = rio_mg*eff_mt_mr_mg*eff_dp_mr_mg*rio_mg*eff_mt_mr_mg*eff_dp_mr_mg*sb_ww_meas_err_data*sb_ww_meas_err_data;
  float e2  = sb_ww_meas_data*rio_mg*eff_mt_mr_mg*sb_ww_meas_data*rio_mg*eff_mt_mr_mg*eff_err_dp_mr_mg*eff_err_dp_mr_mg;
  float e3 = sb_ww_meas_data*eff_mt_mr_mg*eff_dp_mr_mg*sb_ww_meas_data*eff_mt_mr_mg*eff_dp_mr_mg*rio_err_mg*rio_err_mg;
  float e4 = sb_ww_meas_data*rio_mg*eff_dp_mr_mg*sb_ww_meas_data*rio_mg*eff_dp_mr_mg*eff_err_mt_mr_mg*eff_err_mt_mr_mg;
  float sr_ww_meas_err_data = sqrt( e1+e2+e3+e4 );

  if (printAll) {
    pair<float,float> sb_qqww_l = getYield(main_dir+topww_dir+"qqww",  wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_ggww_l = getYield(main_dir+topww_dir+"ggww",  wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_ttbar  = getYield(main_dir+topww_dir+"ttbar", wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_tw     = getYield(main_dir+topww_dir+"tw",    wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_wjets  = getYield(main_dir+topww_dir+"wjets", wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_dymm   = getYield(main_dir+topww_dir+"dymm",  wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_dyee   = getYield(main_dir+topww_dir+"dyee",  wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_dytt   = getYield(main_dir+topww_dir+"dytt",  wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_zz     = getYield(main_dir+topww_dir+"zz_py",    wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    pair<float,float> sb_wz     = getYield(main_dir+topww_dir+"wz",    wwSelection, veto, mass, njets,  "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw);
    cout << "sb data, all mc(no sf): " << sideband_data << " " 
	 << sb_qqww_l.first+sb_ggww_l.first+sb_ttbar.first+sb_tw.first+sb_wjets.first+sb_dyee.first+sb_dymm.first+sb_zz.first+sb_wz.first+sb_dytt.first 
	 << endl;
    cout << "sb ww data, mc: " << sb_ww_meas_data << "+/-" << sb_ww_meas_err_data << " " 
	 << sb_qqww_l.first+sb_ggww_l.first << "+/-" << sqrt(pow(sb_qqww_l.second,2)+pow(sb_ggww_l.second,2)) << endl;
    pair<float,float> wwmc = wwEstimationMC(mass, njets, lumi, applyEff, doPUw);
    cout << "hr ww data, mc: " << sr_ww_meas_data << "+/-" << sr_ww_meas_err_data << " " << wwmc.first << "+/-" << wwmc.second << endl;
    cout << "ratio (signal/side): " << rio_mg << "+/-" << rio_err_mg << endl;
    cout << "sr eff mt, dphi: " << eff_mt_mr_mg << "+/-" << eff_err_mt_mr_mg << " " << eff_dp_mr_mg << "+/-" << eff_err_dp_mr_mg << endl;
  }

  return make_pair<float,float>(sr_ww_meas_data,sr_ww_meas_err_data);
}

void makeWWTable(float lumi=1./*fb-1*/, bool doLatex=false) {

  bool useJson  = true;
  bool applyEff = true;
  bool doPUw    = true;

  pair<float, float> tagEff0j = topVetoEffEstimation(0,0,lumi,"dphijet,minmetvtx,lep2pt15,ptll45,mll20",useJson,applyEff,doPUw);
  pair<float, float> tagEff1j = topVetoEffEstimation(0,1,lumi,"dphijet,minmetvtx,lep2pt15,ptll45,mll20",useJson,applyEff,doPUw);

  //   cout << "topVeto eff 0j: " << tagEff0j.first << "+/-" << tagEff0j.second << endl;
  //   cout << "topVeto eff 1j: " << tagEff1j.first << "+/-" << tagEff1j.second << endl;

  int masses[] = {115,120,130,140,150,160,170,180,190};
  //int masses[] = {115,120,130,140,150};
  //int masses[] = {115,130,150,170,190};
  //int masses[] = {120,140,160,180,200};
  //int masses[] = {115};
  int nmasses = sizeof(masses)/sizeof(int);
  doLatex=false;
  if (!doLatex) {
    cout << "-------------------------------------------------------------------------------------------------------------" << endl;
    cout << "| m_H |    0-j meas    |    0-j exp     |       SF       |    1-j meas    |    1-j exp     |       SF       |" << endl;
    cout << "-------------------------------------------------------------------------------------------------------------" << endl;
  }

  vector<float> vsf0j;
  vector<float> vk0j;
  vector<float> vsf1j;
  vector<float> vk1j;

  for (int j=0;j<nmasses;++j) {

    int mass = masses[j];
    pair<float,float> j0dd = wwEstimationData(mass,0,lumi,tagEff0j.first,tagEff0j.second,useJson,applyEff,doPUw);
    pair<float,float> j0mc = wwEstimationMC  (mass,0,lumi,applyEff,doPUw);
    pair<float,float> j1dd = wwEstimationData(mass,1,lumi,tagEff1j.first,tagEff1j.second,useJson,applyEff,doPUw);
    pair<float,float> j1mc = wwEstimationMC  (mass,1,lumi,applyEff,doPUw);

    float sf0j  = j0dd.first/j0mc.first;
    float sf0je = sqrt(pow(j0dd.second/j0mc.first,2)+pow(j0dd.first*j0mc.second/pow(j0mc.first,2),2));
    float sf1j  = j1dd.first/j1mc.first;
    float sf1je = sqrt(pow(j1dd.second/j1mc.first,2)+pow(j1dd.first*j1mc.second/pow(j1mc.first,2),2));

    TString raw = "| %i | %5.1f +/- %4.1f | %5.1f +/- %4.1f | %5.2f +/- %4.2f | %5.1f +/- %4.1f | %5.1f +/- %4.1f | %5.2f +/- %4.2f |";
    if (doLatex) raw = "%i & %5.1f $\\pm$ %4.1f & %5.1f $\\pm$ %4.1f & %4.2f $\\pm$ %4.2f & %5.1f $\\pm$ %4.1f & %5.1f $\\pm$ %4.1f & %4.2f $\\pm$ %4.2f \\\\";
    cout << Form(raw.Data(),mass,
		 round(j0dd.first*10.)/10.,round(j0dd.second*10.)/10.,
		 round(j0mc.first*10.)/10.,round(j0mc.second*10.)/10.,
		 round(sf0j*100)/100.,
		 round(sf0je*100)/100.,
		 round(j1dd.first*10.)/10.,round(j1dd.second*10.)/10.,
		 round(j1mc.first*10.)/10.,round(j1mc.second*10.)/10.,
		 round(sf1j*100)/100.,
		 round(sf1je*100)/100.) 
	 << endl;
    vsf0j.push_back(sf0j);
    vk0j.push_back(1.+sf0je/sf0j);
    vsf1j.push_back(sf1j);
    vk1j.push_back(1.+sf1je/sf1j);    
  }

  if (!doLatex) cout << "-------------------------------------------------------------------------------------------------------------" << endl;

  if (nmasses==9) {

    ofstream myfile;
    TString fname = "WWBkgScaleFactors.h";
    if (!doMVA) myfile.open(fname);
    else myfile.open(fname,ios::app);
    ostream &out = myfile;

    if (!doMVA) out << Form("Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {\n");
    else out << Form("Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {\n");
    out << Form("assert(jetBin >= 0 && jetBin <= 1);\n");
    out << Form("  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};\n");
    out << Form("  Double_t WWBkgScaleFactorHiggsSelection[2][9] = { \n");
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f},\n",vsf0j[0],vsf0j[1],vsf0j[2],vsf0j[3],vsf0j[4],vsf0j[5],vsf0j[6],vsf0j[7],vsf0j[8]);
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f} };\n",vsf1j[0],vsf1j[1],vsf1j[2],vsf1j[3],vsf1j[4],vsf1j[5],vsf1j[6],vsf1j[7],vsf1j[8]);
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < 9 ; ++m) {\n");
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return 1.0;\n");
    out << Form("  }\n");
    out << Form("}\n");
    
    if (!doMVA) out << Form("Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {\n");
    else out << Form("Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {\n");
    out << Form("assert(jetBin >= 0 && jetBin <= 1);\n");
    out << Form("  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};\n");
    out << Form("  Double_t WWBkgScaleFactorKappaHiggsSelection[2][9] = { \n");
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f},\n",vk0j[0],vk0j[1],vk0j[2],vk0j[3],vk0j[4],vk0j[5],vk0j[6],vk0j[7],vk0j[8]);
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f} };\n",vk1j[0],vk1j[1],vk1j[2],vk1j[3],vk1j[4],vk1j[5],vk1j[6],vk1j[7],vk1j[8]);
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < 9 ; ++m) {\n");
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return 1.0;\n");
    out << Form("  }\n");
    out << Form("}\n");
  }

}

void wwBkg(float lumi=1./*fb-1*/, bool doLatex=false) {
  doMVA = false;
  makeWWTable(lumi, doLatex);
  doMVA = true;
  makeWWTable(lumi, doLatex);
}

void printCutEffic(int mass=160, unsigned int njets=0, float lumi = 1./*fb-1*/, 
		   bool applyEff=true, bool doPUw=true){

  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  float scale_qq = getScale1fb(main_dir+topww_dir+"qqww")*lumi;
  float scale_gg = getScale1fb(main_dir+topww_dir+"ggww")*lumi;

  float sb_qqww_l = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float sb_ggww_l = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float sb_ww_true = sb_qqww_l+sb_ggww_l;
  float sb_ww_true_err = sqrt(scale_qq*sb_qqww_l+scale_gg*sb_ggww_l);

  float mt_qqww_l = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "mtside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float mt_ggww_l = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "mtside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float mt_ww_true = mt_qqww_l+mt_ggww_l;
  float mt_ww_true_err = sqrt(scale_qq*mt_qqww_l+scale_gg*mt_ggww_l);

  float dp_qqww_l = getYield(main_dir+topww_dir+"qqww", wwSelection, veto, mass, njets, "dphiside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float dp_ggww_l = getYield(main_dir+topww_dir+"ggww", wwSelection, veto, mass, njets, "dphiside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, false, applyEff, false, doPUw).first;
  float dp_ww_true = dp_qqww_l+dp_ggww_l;
  float dp_ww_true_err = sqrt(scale_qq*dp_qqww_l+scale_gg*dp_ggww_l);

  pair<float,float> sb_data = getYield(main_dir+topww_dir+"data.root", wwSelection, veto, mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", 0, false, false, false, false);
  pair<float,float> mt_data = getYield(main_dir+topww_dir+"data.root", wwSelection, veto, mass, njets, "mtside,dphijet,minmetvtx,lep2pt15,ptll45,mll20",   0, false, false, false, false);
  pair<float,float> dp_data = getYield(main_dir+topww_dir+"data.root", wwSelection, veto, mass, njets, "dphiside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", 0, false, false, false, false);

  pair<float,float> sb_bkg = getAllBkg(mass, njets, "sideband,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, 0,0, false, applyEff, doPUw);
  pair<float,float> mt_bkg = getAllBkg(mass, njets, "mtside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, 0,0, false, applyEff, doPUw);
  pair<float,float> dp_bkg = getAllBkg(mass, njets, "dphiside,dphijet,minmetvtx,lep2pt15,ptll45,mll20", lumi, 0,0, false, applyEff, doPUw);

  pair<float,float> sb_ww = make_pair<float,float>(sb_data.first-sb_bkg.first,sqrt( pow(sb_data.second,2) + pow(sb_bkg.second,2) ));
  pair<float,float> mt_ww = make_pair<float,float>(mt_data.first-mt_bkg.first,sqrt( pow(mt_data.second,2) + pow(mt_bkg.second,2) ));
  pair<float,float> dp_ww = make_pair<float,float>(dp_data.first-dp_bkg.first,sqrt( pow(dp_data.second,2) + pow(dp_bkg.second,2) ));

  cout << "------------------------------------------------------------------------- mH=" << mass << " ------------------------------------------------------------------------" << endl;
  cout << "|   region   |      yield      |   N-1 effic(%)   ||   region   |      yield      |   N-1 effic(%)   ||   region   |      yield      |   N-1 effic(%)   |" << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << Form("|  ww MC sb  | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", sb_ww_true,sb_ww_true_err, 
	       100*sb_ww_true/sb_ww_true,100*efficiencyErr(sb_ww_true/sb_ww_true,sb_ww_true/scale_qq));//assume error if only qqww
  cout << Form("|   data sb  | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", sb_data.first,sb_data.second,
	       100*sb_data.first/sb_data.first,100*efficiencyErr(sb_data.first/sb_data.first,sb_data.first));
  cout << Form("| ww data sb | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", sb_ww.first,sb_ww.second,
	       100*sb_ww.first/sb_ww.first,100*efficiencyErr(sb_ww.first/sb_ww.first,sb_ww.first)) << endl;
  cout << Form("|  ww MC mt  | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", mt_ww_true,mt_ww_true_err, 
	       100*mt_ww_true/sb_ww_true,100*efficiencyErr(mt_ww_true/sb_ww_true,sb_ww_true/scale_qq));//assume error if only qqww
  cout << Form("|   data mt  | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |",mt_data.first ,mt_data.second, 
	       100*mt_data.first/sb_data.first,100*efficiencyErr(mt_data.first/sb_data.first,sb_data.first));
  cout << Form("| ww data mt | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |",mt_ww.first ,mt_ww.second, 
	       100*mt_ww.first/sb_ww.first,100*efficiencyErr(mt_ww.first/sb_ww.first,sb_ww.first)) << endl;
  cout << Form("|  ww MC dp  | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", dp_ww_true,dp_ww_true_err, 
	       100*dp_ww_true/mt_ww_true,100*efficiencyErr(dp_ww_true/mt_ww_true,mt_ww_true/scale_qq));//assume error if only qqww
  cout << Form("|   data dp  | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", dp_data.first,dp_data.second, 
	       100*dp_data.first/mt_data.first,100*efficiencyErr(dp_data.first/mt_data.first,mt_data.first));
  cout << Form("| ww data dp | %5.1f +/- %5.1f |  %5.1f +/- %5.1f |", dp_ww.first,dp_ww.second, 
	       100*dp_ww.first/mt_ww.first,100*efficiencyErr(dp_ww.first/mt_ww.first,mt_ww.first)) << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

  //cout << scale_qq << endl;
  //cout << efficiencyErr(mt_ww.first/sb_ww.first,sb_ww.first) << endl;
  //cout << ratioPoissErr(mt_ww.first,sqrt(mt_ww.first),sb_ww.first,sqrt(sb_ww.first)) << endl;
  //cout << ratioPoissErr(mt_ww.first,mt_ww.second,sb_ww.first,sb_ww.second) << endl;
  return;

}
