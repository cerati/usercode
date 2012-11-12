#include "common.C"
#include "fakeBg.C"

pair<float, float> evaluateBackground(TString dir, unsigned int cut, unsigned int veto, int mass, int njets, TString myRegion, float lumi, bool useJson=0, 
				      bool applyEff=true, bool doFake=false, bool doPUw=false){


  bool debug = 0;

  //get scale factors for DY
  float dySF = DYBkgScaleFactor(0,njets);

  float qqww = getYield(dir+"qqww",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  float ggww = getYield(dir+"ggww",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  float dyll = dySF*getYield(dir+"dyll",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  //float dyee = dySF*getYield(dir+"dyee",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  float dytt = 0;//getYield(dir+"dytt",  cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  float zz   = getYield(dir+"zz",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  float wz   = getYield(dir+"wz",    cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;
  float wg   = getYield(dir+"wgamma",cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doFake, doPUw).first;

  float mc_other = qqww + ggww + dyll/* + dyee*/ + dytt + zz+ wz+wg;
  pair<float, float> nwj = fakeBgEstimationWithSyst(main_dir+wj_dir,cut, veto, mass, njets, myRegion, lumi, useJson, applyEff, doPUw);

  float background = mc_other+nwj.first;
  float backgr_err = sqrt(pow(0.5*qqww,2) + pow(0.5*ggww,2) + pow(0.5*dyll,2)/* + pow(0.5*dyee,2)*/ + pow(0.5*dytt,2) + pow(0.5*zz,2)+ pow(0.5*wz,2)+pow(0.5*wg,2)+pow(nwj.second,2));//take 50% for MC based
  if (debug) cout << "bkg total, fake, qqww, ggww, dyll"<</*, dyee*/", dytt, zz, wz, wgamma: " << background << " " << nwj.first << " " 
		  << qqww << " " << ggww << " " << dyll/* << " " << dyee*/ << " " << dytt << " " << zz << " " << wz << " " << wg << endl;

  return make_pair<float, float>(background,backgr_err);
}

pair<float, float> topVetoEffEstimation(int mass=160, unsigned int njets=0, float lumi=0., TString region = "=dphijet=minmet40=", 
					bool useJson=0, bool applyEff=true, bool doPUw=false, float sfTW=1.) {

  //mass is not used
  bool debug = 0;
  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  TString region_top    = region+"=btagJet1=";//test noSoftMu= 
  TString region_toptag = region+"=btagJet1=";//test noSoftMu=
  if (njets==1) {
    region_top    = region+"=btagJet2=";
    region_toptag = region+"=btag1and2=";
  }

  //get data yields
  float sideband_data_nj_top     = getYield(main_dir+topww_dir+"data.root", control_top,    veto, mass, nj_top, region_top,    0., useJson, false, false, false).first;
  float sideband_data_nj_top_tag = getYield(main_dir+topww_dir+"data.root", control_toptag, veto, mass, nj_top, region_toptag, 0., useJson, false, false, false).first;

  //subtract contamination
  pair<float, float> sb_nj_top_bkg     = evaluateBackground(main_dir+topww_dir, control_top,    veto, mass, nj_top, region_top,    lumi, false, applyEff, false, doPUw);
  pair<float, float> sb_nj_top_tag_bkg = evaluateBackground(main_dir+topww_dir, control_toptag, veto, mass, nj_top, region_toptag, lumi, false, applyEff, false, doPUw);

  //need to subtract tW as well
  pair<float, float> sb_nj_tw     = getYield(main_dir+topww_dir+"tw",  control_top,    veto, mass, nj_top, region_top,    lumi, false, applyEff, false, doPUw);
  pair<float, float> sb_nj_tw_tag = getYield(main_dir+topww_dir+"tw",  control_toptag, veto, mass, nj_top, region_toptag, lumi, false, applyEff, false, doPUw);

  sb_nj_top_bkg = make_pair<float, float>(sb_nj_top_bkg.first+sfTW*sb_nj_tw.first,sb_nj_top_bkg.second);//FIXME: add uncertainties on tW
  sb_nj_top_tag_bkg = make_pair<float, float>(sb_nj_top_tag_bkg.first+sfTW*sb_nj_tw_tag.first,sb_nj_top_tag_bkg.second);//FIXME: add uncertainties on tW

  sideband_data_nj_top-=sb_nj_top_bkg.first;        //FIXME: add uncertainties on background
  sideband_data_nj_top_tag-=sb_nj_top_tag_bkg.first;//FIXME: add uncertainties on background

  //get tag and veto efficiencies
  float eff_tag_data = sideband_data_nj_top_tag/sideband_data_nj_top;
  float eff_err_tag_data = sqrt(eff_tag_data*(1.-eff_tag_data)/sideband_data_nj_top);//binomial
  //float eff_err_tag_data = eff_tag_data*sqrt(1./sideband_data_nj_top_tag + 1./sideband_data_nj_top);//poissonian
  float eff_veto_data = 0;
  float eff_err_veto_data = 0;
  float fttbar=0.,fttbar_err=0.,twlikett=0.;
  if (njets==0) {
    float wwntv1j_tw  = getYield(main_dir+topww_dir+"tw", wwSelNoMetNoTV,                 veto, mass, 1, region_top,    lumi, false, applyEff, false, doPUw).first;    
    float wwtt1j_tw   = getYield(main_dir+topww_dir+"tw", wwSelNoMetNoTV|TopTagNotInJets, veto, mass, 1, region_toptag, lumi, false, applyEff, false, doPUw).first;    
    twlikett = wwtt1j_tw/wwntv1j_tw;
    //float twlikett_err = 0.;//fixme
    float sideband_ttbar  = getYield(main_dir+topww_dir+"ttbar_powheg", wwSelNoMetNoTV, veto, mass, njets, region, lumi, false, applyEff, false, doPUw).first;
    float sideband_tw     = getYield(main_dir+topww_dir+"tw",    wwSelNoMetNoTV, veto, mass, njets, region, lumi, false, applyEff, false, doPUw).first;
    fttbar = (sideband_ttbar+sideband_tw*twlikett)/(sideband_ttbar+sideband_tw);
    fttbar_err = fttbar*0.17;//from cross section uncertainty
    eff_veto_data = fttbar*(1 - (1-eff_tag_data)*(1-eff_tag_data)) + (1.-fttbar)*eff_tag_data;
    eff_err_veto_data = pow(eff_tag_data-pow(eff_tag_data,2),2)*pow(fttbar_err,2) + pow(fttbar-2*fttbar*eff_tag_data+1,2)*pow(eff_err_tag_data,2);
    eff_err_veto_data = sqrt(eff_err_veto_data);
  } else if (njets==1) {
    eff_veto_data = eff_tag_data;
    eff_err_veto_data = eff_err_tag_data;
  }
  if (debug) {
    cout << "nj cr: " << nj_top << " twsf=" << sfTW << endl;
    if (njets==0) cout << "fttbar, twlikett: " << fttbar << " " << twlikett << endl;
    cout << "uncorr data cr, tag: " << sideband_data_nj_top+sb_nj_top_bkg.first << " " << sideband_data_nj_top_tag+sb_nj_top_tag_bkg.first << endl;
    cout << "backgound cr, tag: " << sb_nj_top_bkg.first << " " << sb_nj_top_tag_bkg.first << endl;
    cout << "tw bkg cr, tag: " << sfTW*sb_nj_tw.first << " " << sfTW*sb_nj_tw_tag.first << endl;
    cout << "corr. data cr, tag: " << sideband_data_nj_top << " " << sideband_data_nj_top_tag << endl;
    cout << "data eff tag,veto: " << eff_tag_data << "+/-" << eff_err_tag_data << " " << eff_veto_data << "+/-" << eff_err_veto_data << endl;
  }
  return make_pair(eff_veto_data, eff_err_veto_data);
}

pair<float, float> topBgEstimation(int mass=160, unsigned int njets=0, float lumi=0., TString region = "=dphijet=minmet40=", 
				   float eff_veto_data=0, float eff_err_veto_data=0, bool useJson=0, 
				   bool applyEff=true, bool doFake=false, bool doPUw=false, float sfTW=1.)  {
  bool debug = 0;
  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  TString toptagreg = "";//test =noSoftMu=
  if (njets==1) toptagreg += "=btagJet1=";

  pair<float, float> sb_data_tag = getYield(main_dir+topww_dir+"data.root", baseline_toptag, veto, mass, njets,  region+toptagreg, 0, useJson, 0, 0, 0);
  float sideband_data_tag  = sb_data_tag.first;

  //compute background contamination
  pair<float, float> sb_tag_bkg = evaluateBackground(main_dir+topww_dir, baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw);
  //subtract background
  sideband_data_tag-=sb_tag_bkg.first;
  float sideband_data_tag_err  = sqrt( pow(sb_data_tag.second,2) + pow(sb_tag_bkg.second,2) );

  //get the efficiency (if needed)
  if (eff_veto_data<1E-5) {
    pair<float, float> tagEff = topVetoEffEstimation(mass,njets,lumi,region,useJson,applyEff,doPUw,sfTW);
    eff_veto_data = tagEff.first;
    eff_err_veto_data = tagEff.second;
  }

  if (debug) cout << "using eff=" << eff_veto_data << " +/- " << eff_err_veto_data << endl;

  //compute the number of top events
  float num_top_data = sideband_data_tag*(1.-eff_veto_data)/eff_veto_data;
  float num_top_err_data = sqrt( pow((1.-eff_veto_data)/eff_veto_data,2)*pow(sideband_data_tag_err,2) + pow(sideband_data_tag,2)*pow(eff_err_veto_data,2)/pow(eff_veto_data,4) );
  if (debug) {
    //float sb_mc_top_tag = getYield(main_dir+topww_dir+"ttbar_powheg", baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first+
    //                      getYield(main_dir+topww_dir+"tw",    baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first;
    //float sb_mc_other_tag = getYield(main_dir+topww_dir+"qqww",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
    //                        getYield(main_dir+topww_dir+"ggww",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first +
    //                        getYield(main_dir+topww_dir+"wjets", baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
    //                        getYield(main_dir+topww_dir+"dyll",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first +
    //                        //getYield(main_dir+topww_dir+"dyee",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
    //                        //getYield(main_dir+topww_dir+"dytt",  baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first +//fixme
    //                        getYield(main_dir+topww_dir+"zz",    baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first + 
    //                        getYield(main_dir+topww_dir+"wz",    baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first;
    //cout << "sb tag mc top, mc other: " << sb_mc_top_tag << " " << sb_mc_other_tag << endl;
    cout << "sb tag data unc, bkg, data corr: " << sideband_data_tag+sb_tag_bkg.first << " " << sb_tag_bkg.first << " " << sideband_data_tag << endl;
    //cout << "sb_tag_purity: " << sb_tag_purity << endl;
    //cout << "sb tag data unc, exp bkg, data corr: " << sideband_data_tag/sb_tag_purity << " " << (sideband_data_tag/sb_tag_purity)*(1.-sb_tag_purity) << " " << sideband_data_tag << endl;
    float sigreg_ttbar  = getYield(main_dir+topww_dir+"ttbar_powheg", wwSelNoMet, veto, mass, njets,  region, lumi, useJson, applyEff, doFake, doPUw).first;
    float sigreg_tw     = getYield(main_dir+topww_dir+"tw",    wwSelNoMet, veto, mass, njets,  region, lumi, useJson, applyEff, doFake, doPUw).first;
    //cout << "data(notag)=" << num_top_data << " mc(notag)=" << sigreg_ttbar+sigreg_tw << " mc(tag)=" << sb_mc_top_tag << endl;
    cout << "data(notag)=" << num_top_data << " mc(notag)=" << sigreg_ttbar+sigreg_tw << endl;
  }
  return make_pair<float, float>(num_top_data, num_top_err_data);
}

void makeTopTable(float lumi) {

  bool useJson  = true;
  bool applyEff = true;
  bool doFake   = false;
  bool doPUw    = true;

  TString anaRegion = "=dphijet=dymvacut=ptll45=";//noSoftMu=lep2pt20allfs=looseVBF=

  int mass = 0;

  bool printAll = 0;

  ///////////////////////////////////////// 2-JET BIN /////////////////////////////////////////
  //signal region
  doVBF=1;
  //pair<float, float> sigreg_ttbar_2j = getYield(main_dir+topww_dir+"ttbar_powheg",    wwSelNoMet, noVeto, mass, 2, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  //pair<float, float> sigreg_tw_2j    = getYield(main_dir+topww_dir+"tw",       wwSelNoMet, noVeto, mass, 2, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_ttbar_2j = getYield(main_dir+topww_dir+"ttbar_powheg",    wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet1=nobJet2=nobJet3=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_2j    = getYield(main_dir+topww_dir+"tw",       wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet1=nobJet2=nobJet3=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_top_2j = make_pair<float, float>(sigreg_ttbar_2j.first+sigreg_tw_2j.first,sqrt(pow(sigreg_ttbar_2j.second,2)+pow(sigreg_tw_2j.second,2)));
  TH2F *data_ctrtag_2j_h = new TH2F("data_ctrtag_2j_h","data_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",data_ctrtag_2j_h, main_dir+topww_dir+"data", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", 0, useJson, 0, 0, 0);
  TH2F *qqww_ctrtag_2j_h = new TH2F("qqww_ctrtag_2j_h","qqww_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",qqww_ctrtag_2j_h, main_dir+topww_dir+"qqww", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *ggww_ctrtag_2j_h = new TH2F("ggww_ctrtag_2j_h","ggww_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",ggww_ctrtag_2j_h, main_dir+topww_dir+"ggww", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *dyll_ctrtag_2j_h = new TH2F("dyll_ctrtag_2j_h","dyll_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",dyll_ctrtag_2j_h, main_dir+topww_dir+"dyll", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", lumi, useJson, applyEff, doFake, doPUw);
  dyll_ctrtag_2j_h->Scale(DYBkgScaleFactor(0,2));
  TH2F *zz_ctrtag_2j_h = new TH2F("zz_ctrtag_2j_h","zz_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",zz_ctrtag_2j_h, main_dir+topww_dir+"zz", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wz_ctrtag_2j_h = new TH2F("wz_ctrtag_2j_h","wz_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",wz_ctrtag_2j_h, main_dir+topww_dir+"wz", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wg_ctrtag_2j_h = new TH2F("wg_ctrtag_2j_h","wg_ctrtag_2j_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",wg_ctrtag_2j_h, main_dir+topww_dir+"wgamma", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F* wj_ctrtag_2j_h = new TH2F("wj_ctrtag_2j_h","wj_ctrtag_2j_h",5,0,2.5,1,0,200);
  //fillPlot("ctrjetetapt",wj_ctrtag_2j_h,main_dir+topww_dir+"data",   wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=", 0, useJson, 0, 1, 0);//fixme
  //   pair<float,float> wj_ctrtag_2j_bin1 = fakeBgEstimationWithSyst(main_dir+wj_dir, wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=ctrjetbin1=", lumi, useJson, applyEff, doPUw);
  //   pair<float,float> wj_ctrtag_2j_bin2 = fakeBgEstimationWithSyst(main_dir+wj_dir, wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=ctrjetbin2=", lumi, useJson, applyEff, doPUw);
  //   pair<float,float> wj_ctrtag_2j_bin3 = fakeBgEstimationWithSyst(main_dir+wj_dir, wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=ctrjetbin3=", lumi, useJson, applyEff, doPUw);
  //   pair<float,float> wj_ctrtag_2j_bin4 = fakeBgEstimationWithSyst(main_dir+wj_dir, wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=ctrjetbin4=", lumi, useJson, applyEff, doPUw);
  //   pair<float,float> wj_ctrtag_2j_bin5 = fakeBgEstimationWithSyst(main_dir+wj_dir, wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bTagNoFwdYesCtr=ctrjetbin5=", lumi, useJson, applyEff, doPUw);
  //   wj_ctrtag_2j_h->SetBinContent(1,1,wj_ctrtag_2j_bin1.first);
  //   wj_ctrtag_2j_h->SetBinContent(2,1,wj_ctrtag_2j_bin1.first);
  //   wj_ctrtag_2j_h->SetBinContent(3,1,wj_ctrtag_2j_bin1.first);
  //   wj_ctrtag_2j_h->SetBinContent(4,1,wj_ctrtag_2j_bin1.first);
  //   wj_ctrtag_2j_h->SetBinContent(5,1,wj_ctrtag_2j_bin1.first);
  data_ctrtag_2j_h->Add(qqww_ctrtag_2j_h,-1.);
  data_ctrtag_2j_h->Add(ggww_ctrtag_2j_h,-1.);
  data_ctrtag_2j_h->Add(dyll_ctrtag_2j_h,-1.);
  data_ctrtag_2j_h->Add(zz_ctrtag_2j_h,-1.);
  data_ctrtag_2j_h->Add(wz_ctrtag_2j_h,-1.);
  data_ctrtag_2j_h->Add(wg_ctrtag_2j_h,-1.);
  data_ctrtag_2j_h->Add(wj_ctrtag_2j_h,-1.);
  //control region
  doVBF=0;
  //these are the background to subtract from data (including tW)
  TH2F *qqww_ctrtag_2j_num_h = new TH2F("qqww_ctrtag_2j_num_h","qqww_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",qqww_ctrtag_2j_num_h, main_dir+topww_dir+"qqww", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *qqww_ctrtag_2j_den_h = new TH2F("qqww_ctrtag_2j_den_h","qqww_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",qqww_ctrtag_2j_den_h, main_dir+topww_dir+"qqww", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *ggww_ctrtag_2j_num_h = new TH2F("ggww_ctrtag_2j_num_h","ggww_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",ggww_ctrtag_2j_num_h, main_dir+topww_dir+"ggww", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *ggww_ctrtag_2j_den_h = new TH2F("ggww_ctrtag_2j_den_h","ggww_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",ggww_ctrtag_2j_den_h, main_dir+topww_dir+"ggww", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *dyll_ctrtag_2j_num_h = new TH2F("dyll_ctrtag_2j_num_h","dyll_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",dyll_ctrtag_2j_num_h, main_dir+topww_dir+"dyll", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  dyll_ctrtag_2j_num_h->Scale(DYBkgScaleFactor(0,2));
  TH2F *dyll_ctrtag_2j_den_h = new TH2F("dyll_ctrtag_2j_den_h","dyll_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",dyll_ctrtag_2j_den_h, main_dir+topww_dir+"dyll", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  dyll_ctrtag_2j_den_h->Scale(DYBkgScaleFactor(0,2));
  TH2F *zz_ctrtag_2j_num_h = new TH2F("zz_ctrtag_2j_num_h","zz_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",zz_ctrtag_2j_num_h, main_dir+topww_dir+"zz", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *zz_ctrtag_2j_den_h = new TH2F("zz_ctrtag_2j_den_h","zz_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",zz_ctrtag_2j_den_h, main_dir+topww_dir+"zz", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wz_ctrtag_2j_num_h = new TH2F("wz_ctrtag_2j_num_h","wz_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",wz_ctrtag_2j_num_h, main_dir+topww_dir+"wz", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wz_ctrtag_2j_den_h = new TH2F("wz_ctrtag_2j_den_h","wz_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",wz_ctrtag_2j_den_h, main_dir+topww_dir+"wz", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *tw_ctrtag_2j_num_h = new TH2F("tw_ctrtag_2j_num_h","tw_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",tw_ctrtag_2j_num_h, main_dir+topww_dir+"tw", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *tw_ctrtag_2j_den_h = new TH2F("tw_ctrtag_2j_den_h","tw_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",tw_ctrtag_2j_den_h, main_dir+topww_dir+"tw", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wg_ctrtag_2j_num_h = new TH2F("wg_ctrtag_2j_num_h","wg_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",wg_ctrtag_2j_num_h, main_dir+topww_dir+"wgamma", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wg_ctrtag_2j_den_h = new TH2F("wg_ctrtag_2j_den_h","wg_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",wg_ctrtag_2j_den_h, main_dir+topww_dir+"wgamma", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
  TH2F *wj_ctrtag_2j_num_h = new TH2F("wj_ctrtag_2j_num_h","wj_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  //fillPlot("ctrjetetapt",wj_ctrtag_2j_num_h, main_dir+topww_dir+"data", wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", 0, useJson, 0, 1, 0);//fixme
  TH2F *wj_ctrtag_2j_den_h = new TH2F("wj_ctrtag_2j_den_h","wj_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  //fillPlot("ctrjetetapt",wj_ctrtag_2j_den_h, main_dir+topww_dir+"data", wwSelNoMetLepTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", 0, useJson, 0, 1, 0);//fixme
  TH2F *other_ctrtag_2j_num_h = new TH2F("other_ctrtag_2j_num_h","other_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  other_ctrtag_2j_num_h->Add(qqww_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(ggww_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(dyll_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(zz_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(wz_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(tw_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(wg_ctrtag_2j_num_h);
  other_ctrtag_2j_num_h->Add(wj_ctrtag_2j_num_h);
  TH2F *other_ctrtag_2j_den_h = new TH2F("other_ctrtag_2j_den_h","other_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  other_ctrtag_2j_den_h->Add(qqww_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(ggww_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(dyll_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(zz_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(wz_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(tw_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(wg_ctrtag_2j_den_h);
  other_ctrtag_2j_den_h->Add(wj_ctrtag_2j_den_h);
  //data..
  TH2F *data_ctrtag_2j_num_h = new TH2F("data_ctrtag_2j_num_h","data_ctrtag_2j_num_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",data_ctrtag_2j_num_h, main_dir+topww_dir+"data", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", 0, useJson, 0, 0, 0);
  TH2F *data_ctrtag_2j_den_h = new TH2F("data_ctrtag_2j_den_h","data_ctrtag_2j_den_h",5,0,2.5,1,0,200);
  fillPlot("ctrjetetapt",data_ctrtag_2j_den_h, main_dir+topww_dir+"data", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", 0, useJson, 0, 0, 0);
  float yield = 0;
  float error = 0;
  for (int i=1;i<6;++i) {
    for (int j=1;j<2;++j) {
      float effBin = (data_ctrtag_2j_num_h->GetBinContent(i,j)-other_ctrtag_2j_num_h->GetBinContent(i,j))/(data_ctrtag_2j_den_h->GetBinContent(i,j)-other_ctrtag_2j_den_h->GetBinContent(i,j));
      float effErrBin = efficiencyErr(effBin, (data_ctrtag_2j_den_h->GetBinContent(i,j)-other_ctrtag_2j_den_h->GetBinContent(i,j)) );
      if (effBin>0) {
	yield += data_ctrtag_2j_h->GetBinContent(i,j)*(1.-effBin)/effBin;
	error += pow(data_ctrtag_2j_h->GetBinContent(i,j)*effErrBin/pow(effBin,2),2) + data_ctrtag_2j_h->GetBinContent(i,j)*pow((1.-effBin)/effBin,2);
      }
    }
  }
  error = sqrt(error);  
  float sf2j = yield/sigreg_top_2j.first;
  float sf2jerr = error/sigreg_top_2j.first;//fixme?
  float k2j = 1.+error/sigreg_top_2j.first;
  pair<float, float> topData2j = make_pair<float, float>(yield,error);
  if (printAll) {
    //MC expected
    cout << "MC after veto: " << sigreg_top_2j.first << " +/- " << sigreg_top_2j.second << endl;
    //MC efficiency calculation
    TH2F *ttbar_ctrtag_2j_num_h = new TH2F("ttbar_ctrtag_2j_num_h","ttbar_ctrtag_2j_num_h",5,0,2.5,1,0,200);
    fillPlot("ctrjetetapt",ttbar_ctrtag_2j_num_h, main_dir+topww_dir+"ttbar_powheg", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=bTagCtr=", lumi, useJson, applyEff, doFake, doPUw);
    TH2F *ttbar_ctrtag_2j_den_h = new TH2F("ttbar_ctrtag_2j_den_h","ttbar_ctrtag_2j_den_h",5,0,2.5,1,0,200);
    fillPlot("ctrjetetapt",ttbar_ctrtag_2j_den_h, main_dir+topww_dir+"ttbar_powheg", wwSelNoMetNoTV, TopTagNotInJets, mass, 2, anaRegion+"=nobJet3=bVetoFwd=", lumi, useJson, applyEff, doFake, doPUw);
    TH2F *top_ctrtag_2j_num_h = new TH2F("top_ctrtag_2j_num_h","top_ctrtag_2j_num_h",5,0,2.5,1,0,200);
    top_ctrtag_2j_num_h->Add(ttbar_ctrtag_2j_num_h);
    TH2F *top_ctrtag_2j_den_h = new TH2F("top_ctrtag_2j_den_h","top_ctrtag_2j_den_h",5,0,2.5,1,0,200);
    top_ctrtag_2j_den_h->Add(ttbar_ctrtag_2j_den_h);
    for (int i=1;i<6;++i) {
      for (int j=1;j<2;++j) {
	float effBin = top_ctrtag_2j_den_h->GetBinContent(i,j)>0 ? top_ctrtag_2j_num_h->GetBinContent(i,j)/top_ctrtag_2j_den_h->GetBinContent(i,j) : 0;
	cout << i << "," << j << " MC: " << top_ctrtag_2j_num_h->GetBinContent(i,j) << " " << top_ctrtag_2j_den_h->GetBinContent(i,j) << " " << effBin << endl;
      }
    }
    //background
    for (int i=1;i<6;++i) {
      for (int j=1;j<2;++j) {
	float effBin = other_ctrtag_2j_den_h->GetBinContent(i,j)>0 ? other_ctrtag_2j_num_h->GetBinContent(i,j)/other_ctrtag_2j_den_h->GetBinContent(i,j) : 0;
	cout << i << "," << j << " OT: " << other_ctrtag_2j_num_h->GetBinContent(i,j) << " " << other_ctrtag_2j_den_h->GetBinContent(i,j) << " " << effBin << endl;
      }
    }
    //data
    for (int i=1;i<6;++i) {
      for (int j=1;j<2;++j) {
	float effBin = (data_ctrtag_2j_num_h->GetBinContent(i,j)-other_ctrtag_2j_num_h->GetBinContent(i,j))/(data_ctrtag_2j_den_h->GetBinContent(i,j)-other_ctrtag_2j_den_h->GetBinContent(i,j));
	//float effErrBin = efficiencyErr(effBin, (data_ctrtag_2j_den_h->GetBinContent(i,j)-other_ctrtag_2j_den_h->GetBinContent(i,j)) );
	cout << i << "," << j << " DA: " << data_ctrtag_2j_h->GetBinContent(i,j) << " " << data_ctrtag_2j_num_h->GetBinContent(i,j)-other_ctrtag_2j_num_h->GetBinContent(i,j) 
	     << " " << data_ctrtag_2j_den_h->GetBinContent(i,j)-other_ctrtag_2j_den_h->GetBinContent(i,j)
	     << " " << effBin << " " << data_ctrtag_2j_h->GetBinContent(i,j)*(1.-effBin)/effBin << endl;
      }
    }
    cout << "data estimation: " << yield << " +/- " << error << endl;
    cout << "sf: " << sf2j << endl;
    cout << "k: " << k2j << endl;
  }

  ///////////////////////////////////////// 1-JET BIN /////////////////////////////////////////
  //   ********* this does not work because of soft muons:
  //   pair<float, float> sigreg_ttbar_0j  = getYield(main_dir+topww_dir+"ttbar_powheg", wwSelNoMet, noVeto, mass, 0, "dphijet,minmet40", lumi, useJson, applyEff, doFake);
  //   pair<float, float> sigreg_tw_0j     = getYield(main_dir+topww_dir+"tw",    wwSelNoMet, noVeto, mass, 0, "dphijet,minmet40", lumi, useJson, applyEff, doFake);
  //   ********* so we need to do this to compare with the data prediction:
  pair<float, float> sigreg_ttbar_1j = getYield(main_dir+topww_dir+"ttbar_powheg",    wwSelNoMetNoTV, TopTagNotInJets, mass, 1, anaRegion+"=nobJet1=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_1j    = getYield(main_dir+topww_dir+"tw",       wwSelNoMetNoTV, TopTagNotInJets, mass, 1, anaRegion+"=nobJet1=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_top_1j = make_pair<float, float>(sigreg_ttbar_1j.first+sigreg_tw_1j.first, sqrt(pow(sigreg_ttbar_1j.second,2)+pow(sigreg_tw_1j.second,2)));
  // here we assume SF for tW is 1
  pair<float, float> vetoEff1j = topVetoEffEstimation(mass, 1, lumi, anaRegion, useJson, applyEff, doPUw, 1.);
  pair<float, float> top_tag_data_1j = getYield(main_dir+topww_dir+"data.root", wwSelNoMetNoTV|OneBJet, TopTagNotInJets, mass, 1, anaRegion, 0, useJson, 0, 0, 0);
  pair<float, float> bkg_exp_cr_1j = evaluateBackground(main_dir+topww_dir, wwSelNoMetNoTV|OneBJet, TopTagNotInJets, mass, 1, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  //this can be replaced with a simple count
  pair<float, float> topData1j = topBgEstimation(mass, 1, lumi, anaRegion, vetoEff1j.first, vetoEff1j.second, useJson, applyEff, doFake, doPUw);
  float sf1j = topData1j.first/sigreg_top_1j.first;
  float sf1jerr = sf1j*sqrt(pow(topData1j.second/topData1j.first,2)+pow(sigreg_top_1j.second/sigreg_top_1j.first,2));//fixme: should not include MC uncertainty?
  //float sf1jpercerr = 100*sqrt(pow(topData1j.second/topData1j.first,2)+pow(sigreg_top_1j.second/sigreg_top_1j.first,2));

  ///////////////////////////////////////// 0-JET BIN /////////////////////////////////////////
  pair<float, float> sigreg_ttbar_0j = getYield(main_dir+topww_dir+"ttbar_powheg",    wwSelNoMet, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_0j    = getYield(main_dir+topww_dir+"tw",       wwSelNoMet, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_top_0j = make_pair<float, float>(sigreg_ttbar_0j.first+sigreg_tw_0j.first,sqrt(pow(sigreg_ttbar_0j.second,2)+pow(sigreg_tw_0j.second,2)));
  // here we use the 1-jet bin SF for tW
  pair<float, float> vetoEff0j = topVetoEffEstimation(mass, 0, lumi, anaRegion, useJson, applyEff, doPUw, sf1j);
  pair<float, float> top_tag_data_0j = getYield(main_dir+topww_dir+"data.root", wwSelNoMetNoTV|TopTagNotInJets, noVeto,  mass, 0, anaRegion, 0, useJson, 0, 0, 0);
  pair<float, float> bkg_exp_cr_0j = evaluateBackground(main_dir+topww_dir, wwSelNoMetNoTV|TopTagNotInJets, noVeto,  mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  //this can be replaced with a simple count
  pair<float, float> topData0j = topBgEstimation(mass, 0, lumi, anaRegion, vetoEff0j.first, vetoEff0j.second, useJson, applyEff, doFake, doPUw);
  float sf0j = topData0j.first/sigreg_top_0j.first;
  float sf0jerr = sf0j*sqrt(pow(topData0j.second/topData0j.first,2)+pow(sigreg_top_0j.second/sigreg_top_0j.first,2));//fixme: should not include MC uncertainty?
  //float sf0jpercerr = 100*sqrt(pow(topData0j.second/topData0j.first,2)+pow(sigreg_top_0j.second/sigreg_top_0j.first,2));

  bool doLatex = false;

  if (!doLatex) {
    cout << "--------------------------------------------------------------------------------------------------" << endl;
    cout << Form("| %40s | %-15s | %-15s | %-15s |","Sample","0-jet","1-jet","2-jet") << endl;
    cout << "--------------------------------------------------------------------------------------------------" << endl;
  } else {
    cout << "\\hline" << endl;
    cout << Form(" %40s & %-15s & %-15s & %-15s \\\\","Sample","0-jet","1-jet","2-jet") << endl;
    cout << "\\hline" << endl;
  }

  TString formstr = "| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |";
  if (doLatex) formstr = " %40s & %5.1f $\\pm$ %-5.1f & %5.1f $\\pm$ %-5.1f & %5.1f $\\pm$ %-5.1f \\\\";
  cout << Form(formstr,
	       "Estimated top events in simulation",
	       round(10.*sigreg_top_0j.first)/10.,round(10.*sigreg_top_0j.second)/10.,
	       round(10.*sigreg_top_1j.first)/10.,round(10.*sigreg_top_1j.second)/10., 
	       round(10.*sigreg_top_2j.first)/10.,round(10.*sigreg_top_2j.second)/10.) 
       << endl;
  cout << Form(formstr,
	       "tagging efficiency (\\%)",
	       round(1000.*vetoEff0j.first)/10.,round(1000.*vetoEff0j.second)/10.,
	       round(1000.*vetoEff1j.first)/10.,round(1000.*vetoEff1j.second)/10.,
	       0.,0.) 
       << endl;
  cout << Form(formstr,
	       "top-tagged events in data",
	       round(10.*top_tag_data_0j.first)/10.,round(10.*top_tag_data_0j.second)/10.,
	       round(10.*top_tag_data_1j.first)/10.,round(10.*top_tag_data_1j.second)/10.,
	       0.,0.)
       << endl;
  cout << Form(formstr,
	       "background events in control region",
	       round(10.*bkg_exp_cr_0j.first)/10.,round(10.*bkg_exp_cr_0j.second)/10.,
	       round(10.*bkg_exp_cr_1j.first)/10.,round(10.*bkg_exp_cr_1j.second)/10.,
	       0.,0.) 
       << endl;
  cout << Form(formstr,
	       "Data-driven top background estimate",
	       round(10.*topData0j.first)/10.,round(10.*topData0j.second)/10.,
	       round(10.*topData1j.first)/10.,round(10.*topData1j.second)/10.,
	       round(10.*topData2j.first)/10.,round(10.*topData2j.second)/10.) 
       << endl;
  formstr = "| %40s | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f |";
  if (doLatex) formstr = " %40s & %5.2f $\\pm$ %-4.2f & %5.2f $\\pm$ %-4.2f & %5.2f $\\pm$ %-4.2f \\\\";
  cout << Form(formstr,
	       "Scale factors",
	       round(100.*sf0j)/100.,round(100.*sf0jerr)/100.,
	       round(100.*sf1j)/100.,round(100.*sf1jerr)/100., 
	       round(100.*sf2j)/100.,round(100.*sf2jerr)/100.) 
       << endl;
  if (!doLatex) cout << "--------------------------------------------------------------------------------------------------" << endl;
  else cout << "\\hline" << endl;

  ofstream myfile;
  TString fname = "TopBkgScaleFactors.h";
  myfile.open(fname);
  ostream &out = myfile;

  out << Form("Double_t TopBkgScaleFactor(Int_t jetBin) {\n");
  out << Form("  assert(jetBin >=0 && jetBin <= 2);\n");
  out << Form("  Double_t TopBkgScaleFactor[3] = { %7.5f,%7.5f,%7.5f   };\n",sf0j,sf1j,sf2j);
  out << Form("  return TopBkgScaleFactor[jetBin];\n");
  out << Form("}\n");
  
  out << Form("Double_t TopBkgScaleFactorKappa(Int_t jetBin) {\n");
  out << Form("  assert(jetBin >=0 && jetBin <= 2);\n");
  out << Form("  Double_t TopBkgScaleFactorKappa[3] = { %7.5f,%7.5f,%7.5f   };\n",1.+sf0jerr/sf0j,1.+sf1jerr/sf1j,k2j);
  out << Form("  return TopBkgScaleFactorKappa[jetBin];\n");
  out << Form("}\n");

}

void topBg(float lumi) {
  makeTopTable(lumi);
}
