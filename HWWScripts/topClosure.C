#include "common.C"
#include "fakeBg.C"

pair<float, float> topVetoEffEstimation(int mass=160, unsigned int njets=0, float lumi=0., TString region = "=dphijet=minmet40=", 
					bool useJson=0, bool applyEff=true, bool doPUw=false, float sfTW=1.) {

  //mass is not used
  bool debug = 1;
  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  TString region_top    = region+"=btagJet1=";//test noSoftMu= 
  TString region_toptag = region+"=btagJet1=";//test noSoftMu=
  if (njets==1) {
    region_top    = region+"=btagJet2=";
    region_toptag = region+"=btag1and2=";
  }

  //get data yields
  float sideband_data_nj_top     = getYield(main_dir+topww_dir+"top-test.root", control_top,    veto, mass, nj_top, region_top,    lumi, false, applyEff, false, doPUw).first;
  float sideband_data_nj_top_tag = getYield(main_dir+topww_dir+"top-test.root", control_toptag, veto, mass, nj_top, region_toptag, lumi, false, applyEff, false, doPUw).first;

  //subtract contamination
  pair<float, float> sb_nj_top_bkg    ;
  pair<float, float> sb_nj_top_tag_bkg;

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
    twlikett = 0*wwtt1j_tw/wwntv1j_tw;
    //float twlikett_err = 0.;//fixme
    float sideband_ttbar  = getYield(main_dir+topww_dir+"ttbar", wwSelNoMetNoTV, veto, mass, njets, region, lumi, false, applyEff, false, doPUw).first;
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
  bool debug = 1;
  unsigned int baseline_toptag=0, control_top=0, control_toptag=0, veto=0, nj_top=0;
  getCutMasks(njets, baseline_toptag, control_top, control_toptag, veto, nj_top);

  TString toptagreg = "";//test =noSoftMu=
  if (njets==1) toptagreg += "=btagJet1=";

  pair<float, float> sb_data_tag = getYield(main_dir+topww_dir+"top-test.root", baseline_toptag, veto, mass, njets,  region+toptagreg, lumi, false, applyEff, false, doPUw);
  float sideband_data_tag  = sb_data_tag.first;

  //compute background contamination
  pair<float, float> sb_tag_bkg ;
  //subtract background
  sideband_data_tag-=sb_tag_bkg.first;
  float sideband_data_tag_err  = sqrt( pow(sb_data_tag.second,2) + pow(sb_tag_bkg.second,2) );

  //get the efficiency (if needed)
  if (eff_veto_data<1E-5) {
    pair<float, float> tagEff = topVetoEffEstimation(mass,njets,lumi,region,useJson,applyEff,doPUw,sfTW);
    eff_veto_data = tagEff.first;
    eff_err_veto_data = tagEff.second;
  }
  //compute the number of top events
  float num_top_data = sideband_data_tag*(1.-eff_veto_data)/eff_veto_data;
  float num_top_err_data = sqrt( pow((1.-eff_veto_data)/eff_veto_data,2)*pow(sideband_data_tag_err,2) + pow(sideband_data_tag,2)*pow(eff_err_veto_data,2)/pow(eff_veto_data,4) );
  if (debug) {
    //float sb_mc_top_tag = getYield(main_dir+topww_dir+"ttbar", baseline_toptag, veto, mass, njets, region+toptagreg, lumi, useJson, applyEff, doFake, doPUw).first+
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
    float sigreg_ttbar  = getYield(main_dir+topww_dir+"ttbar", wwSelNoMet, veto, mass, njets,  region, lumi, useJson, applyEff, doFake, doPUw).first;
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

  TString anaRegion = "=dphijet=dymvacut=ptll45=looseVBF=";//noSoftMu=lep2pt20allfs=

  int mass = 0;

  ///////////////////////////////////////// 2-JET BIN /////////////////////////////////////////
  pair<float, float> sigreg_top_2j;
  float sf2j    = 0;
  float sf2jerr = 0;
  float k2j = 1.;
  pair<float, float> topData2j;

  ///////////////////////////////////////// 1-JET BIN /////////////////////////////////////////
  //   ********* this does not work because of soft muons:
  //   pair<float, float> sigreg_ttbar_0j  = getYield(main_dir+topww_dir+"ttbar", wwSelNoMet, noVeto, mass, 0, "dphijet,minmet40", lumi, useJson, applyEff, doFake);
  //   pair<float, float> sigreg_tw_0j     = getYield(main_dir+topww_dir+"tw",    wwSelNoMet, noVeto, mass, 0, "dphijet,minmet40", lumi, useJson, applyEff, doFake);
  //   ********* so we need to do this to compare with the data prediction:
  pair<float, float> sigreg_ttbar_1j = getYield(main_dir+topww_dir+"ttbar",    wwSelNoMetNoTV, TopTagNotInJets, mass, 1, anaRegion+"=nobJet1=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_1j    = getYield(main_dir+topww_dir+"tw",       wwSelNoMetNoTV, TopTagNotInJets, mass, 1, anaRegion+"=nobJet1=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_top_1j = make_pair<float, float>(sigreg_ttbar_1j.first+sigreg_tw_1j.first, sqrt(pow(sigreg_ttbar_1j.second,2)+pow(sigreg_tw_1j.second,2)));
  // here we assume SF for tW is 1
  pair<float, float> vetoEff1j = topVetoEffEstimation(mass, 1, lumi, anaRegion, useJson, applyEff, doPUw, 1.);
  pair<float, float> top_tag_data_1j = getYield(main_dir+topww_dir+"top-test.root", wwSelNoMetNoTV|OneBJet, TopTagNotInJets, mass, 1, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> bkg_exp_cr_1j;
  //this can be replaced with a simple count
  pair<float, float> topData1j = topBgEstimation(mass, 1, lumi, anaRegion, vetoEff1j.first, vetoEff1j.second, useJson, applyEff, doFake, doPUw);
  float sf1j = topData1j.first/sigreg_top_1j.first;
  float sf1jerr = sf1j*sqrt(pow(topData1j.second/topData1j.first,2)+pow(sigreg_top_1j.second/sigreg_top_1j.first,2));//fixme: should not include MC uncertainty?
  //float sf1jpercerr = 100*sqrt(pow(topData1j.second/topData1j.first,2)+pow(sigreg_top_1j.second/sigreg_top_1j.first,2));

  ///////////////////////////////////////// 0-JET BIN /////////////////////////////////////////
  pair<float, float> sigreg_ttbar_0j = getYield(main_dir+topww_dir+"ttbar",    wwSelNoMet, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_tw_0j    = getYield(main_dir+topww_dir+"tw",       wwSelNoMet, noVeto, mass, 0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> sigreg_top_0j = make_pair<float, float>(sigreg_ttbar_0j.first+sigreg_tw_0j.first,sqrt(pow(sigreg_ttbar_0j.second,2)+pow(sigreg_tw_0j.second,2)));
  // here we use the 1-jet bin SF for tW
  pair<float, float> vetoEff0j = topVetoEffEstimation(mass, 0, lumi, anaRegion, useJson, applyEff, doPUw, sf1j);
  pair<float, float> top_tag_data_0j = getYield(main_dir+topww_dir+"top-test.root", wwSelNoMetNoTV|TopTagNotInJets, noVeto,  mass,  0, anaRegion, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> bkg_exp_cr_0j;
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

void topClosure(float lumi) {
  makeTopTable(lumi);
}
