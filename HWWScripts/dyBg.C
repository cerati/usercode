#include "common.C"
#include "Smurf/Analysis/HWWlvlv/DYRoutinValues.h"

TFile* outRFile;

float getK(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, float lumiSample, 
	   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false){
  float zmmNoMet = getYield(sample, cut, veto, mass, njets, "=zregion=mmfs=", lumiSample, useJson, applyEff, doFake, doPUw).first;//FIXME: should account for pt if needed
  float zeeNoMet = getYield(sample, cut, veto, mass, njets, "=zregion=eefs=", lumiSample, useJson, applyEff, doFake, doPUw).first;
  return sqrt(zeeNoMet/zmmNoMet);//the error is negligible
}

pair<float, float> getDYYieldInData(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString region, 
				    float lumi, float kee, bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false) {
  //region should contain the mass region and the met cut
  bool printAll = 0;
  //sample is assumed to be data
  float lumiSample = 0.;
  //get Z yields after full met
  float zmm = getYield(sample, cut, veto, mass, njets, "=mmfs=mtcut="+region, lumiSample, useJson, false, doFake, false).first;//
  float zme = getYield(sample, cut, veto, mass, njets, "=mefs=mtcut="+region, lumiSample, useJson, false, doFake, false).first;//
  float zem = getYield(sample, cut, veto, mass, njets, "=emfs=mtcut="+region, lumiSample, useJson, false, doFake, false).first;//
  float zee = getYield(sample, cut, veto, mass, njets, "=eefs=mtcut="+region, lumiSample, useJson, false, doFake, false).first;//
  float zmmofs = zmm - 0.5*(zme+zem)/kee;
  float zmmofs_err = sqrt( zmm + 0.25*(zme+zem)/pow(kee,2) );
  float zeeofs = zee - 0.5*(zme+zem)*kee;
  float zeeofs_err = sqrt( zee + 0.25*(zme+zem)*pow(kee,2) );
  float zofs = zmmofs+zeeofs;
  float zofs_err = sqrt( zmm + zee + 0.25*(zme+zem)*pow(kee+1./kee,2) );
  pair<float, float> wzmm_p = getYield(main_dir+dy_dir+"wz",    cut, veto, mass, njets, "=mmfs=mtcut="+region, lumi, useJson, applyEff, doFake, doPUw);//=fromZ
  pair<float, float> zzmm_p = getYield(main_dir+dy_dir+"zz_py", cut, veto, mass, njets, "=mmfs=mtcut="+region, lumi, useJson, applyEff, doFake, doPUw);//=fromZ
  pair<float, float> wzee_p = getYield(main_dir+dy_dir+"wz",    cut, veto, mass, njets, "=eefs=mtcut="+region, lumi, useJson, applyEff, doFake, doPUw);//=fromZ
  pair<float, float> zzee_p = getYield(main_dir+dy_dir+"zz_py", cut, veto, mass, njets, "=eefs=mtcut="+region, lumi, useJson, applyEff, doFake, doPUw);//=fromZ
  float wzmm = wzmm_p.first;
  float wzmm_stat_err = wzmm_p.second;
  float wzmm_syst_err = 0.1*(wzmm);//assume 10% syst
  float wzmm_err = sqrt( pow(wzmm_stat_err,2) + pow(wzmm_syst_err,2) );
  float zzmm = zzmm_p.first;
  float zzmm_stat_err = zzmm_p.second;
  float zzmm_syst_err = 0.1*(zzmm);//assume 10% syst
  float zzmm_err = sqrt( pow(zzmm_stat_err,2) + pow(zzmm_syst_err,2) );
  float wzee = wzee_p.first;
  float wzee_stat_err = wzee_p.second;
  float wzee_syst_err = 0.1*(wzee);//assume 10% syst
  float wzee_err = sqrt( pow(wzee_stat_err,2) + pow(wzee_syst_err,2) );
  float zzee = zzee_p.first;
  float zzee_stat_err = zzee_p.second;
  float zzee_syst_err = 0.1*(zzee);//assume 10% syst
  float zzee_err = sqrt( pow(zzee_stat_err,2) + pow(zzee_syst_err,2) );
  float vzmm_syst_err = 0.1*(wzmm+zzmm);//assume 10% syst
  float vzee_syst_err = 0.1*(wzee+zzee);//assume 10% syst
  float vz_syst_err = 0.1*(wzmm+zzmm+wzee+zzee);//assume 10% syst
  float zmmofs_corr = zmmofs-wzmm-zzmm;
  float zmmofs_corr_err = sqrt( pow(zmmofs_err,2) + pow(wzmm_stat_err,2) + pow(zzmm_stat_err,2) + pow(vzmm_syst_err,2) );
  float zeeofs_corr = zeeofs-wzee-zzee;
  float zeeofs_corr_err = sqrt( pow(zeeofs_err,2) + pow(wzee_stat_err,2) + pow(zzee_stat_err,2) + pow(vzee_syst_err,2) );
  float zofs_corr = zofs-wzmm-zzmm-wzee-zzee;
  float zofs_corr_err = sqrt( pow(zofs_err,2) + pow(wzmm_stat_err,2) + pow(zzmm_stat_err,2) + pow(wzee_stat_err,2) + pow(zzee_stat_err,2) + pow(vz_syst_err,2) );
  if (printAll){
    cout << "k: " << kee << endl;
    cout << "Z(after full met): " << zmm << " " << zme  << " " << zem  << " " << zee << endl;
    cout << "OF corr mm, ee, all: " << zmmofs << "+/-" << zmmofs_err << " " << zeeofs << "+/-" << zeeofs_err << " " << zofs << "+/-" << zofs_err << endl;
    cout << "WZ mm, ee: " << wzmm << "+/-" << wzmm_err << " " << wzee << "+/-" << wzee_err << endl;
    cout << "ZZ mm, ee: " << zzmm << "+/-" << zzmm_err << " " << zzee << "+/-" << zzee_err << endl;
    cout << "OF+VV corr mm, ee, all: " << zmmofs_corr << "+/-" << zmmofs_corr_err << " " << zeeofs_corr << "+/-" << zeeofs_corr_err << " " << zofs_corr << "+/-" << zofs_corr_err << endl;
  }
  return make_pair<float, float>(zofs_corr,zofs_corr_err);
}

pair<float, float> computeRoutinData(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, TString regionOut, TString metcut, 
				     float lumi, float kee, bool useJson=0, bool applyEff=false, bool doFake=false, bool doPUw=false)  {
  bool printAll = 0;
  pair<float, float> outYield = getDYYieldInData(main_dir+dy_dir+"data.root", cut, veto, mass, njets, regionOut+metcut, lumi, kee, useJson, applyEff, doFake, doPUw);
  pair<float, float> inzYield = getDYYieldInData(main_dir+dy_dir+"data.root", cut, veto, mass, njets, regionIn +metcut, lumi, kee, useJson, applyEff, doFake, doPUw);
  float r_all = outYield.first/inzYield.first;
  float r_all_err = sqrt( pow(outYield.second,2)/pow(inzYield.first,2) + outYield.first*pow(inzYield.second,2)/pow(inzYield.first,4) );
  if (printAll) {
    cout << "computing Rout/in for metcut: " << metcut << endl;
    cout << "Nin: " << inzYield.first << "+/-" << inzYield.second << endl;
    cout << "Nout: " << outYield.first << "+/-" << outYield.second << endl;
    cout << "R(out/in): " << r_all << "+/-" << r_all_err << endl;
  }
  return make_pair<float, float>(r_all,r_all_err);
}

pair<float, float> computeRoutinDatawithSyst(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, TString regionOut, float lumi, 
					   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false)  {
  bool printAll = 0;
  float kee = getK(main_dir+dy_dir+"data.root", cut, veto, mass, njets, 0, useJson);
  pair<float, float> rbin1 = computeRoutinData(cut, veto, mass, njets, regionIn+"=zregion=", regionOut, "=met2025=", lumi, kee, useJson, applyEff, doFake, doPUw);
  pair<float, float> rbin2 = computeRoutinData(cut, veto, mass, njets, regionIn+"=zregion=", regionOut, "=met2530=", lumi, kee, useJson, applyEff, doFake, doPUw);
  pair<float, float> rbin3 = computeRoutinData(cut, veto, mass, njets, regionIn+"=zregion=", regionOut, "=met3037=", lumi, kee, useJson, applyEff, doFake, doPUw);
  float r_all = rbin3.first;
  float r_all_stat_err = rbin3.second;
  float r_all_syst_err = max(fabs(r_all-rbin2.first),fabs(r_all-rbin1.first));
  float r_all_err = sqrt( pow(r_all_stat_err,2) + pow(r_all_syst_err,2) );

  pair<float, float> rbin4 = computeRoutinData(cut, veto, mass, njets, regionIn+"=zregion=", regionOut, "=met37up=", lumi, kee, useJson);

  outRFile->cd();
  TH1F* rplot = new TH1F(Form("data_mh%i_%ij",mass,njets),Form("data_mh%i_%ij",mass,njets),4,0,4);
  rplot->SetBinContent(1,rbin1.first);rplot->SetBinError(1,rbin1.second);
  rplot->SetBinContent(2,rbin2.first);rplot->SetBinError(2,rbin2.second);
  rplot->SetBinContent(3,rbin3.first);rplot->SetBinError(3,rbin3.second);
  rplot->SetBinContent(4,rbin4.first);rplot->SetBinError(4,rbin4.second);
  rplot->Write();
  delete rplot;

  if (printAll) {
    cout << "r values met bins (data): " << rbin1.first << " " << rbin2.first << " " << rbin3.first << endl;
    cout << "r_all(data): " << r_all << "+/-" << r_all_stat_err << "+/-" << r_all_syst_err << endl;
  }
  return make_pair<float, float>(r_all,r_all_err);
}

pair<float, float> computeRoutinMC(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, TString regionOut, TString metcut, float lumi, 
				   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false)  {
  bool printAll = 0;
  pair<float, float> in_mm = getYield(main_dir+dy_dir+"dymm", wwSelLepOnly|cut, veto, mass, njets, regionIn+"=zregion=mmfs="+metcut, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> in_ee = getYield(main_dir+dy_dir+"dyee", wwSelLepOnly|cut, veto, mass, njets, regionIn+"=zregion=eefs="+metcut, lumi, useJson, applyEff, doFake, doPUw);
  float zmm_in = in_mm.first;
  float zee_in = in_ee.first;
  float zmm_in_err = in_mm.second;
  float zee_in_err = in_ee.second;
  float z_in = zee_in + zmm_in;
  float z_in_err = sqrt( pow(zee_in_err,2) + pow(zmm_in_err,2) );
  pair<float, float> out_mm = getYield(main_dir+dy_dir+"dymm", wwSelLepOnly|ZVeto|cut, veto, mass, njets, regionOut+"=mmfs="+metcut, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> out_ee = getYield(main_dir+dy_dir+"dyee", wwSelLepOnly|ZVeto|cut, veto, mass, njets, regionOut+"=eefs="+metcut, lumi, useJson, applyEff, doFake, doPUw);
  float zmm_out = out_mm.first;
  float zee_out = out_ee.first;
  float zmm_out_err = out_mm.second;
  float zee_out_err = out_ee.second;
  float z_out = zee_out + zmm_out;
  float z_out_err = sqrt( pow(zee_out_err,2) + pow(zmm_out_err,2) );
  float r_all = z_out/z_in;
  float r_all_err = r_all*sqrt( pow(z_out_err/z_out,2) + pow(z_in_err/z_in,2) );
  if (printAll) {
    cout << "computing Rout/in for metcut: " << metcut << endl;
    cout << "Nin mm, ee, all: " << zmm_in << "+/-" << zmm_in_err << " " << zee_in << "+/-" << zee_in_err << " " << z_in << "+/-" << z_in_err << endl;
    cout << "Nout mm, ee, all: " << zmm_out << "+/-" << zmm_out_err << " " << zee_out << "+/-" << zee_out_err << " " << z_out << "+/-" << z_out_err << endl;
    float rmm_all = zmm_out/zmm_in;
    float rmm_all_err = rmm_all*sqrt( pow(zmm_out_err/zmm_out,2) + pow(zmm_in_err/zmm_in,2) );
    float ree_all = zee_out/zee_in;
    float ree_all_err = ree_all*sqrt( pow(zee_out_err/zee_out,2) + pow(zee_in_err/zee_in,2) );
    cout << "R(out/in) mm, ee, all: " << rmm_all << "+/-" << rmm_all_err << " " << ree_all << "+/-" << ree_all_err << " " << r_all << "+/-" << r_all_err << endl;
  }
  return make_pair<float, float>(r_all,r_all_err);
}

pair<float, float> computeRoutinMCwithSyst(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, TString regionOut, float lumi, 
					   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false)  {
  bool printAll = 0;
  pair<float, float> rbin1 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "=met2025=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> rbin2 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "=met2530=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> rbin3 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "=met3037=", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> rbin4 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "=met37up=", lumi, useJson, applyEff, doFake, doPUw);
  if (rbin4.second/rbin4.first>0.40 || !isfinite(rbin4.first)) {
    //do not consider the last bin
    rbin4 = make_pair<float, float>(rbin3.first,rbin3.second);
  }
  float r_all = rbin3.first;
  float r_all_stat_err = rbin3.second;
  float r_all_syst_err = max(fabs(r_all-rbin4.first),max(fabs(r_all-rbin2.first),fabs(r_all-rbin1.first)));
  float r_all_err = sqrt( pow(r_all_stat_err,2) + pow(r_all_syst_err,2) );

  outRFile->cd();
  TH1F* rplot = new TH1F(Form("mc_mh%i_%ij",mass,njets),Form("mc_mh%i_%ij",mass,njets),4,0,4);
  rplot->SetBinContent(1,rbin1.first);rplot->SetBinError(1,rbin1.second);
  rplot->SetBinContent(2,rbin2.first);rplot->SetBinError(2,rbin2.second);
  rplot->SetBinContent(3,rbin3.first);rplot->SetBinError(3,rbin3.second);
  rplot->SetBinContent(4,rbin4.first);rplot->SetBinError(4,rbin4.second);
  rplot->Write();
  delete rplot;

  if (printAll) {
    cout << "r values in met bins (MC): " << rbin1.first << " " << rbin2.first << " " << rbin3.first << " " << rbin4.first << endl;
    cout << "r_all(MC): " << r_all << "+/-" << r_all_stat_err << "+/-" << r_all_syst_err << endl;
  }
  return make_pair<float, float>(r_all,r_all_err);
}

pair<float, float> getZYieldInData(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, 
				   float lumi, bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false) {

  bool printAll = 0;
  //sample is assumed to be data
  float lumiSample = 0.;

  //compute k
  float zmmNoMet = getYield(sample, cut, veto, mass, njets, "=zregion=mmfs=", lumiSample, useJson, false, doFake, false).first;//FIXME: should account for pt if needed
  float zeeNoMet = getYield(sample, cut, veto, mass, njets, "=zregion=eefs=", lumiSample, useJson, false, doFake, false).first;
  float kee = sqrt(zeeNoMet/zmmNoMet);//the error is negligible

  //get Z yields after full met
  float zmm = getYield(sample, cut, veto, mass, njets, "=zregion=mmfs=mmvtxallfs=mtcut="+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zme = getYield(sample, cut, veto, mass, njets, "=zregion=mefs=mmvtxallfs=mtcut="+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zem = getYield(sample, cut, veto, mass, njets, "=zregion=emfs=mmvtxallfs=mtcut="+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zee = getYield(sample, cut, veto, mass, njets, "=zregion=eefs=mmvtxallfs=mtcut="+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zmmofs = zmm - 0.5*(zme+zem)/kee;
  float zmmofs_err = sqrt( zmm + 0.25*(zme+zem)/pow(kee,2) );
  float zeeofs = zee - 0.5*(zme+zem)*kee;
  float zeeofs_err = sqrt( zee + 0.25*(zme+zem)*pow(kee,2) );
  float zofs = zmmofs+zeeofs;
  float zofs_err = sqrt( zmm + zee + 0.25*(zme+zem)*pow(kee+1./kee,2) );
  pair<float, float> wzmm_p = getYield(main_dir+dy_dir+"wz",    cut, veto, mass, njets, "=zregion=mmfs=minmetvtx=mtcut=fromZ="+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zzmm_p = getYield(main_dir+dy_dir+"zz_py", cut, veto, mass, njets, "=zregion=mmfs=minmetvtx=mtcut=fromZ="+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wzee_p = getYield(main_dir+dy_dir+"wz",    cut, veto, mass, njets, "=zregion=eefs=minmetvtx=mtcut=fromZ="+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zzee_p = getYield(main_dir+dy_dir+"zz_py", cut, veto, mass, njets, "=zregion=eefs=minmetvtx=mtcut=fromZ="+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  float wzmm = wzmm_p.first;
  float wzmm_stat_err = wzmm_p.second;
  float wzmm_syst_err = 0.1*(wzmm);//assume 10% syst
  float wzmm_err = sqrt( pow(wzmm_stat_err,2) + pow(wzmm_syst_err,2) );
  float zzmm = zzmm_p.first;
  float zzmm_stat_err = zzmm_p.second;
  float zzmm_syst_err = 0.1*(zzmm);//assume 10% syst
  float zzmm_err = sqrt( pow(zzmm_stat_err,2) + pow(zzmm_syst_err,2) );
  float wzee = wzee_p.first;
  float wzee_stat_err = wzee_p.second;
  float wzee_syst_err = 0.1*(wzee);//assume 10% syst
  float wzee_err = sqrt( pow(wzee_stat_err,2) + pow(wzee_syst_err,2) );
  float zzee = zzee_p.first;
  float zzee_stat_err = zzee_p.second;
  float zzee_syst_err = 0.1*(zzee);//assume 10% syst
  float zzee_err = sqrt( pow(zzee_stat_err,2) + pow(zzee_syst_err,2) );
  float vzmm_syst_err = 0.1*(wzmm+zzmm);//assume 10% syst
  float vzee_syst_err = 0.1*(wzee+zzee);//assume 10% syst
  float vz_syst_err = 0.1*(wzmm+zzmm+wzee+zzee);//assume 10% syst
  float zmmofs_corr = zmmofs-wzmm-zzmm;
  float zmmofs_corr_err = sqrt( pow(zmmofs_err,2) + pow(wzmm_stat_err,2) + pow(zzmm_stat_err,2) + pow(vzmm_syst_err,2) );
  float zeeofs_corr = zeeofs-wzee-zzee;
  float zeeofs_corr_err = sqrt( pow(zeeofs_err,2) + pow(wzee_stat_err,2) + pow(zzee_stat_err,2) + pow(vzee_syst_err,2) );
  float zofs_corr = zofs-wzmm-zzmm-wzee-zzee;
  float zofs_corr_err = sqrt( pow(zofs_err,2) + pow(wzmm_stat_err,2) + pow(zzmm_stat_err,2) + pow(wzee_stat_err,2) + pow(zzee_stat_err,2) + pow(vz_syst_err,2) );
  if (printAll){
    cout << "k: " << kee << endl;
    cout << "Z(after full met): " << zmm << " " << zme  << " " << zem  << " " << zee << endl;
    cout << "OF corr mm, ee, all: " << zmmofs << "+/-" << zmmofs_err << " " << zeeofs << "+/-" << zeeofs_err << " " << zofs << "+/-" << zofs_err << endl;
    cout << "WZ mm, ee: " << wzmm << "+/-" << wzmm_err << " " << wzee << "+/-" << wzee_err << endl;
    cout << "ZZ mm, ee: " << zzmm << "+/-" << zzmm_err << " " << zzee << "+/-" << zzee_err << endl;
    cout << "OF+VV corr mm, ee, all: " << zmmofs_corr << "+/-" << zmmofs_corr_err << " " << zeeofs_corr << "+/-" << zeeofs_corr_err << " " << zofs_corr << "+/-" << zofs_corr_err << endl;
  }
  return make_pair<float, float>(zofs_corr,zofs_corr_err);

}

pair<float, float> dyBkgEstimation(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, TString regionOut, 
				   float lumi, float r_all=0., float r_all_err=0., float zofs_corr=0., float zofs_corr_err=0., 
				   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false)  {

  //warning: cut and veto are common for in and out and 
  //should not contain met, mT nor Z veto cuts
  //(they are properly applied in this method) 
  //also regions should not contain MET cuts
  //the regionIn cut has to be symmetric for all fs

  bool printAll = 0;

  //get yield under z peak
  if (fabs(zofs_corr)<1E-5) {
    pair<float, float> z = getZYieldInData(sample, cut, veto, mass, njets, regionIn,lumi,  useJson, applyEff, doFake, doPUw);
    zofs_corr = z.first;
    zofs_corr_err = z.second;
  }

  //now compute R...
  if (fabs(r_all)<1E-5) {
    pair<float, float> r = computeRoutinMCwithSyst(cut, veto, mass, njets, regionIn, regionOut, lumi, useJson, applyEff, doFake, doPUw);
    r_all     = r.first;
    r_all_err = r.second;
  }

  //get the estimate
  float dy_est = zofs_corr*r_all;
  float dy_est_err = sqrt( pow(zofs_corr,2)*pow(r_all_err,2) + pow(zofs_corr_err,2)*pow(r_all,2) );
  if (printAll) {
    cout << "r_all: " << r_all << "+/-" << r_all_err << endl;
    cout << "OF+VV corr all: " << zofs_corr << "+/-" << zofs_corr_err << endl;
    cout << "data driven DY estimate: " << dy_est << "+/-" << dy_est_err << endl;
    cout << "MC expect.: " << getYield(main_dir+dy_dir+"dymm", ZVeto|cut, veto, mass, njets, regionOut+"=mmfs=minmetvtx=", lumi, useJson, applyEff, doFake, doPUw).first+
                              getYield(main_dir+dy_dir+"dyee", ZVeto|cut, veto, mass, njets, regionOut+"=eefs=minmetvtx=", lumi, useJson, applyEff, doFake, doPUw).first << endl;
  }
  return make_pair<float, float>(dy_est,dy_est_err);
}

void makeDYTable(float lumi) {

  bool useJson  = true;
  bool applyEff = true;
  bool doFake   = false;
  bool doPUw    = true;

  outRFile = TFile::Open("outRFile.root","RECREATE");

  TString regionIn  = "=dpjallfs=leppts=dphicut=ptll45=lep2pt15allfs=";
  TString regionOut = "=dphijet=leppts=dphicut=masscut=ptll45=lep2pt15=mll20=zvetoall=";
  //TString regionIn  = "leppts,dphicut,ptll45";
  //TString regionOut = "leppts,dphicut,masscut,ptll45,mll20,zvetoall";

  //int jetbins[] = {2};
  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  //int masses[] = {0};
  //int masses[] = {0,120,140,160,180,200};
  int masses[] = {0,115,120,130,140,150,160,170,180,190,200,250,300};
  int nmasses = sizeof(masses)/sizeof(int);

  bool doLatex = false;

  vector<float> vsf0j;
  vector<float> vk0j;
  vector<float> vsf1j;
  vector<float> vk1j;
  vector<float> vsf2j;
  vector<float> vk2j;

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];
    if (!doLatex) {
      cout << "----------------------------------------------- " << njets << "-jet bin -----------------------------------------------" << endl;
      cout << Form("| %10s | %-16s | %-15s | %-15s | %-15s | %-14s  |","mass","Nin(data)","R_out/in","Nout(data)","Nout(MC)","SF(Data/MC)") << endl;
    } else {
      cout << "\\hline" << endl;
      cout << Form("\\multicolumn{6}{c}{%i-jet bin} \\\\",njets) << endl;
      cout << "\\hline" << endl;
      cout << Form(" %10s & %-16s & %-15s & %-15s & %-15s & %-14s  \\\\","mass","$N_{in}$(data)","$R_{out/in}$","$N_{out}$(data)","$N_{out}$(MC)","SF(Data/MC)") << endl;
    }

    for (int jj=0;jj<nmasses;++jj) {

      int mass = masses[jj];
      if (njets==2 && mass>0) continue;
      if (njets==2) {
	regionIn+=",looseVBF,";
	regionOut+=",looseVBF,";
	doVBF=1;
      }

      pair<float, float> dymmMC   = getYield(main_dir+dy_dir+"dymm",  wwSelNoMet, noVeto, mass, njets, "=mmfs=minmetvtx=mtcut="+regionOut, lumi, false, applyEff, doFake, doPUw);
      pair<float, float> dyeeMC   = getYield(main_dir+dy_dir+"dyee",  wwSelNoMet, noVeto, mass, njets, "=eefs=minmetvtx=mtcut="+regionOut, lumi, false, applyEff, doFake, doPUw);
      
      pair<float, float> r  = computeRoutinMCwithSyst(wwSelNoZVNoMet, noVeto, mass, njets, regionIn, regionOut, lumi, useJson, applyEff, doFake, doPUw);
      //pair<float, float> rd = computeRoutinDatawithSyst(wwSelNoZVNoMet, noVeto, mass, njets, regionIn, regionOut, lumi, useJson, applyEff, doFake, doPUw);
      //take it from frozen values
      //pair<float, float> r = make_pair<float, float>(RoutinValue(mass,njets),sqrt(pow(RoutinStatError(mass,njets),2)+pow(RoutinSystError(mass,njets),2)));

      pair<float, float> z = getZYieldInData(main_dir+dy_dir+"data.root", wwSelNoZVNoMet, noVeto, mass, njets, regionIn, lumi, useJson, applyEff, doFake, doPUw);
      //avoid negative results putting a higher bound
      if (z.first<0) z = make_pair<float, float>(1.0,1.0);
      
      //Z veto, full met and SF cuts are already applied in the method
      pair<float, float> dyData = dyBkgEstimation(main_dir+dy_dir+"data.root", wwSelNoZVNoMet, noVeto, mass, njets, regionIn, regionOut, lumi, 
						  r.first, r.second, z.first, z.second, useJson, applyEff, doFake, doPUw);

      float sf = 1.0;
      float sf_err = 1.0;
      if ((dymmMC.first+dyeeMC.first)>0.) {
	sf = dyData.first/(dymmMC.first+dyeeMC.first);
	//mc uncertainty already included in the cards
	sf_err = dyData.second/(dymmMC.first+dyeeMC.first);
      }

      if (njets==0) {
	vsf0j.push_back(sf);
	vk0j.push_back(1.+sf_err/sf);
      } else if (njets==1) {
	vsf1j.push_back(sf);
	vk1j.push_back(1.+sf_err/sf);
      } else {
	vsf2j.push_back(sf);
	vk2j.push_back(1.+sf_err/sf);
      }

      if (mass==0) {
	TString formstr = "| %10s | %6.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f |";
	if (doLatex) formstr = " %10s & %6.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f \\\\";
	cout << Form(formstr,
		     "WW",
		     round(100.*z.first)/100.,round(100.*z.second)/100.,
		     round(100.*r.first)/100.,round(100.*r.second)/100.,
		     round(100.*dyData.first)/100.,round(100.*dyData.second)/100.,
		     round(100.*(dymmMC.first+dyeeMC.first))/100.,round(100.*sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2)))/100.,
		     //round(100.*sf)/100.,round(10.*sf_percerr)/10.)
		     round(100.*sf)/100.,round(100.*sf_err)/100.)
	     << endl;
      } else {
	TString formstr = "| %6i GeV | %6.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f |";
	if (doLatex) formstr = " %6i \\GeVcc & %6.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f \\\\";
	cout << Form(formstr,
		     mass,
		     round(100.*z.first)/100.,round(100.*z.second)/100.,
		     round(100.*r.first)/100.,round(100.*r.second)/100.,
		     round(100.*dyData.first)/100.,round(100.*dyData.second)/100.,
		     round(100.*(dymmMC.first+dyeeMC.first))/100.,round(100.*sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2)))/100.,
		     //round(100.*sf)/100.,round(10.*sf_percerr)/10.)
		     round(100.*sf)/100.,round(100.*sf_err)/100.)
	     << endl;
      }
    }
    if (!doLatex) cout << "---------------------------------------------------------------------------------------------------------" << endl;
    else cout << "\\hline" << endl;
  }

  if (nmasses==13) {
    ofstream myfile;
    TString fname = "DYBkgScaleFactors.h";
    myfile.open(fname);
    ostream &out = myfile;

    out << Form("Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {\n");
    out << Form("  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};\n");
    out << Form("  Double_t DYBkgScaleFactorWWPreselection[3] = { %7.5f,%7.5f,%7.5f  };\n",vsf0j[0],vsf1j[0],vsf2j[0]);
    out << Form("  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { \n");
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f },\n",
	   vsf0j[1],vsf0j[2],vsf0j[3],vsf0j[4],vsf0j[5],vsf0j[6],vsf0j[7],vsf0j[8],vsf0j[9],vsf0j[10],vsf0j[11],vsf0j[12]);
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f },\n",
	   vsf1j[1],vsf1j[2],vsf1j[3],vsf1j[4],vsf1j[5],vsf1j[6],vsf1j[7],vsf1j[8],vsf1j[9],vsf1j[10],vsf1j[11],vsf1j[12]);
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f } };\n",
	   vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0],vsf2j[0]);
    out << Form("  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];\n");
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < 12 ; ++m) {\n");
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return DYBkgScaleFactorWWPreselection[jetBin];\n");
    out << Form("  }\n");
    out << Form("}\n");

    out << Form("Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {\n");
    out << Form("  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};\n");
    out << Form("  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { %7.5f,%7.5f,%7.5f  };\n",vk0j[0],vk1j[0],vk2j[0]);
    out << Form("  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { \n");
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f },\n",
	   vk0j[1],vk0j[2],vk0j[3],vk0j[4],vk0j[5],vk0j[6],vk0j[7],vk0j[8],vk0j[9],vk0j[10],vk0j[11],vk0j[12]);
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f },\n",
	   vk1j[1],vk1j[2],vk1j[3],vk1j[4],vk1j[5],vk1j[6],vk1j[7],vk1j[8],vk1j[9],vk1j[10],vk1j[11],vk1j[12]);
    out << Form("    { %7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f } };\n",
	   vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0],vk2j[0]);
    out << Form("  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];\n");
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < 12 ; ++m) {\n");
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return DYBkgScaleFactorWWPreselectionKappa[jetBin];\n");
    out << Form("  }\n");
    out << Form("}\n");
  }

  outRFile->Close();

}

void dyBg(float lumi) {
  makeDYTable(lumi);
}
