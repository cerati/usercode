#include "common.C"

//fixme use of eff SF

pair<float, float> computeRoutinData(TString sample, unsigned int veto, int mass, unsigned int njets, TString region, TString metcut, float kee, bool useJson=0)  {
  cout << "computeRoutinData: WARNING, THIS METHOD IS OUTDATED" << endl;
  bool printAll = 0;
  float zmm_in = getYield(sample, wwSelNoZVNoMet, veto, mass, njets, "zregion,mmfs,dpjallfs,"+metcut, 0, useJson).first;
  float zme_in = getYield(sample, wwSelNoZVNoMet, veto, mass, njets, "zregion,mefs,dpjallfs,"+metcut, 0, useJson).first;
  float zem_in = getYield(sample, wwSelNoZVNoMet, veto, mass, njets, "zregion,emfs,dpjallfs,"+metcut, 0, useJson).first;
  float zee_in = getYield(sample, wwSelNoZVNoMet, veto, mass, njets, "zregion,eefs,dpjallfs,"+metcut, 0, useJson).first;
  if (printAll) cout << "uncorr in: " << zmm_in << " " << zme_in  << " " << zem_in  << " " << zee_in << endl;
  float zmmofs_in = zmm_in - 0.5*(zme_in+zem_in)/kee;
  float zmmofs_in_err = sqrt( zmm_in + 0.25*(zme_in+zem_in)/pow(kee,2) );
  float zeeofs_in = zee_in - 0.5*(zme_in+zem_in)*kee;
  float zeeofs_in_err = sqrt( zee_in + 0.25*(zme_in+zem_in)*pow(kee,2) );
  float zofs_in = zee_in + pow(kee,2)*zmm_in - kee*(zme_in+zem_in);
  float zofs_in_err = sqrt( zee_in + pow(kee,4)*zmm_in + pow(kee,2)*(zme_in+zem_in) );
  if (printAll) cout << "Nin mm, ee, all: " << zmmofs_in << "+/-" << zmmofs_in_err << " " << zeeofs_in << "+/-" << zeeofs_in_err << " " << zofs_in << "+/-" << zofs_in_err << endl;
  float zmm_out = getYield(sample, wwSelNoMet, veto, mass, njets, region+",mmfs,dpjallfs,"+metcut, 0, useJson).first;
  float zme_out = getYield(sample, wwSelNoMet, veto, mass, njets, region+",mefs,dpjallfs,zvetoall,"+metcut, 0, useJson).first;
  float zem_out = getYield(sample, wwSelNoMet, veto, mass, njets, region+",emfs,dpjallfs,zvetoall,"+metcut, 0, useJson).first;
  float zee_out = getYield(sample, wwSelNoMet, veto, mass, njets, region+",eefs,dpjallfs,"+metcut, 0, useJson).first;
  if (printAll) cout << "uncorr out: " << zmm_out << " " << zme_out  << " " << zem_out  << " " << zee_out << endl;
  float zmmofs_out = zmm_out - 0.5*(zme_out+zem_out)/kee;
  float zmmofs_out_err = sqrt( zmm_out + 0.25*(zme_out+zem_out)/pow(kee,2) );
  float zeeofs_out = zee_out - 0.5*(zme_out+zem_out)*kee;
  float zeeofs_out_err = sqrt( zee_out + 0.25*(zme_out+zem_out)*pow(kee,2) );
  float zofs_out = zee_out + pow(kee,2)*zmm_out - kee*(zme_out+zem_out);
  float zofs_out_err = sqrt( zee_out + pow(kee,4)*zmm_out + pow(kee,2)*(zme_out+zem_out) );
  if (printAll) cout << "Nout mm, ee, all: " << zmmofs_out << "+/-" << zmmofs_out_err << " " << zeeofs_out << "+/-" << zeeofs_out_err << " " << zofs_out << "+/-" << zofs_out_err << endl;
  float r_all = zofs_out/zofs_in;
  float r_all_err = sqrt( pow(zofs_out_err,2)/pow(zofs_in,2) + zofs_out*pow(zofs_in_err,2)/pow(zofs_in,4) );
  return make_pair<float, float>(r_all,r_all_err);
}

pair<float, float> computeRoutinMC(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, TString regionOut, TString metcut, float lumi, 
				   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false)  {
  bool printAll = 0;
  pair<float, float> in_mm = getYield(dir_mc+"dymm", wwSelLepOnly|cut, veto, mass, njets, regionIn+",zregion,mmfs,"+metcut, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> in_ee = getYield(dir_mc+"dyee", wwSelLepOnly|cut, veto, mass, njets, regionIn+",zregion,eefs,"+metcut, lumi, useJson, applyEff, doFake, doPUw);
  float zmm_in = in_mm.first;
  float zee_in = in_ee.first;
  float zmm_in_err = in_mm.second;
  float zee_in_err = in_ee.second;
  float z_in = zee_in + zmm_in;
  float z_in_err = sqrt( pow(zee_in_err,2) + pow(zmm_in_err,2) );
  pair<float, float> out_mm = getYield(dir_mc+"dymm", wwSelLepOnly|ZVeto|cut, veto, mass, njets, regionOut+",mmfs,"+metcut, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> out_ee = getYield(dir_mc+"dyee", wwSelLepOnly|ZVeto|cut, veto, mass, njets, regionOut+",eefs,"+metcut, lumi, useJson, applyEff, doFake, doPUw);
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
  pair<float, float> r2025 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "met2025", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> r2530 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "met2530", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> r3037 = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "met3037", lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> r37up = computeRoutinMC(cut, veto, mass, njets, regionIn, regionOut, "met37up", lumi, useJson, applyEff, doFake, doPUw);
  float r_all = r37up.first;
  float r_all_stat_err = r37up.second;
  if (r_all_stat_err/r_all>0.65) {
    r_all = r3037.first;
    r_all_stat_err = r3037.second;
  }
  float r_all_syst_err = max(fabs(r_all-r3037.first),max(fabs(r_all-r2530.first),fabs(r_all-r2025.first)));
  float r_all_err = sqrt( pow(r_all_stat_err,2) + pow(r_all_syst_err,2) );
  if (printAll) {
    cout << "r values in met bins: " << r2025.first << " " << r2530.first << " " << r3037.first << " " << r37up.first << endl;
    cout << "r_all: " << r_all << "+/-" << r_all_stat_err << "+/-" << r_all_syst_err << endl;
  }
  return make_pair<float, float>(r_all,r_all_err);
}

pair<float, float> getZYieldInData(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString regionIn, 
				   float lumi, bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false) {

  bool printAll = 0;
  //sample is assumed to be data
  float lumiSample = 0.;

  //compute k
  float zmmNoMet = getYield(sample, cut, veto, mass, njets, "zregion,mmfs,", lumiSample, useJson, false, doFake, false).first;//FIXME: should account for pt if needed
  float zeeNoMet = getYield(sample, cut, veto, mass, njets, "zregion,eefs,", lumiSample, useJson, false, doFake, false).first;
  float kee = sqrt(zeeNoMet/zmmNoMet);//the error is negligible

  //get Z yields after full met
  float zmm = getYield(sample, FullMET|cut, veto, mass, njets, "zregion,mmfs,mmvtxallfs,mtcut,"+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zme = getYield(sample, FullMET|cut, veto, mass, njets, "zregion,mefs,mmvtxallfs,mtcut,"+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zem = getYield(sample, FullMET|cut, veto, mass, njets, "zregion,emfs,mmvtxallfs,mtcut,"+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zee = getYield(sample, FullMET|cut, veto, mass, njets, "zregion,eefs,mmvtxallfs,mtcut,"+regionIn, lumiSample, useJson, false, doFake, false).first;
  float zmmofs = zmm - 0.5*(zme+zem)/kee;
  float zmmofs_err = sqrt( zmm + 0.25*(zme+zem)/pow(kee,2) );
  float zeeofs = zee - 0.5*(zme+zem)*kee;
  float zeeofs_err = sqrt( zee + 0.25*(zme+zem)*pow(kee,2) );
  float zofs = zmmofs+zeeofs;
  float zofs_err = sqrt( zmm + zee + 0.25*(zme+zem)*pow(kee+1./kee,2) );
  pair<float, float> wzmm_p = getYield(dir_mc+"wz", FullMET|cut, veto, mass, njets, "zregion,mmfs,minmetvtx,mtcut,fromZ,"+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zzmm_p = getYield(dir_mc+"zz", FullMET|cut, veto, mass, njets, "zregion,mmfs,minmetvtx,mtcut,fromZ,"+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wzee_p = getYield(dir_mc+"wz", FullMET|cut, veto, mass, njets, "zregion,eefs,minmetvtx,mtcut,fromZ,"+regionIn, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zzee_p = getYield(dir_mc+"zz", FullMET|cut, veto, mass, njets, "zregion,eefs,minmetvtx,mtcut,fromZ,"+regionIn, lumi, useJson, applyEff, doFake, doPUw);
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
    cout << "MC expect.: " << getYield(dir_mc+"dymm", ZVeto|FullMET|cut, veto, mass, njets, regionOut+",mmfs,minmetvtx", lumi, useJson, applyEff, doFake, doPUw).first+
                              getYield(dir_mc+"dyee", ZVeto|FullMET|cut, veto, mass, njets, regionOut+",eefs,minmetvtx", lumi, useJson, applyEff, doFake, doPUw).first << endl;
  }
  return make_pair<float, float>(dy_est,dy_est_err);
}

void makeDYTable(float lumi) {

  bool useJson  = false;
  bool applyEff = true;
  bool doFake   = false;
  bool doPUw    = true;

  TString regionIn  = "dpjallfs,leppts,dphicut,ptll45,lep2pt15allfs";
  TString regionOut = "dphijet,leppts,dphicut,masscut,ptll45,lep2pt15";

  //int jetbins[] = {0};
  int jetbins[] = {0,1};
  int njetbins = sizeof(jetbins)/sizeof(int);

  //int masses[] = {160};
  int masses[] = {0,120,140,160,180,200};
  int nmasses = sizeof(masses)/sizeof(int);

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];
    cout << "----------------------------------------------- " << njets << "-jet bin -----------------------------------------------" << endl;
    cout << Form("| %10s | %-15s | %-15s | %-15s | %-15s | %-15s  |","mass","Nin(data)","R_out/in","Nout(data)","Nout(MC)","SF(Data/MC)") << endl;

    for (int jj=0;jj<nmasses;++jj) {

      int mass = masses[jj];

      pair<float, float> dymmMC   = getYield(dir_mc+"dymm",  wwSelection, noVeto, mass, njets, "mmfs,minmetvtx,mtcut"+regionOut, lumi, false, applyEff, doFake, doPUw);
      pair<float, float> dyeeMC   = getYield(dir_mc+"dyee",  wwSelection, noVeto, mass, njets, "eefs,minmetvtx,mtcut"+regionOut, lumi, false, applyEff, doFake, doPUw);
      
      pair<float, float> r = computeRoutinMCwithSyst(wwSelNoZVNoMet, noVeto, mass, njets, regionIn, regionOut, lumi, useJson, applyEff, doFake, doPUw);
      
      pair<float, float> z = getZYieldInData(data_file, wwSelNoZVNoMet, noVeto, mass, njets, regionIn, lumi, useJson, applyEff, doFake, doPUw);
      
      //Z veto, full met and SF cuts are already applied in the method
      pair<float, float> dyData = dyBkgEstimation(data_file, wwSelNoZVNoMet, noVeto, mass, njets, regionIn, regionOut, lumi, r.first, r.second, z.first, z.second, useJson, applyEff, doFake, doPUw);

      float sf = dyData.first/(dymmMC.first+dyeeMC.first);
      float sf_err = sf*sqrt( pow(dyData.second/dyData.first,2) + pow(sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2))/(dymmMC.first+dyeeMC.first),2) );
      float sf_percerr = 100.*sf_err/sf;

      if (mass==0) {
	cout << Form("| %10s | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.1f%% |",
		     "WW",
		     round(100.*z.first)/100.,round(100.*z.second)/100.,
		     round(100.*r.first)/100.,round(100.*r.second)/100.,
		     round(100.*dyData.first)/100.,round(100.*dyData.second)/100.,
		     round(100.*(dymmMC.first+dyeeMC.first))/100.,round(100.*sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2)))/100.,
		     round(100.*sf)/100.,round(10.*sf_percerr)/10.)
	     << endl;
      } else {
	cout << Form("| %6i GeV | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.1f%% |",
		     mass,
		     round(100.*z.first)/100.,round(100.*z.second)/100.,
		     round(100.*r.first)/100.,round(100.*r.second)/100.,
		     round(100.*dyData.first)/100.,round(100.*dyData.second)/100.,
		     round(100.*(dymmMC.first+dyeeMC.first))/100.,round(100.*sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2)))/100.,
		     round(100.*sf)/100.,round(10.*sf_percerr)/10.)
	     << endl;
      }
    }
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
  }

}
