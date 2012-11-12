#include "common.C"

float getK(TString sample, unsigned int cut, unsigned int veto, int mass, unsigned int njets, float lumiSample, 
	   bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false){
  float zmmNoMet = getYield(sample, cut, veto, mass, njets, "=zregion=mmfs=", lumiSample, useJson, applyEff, doFake, doPUw).first;//FIXME: should account for pt if needed
  float zeeNoMet = getYield(sample, cut, veto, mass, njets, "=zregion=eefs=", lumiSample, useJson, applyEff, doFake, doPUw).first;
  return sqrt(zeeNoMet/zmmNoMet);//the error is negligible
}

pair<float, float> getZYieldInData(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString region, 
				   float lumi, float kee, TString mmee, bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false) {

  bool printAll = 0;
  //sample is assumed to be data
  float lumiSample = 0.;

  //get Z yields after full met
  float zmm = getYield(main_dir+dy_dir+"data", cut, veto, mass, njets, "=zregion=mmfs=dymvaallfs="+region, lumiSample, useJson, false, doFake, false).first;
  float zme = getYield(main_dir+dy_dir+"data", cut, veto, mass, njets, "=zregion=mefs=dymvaallfs="+region, lumiSample, useJson, false, doFake, false).first;
  float zem = getYield(main_dir+dy_dir+"data", cut, veto, mass, njets, "=zregion=emfs=dymvaallfs="+region, lumiSample, useJson, false, doFake, false).first;
  float zee = getYield(main_dir+dy_dir+"data", cut, veto, mass, njets, "=zregion=eefs=dymvaallfs="+region, lumiSample, useJson, false, doFake, false).first;
  float zmmofs = zmm - 0.5*(zme+zem)/kee;
  float zmmofs_err = sqrt( zmm + 0.25*(zme+zem)/pow(kee,2) );
  float zeeofs = zee - 0.5*(zme+zem)*kee;
  float zeeofs_err = sqrt( zee + 0.25*(zme+zem)*pow(kee,2) );
  float zofs = zmmofs+zeeofs;
  float zofs_err = sqrt( zmm + zee + 0.25*(zme+zem)*pow(kee+1./kee,2) );
  pair<float, float> wzmm_p = getYield(main_dir+dy_dir+"wz",    cut, veto, mass, njets, "=zregion=mmfs=dymvaallfs=fromZ="+region, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zzmm_p = getYield(main_dir+dy_dir+"zz", cut, veto, mass, njets, "=zregion=mmfs=dymvaallfs=fromZ="+region, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> wzee_p = getYield(main_dir+dy_dir+"wz",    cut, veto, mass, njets, "=zregion=eefs=dymvaallfs=fromZ="+region, lumi, useJson, applyEff, doFake, doPUw);
  pair<float, float> zzee_p = getYield(main_dir+dy_dir+"zz", cut, veto, mass, njets, "=zregion=eefs=dymvaallfs=fromZ="+region, lumi, useJson, applyEff, doFake, doPUw);
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
  if (mmee=="mm") return make_pair<float, float>(zmmofs_corr,zmmofs_corr_err);
  else if (mmee=="ee") return make_pair<float, float>(zeeofs_corr,zeeofs_corr_err);
  else return make_pair<float, float>(zofs_corr,zofs_corr_err);

}

pair<float, float> dyBkgEstimation(unsigned int cut, unsigned int veto, int mass, unsigned int njets, TString region, float lumi,
				   float kee, TString mmee, bool useJson=false, bool applyEff=false, bool doFake=false, bool doPUw=false)  {

  //warning: region and cut should not contain MET cuts (only min(pmet,pTrackMet)>20.)!

  bool printAll = 0;

  TString sffs = "=sffs=";
  TString offs = "=offs=";
  //kee=0.770;//fixme
  float lumicorr = (1.+kee*kee)/(2.*kee);
  if (mmee=="mm") {
    sffs = "=mmfs=";
    offs = "=offs=";
    lumicorr = 0.5/kee;
  } else if(mmee=="ee") {
    sffs = "=eefs=";
    offs = "=offs=";
    lumicorr = 0.5*kee;
  }  

  float yield = 0;
  float error2 = 0;

  TFile* fzeta = TFile::Open("./zeta_allmasses.root");
  TH1F* zeta = 0;
  if (njets<2) {
    if (doMVA) zeta = (TH1F*) fzeta->Get(Form("zeta_dymva_mass0_%ij",njets)); 
    else zeta = (TH1F*) fzeta->Get(Form("zeta_dymva_mass%i_%ij",mass,njets)); 
  } else {
    if (doMVA) zeta = (TH1F*) fzeta->Get(Form("zeta_cut_mass0_%ij",njets)); 
    else zeta = (TH1F*) fzeta->Get(Form("zeta_cut_mass%i_%ij",mass,njets)); 
  }
  assert(zeta);
  //data SF
  TH1F* lomet_sf = new TH1F("lomet_sf","lomet_sf",zeta->GetNbinsX(),zeta->GetXaxis()->GetXmin(),zeta->GetXaxis()->GetXmax());
  fillPlot("ptll",lomet_sf, main_dir+dy_dir+"data", cut, veto, mass, njets, region+"=dymvazloose="+sffs, 0, useJson, applyEff, doFake, doPUw);
  //data OF
  TH1F* lomet_of = new TH1F("lomet_of","lomet_of",zeta->GetNbinsX(),zeta->GetXaxis()->GetXmin(),zeta->GetXaxis()->GetXmax());
  fillPlot("ptll",lomet_of, main_dir+dy_dir+"data", cut, veto, mass, njets, region+"=dymvazloose="+offs, 0, useJson, applyEff, doFake, doPUw);
  lomet_of->Scale(-1.*lumicorr);
  //MC WZ
  TH1F* lomet_wz = new TH1F("lomet_wz","lomet_wz",zeta->GetNbinsX(),zeta->GetXaxis()->GetXmin(),zeta->GetXaxis()->GetXmax());
  fillPlot("ptll",lomet_wz, main_dir+dy_dir+"wz", cut, veto, mass, njets, region+"=dymvazloose="+sffs, lumi, useJson, applyEff, doFake, doPUw);
  lomet_wz->Scale(-1.);
  //MC ZZ
  TH1F* lomet_zz = new TH1F("lomet_zz","lomet_zz",zeta->GetNbinsX(),zeta->GetXaxis()->GetXmin(),zeta->GetXaxis()->GetXmax());
  fillPlot("ptll",lomet_zz, main_dir+dy_dir+"zz", cut, veto, mass, njets, region+"=dymvazloose="+sffs, lumi, useJson, applyEff, doFake, doPUw);
  lomet_zz->Scale(-1.);
  TH1F* lomet = new TH1F("lomet","lomet",zeta->GetNbinsX(),zeta->GetXaxis()->GetXmin(),zeta->GetXaxis()->GetXmax());
  lomet->Add(lomet_sf);
  lomet->Add(lomet_of);
  lomet->Add(lomet_wz);
  lomet->Add(lomet_zz);
  for (int bin=1;bin<=lomet->GetNbinsX();++bin) {
    yield += lomet->GetBinContent(bin)*zeta->GetBinContent(bin);
    error2 += pow(lomet->GetBinContent(bin)*zeta->GetBinError(bin),2)+pow(lomet->GetBinError(bin)*zeta->GetBinContent(bin),2);
  }
  delete lomet;
  delete lomet_sf;
  delete lomet_of;
  delete lomet_wz;
  delete lomet_zz;

  fzeta->Close();
  //get the estimate
  float dy_est = TMath::Max(yield,float(0.));;
  float dy_est_err = sqrt(error2);//only stat for now
  if (printAll) {}
  return make_pair<float, float>(dy_est,dy_est_err);
}

/*
void makeDYTable(float lumi) {

  bool useJson  = true;
  bool applyEff = true;
  bool doFake   = false;
  bool doPUw    = true;

  TString zregion  = "=leppts=ptll45=mtcut=dphicut=zregion=dpjallfs=";
  TString dyregion = "=leppts=ptll45=mtcut=dphicut=masscut=zvetoall=dpjallfs=";

  //int jetbins[] = {0};
  //int jetbins[] = {0,1};
  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  //int masses[] = {0};
  //int masses[] = {0,120};
  int masses[] = {0,120,140,160,180,200};
  //int masses[] = {0,115,120,130,140,150,160,170,180,190,200,250,300};
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
      cout << "---------------------------------------------------------------- " << njets << "-jet bin ---------------------------------------------------------------" << endl;
      cout << 
	Form("| %8s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s |","mass","Z-peak MC","Z-peak Exp.","Z-peak Pred.","Z-peak Bias","Sig. Reg. MC","Sig. Reg. Pred.","Scale Factor")
	   << endl;
    } else {
      cout << "\\hline" << endl;
      cout << Form("\\multicolumn{6}{c}{%i-jet bin} \\\\",njets) << endl;
      cout << "\\hline" << endl;
      cout << Form(" %10s & %-16s & %-15s & %-15s & %-15s & %-14s  \\\\","mass","$N_{in}$(data)","$R_{out/in}$","$N_{out}$(data)","$N_{out}$(MC)","SF(Data/MC)") << endl;
    }

    //compute k
    float zmmNoMet = getYield(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, "=zregion=mmfs=", 0., useJson, false, doFake, false).first;//FIXME: should account for pt if needed
    float zeeNoMet = getYield(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, "=zregion=eefs=", 0., useJson, false, doFake, false).first;
    float kee = sqrt(zeeNoMet/zmmNoMet);//the error is negligible
    //cout << "kee: " << kee << endl;      

    for (int jj=0;jj<nmasses;++jj) {

      float sf = 1.0;
      float sf_err = 1.0;

      int mass = masses[jj];
      if (njets==2 && mass>0) continue;

      if (njets!=2) {
	//cout << "mass: " << mass << endl;      
	pair<float, float> zmmMC   = getYield(main_dir+dy_dir+"dyll",  wwSelNoZVNoMet, noVeto, mass, njets, "=mmfs=dymvacut="+zregion, lumi, false, applyEff, doFake, doPUw);
	pair<float, float> zeeMC   = getYield(main_dir+dy_dir+"dyll",  wwSelNoZVNoMet, noVeto, mass, njets, "=eefs=dymvacut="+zregion, lumi, false, applyEff, doFake, doPUw);	
	zregion  += "evenevt,";	
	pair<float, float> zmmExp = getZYieldInData(wwSelNoZVNoMet, noVeto, mass, njets, zregion, lumi, kee, "mm", useJson, applyEff, doFake, doPUw);
	pair<float, float> zeeExp = getZYieldInData(wwSelNoZVNoMet, noVeto, mass, njets, zregion, lumi, kee, "ee", useJson, applyEff, doFake, doPUw);	
	pair<float, float> zmmPred = dyBkgEstimation(wwSelNoZVNoMet, noVeto, mass, njets, zregion, lumi, kee, "mm", useJson, applyEff, doFake, doPUw);
	pair<float, float> zeePred = dyBkgEstimation(wwSelNoZVNoMet, noVeto, mass, njets, zregion, lumi, kee, "ee", useJson, applyEff, doFake, doPUw);	
	pair<float, float> biasMM = make_pair<float, float>(zmmPred.first/zmmExp.first-1., ratioPoissErr(zmmPred.first, zmmPred.second, zmmExp.first, zmmExp.second));
	pair<float, float> biasEE = make_pair<float, float>(zeePred.first/zeeExp.first-1., ratioPoissErr(zeePred.first, zeePred.second, zeeExp.first, zeeExp.second));	
	pair<float, float> dymmMC   = getYield(main_dir+dy_dir+"dyll",  wwSelNoZVNoMet, noVeto, mass, njets, "=mmfs=dymvacut="+dyregion, lumi, false, applyEff, doFake, doPUw);
	pair<float, float> dyeeMC   = getYield(main_dir+dy_dir+"dyll",  wwSelNoZVNoMet, noVeto, mass, njets, "=eefs=dymvacut="+dyregion, lumi, false, applyEff, doFake, doPUw);
	pair<float, float> dymmPred = dyBkgEstimation(wwSelNoZVNoMet, noVeto, mass, njets, dyregion, lumi, kee, "mm", useJson, applyEff, doFake, doPUw);
	pair<float, float> dyeePred = dyBkgEstimation(wwSelNoZVNoMet, noVeto, mass, njets, dyregion, lumi, kee, "ee", useJson, applyEff, doFake, doPUw);	
	//Z yields are divided by 2 because of even/odd
	cout << Form("| %3i - mm | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.2f +/- %-5.2f | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.2f +/- %-5.2f |",
		     mass,zmmMC.first/2.,zmmMC.second/2.,zmmExp.first,zmmExp.second,zmmPred.first,zmmPred.second,biasMM.first,biasMM.second,
		     dymmMC.first,dymmMC.second,dymmPred.first,dymmPred.second,dymmPred.first/dymmMC.first,dymmPred.second/dymmMC.first) << endl;	
	cout << Form("| %3i - ee | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.2f +/- %-5.2f | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.2f +/- %-5.2f |",
		     mass,zeeMC.first/2.,zeeMC.second/2.,zeeExp.first,zeeExp.second,zeePred.first,zeePred.second,biasEE.first,biasEE.second,
		     dyeeMC.first,dyeeMC.second,dyeePred.first,dyeePred.second,dyeePred.first/dyeeMC.first,dyeePred.second/dyeeMC.first) << endl;	
	if ((dymmMC.first+dyeeMC.first)>0. && (dymmPred.first+dyeePred.first)>0. && dymmPred.second>0. && dyeePred.second>0.) {
	  sf = (dymmPred.first+dyeePred.first)/(dymmMC.first+dyeeMC.first);
	  //mc uncertainty already included in the cards
	  sf_err = sqrt(pow(dymmPred.second,2)+pow(dyeePred.second,2))/(dymmMC.first+dyeeMC.first);
	}
      }
      //cout << sf << " +/- " << sf_err << endl;

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

//       if (mass==0) {
// 	TString formstr = "| %10s | %6.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f |";
// 	if (doLatex) formstr = " %10s & %6.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f \\\\";
// 	cout << Form(formstr,
// 		     "WW",
// 		     round(100.*z.first)/100.,round(100.*z.second)/100.,
// 		     round(100.*r.first)/100.,round(100.*r.second)/100.,
// 		     round(100.*dyData.first)/100.,round(100.*dyData.second)/100.,
// 		     round(100.*(dymmMC.first+dyeeMC.first))/100.,round(100.*sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2)))/100.,
// 		     //round(100.*sf)/100.,round(10.*sf_percerr)/10.)
// 		     round(100.*sf)/100.,round(100.*sf_err)/100.)
// 	     << endl;
//       } else {
// 	TString formstr = "| %6i GeV | %6.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f |";
// 	if (doLatex) formstr = " %6i \\GeVcc & %6.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f & %5.2f $\\pm$ %-5.2f \\\\";
// 	cout << Form(formstr,
// 		     mass,
// 		     round(100.*z.first)/100.,round(100.*z.second)/100.,
// 		     round(100.*r.first)/100.,round(100.*r.second)/100.,
// 		     round(100.*dyData.first)/100.,round(100.*dyData.second)/100.,
// 		     round(100.*(dymmMC.first+dyeeMC.first))/100.,round(100.*sqrt(pow(dymmMC.second,2)+pow(dyeeMC.second,2)))/100.,
// 		     //round(100.*sf)/100.,round(10.*sf_percerr)/10.)
// 		     round(100.*sf)/100.,round(100.*sf_err)/100.)
// 	     << endl;
//       }

    }
    if (!doLatex) cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
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

}
*/

void makeDYTableFast(float lumi) {

  bool useJson  = true;
  bool applyEff = true;
  bool doFake   = false;
  bool doPUw    = true;


  TString dyregion = "=leppts=ptll45=mtcut=dphicut=masscut=zvetoall=dpjallfs=mt80=";

  //int jetbins[] = {0};
  //int jetbins[] = {0,1};
  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  //int masses[] = {0};
  //int masses[] = {140,160};
  int masses[] = {0,120,140,160,180,200};
  //int masses[] = {0,115,120,125,130,140,145,150,160,170,180,190,200,250,300};
  int nmasses = sizeof(masses)/sizeof(int);

  vector<float> vsf0j;
  vector<float> vk0j;
  vector<float> vsf1j;
  vector<float> vk1j;
  vector<float> vsf2j;
  vector<float> vk2j;

  for (int j=0;j<njetbins;++j) {

    int njets = jetbins[j];
    cout << "---------------------------------------------------------------- " << njets << "-jet bin --------------------------------------------------------------" << endl;
    cout << 
      Form("| %7s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s |","mass","MC - MM","Data - MM","MC - EE","Data - EE","MM/EE","Data - LL","Scale Factor")
	 << endl;

    //compute k
    float zmmNoMet = getYield(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, "=zregion=mmfs=", 0., useJson, false, doFake, false).first;//FIXME: should account for pt if needed
    float zeeNoMet = getYield(main_dir+dy_dir+"data", wwSelNoZVNoMet, noVeto, 0, njets, "=zregion=eefs=", 0., useJson, false, doFake, false).first;
    float kee = sqrt(zeeNoMet/zmmNoMet);//the error is negligible
    //cout << "kee: " << kee << endl;      

    for (int jj=0;jj<nmasses;++jj) {

      float sf = 1.0;
      float sf_err = 1.0;

      int mass = masses[jj];
      if (njets==2 && doMVA) continue;
      doVBF=0;
      if (njets==2) {
	doVBF=1;
	dyregion.ReplaceAll("mt80","mt30");
	if (mass>100) {
	  dyregion.ReplaceAll("looseVBF=","");
	} else {
	  dyregion+="=looseVBF=";
	}
      }

      pair<float, float> dymmMC   = getYield(main_dir+dy_dir+"dyll",  wwSelNoZVNoMet, noVeto, mass, njets, "=mmfs=dymvacut="+dyregion, lumi, false, applyEff, doFake, doPUw);
      pair<float, float> dyeeMC   = getYield(main_dir+dy_dir+"dyll",  wwSelNoZVNoMet, noVeto, mass, njets, "=eefs=dymvacut="+dyregion, lumi, false, applyEff, doFake, doPUw);
      pair<float, float> dymmPred = dyBkgEstimation(wwSelNoZVNoMet, noVeto, mass, njets, dyregion, lumi, kee, "mm", useJson, applyEff, doFake, doPUw);
      pair<float, float> dyeePred = dyBkgEstimation(wwSelNoZVNoMet, noVeto, mass, njets, dyregion, lumi, kee, "ee", useJson, applyEff, doFake, doPUw);	
      if ((dymmMC.first+dyeeMC.first)>0. && (dymmPred.first+dyeePred.first)>0. && dymmPred.second>0. && dyeePred.second>0.) {
	sf = (dymmPred.first+dyeePred.first)/(dymmMC.first+dyeeMC.first);
	//mc uncertainty already included in the cards
	sf_err = sqrt(pow(dymmPred.second,2)+pow(dyeePred.second,2))/(dymmMC.first+dyeeMC.first);
      }
      cout << Form("| %3i GeV | %5.2f +/- %-5.2f | %5.1f +/- %-5.1f | %5.2f +/- %-5.2f | %5.1f +/- %-5.1f | %5.2f +/- %-5.2f |  %5.1f +/- %-5.1f |%5.2f +/- %-5.2f |",
		   mass,dymmMC.first,dymmMC.second,dymmPred.first,dymmPred.second,dyeeMC.first,dyeeMC.second,dyeePred.first,dyeePred.second,
		   dymmPred.first/dyeePred.first,ratioPoissErr(dymmPred.first, dymmPred.second, dyeePred.first, dyeePred.second),
		   dymmPred.first+dyeePred.first,sqrt(pow(dymmPred.second,2)+pow(dyeePred.second,2)+pow(0.4*dymmPred.first+0.4*dyeePred.first,2)),sf,sf_err) << endl;
      
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
    }

    cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;

  }

  /*
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
  */

}

void dyBg_zeta(float lumi) {
  makeDYTableFast(lumi);
//   doMVA=1;
//   makeDYTableFast(lumi);
//   doMVA=0;
}
