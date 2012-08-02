#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>

void printRValues(bool doMVA=0, bool useMC=1){

  TFile* f = doMVA ? TFile::Open("outRFile_shape.root") : TFile::Open("outRFile_cut.root");

  int jetbins[] = {0,1,2};
  int njetbins = sizeof(jetbins)/sizeof(int);

  int masses[] = {0,115,120,125,130,140,145,150,160,170,180,190,200,250,300};
  int nmasses = sizeof(masses)/sizeof(int);

  vector<float> vr0j;
  vector<float> vst0j;
  vector<float> vsy0j;

  vector<float> vr1j;
  vector<float> vst1j;
  vector<float> vsy1j;

  vector<float> vr2j;
  vector<float> vst2j;
  vector<float> vsy2j;

  for (int j=0;j<njetbins;++j) {
    int njets = jetbins[j];

    cout << endl;
    cout << "njets=" << njets << endl;

    for (int jj=0;jj<nmasses;++jj) {
      int mass = masses[jj];

      TH1F* h = 0;
      h = useMC ? (TH1F*) f->Get(Form("mc_mh%i_%ij",mass,njets)) :  (TH1F*) f->Get(Form("data_mh%i_%ij",mass,njets));

      pair<float, float> rbin1;
      pair<float, float> rbin2;
      pair<float, float> rbin3;
      pair<float, float> rbin4;

      if (h) {
	rbin1 = make_pair<float, float>(h->GetBinContent(1),h->GetBinError(1));
	rbin2 = make_pair<float, float>(h->GetBinContent(2),h->GetBinError(2));
	rbin3 = make_pair<float, float>(h->GetBinContent(3),h->GetBinError(3));
	rbin4 = make_pair<float, float>(h->GetBinContent(4),h->GetBinError(4));
      }

      if (rbin4.second/rbin4.first>0.40 || !isfinite(rbin4.first)) {
	//do not consider the last bin
	rbin4 = make_pair<float, float>(rbin3.first,rbin3.second);
      }
      float r_all = rbin3.first;
      float r_all_stat_err = rbin3.second;
      float r_all_syst_err = max(fabs(r_all-rbin4.first),max(fabs(r_all-rbin2.first),fabs(r_all-rbin1.first)));
      float r_all_err = sqrt( pow(r_all_stat_err,2) + pow(r_all_syst_err,2) );
      
      cout << Form("%3i %4.2f +/- %4.2f +/- %4.2f = %4.2f +/- %4.2f",mass,r_all,r_all_stat_err,r_all_syst_err,r_all,r_all_err) << endl;

      if (njets==0) {
	vr0j.push_back(r_all);
	vst0j.push_back(r_all_stat_err);
	vsy0j.push_back(r_all_syst_err);
      }
      if (njets==1) {
	vr1j.push_back(r_all);
	vst1j.push_back(r_all_stat_err);
	vsy1j.push_back(r_all_syst_err);
      }
      if (njets==2) {
	vr2j.push_back(r_all);
	vst2j.push_back(r_all_stat_err);
	vsy2j.push_back(r_all_syst_err);
      }

    }
  }

  if (nmasses>6 && njetbins==3 && masses[0]==0) {
    ofstream myfile;
    TString fname = "DYRoutinValues.h";
    if (!doMVA) {
      myfile.open(fname);
    } else {
      myfile.open(fname,ios::app);
    }
    ostream &out = myfile;
    
    if (!doMVA) out << Form("Double_t RoutinValue(Int_t mH, Int_t jetBin) {\n");
    else  out << Form("Double_t RoutinValueBDT(Int_t mH, Int_t jetBin) {\n");
    out << Form("  Int_t mHiggs[%i] = {",nmasses-1);
    for (int j=1;j<nmasses-1;++j) out << Form("%i,",masses[j]);
    out << Form("%i};\n",masses[nmasses-1]);
    out << Form("  Double_t RoutinValueWWPreselection[3] = { %7.5f,%7.5f,%7.5f  };\n",vr0j[0],vr1j[0],vr2j[0]);
    out << Form("  Double_t RoutinValueHiggsSelection[3][%i] = { \n",nmasses-1);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vr0j[j]);
    out << Form("%7.5f}, \n",vr0j[nmasses-1]);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vr1j[j]);
    out << Form("%7.5f}, \n",vr1j[nmasses-1]);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vr2j[j]);
    out << Form("%7.5f} }; \n",vr2j[nmasses-1]);
    out << Form("  if(mH == 0) return RoutinValueWWPreselection[jetBin];\n");
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < %i ; ++m) {\n",nmasses);
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return RoutinValueHiggsSelection[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return RoutinValueWWPreselection[jetBin];\n");
    out << Form("  }\n");
    out << Form("}\n");
    
    if (!doMVA) out << Form("Double_t RoutinValueStatError(Int_t mH, Int_t jetBin) {\n");
    else  out << Form("Double_t RoutinValueBDTStatError(Int_t mH, Int_t jetBin) {\n");
    out << Form("  Int_t mHiggs[%i] = {",nmasses-1);
    for (int j=1;j<nmasses-1;++j) out << Form("%i,",masses[j]);
    out << Form("%i};\n",masses[nmasses-1]);
    out << Form("  Double_t RoutinValueWWPreselectionStatError[3] = { %7.5f,%7.5f,%7.5f  };\n",vst0j[0],vst1j[0],vst2j[0]);
    out << Form("  Double_t RoutinValueHiggsSelectionStatError[3][%i] = { \n",nmasses-1);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vst0j[j]);
    out << Form("%7.5f}, \n",vst0j[nmasses-1]);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vst1j[j]);
    out << Form("%7.5f}, \n",vst1j[nmasses-1]);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vst2j[j]);
    out << Form("%7.5f} }; \n",vst2j[nmasses-1]);
    out << Form("  if(mH == 0) return RoutinValueWWPreselectionStatError[jetBin];\n");
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < %i ; ++m) {\n",nmasses);
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return RoutinValueHiggsSelectionStatError[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return RoutinValueWWPreselectionStatError[jetBin];\n");
    out << Form("  }\n");
    out << Form("}\n");
    
    if (!doMVA) out << Form("Double_t RoutinValueSystError(Int_t mH, Int_t jetBin) {\n");
    else  out << Form("Double_t RoutinValueBDTSystError(Int_t mH, Int_t jetBin) {\n");
    out << Form("  Int_t mHiggs[%i] = {",nmasses-1);
    for (int j=1;j<nmasses-1;++j) out << Form("%i,",masses[j]);
    out << Form("%i};\n",masses[nmasses-1]);
    out << Form("  Double_t RoutinValueWWPreselectionSystError[3] = { %7.5f,%7.5f,%7.5f  };\n",vsy0j[0],vsy1j[0],vsy2j[0]);
    out << Form("  Double_t RoutinValueHiggsSelectionSystError[3][%i] = { \n",nmasses-1);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vsy0j[j]);
    out << Form("%7.5f}, \n",vsy0j[nmasses-1]);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vsy1j[j]);
    out << Form("%7.5f}, \n",vsy1j[nmasses-1]);
    out << Form("    { ");
    for (int j=1;j<nmasses-1;++j) out << Form("%7.5f,",vsy2j[j]);
    out << Form("%7.5f} }; \n",vsy2j[nmasses-1]);
    out << Form("  if(mH == 0) return RoutinValueWWPreselectionSystError[jetBin];\n");
    out << Form("  Int_t massIndex = -1;\n");
    out << Form("  for (UInt_t m=0; m < %i ; ++m) {\n",nmasses);
    out << Form("    if (mH == mHiggs[m]) massIndex = m;\n");
    out << Form("  }\n");
    out << Form("  if (massIndex >= 0) {\n");
    out << Form("    return RoutinValueHiggsSelectionSystError[jetBin][massIndex];\n");
    out << Form("  } else {\n");
    out << Form("    return RoutinValueWWPreselectionSystError[jetBin];\n");
    out << Form("  }\n");
    out << Form("}\n");
    
    myfile.close();
  }

}
