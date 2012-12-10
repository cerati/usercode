#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <sstream>
double string_to_double( const std::string& s ) {
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return 0;
  return x;
} 

float devstd(vector<float> vec) {
  float avg = 0;
  for (unsigned int i=0;i<vec.size();++i) avg+=vec[i];
  avg/=float(vec.size());
  float dev = 0;
  for (unsigned int i=0;i<vec.size();++i) dev+=pow(vec[i]-avg,2);
  dev=sqrt(dev/float(vec.size()-1.));
  return dev;
}

void averageResultsNew(int mass,TString mode,TString inj,TString dir) {

  int ng = 0;

  vector<float> med_obs;
  vector<float> med_exp;
  vector<float> med_l1s;
  vector<float> med_h1s;
  vector<float> med_l2s;
  vector<float> med_h2s;
  vector<float> med_stro;
  vector<float> med_strl;
  vector<float> med_strh;
  vector<float> med_sigo;
  vector<float> med_sige;
  ifstream indump(Form("%s/logs/%i/log_%s_%s_%i.log",dir.Data(),mass,mode.Data(),inj.Data(),mass));
  string line;
  if (indump.is_open()) {
    while ( indump.good() ) {
      getline (indump,line);
      //cout << line << endl;
      size_t found=line.find(Form("limit: %i ",mass));
      if (found!=string::npos) {
	TString myline(line);
	myline.ReplaceAll('[',"");
	myline.ReplaceAll(']',"");
	myline.ReplaceAll(',',' ');
	myline.ReplaceAll("---",' ');
	myline.ReplaceAll(Form("limit: %i ",mass),' ');
	myline.ReplaceAll("strength:",' ');
	myline.ReplaceAll("significance:",' ');
	//cout << myline << endl;
	if ((*myline.Tokenize(' ')).GetEntries()<11) continue;
	float obs = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[0])->GetString().Data() );
	float exp = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[1])->GetString().Data() );
	float l1s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[2])->GetString().Data() );
	float h1s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[3])->GetString().Data() );
	float l2s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[4])->GetString().Data() );
	float h2s = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[5])->GetString().Data() );
	float stro = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[6])->GetString().Data() );
	float strl = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[7])->GetString().Data() );
	float strh = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[8])->GetString().Data() );
	float sigo = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[9])->GetString().Data() );
	float sige = string_to_double( ((TObjString*) (*myline.Tokenize(' '))[10])->GetString().Data() );
	//if ( exp<10. && obs<20. && l1s<50.  && l2s<50. && h1s<50.  && h2s<50.1) {
	  ng++;
	  med_obs.push_back(obs);
	  med_exp.push_back(exp);
	  med_l1s.push_back(l1s);
	  med_h1s.push_back(h1s);
	  med_l2s.push_back(l2s);
	  med_h2s.push_back(h2s);
	  med_stro.push_back(stro);
	  med_strl.push_back(strl);
	  med_strh.push_back(strh);
	  med_sigo.push_back(sigo);
	  med_sige.push_back(sige);
	  //}
      }
    }
    indump.close();
  }

  vector<int> idx_lim;
  vector<int> idx_str;
  vector<int> idx_sig;
  if (ng==0) return;
  for (int i=0;i<ng;++i) {
    idx_lim.push_back(0);
    idx_str.push_back(0);
    idx_sig.push_back(0);
  }
  const float* limv = &med_obs[0];
  TMath::Sort(ng,limv,&idx_lim[0],0);
  const float* strv = &med_stro[0];
  TMath::Sort(ng,strv,&idx_str[0],0);
  const float* sigv = &med_sigo[0];
  TMath::Sort(ng,sigv,&idx_sig[0],0);

  cout << Form("%i  %4.2f - %4.2f + %4.2f  %4.2f  [%4.2f,  %4.2f]  [%4.2f,  %4.2f] --- %4.2f - %4.2f + %4.2f  %4.2f  %4.2f  ---  %4.2f - %4.2f + %4.2f  %4.2f ",
	       mass,
	       TMath::Median(ng,&med_obs[0]), (med_obs[idx_lim[int(0.5*ng)]]-med_obs[idx_lim[int(0.16*ng)]]), (med_obs[idx_lim[int(0.84*ng)]]-med_obs[idx_lim[int(0.5*ng)]]),
	       TMath::Median(ng,&med_exp[0]),
	       TMath::Median(ng,&med_l1s[0]),
	       TMath::Median(ng,&med_h1s[0]),
	       TMath::Median(ng,&med_l2s[0]),
	       TMath::Median(ng,&med_h2s[0]),
	       TMath::Median(ng,&med_stro[0]), (med_stro[idx_str[int(0.5*ng)]]-med_stro[idx_str[int(0.16*ng)]]), (med_stro[idx_str[int(0.84*ng)]]-med_stro[idx_str[int(0.5*ng)]]),
	       TMath::Median(ng,&med_strl[0]),
	       TMath::Median(ng,&med_strh[0]),
	       TMath::Median(ng,&med_sigo[0]), (med_sigo[idx_sig[int(0.5*ng)]]-med_sigo[idx_sig[int(0.16*ng)]]), (med_sigo[idx_sig[int(0.84*ng)]]-med_sigo[idx_sig[int(0.5*ng)]]),
	       TMath::Median(ng,&med_sige[0])
	       ) 
       << endl;

}



void averageResultsNew(TString mode,TString inj,TString dir) {

  //int masses[] = {200,250,300};
  int masses[] = {110,115,120,125,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  int nmasses = sizeof(masses)/sizeof(int);
  for (int j=0;j<nmasses;++j) {
    int mass = masses[j];
    averageResultsNew(mass,mode,inj,dir);
  }


}

/*
root -b -q averageResultsNew.C+\(\"ALLC\",\"125\",\".\"\)
root -b -q averageResultsNew.C+\(\"ALLC\",\"125\",\"systAndStatResults\"\)
*/
