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


void removeSyst(TString fname) {

  ifstream indump(fname);
  ofstream outdump("tmp.txt");
  string line;
  bool goon = true;
  if (indump.is_open()) {
    while ( indump.good() && goon ) {
      getline (indump,line);
      //cout << line << endl;
      size_t found=line.find("rate");
      outdump << line << "\n";
      if (found!=string::npos) {
	goon = false;
      }
    }
  }
  gSystem->Exec(TString("mv tmp.txt "+fname).Data());

}

// for file in cards_inj_stat/*/*.txt; do root -b -q removeSyst.C+\(\"${file}\"\); done
