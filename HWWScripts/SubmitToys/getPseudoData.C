#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>

using namespace std;

// ----------------------------
// First specify these inputs!
// ----------------------------

// directory of the cards where you have the original cards 
// make sure there is no hard coded locations of the input root files since we will use combine
TString dir_cards = "/afs/cern.ch/user/c/cerati/scratch0/Moriond_Injection/cards_def/";
// location of the histograms from the toy creation 
// from this command done with lands
// lands.exe -d card.txt -M Hybrid  -m 125 --minuitSTRATEGY 0 --bWriteToys 1 -n Name --nToysForCLsb 1000 --nToysForCLb 1 --singlePoint 1 --seed 12344  -rMin 0 -rMax 5
  
// ----------------------------

void getPseudoData(TString dir_toys, int mass, TString cms, TString mode, int njets, TString fs, TString ana,int min=0, int max=1) {

  bool noinjection = false;
  TString injmh = "125";

  TString path = "testcards/cards_"+injmh;
  if (noinjection) path = "testcards/cards_def";

  for (int i=min;i<=max;i++) {

    gSystem->Exec(Form("mkdir -p %s_N%i/%i",path.Data(),i,mass));    

    TFile* outf = 0;
    TFile* inf = 0;
    if (mode=="shape") {
      //first copy the input histo except histo_Data
      outf = TFile::Open(Form("%s_N%i/%i/%s%s_%ij.input_%s.root",path.Data(),i,mass,ana.Data(),fs.Data(),njets,cms.Data()),"RECREATE");
      inf = TFile::Open(Form("%s/%i/%s%s_%ij.input_%s.root",dir_cards.Data(),mass,ana.Data(),fs.Data(),njets,cms.Data()));
      TIter next(inf->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();
	if (TString(h->GetName()).Contains("histo_Data")) continue;
	//cout << h->GetName() << endl;
	outf->cd();
	h->Write();
      }
      inf->Close();
      outf->Close();
    }

    //now get the data from mingshui's toys    
    TFile* msf = TFile::Open(Form("%s/%i/%s%s_%ij_%s_%s_PseudoData_sb.root",dir_toys.Data(),mass,ana.Data(),fs.Data(),njets,mode.Data(),cms.Data()));
    TH1F* msh = (TH1F*) msf->Get(Form("j%i%s_%i",njets,fs=="" ? "ll" : fs.Data(),i));

    if (mode=="shape") {
      outf = TFile::Open(Form("%s_N%i/%i/%s%s_%ij.input_%s.root",path.Data(),i,mass,ana.Data(),fs.Data(),njets,cms.Data()),"UPDATE");
    }
    TH1F* data =  (TH1F*) msh->Clone("histo_Data");

    ofstream outcard;
    outcard.open(Form("%s_N%i/%i/%s%s_%ij_%s_%s.txt",path.Data(),i,mass,ana.Data(),fs.Data(),njets,mode.Data(),cms.Data()));
    ifstream incard (Form("%s/%i/%s%s_%ij_%s_%s.txt",dir_cards.Data(),mass,ana.Data(),fs.Data(),njets,mode.Data(),cms.Data()));
    string line;
    if (incard.is_open()) {
      while ( incard.good() ) {
	getline (incard,line);
	size_t found=line.find("Observation");
	if (found!=string::npos) outcard << "Observation " << int(data->Integral()) << endl;
	//else if (line.find("FakeRate")!=string::npos) outcard << TString(line).ReplaceAll("1.360","9.999") << endl;
	else outcard << line << endl;
      }
      incard.close();
    }
    outcard.close();
    
    if (mode=="shape") {
      outf->cd();
      data->Write();
      outf->Close();
    }
  }

}

/*
root -b -q getPseudoData.C\(110,\"8TeV\",\"shape\",0,\"of\",0,1\)
*/
