#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TKey.h"
#include "TClass.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "histTools.C"
#include <iostream>
#include <fstream>
#include <sstream>

void mixCards(){

  gStyle->SetOptStat(0);

  TCanvas c1;

  TFile *_file_G = TFile::Open("hwwof_0j.input_8TeV.root");
  TFile *_file_g = TFile::Open("hwwof_0j.input.root");

  TFile *_file_new = TFile::Open("hwwof_0j.input_mixtest.root","RECREATE");
  _file_new->cd();

  TString discriminator = "MVAMETResBounding";

  TIter next_G(_file_G->GetListOfKeys());
  TKey *key_G;
  while ((key_G = (TKey*)next_G())) {
    TClass *cl = gROOT->GetClass(key_G->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_G = (TH1F*)key_G->ReadObj();
    if (TString(h_G->GetName()).Contains(discriminator.Data())==0) continue;
    TH1F* h_G_clone = new TH1F();
    h_G->Copy(*h_G_clone);
    h_G_clone->Write();
  }

  TIter next_g(_file_g->GetListOfKeys());
  TKey *key_g;
  while ((key_g = (TKey*)next_g())) {
    TClass *cl = gROOT->GetClass(key_g->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_g = (TH1F*)key_g->ReadObj();
    if (TString(h_g->GetName()).Contains(discriminator.Data())) continue;
    TH1F* h_g_clone = new TH1F();
    h_g->Copy(*h_g_clone);
    h_g_clone->Write();
  }

  _file_new->Close();

}



//void hackCards(TString card, TString alter, TString central){
void alterShapeToCentral(TString card, TString alter, TString central){

  gStyle->SetOptStat(0);

  TCanvas c1;

  gSystem->Exec("mv "+card+" tmp.root");

  TFile *_file_old = TFile::Open("tmp.root");

  TFile *_file_new = TFile::Open(card,"RECREATE");
  _file_new->cd();

  TH1F* h_central = new TH1F();
  TH1F* h_alter = new TH1F();

  TIter next_old(_file_old->GetListOfKeys());
  TKey *key_old;
  while ((key_old = (TKey*)next_old())) {
    TClass *cl = gROOT->GetClass(key_old->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_old = (TH1F*)key_old->ReadObj();
    if (TString(h_old->GetName())==central.Data()){
      continue;
    }
    if (TString(h_old->GetName())==alter.Data()){
      h_old->Copy(*h_alter); 
      continue;
    }
    TH1F* h_old_clone = new TH1F();
    h_old->Copy(*h_old_clone);
    h_old_clone->Write();
  }

  h_central = (TH1F*) h_alter->Clone(central.Data());
  h_central->SetTitle(central.Data());
  h_central->Write();
  h_alter->Write();

  _file_new->Close();

  gSystem->Exec("rm tmp.root");

}


//void hackCards(TString card, TString alter, TString central){
void alterShapeToCentralWithSyst(TString card, TString alter, TString central){

  //fixme: need to account for more than one alternative shape (for now it works only in statonly mode)

  gStyle->SetOptStat(0);

  TCanvas c1;

  gSystem->Exec("mv "+card+" tmp.root");

  TFile *_file_old = TFile::Open("tmp.root");

  TFile *_file_new = TFile::Open(card,"RECREATE");
  _file_new->cd();

  TH1F* h_central = new TH1F();
  TH1F* h_alterUp = new TH1F();
  TH1F* h_alterDn = new TH1F();

  TIter next_old(_file_old->GetListOfKeys());
  TKey *key_old;
  while ((key_old = (TKey*)next_old())) {
    TClass *cl = gROOT->GetClass(key_old->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_old = (TH1F*)key_old->ReadObj();
    if (TString(h_old->GetName())==central.Data()){
      h_old->Copy(*h_central);
      continue;
    }
    if (TString(h_old->GetName()).Contains(alter.Data())){
      if (TString(h_old->GetName()).Contains("Up"))   h_old->Copy(*h_alterUp); 
      if (TString(h_old->GetName()).Contains("Down")) h_old->Copy(*h_alterDn); 
      continue;
    }
    TH1F* h_old_clone = new TH1F();
    h_old->Copy(*h_old_clone);
    h_old_clone->Write();
  }

  TH1F* h_central_new = (TH1F*) h_alterUp->Clone(central.Data());
  h_central_new->SetTitle(central.Data());
  h_central_new->Write();

  TH1F* h_alterUp_new = (TH1F*) h_central->Clone(Form("%s%s",alter.Data(),"Up"));
  h_alterUp_new->SetTitle(Form("%s%s",alter.Data(),"Up"));
  h_alterUp_new->Write();

  TH1F* h_alterDn_new = (TH1F*) h_central->Clone(Form("%s%s",alter.Data(),"Down"));
  h_alterDn_new->SetTitle(Form("%s%s",alter.Data(),"Down"));
  h_alterDn_new->Reset();
  fillDownMirrorUp(h_central_new,h_alterUp_new,h_alterDn_new);
  h_alterDn_new->Write();

  _file_new->Close();

  gSystem->Exec("rm tmp.root");

}

//void scaleProcess(TString cardtxt, TString process, float scale) {
void hackCards(TString cardtxt, TString process, float scale) {

  //scale the yield
  gSystem->Exec("mv "+cardtxt+" tmp.txt");
  ofstream outcard;
  outcard.open(cardtxt);
  ifstream incard("tmp.txt");
  string line;
  int it_proc = -1;
  if (incard.is_open()) {
    while ( incard.good() ) {
      getline (incard,line);
      //cout << line << endl;
      if (line.find("rate")!=string::npos) {
        TString myline(line);
	TObjArray tok_str = *myline.Tokenize(' ');
	for (int i=0;i<tok_str.GetEntries();++i) {
	  TObjString* oldrate = (TObjString*) tok_str.At(i);
	  if (it_proc==i) {
	    float newrate = oldrate->GetString().Atof()*scale;
	    outcard << newrate << " ";
	  }
	  else outcard << oldrate->GetString() << " ";
	}
	outcard <<  endl;
      } else {
	if (line.find("process")!=string::npos) {
	  TString myline(line);
	  TObjArray tok_str = *myline.Tokenize(' ');
	  for (int i=0;i<tok_str.GetEntries();++i) {
	    TObjString* token = (TObjString*) tok_str.At(i);
	    //cout << token->GetString() << endl;
	    if (token->GetString()==process) {
	      it_proc=i;
	      //cout << "found proc=" << it_proc << endl;
	      break;
	    }
	  }
	}
	outcard << line << endl;
      }
    }
    incard.close();
  }
  outcard.close();
  gSystem->Exec("rm tmp.txt");

  //scale the shapes
  if (cardtxt.Contains("shape")) {
    TString card = cardtxt;
    card.ReplaceAll("_shape",".input");
    card.ReplaceAll("txt","root");
    gSystem->Exec("mv "+card+" tmp.root");
    TFile *_file_old = TFile::Open("tmp.root");
    TFile *_file_new = TFile::Open(card,"RECREATE");
    _file_new->cd();
    TIter next_old(_file_old->GetListOfKeys());
    TKey *key_old;
    while ((key_old = (TKey*)next_old())) {
      TClass *cl = gROOT->GetClass(key_old->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TH1F *h_old = (TH1F*)key_old->ReadObj();
      TH1F* h_old_clone = new TH1F();
      h_old->Copy(*h_old_clone);
      if (TString(h_old->GetName()).Contains(process.Data())){
	h_old_clone->Scale(scale); 
      }
      h_old_clone->Write();
    }
    _file_new->Close();
    gSystem->Exec("rm tmp.root");
  }

}

//for f in cards_inj_stat_wj_alter/*/hww*.root; do root -b -q hackCards.C+\(\"${f}\",\"histo_Wjets_CMS_hww_MVAWBoundingUp\",\"histo_Wjets\"\); done
//for f in cards_inj_nuin_wj_alter/*/hww*.root; do root -b -q hackCards.C+\(\"${f}\",\"histo_Wjets_CMS_hww_MVAWBounding\",\"histo_Wjets\"\); done
//for f in cards_inj_nuin_wj_scale/*/hww*.txt; do root -b -q hackCards.C+\(\"${f}\",\"Wjets\",1.360\); done
