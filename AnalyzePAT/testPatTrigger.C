#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatTrigger() {

  gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patTrigger.root","RECREATE");

  TH1D * h_hlt = new TH1D("hlt", "HLT bits", 150, -0.5, 149.5 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {

    std::vector<bool> trigBits;
    fwlite::Handle<edm::TriggerResults> hltresults;
    hltresults.getByLabel(ev,"TriggerResults","","HLT");
    edm::TriggerNames triggerNames_ = ev.triggerNames(*hltresults);
    int ntrigs=hltresults->size();
    vector<string> triggernames = triggerNames_.triggerNames();
    //int trigjolly = 0;
    for (int itrig = 0; itrig != ntrigs; ++itrig){
      string trigName=triggerNames_.triggerName(itrig);
      bool accept = hltresults->accept(itrig);
      if (accept) h_hlt->Fill(itrig);
      //cout << trigName << " " << accept << endl;
    }

  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
