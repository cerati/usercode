#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Muon.h"
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatMuons() {

  gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patMuons.root","RECREATE");

  TH1D * h_pt = new TH1D("pt", "Muon p_{T}", 20, 0, 100 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {
    fwlite::Handle<std::vector<pat::Muon> > handle;
    handle.getByLabel(ev,"patMuons");

    if (!handle.isValid() ) {cout << "could not get the collection!" << endl;return;}
    vector<pat::Muon> const & muons = *handle;
    
    //loop through each Muon
    for (unsigned int i=0; i!=muons.size(); ++i) {
      h_pt->Fill( muons[i].pt() );
    }   //end Muon loop   

  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
