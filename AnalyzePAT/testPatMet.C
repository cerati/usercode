#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Met.h"
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatMet() {

  gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patMet.root","RECREATE");

  TH1D * h_pt = new TH1D("pt", "Met p_{T}", 20, 0, 100 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {
    fwlite::Handle<std::vector<pat::MET> > handle;
    handle.getByLabel(ev,"patPFMETs");

    if (!handle.isValid() ) {cout << "could not get the collection!" << endl;return;}
    vector<pat::MET> const & mets = *handle;
    
    //loop through each Met
    for (unsigned int i=0; i!=mets.size(); ++i) {
      h_pt->Fill( mets[i].pt() );
    }   //end Met loop   

  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
