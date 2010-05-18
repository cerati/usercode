#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// #include "DataFormats/PatCandidates/interface/Electron.h"
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatElectrons() {

  gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patElectrons.root","RECREATE");

  TH1D * h_pt = new TH1D("pt", "Electron p_{T}", 20, 0, 100 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {
    fwlite::Handle<std::vector<pat::Electron> > handle;
    handle.getByLabel(ev,"patElectrons");

    if (!handle.isValid() ) {cout << "could not get the collection!" << endl;return;}
    vector<pat::Electron> const & electrons = *handle;
    
    //loop through each Electron
    for (unsigned int i=0; i!=electrons.size(); ++i) {
      h_pt->Fill( electrons[i].pt() );
    }   //end Electron loop   

  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
