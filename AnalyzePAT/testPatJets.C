#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Jet.h"
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatJets() {

  gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patJets.root","RECREATE");

  TH1D * h_pt = new TH1D("pt", "Jet p_{T}", 20, 0, 100 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {
    fwlite::Handle<std::vector<pat::Jet> > handle;
    handle.getByLabel(ev,"patJets");

    if (!handle.isValid() ) {cout << "could not get the collection!" << endl;return;}
    vector<pat::Jet> const & jets = *handle;
    
    //loop through each Jet
    for (unsigned int i=0; i!=jets.size(); ++i) {
      h_pt->Fill( jets[i].pt() );
    }   //end Jet loop   

  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
