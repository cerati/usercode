#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// #include "DataFormats/PatCandidates/interface/Vertex.h"
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatVertex() {

  gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patVertex.root","RECREATE");

  TH1D * h_ndof = new TH1D("ndof", "ndof", 20, 0, 100 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {
    fwlite::Handle<std::vector<reco::Vertex> > handle;
    handle.getByLabel(ev,"offlinePrimaryVertices");

    if (!handle.isValid() ) {cout << "could not get the collection!" << endl;return;}
    vector<reco::Vertex> const & vertices = *handle;
    
    //loop through each Vertex
    for (unsigned int i=0; i!=vertices.size(); ++i) {
      h_ndof->Fill( vertices[i].ndof() );
    }   //end Vertex loop   

  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
