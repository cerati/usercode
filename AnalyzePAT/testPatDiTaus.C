#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#endif

#include <iostream>
#include <cmath>      //necessary for absolute function fabs()

using namespace std;

void testPatDiTaus() {

//   gROOT-> ProcessLine(".x rootlogon.C");

  TFile  * file = new TFile("/tmp/cerati/tauAnalysisPatTupleElecMu.root");

  TFile  * outfile = new TFile("patDiTaus.root","RECREATE");

  TH1D * h_vismass = new TH1D("vismass", "vismass", 20, 0, 100 );
  fwlite::Event ev(file);

  //loop through each event
  for( ev.toBegin(); !ev.atEnd(); ++ev) {

    fwlite::Handle<std::vector<pat::Electron> > handleE;
    handleE.getByLabel(ev,"patElectrons");
    fwlite::Handle<std::vector<pat::Muon> > handleM;
    handleM.getByLabel(ev,"patMuons");

    if (!handleE.isValid() ) {cout << "could not get the collection!" << endl;return;}
    if (!handleM.isValid() ) {cout << "could not get the collection!" << endl;return;}
    vector<pat::Electron> const & electrons = *handleE;
    vector<pat::Muon> const & muons = *handleM;
    
    //loop through each Electron
    for (unsigned int i=0; i!=electrons.size(); ++i) {
      double me=0.00051;
      TLorentzVector p4e(electrons[i].px(),electrons[i].py(),electrons[i].pz(),sqrt(electrons[i].p()*electrons[i].p()+me*me));
      //loop through each Muon
      for (unsigned int j=0; j!=muons.size(); ++j) {
	double mm=0.10565;
	TLorentzVector p4m(muons[j].px(),muons[j].py(),muons[j].pz(),sqrt(muons[j].p()*muons[j].p()+mm*mm));
	TLorentzVector p4vis = p4e+p4m;
	h_vismass->Fill( p4vis.M() );
      } //end Muon loop
    }  //end Electron loop   
  }   //end event loop
  
  outfile->Write();
  outfile->Close();

}
