// C++
#include <iostream>
#include <vector>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TMath.h"

#include "./looper.h"

// CMS2
#include "./CORE/conversionTools.h"
#include "./CORE/electronSelections.h"
#include "./CORE/eventSelections.h"
#include "./CORE/jetSelections.h"
#include "./CORE/mcSelections.h"
#include "./CORE/metSelections.h"
#include "./CORE/MITConversionUtilities.h"
#include "./CORE/muonSelections.h"
#include "./CORE/trackSelections.h"
#include "./CORE/triggerUtils.h"
#include "./CORE/utilities.h"

using namespace tas;
using namespace std;

int looper::ScanChain( TChain* chain, const char* prefix, bool isData, int nEvents) {

  makebaby       = false;
  makehist       = false;
  maketext       = false;
  
  if (makebaby) MakeBabyNtuple( Form( "%s_baby.root", prefix ) );
  if (makehist) BookHistos( Form( "%s_histos.root", prefix ) );

  // File Loop
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    // Event Loop
    unsigned int nEvents = tree->GetEntries();

    for( unsigned int event = 0; event < nEvents; ++event) {
    
      // Get Event Content
      cms2.GetEntry(event);
      ++nEventsTotal;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
	  fflush(stdout);
	}
      }

      if (makebaby) InitBabyNtuple();

      //-------- DO THE ANALYSIS HERE --------//

      //print info
      if (maketext) printEvent();

      //fill baby
      run_   = evt_run();
      ls_    = evt_lumiBlock();
      evt_   = evt_event();
      weight_ = isData ? 1. : evt_scale1fb();

      for (unsigned int i_hyp = 0; i_hyp<hyp_p4().size(); i_hyp++) {

	*p4_ = hyp_p4().at(i_hyp);
	if (makebaby) FillBabyNtuple();
	
      }

      //fill histos
      if (makehist) h_dummy->Fill(1);

    }

    delete tree;
    f.Close();
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  if (makebaby) CloseBabyNtuple();
  if (makehist) SaveHistos();

  return 0;
}

void looper::printEvent(  ostream& ostr ){
  ostr << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl; 
}

// Book the baby ntuple
void looper::MakeBabyNtuple(const char *babyFilename) {
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("tree", "A Baby Ntuple");
  
  babyTree_->Branch("run", &run_, "run/I");
  babyTree_->Branch("ls", &ls_, "ls/I");
  babyTree_->Branch("evt", &evt_, "evt/I");
  babyTree_->Branch("weight", &weight_, "weight/F");

  babyTree_->Branch("p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &p4_);
}

// Init the baby
void looper::InitBabyNtuple() { 
  run_ = -9999;
  ls_ = -9999;
  evt_ = -9999;
  weight_ = -9999.;

  babyTree_->SetBranchAddress("p4",          &p4_);
  *p4_ = LorentzVector();

}

void looper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}
