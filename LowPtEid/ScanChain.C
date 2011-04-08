/* Usage:
   root [0] .L ScanChain.C++
   root [1] TChain *chain = new TChain("Events")
   root [2] chain->Add("merged_ntuple.root")
   root [3] ScanChain(chain)
*/

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TMath.h"

// CMS2
#include "../CORE/CMS2.h"
#include "../CORE/electronSelections.h"
#include "../CORE/MITConversionUtilities.h"

using namespace tas;


int ScanChain( TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

  // File Loop
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  int   mitPar1[] = {1,1,1,1,1,1,2,999,0,0,0,0,0,0};
  float mitPar2[] = {1e-6,1e-6,1e-6,1e-6,1e-6,1e-10,1e-10,1e-10,1e-6,1e-6,1e-6,1e-6,1e-6,1e-10};
  float mitPar3[] = {2.0,2.0,2.0,-1.0,-999.9,-999.9,-999.9,-999.9,2.0,2.0,2.0,-1.0,-999.9,-999.9};
  int   mitPar4[] = {0,0,1,1,1,1,1,1,0,0,1,1,1,1};
  int   mitPar5[] = {1,0,0,0,0,0,0,0,1,0,0,0,0,0};
  int nmit = sizeof(mitPar1)/sizeof(int);

  float sntDist[] = {0.02,0.02,0.05,0.05,0.02,0.02,0.05,0.05,0.1};
  float sntDcot[] = {0.02,0.02,0.02,0.02,0.05,0.05,0.05,0.05,0.1};
  int   sntDmht[] = {2   ,1   ,2   ,1   ,2   ,1   ,2   ,1   ,2  };
  int nsnt = sizeof(sntDmht)/sizeof(int);

  vector<int> v_nEls;
  v_nEls.reserve(6);
  for(unsigned int i = 0; i < 6; i++) v_nEls.push_back(0);
  vector<int> nEls_numMIT;
  vector<int> nEls_numSnT;
  vector<int> nEls_numSnT_old;
  for (int i=0;i<nmit;++i) nEls_numMIT.push_back(0);
  for (int i=0;i<nsnt;++i) nEls_numSnT.push_back(0);
  for (int i=0;i<nsnt;++i) nEls_numSnT_old.push_back(0);

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

      for(unsigned int i_els = 0; i_els < els_p4().size(); i_els++) {

	int cut_counter = 0;
	if(els_p4()[i_els].Pt() < 20.) continue;
	v_nEls[cut_counter]++;
	cut_counter++;

	//if(abs(els_mc_id()[i_els]) !=22) continue;
	//if(abs(els_mc3_id()[i_els]) != 11) continue;
	v_nEls[cut_counter]++;
	cut_counter++;

	
	if(electronIsolation_rel_v1(i_els, true) > 0.10) continue;
	//if(fabs(els_etaSC()[i_els])<1.479 && electronIsolation_rel_v1(i_els, true) > 0.07) continue;
	//else if (fabs(els_etaSC()[i_els])>1.479 && electronIsolation_rel_v1(i_els, true) > 0.06) continue;
	v_nEls[cut_counter]++;
	cut_counter++;

	unsigned int answer_vbtf = electronId_VBTF(i_els, VBTF_35X_80, false, false);
	if( ( answer_vbtf & (1ll<<ELEID_ID) ) == (1ll<<ELEID_ID) != true) continue;

	v_nEls[cut_counter]++;
	cut_counter++;

	if(electron_d0PV_mindz(i_els) > 0.02) continue;

	v_nEls[cut_counter]++;
	cut_counter++;


	if(els_exp_innerlayers()[i_els] != 0) continue;
	v_nEls[cut_counter]++;
	cut_counter++;
	
	for (int isnt=0;isnt<nsnt;isnt++) {
	  if(els_conv_flag()[i_els] < 0) continue;

	  //compute delta z0 between electron and partner track
	  //version 1, take it from the tracks that are paired by the convInfo algo
	  float dz0 = 999.;
	  if (els_conv_flag()[i_els]==0) dz0 = fabs(trks_z0().at(els_conv_tkidx()[i_els]) - trks_z0().at(els_trkidx()[i_els]));
	  else if (els_conv_flag()[i_els]==1) dz0 = fabs(trks_z0().at(els_conv_tkidx()[i_els]) - gsftrks_z0().at(els_gsftrkidx()[i_els]));
	  else if (els_conv_flag()[i_els]==2) dz0 = fabs(gsftrks_z0().at(els_conv_gsftkidx()[i_els]) - trks_z0().at(els_trkidx()[i_els]));
	  else if (els_conv_flag()[i_els]==3) dz0 = fabs(gsftrks_z0().at(els_conv_gsftkidx()[i_els]) - gsftrks_z0().at(els_gsftrkidx()[i_els]));

	  //compute delta z0 between electron and partner track
	  //version 2, take it the min delta z0 between all available tracks
	  //(the idea is that dz might be badly computed for the ctf electron track and thus we make the cut a bit looser) 
	  float mindz0 = 999., dz01 = 999., dz02 = 999., dz03 = 999., dz04 = 999.;
	  if (els_conv_tkidx()[i_els]>-1 && els_trkidx()[i_els]>-1) 
	    dz01 = fabs(trks_z0().at(els_conv_tkidx()[i_els]) - trks_z0().at(els_trkidx()[i_els]));
	  if (els_conv_tkidx()[i_els]>-1 && els_gsftrkidx()[i_els]>-1) 
	    dz02 = fabs(trks_z0().at(els_conv_tkidx()[i_els]) - gsftrks_z0().at(els_gsftrkidx()[i_els]));
	  if (els_conv_gsftkidx()[i_els]>-1 && els_trkidx()[i_els]>-1) 
	    dz03 = fabs(gsftrks_z0().at(els_conv_gsftkidx()[i_els]) - trks_z0().at(els_trkidx()[i_els]));
	  if (els_conv_gsftkidx()[i_els]>-1 && els_gsftrkidx()[i_els]>-1) 
	    dz04 = fabs(gsftrks_z0().at(els_conv_gsftkidx()[i_els]) - gsftrks_z0().at(els_gsftrkidx()[i_els]));
	  mindz0 = min(dz01,min(dz02,min(dz03,dz04)));

	  bool domin = true;
	  if (domin) dz0 = mindz0; 

	  if (fabs(els_conv_dist()[i_els]) < sntDist[isnt]
	      && fabs(els_conv_dcot()[i_els]) < sntDcot[isnt] 
	      && abs(els_conv_delMissHits()[i_els]) < sntDmht[isnt]
	      && dz0 < 0.2 )
	    nEls_numSnT[isnt]++;

	  bool tests = false;
	  if (tests) {
	    if (dz0-mindz0>0.1)
	      cout << dz0 << " " << mindz0 << " " << dz01 << " " << dz02 << " " << dz03 << " " << dz04 << " " 
		   << els_conv_flag()[i_els] << " " << trks_z0().at(els_trkidx()[i_els]) << " " << gsftrks_z0().at(els_gsftrkidx()[i_els]) 
		   << " " << trks_trk_p4().at(els_trkidx()[i_els]).pt() << " " << gsftrks_p4().at(els_gsftrkidx()[i_els]).pt()
		   << " " << trks_trk_p4().at(els_trkidx()[i_els]).eta() << " " << gsftrks_p4().at(els_gsftrkidx()[i_els]).eta()
		   << " " << trks_trk_p4().at(els_trkidx()[i_els]).phi() << " " << gsftrks_p4().at(els_gsftrkidx()[i_els]).phi()
		   << endl;
	  }
	}

	for (int isnt=0;isnt<nsnt;isnt++) {
	  if(els_conv_flag()[i_els] > -1 
	     && fabs(els_conv_old_dist()[i_els]) < sntDist[isnt]
	     && fabs(els_conv_old_dcot()[i_els]) < sntDcot[isnt] 
	     && abs(els_conv_old_delMissHits()[i_els]) < sntDmht[isnt])
	    nEls_numSnT_old[isnt]++;
	}

	//MIT stuff
	for (int imit=0;imit<nmit;imit++) {
	  if(isMITConversion(i_els, mitPar1[imit], mitPar2[imit], mitPar3[imit], mitPar4[imit], mitPar5[imit])) 
	    nEls_numMIT[imit]++;
	}

      }

    
    }

    delete tree;
    f.Close();
  }
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  cout << "After pt:" << v_nEls[0] << endl;
  cout << "After gen match:" << v_nEls[1] << endl;
  cout << "After iso:" << v_nEls[2] << endl;    
  cout << "After id:" << v_nEls[3] << endl;    
  cout << "After d0:" << v_nEls[4] << endl;    
  cout << "After exp hits:" << v_nEls[5] << endl;    

  cout << endl;


  for (int isnt=0;isnt<nsnt;isnt++) {
    cout << "Number not ided as conversions, old_SnTv" << isnt+1 << ":" << v_nEls[5] - nEls_numSnT_old[isnt]
	 << "(" << (float)(v_nEls[5] - nEls_numSnT_old[isnt])/v_nEls[5] << ")" << endl;
  }

  cout << endl;

  for (int isnt=0;isnt<nsnt;isnt++) {
    cout << "Number not ided as conversions, SnTv" << isnt+1 << ":" << v_nEls[5] - nEls_numSnT[isnt]
	 << "(" << (float)(v_nEls[5] - nEls_numSnT[isnt])/v_nEls[5] << ")" << endl;
  }

  cout << endl;
 
  for (int imit=0;imit<nmit;imit++) {
    cout << "Number not ided as conversions, MITv" << imit+1 << ":" << v_nEls[5] - nEls_numMIT[imit]
	 << "(" << (float)(v_nEls[5] - nEls_numMIT[imit])/v_nEls[5] << ")" << endl;
  }

  return 0;
}
