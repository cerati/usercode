#include <iostream>
#include <exception>
#include <set>

// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

// TAS includes
#include "./CORE/CMS2.cc"

#include "./CORE/trackSelections.cc"
#include "./CORE/eventSelections.cc"

#include "./CORE/muonSelections.cc"
#include "./CORE/electronSelections.cc"
#include "./CORE/electronSelectionsParameters.cc"
#include "./CORE/metSelections.cc"
#include "./CORE/utilities.cc"
#include "./CORE/jetSelections.cc"
#include "./Tools/goodrun.cc"
#include "./Tools/EgammaAnalysisTools/include/LikelihoodUtil.h"
#include "./CORE/conversionTools.cc"
#include "./CORE/MITConversionUtilities.cc"

#include "./CORE/triggerUtils.cc"
#include "./CORE/triggerSuperModel.cc"
#include "./myBabyMaker.h"

using namespace std;
using namespace tas;

// function for dR matching offline letpon to trigger object 
Int_t TriggerMatch( LorentzVector lepton_p4, const char* trigString, double dR_cut = 0.4 ){
	float dR_min = numeric_limits<float>::max();
	dR_min = 99.0;
	Int_t nTrig = nHLTObjects( trigString );
	//cout << trigString << " before : " <<   nTrig << ", " << dR_min << endl;
	if (nTrig > 0) {
		bool match = false;
		for (int itrg=0; itrg<nTrig; itrg++) {
			LorentzVector p4tr = p4HLTObject( trigString, itrg );
			double dr = ROOT::Math::VectorUtil::DeltaR( lepton_p4, p4tr);
			//cout << "dr :" << dr << endl;
			if (dr < dR_cut) match = true;
			if (dr < dR_min) dR_min = dr;
		}
		if (match) {
			nTrig = 2;
		} 
		else {
			nTrig = 1;
		}
	}
	//cout << trigString << " after : " <<   nTrig << ", " << dR_min << endl;
	return nTrig;
}


float globalJESRescale = 1.;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

/* THIS NEEDS TO BE IN CORE */

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

// transverse mass
float Mt( LorentzVector p4, float met, float met_phi ){
  return sqrt( 2*met*( p4.pt() - ( p4.Px()*cos(met_phi) + p4.Py()*sin(met_phi) ) ) );
}

void myBabyMaker::ScanChain( TChain* chain, const char *babyFilename, bool isData) {

  //cout << "1" << endl;
	
  already_seen.clear();

  // Make a baby ntuple
  MakeBabyNtuple(babyFilename);

  // Set the JSON file
  if(isData){
    set_goodrun_file("goodruns.txt");
  }

  // initiate likelihood object with pdf location
  LikelihoodUtil likelihoodUtil("./Tools/EgammaAnalysisTools/PDFs/pdfs_MC.root");

  //cout << "2" << endl;

  //--------------------------
  // File and Event Loop
  //---------------------------
  int i_permilleOld = 0;
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = 0;
  int nEvents = -1;
  if (nEvents==-1){
    nEventsChain = chain->GetEntries();
  } else {
    nEventsChain = nEvents;
  }

  //cout << "3" << endl;

  nEventsChain = chain->GetEntries();
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  map<int,int> m_events;
  while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
    TString filename = currentFile->GetTitle();

  //cout << "4" << endl;

    TFile f(filename.Data());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;
    unsigned int z;
    for( z = 0; z < nLoop; z++) 
    { // Event Loop
      cms2.GetEntry(z);

  //cout << "5" << endl;
      if(isData){
        // Good  Runs
        if(!goodrun( evt_run(), evt_lumiBlock() )) continue;

        // check for duplicated
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){ 
          cout << "\t! ERROR: found duplicate." << endl;
          continue;
        }
      }

      // looper progress
      ++nEventsTotal;
      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
      if (i_permille != i_permilleOld) {
        printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
        fflush(stdout);
        i_permilleOld = i_permille;
      }
    
  //cout << "6" << endl;

      // Event cleaning (careful, it requires technical bits)
      if (!cleaning_standardAugust2010(isData)) continue;

      // define eleId bits
      static const cuts_t electronSelection_VBTF95 = (1ll<<ELEID_VBTF_35X_95);
      static const cuts_t electronSelection_VBTF90 = (1ll<<ELEID_VBTF_35X_90);
      static const cuts_t electronSelection_VBTF80 = (1ll<<ELEID_VBTF_35X_80);
      static const cuts_t electronSelection_VBTF70 = (1ll<<ELEID_VBTF_35X_70);

  //cout << "7" << endl;
      InitBabyNtuple();  
      // loop over hyp
      for(unsigned int iEl = 0; iEl < els_p4().size(); ++iEl) 
      {


	if (els_gsftrkidx().at(iEl)<0) continue;
  //cout << "8" << endl;
	      
	      //  10 < pt < 30 GeV
	      if (els_p4().at(iEl).pt()<10) continue;
	      //if (els_p4().at(iEl).pt()>30) continue;

	      // ww electron selection
	      //if(!pass_electronSelection(iEl, electronSelection_wwV1_isoLoose, false, false)) continue;

	      // Z veto flag
	      for(unsigned int jEl = 0; jEl < iEl; jEl++) 
	      {
		      LorentzVector dil(els_p4().at(iEl)+els_p4().at(jEl));
		      float MZ = dil.mass();                           
		      if(MZ>76 && MZ<106) el_MZ_ = MZ;
	      }

	      //
	      // store in ntuples
	      //
	      //scale1fb_ = evt_scale1fb();
	      run_ = evt_run();
              ls_    = evt_lumiBlock();
	      evt_ = evt_event();

	      nvtx_ = 0;
	      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v) {
		if(isGoodVertex(v)) nvtx_++;
	      }

	      met_		= evt_pfmet();
	      met_phi_		= evt_pfmetPhi();

	      float  leadJetPt = -1.;
	      int iJet = -1;
	      for (int ij=0;ij<pfjets_p4().size();++ij){
		if (ROOT::Math::VectorUtil::DeltaR( pfjets_p4().at(ij), els_p4().at(iEl) )<1.0) continue;
		float jetPt = pfjets_p4().at(ij).pt()*pfjets_corL1FastL2L3().at(ij);
		if (jetPt>leadJetPt){
		  leadJetPt=jetPt;
		  iJet = ij;
		}
	      }
	      if (iJet>=0){
		jet_pt_ 	= pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet);
		jet_eta_ 	= pfjets_p4().at(iJet).eta();
		jet_phi_ 	= pfjets_p4().at(iJet).phi();
		jet_dR_         = ROOT::Math::VectorUtil::DeltaR( pfjets_p4().at(iJet), els_p4().at(iEl) );
	      } else {
		jet_pt_ 	= -999.;
		jet_eta_ 	= -999.;
		jet_phi_ 	= -999.;
		jet_dR_         = -999.;
	      }

	      el_q_ 	= els_charge().at(iEl);
	      el_pt_ 	= els_p4().at(iEl).pt();
	      el_px_ 	= els_p4().at(iEl).px();
	      el_py_ 	= els_p4().at(iEl).py();
	      el_pz_ 	= els_p4().at(iEl).pz();
	      el_eta_ 	= els_p4().at(iEl).eta();
	      el_phi_ 	= els_p4().at(iEl).phi();
	      el_mva_   = els_mva().at(iEl);
	      el_hOverE_ = els_hOverE().at(iEl); 
	      el_dEtaIn_ = els_dEtaIn().at(iEl); 
	      el_dPhiIn_ = els_dPhiIn().at(iEl); 
	      el_sigmaIEtaIEta_ = els_sigmaIEtaIEta().at(iEl); 
	      el_sigmaIPhiIPhi_ = els_sigmaIPhiIPhi().at(iEl); 
	      el_e2x5Max_ 	= els_e2x5Max().at(iEl); 
	      el_e1x5_ 		= els_e1x5().at(iEl); 
	      el_e5x5_ 		= els_e5x5().at(iEl); 
	      el_eSC_ 		= els_eSC().at(iEl);
	      el_etaSC_ 	= els_etaSC().at(iEl);
	      el_eOverPIn_ 	= els_eOverPIn().at(iEl); 
	      el_eOverPOut_ 	= els_eOverPOut().at(iEl);
	      el_fbrem_ 	= els_fbrem().at(iEl); 
	      el_nbrem_ 	= els_nSeed().at(iEl); 
	      el_Mt_		= Mt(els_p4().at(iEl), met_, met_phi_);
	      el_relIso_ 	= electronIsolation_rel(iEl, true);
	      el_ecalIso_ 	= els_ecalIso().at(iEl);
	      el_hcalIso_ 	= els_hcalIso().at(iEl);
	      el_tkIso_  	= els_tkIso().at(iEl);
	      el_ecalIso04_ 	= els_ecalIso04().at(iEl);
	      el_hcalIso04_ 	= els_hcalIso04().at(iEl);
	      el_tkIso04_  	= els_tkIso04().at(iEl);
	      el_pfIso03_ = electronIsoValuePF(iEl,0,0.3,1.0,0.1,0.07,0.025,0.025);
	      el_pfIso04_ = electronIsoValuePF(iEl,0,0.4,1.0,0.1,0.07,0.025,0.025);
	      el_innerlayer_	= els_exp_innerlayers().at(iEl);
	      el_conv_dist_	= els_conv_dist().at(iEl);
	      el_conv_dcot_	= els_conv_dcot().at(iEl);
	      //el_conv_rad_	= els_conv_radius().at(iEl);
	      el_conv_ddz_	= 0;
	      if (els_conv_tkidx().at(iEl)>=0) el_conv_ddz_	= trks_z0().at(els_conv_tkidx().at(iEl))-els_z0().at(iEl);
	      el_inner_posrho_	= els_inner_position().at(iEl).rho();
	      el_LL_		= likelihoodUtil.getValue(iEl);
	      if (el_LL_<=0) el_lh_=-20.;
	      else if (el_LL_==1) el_lh_=20.; 
	      else el_lh_            = log(el_LL_/(1.0-el_LL_));
	      el_type_		= els_type().at(iEl);

	      el_mitconv_ = isFromConversionMIT(iEl);

	      el_CICt_ 	        = electronId_CIC(iEl, 6, CIC_TIGHT, false, false); // 3? CIC_MEDIUM?
	      el_CICst_ 	= electronId_CIC(iEl, 6, CIC_SUPERTIGHT, false, false); // 3? CIC_MEDIUM?
	      el_CICht1_ 	= electronId_CIC(iEl, 6, CIC_HYPERTIGHT1, false, false); // 3? CIC_MEDIUM?
	      el_CICht2_ 	= electronId_CIC(iEl, 6, CIC_HYPERTIGHT2, false, false); // 3? CIC_MEDIUM?
	      el_CICht3_ 	= electronId_CIC(iEl, 6, CIC_HYPERTIGHT3, false, false); // 3? CIC_MEDIUM?
	      el_CICht4_ 	= electronId_CIC(iEl, 6, CIC_HYPERTIGHT4, false, false); // 3? CIC_MEDIUM?

	      if(el_pt_>10 && el_pt_<20 && (el_fbrem_>0.2 || (el_fbrem_<0.2 && abs(el_etaSC_)<1.479 && el_eOverPIn_>0.95)))      
  	      el_fbrem_eOverP_ = true;

	      el_class_ 	= els_class().at(iEl);
	      el_category_ 	= els_category().at(iEl);
	      //el_eg_roblooseid_ = els_egamma_robustLooseId().at(iEl);
	      //el_eg_robtightid_ = els_egamma_robustTightId().at(iEl);
	      //el_eg_looseid_ 	= els_egamma_looseId().at(iEl);
	      //el_eg_tightid_ 	= els_egamma_tightId().at(iEl);
	      //el_eg_robhighE_	= els_egamma_robustHighEnergy().at(iEl);

	      // id
	      el_VBTF70_ = pass_electronSelection(iEl, electronSelection_VBTF70, false, false);
	      el_VBTF80_ = pass_electronSelection(iEl, electronSelection_VBTF80, false, false);
	      el_VBTF90_ = pass_electronSelection(iEl, electronSelection_VBTF90, false, false);
	      el_VBTF95_ = pass_electronSelection(iEl, electronSelection_VBTF95, false, false);
   
	      // selection 
	      el_wwV0_ 	= pass_electronSelection(iEl, electronSelection_wwV0, false, false); 
	      el_wwV0b_ = pass_electronSelection(iEl, electronSelection_wwV0b, false, false); 
	      el_wwV1_ 	= pass_electronSelection(iEl, electronSelection_wwV1, false, false); 


	      // HLT mathcing dR<0.4. returns 2 if matching
	      el8_v1_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_v1");
	      el8_v2_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_v2");
	      el8_v3_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_v3");
	      el8idisojet40_v1_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1");
	      el8idisojet40_v2_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2");
	      el8idisojet40_v3_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3");
	      el8idiso_v1_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_v1");
	      el8idiso_v2_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_v2");
	      el8idiso_v3_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_v3");
	      el17idiso_v1_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_CaloIdL_CaloIsoVL_v1");
	      el17idiso_v2_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_CaloIdL_CaloIsoVL_v2");
	      el17idiso_v3_  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_CaloIdL_CaloIsoVL_v3");
	      el8pho20_v1_  = TriggerMatch( els_p4().at(iEl), "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1");
	      el8pho20_v2_  = TriggerMatch( els_p4().at(iEl), "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2");
	      el8pho20_v3_  = TriggerMatch( els_p4().at(iEl), "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3");

	      vector<ConversionInfo> v_convInfos = getConversionInfos(iEl, evt_bField(), 0.45);   
	      ConversionInfo bestConv = findBestConversionMatch(v_convInfos);
	      el_newconv_dist_=bestConv.dist();
	      el_newconv_dcot_=bestConv.dcot();
	      el_newconv_rad_=bestConv.radiusOfConversion();
	      el_newconv_dmh_=bestConv.deltaMissingHits();

	      el_validhits_ = els_validHits().at(iEl);
	      el_invalidhits_ = els_lostHits().at(iEl);
	      el_misshitsout_ = els_exp_outerlayers().at(iEl);
	      el_nchi2_ = els_chi2().at(iEl)/els_ndof().at(iEl);
	      el_d0corr_ = els_d0corr().at(iEl);
	      el_d0Err_ = els_d0Err().at(iEl);

	      el_d0pv2d_ = electron_d0PV_wwV1(iEl);
	      //el_d0pv3d_ = electron_d0PV3d_wwV1(iEl);

	      //el_d0pv2dErr_ = electron_d0PVErr_wwV1(iEl);
	      //el_d0pv3dErr_ = electron_d0PV3dErr_wwV1(iEl);
	      el_dzpv_ = gsftrks_dz_pv(els_gsftrkidx().at(iEl), 0, true).first;

	      FillBabyNtuple(); 

      }

    }// closes loop over events

  }  // closes loop over files

 
  CloseBabyNtuple();
  return;

} // closes myLooper function  


//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void myBabyMaker::InitBabyNtuple () 
{
  run_ = -1;
  ls_  = -1;
  evt_ = -1;
  scale1fb_ = -999.;
  nvtx_ = -1;

  met_ = -999.;
  met_phi_ = -999.;

  jet_pt_ = -999.;
  jet_eta_ = -999.;
  jet_phi_ = -999.;
  jet_dR_ = -999.;

  el_q_ = -1;
  el_pt_ = -999.;
  el_px_ = -999.;
  el_py_ = -999.;
  el_pz_ = -999.;
  el_eta_ = -999.;
  el_phi_ = -999.;
  el_mva_ = -999.;
  el_hOverE_ = -999.;
  el_dEtaIn_ = -999.;
  el_dPhiIn_ = -999.;
  el_sigmaIEtaIEta_ = -999.;
  el_sigmaIPhiIPhi_ = -999.;
  el_e2x5Max_ = -999.;
  el_e1x5_ = -999.;
  el_e5x5_ = -999.;
  el_eSC_ = -999.;
  el_etaSC_ = -999.;
  el_eOverPIn_ = -999.;
  el_eOverPOut_ = -999.;
  el_fbrem_ = -999.;
  el_nbrem_ = -999;
  el_Mt_= -999.;
  el_relIso_ = -999.;
  el_tkIso_ = -999.;
  el_ecalIso_ = -999.;
  el_hcalIso_ = -999.;
  el_tkIso04_ = -999.;
  el_ecalIso04_ = -999.;
  el_hcalIso04_ = -999.;
  el_innerlayer_ = -999.;
  el_conv_dist_ = -999.; 
  el_conv_dcot_ = -999.; 
  el_conv_rad_ = -999.; 
  el_conv_ddz_ = -999.; 
  el_inner_posrho_ = -999.; 
  el_LL_ = -999.; 
  el_lh_ = -999.; 
  el_type_ = -1; 
  el_CICt_ = -1; 
  el_CICst_ = -1; 
  el_CICht1_ = -1; 
  el_CICht2_ = -1; 
  el_CICht3_ = -1; 
  el_CICht4_ = -1; 

  el_mitconv_ = false;

  el_VBTF70_ = false;
  el_VBTF80_ = false;
  el_VBTF90_ = false;
  el_VBTF95_ = false;
  el_wwV0_   = false;
  el_wwV0b_  = false;
  el_wwV1_   = false;

  el_MZ_ = -999.;
 
  el_fbrem_eOverP_ = false;

  el_class_ = -1;
  el_category_ = -1;
  el_eg_roblooseid_ = -999. ;
  el_eg_robtightid_ = -999.;
  el_eg_looseid_ = -999.;
  el_eg_tightid_ = -999.;
  el_eg_robhighE_= -999.;
  
  el8_v1_ = -1; // HLT_Ele8_v1
  el8_v2_ = -1; // HLT_Ele8_v2
  el8_v3_ = -1; // HLT_Ele8_v3
  el8idisojet40_v1_ = -1; // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1
  el8idisojet40_v2_ = -1; // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2
  el8idisojet40_v3_ = -1; // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3
  el8idiso_v1_ = -1; // HLT_Ele8_CaloIdL_CaloIsoVL_v1
  el8idiso_v2_ = -1; // HLT_Ele8_CaloIdL_CaloIsoVL_v2
  el8idiso_v3_ = -1; // HLT_Ele8_CaloIdL_CaloIsoVL_v3
  el17idiso_v1_ = -1; // HLT_Ele17_CaloIdL_CaloIsoVL_v1
  el17idiso_v2_ = -1; // HLT_Ele17_CaloIdL_CaloIsoVL_v2
  el17idiso_v3_ = -1; // HLT_Ele17_CaloIdL_CaloIsoVL_v3
  el8pho20_v1_ = -1; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1
  el8pho20_v2_ = -1; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2
  el8pho20_v3_ = -1; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3

  el_newconv_dist_ = -999.; 
  el_newconv_dcot_ = -999.; 
  el_newconv_rad_ = -999.; 
  el_newconv_dmh_ = -999.; 

  el_validhits_ = -999;
  el_invalidhits_ = -999;
  el_misshitsout_ = -999;
  el_nchi2_ = -999.;
  el_d0corr_ = -999.;
  el_d0Err_ = -999.;
  el_d0pv2d_ = -999.;
  el_d0pv3d_ = -999.;
  el_d0pv2dErr_ = -999.;
  el_d0pv3dErr_ = -999.;
  el_dzpv_ = -999.;

  el_pfIso03_ = -999.;
  el_pfIso04_ = -999.;

}

//-------------------------------------
// Book the baby ntuple
//-------------------------------------
void myBabyMaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    babyTree_->Branch("run",          	&run_,         	"run/I"         );
    babyTree_->Branch("ls",           	&ls_,          	"ls/I"          );
    babyTree_->Branch("evt",          	&evt_,         	"evt/I"         );
    babyTree_->Branch("scale1fb",     	&scale1fb_,	"scale1fb/F"    );
    babyTree_->Branch("nvtx",          	&nvtx_,        	"nvtx/I"         );

    babyTree_->Branch("met",     		&met_,			"met/F");
    babyTree_->Branch("met_phi",     		&met_phi_,		"met_phi/F");

    babyTree_->Branch("jet_pt",     		&jet_pt_,		"jet_pt/F");
    babyTree_->Branch("jet_eta",     		&jet_eta_,		"jet_eta/F");
    babyTree_->Branch("jet_phi",     		&jet_phi_,		"jet_phi/F");
    babyTree_->Branch("jet_dR",     		&jet_dR_,		"jet_dR/F");

    babyTree_->Branch("el_q",     		&el_q_,			"el_q/I");
    babyTree_->Branch("el_pt",     		&el_pt_,		"el_pt/F");
    babyTree_->Branch("el_px",     		&el_px_,		"el_px/F");
    babyTree_->Branch("el_py",     		&el_py_,		"el_py/F");
    babyTree_->Branch("el_pz",     		&el_pz_,		"el_pz/F");
    babyTree_->Branch("el_eta",     		&el_eta_,		"el_eta/F");
    babyTree_->Branch("el_phi",     		&el_phi_,		"el_phi/F");
    babyTree_->Branch("el_mva",     		&el_mva_,		"el_mva/F");
    babyTree_->Branch("el_hOverE",     		&el_hOverE_,		"el_hOverE/F");
    babyTree_->Branch("el_dEtaIn",     		&el_dEtaIn_,		"el_dEtaIn/F");
    babyTree_->Branch("el_dPhiIn",     		&el_dPhiIn_,		"el_dPhiIn/F");
    babyTree_->Branch("el_sigmaIEtaIEta",     	&el_sigmaIEtaIEta_,	"el_sigmaIEtaIEta/F");
    babyTree_->Branch("el_sigmaIPhiIPhi",     	&el_sigmaIPhiIPhi_,	"el_sigmaIPhiIPhi/F");
    babyTree_->Branch("el_e2x5Max",     	&el_e2x5Max_,		"el_e2x5Max/F");
    babyTree_->Branch("el_e1x5",     		&el_e1x5_,		"el_e1x5/F");
    babyTree_->Branch("el_e5x5",     		&el_e5x5_,		"el_e5x5/F");
    babyTree_->Branch("el_eSC",     		&el_eSC_,		"el_eSC/F");
    babyTree_->Branch("el_etaSC",     		&el_etaSC_,		"el_etaSC/F");
    babyTree_->Branch("el_eOverPIn",     	&el_eOverPIn_,		"el_eOverPIn/F");
    babyTree_->Branch("el_eOverPOut",   	&el_eOverPOut_,		"el_eOverPOut/F");
    babyTree_->Branch("el_fbrem",     		&el_fbrem_,		"el_fbrem/F");
    babyTree_->Branch("el_nbrem",     		&el_nbrem_,		"el_nbrem/I");
    babyTree_->Branch("el_Mt",     		&el_Mt_,		"el_Mt/F");
    babyTree_->Branch("el_relIso",     		&el_relIso_,		"el_relIso/F");
    babyTree_->Branch("el_tkIso",     		&el_tkIso_,		"el_tkIso/F");
    babyTree_->Branch("el_ecalIso",    		&el_ecalIso_,		"el_ecalIso/F");
    babyTree_->Branch("el_hcalIso",    		&el_hcalIso_,		"el_hcalIso/F");
    babyTree_->Branch("el_tkIso04",     	&el_tkIso04_,		"el_tkIso04/F");
    babyTree_->Branch("el_ecalIso04",    	&el_ecalIso04_,		"el_ecalIso04/F");
    babyTree_->Branch("el_hcalIso04",   	&el_hcalIso04_,		"el_hcalIso04/F");
    babyTree_->Branch("el_innerlayer",   	&el_innerlayer_,	"el_innerlayer/F");
    babyTree_->Branch("el_conv_dist",   	&el_conv_dist_,		"el_conv_dist/F");
    babyTree_->Branch("el_conv_dcot",   	&el_conv_dcot_,		"el_conv_dcot/F");
    babyTree_->Branch("el_conv_rad",	   	&el_conv_rad_,		"el_conv_rad/F");
    babyTree_->Branch("el_conv_ddz",	   	&el_conv_ddz_,		"el_conv_ddz/F");
    babyTree_->Branch("el_inner_posrho",	&el_inner_posrho_,	"el_inner_posrho/F");
    babyTree_->Branch("el_LL",   		&el_LL_,		"el_LL/F");
    babyTree_->Branch("el_lh",   		&el_lh_,		"el_lh/F");
    babyTree_->Branch("el_type",   		&el_type_,		"el_type/I");
    babyTree_->Branch("el_CICt",   		&el_CICt_,		"el_CICt/I");
    babyTree_->Branch("el_CICst",   		&el_CICst_,		"el_CICst/I");
    babyTree_->Branch("el_CICht1",  		&el_CICht1_,		"el_CICht1/I");
    babyTree_->Branch("el_CICht2",  		&el_CICht2_,		"el_CICht2/I");
    babyTree_->Branch("el_CICht3",  		&el_CICht3_,		"el_CICht3/I");
    babyTree_->Branch("el_CICht4",  		&el_CICht4_,		"el_CICht4/I");

    babyTree_->Branch("el_VBTF70",          	&el_VBTF70_,        	"el_VBTF70/O"   );
    babyTree_->Branch("el_VBTF80",          	&el_VBTF80_,        	"el_VBTF80/O"   );
    babyTree_->Branch("el_VBTF90",          	&el_VBTF90_,        	"el_VBTF90/O"   );
    babyTree_->Branch("el_VBTF95",          	&el_VBTF95_,        	"el_VBTF95/O"   );
    babyTree_->Branch("el_wwV0",          	&el_wwV0_,        	"el_wwV0/O"	);
    babyTree_->Branch("el_wwV0b",          	&el_wwV0b_,        	"el_wwV0b/O"    );
    babyTree_->Branch("el_wwV1",          	&el_wwV1_,        	"el_wwV1/O"     );
    
    babyTree_->Branch("el_MZ",          	&el_MZ_,        	"el_MZ/F"     );
   
    babyTree_->Branch("el_fbrem_eOverP",   	&el_fbrem_eOverP_,  	"el_fbrem_eOverP/O"     );
   
    babyTree_->Branch("el_class",          	&el_class_,        	"el_class/I"     );
    babyTree_->Branch("el_category",          	&el_category_,        	"el_category/I"     );
    babyTree_->Branch("el_eg_roblooseid",      	&el_eg_roblooseid_,    	"el_eg_roblooseid/F"     );
    babyTree_->Branch("el_eg_robtightid",      	&el_eg_robtightid_,    	"el_eg_robtightid/F"     );
    babyTree_->Branch("el_eg_looseid",      	&el_eg_looseid_,    	"el_eg_looseid/F"     );
    babyTree_->Branch("el_eg_tightid",      	&el_eg_tightid_,    	"el_eg_tightid/F"     );
    babyTree_->Branch("el_eg_robhighE",      	&el_eg_robhighE_,    	"el_eg_robhighE/F"     );

    babyTree_->Branch("el8_v1",          	&el8_v1_,        	"el8_v1/I"	);
    babyTree_->Branch("el8_v2",          	&el8_v2_,        	"el8_v2/I"	);
    babyTree_->Branch("el8_v3",          	&el8_v3_,        	"el8_v3/I"	);
    babyTree_->Branch("el8idisojet40_v1",      	&el8idisojet40_v1_,     "el8idisojet40_v1/I"	);
    babyTree_->Branch("el8idisojet40_v2",     	&el8idisojet40_v2_,     "el8idisojet40_v2/I"	);
    babyTree_->Branch("el8idisojet40_v3",     	&el8idisojet40_v3_,     "el8idisojet40_v3/I"	);
    babyTree_->Branch("el8idiso_v1",          	&el8idiso_v1_,        	"el8idiso_v1/I"	);
    babyTree_->Branch("el8idiso_v2",          	&el8idiso_v2_,        	"el8idiso_v2/I"	);
    babyTree_->Branch("el8idiso_v3",          	&el8idiso_v3_,        	"el8idiso_v3/I"	);
    babyTree_->Branch("el17idiso_v1",          	&el17idiso_v1_,        	"el17idiso_v1/I"	);
    babyTree_->Branch("el17idiso_v2",          	&el17idiso_v2_,        	"el17idiso_v2/I"	);
    babyTree_->Branch("el17idiso_v3",          	&el17idiso_v3_,        	"el17idiso_v3/I"	);
    babyTree_->Branch("el8pho20_v1",          	&el8pho20_v1_,        	"el8pho20_v1/I"	);
    babyTree_->Branch("el8pho20_v2",          	&el8pho20_v2_,        	"el8pho20_v2/I"	);
    babyTree_->Branch("el8pho20_v3",          	&el8pho20_v3_,        	"el8pho20_v3/I"	);

    babyTree_->Branch("el_newconv_dist",   	&el_newconv_dist_,		"el_newconv_dist/F");
    babyTree_->Branch("el_newconv_dcot",   	&el_newconv_dcot_,		"el_newconv_dcot/F");
    babyTree_->Branch("el_newconv_rad",	   	&el_newconv_rad_,		"el_newconv_rad/F");
    babyTree_->Branch("el_newconv_dmh",	   	&el_newconv_dmh_,		"el_newconv_dmh/F");

    babyTree_->Branch("el_mitconv",          	&el_mitconv_,        	"el_mitconv/O"   );

    babyTree_->Branch("el_validhits",	   	&el_validhits_,		"el_validhits/I");
    babyTree_->Branch("el_invalidhits",	   	&el_invalidhits_,		"el_invalidhits/I");
    babyTree_->Branch("el_misshitsout",	   	&el_misshitsout_,		"el_misshitsout/I");
    babyTree_->Branch("el_nchi2",	   	&el_nchi2_,		"el_nchi2/F");
    babyTree_->Branch("el_d0corr",	   	&el_d0corr_,		"el_d0corr/F");
    babyTree_->Branch("el_d0Err",	   	&el_d0Err_,		"el_d0Err/F");
    babyTree_->Branch("el_dzpv",	   	&el_dzpv_,		"el_dzpv/F");

    babyTree_->Branch("el_d0pv2d",	   	&el_d0pv2d_,		"el_d0pv2d/F");
    babyTree_->Branch("el_d0pv3d",	   	&el_d0pv3d_,		"el_d0pv3d/F");
    babyTree_->Branch("el_d0pv2dErr",	   	&el_d0pv2dErr_,		"el_d0pv2dErr/F");
    babyTree_->Branch("el_d0pv3dErr",	   	&el_d0pv3dErr_,		"el_d0pv3dErr/F");

    babyTree_->Branch("el_pfIso03",     	&el_pfIso03_,		"el_pfIso03/F");
    babyTree_->Branch("el_pfIso04",     	&el_pfIso04_,		"el_pfIso04/F");

}
//----------------------------------
// Fill the baby
//----------------------------------
void myBabyMaker::FillBabyNtuple()
{
    babyTree_->Fill();
}
//--------------------------------
// Close the baby
//--------------------------------
void myBabyMaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}

