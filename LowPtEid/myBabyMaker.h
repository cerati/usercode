#ifndef myBabyMaker_h
#define myBabyMaker_h

#include <stdint.h>
#include "TFile.h"
#include "TTree.h"

typedef ULong64_t   uint64;
typedef uint64      electronIdComponent_t;

class myBabyMaker
{
    public:
        myBabyMaker() {};
        ~myBabyMaker() {
            delete babyFile_;
            delete babyTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain *, const char *, bool);

    private:
        //
        // BABY NTUPLE VARIABLES
        //
        TFile *babyFile_;
	TTree *babyTree_;

	// event var 
	Int_t   run_;
	Int_t   ls_;
	Int_t   evt_;
	Float_t scale1fb_;
	Int_t   nvtx_;

	//jet
	Float_t jet_pt_;
	Float_t jet_eta_;
	Float_t jet_phi_;
	Float_t jet_dR_;

	// charge
	Int_t   el_q_;

	// kinematic 
	Float_t el_px_;
	Float_t el_py_;
	Float_t el_pz_;
	Float_t el_pt_;
	Float_t el_eta_;
	Float_t el_phi_;

	// multi-variate
	Float_t el_mva_;

	// id var
	Float_t el_hOverE_;
	Float_t el_dEtaIn_;
	Float_t el_dPhiIn_;
	Float_t el_sigmaIEtaIEta_;
	Float_t el_sigmaIPhiIPhi_;

	// cluster var
	Float_t el_e2x5Max_;
	Float_t el_e1x5_;
	Float_t el_e5x5_;
	Float_t el_eSC_;
	Float_t el_etaSC_;

	// e/p
	Float_t el_eOverPIn_;
	Float_t el_eOverPOut_;

	// brem
	Float_t el_fbrem_;
	Int_t el_nbrem_;

	// transverse mass
	Float_t el_Mt_;

	// isolation
	Float_t el_relIso_;
	Float_t el_tkIso_;
	Float_t el_ecalIso_;
	Float_t el_hcalIso_;

	Float_t el_tkIso04_;
	Float_t el_ecalIso04_;
	Float_t el_hcalIso04_;

	Float_t el_pfIso03_;
	Float_t el_pfIso04_;

	// conversion
	Float_t el_innerlayer_;
	Float_t el_conv_dist_;
	Float_t el_conv_dcot_;
	Float_t el_conv_rad_;
	Float_t el_conv_ddz_;
	Float_t el_inner_posrho_;
	// new conversion
	Float_t el_newconv_dist_;
	Float_t el_newconv_dcot_;
	Float_t el_newconv_rad_;
	Float_t el_newconv_dmh_;
	// mit conversion
	Bool_t el_mitconv_;

	// hits and tk additional stuff
	Int_t el_validhits_;
	Int_t el_invalidhits_;
	Int_t el_misshitsout_;
	Float_t el_nchi2_;

	Float_t el_dzpv_;

	Float_t el_d0corr_;
	Float_t el_d0Err_;
	Float_t el_d0pv2d_;
	Float_t el_d0pv3d_;
	Float_t el_d0pv2dErr_;
	Float_t el_d0pv3dErr_;

	// likelihood
	Float_t el_LL_;
	Float_t el_lh_;

        // seeding
	Int_t el_type_;

	// electron CIC
	electronIdComponent_t el_CICt_;
	electronIdComponent_t el_CICst_;
	electronIdComponent_t el_CICht1_;
	electronIdComponent_t el_CICht2_;
	electronIdComponent_t el_CICht3_;
	electronIdComponent_t el_CICht4_;

	// met
	Float_t met_;
	Float_t met_phi_;

	// id and selection
	Bool_t el_VBTF70_;
	Bool_t el_VBTF80_;
	Bool_t el_VBTF90_;
	Bool_t el_VBTF95_;
	Bool_t el_wwV0_;
	Bool_t el_wwV0b_;
	Bool_t el_wwV1_;

	// dilep mass
	Float_t el_MZ_;

	// new cut 10<pt<20 & (fbrem>0.2 || (fbrem<0.2 & EB & e/p>0.95))
	Bool_t el_fbrem_eOverP_;

	// more var
	Int_t el_class_;
	Int_t el_category_;
	Float_t el_eg_roblooseid_;
	Float_t el_eg_robtightid_;
	Float_t el_eg_looseid_;
	Float_t el_eg_tightid_;
	Float_t el_eg_robhighE_;


	// trigger matching. 2 means matching
	Int_t el8_v1_; // HLT_Ele8_v1
	Int_t el8_v2_; // HLT_Ele8_v2
	Int_t el8_v3_; // HLT_Ele8_v3
	Int_t el8idisojet40_v1_; // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1
	Int_t el8idisojet40_v2_; // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2
	Int_t el8idisojet40_v3_; // HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3
	Int_t el8idiso_v1_; // HLT_Ele8_CaloIdL_CaloIsoVL_v1
	Int_t el8idiso_v2_; // HLT_Ele8_CaloIdL_CaloIsoVL_v2
	Int_t el8idiso_v3_; // HLT_Ele8_CaloIdL_CaloIsoVL_v3
	Int_t el17idiso_v1_; // HLT_Ele17_CaloIdL_CaloIsoVL_v1
	Int_t el17idiso_v2_; // HLT_Ele17_CaloIdL_CaloIsoVL_v2
	Int_t el17idiso_v3_; // HLT_Ele17_CaloIdL_CaloIsoVL_v3
	Int_t el8pho20_v1_; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1
	Int_t el8pho20_v2_; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2
	Int_t el8pho20_v3_; // HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3

};

float deltaPhi( float phi1 , float phi2);


#endif
