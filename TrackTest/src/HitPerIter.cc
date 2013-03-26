// -*- C++ -*-
//
// Package:    HitPerIter
// Class:      HitPerIter
// 
/**\class HitPerIter HitPerIter.cc Test/TrackTest/src/HitPerIter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giuseppe Cerati,28 S-012,+41227678302,
//         Created:  Thu Dec 13 12:06:32 CET 2012
// $Id: HitPerIter.cc,v 1.3 2013/03/14 12:49:26 cerati Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/ContainerMask.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"

//
// class declaration
//

class HitPerIter : public edm::EDAnalyzer {
   public:
      explicit HitPerIter(const edm::ParameterSet&);
      ~HitPerIter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      template<class T, class U> void makeFillHisto1D(const char* name,const char* title,int nbins,U minx,U maxx,U value);

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag trackTags_; //used to select what tracks to read from configuration file
      std::string builderName;
      const TransientTrackingRecHitBuilder* builder;
      bool saveHists;
      bool printInfo;
};

//
// constructors and destructor
//
HitPerIter::HitPerIter(const edm::ParameterSet& iConfig)
  :trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
   builderName(iConfig.getParameter<std::string>("TTRHBuilder")), 
   saveHists(iConfig.getParameter<bool>("saveHists")), 
   printInfo(iConfig.getParameter<bool>("printInfo")) 
{
  if (saveHists) {
    edm::Service<TFileService> fs;
    fs->make<TDirectory>("dummy","dummy");//do this so that the file is properly created
  }
}


HitPerIter::~HitPerIter() {}

template<class T, class U> void HitPerIter::makeFillHisto1D(const char* name,const char* title,int nbins,U minx,U maxx,U value) {
  edm::Service<TFileService> fs;
  T* h = 0;
  try {
    h = fs->getObject<T>(name);
  } catch (cms::Exception e) {
    if (e.category()=="ObjectNotFound") {
      //std::cout << name <<" "<< title <<" "<< nbins <<" "<< minx <<" "<< maxx << std::endl;
      h = fs->make<T>(name, title, nbins, minx, maxx);
    } else throw e;
  }
  h->Fill(value);
}

//
// member functions
//

// ------------ method called for each event  ------------
void HitPerIter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;

  //tracks
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
  int totTk = 0, totHp = 0;
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    totTk++;
    if (saveHists) {
      //track basic plots
      makeFillHisto1D<TH1F,int>("algo","algo",10,0,10,itTrack->algo()-4);
      if (itTrack->quality(TrackBase::highPurity)) {
	totHp++;
	makeFillHisto1D<TH1F,int>("algo_hp","algo_hp",10,0,10,itTrack->algo()-4);
      }
    }
  }
  makeFillHisto1D<TH1F,int>("ntk","ntk",100,0,5000,totTk);
  makeFillHisto1D<TH1F,int>("nhp","nhp",100,0,5000,totHp);

  //hits
  edm::Handle<SiPixelRecHitCollection> pixelHits;
  iEvent.getByLabel("siPixelRecHits", pixelHits);

  edm::Handle<SiStripMatchedRecHit2DCollection> matchedHits;
  iEvent.getByLabel("siStripMatchedRecHits","matchedRecHit", matchedHits);

  edm::Handle<SiStripRecHit2DCollection> rphiHits;
  iEvent.getByLabel("siStripMatchedRecHits","rphiRecHit", rphiHits);

  edm::Handle<SiStripRecHit2DCollection> rphiuHits;
  iEvent.getByLabel("siStripMatchedRecHits","rphiRecHitUnmatched", rphiuHits);

  edm::Handle<SiStripRecHit2DCollection> stereoHits;
  iEvent.getByLabel("siStripMatchedRecHits","stereoRecHit", stereoHits);

  edm::Handle<SiStripRecHit2DCollection> stereouHits;
  iEvent.getByLabel("siStripMatchedRecHits","stereoRecHitUnmatched", stereouHits);

  string labels[] = {"","lowPtTripletStepClusters","pixelPairStepClusters","detachedTripletStepClusters","mixedTripletStepClusters","pixelLessStepClusters","tobTecStepClusters"};

  for (int iter=0;iter<7;iter++) {

    edm::Handle<edm::ContainerMask<edmNew::DetSetVector<SiStripCluster> > > stripClusterMask;
    if (labels[iter]!="") iEvent.getByLabel(labels[iter],stripClusterMask);

    edm::Handle<edm::ContainerMask<edmNew::DetSetVector<SiPixelCluster> > > pixelClusterMask;
    if (labels[iter]!="") iEvent.getByLabel(labels[iter],pixelClusterMask);
    
    int mtchTot = 0;
    int mtchUnm = 0;
    for (auto it = matchedHits->begin(); it!=matchedHits->end(); it++ ) {
      DetId hitId = it->detId();
      for (auto hit = it->begin(); hit!=it->end(); hit++ ) {
	if (labels[iter]=="" || (stripClusterMask->mask(hit->monoClusterRef().key())==0 || stripClusterMask->mask(hit->stereoClusterRef().key())==0) ) {//mask means skip
	  makeFillHisto1D<TH1F,int>("mtch","mtch",8,0,8,iter);
	  if (hitId.subdetId() == StripSubdetector::TIB){
	    if (TIBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("mtch_tib1","mtch_tib1",8,0,8,iter);
	    else if (TIBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("mtch_tib2","mtch_tib2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TOB){
	    if (TOBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("mtch_tob1","mtch_tob1",8,0,8,iter);
	    else if (TOBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("mtch_tob2","mtch_tob2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TID) {
	    if (TIDDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("mtch_tid1","mtch_tid1",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("mtch_tid2","mtch_tid2",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("mtch_tid3","mtch_tid3",8,0,8,iter);	      
	  } else if (hitId.subdetId() == StripSubdetector::TEC) {
	    if (TECDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("mtch_tec1","mtch_tec1",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("mtch_tec2","mtch_tec2",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("mtch_tec3","mtch_tec3",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==4) makeFillHisto1D<TH1F,int>("mtch_tec4","mtch_tec4",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==5) makeFillHisto1D<TH1F,int>("mtch_tec5","mtch_tec5",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==6) makeFillHisto1D<TH1F,int>("mtch_tec6","mtch_tec6",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==7) makeFillHisto1D<TH1F,int>("mtch_tec7","mtch_tec7",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==8) makeFillHisto1D<TH1F,int>("mtch_tec8","mtch_tec8",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==9) makeFillHisto1D<TH1F,int>("mtch_tec9","mtch_tec9",8,0,8,iter);	      
	  }
	  mtchUnm++;
	} 
	mtchTot++;
      }
    }
    //std::cout << "strip matched hits:" << mtchTot  << " " << mtchUnm << std::endl;
    makeFillHisto1D<TH1F,int>(Form("nmtch_tot_it%i",iter),Form("nmtch_tot_it%i",iter),100,0,100000,mtchTot);
    makeFillHisto1D<TH1F,int>(Form("nmtch_unm_it%i",iter),Form("nmtch_unm_it%i",iter),100,0,100000,mtchUnm);
  
    int rphiTot = 0;
    int rphiUnm = 0;
    for (auto it = rphiHits->begin(); it!=rphiHits->end(); it++ ) {
      DetId hitId = it->detId();
      for (auto hit = it->begin(); hit!=it->end(); hit++ ) {
	if (labels[iter]=="" || stripClusterMask->mask(hit->cluster().key())==0) {
	  makeFillHisto1D<TH1F,int>("rphi","rphi",8,0,8,iter);
	  if (hitId.subdetId() == StripSubdetector::TIB){
	    if (TIBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("rphi_tib1","rphi_tib1",8,0,8,iter);
	    else if (TIBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("rphi_tib2","rphi_tib2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TOB){
	    if (TOBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("rphi_tob1","rphi_tob1",8,0,8,iter);
	    else if (TOBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("rphi_tob2","rphi_tob2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TID) {
	    if (TIDDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("rphi_tid1","rphi_tid1",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("rphi_tid2","rphi_tid2",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("rphi_tid3","rphi_tid3",8,0,8,iter);	      
	  } else if (hitId.subdetId() == StripSubdetector::TEC) {
	    if (TECDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("rphi_tec1","rphi_tec1",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("rphi_tec2","rphi_tec2",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("rphi_tec3","rphi_tec3",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==4) makeFillHisto1D<TH1F,int>("rphi_tec4","rphi_tec4",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==5) makeFillHisto1D<TH1F,int>("rphi_tec5","rphi_tec5",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==6) makeFillHisto1D<TH1F,int>("rphi_tec6","rphi_tec6",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==7) makeFillHisto1D<TH1F,int>("rphi_tec7","rphi_tec7",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==8) makeFillHisto1D<TH1F,int>("rphi_tec8","rphi_tec8",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==9) makeFillHisto1D<TH1F,int>("rphi_tec9","rphi_tec9",8,0,8,iter);	      
	  }
	  rphiUnm++;
	}
	rphiTot++;
      }
    }
    //std::cout << "strip rphi hits:" << rphiTot  << " " << rphiUnm << std::endl;
    makeFillHisto1D<TH1F,int>(Form("nrphi_tot_it%i",iter),Form("nrphi_tot_it%i",iter),100,0,100000,rphiTot);
    makeFillHisto1D<TH1F,int>(Form("nrphi_unm_it%i",iter),Form("nrphi_unm_it%i",iter),100,0,100000,rphiUnm);
    
    int rphiuTot = 0;
    int rphiuUnm = 0;
    for (auto it = rphiuHits->begin(); it!=rphiuHits->end(); it++ ) {
      DetId hitId = it->detId();
      for (auto hit = it->begin(); hit!=it->end(); hit++ ) {
	if (labels[iter]=="" || stripClusterMask->mask(hit->cluster().key())==0) {
	  makeFillHisto1D<TH1F,int>("rphiu","rphiu",8,0,8,iter);
	  if (hitId.subdetId() == StripSubdetector::TIB){
	    if (TIBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("rphiu_tib1","rphiu_tib1",8,0,8,iter);
	    else if (TIBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("rphiu_tib2","rphiu_tib2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TOB){
	    if (TOBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("rphiu_tob1","rphiu_tob1",8,0,8,iter);
	    else if (TOBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("rphiu_tob2","rphiu_tob2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TID) {
	    if (TIDDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("rphiu_tid1","rphiu_tid1",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("rphiu_tid2","rphiu_tid2",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("rphiu_tid3","rphiu_tid3",8,0,8,iter);	      
	  } else if (hitId.subdetId() == StripSubdetector::TEC) {
	    if (TECDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("rphiu_tec1","rphiu_tec1",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("rphiu_tec2","rphiu_tec2",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("rphiu_tec3","rphiu_tec3",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==4) makeFillHisto1D<TH1F,int>("rphiu_tec4","rphiu_tec4",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==5) makeFillHisto1D<TH1F,int>("rphiu_tec5","rphiu_tec5",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==6) makeFillHisto1D<TH1F,int>("rphiu_tec6","rphiu_tec6",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==7) makeFillHisto1D<TH1F,int>("rphiu_tec7","rphiu_tec7",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==8) makeFillHisto1D<TH1F,int>("rphiu_tec8","rphiu_tec8",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==9) makeFillHisto1D<TH1F,int>("rphiu_tec9","rphiu_tec9",8,0,8,iter);	      
	  }
	  rphiuUnm++;
	}
	rphiuTot++;
      }
    }
    //std::cout << "strip rphiu hits:" << rphiuTot  << " " << rphiuUnm << std::endl;
    makeFillHisto1D<TH1F,int>(Form("nrphiu_tot_it%i",iter),Form("nrphiu_tot_it%i",iter),100,0,100000,rphiuTot);
    makeFillHisto1D<TH1F,int>(Form("nrphiu_unm_it%i",iter),Form("nrphiu_unm_it%i",iter),100,0,100000,rphiuUnm);

    int stereoTot = 0;
    int stereoUnm = 0;
    for (auto it = stereoHits->begin(); it!=stereoHits->end(); it++ ) {
      DetId hitId = it->detId();
      for (auto hit = it->begin(); hit!=it->end(); hit++ ) {
	if (labels[iter]=="" || stripClusterMask->mask(hit->cluster().key())==0) {
	  makeFillHisto1D<TH1F,int>("stereo","stereo",8,0,8,iter);
	  if (hitId.subdetId() == StripSubdetector::TIB){
	    if (TIBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("stereo_tib1","stereo_tib1",8,0,8,iter);
	    else if (TIBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("stereo_tib2","stereo_tib2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TOB){
	    if (TOBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("stereo_tob1","stereo_tob1",8,0,8,iter);
	    else if (TOBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("stereo_tob2","stereo_tob2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TID) {
	    if (TIDDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("stereo_tid1","stereo_tid1",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("stereo_tid2","stereo_tid2",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("stereo_tid3","stereo_tid3",8,0,8,iter);	      
	  } else if (hitId.subdetId() == StripSubdetector::TEC) {
	    if (TECDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("stereo_tec1","stereo_tec1",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("stereo_tec2","stereo_tec2",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("stereo_tec3","stereo_tec3",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==4) makeFillHisto1D<TH1F,int>("stereo_tec4","stereo_tec4",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==5) makeFillHisto1D<TH1F,int>("stereo_tec5","stereo_tec5",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==6) makeFillHisto1D<TH1F,int>("stereo_tec6","stereo_tec6",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==7) makeFillHisto1D<TH1F,int>("stereo_tec7","stereo_tec7",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==8) makeFillHisto1D<TH1F,int>("stereo_tec8","stereo_tec8",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==9) makeFillHisto1D<TH1F,int>("stereo_tec9","stereo_tec9",8,0,8,iter);	      
	  }
	  stereoUnm++;
	}
	stereoTot++;
      }
    }
    //std::cout << "strip stereo hits:" << stereoTot  << " " << stereoUnm << std::endl;
    makeFillHisto1D<TH1F,int>(Form("nstereo_tot_it%i",iter),Form("nstereo_tot_it%i",iter),100,0,100000,stereoTot);
    makeFillHisto1D<TH1F,int>(Form("nstereo_unm_it%i",iter),Form("nstereo_unm_it%i",iter),100,0,100000,stereoUnm);
    
    int stereouTot = 0;
    int stereouUnm = 0;
    for (auto it = stereouHits->begin(); it!=stereouHits->end(); it++ ) {
      DetId hitId = it->detId();
      for (auto hit = it->begin(); hit!=it->end(); hit++ ) {
	if (labels[iter]=="" || stripClusterMask->mask(hit->cluster().key())==0) {
	  makeFillHisto1D<TH1F,int>("stereou","stereou",8,0,8,iter);
	  if (hitId.subdetId() == StripSubdetector::TIB){
	    if (TIBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("stereou_tib1","stereou_tib1",8,0,8,iter);
	    else if (TIBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("stereou_tib2","stereou_tib2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TOB){
	    if (TOBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("stereou_tob1","stereou_tob1",8,0,8,iter);
	    else if (TOBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("stereou_tob2","stereou_tob2",8,0,8,iter);
	  } else if (hitId.subdetId() == StripSubdetector::TID) {
	    if (TIDDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("stereou_tid1","stereou_tid1",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("stereou_tid2","stereou_tid2",8,0,8,iter);
	    else if (TIDDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("stereou_tid3","stereou_tid3",8,0,8,iter);	      
	  } else if (hitId.subdetId() == StripSubdetector::TEC) {
	    if (TECDetId(hitId).wheel()==1) makeFillHisto1D<TH1F,int>("stereou_tec1","stereou_tec1",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==2) makeFillHisto1D<TH1F,int>("stereou_tec2","stereou_tec2",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==3) makeFillHisto1D<TH1F,int>("stereou_tec3","stereou_tec3",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==4) makeFillHisto1D<TH1F,int>("stereou_tec4","stereou_tec4",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==5) makeFillHisto1D<TH1F,int>("stereou_tec5","stereou_tec5",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==6) makeFillHisto1D<TH1F,int>("stereou_tec6","stereou_tec6",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==7) makeFillHisto1D<TH1F,int>("stereou_tec7","stereou_tec7",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==8) makeFillHisto1D<TH1F,int>("stereou_tec8","stereou_tec8",8,0,8,iter);
	    else if (TECDetId(hitId).wheel()==9) makeFillHisto1D<TH1F,int>("stereou_tec9","stereou_tec9",8,0,8,iter);	      
	  }
	  stereouUnm++;
	}
	stereouTot++;
      }
    }
    //std::cout << "strip stereou hits:" << stereouTot  << " " << stereouUnm << std::endl;
    makeFillHisto1D<TH1F,int>(Form("nstereou_tot_it%i",iter),Form("nstereou_tot_it%i",iter),100,0,100000,stereouTot);
    makeFillHisto1D<TH1F,int>(Form("nstereou_unm_it%i",iter),Form("nstereou_unm_it%i",iter),100,0,100000,stereouUnm);

    int pixelTot = 0;
    int pixelUnm = 0;
    for (auto it = pixelHits->begin(); it!=pixelHits->end(); it++ ) {
      DetId hitId = it->detId();
      for (auto hit = it->begin(); hit!=it->end(); hit++ ) {
	//TransientTrackingRecHit::RecHitPointer ttrh = builder->build(&*hit);
	//std::cout << ttrh->globalPosition() << std::endl;
	if (labels[iter]=="" || pixelClusterMask->mask(hit->cluster().key())==0) {
	  makeFillHisto1D<TH1F,int>("pixel","pixel",8,0,8,iter);
	  if (hitId.subdetId() == (int) PixelSubdetector::PixelBarrel ) {
	    if (PXBDetId(hitId).layer()==1) makeFillHisto1D<TH1F,int>("pixel_pxb1","pixel_pxb1",8,0,8,iter);
	    else if (PXBDetId(hitId).layer()==2) makeFillHisto1D<TH1F,int>("pixel_pxb2","pixel_pxb2",8,0,8,iter);
	    else if (PXBDetId(hitId).layer()==3) makeFillHisto1D<TH1F,int>("pixel_pxb3","pixel_pxb3",8,0,8,iter);
	  } else if (hitId.subdetId() == (int) PixelSubdetector::PixelEndcap ) {
	    if (PXFDetId(hitId).disk()==1) makeFillHisto1D<TH1F,int>("pixel_pxf1","pixel_pxf1",8,0,8,iter);
	    else if (PXFDetId(hitId).disk()==2) makeFillHisto1D<TH1F,int>("pixel_pxf2","pixel_pxf2",8,0,8,iter);
	  }
	  pixelUnm++;
	}
	pixelTot++;
      }
    }
    //std::cout << "pixel hits:" << pixelTot  << " " << pixelUnm << std::endl;
    makeFillHisto1D<TH1F,int>(Form("npixel_tot_it%i",iter),Form("npixel_tot_it%i",iter),100,0,10000,pixelTot);
    makeFillHisto1D<TH1F,int>(Form("npixel_unm_it%i",iter),Form("npixel_unm_it%i",iter),100,0,10000,pixelUnm);
    

  }

}


// ------------ method called once each job just before starting event loop  ------------
void HitPerIter::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void HitPerIter::endJob() {}

// ------------ method called when starting to processes a run  ------------
void HitPerIter::beginRun(edm::Run const& run, edm::EventSetup const& setup) {
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName,theBuilder);
  builder=theBuilder.product();
}

// ------------ method called when ending the processing of a run  ------------
void HitPerIter::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void HitPerIter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void HitPerIter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HitPerIter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

 //Specify that only 'tracks' is allowed
 //To use, remove the default given above and uncomment below
 //ParameterSetDescription desc;
 //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
 //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitPerIter);
