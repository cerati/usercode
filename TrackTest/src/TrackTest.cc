// -*- C++ -*-
//
// Package:    TrackTest
// Class:      TrackTest
// 
/**\class TrackTest TrackTest.cc Test/TrackTest/src/TrackTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giuseppe Cerati,28 S-012,+41227678302,
//         Created:  Thu Dec 13 12:06:32 CET 2012
// $Id: TrackTest.cc,v 1.1 2012/12/13 17:11:28 cerati Exp $
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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//
// class declaration
//

class TrackTest : public edm::EDAnalyzer {
   public:
      explicit TrackTest(const edm::ParameterSet&);
      ~TrackTest();

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
TrackTest::TrackTest(const edm::ParameterSet& iConfig)
  :trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
   builderName(iConfig.getParameter<std::string>("TTRHBuilder")), 
   saveHists(iConfig.getParameter<bool>("saveHists")), 
   printInfo(iConfig.getParameter<bool>("printInfo")) 
{
  edm::Service<TFileService> fs;
  fs->make<TDirectory>("dummy","dummy");//do this so that the file is properly created
}


TrackTest::~TrackTest() {}

template<class T, class U> void TrackTest::makeFillHisto1D(const char* name,const char* title,int nbins,U minx,U maxx,U value) {
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
void TrackTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot", beamSpot);
  
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
  int ntrk = 0 ;
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {

    if (saveHists) {
      //track basic plots
      makeFillHisto1D<TH1F,int>("algo","algo",20,0,20,itTrack->algo());
      makeFillHisto1D<TH1F,int>("isHighPurity","isHighPurity",3,0,3,itTrack->quality(TrackBase::highPurity));      
      makeFillHisto1D<TH1F,float>("pT","p_{T}",100,0,100,itTrack->pt());
      makeFillHisto1D<TH1F,float>("eta","#eta",60,-3,3,itTrack->eta());
      makeFillHisto1D<TH1F,float>("phi","#phi",64,-M_PI,M_PI,itTrack->phi());
      makeFillHisto1D<TH1F,float>("dxy","dxy",100,-50,50,itTrack->dxy(beamSpot->position()));
      makeFillHisto1D<TH1F,float>("dz","dz",100,-50,50,itTrack->dz(beamSpot->position()));      
      makeFillHisto1D<TH1F,float>("pT_err","#Delta p_{T}",100,0,10,itTrack->ptError());
      makeFillHisto1D<TH1F,float>("eta_err","#Delta #eta",100,0,0.1,itTrack->etaError());
      makeFillHisto1D<TH1F,float>("phi_err","#Delta #phi",100,0,0.1,itTrack->phiError());
      makeFillHisto1D<TH1F,float>("dxy_err","#Delta dxy",100,0,0.5,itTrack->dxyError());
      makeFillHisto1D<TH1F,float>("dz_err","#Delta dz",100,0,0.5,itTrack->dzError());      
      makeFillHisto1D<TH1F,float>("pT_sig","#Delta p_{T}/p_{T}",100,0,1,itTrack->ptError()/itTrack->pt());
      makeFillHisto1D<TH1F,float>("eta_sig","#Delta #eta/#eta",100,0,1,itTrack->etaError()/fabs(itTrack->eta()));
      makeFillHisto1D<TH1F,float>("phi_sig","#Delta #phi/#phi",100,0,1,itTrack->phiError()/fabs(itTrack->phi()));
      makeFillHisto1D<TH1F,float>("dxy_sig","#sigma dxy",100,0,50,fabs(itTrack->dxy(beamSpot->position()))/itTrack->dxyError());
      makeFillHisto1D<TH1F,float>("dz_sig","#sigma dz",100,0,50,fabs(itTrack->dz(beamSpot->position()))/itTrack->dzError());
      makeFillHisto1D<TH1F,int>("nhits","nhits",50,0,50,itTrack->found());
      //hit level plots
      HitPattern hp = itTrack->hitPattern();
      makeFillHisto1D<TH1F,int>("nlayer_pixel","nlayer_pixel",10,0,10,hp.pixelLayersWithMeasurement());
      makeFillHisto1D<TH1F,int>("nlayer_strip","nlayer_strip",30,0,30,hp.stripLayersWithMeasurement());
      makeFillHisto1D<TH1F,int>("nlayer_strip_3d","nlayer_strip_3d",10,0,10,hp.numberOfValidStripLayersWithMonoAndStereo());
      makeFillHisto1D<TH1F,int>("nlayer_tot_3d","nlayer_tot_3d",15,0,15,hp.pixelLayersWithMeasurement()+hp.numberOfValidStripLayersWithMonoAndStereo());
    }
    
    if (printInfo) {
      cout << "Track #" << ntrk << " with q=" << itTrack->charge() 
	   << ", pT=" << itTrack->pt() << " GeV, eta=" << itTrack->eta() 
	   << ", Nhits=" << itTrack->recHitsSize() 
	   << ", algo=" << itTrack->algoName(itTrack->algo()).c_str() << endl;
      int nhit = 0;
      for (trackingRecHit_iterator i=itTrack->recHitsBegin(); i!=itTrack->recHitsEnd(); i++){
	cout << "hit #" << nhit;
	TransientTrackingRecHit::RecHitPointer hit=builder->build(&**i );
	DetId hitId = hit->geographicalId();
	if(hitId.det() == DetId::Tracker) {
	  if (hitId.subdetId() == StripSubdetector::TIB )  
	    cout << " - TIB " << TIBDetId(hitId).layer();
	  else if (hitId.subdetId() == StripSubdetector::TOB ) 
	    cout << " - TOB " << TOBDetId(hitId).layer();
	  else if (hitId.subdetId() == StripSubdetector::TEC ) 
	    cout << " - TEC " << TECDetId(hitId).wheel();
	  else if (hitId.subdetId() == StripSubdetector::TID ) 
	    cout << " - TID " << TIDDetId(hitId).wheel();
	  else if (hitId.subdetId() == StripSubdetector::TID ) 
	    cout << " - TID " << TIDDetId(hitId).wheel();
	  else if (hitId.subdetId() == (int) PixelSubdetector::PixelBarrel ) 
	    cout << " - PixBar " << PXBDetId(hitId).layer();
	  else if (hitId.subdetId() == (int) PixelSubdetector::PixelEndcap )
	    cout << " - PixFwd " << PXFDetId(hitId).disk();
	  else 
	    cout << " UNKNOWN TRACKER HIT TYPE ";
	}
	if (hit->isValid()) cout << " - globalPos =" << hit->globalPosition() << endl;
	else cout << " - invalid hit" << endl;
	nhit++;
      }    
      cout << endl;
    }
    ntrk++;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void TrackTest::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void TrackTest::endJob() {}

// ------------ method called when starting to processes a run  ------------
void TrackTest::beginRun(edm::Run const& run, edm::EventSetup const& setup) {
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName,theBuilder);
  builder=theBuilder.product();
}

// ------------ method called when ending the processing of a run  ------------
void TrackTest::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void TrackTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void TrackTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(TrackTest);
