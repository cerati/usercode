// -*- C++ -*-
//
// Package:    SeedAndTrackAnalyzer
// Class:      SeedAndTrackAnalyzer
// 
/**\class SeedAndTrackAnalyzer SeedAndTrackAnalyzer.cc Demo/SeedAndTrackAnalyzer/src/SeedAndTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giuseppe Cerati,28 S-012,+41227678302,
//         Created:  Tue Apr 16 13:22:56 CEST 2013
// $Id$
//
//

// skelethon generated with command: mkedanlzr -track SeedAndTrackAnalyzer

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

//
// class declaration
//

class SeedAndTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SeedAndTrackAnalyzer(const edm::ParameterSet&);
      ~SeedAndTrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag seedTag_; //used to select what seeds to read from configuration file
      edm::InputTag trackTag_; //used to select what tracks to read from configuration file
      std::string builderName_;
      edm::ESHandle<MagneticField> theMF;
      edm::ESHandle<TransientTrackingRecHitBuilder> theTTRHBuilder;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SeedAndTrackAnalyzer::SeedAndTrackAnalyzer(const edm::ParameterSet& iConfig):
  seedTag_(iConfig.getUntrackedParameter<edm::InputTag>("seeds")),
  trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
  builderName_(iConfig.getParameter<std::string>("TTRHBuilder"))
{

}


SeedAndTrackAnalyzer::~SeedAndTrackAnalyzer() {}


//
// member functions
//

// ------------ method called for each event  ------------
void
SeedAndTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  Handle<TrajectorySeedCollection> seeds;
  iEvent.getByLabel(seedTag_,seeds);
  for(TrajectorySeedCollection::const_iterator itSeed = seeds->begin(); itSeed != seeds->end(); ++itSeed) {
    //find seed class at: http://cmssdt.cern.ch/SDT/lxr/source/DataFormats/TrajectorySeed/interface/TrajectorySeed.h?v=CMSSW_6_2_0_pre2
    int nHits = itSeed->nHits();  
    if (nHits==3) {
      TransientTrackingRecHit::RecHitPointer recHit0 = theTTRHBuilder->build(&*(itSeed->recHits().first));
      TransientTrackingRecHit::RecHitPointer recHit1 = theTTRHBuilder->build(&*(itSeed->recHits().first+1));
      TransientTrackingRecHit::RecHitPointer recHit2 = theTTRHBuilder->build(&*(itSeed->recHits().first+2));

      TrajectoryStateOnSurface state = trajectoryStateTransform::transientState( itSeed->startingState(), recHit2->surface(), theMF.product());

      cout << "seed number: " <<  itSeed-seeds->begin()
	   << " pt=" << state.globalMomentum().perp()
	   << " hit0=" << recHit0->globalPosition()
	   << " hit1=" << recHit1->globalPosition()
	   << " hit2=" << recHit2->globalPosition() << endl;
    }

  }

  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTag_,tracks);
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    //find track class at: http://cmssdt.cern.ch/SDT/lxr/source/DataFormats/TrackReco/interface/Track.h?v=CMSSW_6_2_0_pre2

    //use only highPurity tracks
    if (itTrack->quality(TrackBase::highPurity)==0) continue;

    int charge = itTrack->charge();  
    float pt = itTrack->pt();  
    float eta = itTrack->eta();  
    float phi = itTrack->phi();  
    int nHits = itTrack->numberOfValidHits();  
    edm::RefToBase<TrajectorySeed> trackSeed = itTrack->seedRef();
    if (trackSeed.isNull()==0) {
      TransientTrackingRecHit::RecHitPointer recHit = theTTRHBuilder->build(&*(trackSeed->recHits().second-1));
      TrajectoryStateOnSurface state = trajectoryStateTransform::transientState( trackSeed->startingState(), recHit->surface(), theMF.product());
      cout << "track q/pt/eta/phi/nhits: " << charge << " / " << pt << " / " << eta << " / " << phi << " / " << nHits;
      cout << " --- seed number: " << trackSeed.key() << " pt=" << state.globalMomentum().perp() << endl;
    }
  }
  
}


// ------------ method called once each job just before starting event loop  ------------
void SeedAndTrackAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SeedAndTrackAnalyzer::endJob() {}

// ------------ method called when starting to processes a run  ------------
void SeedAndTrackAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
  iSetup.get<IdealMagneticFieldRecord>().get(theMF);  
  iSetup.get<TransientRecHitRecord>().get(builderName_,theTTRHBuilder);
}

// ------------ method called when ending the processing of a run  ------------
//void SeedAndTrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
//void SeedAndTrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

// ------------ method called when ending the processing of a luminosity block  ------------
//void SeedAndTrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SeedAndTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SeedAndTrackAnalyzer);
