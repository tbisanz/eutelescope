/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelHitProjectionCALICE.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "marlin/AIDAProcessor.h"

#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>

// lcio includes <.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;
using namespace marlin;

EUTelHitProjectionCALICE::EUTelHitProjectionCALICE()
  :Processor("EUTelHitProjectionCALICE"),
   _trackCollectionNameInput(), _caliceRawCollectionNameInput() {

  _description = "EUTelHitProjectionCALICE etrapolates tracks onto the CALICE layer and performs a preliminary analysis of the data";
  
  registerInputCollection(LCIO::TRACK,
			   "trackCollectionNameInput", "Input track collection name, for tracks from the tracking processor",
			   _trackCollectionNameInput, std::string ("track"));
	
  registerInputCollection(LCIO::TRACKERDATA,
			   "caliceCollectionNameInput", "Input collection name of the calice raw data",
			   _caliceRawCollectionNameInput, std::string ("calice_raw"));		   
}

void EUTelHitProjectionCALICE::init() {
  geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, 
					     EUTELESCOPE::DUMPGEOROOT);
  bookHistos();
}

void EUTelHitProjectionCALICE::processEvent(LCEvent* event) {
  auto evt = static_cast<EUTelEventImpl*>(event);				
  if( evt->getEventType() == kEORE ) {
    streamlog_out( MESSAGE5 ) << "EORE found: nothing else to do." << std::endl;
    return;
  }

  LCCollection* trackCollection = nullptr;
  try {
    trackCollection = evt->getCollection(_trackCollectionNameInput);
  } catch(DataNotAvailableException& e) {
    streamlog_out( WARNING2 ) << _trackCollectionNameInput 
			      << " collection not available" << std::endl;
    return;
  }

  bool hasMeasHitBeforeLast = false;
  bool hasMeasHitLast = false;
  std::array<double,3> measHitBeforeLast;
  std::array<double,3> measHitLast {0,0,0};
  std::array<double,3> predHitOnCalo {0,0,0};

  //prepare some decoders
  CellIDDecoder<TrackerHit> trackerHitDecoder(EUTELESCOPE::HITENCODING);
  CellIDDecoder<TrackerDataImpl> trackerDataDecoder(EUTELESCOPE::ZSDATADEFAULTENCODING);

  //We only want to process events with a single track, other events are discarded and the single track rejection counter is incremented
  auto noTracks = trackCollection->getNumberOfElements();
  if(noTracks == 1) {
    auto track = static_cast<TrackImpl*>(trackCollection->getElementAt(0));
    auto trackerHits = track->getTrackerHits();
    for(auto& hit: trackerHits) {
      auto sensorID = trackerHitDecoder(hit)["sensorID"];
      //We're interested in sensor 4 and 5 as they are the last telescope sensors and sensor 8 as it is the sensor we're investigating
      if(sensorID == 4 or sensorID == 5 or sensorID == 8){
        auto isSimulated = (trackerHitDecoder(hit)["properties"] & HitProperties::kFittedHit) == HitProperties::kFittedHit;
        std::array<double,3> position { hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] };
	//a non simulated hit is a measured one, we only want tracks with measurements on the next to last and last sensor
        if(!isSimulated) {
          if(sensorID == 4){
            measHitBeforeLast = position;
            hasMeasHitBeforeLast = true;
          } else if(sensorID == 5) {
            measHitLast = position;
            hasMeasHitLast = true;
          }
        } else {
	  //if it is simulated and on sensor 8, it's the track prediction onto our device of interest
          if(sensorID == 8){
            predHitOnCalo = position;
          }
        }
      }   	
    }
  } else {
    singleTrackRejects++;
    return;
  }

  //events without hits on the last two telescope layers are also rejected
  if(!(hasMeasHitBeforeLast && hasMeasHitLast)){
    twoHitRejects++;
    return;
  }

  LCCollection* caliceCollection = nullptr;
  try {
    caliceCollection = evt->getCollection(_caliceRawCollectionNameInput);
  } catch(DataNotAvailableException& e) {
    streamlog_out( WARNING2 ) << _caliceRawCollectionNameInput 
			      << " collection not available" << std::endl;
    return;
  }

  //This is a somewhat unelegant solution to store something like tuples of {x,y,charge} in a dynamic container
  std::vector<std::pair<int,int>> hitCaloCells; //holds {x,y}
  std::vector<int> hitCaloCellsADC; // holds {charge} in the same order

  for(int iHit = 0; iHit < caliceCollection->getNumberOfElements(); iHit++) {
    auto rawData = dynamic_cast<TrackerDataImpl*>(caliceCollection->getElementAt(iHit));
    auto type = static_cast<SparsePixelType>(static_cast<int>(trackerDataDecoder(rawData)["sparsePixelType"]));
    int sensorID = static_cast<int>(trackerDataDecoder(rawData)["sensorID"]);

    //we're only interested in hits on the first layer, this is layer 8 in the august test beam
    if(sensorID == 8) {
      auto sparseData = Utility::getSparseData(rawData, type);
      auto hitPixelVec = sparseData->getPixels();
      for(auto& pixel: hitPixelVec){
        hitCaloCells.emplace_back(std::make_pair(pixel.get().getXCoord(), pixel.get().getYCoord()));
	hitCaloCellsADC.emplace_back(pixel.get().getSignal());
     }
    }
  }

  //In order to process the predicted hits on the calorimeter, we need to transform them into
  //the local coordinate system on that detector
  auto localPosEUTelCoord = std::array<double,3>{0.,0.,0.};
  auto sensorID = 8;
  geo::gGeometry().master2Local(sensorID, predHitOnCalo, localPosEUTelCoord);
 
  //This is a hack to change and flix axese, this we should fix!
  auto localPos = std::array<double,3>{-localPosEUTelCoord[1],-localPosEUTelCoord[0],localPosEUTelCoord[2]};
 
  //Counting the spatially resolved tracks for single track events where the last two detector
  //layers also registered a hit
  allTrackCountHist->fill(localPos[0],localPos[1]);  

  //If we did't detect a hit on the calo layer, increment inefficiency map
  if(hitCaloCells.empty()){
    ineffHist->fill(localPos[0],localPos[1]);
  }

  //If we didn't detect a single hit, but zero or more than one hit we don't resolve the unambiguity
  if(hitCaloCells.size() != 1){
    singleCaloHitRejects++;
    return;
  }

  ADCHist->fill(localPos[0],localPos[1],hitCaloCellsADC.front());  
  ADCCountHist->fill(localPos[0],localPos[1]);  
  auto caloCell = hitCaloCells.front();
  histoMap[caloCell.first][caloCell.second]->fill(measHitLast[0],measHitLast[1]);
}

void EUTelHitProjectionCALICE::bookHistos() {
  ineffHist = AIDAProcessor::histogramFactory(this)-> createHistogram2D("missinghit_det8", 2500, -250, 250, 2500, -250, 250);
  ineffHist->setTitle(";position-x [mm];position-y [mm]");

  ADCHist = AIDAProcessor::histogramFactory(this)-> createHistogram2D("adc_hist", 2500, -250, 250, 2500, -250, 250);
  ADCHist->setTitle(";position-x [mm];position-y [mm]");

  ADCCountHist = AIDAProcessor::histogramFactory(this)-> createHistogram2D("count_hist", 2500, -250, 250, 2500, -250, 250);
  ADCCountHist->setTitle(";position-x [mm];position-y [mm]");

  allTrackCountHist = AIDAProcessor::histogramFactory(this)-> createHistogram2D("all_count_hist", 2500, -250, 250, 2500, -250, 250);
  allTrackCountHist->setTitle(";position-x [mm];position-y [mm]");

  AIDAProcessor::tree(this)->mkdir("individualhitmaps");
  for(size_t ix = 0; ix < 24; ix++){
	for(size_t jx = 0; jx < 24; jx++){
		histoMap[ix][jx] = AIDAProcessor::histogramFactory(this)->
    		  createHistogram2D("individualhitmaps/hitpos_"+std::to_string(ix)+"_"+std::to_string(jx), 180, -18, 18, 80, -8, 8);
		histoMap[ix][jx]->setTitle(";position-x [mm];position-y [mm]");
        }
  }
}

void EUTelHitProjectionCALICE::end()
{
  std::cout << "sing track rejetcs: " << singleTrackRejects << " two hit rejetcs: " << twoHitRejects << " single calo hit rejects: " << singleCaloHitRejects << '\n';
  AIDAProcessor::histogramFactory(this)->divide("ineff_map", *ineffHist, *allTrackCountHist);
  AIDAProcessor::histogramFactory(this)->divide("avg_adc_map", *ADCCountHist, *allTrackCountHist);
}
