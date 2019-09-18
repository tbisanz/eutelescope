/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelCALICE.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"

#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCGenericObject.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <numeric>
#include <utility>

using namespace eutelescope;
using namespace marlin;

EUTelCALICE::EUTelCALICE()
  :Processor("EUTelCALICE"),
   _CALICERawDataCollectionName() {

  _description = "EUTelCALICE processes raw CALICE AHCAL data and provides them in a EUTelescope/LCIO consistent format.";

  registerInputCollection(LCIO::LCGENERICOBJECT,
			  "CALICECollectionNameInput", "Input hit collection name",
			  _CALICERawDataCollectionName, std::string("calice_daq_data"));
}

void EUTelCALICE::init() {
  bookHistos();
}

/** Maps the chip and index of a DAQ channel to the index of the large layer or HBU */
std::pair<int,int> EUTelCALICE::mapChipSiPM(int chip, int sipm){
	auto mod6 = sipm%6;
	auto div6 = sipm/6;

	auto x = 0;
	auto y = 0;

//This is the mapping to an HBU
/* may data
	if(chip==0){
		x = 6+div6;
		y = 6+mod6;
	} else if(chip==1){
		x = 6+div6;
		y = mod6;
	} else if(chip==2){
		x = div6;
		y = 11-mod6;
	} else if(chip==3){
		x=div6;
		y=5-mod6;
	}
*/
	if(chip==0){
		x = 18+div6;
		y = 18+mod6;
	} else if (chip==1) {
		x = 18+div6;
		y = 12+mod6;
	} else if (chip==2) {
		x = 12+div6;
		y = 23-mod6;
	} else if (chip==3) {
		x = 12+div6;
		y = 17-mod6;
	} else if (chip==4) {
		x = 18+div6;
		y = 6+mod6;
	} else if (chip==5) {
		x = 18+div6;
		y = mod6;
	} else if (chip==6) {
		x = 12+div6;
		y = 11-mod6;
	} else if (chip==7) {
		x = 12+div6;
		y = 5-mod6;
	} else if (chip==8) {
		x = 6+div6;
		y = 18+mod6;
	} else if (chip==9) {
		x = 6+div6;
		y = 12+mod6;
	} else if (chip==10) {
		x = div6;
		y = 23-mod6;
	} else if (chip==11) {
		x = div6;
		y = 17-mod6;
	} else if (chip==12) {
		x = 6+div6;
		y = 6+mod6;
	} else if (chip==13) {
		x = 6+div6;
		y = mod6;
	} else if (chip==14) {
		x = div6;
		y = 11-mod6;
	} else if (chip==15) {
		x = div6;
		y = 5-mod6;
	}

	return std::make_pair(x,y);
}

void EUTelCALICE::processEvent(LCEvent* event)
{
  //check event type and if it is the last event
  auto evt = static_cast<EUTelEventImpl*>(event);				
  if( evt->getEventType() == kEORE ) {
    streamlog_out( MESSAGE5 ) << "EORE found: nothing else to do." << std::endl;
    return;
  }

  //opens the CALICE raw data collection into which the DAQ wrote the data
  LCCollection* caliceCollection = nullptr;
  try {
    caliceCollection = evt->getCollection(_CALICERawDataCollectionName);
  } catch(DataNotAvailableException& e) {
    streamlog_out( WARNING2 ) << _CALICERawDataCollectionName 
			      << " collection not available" << std::endl;
    return;
  }

  //This collection is created by this processor, it stores an LCIO hit for every AHCAL hit
  LCCollectionVec* hitCollection = nullptr;
  std::string _hitCollectionName = "calice_direct_hit";
  try {
    hitCollection = static_cast<LCCollectionVec *>(
        event->getCollection(_hitCollectionName));
  } catch(...) {
    hitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
  }


  //This collection is created by this processor, it stores an LCIO raw hit for every AHCAL hit
  LCCollectionVec* rawCollection = nullptr;
  std::string _rawCollectionName = "calice_raw";
  try {
    rawCollection = static_cast<LCCollectionVec *>(
        event->getCollection(_rawCollectionName));
  } catch(...) {
    rawCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  }

  CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING, hitCollection);
  CellIDEncoder<TrackerDataImpl> idRawEncoder(EUTELESCOPE::ZSDATADEFAULTENCODING, rawCollection);

  //TODO: these are the AHCAL detector layers, hard coded as it is right now!
  std::map<int, std::vector<float>> rawData = {	{8, std::vector<float>()},
						{9, std::vector<float>()},
						{10,std::vector<float>()},
						{11,std::vector<float>()},
						{13,std::vector<float>()} };

    //[START] loop over hits
    for(int iEntry = 0; iEntry < caliceCollection->getNumberOfElements(); ++iEntry) {  
	auto obj = static_cast<EVENT::LCGenericObject*>(caliceCollection->getElementAt(iEntry));
	auto chipIDRaw = obj->getIntVal(3);

	auto moduleNumber = chipIDRaw >> 8;
	auto chipID = chipIDRaw & 0x00FF;

	std::vector<int> dacs;
	for(size_t ix = 0; ix < 36; ix++){
		/*There are 5 generic entries, then 36 TDC and 36 ACD values,
		 *we're interested in the ADCs */
		dacs.emplace_back( obj->getIntVal(5+36+ix) );
	}
	auto average = std::accumulate(dacs.begin(), dacs.end(), 0)/dacs.size();

	for(auto ix = dacs.begin(); ix != dacs.end(); ix++){
		//This is the hit bit and indicates if an ADC reading corresponds to a hit
		if(*ix & 0x1000) {
			auto hit = mapChipSiPM(chipID, ix-dacs.begin());
			_hitsVsTracksHistos[moduleNumber]->fill(hit.first,hit.second);
			auto lchit = new TrackerHitImpl;
			std::array<double,3> pos = { (hit.first - 6 + 0.5)*31.5, (hit.second - 6 + 0.5)*31.5, 0.0 };
			lchit->setPosition(pos.data());

			idHitEncoder["sensorID"] = moduleNumber;
			idHitEncoder["properties"] = 0;
    			idHitEncoder.setCellID(lchit);
    			hitCollection->push_back(lchit);

			/**This is a makeshift hack and exploits the fact, that EUTelGenericSparsePixel has
                          * four entries per hit: X, Y, charge, time */
			rawData[moduleNumber].emplace_back(hit.first);
			rawData[moduleNumber].emplace_back(hit.second);
			rawData[moduleNumber].emplace_back(*ix-average);
			rawData[moduleNumber].emplace_back(0);
		}
	}
    }//[END] loop over hits

    //For each module, we write out the hit calo cells as raw hits
    for(auto& module: rawData) {
	auto moduleID = module.first;
	auto data = new TrackerDataImpl;
	data->setChargeValues(module.second);
	idRawEncoder["sensorID"] = moduleID;
	idRawEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;
	idRawEncoder.setCellID(data);
	rawCollection->push_back(data);
    } 
  
  //This checks if the collection already exists and adds it to the LCIO file if it does not yet
  try {
    event->getCollection(_hitCollectionName);
  } catch(...) {
    event->addCollection(hitCollection, _hitCollectionName);
  }
  try {
    event->getCollection(_rawCollectionName);
  } catch(...) {
    event->addCollection(rawCollection, _rawCollectionName);
  }
}

void EUTelCALICE::bookHistos() {
	std::string histName_hitsVsTracksAll = "/all";
	hitsVsTracksAll = AIDAProcessor::histogramFactory(this)->
	  createHistogram2D(histName_hitsVsTracksAll, 4, -0.5, 3.5, 25, -2.5, 22.5);
	hitsVsTracksAll->setTitle(";Tracks;Hits on all detectors");
	nTracksHist = AIDAProcessor::histogramFactory(this)->
	  createHistogram1D("nTracks", 5, -0.5, 4.5);

	std::vector<int> HBUmodules = {8,9,10,11,13};
	for(auto module: HBUmodules) {
		//create folder for current detector
		std::string basePath = "detector_" + std::to_string(module);
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		std::string histName_hitsVsTracks = basePath+"/hitmap_module"+ std::to_string(module);
		AIDA::IHistogram2D* hitsVsTracks = AIDAProcessor::histogramFactory(this)->
		  createHistogram2D(histName_hitsVsTracks, 24, -0.5, 23.5, 24, -0.5, 23.5);
		hitsVsTracks->setTitle("hits vs tracks");
		_hitsVsTracksHistos.insert(std::make_pair(module, hitsVsTracks));
	}
}

void EUTelCALICE::end() {
  streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}
