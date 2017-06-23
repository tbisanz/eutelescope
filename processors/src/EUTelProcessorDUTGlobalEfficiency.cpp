//Writen by Alexander Morton
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorDUTGlobalEfficiency.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;

EUTelProcessorDUTGlobalEfficiency::EUTelProcessorDUTGlobalEfficiency():
Processor("EUTelProcessorDUTGlobalEfficiency"),
_trackCollectionName(), 
_hitCollectionName(),
_totalMatchedTracks(0)

{
		_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the EUTelGeometryClass";

		registerInputCollection(LCIO::TRACK, "trackCollectionName", "Local input hit collection name", _trackCollectionName, std::string("local_hit"));

		registerInputCollection(LCIO::TRACKERHIT, "hitCollectionName", "Local input hit collection name", _hitCollectionName, std::string("local_hit"));

		//registerOptionalParameter("Undo Alignment (boolean)", "Set to true to undo the alignment instead", _undoAlignment, bool(false));
}

void EUTelProcessorDUTGlobalEfficiency::init() {
		geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

	_DUTIDs.push_back(21);
	_DUTIDs.push_back(22);
	_DUTIDs.push_back(23);
		
    for(auto& id: _DUTIDs) {
		 _DUTHits[id] = std::vector<const double*>();
	}
    for(auto& id: _DUTIDs) {
		 _DUTMatchedTrackes[id] = 0; 
	}

}

void EUTelProcessorDUTGlobalEfficiency::processRunHeader(LCRunHeader* rdr)
{
	return;
}

void EUTelProcessorDUTGlobalEfficiency::processEvent(LCEvent* event)
{
		//Check the event type and if it is the last event.
		EUTelEventImpl* evt	= static_cast<EUTelEventImpl*>(event);				
		if( evt->getEventType() == kEORE )
		{
				streamlog_out( MESSAGE5 ) << "EORE found: nothing else to do." << std::endl;
				return;
		} 
		else if( evt->getEventType() == kUNKNOWN )
		{
				streamlog_out( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
						<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}

		//Opens collection for input.
		LCCollection* inputTrackCollection = nullptr;
		LCCollection* inputHitCollection = nullptr;

		try
		{
				inputTrackCollection = evt->getCollection(_trackCollectionName);
				inputHitCollection = evt->getCollection(_hitCollectionName);
		}
		catch (DataNotAvailableException& e)
		{
				streamlog_out( WARNING2 ) << "Input collection not available" << std::endl;
				return;
		}

				//std::cout << "----------------- next event -------------------" << std::endl;

		auto hitEncoding = inputHitCollection->getParameters().getStringVal( LCIO::CellIDEncoding );
		if(hitEncoding.empty()) hitEncoding = EUTELESCOPE::HITENCODING;
		lcio::CellIDDecoder<TrackerHitImpl> hitDecoder ( hitEncoding );

		//Now get each individual hit LOOP OVER!
		for(int iHit = 0; iHit < inputHitCollection->getNumberOfElements(); ++iHit)
		{  
			auto inputHit = static_cast<TrackerHitImpl*>(inputHitCollection->getElementAt(iHit));

			//Call the local2masterHit/master2localHit function defined int EUTelGeometryTelescopeDescription
			//int properties = hitDecoder(inputHit)["properties"];
			int sensorID = hitDecoder(inputHit)["sensorID"];
			if(std::find(_DUTIDs.begin(), _DUTIDs.end(), sensorID) != _DUTIDs.end() ){
				auto inputPos = inputHit->getPosition();
				_DUTHits[sensorID].emplace_back(inputPos);
				//std::cout << "Pushed back " << inputPos[0] << ", " << inputPos[1] << ", " << inputPos[2] << " on plane " << sensorID << std::endl;
			}	
		}
	
		for(int iTrack = 0; iTrack < inputTrackCollection->getNumberOfElements(); ++iTrack)
		{  
			auto inputTrack = static_cast<TrackImpl*>(inputTrackCollection->getElementAt(iTrack));
			auto trackHits = inputTrack->getTrackerHits();
			//std::cout << "Track has " << trackHits.size() << " hits" << std::endl;
			bool matchedTrack = false;

			for(auto& trackHitRaw: trackHits) {
				auto trackHit = static_cast<TrackerHitImpl*>(trackHitRaw);
				int sensorID = hitDecoder(trackHit)["sensorID"];
				if(std::find(_DUTIDs.begin(), _DUTIDs.end(), sensorID) != _DUTIDs.end() ){
					auto predPos = trackHit->getPosition();
					//std::cout << "Track has " << trackHit->getPosition()[0] << ", " << trackHit->getPosition()[1] << ", " << trackHit->getPosition()[2] << " on plane " << sensorID << std::endl;
					for( auto& measPos: _DUTHits[sensorID]) {
						if( fabs( measPos[0] - predPos[0] ) < 0.125*1.5 && fabs( measPos[1] - predPos[1] ) < 0.05*1.5  ){
							//std::cout << "Track matched to DUT: " << sensorID << std::endl;
							matchedTrack = true;
							_DUTMatchedTrackes[sensorID]++;
							break;
						}
					}
				}
			}
			if(matchedTrack) _totalMatchedTracks++;	
		}		

	//Call the local2masterHit/master2localHit function defined int EUTelGeometryTelescopeDescription
			//int properties = hitDecoder(inputHit)["properties"];
		//	int sensorID = hitDecoder(inputHit)["sensorID"];
		//	if(std::find(_DUTIDs.begin(), _DUTIDs.end(), sensorID) != _DUTIDs.end() ){
		//		const double* inputPos = inputHit->getPosition();
		//		_DUTHits[sensorID].emplace_back(sensorID);		
			//}	

    for(auto& id: _DUTIDs) {
		 _DUTHits[id].clear();
	}

}

void EUTelProcessorDUTGlobalEfficiency::end()
{
	std::cout << "Total matched tracks: " << _totalMatchedTracks << std::endl;
	for(auto& pair: _DUTMatchedTrackes) {
		std::cout << "Sensor " << pair.first << " saw " << pair.second << " tracks!" << std::endl;
		
		double N = _totalMatchedTracks;
		double k = pair.second;
		double binomialUnc = (1/N)*sqrt(k*(1-k/N));

		std::cout << "Yielding an efficiency of: " << k/N << " +/- " << binomialUnc << std::endl;
	}
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}
