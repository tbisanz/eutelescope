/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorSelectHitRegion.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "CellIDReencoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;

EUTelProcessorSelectHitRegion::EUTelProcessorSelectHitRegion():
Processor("EUTelProcessorSelectHitRegion"),
_hitCollectionNameInput(), 
_hitCollectionNameOutput()
{
		_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the EUTelGeometryClass";

		registerInputCollection(LCIO::TRACKERHIT, "hitCollectionNameInput", "Input hit collection name", _hitCollectionNameInput, std::string("inputhit"));
		registerOutputCollection(LCIO::TRACKERHIT,"hitCollectionNameOutput", "Output hit collection name", _hitCollectionNameOutput, std::string ("outputhit"));

		registerOptionalParameter("xmax", "Set x-max value of selected area", _xMax,  static_cast<float>(0));
		registerOptionalParameter("xmin", "Set x-min value of selected area", _xMin,  static_cast<float>(0));
		registerOptionalParameter("ymax", "Set y-max value of selected area", _yMax,  static_cast<float>(0));
		registerOptionalParameter("ymin", "Set y-min value of selected area", _yMin,  static_cast<float>(0));
  		registerOptionalParameter("PlaneVec", "Planes on which the cuts shall be applied", _planes, std::vector<int>() );

}

void EUTelProcessorSelectHitRegion::processEvent(LCEvent* event) {

		//Check the event type and if it is the last event.
		EUTelEventImpl* evt	= static_cast<EUTelEventImpl*>(event);				
		if( evt->getEventType() == kEORE ) {
				streamlog_out( MESSAGE5 ) << "EORE found: nothing else to do." << std::endl;
				return;
		} else if( evt->getEventType() == kUNKNOWN ) {
				streamlog_out( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
						<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}

		//Opens collection for input.
		LCCollection* inputCollection = nullptr;
		try {
				inputCollection = evt->getCollection(_hitCollectionNameInput);
		} catch (DataNotAvailableException& e) {
				streamlog_out( WARNING2 ) << _hitCollectionNameInput << " collection not available" << std::endl;
				return;
		}

		LCCollectionVec* outputCollection = nullptr;
		try {
				outputCollection  = static_cast<LCCollectionVec*> (event->getCollection( _hitCollectionNameOutput ));
		} catch(...) {
				outputCollection = new LCCollectionVec(LCIO::TRACKERHIT);
		}

		std::string encoding = inputCollection->getParameters().getStringVal( LCIO::CellIDEncoding );

		if(encoding.empty()) {
			encoding = EUTELESCOPE::HITENCODING;
		}

		lcio::CellIDDecoder<TrackerHitImpl> hitDecoder ( encoding );

		//Now get each individual hit LOOP OVER!
		for(int iHit = 0; iHit < inputCollection->getNumberOfElements(); ++iHit) {  
			TrackerHitImpl*	inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));
			TrackerHitImpl* outputHit = new IMPL::TrackerHitImpl(); 

			int sensorID = hitDecoder(inputHit)["sensorID"];
			
			const double* inputPos = inputHit->getPosition();

			auto it = std::find(_planes.begin(), _planes.end(), sensorID);
			if(	 	it == _planes.end() || 
					(inputPos[0] >= _xMin &&
 					inputPos[0] <= _xMax &&
					inputPos[1] >= _yMin &&
					inputPos[1] <= _yMax)) {	
					
					outputHit->setPosition( inputPos );
					outputHit->setCovMatrix( inputHit->getCovMatrix());
					outputHit->setType( inputHit->getType() );
					outputHit->setTime( inputHit->getTime() );
					outputHit->setCellID0( inputHit->getCellID0() );
					outputHit->setCellID1( inputHit->getCellID1() );
					outputHit->setQuality( inputHit->getQuality() );
					outputHit->rawHits() = inputHit->getRawHits();

					outputCollection->push_back(outputHit);
				}
		}
	
		//Now push the hit for this event onto the collection
		try {	
				event->addCollection(outputCollection, _hitCollectionNameOutput );
		} catch(...) {
				streamlog_out ( WARNING5 )  << "Problem with pushing collection onto event" << std::endl;
		}
}

void EUTelProcessorSelectHitRegion::end() {
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}
