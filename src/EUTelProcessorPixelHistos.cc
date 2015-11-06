/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorPixelHistos.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// eutelescope geometry
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

#include "marlin/AIDAProcessor.h"
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogramFactory.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <EVENT/LCEvent.h>
#include <Exceptions.h>

// system includes <>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <iterator>
#include <algorithm>

using namespace marlin;
using namespace eutelescope;

EUTelProcessorPixelHistos::EUTelProcessorPixelHistos(): 
  Processor("EUTelProcessorPixelHistos"),
  _zsDataCollectionName(""),
  _iRun(0),
  _iEvt(0),
  _chosenSensorID(0)
{
  //processor description
  _description = "EUTelProcessorPixelHistos computes the firing frequency of pixels and applies a cut on this value to mask (NOT remove) hot pixels.";
  registerInputCollection(LCIO::TRACKERDATA, "ZSDataCollectionName", "Zero suppressed data input collection name", _zsDataCollectionName, std::string("zsdata") );
  registerOptionalParameter("MinX", "Min X pixel index to be included.", _minX, -1);
  registerOptionalParameter("MaxX", "Max X pixel index to be included.", _maxX, -1);
  registerOptionalParameter("MinY", "Min Y pixel index to be included.", _minY, -1);
  registerOptionalParameter("MaxY", "Max Y pixel index to be included.", _maxY, -1);
  registerProcessorParameter("SensorID", "Sensor ID of the investigated sensor", _chosenSensorID, 0);
}

std::pair<int, int> EUTelProcessorPixelHistos::mapPixels( int x, int y) {
	return std::make_pair(x-_minX, y-_minY);
}


void EUTelProcessorPixelHistos::init() {
	// this method is called only once even when the rewind is active usually a good idea to
	printParameters ();

	// set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;
	_noOfEvents = 0;

	//init new geometry
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,true);

	if( _minX == -1 ) _minX = 0;
	if( _maxX == -1 ) _maxX = geo::gGeometry().siPlaneXNpixels(_chosenSensorID)-1;
	if( _minY == -1 ) _minY = 0;
	if( _maxY == -1 ) _maxY = geo::gGeometry().siPlaneYNpixels(_chosenSensorID)-1;

	//_sensorIDVec = geo::gGeometry().sensorIDsVec();
	bookHistos();
}

void EUTelProcessorPixelHistos::processRunHeader(LCRunHeader* rdr) {
	std::unique_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl(rdr) );
	runHeader->addProcessor(type());
	// increment the run counter
	++_iRun;
	// reset the event counter
	_iEvt = 0;
}

void EUTelProcessorPixelHistos::processEvent(LCEvent* event) {
	if( event == nullptr ) {
		streamlog_out ( WARNING2 ) <<  "Event does not exist! Skipping!" <<  std::endl;       
		return;
	}

	EUTelEventImpl* evt = static_cast<EUTelEventImpl*>(event);
	if ( evt->getEventType() == kEORE ) {
		streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  std::endl;
		return;
	} else if ( evt->getEventType() == kUNKNOWN ) {
		streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	try {
		LCCollectionVec* zsInputCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );

		for( size_t iDetector = 0 ; iDetector < zsInputCollectionVec->size(); iDetector++ ) {

			// get the TrackerData and guess which kind of sparsified data it contains.
			TrackerDataImpl* zsData = dynamic_cast< TrackerDataImpl* > ( zsInputCollectionVec->getElementAt( iDetector ) );
			int sensorID            = static_cast<int > ( cellDecoder( zsData )["sensorID"] );

			if( sensorID != _chosenSensorID ) continue;

			// now prepare the EUTelescope interface to sparsified data.  
			std::unique_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData (new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ));
			EUTelGenericSparsePixel* genericPixel =  new EUTelGenericSparsePixel();

			// loop over all pixels in the sparseData object, these are the hit pixels!
			for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) {
				//get the pixel
				sparseData->getSparsePixelAt( iPixel, genericPixel );
				int xCo = genericPixel->getXCoord(); 
				int yCo = genericPixel->getYCoord();
				if( xCo <= _maxX && xCo >= _minX && yCo <= _maxY && yCo >= _minY) {
					auto indices = mapPixels(xCo, yCo);
					_arrSignal.at(indices.first).at(indices.second)->fill(genericPixel->getSignal());
					_arrTime.at(indices.first).at(indices.second)->fill(genericPixel->getTime());
				}
			}
			delete genericPixel;
		}

	} catch (lcio::DataNotAvailableException& e ) {
		streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << e.what() << std::endl;
		return;    
	}
	//don't forget to increment the event counter
	_iEvt++;
}

void EUTelProcessorPixelHistos::end() {
}

void EUTelProcessorPixelHistos::bookHistos() {
		auto basePath = "detector_" + to_string(_chosenSensorID);
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		auto const chargeHistoName = "raw_hits_charge_d"+to_string(_chosenSensorID)+"_px";
		auto const timeHistoName = "raw_hits_time_d"+to_string(_chosenSensorID)+"_px";
		
		int sizeX = _maxX-_minX+1;
		int sizeY = _maxY-_minY+1;

		//this is the 2-dimensional array used to store all the histograms
		_arrSignal = std::vector<std::vector<AIDA::IHistogram1D*>>(sizeX);
		_arrTime = std::vector<std::vector<AIDA::IHistogram1D*>>(sizeX);
		for(size_t i = 0; i < _arrSignal.size(); i++ ) {
			_arrSignal.at(i).resize(sizeY, nullptr);
			_arrTime.at(i).resize(sizeY, nullptr);
			for(int j = 0; j < sizeY; j++) {
				auto chargeName = basePath+chargeHistoName+to_string(i+_minX)+"_"+to_string(j+_minY);
				auto timeName = basePath+timeHistoName+to_string(i+_minX)+"_"+to_string(j+_minY);
				_arrSignal.at(i).at(j) = AIDAProcessor::histogramFactory(this)->createHistogram1D( chargeName.c_str(),  20, -0.5, 19.5);
				_arrTime.at(i).at(j) = AIDAProcessor::histogramFactory(this)->createHistogram1D( timeName.c_str(),  20, -0.5, 19.5);
			}	
		}
}
