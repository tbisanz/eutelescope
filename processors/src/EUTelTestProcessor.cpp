/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelTestProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "CellIDReencoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

#include "marlin/AIDAProcessor.h"
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;

EUTelTestProcessor::EUTelTestProcessor():
Processor("EUTelTestProcessor"),
_hitCollection1NameInput(),
_hitCollection2NameInput()	
{
		_description ="This is a test to learn what we're doing";
		registerInputCollection(LCIO::TRACKERHIT, "hitCollection1NameInput", "Test Input Collection", _hitCollection1NameInput, std::string("local_hit"));

                _description ="This is a second Test to learn what we're doing";
                registerInputCollection(LCIO::TRACKERHIT, "hitCollection2NameInput", "Test Input Collection2", _hitCollection2NameInput, std::string("hit"));
}


void EUTelTestProcessor::init() {
		geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

		_sensorIDVec = geo::gGeometry().sensorIDsVec();
		bookHistos();

		std::cout << "EUTelTestProcessor::init()" << std::endl;
}

void EUTelTestProcessor::processRunHeader(LCRunHeader* rdr)
{
	//maybe put some stuff here?
}

void EUTelTestProcessor::processEvent(LCEvent* event)
{
		//Opens collection for input.
		LCCollection* inputCollection = nullptr;
		EUTelEventImpl *evt = static_cast<EUTelEventImpl *>(event);
	//	try {
				inputCollection = evt->getCollection(_hitCollection1NameInput);
	//	} catch (DataNotAvailableException& e) {
	//			streamlog_out( WARNING2 ) << _hitCollectionNameInput << " collection not available" << std::endl;
	//			return;
	//	}
		LCCollection* inputCollection2 = nullptr;
		
				inputCollection2 = evt->getCollection(_hitCollection2NameInput);
		/*std::string encoding = inputCollection->getParameters().getStringVal( LCIO::CellIDEncoding );
		if(encoding.empty()) {
			encoding = EUTELESCOPE::HITENCODING;
		}*/
		lcio::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );

		//Now get each individual hit LOOP OVER!
		for(int iHit = 0; iHit < inputCollection->getNumberOfElements(); ++iHit) {  
			TrackerHitImpl*	inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));
			TrackerHitImpl* inputHit2 = static_cast<TrackerHitImpl*>(inputCollection2->getElementAt(iHit));
			//Call the local2masterHit/master2localHit function defined int EUTelGeometryTelescopeDescription
			//int properties = hitDecoder(inputHit)["properties"];
			int sensorID = hitDecoder(inputHit)["sensorID"];
			const double* inputPos = inputHit->getPosition();
			const double* inputPos2 = inputHit2->getPosition();
			std::array<double, 3> inputArray1 = { inputPos[0], inputPos[1], inputPos[2] };
			std::array<double, 3> inputArray2 = { inputPos2[0], inputPos2[1], inputPos2[2] };
			double X=0;
			double X2=0;
			double Y=0;
			double Y2=0;
			//std::cout << "local x: " << *inputPos << " local Y: " << *(inputPos+1) << " local Z: " << *(inputPos+2) << std::endl;
                        //std::cout << "global X: " << *inputPos2 << " global Y: " << *(inputPos2+1) << " global Z: " << *(inputPos2+2) << std::endl;
			if(*(inputPos)<0){
				X=*(inputPos)*-1;
			}
			else{
				X=*(inputPos);
			}
			if(*(inputPos2)<0){
                                X2=*(inputPos2)*-1;
                        }
			else{
				X2=*(inputPos2);
			}
			if(*(inputPos+1)<0){
                                Y=*(inputPos+1)*-1;
                        }
                        else{
                                Y=*(inputPos+1);
                        }
                        if(*(inputPos2+1)<0){
                                Y2=*(inputPos2+1)*-1;
                        }
                        else{
                                Y2=*(inputPos2+1);
                        }

			double difX = X - X2;
			double difY = Y - Y2;
			//double difZ = *(Pos+2) - *(Pos2+2);
			double dif = sqrt(difX * difX+difY * difY);
			//std::cout << "Difference between local and global coordinates without Z: " << dif << std::endl;
			_HitDifHisto.at(sensorID)->fill(dif);

		}
}

void EUTelTestProcessor::end()
{
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}

void EUTelTestProcessor::bookHistos() {
       	streamlog_out(MESSAGE1) << "Booking histograms " << std::endl;

 	 for (auto sensorID : _sensorIDVec) {
   		 auto basePath = "detector_" + to_string(sensorID);
		 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
  		 basePath.append("/");

		 auto const countHistoName = "difference_in_local_and_global_coordinates" + to_string(sensorID);

		 auto HitDifHisto =
       			 marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
           		 (basePath + countHistoName).c_str(), 100, -0.02, 0.02);

		 _HitDifHisto[sensorID] = HitDifHisto;
	}
}

	
