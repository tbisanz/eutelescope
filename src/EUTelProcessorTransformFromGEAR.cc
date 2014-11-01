// eutelescope includes ".h"
#include "EUTelProcessorTransformFromGEAR.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// gear includes <.h>
#include "marlin/Global.h"
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// system includes <>
#include <memory>
using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelProcessorTransformFromGEAR::EUTelProcessorTransformFromGEAR () :Processor("EUTelProcessorTransformFromGEAR")
{
  _description = "Apply alignment constants to hit collection";

  registerInputCollection(LCIO::TRACKERHIT, "InputHitCollectionName","The name of the input hit collection in the local frame", _inputHitCollectionName, std::string("inhit"));
  
  registerOutputCollection(LCIO::TRACKERHIT, "OutputHitCollectionName","The name of the output hit collection in the global frame", _outputHitCollectionName, std::string("outhit"));
}


void EUTelProcessorTransformFromGEAR::init()
{
	//for info, good idea to:
	printParameters();

	// set to zero the run and event counters
	_iRun = 0;  
	_iEvt = 0;

	//get parameters from GEAR
	_siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

	for( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
	{
		GEAREntries entries;

		int sensorID = _siPlanesLayerLayout->getID(iPlane);
		
		entries.offX = _siPlanesLayerLayout->getLayerPositionX(iPlane);
		entries.offY = _siPlanesLayerLayout->getLayerPositionY(iPlane);
		entries.offZ = _siPlanesLayerLayout->getLayerPositionZ(iPlane);

		entries.r1 = _siPlanesLayerLayout->getSensitiveRotation1(iPlane);
		entries.r2 = _siPlanesLayerLayout->getSensitiveRotation2(iPlane);
		entries.r3 = _siPlanesLayerLayout->getSensitiveRotation3(iPlane);
		entries.r4 = _siPlanesLayerLayout->getSensitiveRotation4(iPlane);

		entries.alpha = _siPlanesLayerLayout->getLayerRotationZY(iPlane);
		entries.beta = _siPlanesLayerLayout->getLayerRotationZX(iPlane);
		entries.gamma = _siPlanesLayerLayout->getLayerRotationXY(iPlane);
	
		_GEAREntriesMap[sensorID] = entries;

		std::cout << "For sensor: " << sensorID << "Offset vaules: " << entries.offX << "," << entries.offY << "," << entries.offZ << std::endl;
	}

	//Prepare Rotation Matrix!
}

void EUTelProcessorTransformFromGEAR::processRunHeader(LCRunHeader* rdr)
{
	auto_ptr<EUTelRunHeaderImpl>runHeader(new EUTelRunHeaderImpl( rdr ));
	runHeader->addProcessor(type());
	++_iRun;
}

void EUTelProcessorTransformFromGEAR::processEvent(LCEvent* event)
{
	++_iEvt;
 	LCCollectionVec* inputCollectionVec = NULL;
	
	try
	{
		inputCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_inputHitCollectionName));

  	}

	catch(DataNotAvailableException& e)
	{
	       return;	
	}
	
	UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder( inputCollectionVec );
	
	//now prepare output collection
	LCCollectionVec* outputCollectionVec;
	bool outputCollectionVecExists = false;
  	_initialOutputCollectionSize = 0;

  	try 
  	{
   		outputCollectionVec = dynamic_cast< LCCollectionVec* > ( event->getCollection( _outputHitCollectionName ) );
    		outputCollectionVecExists = true;
    		_initialOutputCollectionSize = outputCollectionVec->size();
  	} 
  	catch ( lcio::DataNotAvailableException& e ) 
  	{
    		outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
  	}

        //read the encoding string from the input collection
	std::string encodingString = inputCollectionVec->getParameters().getStringVal( LCIO::CellIDEncoding );	
	//and the encoder for the output data
	CellIDEncoder<TrackerHitImpl> idZSGenDataEncoder( encodingString , outputCollectionVec);


	//Loop over all hits and determine the hit on the fixed plane:
	for(size_t hitNo = 0; hitNo < inputCollectionVec->size(); hitNo++)
	{
		TrackerHitImpl* inputHit = dynamic_cast<TrackerHitImpl*>(inputCollectionVec->getElementAt(hitNo));
		
		int sensorID = hitDecoder(inputHit)["sensorID"];
		const double* hitPos = inputHit->getPosition();
		
		double newPos[3];

		newPos[0] = hitPos[0] -  _GEAREntriesMap[sensorID].offX; 
		newPos[1] = hitPos[1] -  _GEAREntriesMap[sensorID].offY; 
		newPos[2] = hitPos[2]; 

		//TrackerHitImpl for the output collection
		auto_ptr<TrackerHitImpl> outputHit ( new TrackerHitImpl );

		//copy the information which is the same
		outputHit->setCellID0( inputHit->getCellID0() );
		outputHit->setCellID1( inputHit->getCellID1() );
		outputHit->setTime( inputHit->getTime() );
		outputHit->setCovMatrix( inputHit->getCovMatrix() );
		outputHit->setQuality( inputHit->getQuality() );
		outputHit->rawHits() =  inputHit->getRawHits();
		outputHit->setPosition( &newPos[0] );

		outputCollectionVec->push_back( outputHit.release() );
	}


	//add the collection if necessary
	if ( !outputCollectionVecExists && ( outputCollectionVec->size() != _initialOutputCollectionSize )) 
	{
		event->addCollection( outputCollectionVec, _outputHitCollectionName );
	}

	if ( !outputCollectionVecExists && ( outputCollectionVec->size() == _initialOutputCollectionSize ) ) 
	{
		delete outputCollectionVec;
	}	
	//rest of memory cleaned up by auto_ptrs
}

void EUTelProcessorTransformFromGEAR::end()
{
	//NOP NOP NOP
}
