// Version: $Id$
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelPreAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// gear includes <.h>
#include "marlin/Global.h"
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gearxml/GearXML.h>

// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <memory>
#include <cstdio>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelPreAlign::EUTelPreAlign () :Processor("EUTelPreAlign")
{
  _description = "Apply alignment constants to hit collection";

  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName","The name of the input hit collection",
                           _inputHitCollectionName, string ("hit"));

  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedID, 0);

  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored", 
			    _alignmentConstantLCIOFile, std::string( "alignment.slcio" ) );

  registerProcessorParameter ("Events", "How many events should be used for an approximation to the X,Y shifts (pre-alignment)? (default=50000)",
                              _events, static_cast <int> (50000) );

  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax, std::vector<float > (6,  10.) );

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax, std::vector<float > (6,  10.) );

  registerOptionalParameter ("MinNumberOfCorrelatedHits", "If there are more then this number of correlated hits (planes->track candidate) (default=5)",
			     _minNumberOfCorrelatedHits, static_cast <int> (5) );

  registerOptionalParameter("HistogramFilling","Switch on or off the histogram filling",_fillHistos, static_cast< bool > ( true ) );

}


void EUTelPreAlign::init () {
  //for info, good idea to:
  printParameters();

  // set to zero the run and event counters
  _iRun = 0;  
  _iEvt = 0;

  // clear the sensor ID vector
  _sensorIDVec.clear();

  // clear the sensor ID map
  _sensorIDVecMap.clear();
  _sensorIDtoZOrderMap.clear();

  // clear the sensor ID vector (z-axis order)
  _sensorIDVecZOrder.clear();

  //get parameters from GEAR
  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  for( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
  {

    if(_siPlanesLayerLayout->getID(iPlane) == _fixedID ) { 
      //Get Zpos of ref plane
      _fixedZ = _siPlanesLayerLayout->getSensitivePositionZ(iPlane); 
    } else {
      //Get 
      _preAligners.push_back( PreAligner( _siPlanesLayerLayout->getSensitivePitchX(iPlane) /10.,
					  _siPlanesLayerLayout->getSensitivePitchY(iPlane) /10.,
					  _siPlanesLayerLayout->getSensitivePositionZ(iPlane),
					  _siPlanesLayerLayout->getID(iPlane)) );
    }
  }


  _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
    {

      _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
      int sensorID = _siPlanesLayerLayout->getID( iPlane );

      _sensorIDVec.push_back( sensorID );
      _sensorIDVecMap.insert( make_pair( sensorID, iPlane ) );

      // count number of the sensors to the left of the current one:
      int _sensors_to_the_left = 0;
      for ( int jPlane = 0 ; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++ ) 
	{
	  if( _siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 < _siPlaneZPosition[ iPlane ] )
	    {
	      _sensors_to_the_left++;
	    }
	}
      streamlog_out ( DEBUG5 ) << " iPlane " << iPlane << " sensor_#_along_Z_axis " << _sensors_to_the_left << "[z= " << setprecision (3) << _siPlaneZPosition[iPlane] << " ] [sensorID " << sensorID << " ]  " << endl;

      _sensorIDVecZOrder.push_back( _sensors_to_the_left );
      _sensorIDtoZOrderMap.insert( make_pair( sensorID, _sensors_to_the_left ) );
    }

  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
    {
      _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
      int sensorID = _siPlanesLayerLayout->getID( iPlane );
      _sensorIDinZordered.insert( make_pair( _sensorIDtoZOrderMap[ sensorID ], sensorID ) );
    }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  string tempHistoName = "";
  string basePath; 

  if( _fillHistos ) {

    // Allow any plane to be the fixed reference:
    for(unsigned int i = 0; i < _sensorIDVecZOrder.size(); i++)
      {
	int sensorID = _sensorIDinZordered[i];

	basePath = "plane_" + to_string( sensorID );
	AIDAProcessor::tree(this)->mkdir(basePath.c_str());
	basePath.append("/");
 
	tempHistoName = "hitXCorr_fixed_to_" + to_string( sensorID ) ;
	AIDA::IHistogram1D * histo1Da = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.);
	_hitXCorr.insert( make_pair( sensorID, histo1Da) );
 
	tempHistoName = "hitYCorr_fixed_to_" + to_string( sensorID) ;
	AIDA::IHistogram1D * histo1Db = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.) ;
	_hitYCorr.insert( make_pair( sensorID, histo1Db) );
      }
  }
#endif

}

void EUTelPreAlign::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
  runHeader->addProcessor( type() );
  ++_iRun;
}

void EUTelPreAlign::processEvent (LCEvent * event) {

  ++_iEvt;

  if(_iEvt > _events) return;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  
  if(  evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if(  evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  try {
    LCCollectionVec * inputCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
    
    std::vector<float> residX;
    std::vector<float> residY;
    std::vector<PreAligner*> prealign;

    //Loop over all hits and determine the hit on the fixed plane:
    for( size_t ref = 0; ref < inputCollectionVec->size(); ref++ )
    {
      TrackerHitImpl* refHit = dynamic_cast<TrackerHitImpl*>( inputCollectionVec->getElementAt(ref) );
      const double* refPos = refHit->getPosition();

      //if the hit is the one on the fixed plane, we go further, otherwise we try the next hit
      int sensorID = hitDecoder(refHit)["sensorID"];
      if( sensorID != _fixedID ) continue;

      residX.clear();
      residY.clear();
      prealign.clear();

      //Now we again loop over all hits and compute the residuals in comparison
      for( size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++)
      {
	TrackerHitImpl* hit = dynamic_cast<TrackerHitImpl*>( inputCollectionVec->getElementAt(iHit) );
        
        const double * pos = hit->getPosition();

	//we dont consider the case that they equal
        int iHitID = hitDecoder(hit)["sensorID"]; 
        if( iHitID == _fixedID ) continue;
        
	bool gotIt = false;

	for(size_t ii = 0; ii < _preAligners.size(); ii++)
	{
	    PreAligner& pa = _preAligners.at(ii);

	    if( pa.getIden() != iHitID  ) continue;

	    gotIt = true;

	    double correlationX =  refPos[0] - pos[0] ;
	    double correlationY =  refPos[1] - pos[1] ;

	    int idZ = _sensorIDtoZOrderMap[ iHitID ];

	    //check if the residuals are inside the defined window
	    if( 
	       (_residualsXMin[idZ] < correlationX ) && ( correlationX < _residualsXMax[idZ]) &&
	       (_residualsYMin[idZ] < correlationY ) && ( correlationY < _residualsYMax[idZ]) 
		) 
	    {
	      residX.push_back( correlationX );
	      residY.push_back( correlationY );
	      prealign.push_back(&pa);
	    }
	    break;
        }
	if( !gotIt ) 
	{
	    streamlog_out ( ERROR5 ) << "Mismatched hit at " << pos[2] << endl;
	}
      }

      if( prealign.size() > static_cast< unsigned int >(_minNumberOfCorrelatedHits) && residX.size() == residY.size() ) {
	for( unsigned int ii = 0 ;ii < prealign.size(); ii++ ) {

	  prealign[ii]->addPoint( residX[ii], residY[ii] );

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	  if( _fillHistos ) {
	    ( dynamic_cast<AIDA::IHistogram1D*> (_hitXCorr[ prealign[ii]->getIden() ] ) )->fill( residX[ii] );
	    ( dynamic_cast<AIDA::IHistogram1D*> (_hitYCorr[ prealign[ii]->getIden() ] ) )->fill( residY[ii] );
	  }
#endif
	}
      }
    }
  }
  catch( DataNotAvailableException& e) { 
    streamlog_out  ( WARNING2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

  if( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelPreAlign::end() {
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW    );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  streamlog_out ( MESSAGE5 ) << "Writing to " << _alignmentConstantLCIOFile << endl;

  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;
  LCEventImpl * event = new LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );
  LCTime * now = new LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );


  for(size_t ii=0; ii<_sensorIDVec.size(); ii++)
    {
      bool ifound = false;
      for(size_t jj=0; jj< _preAligners.size(); jj++)
	{
	  int sensorID = _preAligners.at(jj).getIden();
	  if( _sensorIDVec[ii] == sensorID ) { ifound = true; break; }
	}
      if( ifound == false)
	{
	  EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
	  constant->setXOffset( 0.0 );
	  constant->setYOffset( 0.0 );
	  constant->setSensorID( _sensorIDVec[ii] );
	  constantsCollection->push_back( constant );
	  streamlog_out ( MESSAGE5 ) << (*constant) << endl;
	  continue; 
	}
    }


  for(size_t ii = 0 ; ii < _preAligners.size(); ii++){
    EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();

    if( abs( _preAligners.at(ii).getPeakX() ) < 1000. )
      constant->setXOffset( -1.0 * _preAligners.at(ii).getPeakX() );
    else
      constant->setXOffset( -1.0 * 0.0                           );
 
    if( abs( _preAligners.at(ii).getPeakY() ) < 1000. )
      constant->setYOffset( -1.0 * _preAligners.at(ii).getPeakY() );
    else
      constant->setYOffset( -1.0 * 0.0                           );
 
    int sensorID = _preAligners.at(ii).getIden();
    constant->setSensorID( sensorID );
    constantsCollection->push_back( constant );
    streamlog_out ( MESSAGE5 ) << (*constant) << endl;
  
    int GEARID = -1;
    for(int GEARLayers = 0; GEARLayers < _siPlanesLayerLayout->getNLayers(); GEARLayers++) 
    {
	if(_siPlanesLayerLayout->getID(GEARLayers) == sensorID )
	{
		GEARID = GEARLayers;
		break;
	}
    }
    
    double sensX = _siPlanesLayerLayout->getLayerPositionX(GEARID);
    double sensY = _siPlanesLayerLayout->getLayerPositionY(GEARID);
  
    sensX += _preAligners.at(ii).getPeakX();
    sensY += _preAligners.at(ii).getPeakY();

    _siPlanesLayerLayout->setLayerPositionX(GEARID, sensX);
    _siPlanesLayerLayout->setLayerPositionY(GEARID, sensY);

  }

  streamlog_out( DEBUG5 ) << " adding Collection " << "alignment " << endl;
 
  event->addCollection( constantsCollection, "alignment" );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();

  GearXML::createXMLFile(Global::GEAR,"updatedGEAR.xml");
}
#endif
