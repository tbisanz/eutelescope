// Contact: Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
//
// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTelMille2.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "mille/Mille.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <Eigen/Core>
#include <Eigen/QR>

// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TMinuit.h>
#include <TSystem.h>
#include <TMath.h>
#include <TVector3.h>
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

//evil global variables ...
EUTelMille2::hit2* hitsarray2;
unsigned int number_of_datapoints2;
void fcn_wrapper2(int& /*npar*/, double* /*gin*/, double &f, double *par, int /*iflag*/)
{
	EUTelMille2::trackfitter2 fobj(hitsarray2,number_of_datapoints2);
	f = fobj.fit(par);
}


EUTelMille2::EUTelMille2 () : Processor("EUTelMille2") {

  //some default values
  FloatVec MinimalResidualsX;
  FloatVec MinimalResidualsY;
  FloatVec MaximalResidualsX;
  FloatVec MaximalResidualsY;

  FloatVec PedeUserStartValuesX;
  FloatVec PedeUserStartValuesY;

  FloatVec PedeUserStartValuesGamma;

  FloatVec SensorZPositions;

  FloatVec SensorXShifts;
  FloatVec SensorYShifts;

  FloatVec SensorGamma;
  FloatVec SensorAlpha;
  FloatVec SensorBeta;

  //maybe one has to chose a larger value than 6?
  for(int i =0; i<6;i++)
    {
      MinimalResidualsX.push_back(0.0);
      MinimalResidualsY.push_back(0.0);
      MaximalResidualsX.push_back(0.0);
      MaximalResidualsY.push_back(0.0);

      PedeUserStartValuesX.push_back(0.0);
      PedeUserStartValuesY.push_back(0.0);

      PedeUserStartValuesGamma.push_back(0.0);

      float zpos = 20000.0 +  20000.0 * static_cast< float >(i);
      SensorZPositions.push_back(zpos);

      SensorXShifts.push_back(0.0);
      SensorYShifts.push_back(0.0);

      SensorGamma.push_back(0.0);
      SensorAlpha.push_back(0.0);
      SensorBeta.push_back(0.0);
    }

  // modify processor description
  _description = "EUTelMille2 uses the MILLE program to write data files for MILLEPEDE II.";

  // choose input mode
  registerOptionalParameter("InputMode","Selects the source of input hits."
                            "\n0 - hits read from hitfile and simple straight line trackfinding will be  performed internally."
                            "\n1 - hits read from output of tracking processor."
                            "\n(DEFUNCT!!) 2 - Test mode. Simple internal simulation and simple trackfinding."
                            "\n(DEFUNCT) 3 - Mixture of a track collection from the telescope and hit collections for the DUT (only one DUT layer can be used unfortunately)",
                            _inputMode, static_cast <int> (0));

  registerOptionalParameter("AllowedMissingHits","Set how many hits (=planes) can be missing on a track candidate.",  _allowedMissingHits, static_cast <int> (0));

  // input collections
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("corrhits");

  registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName", "Hit collections name", _hitCollectionName,HitCollectionNameVecExample);

  registerInputCollection(LCIO::TRACK,"TrackCollectionName", "Track collection name. This is only relevant if InputMode is set to larger to 1", _trackCollectionName,std::string("fittracks"));

  // parameters
  registerOptionalParameter("DistanceMax","Maximal allowed distance between hits entering the fit per 10 cm space between the planes.", _distanceMax, static_cast <float> (2000.0));

  registerOptionalParameter("DistanceMaxVec","Maximal allowed distance between hits entering the fit per 10 cm space between the planes. One value for each neighbor planes. DistanceMax will be used for each pair if this vector is empty.", _distanceMaxVec, FloatVec ());

  registerOptionalParameter("ExcludePlanes","Exclude planes from fit according to their sensor ids.",_excludePlanes_sensorIDs ,IntVec());

  registerOptionalParameter("FixedPlanes","Fix sensor planes in the fit according to their sensor ids.",_FixedPlanes_sensorIDs ,IntVec());

  registerOptionalParameter("MaxTrackCandidatesTotal","Stop processor after this maximum number of track candidates (Total) is reached.",_maxTrackCandidatesTotal, static_cast <int> (10000000));

  registerOptionalParameter("MaxTrackCandidates","Maximal number of track candidates in a event.",_maxTrackCandidates, static_cast <int> (2000));

  registerOptionalParameter("BinaryFilename","Name of the Millepede binary file.",_binaryFilename, string ("mille.bin"));

  registerOptionalParameter("TelescopeResolution","(default) Resolution of the telescope for Millepede (sigma_x=sigma_y) used only if plane dependent resolution is set inconsistently.",_telescopeResolution, static_cast <float> (3.0));

  registerOptionalParameter("OnlySingleHitEvents","Use only events with one hit in every plane.",_onlySingleHitEvents, static_cast <bool> (false));

  registerOptionalParameter("OnlySingleTrackEvents","Use only events with one track candidate.",_onlySingleTrackEvents, static_cast <bool> (false));

  registerOptionalParameter("AlignMode","Number of alignment constants used. Available mode are: "
                            "\n1 - shifts in the X and Y directions and a rotation around the Z axis,"
                            "\n2 - only shifts in the X and Y directions"
                            "\n3 - (EXPERIMENTAL) shifts in the X,Y and Z directions and rotations around all three axis",
                            _alignMode, static_cast <int> (1));

  registerOptionalParameter("UseResidualCuts","Use cuts on the residuals to reduce the combinatorial background.",_useResidualCuts, static_cast <bool> (false));

  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin,MinimalResidualsX);

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin,MinimalResidualsY);

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax,MaximalResidualsX);

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax,MaximalResidualsY);

  registerOptionalParameter("ResolutionX","X resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionX,  FloatVec (static_cast <int> (6), 10.));

  registerOptionalParameter("ResolutionY","Y resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionY,FloatVec (static_cast <int> (6), 10.));

  registerOptionalParameter("ResolutionZ","Z resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionZ,FloatVec (static_cast <int> (6), 10.));

  registerOptionalParameter("FixParameter","Fixes the given alignment parameters in the fit if alignMode==3 is used. For each sensor an integer must be specified (If no value is given, then all parameters will be free). bit 0 = x shift, bit 1 = y shift, bit 2 = z shift, bit 3 = alpha, bit 4 = beta, bit 5 = gamma. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_FixParameter, IntVec (static_cast <int> (6), 24));

  registerOptionalParameter("GeneratePedeSteerfile","Generate a steering file for the pede program.",_generatePedeSteerfile, static_cast <int> (0));

  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));

  registerOptionalParameter("PedeSteeringAdditionalCmds","FOR EXPERTS: List of commands that should be included in the pede steering file. Use '\\' to seperate options and introduce a line break.",_pedeSteerAddCmds, StringVec());

  registerOptionalParameter("UsePedeUserStartValues","Give start values for pede by hand (0 - automatic calculation of start values, 1 - start values defined by user).", _usePedeUserStartValues, static_cast <int> (0));

  registerOptionalParameter("PedeUserStartValuesX","Start values for the alignment for shifts in the X direction.",_pedeUserStartValuesX,PedeUserStartValuesX);

  registerOptionalParameter("PedeUserStartValuesY","Start values for the alignment for shifts in the Y direction.",_pedeUserStartValuesY,PedeUserStartValuesY);

  registerOptionalParameter("PedeUserStartValuesZ","Start values for the alignment for shifts in the Z direction.",_pedeUserStartValuesZ,FloatVec (static_cast <int> (6), 0.0));

  registerOptionalParameter("PedeUserStartValuesAlpha","Start values for the alignment for the angle alpha.",_pedeUserStartValuesAlpha,FloatVec (static_cast <int> (6), 0.0));
  
  registerOptionalParameter("PedeUserStartValuesBeta","Start values for the alignment for the angle beta.",_pedeUserStartValuesBeta,FloatVec (static_cast <int> (6), 0.0));
  
  registerOptionalParameter("PedeUserStartValuesGamma","Start values for the alignment for the angle gamma.",_pedeUserStartValuesGamma,PedeUserStartValuesGamma);
}

void EUTelMille2::init()
{
	// Getting access to geometry description
	std::string name("test.root");
	geo::gGeometry().initializeTGeoDescription(name,false);


	// check if the GEAR manager pointer is not null!
	if( Global::GEAR == 0x0 )
	{
		streamlog_out( ERROR2 ) << "The GearMgr is not available, for an unknown reason." << endl;
		throw InvalidGeometryException("GEAR manager is not initialised");
	}

	//  sensor-planes in geometry navigation:
	_siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

	// clear the sensor ID vector
	_sensorIDVec.clear();

	// clear the sensor ID map
	_sensorIDVecMap.clear();
	_sensorIDtoZOrderMap.clear();

	// clear the sensor ID vector (z-axis order)
	_sensorIDVecZOrder.clear();

	// copy-paste from another class (should be ideally part of GEAR!)
	double* keepZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
	for( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
	{
		int sensorID = _siPlanesLayerLayout->getID( iPlane );
		keepZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);

		_sensorIDVec.push_back( sensorID );
		_sensorIDVecMap.insert( make_pair( sensorID, iPlane ) );

		// count number of the sensors to the left of the current one:
		int _sensors_to_the_left = 0;
		for( int jPlane = 0 ; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++ ) 
		{
			if( _siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 <     keepZPosition[ iPlane ] )
			{
				_sensors_to_the_left++;
			}
		}

		_sensorIDVecZOrder.push_back( _sensors_to_the_left );
		_sensorIDtoZOrderMap.insert(make_pair( sensorID, _sensors_to_the_left));
	}

	delete[] keepZPosition;

	//lets guess the number of planes
	_nPlanes = geo::gGeometry().nPlanes();

	// an associative map for getting also the sensorID ordered
	map< double, int > sensorIDMap;
	//lets create an array with the z positions of each layer
	for( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
		_siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
		sensorIDMap.insert( make_pair( _siPlanesLayerLayout->getLayerPositionZ(iPlane), _siPlanesLayerLayout->getID(iPlane) ) );
	}

	//lets sort the array with increasing z
	sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

	//the user is giving sensor ids for the planes to be excluded. this
	//sensor ids have to be converted to a local index according to the
	//planes positions along the z axis.
	for(size_t i = 0; i < _FixedPlanes_sensorIDs.size(); i++)
	{
		map< double, int >::iterator iter = sensorIDMap.begin();
		int counter = 0;
		while( iter != sensorIDMap.end() ) {
			if( iter->second == _FixedPlanes_sensorIDs[i])
			{
				_FixedPlanes.push_back(counter);
				break;
			}
			++iter;
			++counter;
		}
	}
	for(size_t i = 0; i < _excludePlanes_sensorIDs.size(); i++)
	{
		map< double, int >::iterator iter = sensorIDMap.begin();
		int counter = 0;
		while( iter != sensorIDMap.end() ) {
			if( iter->second == _excludePlanes_sensorIDs[i])
			{
				_excludePlanes.push_back(counter);
				break;
			}
			++iter;
			++counter;
		}
	}

	// strip from the map the sensor id already sorted.
	map< double, int >::iterator iter = sensorIDMap.begin();
	unsigned int counter = 0;
	while( iter != sensorIDMap.end() ) {
		bool excluded = false;
		for(size_t i = 0; i < _excludePlanes.size(); i++)
		{
			if(_excludePlanes[i] == counter)
			{
				excluded = true;
				break;
			}
		}
		if(!excluded)
			_orderedSensorID_wo_excluded.push_back( iter->second );
		_orderedSensorID.push_back( iter->second );

		++iter;
		++counter;
	}

	//consistency
	if(_siPlaneZPosition.size() != _nPlanes)
	{
		streamlog_out( ERROR2 ) << "the number of detected planes is " << _nPlanes << " but only " << _siPlaneZPosition.size() << " layer z positions were found!"  << endl;
		throw InvalidParameterException("number of layers and layer z positions mismatch");
	}

	// this method is called only once even when the rewind is active
	// usually a good idea to
	printParameters ();

	// set to zero the run and event counters
	_iRun = 0;
	_iEvt = 0;

	// Initialize number of excluded planes
	_nExcludePlanes = _excludePlanes.size();

	streamlog_out( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

	// Initialise Mille statistics
	_nMilleDataPoints = 0;
	_nMilleTracks = 0;

	_waferResidX = new double[_nPlanes];
	_waferResidY = new double[_nPlanes];
	_waferResidZ = new double[_nPlanes];

	_xFitPos = new double[_nPlanes];
	_yFitPos = new double[_nPlanes];

	//check the consistency of the resolution parameters
	if(_alignMode == 3)
	{
		if( _resolutionX.size() != _resolutionY.size() )
		{
			throw InvalidParameterException("WARNING, length of resolution X and Y is not the same \n");
		}
		if( _resolutionY.size() != _resolutionZ.size() )
		{
			throw InvalidParameterException("WARNING, length of resolution Y and Z is not the same \n");
		}

		if(
				_resolutionX.size() != static_cast<unsigned int>(_nPlanes ) ||
				_resolutionY.size() != static_cast<unsigned int>(_nPlanes ) ||
				_resolutionZ.size() != static_cast<unsigned int>(_nPlanes )
		  )
		{
			streamlog_out( WARNING2 ) << "Consistency check of the resolution parameters failed. The array size is different than the number of found planes! "
				"The resolution parameters are set to default values now (see variable TelescopeResolution). "
				"This introduces a bias if the real values for X,Y and Z are rather different." << endl;
			_resolutionX.clear();
			_resolutionY.clear();
			_resolutionZ.clear();
			for(size_t i = 0; i < _nPlanes; i++)
			{
				_resolutionX.push_back(_telescopeResolution);
				_resolutionY.push_back(_telescopeResolution);
				_resolutionZ.push_back(_telescopeResolution);
			}
		}

		if(_FixParameter.size() != static_cast<unsigned int>(_nPlanes ) && !_FixParameter.empty())
		{       
			streamlog_out( WARNING2 ) << "Consistency check of the fixed parameters array failed. The array size is different than the number of found planes! The array is now set to default values, which means that all parameters are free in the fit." << endl;
			_FixParameter.clear();
			for(size_t i = 0; i < _nPlanes; i++)
			{
				_FixParameter.push_back(0);
			}
		}
		if(_FixParameter.empty())
		{
			streamlog_out( WARNING2 ) << "The fixed parameters array was found to be empty. It will be filled with default values. All parameters are free in the fit now." << endl;
			_FixParameter.clear();
			for(size_t i = 0; i < _nPlanes; i++)
			{
				_FixParameter.push_back(0);
			}
		}
		
		number_of_datapoints2 = _nPlanes;
		hitsarray2 = new hit2[number_of_datapoints2];
	}

	streamlog_out( MESSAGE5 ) << "Initialising Mille..." << endl;

	_mille = new Mille(_binaryFilename.c_str());

	_xPos.clear();
	_yPos.clear();
	_zPos.clear();

	_trackResidX.clear();
	_trackResidY.clear();
	_trackResidZ.clear();

	for(int i = 0; i < _maxTrackCandidates; i++)
	{
		_xPos.push_back(DoubleVec(_nPlanes,0.0));
		_yPos.push_back(DoubleVec(_nPlanes,0.0));
		_zPos.push_back(DoubleVec(_nPlanes,0.0));

		_trackResidX.push_back(DoubleVec(_nPlanes,0.0));
		_trackResidY.push_back(DoubleVec(_nPlanes,0.0));
		_trackResidZ.push_back(DoubleVec(_nPlanes,0.0));
	}

	if(!_distanceMaxVec.empty())
	{
		if(_distanceMaxVec.size() !=  static_cast<unsigned int>(_nPlanes ) )
		{
			streamlog_out( WARNING2 ) << "Consistency check of the DistanceMaxVec array failed. Its size is different compared to the number of planes! Will now use _distanceMax for each pair of planes." << endl;
			_distanceMaxVec.clear();
			for(size_t i = 0; i < _nPlanes-1; i++)
			{
				_distanceMaxVec.push_back(_distanceMax);
			}
		}
	}
	else
	{
		_distanceMaxVec.clear();
		for(size_t i = 0; i < _nPlanes-1; i++)
		{
			_distanceMaxVec.push_back(_distanceMax);
		}
	}

	streamlog_out( MESSAGE4 ) << "end of initialisation" << endl;
}

void EUTelMille2::processRunHeader(LCRunHeader* rdr) {

  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl(rdr) );
  header->addProcessor( type() );

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if( header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID() ) {
    streamlog_out( ERROR2 ) << "Error during the geometry consistency check: " << endl;
    streamlog_out( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << endl;
    streamlog_out( ERROR2 ) << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << endl;
  }

  // increment the run counter
  ++_iRun;
}

void EUTelMille2::findtracks2(
                            int missinghits,
                            std::vector<IntVec > &indexarray,
                            IntVec vec,
                            std::vector<std::vector<EUTelMille2::HitsInPlane2> > &_allHitsArray,
                            unsigned int i,
                            int y
                            )
{
	if(y==-1) missinghits++;

	if( missinghits > getAllowedMissingHits() ) 
	{
		// recursive chain is dropped here;
		streamlog_out(DEBUG9) << "indexarray size:" << indexarray.size() << std::endl;
		return;
	}

	if(i>0)
	{ 
		vec.push_back(y); // recall hit id from the plane (i-1)
	}


	if( (_allHitsArray[i].size() == 0) && (i <_allHitsArray.size()-1) )
	{
		findtracks2(missinghits,indexarray,vec, _allHitsArray, i+1, -1 ); 
	} 

	for(size_t j =0; j < _allHitsArray[i].size(); j++)
	{
		int ihit = static_cast< int >(j);
		streamlog_out(DEBUG5) << "ihit:" << ihit << std::endl;

		//if we are not in the last plane, call this method again
		if(i < _allHitsArray.size()-1)
		{
			vec.push_back( ihit); //index of the cluster in the last plane

			//track candidate requirements
			bool taketrack = true;
			const int e = vec.size()-2;
			if(e >= 0)
			{
				double residualX  = -999999.;
				double residualY  = -999999.;
				double residualZ  = -999999.;

				// ACTIVE
				// stop on the last non-zero hit

				// now loop through all hits on a track candidate "vec"
				// start at the end, stop on the first non-zero hit
				for(int ivec=e; ivec>=e; --ivec)                  //     <-> OFF
				{
					if(vec[ivec]>=0) // non zero hit has id vec[ivec]>=0 {otherwise -1}
					{
						double x = _allHitsArray[ivec][vec[ivec]].measuredX;
						double y = _allHitsArray[ivec][vec[ivec]].measuredY;
						double z = _allHitsArray[ivec][vec[ivec]].measuredZ;
						residualX  = abs(x - _allHitsArray[e+1][vec[e+1]].measuredX);
						residualY  = abs(y - _allHitsArray[e+1][vec[e+1]].measuredY);
						residualZ  = abs(z - _allHitsArray[e+1][vec[e+1]].measuredZ);
						break; 
					}   
				}

				if( 
						residualX < _residualsXMin[e] || residualX > _residualsXMax[e] ||
						residualY < _residualsYMin[e] || residualY > _residualsYMax[e] 
				  )
					taketrack = false;

				if( taketrack == false )
				{
					taketrack = true; 
					ihit=-1;
				} 
			}
			vec.pop_back(); 

			if(taketrack)
			{ 
				findtracks2(missinghits, indexarray, vec, _allHitsArray, i+1, ihit );
			}
		}
		else
		{
			//we are in the last plane
			vec.push_back( ihit ); //index of the cluster in the last plane

			//track candidate requirements
			bool taketrack = true;
			const int e = vec.size()-2;
			if(e >= 0)
			{
				double residualX  = -999999.;
				double residualY  = -999999.;
				//double residualZ  = -999999.;

				// now loop through all hits on a track candidate "vec"
				// start at the end, stop on the first non-zero hit
				for(int ivec=e; ivec>=e; --ivec)                        //   <-> OFF
				{
					if(vec[ivec]>=0) // non zero hit has id vec[ivec]>=0 {otherwise -1}
					{
						double x = _allHitsArray[ivec][vec[ivec]].measuredX;
						double y = _allHitsArray[ivec][vec[ivec]].measuredY;
						//double z = _allHitsArray[ivec][vec[ivec]].measuredZ;
						residualX  = abs(x - _allHitsArray[e+1][vec[e+1]].measuredX);
						residualY  = abs(y - _allHitsArray[e+1][vec[e+1]].measuredY);
						//residualZ  = abs(z - _allHitsArray[e+1][vec[e+1]].measuredZ);
						break; 
					}   
				}

				if( 
						residualX < _residualsXMin[e] || residualX > _residualsXMax[e] ||
						residualY < _residualsYMin[e] || residualY > _residualsYMax[e] 
				  )
					taketrack = false;

				if( taketrack == false )
				{
					taketrack = true; 
					ihit=-1;
				} 
			}

			if(static_cast< int >(indexarray.size()) >= _maxTrackCandidates)
				taketrack = false;

			if(taketrack)
			{
				indexarray.push_back(vec);
				streamlog_out(DEBUG9) << "indexarray size at last plane:" << indexarray.size() << std::endl;
			}
			vec.pop_back(); //last element must be removed because the
			//vector is still used -> we are in a last plane hit loop!

		}
	}

	if( (_allHitsArray[i].size() == 0) && (i >= _allHitsArray.size()-1) )
	{
		indexarray.push_back(vec);
	} 
}

void  EUTelMille2::findMatchedHits(int& _ntrack, Track* TrackHere)
{
	// hit list assigned to track
	std::vector<EVENT::TrackerHit*> TrackHitsHere = TrackHere->getTrackerHits();

	// check for a hit in every plane
	streamlog_out( MESSAGE1 ) << "track " << _ntrack << " has " << TrackHitsHere.size() << " hits " << endl;

	// assume hits are ordered in z! start counting from 0
	int nPlaneHere = 0;

	// setup cellIdDecoder to decode the hit properties
	CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

	std::vector<TrackerHit*> hit;
	std::vector<TrackerHit*> fit;

	// loop over all hits and fill arrays
	for(int nHits = 0; nHits < int(TrackHitsHere.size()); nHits++)
	{
		TrackerHit *HitHere = TrackHitsHere.at(nHits);
		int sensorID = Utility::getSensorIDfromHit ( HitHere ) ;

		// check if this is a measured hit or a fitted hit, want measured hit
		streamlog_out( MESSAGE0 ) << "hit on plane [" << sensorID << "] properties : " << ( hitCellDecoder(HitHere)["properties"] & kFittedHit ) << std::endl;

		if( ((hitCellDecoder(HitHere)["properties"] & kFittedHit) >> 1) == 0 ){  hit.push_back(HitHere); }
		if( ((hitCellDecoder(HitHere)["properties"] & kFittedHit) >> 1) == 1 ){  fit.push_back(HitHere); }

	}

	nPlaneHere = 0;

	for(std::vector<TrackerHit*>::iterator ihit = hit.begin(); ihit!= hit.end(); ihit++)
	{
		int hitID = Utility::getSensorIDfromHit ( (*ihit) );

		for(std::vector<TrackerHit*>::iterator ifit = fit.begin(); ifit!= fit.end(); ifit++)
		{
			int fitID = Utility::getSensorIDfromHit ( (*ifit) );          
			if( fitID != hitID ) continue;

			// hit positions
			const double *hitPosition = (*ihit)->getPosition();
			const double *fitPosition = (*ifit)->getPosition();

			// fill hits to arrays
			_xPos[_ntrack][nPlaneHere] = hitPosition[0] * 1000.;
			_yPos[_ntrack][nPlaneHere] = hitPosition[1] * 1000.;
			_zPos[_ntrack][nPlaneHere] = hitPosition[2] * 1000.;

			_trackResidX[_ntrack][nPlaneHere] = ( fitPosition[0] - hitPosition[0] ) * 1000. ;
			_trackResidY[_ntrack][nPlaneHere] = ( fitPosition[1] - hitPosition[1] ) * 1000. ;
			_trackResidZ[_ntrack][nPlaneHere] = ( fitPosition[2] - hitPosition[2] ) * 1000. ;

			nPlaneHere++;
		}
	}
	_ntrack++;
}

void EUTelMille2::processEvent(LCEvent * event)
{

	CellIDDecoder<TrackerHit>  hitDecoder(EUTELESCOPE::HITENCODING);

	if(_iEvt % 1000 == 0) 
	{
		streamlog_out( MESSAGE5 ) << "Currently having " << _nMilleDataPoints << " data points in " << _nMilleTracks << " tracks " << endl;
	}

	if( _nMilleTracks > _maxTrackCandidatesTotal )
	{
		return; //throw StopProcessingException(this);
	}
	EUTelEventImpl* evt = static_cast<EUTelEventImpl*>(event) ;

	if( evt->getEventType() == kEORE )
	{
		streamlog_out( DEBUG2 ) << "EORE found: nothing else to do." << std::endl;
		return;
	}

	IntVec indexconverter(_nPlanes,-1);
	std::vector<std::vector<EUTelMille2::HitsInPlane2> > _allHitsArray(_nPlanes, std::vector<EUTelMille2::HitsInPlane2>());

	int icounter = 0;

	for(size_t plane = 0; plane < _nPlanes; plane++)
	{
		int excluded = 0; //0 - not excluded, 1 - excluded
		if( _nExcludePlanes > 0 )
		{
			for(int exPlane = 0; exPlane < _nExcludePlanes; exPlane++)
			{
				if(plane == _excludePlanes[exPlane] ) {
					excluded = 1;
					break;//leave the for loop
				}
			}
		}
		if(excluded == 1)
		{
			indexconverter[plane] = -1;
		}
		else
		{
			indexconverter[plane] = icounter;
			icounter++;
		}
	}

	int _nTracks = 0;
	int _nGoodTracks = 0;

	if(_inputMode == 0)
	{
		for(size_t i =0;i < _hitCollectionName.size();i++)
		{
			LCCollection* collection;
			try
			{
				collection = event->getCollection(_hitCollectionName[i]);
			} 
			catch(DataNotAvailableException& e)
			{
				streamlog_out( WARNING2 ) << "No input collection " << _hitCollectionName[i] << " found for event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
				throw SkipEventException(this);
			}

			int layerIndex = -1;
			HitsInPlane2 hitsInPlane;

			// loop over all hits in collection
			for( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ )
			{
				TrackerHitImpl* hit = static_cast<TrackerHitImpl*>( collection->getElementAt(iHit) );

				int localSensorID = hitDecoder(hit)["sensorID"]; 
				layerIndex = _sensorIDVecMap[localSensorID] ;

				// Getting positions of the hits.
				hitsInPlane.measuredX = 1000. * hit->getPosition()[0];
				hitsInPlane.measuredY = 1000. * hit->getPosition()[1];
				hitsInPlane.measuredZ = 1000. * hit->getPosition()[2];

				_allHitsArray[layerIndex].push_back(hitsInPlane);
			} // end loop over all hits in collection
		}

		//perform simple track finding
		std::vector<IntVec > indexarray;

		findtracks2(0, indexarray, IntVec(), _allHitsArray, 0, 0);

		for(size_t i = 0; i < indexarray.size(); i++)
		{
			for(size_t j = 0; j <  _nPlanes; j++)
			{
				if( _allHitsArray[j].size()>0 &&  indexarray[i][j] >= 0 )
				{              
					_xPos[i][j] = _allHitsArray[j][indexarray[i][j]].measuredX;
					_yPos[i][j] = _allHitsArray[j][indexarray[i][j]].measuredY;
					_zPos[i][j] = _allHitsArray[j][indexarray[i][j]].measuredZ;
				}
				else
				{
					_xPos[i][j] = 0.;
					_yPos[i][j] = 0.;
					_zPos[i][j] = 0.;
				}  
			}
		}	
		_nTracks = static_cast<int>(indexarray.size());
	}

	//otherwise we take provided tracks!
	else if(_inputMode == 1)
	{

		LCCollection* collection;
		try
		{
			collection = event->getCollection(_trackCollectionName);
		} 
		catch (DataNotAvailableException& e)
		{
			streamlog_out( WARNING2 ) << "No input track collection " << _trackCollectionName  << " found for event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
			throw SkipEventException(this);
		}
		const int nTracksHere = collection->getNumberOfElements();

		// loop over all tracks
		for(int nTracksEvent = 0; nTracksEvent < nTracksHere && nTracksEvent < _maxTrackCandidates; nTracksEvent++)
		{
			Track *TrackHere = dynamic_cast<Track*>(collection->getElementAt(nTracksEvent));
			findMatchedHits( _nTracks, TrackHere );
		}
		streamlog_out( MESSAGE1 ) << "Number of tracks available in track collection: " << nTracksHere << " tracks selected for Mille: " << _nTracks << std::endl;
	}
       	
	if( _inputMode != 1 )
	{
		streamlog_out( MESSAGE1 ) << "Number of hits in the individual planes: ";
		for(size_t i = 0; i < _allHitsArray.size(); i++)
			streamlog_out( MESSAGE1 ) << _allHitsArray[i].size() << " ";
		streamlog_out( MESSAGE1 ) << endl;
	}

	streamlog_out( MESSAGE1 ) << "Number of track candidates found: " << _iEvt << ": " << _nTracks << endl;

	DoubleVec lambda;
	lambda.reserve(_nPlanes);
	bool validminuittrack = false;
	double Chiquare[2] = {0,0};
	double angle[2] = {0,0};

	// loop over all track candidates
	for(int track = 0; track < _nTracks; track++)
	{
		_xPosHere = new double[_nPlanes];
		_yPosHere = new double[_nPlanes];
		_zPosHere = new double[_nPlanes];

		for(unsigned int help = 0; help < _nPlanes; help++)
		{
			_xPosHere[help] = _xPos[track][help];
			_yPosHere[help] = _yPos[track][help];
			_zPosHere[help] = _zPos[track][help];

			if( _inputMode == 1 )
			{
				_waferResidX[help] = _trackResidX[track][help];
				_waferResidY[help] = _trackResidY[track][help];
				_waferResidZ[help] = _trackResidZ[track][help];
			} 
		}
		Chiquare[0] = 0.0;
		Chiquare[1] = 0.0;

		streamlog_out( MESSAGE1 ) << "Adding track using the following coordinates: ";

		{ 
		if(_alignMode == 3)
		{
				streamlog_out(MESSAGE1) << " AlignMode = " << _alignMode << " _inputMode = " << _inputMode << std::endl;

				//use minuit to find tracks

				size_t mean_n = 0  ;
				double mean_x = 0.0;
				double mean_y = 0.0;
				double mean_z = 0.0;
				int index_hitsarray=0;
				double x0 = -1.;
				double y0 = -1.;
				//double z0 = -1.;
				for(unsigned int help = 0; help < _nPlanes; help++) 
				{
					bool excluded = false;
					// check if actual plane is excluded
					if(_nExcludePlanes > 0) 
					{
						for(int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
						{
							if(help == _excludePlanes[helphelp]) 
							{
								excluded = true;
							}
						}
					}
					const double x = _xPos[track][help];
					const double y = _yPos[track][help];
					const double z = _zPos[track][help];
					if( abs(x)>1e-06 && abs(y)>1e-06  )
					{
						x0 = _xPos[track][help];
						y0 = _yPos[track][help];
						//z0 = _zPos[track][help];
					}
					const double xresid = x0 - x;
					const double yresid = y0 - y;
					streamlog_out( MESSAGE1 ) << " x0 = " << x0 << " x= " << x << " ;; y0 = " << y0 << " y = " << y << std::endl;

					if( xresid < _residualsXMin[help] || xresid > _residualsXMax[help]) 
					{
						continue;
					}
					if( yresid < _residualsYMin[help] || yresid > _residualsYMax[help]) 
					{
						continue;
					}

					if(!excluded)
					{
						double sigmax  = _resolutionX[help];
						double sigmay  = _resolutionY[help];
						double sigmaz  = _resolutionZ[help];

						if( !( abs(x)<1e-06 && abs(y)<1e-06) )
						{
							mean_z += z;
							mean_x += x;
							mean_y += y;
							mean_n++;  
						} 

						if( abs(x)<1e-06 && abs(y)<1e-06  )
						{
							sigmax = 1000000.;
							sigmay = 1000000.;
							sigmaz = 1000000.;
						}

						hitsarray2[index_hitsarray] = (hit2(x, y, z, sigmax, sigmay, sigmaz, help));
						index_hitsarray++;
					}
				}

				mean_z = mean_z / static_cast<double>(mean_n);
				mean_x = mean_x / static_cast<double>(mean_n);
				mean_y = mean_y / static_cast<double>(mean_n);

				int diff_mean = _nPlanes - mean_n;
				streamlog_out( MESSAGE0 ) << " diff_mean: " << diff_mean << " _nPlanes = " << _nPlanes << " mean_n = " << mean_n << std::endl;

				if( diff_mean > getAllowedMissingHits() ) 
				{
					continue;
				}

				static bool firstminuitcall = true;
				if(firstminuitcall)
				{
					gSystem->Load("libMinuit");//is this really needed?
					firstminuitcall = false;
				}

				TMinuit* gMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 4 params

				//set print level (-1 = quiet, 0 = normal, 1 = verbose)
				gMinuit->SetPrintLevel(-1);
				gMinuit->SetFCN(fcn_wrapper2);

				double arglist[10];
				int ierflg = 0;

				//minimization strategy (1 = standard, 2 = slower)
				arglist[0] = 2;
				gMinuit->mnexcm("SET STR",arglist,2,ierflg);

				//set error definition (1 = for chi square)
				arglist[0] = 1;
				gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

				//analytic track fit to guess the starting parameters
				double sxx = 0.0;
				double syy = 0.0;
				double szz = 0.0;

				double szx = 0.0;
				double szy = 0.0;

				for(size_t i = 0; i< number_of_datapoints2; i++)
				{
					const double x = hitsarray2[i].x;
					const double y = hitsarray2[i].y;
					const double z = hitsarray2[i].z;
					if( !(abs(x)<1e-06 && abs(y)<1e-06) )
					{
						sxx += pow(x-mean_x,2);
						syy += pow(y-mean_y,2);
						szz += pow(z-mean_z,2);

						szx += (x-mean_x)*(z-mean_z);
						szy += (y-mean_y)*(z-mean_z);
					}
				}
				double linfit_x_a1 = szx/szz; //slope
				double linfit_y_a1 = szy/szz; //slope

				double linfit_x_a0 = mean_x - linfit_x_a1 * mean_z; //offset
				double linfit_y_a0 = mean_y - linfit_y_a1 * mean_z; //offset

				double del= -1.0*atan(linfit_y_a1);//guess of delta
				double ps = atan(linfit_x_a1/sqrt(1.0+linfit_y_a1*linfit_y_a1));//guess of psi

				//Set starting values and step sizes for parameters
				Double_t vstart[4] = {linfit_x_a0, linfit_y_a0, del, ps};
				//duble vstart[4] = {0.0, 0.0, 0.0, 0.0};
				double step[4] = {0.01, 0.01, 0.01, 0.01};

				gMinuit->mnparm(0, "b0", vstart[0], step[0], 0,0,ierflg);
				gMinuit->mnparm(1, "b1", vstart[1], step[1], 0,0,ierflg);
				gMinuit->mnparm(2, "delta", vstart[2], step[2], -1.0*TMath::Pi(), 1.0*TMath::Pi(),ierflg);
				gMinuit->mnparm(3, "psi", vstart[3], step[3], -1.0*TMath::Pi(), 1.0*TMath::Pi(),ierflg);

				//Now ready for minimization step
				arglist[0] = 2000;
				arglist[1] = 0.01;
				gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);

				bool ok = true;

				if(ierflg != 0)
				{
					ok = false;            
				}

				//get results from migrad
				double b0 = 0.0;
				double b1 = 0.0;
				double delta = 0.0;
				double psi = 0.0;
				double b0_error = 0.0;
				double b1_error = 0.0;
				double delta_error = 0.0;
				double psi_error = 0.0;

				gMinuit->GetParameter(0,b0,b0_error);
				gMinuit->GetParameter(1,b1,b1_error);
				gMinuit->GetParameter(2,delta,delta_error);
				gMinuit->GetParameter(3,psi,psi_error);

				double c0 = 1.0;
				double c1 = 1.0;
				double c2 = 1.0;

				if(ok)
				{
					c0 = TMath::Sin(psi);
					c1 = -1.0*TMath::Cos(psi) * TMath::Sin(delta);
					c2 = TMath::Cos(delta) * TMath::Cos(psi);
					//cout << " b0: " << b0 << ", b1: " << b1 << ", c2: " << c2 << endl;
					validminuittrack = true;

					for(unsigned int help =0; help < _nPlanes; help++)
					{
						const double x = _xPos[track][help];
						const double y = _yPos[track][help];
						const double z = _zPos[track][help];

						//calculate the lambda parameter
						const double la = -1.0*b0*c0-b1*c1+c0*x+c1*y+sqrt(1-c0*c0-c1*c1)*z;
						lambda.push_back(la);

						//determine the residuals without reference vector
						_waferResidX[help] = b0 + la*c0 - x;
						_waferResidY[help] = b1 + la*c1 - y;
						_waferResidZ[help] = la*sqrt(1.0 - c0*c0 - c1*c1) - z;
						
						int sensorID =  _sensorIDVec[help];
						//std::cout << "Getting sensorID: " << sensorID << std::endl;
						Eigen::Vector3d xPrime = geo::gGeometry().globalXAxis(sensorID);
						Eigen::Vector3d yPrime = geo::gGeometry().globalYAxis(sensorID);
						Eigen::Vector3d offV = geo::gGeometry().getOffsetVector(sensorID)*1000;
					

						Eigen::Matrix3d leq;
						leq << 	xPrime(0), yPrime(0), -c0,
							xPrime(1), yPrime(1), -c1,
							xPrime(2), yPrime(2), -c2;
						
						Eigen::Vector3d leqV;
						leqV << b0-offV(0), b1-offV(1), -offV(2);

						Eigen::Vector3d result = leq.colPivHouseholderQr().solve(leqV);

						double resX = b0 + result(2)*c0 - x;
						//double resX1 = result(0)*xPrime(0)+result(1)*yPrime(0)+offV(0) -x;
						double resY = b1 + result(2)*c1 - y;
						//double resY1 = result(0)*xPrime(1)+result(1)*yPrime(1)+offV(1)-y;
						double resZ =  result(2)*c2-z;
						//double resZ1 = result(0)*xPrime(2)+result(1)*yPrime(2)+offV(2)-z;
				
						/*
						std::cout << "Residuals: " << std::endl;
						std::cout << "HitPos: " << x << ", " << y << ", " << z << std::endl;
						std::cout << "Old residual: " <<  _waferResidX[help] << ", " <<  _waferResidY[help] << ", " <<  _waferResidZ[help] << std::endl;
						std::cout << "New residual: " <<  resX << ", " <<  resY << ", " <<  resZ << std::endl;
						*/
						//std::cout << "New residual1: " <<  resX1 << ", " <<  resY1 << ", " <<  resZ1 << std::endl;
					
						_waferResidX[help] = resX;
						_waferResidY[help] = resY;
						_waferResidZ[help] = resZ;
						      
					}
				}
				delete gMinuit;
			}
			else
			{
				streamlog_out(MESSAGE1) << " AlignMode = " << _alignMode << "not supported" << std::endl;
			}
		}

		//Just printing out residual info if verbosity level is low enough
		streamlog_out( MESSAGE1 ) << "Residuals X: ";
		for(unsigned int help = 0; help < _nPlanes; help++)
		{
			streamlog_out( MESSAGE1 ) << _waferResidX[help] << " ";
		}
		streamlog_out( MESSAGE1 ) << endl;

		streamlog_out( MESSAGE1 ) << "Residuals Y: ";
		for(unsigned int help = 0; help < _nPlanes; help++) 
		{
			streamlog_out( MESSAGE1 ) << _waferResidY[help] << " ";
		}
		streamlog_out( MESSAGE1 ) << endl;

		streamlog_out( MESSAGE1 ) << "Residuals Z: ";
		for(unsigned int help = 0; help < _nPlanes; help++)
		{
			streamlog_out( MESSAGE1 ) << _waferResidZ[help] << " ";
		}
		streamlog_out( MESSAGE1 ) << endl;


		int residualsXOkay = 1;
		int residualsYOkay = 1;

		// check if residal cuts are used
		if(_useResidualCuts != 0) 
		{
			// loop over all sensors
			for(unsigned int help = 0; help < _nPlanes; help++) 
			{
				int excluded = 0; //0 not excluded, 1 excluded
				if(_nExcludePlanes > 0) 
				{
					for(int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
					{
						if(help == _excludePlanes[helphelp]) 
						{
							excluded = 1;
						}
					}
				}
				if(excluded == 0)
				{
					if(_waferResidX[help] < _residualsXMin[help] || _waferResidX[help] > _residualsXMax[help]) 
					{
						residualsXOkay = 0;
					}
					if(_waferResidY[help] < _residualsYMin[help] || _waferResidY[help] > _residualsYMax[help]) 
					{
						residualsYOkay = 0;
					}
				}

			} // end loop over all sensors

		} // end check if residual cuts are used

		if(_useResidualCuts != 0 && (residualsXOkay == 0 || residualsYOkay == 0)) {
			streamlog_out( MESSAGE1 ) << "Track did not pass the residual cuts." << endl;
		}

		// apply track cuts (at the moment only residuals)
		if(_useResidualCuts == 0 || (residualsXOkay == 1 && residualsYOkay == 1))
		{
			// Add track to Millepede
			if(_alignMode == 3)
			{
				if(validminuittrack || _inputMode == 1 )
				{
					const int nLC = 4; // number of local parameters
					const int nGL = _nPlanes * 6; // number of global parameters

					float* derLC = new float[nLC]; // array of derivatives for local parameters
					float* derGL = new float[nGL]; // array of derivatives for global parameters

					int* label = new int[nGL]; // array of labels

					float residual;

					// create labels
					for(int gPara = 0; gPara < nGL; gPara++) 
					{
						label[gPara] = gPara + 1;
						derGL[gPara] = 0;
					}

					for(int lPara = 0; lPara < nLC; lPara++) 
					{
						derLC[lPara] = 0;
					}

					int nExcluded = 0;

					// loop over all planes
					for(unsigned int help = 0; help < _nPlanes; help++) 
					{
						int excluded = 0;

						double sigmax = _resolutionX[help];
						double sigmay = _resolutionY[help];
						double sigmaz = _resolutionZ[help];

						if( 
								abs(_xPosHere[help]) < 1e-06 &&
								abs(_yPosHere[help]) < 1e-06 &&
								abs(_zPosHere[help]) < 1e-06   
						  )
						{
							sigmax *= 1000000.;
							sigmay *= 1000000.;
							sigmaz *= 1000000.;
							continue;
						}

						if(excluded == 0)
						{
							//     cout << "--" << endl;
							int helphelp = help - nExcluded; // index of plane after
							// excluded planes have
							// been removed

							//local parameters: b0, b1, c0, c1
							//                  const double la = lambda[help];

							//                  double z_sensor = _siPlanesLayerLayout -> getSensitivePositionZ(help) + 0.5 * _siPlanesLayerLayout->getSensitiveThickness( help );
							//                  z_sensor *= 1000;		// in microns
							// reset all derivatives to zero!
							for(int i = 0; i < nGL; i++ ) 
							{
								derGL[i] = 0.000;
							}

							for(int i = 0; i < nLC; i++ ) 
							{
								derLC[i] = 0.000;
							}

							double x_sensor = 0.;
							double y_sensor = 0.;
							double z_sensor = 0.;

							x_sensor *= 1000.;
							y_sensor *= 1000.;
							z_sensor *= 1000.;

							// track model : fit-reco => 
							//   (a_X*x+b_X, a_Y*y+b_Y)  :   /  1   -g    b \   / x          \   :: shouldn't it be x-xcenter-of-the-sensor ??
							//                           : - |  g    1   -a |   | y          |   ::  and y-ycenter-of-the-sensor ??   
							//                           :   \ -b    a    1 /   \ z-z_sensor /   ::  == z-zcenter-of-the-sensor  ?? (already)
							// make angles sings consistent with X->Y->Z->X rotations.
							// correct likewise all matrices in ApplyAlignment processor
							// Igor Rubinsky 09-10-2011
							
							
							//Track can be computed as following: r_global = R*r_local + o(ffset)
							//Where R = R_y(beta)*R_x(alpha)*R_z(gamma)
							//This yields, using small angle approx.:
							//  r_x = -x-gamma*y+beta*z+o_x
							//  r_y = gamma*x+y-alpha*z+o_y
							//  r_z = -beta*x+alpha*y+z+o_z

							//GLOBAL DERIVATIVES
							//r_x derivatives (o_x, o_y, o_z, alpha, beta, gamma)
							derGL[((helphelp * 6) + 0)] = -1.0;				//d(r_x)/d(o_x)
							derGL[((helphelp * 6) + 1)] =  0.0;				//d(r_x)/d(o_y)
							derGL[((helphelp * 6) + 2)] =  0.0;				//d(r_x)/d(o_z)
							derGL[((helphelp * 6) + 3)] =  0.0;				//d(r_x)/d(alpha)
							derGL[((helphelp * 6) + 4)] = -1.0*(_zPosHere[help]);		//d(r_x)/d(beta)
							derGL[((helphelp * 6) + 5)] =  1.0*(_yPosHere[help]);		//d(r_x)/d(gamma)
							//LOCAL DERIVATIVES
							//r_x = b1+tan(c1)*z ~~ b1+c1*z
							derLC[0] = 1.0;							//d(r_x)/d(b0)
							derLC[1] = 0.0;							//d(r_x)/d(b1)
							derLC[2] = _zPosHere[help] + _waferResidZ[help];		//d(r_x)/d(c0)
							derLC[3] = 0.0;							//d(r_x)/d(c1)

							residual = _waferResidX[help];
							_mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmax);

							//same for r_y
							derGL[((helphelp * 6) + 0)] =  0.0;
							derGL[((helphelp * 6) + 1)] = -1.0; 
							derGL[((helphelp * 6) + 2)] =  0.0;
							derGL[((helphelp * 6) + 3)] =  1.0*(_zPosHere[help] - z_sensor); 
							derGL[((helphelp * 6) + 4)] =  0.0; 
							derGL[((helphelp * 6) + 5)] = -1.0*(_xPosHere[help] - x_sensor);
				
							derLC[0] = 0.0;
							derLC[1] = 1.0;
							derLC[2] = 0.0;
							derLC[3] = _zPosHere[help] + _waferResidZ[help];

							residual = _waferResidY[help];
							_mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmay);

							//same for r_z
							derGL[((helphelp * 6) + 0)] =  0.0;
							derGL[((helphelp * 6) + 1)] =  0.0;
							derGL[((helphelp * 6) + 2)] = -1.0; 
							derGL[((helphelp * 6) + 3)] = -1.0*(_yPosHere[help]-y_sensor); 
							derGL[((helphelp * 6) + 4)] =  1.0*(_xPosHere[help]-x_sensor);
							derGL[((helphelp * 6) + 5)] =  0.0;
							//r_z = c1*x + c2*y ?? TODO
							derLC[0] = 0.0;
							derLC[1] = 0.0;
							derLC[2] = _xPosHere[help] + _waferResidX[help];
							derLC[3] = _yPosHere[help] + _waferResidY[help];

							residual = _waferResidZ[help];

							_mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmaz);
							_nMilleDataPoints++;
						}// end if plane is not excluded
					} // end loop over all planes
					
					// clean up
					delete[] derLC;
					delete[] derGL;
					delete[] label;
				}
			} 
			else	
			{
				streamlog_out( ERROR2 ) << _alignMode << " is not a valid mode. Please choose 1,2 or 3." << endl;
			}

			_nGoodTracks++;

			// end local fit
			_mille->end();
			_nMilleTracks++;
		} // end if apply track cuts

		// clean up
		delete[] _zPosHere;
		delete[] _yPosHere;
		delete[] _xPosHere;

	}//end loop over all track candidates

	streamlog_out( MESSAGE1 ) << "Finished fitting tracks in event " << _iEvt << endl;

	// count events
	++_iEvt;
	if( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelMille2::end()
{
	delete[] _yFitPos;
	delete[] _xFitPos;
	delete[] _waferResidY;
	delete[] _waferResidX;
	delete[] _waferResidZ;

	if(_alignMode == 3)
	{
		delete[] hitsarray2;
	}

	// close the output file
	delete _mille;

	// if write the pede steering file
	if(_generatePedeSteerfile)
	{
		streamlog_out( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

		double* meanX = new double[_nPlanes];
		double* meanY = new double[_nPlanes];
		double* meanZ = new double[_nPlanes];

		ofstream steerFile;
		steerFile.open(_pedeSteerfileName.c_str());

		if(steerFile.is_open())
		{

			// find first and last excluded plane
			unsigned int firstnotexcl = _nPlanes;
			unsigned int lastnotexcl = 0;

			// loop over all planes
			for(unsigned int help = 0; help < _nPlanes; help++) 
			{

				int excluded = 0;

				// loop over all excluded planes
				for(int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
				{
					if(help == _excludePlanes[helphelp]) 
					{
						//            excluded = 1;
					}
				} // end loop over all excluded planes

				if(excluded == 0 && firstnotexcl > help) 
				{
					firstnotexcl = help;
				}

				if(excluded == 0 && lastnotexcl < help) 
				{
					lastnotexcl = help;
				}
			} // end loop over all planes

			// calculate average
			double averageX = (meanX[firstnotexcl] + meanX[lastnotexcl]) / 2;
			double averageY = (meanY[firstnotexcl] + meanY[lastnotexcl]) / 2;
			double averageZ = (meanZ[firstnotexcl] + meanZ[lastnotexcl]) / 2;

			steerFile << "Cfiles" << endl;
			steerFile << _binaryFilename << endl;
			steerFile << endl;

			steerFile << "Parameter" << endl;

			int counter = 0;

			// loop over all planes
			for(unsigned int help = 0; help < _nPlanes; help++) {

				int excluded = 0; // flag for excluded planes

				// loop over all excluded planes
				for(int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {

					if(help == _excludePlanes[helphelp]) {
						//            excluded = 1;
					}

				} // end loop over all excluded planes

				// if plane not excluded
				if(excluded == 0) {

					bool fixed = false;
					for(size_t i = 0;i< _FixedPlanes.size(); i++)
					{
						if(_FixedPlanes[i] == static_cast< int >(help))
							fixed = true;
					}

					if( fixed || (_FixedPlanes.empty() && (help == firstnotexcl || help == lastnotexcl) ) )
					{
						if(_alignMode == 3) {
							steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
							steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
							steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
							steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
							steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
							steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
						}

					} else {

						if(_alignMode == 3) {
							if(_usePedeUserStartValues == 0)
							{
								if(_FixParameter[help] & (1 << 0))
									steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 1) << " " << (averageX - meanX[help]) << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 1))
									steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 2) << " " << (averageY - meanY[help]) << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 2))
									steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 3) << " " << (averageZ - meanZ[help]) << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 3))
									steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 4) << " 0.0 0.0" << endl;

								if(_FixParameter[help] & (1 << 4))
									steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 5) << " 0.0 0.0" << endl;

								if(_FixParameter[help] & (1 << 5))
									steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 6) << " 0.0 0.0" << endl;
							}
							else
							{
								if(_FixParameter[help] & (1 << 0))
									steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 1))
									steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 2))
									steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 3) << " " << _pedeUserStartValuesZ[help] << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 3))
									steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 4) << " " << _pedeUserStartValuesAlpha[help] << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 4))
									steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 5) << " " << _pedeUserStartValuesBeta[help] << " 0.0" << endl;

								if(_FixParameter[help] & (1 << 5))
									steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
								else
									steerFile << (counter * 6 + 6) << " " << _pedeUserStartValuesGamma[help] << " 0.0" << endl;
							}
						}
					}

					counter++;

				} // end if plane not excluded

			} // end loop over all planes

			steerFile << endl;
			for(StringVec::iterator it = _pedeSteerAddCmds.begin(); it != _pedeSteerAddCmds.end(); ++it)
			{
				// two backslashes will be interpreted as newline
				if(*it == "\\\\")
					steerFile << endl;
				else
					steerFile << *it << " ";
			}
			steerFile << endl;
			steerFile << "method inversion 10 0.001" << endl;
			steerFile << endl;
			steerFile << "histprint" << endl;
			steerFile << endl;
			steerFile << "end" << endl;

			steerFile.close();

			streamlog_out( MESSAGE5 ) << "File " << _pedeSteerfileName << " written." << endl;

		}
		else
		{
			streamlog_out( ERROR2 ) << "Could not open steering file." << endl;
		}

		//cleaning up
		delete[] meanX;
		delete[] meanY;
		delete[] meanZ;
	} // end if write the pede steering file

	streamlog_out( MESSAGE7 ) << "Number of data points used: " << _nMilleDataPoints << endl;
	streamlog_out( MESSAGE7 ) << "Number of tracks used: " << _nMilleTracks << endl;

	streamlog_out( MESSAGE2 ) << endl;
	streamlog_out( MESSAGE2 ) << "Successfully finished" << endl;
}
#endif // USE_GEAR
