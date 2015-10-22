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
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelReferenceHit.h"
#include "EUTelCDashMeasurement.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "mille/Mille.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TRandom.h>
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
//std::vector<EUTelMille2::hit> hitsarray;
EUTelMille2::hit *hitsarray2;
unsigned int number_of_datapoints2;
void fcn_wrapper2(int& /*npar*/, double* /*gin*/, double &f, double *par, int /*iflag*/)
{
  EUTelMille2::trackfitter fobj(hitsarray2,number_of_datapoints2);
  f = fobj.fit(par);
}


// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMille2::_numberTracksLocalname   = "NumberTracks";
std::string EUTelMille2::_chi2XLocalname          = "Chi2X";
std::string EUTelMille2::_chi2YLocalname          = "Chi2Y";
std::string EUTelMille2::_residualXLocalname      = "ResidualX";
std::string EUTelMille2::_residualYLocalname      = "ResidualY";
std::string EUTelMille2::_residualZLocalname      = "ResidualZ";

std::string EUTelMille2::_residualXvsXLocalname      = "ResidualXvsX";
std::string EUTelMille2::_residualXvsYLocalname      = "ResidualXvsY";
std::string EUTelMille2::_residualYvsXLocalname      = "ResidualYvsX";
std::string EUTelMille2::_residualYvsYLocalname      = "ResidualYvsY";
std::string EUTelMille2::_residualZvsXLocalname      = "ResidualZvsX";
std::string EUTelMille2::_residualZvsYLocalname      = "ResidualZvsY";


#endif




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
  _description =
    "EUTelMille2 uses the MILLE program to write data files for MILLEPEDE II.";

  // choose input mode
  registerOptionalParameter("InputMode","Selects the source of input hits."
                            "\n0 - hits read from hitfile and simple straight line trackfinding will be  performed internally. "
                            "\n1 - hits read from output of tracking processor. "
                            "\n2 - Test mode. Simple internal simulation and simple trackfinding. "
                            "\n3 - Mixture of a track collection from the telescope and hit collections for the DUT (only one DUT layer can be used unfortunately)",
                            _inputMode, static_cast <int> (0));

  registerOptionalParameter("AllowedMissingHits","Set how many hits (=planes) can be missing on a track candidate.",
                            _allowedMissingHits, static_cast <int> (0));

  registerOptionalParameter("MimosaClusterChargeMin","Remove Mimosa26 clusters with a charge (i.e. number of fired pixels in cluster) below or equal to this value",
                            _mimosa26ClusterChargeMin,  static_cast <int> (1) );


  // input collections
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("corrhits");

  registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName",
                           "Hit collections name",
                           _hitCollectionName,HitCollectionNameVecExample);

  registerInputCollection(LCIO::TRACK,"TrackCollectionName",
                          "Track collection name. This is only relevant if InputMode is set to larger to 1",
                          _trackCollectionName,std::string("fittracks"));

  // parameters

  registerOptionalParameter("DistanceMax","Maximal allowed distance between hits entering the fit per 10 cm space between the planes.",
                            _distanceMax, static_cast <float> (2000.0));

  registerOptionalParameter("DistanceMaxVec","Maximal allowed distance between hits entering the fit per 10 cm space between the planes. One value for each neighbor planes. DistanceMax will be used for each pair if this vector is empty.",
                            _distanceMaxVec, FloatVec ());

  

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

  registerOptionalParameter("UseResidualCuts","Use cuts on the residuals to reduce the combinatorial background.",_useResidualCuts,
                            static_cast <bool> (false));

  registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
                            "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );

  registerOptionalParameter("AlignmentConstantCollectionName", "This is the name of the alignment collection to be saved into the slcio file",
                            _alignmentConstantCollectionName, static_cast< string > ( "alignment" ));

  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin,MinimalResidualsX);

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin,MinimalResidualsY);

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax,MaximalResidualsX);

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax,MaximalResidualsY);



  registerOptionalParameter("ResolutionX","X resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionX,  FloatVec (static_cast <int> (6), 10.));

  registerOptionalParameter("ResolutionY","Y resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionY,FloatVec (static_cast <int> (6), 10.));

  registerOptionalParameter("ResolutionZ","Z resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionZ,FloatVec (static_cast <int> (6), 10.));

  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("referenceHit") );
 
  registerOptionalParameter("UseReferenceCollection","Do you want the reference hit collection to be used for coordinate transformations?",  _useReferenceHitCollection, static_cast< bool   > ( true ));

  registerOptionalParameter("FixParameter","Fixes the given alignment parameters in the fit if alignMode==3 is used. For each sensor an integer must be specified (If no value is given, then all parameters will be free). bit 0 = x shift, bit 1 = y shift, bit 2 = z shift, bit 3 = alpha, bit 4 = beta, bit 5 = gamma. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_FixParameter, IntVec (static_cast <int> (6), 24));

  registerOptionalParameter("GeneratePedeSteerfile","Generate a steering file for the pede program.",_generatePedeSteerfile, static_cast <int> (0));

  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));

  registerOptionalParameter("PedeSteeringAdditionalCmds","FOR EXPERTS: List of commands that should be included in the pede steering file. Use '\\' to seperate options and introduce a line break.",_pedeSteerAddCmds, StringVec());

  registerOptionalParameter("RunPede","Execute the pede program using the generated steering file.",_runPede, static_cast <bool> (true));

  registerOptionalParameter("UsePedeUserStartValues","Give start values for pede by hand (0 - automatic calculation of start values, 1 - start values defined by user).", _usePedeUserStartValues, static_cast <int> (0));

  registerOptionalParameter("PedeUserStartValuesX","Start values for the alignment for shifts in the X direction.",_pedeUserStartValuesX,PedeUserStartValuesX);

  registerOptionalParameter("PedeUserStartValuesY","Start values for the alignment for shifts in the Y direction.",_pedeUserStartValuesY,PedeUserStartValuesY);

  registerOptionalParameter("PedeUserStartValuesZ","Start values for the alignment for shifts in the Z direction.",_pedeUserStartValuesZ,FloatVec (static_cast <int> (6), 0.0));

  registerOptionalParameter("PedeUserStartValuesAlpha","Start values for the alignment for the angle alpha.",_pedeUserStartValuesAlpha,FloatVec (static_cast <int> (6), 0.0));
  
  registerOptionalParameter("PedeUserStartValuesBeta","Start values for the alignment for the angle beta.",_pedeUserStartValuesBeta,FloatVec (static_cast <int> (6), 0.0));
  
  registerOptionalParameter("PedeUserStartValuesGamma","Start values for the alignment for the angle gamma.",_pedeUserStartValuesGamma,PedeUserStartValuesGamma);

  registerOptionalParameter("TestModeSensorResolution","Resolution assumed for the sensors in test mode.",_testModeSensorResolution, static_cast <float> (3.0));

  registerOptionalParameter("TestModeXTrackSlope","Width of the track slope distribution in the x direction",_testModeXTrackSlope, static_cast <float> (0.0005));

  registerOptionalParameter("TestModeYTrackSlope","Width of the track slope distribution in the y direction",_testModeYTrackSlope, static_cast <float> (0.0005));

  registerOptionalParameter("TestModeSensorZPositions","Z positions of the sensors in test mode.",_testModeSensorZPositions,SensorZPositions);

  registerOptionalParameter("TestModeSensorXShifts","X shifts of the sensors in test mode (to be determined by the alignment).",
                            _testModeSensorXShifts,SensorXShifts);

  registerOptionalParameter("TestModeSensorYShifts","Y shifts of the sensors in test mode (to be determined by the alignment).",
                            _testModeSensorYShifts,SensorYShifts);


  registerOptionalParameter("TestModeSensorGamma","Rotation around the z axis of the sensors in test mode (to be determined by the alignment).",
                            _testModeSensorGamma,SensorGamma);


  registerOptionalParameter("TestModeSensorAlpha","Rotation around the x axis of the sensors in test mode (to be determined by the alignment).",
                            _testModeSensorAlpha,SensorAlpha);


  registerOptionalParameter("TestModeSensorBeta","Rotation around the y axis of the sensors in test mode (to be determined by the alignment).",
                            _testModeSensorBeta,SensorBeta);

  IntVec initRect;
  registerOptionalParameter("UseSensorRectangular","Do not use all pixels for alignment, only these in the rectangular (A|B) e.g. (0,0) and (C|D) e.g. (100|100) of sensor S. Type in the way S1 A1 B1 C1 D1 S2 A2 B2 C2 D2 ...",
                            _useSensorRectangular,initRect);

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection to be saved into the output slcio file",
                             _hotPixelCollectionName, static_cast< string > ( "" ));

}

void EUTelMille2::init() {

    // Getting access to geometry description
    std::string name("test.root");
    geo::gGeometry().initializeTGeoDescription(name,false);

  _sensorIDVec.clear();
  _sensorIDVecMap.clear();
	//TODO: get this directly

  // an associative map for getting also the sensorID ordered
  map< double, int > sensorIDMap;

  for ( size_t i = 0; i < geo::gGeometry().sensorIDsVec().size(); i++) 
   {
	int sensorID = geo::gGeometry().sensorIDsVec().at(i);
    	_sensorIDVec.push_back( sensorID );
    	_sensorIDVecMap.insert( make_pair( sensorID, i ) );
    	sensorIDMap.insert( make_pair( geo::gGeometry().siPlaneZPosition(sensorID), sensorID ) );
   }
     
  _histogramSwitch = true;
  _referenceHitVec = NULL;
  
  _nPlanes = geo::gGeometry().nPlanes();

  //the user is giving sensor ids for the planes to be excluded. this
  //sensor ids have to be converted to a local index according to the
  //planes positions along the z axis.
  for (size_t i = 0; i < _FixedPlanes_sensorIDs.size(); i++)
    {
      map< double, int >::iterator iter = sensorIDMap.begin();
      int counter = 0;
      while ( iter != sensorIDMap.end() ) {
        if( iter->second == _FixedPlanes_sensorIDs[i])
          {
            _FixedPlanes.push_back(counter);
            break;
          }
        ++iter;
        ++counter;
      }
    }
  for (size_t i = 0; i < _excludePlanes_sensorIDs.size(); i++)
    {
      map< double, int >::iterator iter = sensorIDMap.begin();
      int counter = 0;
      while ( iter != sensorIDMap.end() ) {
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
  while ( iter != sensorIDMap.end() ) {
    bool excluded = false;
    for (size_t i = 0; i < _excludePlanes.size(); i++)
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

  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // Initialize number of excluded planes
  _nExcludePlanes = _excludePlanes.size();

  streamlog_out ( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

  // Initialise Mille statistics
  _nMilleDataPoints = 0;
  _nMilleTracks = 0;

  _waferResidX = new double[_nPlanes];
  _waferResidY = new double[_nPlanes];
  _waferResidZ = new double[_nPlanes];
  
  _xFitPos = new double[_nPlanes];
  _yFitPos = new double[_nPlanes];

  _telescopeResolX = new double[_nPlanes];
  _telescopeResolY = new double[_nPlanes];
  _telescopeResolZ = new double[_nPlanes];

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
          streamlog_out ( WARNING2 ) << "Consistency check of the resolution parameters failed. The array size is different than the number of found planes! "
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
         streamlog_out ( WARNING2 ) << "Consistency check of the fixed parameters array failed. The array size is different than the number of found planes! The array is now set to default values, which means that all parameters are free in the fit." << endl;
          _FixParameter.clear();
          for(size_t i = 0; i < _nPlanes; i++)
            {
              _FixParameter.push_back(0);
            }
        }
      if(_FixParameter.empty())
        {
          streamlog_out ( WARNING2 ) << "The fixed parameters array was found to be empty. It will be filled with default values. All parameters are free in the fit now." << endl;
          _FixParameter.clear();
          for(size_t i = 0; i < _nPlanes; i++)
            {
              _FixParameter.push_back(0);
            }
        }
    }

  if(_alignMode == 3)
    {
      number_of_datapoints2 = _nPlanes;
      hitsarray2 = new hit[number_of_datapoints2];
    }

  // booking histograms
  bookHistos();

  streamlog_out ( MESSAGE5 ) << "Initialising Mille..." << endl;
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
          streamlog_out ( WARNING2 ) << "Consistency check of the DistanceMaxVec array failed. Its size is different compared to the number of planes! Will now use _distanceMax for each pair of planes." << endl;
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

  streamlog_out ( MESSAGE4 ) << "end of initialisation" << endl;
}

void EUTelMille2::processRunHeader (LCRunHeader * rdr) {

  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if ( (unsigned int)header->getGeoID() != geo::gGeometry().getSiPlanesLayoutID() ) {
    streamlog_out ( ERROR2 ) << "Error during the geometry consistency check: " << endl;
    streamlog_out ( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << endl;
    streamlog_out ( ERROR2 ) << "The GEAR description says is     " << geo::gGeometry().getSiPlanesLayoutID() << endl;
  }

  // increment the run counter
  ++_iRun;
}

void  EUTelMille2::findMatchedHits(int& _ntrack, Track* TrackHere) {
 
      // hit list assigned to track
      std::vector<EVENT::TrackerHit*> TrackHitsHere = TrackHere->getTrackerHits();

      // check for a hit in every plane
      streamlog_out ( MESSAGE1 ) << "track " << _ntrack << " has " << TrackHitsHere.size() << " hits " << endl;

      // assume hits are ordered in z! start counting from 0
      int nPlaneHere = 0;

      // setup cellIdDecoder to decode the hit properties
      CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);



      std::vector <TrackerHit* > hit;
      std::vector <TrackerHit* > fit;

      // loop over all hits and fill arrays
      for (int nHits = 0; nHits < int(TrackHitsHere.size()); nHits++) {

          TrackerHit *HitHere = TrackHitsHere.at(nHits);
          int sensorID = Utility::getSensorIDfromHit ( HitHere ) ;

          // check if this is a measured hit or a fitted hit, want measured hit
          streamlog_out( MESSAGE0 ) << "hit on plane [" << sensorID << "] properties : " << ( hitCellDecoder(HitHere)["properties"] & kFittedHit ) << std::endl;

          if ( ((hitCellDecoder(HitHere)["properties"] & kFittedHit) >> 1) == 0 ){  hit.push_back(HitHere); }

          if ( ((hitCellDecoder(HitHere)["properties"] & kFittedHit) >> 1) == 1 ){  fit.push_back(HitHere); }

      }

      nPlaneHere = 0;

      for( std::vector<TrackerHit*>::iterator ihit = hit.begin(); ihit!= hit.end(); ihit++){
         int hitID = Utility::getSensorIDfromHit ( (*ihit) );

         streamlog_out ( MESSAGE1 ) << "hit: @ " << hitID << " " << std::endl;

         for( std::vector<TrackerHit*>::iterator ifit = fit.begin(); ifit!= fit.end(); ifit++){
             int fitID = Utility::getSensorIDfromHit ( (*ifit) );          
             streamlog_out ( MESSAGE1 ) << "fit: @ " << fitID << " " << std::endl;
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
   
            streamlog_out ( MESSAGE1 ) << "hit: @ " << hitID << " " 
                          << _xPos[_ntrack][nPlaneHere] << " " 
                          << _yPos[_ntrack][nPlaneHere] << " " 
                          << _zPos[_ntrack][nPlaneHere] << " type: " << (*ihit)->getType() ;

            streamlog_out ( MESSAGE1 ) << "res " 
                         << _trackResidX[_ntrack][nPlaneHere] << " " 
                         << _trackResidY[_ntrack][nPlaneHere] << " " 
                         << _trackResidZ[_ntrack][nPlaneHere] ;
            streamlog_out ( MESSAGE1 ) << endl;


            nPlaneHere++;
             

         }
      }

      _ntrack++;

}

void EUTelMille2::processEvent (LCEvent * event) {

  CellIDDecoder<TrackerHit>  hitDecoder(EUTELESCOPE::HITENCODING);

  if ( _useReferenceHitCollection ){
    try {
    _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
    }
    catch (...){
      streamlog_out ( ERROR5 ) <<  "Reference Hit Collection " << _referenceHitCollectionName.c_str() << " could not be retrieved for event " << event->getEventNumber()<< "! Please check your steering files! " << endl;
    }
  }
  
  if (_iEvt % 1000 == 0) 
  {
    streamlog_out( MESSAGE5 ) << "Currently having " << _nMilleDataPoints << " data points in "
                              << _nMilleTracks << " tracks " << endl;
  }
  
  if( _nMilleTracks > _maxTrackCandidatesTotal )
  {
	return; //throw StopProcessingException(this);
  }
  
  // fill resolution arrays
  for (size_t help = 0; help < _nPlanes; help++) {
    _telescopeResolX[help] = _telescopeResolution;
    _telescopeResolY[help] = _telescopeResolution;
  }

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  std::vector<std::vector<EUTelMille2::HitsInPlane> > _hitsArray(_nPlanes - _nExcludePlanes, std::vector<EUTelMille2::HitsInPlane>());
  IntVec indexconverter (_nPlanes,-1);
 
  std::vector<std::vector<EUTelMille2::HitsInPlane> > _allHitsArray(_nPlanes, std::vector<EUTelMille2::HitsInPlane>());
  
  {
 
    
      int icounter = 0;
      for(size_t i = 0; i < _nPlanes; i++)
      {
          int excluded = 0; //0 - not excluded, 1 - excluded
          if ( _nExcludePlanes > 0 )
          {
              for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
                  if (i == _excludePlanes[helphelp] ) {
                      excluded = 1;
                      break;//leave the for loop
                  }
              }
          }
          if(excluded == 1)
              indexconverter[i] = -1;
          else
          {
              indexconverter[i] = icounter;
              icounter++;
          }
       }
  }

  int _nTracks = 0;

  int _nGoodTracks = 0;

  if (_inputMode == 1) {

    LCCollection* collection;
    try {
      collection = event->getCollection(_trackCollectionName);
    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input track collection " << _trackCollectionName  << " found for event " << event->getEventNumber()
                                 << " in run " << event->getRunNumber() << endl;
      throw SkipEventException(this);
    }
    const int nTracksHere = collection->getNumberOfElements();


    // loop over all tracks
    for (int nTracksEvent = 0; nTracksEvent < nTracksHere && nTracksEvent < _maxTrackCandidates; nTracksEvent++) {

      Track *TrackHere = dynamic_cast<Track*> (collection->getElementAt(nTracksEvent));

      findMatchedHits( _nTracks, TrackHere );

    } // end loop over all tracks

    streamlog_out ( MESSAGE1 ) << "Number of tracks available in track collection: " << nTracksHere << " tracks selected for Mille: " << _nTracks << std::endl;

  } 
  
  streamlog_out ( MESSAGE1 ) << "Number of track candidates found: " << _iEvt << ": " << _nTracks << endl;

  // Perform fit for all found track candidates
  // ------------------------------------------

  // only one track or no single track event
//  if (_nTracks == 1 || _onlySingleTrackEvents == 0) 
  {

    DoubleVec lambda;
    lambda.reserve(_nPlanes);
    bool validminuittrack = false;

    double Chiquare[2] = {0,0};
    double angle[2] = {0,0};

    // loop over all track candidates
    for (int track = 0; track < _nTracks; track++) {

      _xPosHere = new double[_nPlanes];
      _yPosHere = new double[_nPlanes];
      _zPosHere = new double[_nPlanes];


      for (unsigned int help = 0; help < _nPlanes; help++) {
        _xPosHere[help] = _xPos[track][help];
        _yPosHere[help] = _yPos[track][help];
        _zPosHere[help] = _zPos[track][help];

        if( _inputMode == 1 ) {
          _waferResidX[help] = _trackResidX[track][help];
          _waferResidY[help] = _trackResidY[track][help];
          _waferResidZ[help] = _trackResidZ[track][help];
        } 
      }

      Chiquare[0] = 0.0;
      Chiquare[1] = 0.0;

      streamlog_out ( MESSAGE1 ) << "Adding track using the following coordinates: ";

      // loop over all planes
      for (unsigned int help = 0; help < _nPlanes; help++) 
      {
        int excluded = 0;

        // check if actual plane is excluded
        if (_nExcludePlanes > 0) {
          for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {
            if (help == _excludePlanes[helphelp]) {
              excluded = 1;
            }
          }
        }

        if (excluded == 0) {
          streamlog_out ( MESSAGE1 ) << 
                std::endl << " not Excluded @ " << help  << "["<<_nPlanes << "] " <<_xPosHere[help] << " " << _yPosHere[help] << " " << _zPosHere[help] ;
        }
        streamlog_out ( MESSAGE1 ) << std::endl;

      } // end loop over all planes

      streamlog_out ( MESSAGE1 ) << endl;
     
//    if( _inputMode == 1 ) {
//    }else
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
          for (unsigned int help = 0; help < _nPlanes; help++) 
          {
            bool excluded = false;
            // check if actual plane is excluded
            if (_nExcludePlanes > 0) 
            {
              for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
              {
                if (help == _excludePlanes[helphelp]) 
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
            streamlog_out ( MESSAGE1 ) << " x0 = " << x0 << " x= " << x << " ;; y0 = " << y0 << " y = " << y << std::endl;

            if ( xresid < _residualsXMin[help] || xresid > _residualsXMax[help]) 
              {
                continue;
              }
            if ( yresid < _residualsYMin[help] || yresid > _residualsYMax[help]) 
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

                hitsarray2[index_hitsarray] = (hit(
                                                  x, y, z,
                                                  sigmax, sigmay, sigmaz,
                                                  help
                                                  ));
                
                index_hitsarray++;
              }
          }
          mean_z = mean_z / static_cast< double >(mean_n);
          mean_x = mean_x / static_cast< double >(mean_n);
          mean_y = mean_y / static_cast< double >(mean_n);
	  
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
          TMinuit *gMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 4 params
          
          //  set print level (-1 = quiet, 0 = normal, 1 = verbose)
          gMinuit->SetPrintLevel(-1);
  
          gMinuit->SetFCN(fcn_wrapper2);
          
          double arglist[10];
          int ierflg = 0;

   
          // minimization strategy (1 = standard, 2 = slower)
          arglist[0] = 2;
          gMinuit->mnexcm("SET STR",arglist,2,ierflg);
          
          // set error definition (1 = for chi square)
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
          double ps = atan(linfit_x_a1/sqrt(1.0+linfit_y_a1*linfit_y_a1));//guess
                                                                          //of psi
            
          //  Set starting values and step sizes for parameters
          Double_t vstart[4] = {linfit_x_a0, linfit_y_a0, del, ps};
          //duble vstart[4] = {0.0, 0.0, 0.0, 0.0};
          double step[4] = {0.01, 0.01, 0.01, 0.01};
            
          gMinuit->mnparm(0, "b0", vstart[0], step[0], 0,0,ierflg);
          gMinuit->mnparm(1, "b1", vstart[1], step[1], 0,0,ierflg);
          gMinuit->mnparm(2, "delta", vstart[2], step[2], -1.0*TMath::Pi(), 1.0*TMath::Pi(),ierflg);
          gMinuit->mnparm(3, "psi", vstart[3], step[3], -1.0*TMath::Pi(), 1.0*TMath::Pi(),ierflg);
            

          //  Now ready for minimization step
          arglist[0] = 2000;
          arglist[1] = 0.01;
          gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
           
          bool ok = true;
            
          if(ierflg != 0)
          {
            ok = false;            
          }
                    
          //   get results from migrad
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
              
              for (unsigned int help =0; help < _nPlanes; help++)
                {
                  const double x = _xPos[track][help];
                  const double y = _yPos[track][help];
                  const double z = _zPos[track][help];
                 
                  //calculate the lambda parameter
                  const double la = -1.0*b0*c0-b1*c1+c0*x+c1*y+sqrt(1-c0*c0-c1*c1)*z;
                  lambda.push_back(la);
/*
		  if (_referenceHitVec == 0){                  
                  //determine the residuals without reference vector
                  _waferResidX[help] = b0 + la*c0 - x;
                  _waferResidY[help] = b1 + la*c1 - y;
                  _waferResidZ[help] = la*sqrt(1.0 - c0*c0 - c1*c1) - z;

		  } else {
*/
		    // use reference vector
		    TVector3 vpoint(b0,b1,0.);
		    TVector3 vvector(c0,c1,c2);
		    TVector3 point=Line2Plane(help, vpoint,vvector);

		    //determine the residuals
		    if( !(abs(x)<1e-06 && abs(y)<1e-06)  )
		      {
			_waferResidX[help] = point[0] - x;
			_waferResidY[help] = point[1] - y;
			_waferResidZ[help] = point[2] - z;
			// TODO: THIS LOOKS LIKE A BUG!!
			//cout << "point[2] = " << point[2] << ", z = " << z << endl;
			//cout << "z residual w/ ref vec: " << _waferResidZ[help] << ", w/o ref vec: " << la*sqrt(1.0 - c0*c0 - c1*c1) - z << endl;
		      }else{
		      _waferResidX[help] = 0.;
		      _waferResidY[help] = 0.;
		      _waferResidZ[help] = 0.;
		    }
		/*  
		}
		*/
                }
            }
          delete gMinuit;
        }
      }
 
      streamlog_out ( MESSAGE1 ) << "Residuals X: ";

      for (unsigned int help = 0; help < _nPlanes; help++) {
        streamlog_out ( MESSAGE1 ) << _waferResidX[help] << " ";
      }

      streamlog_out ( MESSAGE1 ) << endl;

      streamlog_out ( MESSAGE1 ) << "Residuals Y: ";

      for (unsigned int help = 0; help < _nPlanes; help++) {
        streamlog_out ( MESSAGE1 ) << _waferResidY[help] << " ";
      }

      streamlog_out ( MESSAGE1 ) << endl;

      streamlog_out ( MESSAGE1 ) << "Residuals Z: ";

      for (unsigned int help = 0; help < _nPlanes; help++) {
        streamlog_out ( MESSAGE1 ) << _waferResidZ[help] << " ";
      }

      streamlog_out ( MESSAGE1 ) << endl;


      int residualsXOkay = 1;
      int residualsYOkay = 1;

      // check if residal cuts are used
      if (_useResidualCuts != 0) 
      {

        // loop over all sensors
        for (unsigned int help = 0; help < _nPlanes; help++) 
        {
          int excluded = 0; //0 not excluded, 1 excluded
          if (_nExcludePlanes > 0) 
          {
            for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
            {
              if (help == _excludePlanes[helphelp]) 
              {
                excluded = 1;
              }
            }
          }
          if(excluded == 0)
          {
              if (_waferResidX[help] < _residualsXMin[help] || _waferResidX[help] > _residualsXMax[help]) 
              {
                residualsXOkay = 0;
              }
              if (_waferResidY[help] < _residualsYMin[help] || _waferResidY[help] > _residualsYMax[help]) 
              {
                residualsYOkay = 0;
              }
          }

        } // end loop over all sensors

      } // end check if residual cuts are used

      if (_useResidualCuts != 0 && (residualsXOkay == 0 || residualsYOkay == 0)) {
        streamlog_out ( MESSAGE1 ) << "Track did not pass the residual cuts." << endl;
      }

      // apply track cuts (at the moment only residuals)
      if (_useResidualCuts == 0 || (residualsXOkay == 1 && residualsYOkay == 1)) {

        // Add track to Millepede
        // ---------------------------

	if (_alignMode == 3) {
          if(validminuittrack || _inputMode == 1 )
            {
              const int nLC = 4; // number of local parameters
              const int nGL = _nPlanes * 6; // number of global parameters

              float *derLC = new float[nLC]; // array of derivatives for local parameters
              float *derGL = new float[nGL]; // array of derivatives for global parameters

              int *label = new int[nGL]; // array of labels

              float residual;

              // create labels
              for (int help = 0; help < nGL; help++) 
              {
                label[help] = help + 1;
              }

              for (int help = 0; help < nGL; help++) 
              {
                derGL[help] = 0;
              }

              for (int help = 0; help < nLC; help++) 
              {
                derLC[help] = 0;
              }

              int nExcluded = 0;

              // loop over all planes
             
              for (unsigned int help = 0; help < _nPlanes; help++) 
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

                if (excluded == 0) {
                  int helphelp = help - nExcluded; // index of plane after

		  for (int i = 0; i < nGL; i++ ) 
                  {
		    derGL[i] = 0.000;
		  }

		  for (int i = 0; i < nLC; i++ ) 
                  {
		    derLC[i] = 0.000;
		  }

                  double x_sensor = 0.;
                  double y_sensor = 0.;
                  double z_sensor = 0.;

		  if (_referenceHitVec != 0)
{
		    for(int ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
		      {
			EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
			if( _sensorIDVec[help] == refhit->getSensorID() )
			  {
			    x_sensor =  refhit->getXOffset();
			    y_sensor =  refhit->getYOffset();
			    z_sensor =  refhit->getZOffset();
			  } 
		      }
}
else
{
	int sensorID = _sensorIDVec[help];
	x_sensor = geo::gGeometry().siPlaneXPosition(sensorID); 
	y_sensor = geo::gGeometry().siPlaneYPosition(sensorID); 
	z_sensor = geo::gGeometry().siPlaneZPosition(sensorID); 
//std::cout << "Retrived: " << x_sensor << ", " << y_sensor << ", " << z_sensor << std::endl;
}
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
//
                  // shift in X
                  derGL[((helphelp * 6) + 0)] = -1.0;                                 // dx
                  derGL[((helphelp * 6) + 1)] =  0.0;                                 // dy
                  derGL[((helphelp * 6) + 2)] =  0.0;                                 // dz
                  // rotation in ZY ( alignment->getAlpfa() )
                  derGL[((helphelp * 6) + 3)] =      0.0;                             // alfa  - ZY :: Y->Z
                  derGL[((helphelp * 6) + 4)] = -1.0*(_zPosHere[help] - z_sensor);    // beta  - ZX :: Z->X
                  derGL[((helphelp * 6) + 5)] =  1.0*(_yPosHere[help] - y_sensor);                 // gamma - XY :: X->Y

                  derLC[0] = 1.0;
                  derLC[1] = 0.0;
                  derLC[2] = _zPosHere[help] + _waferResidZ[help];
                  derLC[3] = 0.0;
            
                  residual = _waferResidX[help];
                 
                  _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmax);

             
                  // shift in Y
                  derGL[((helphelp * 6) + 0)] =  0.0;
                  derGL[((helphelp * 6) + 1)] = -1.0; 
                  derGL[((helphelp * 6) + 2)] =  0.0;
                  // rotation in ZX
                  derGL[((helphelp * 6) + 3)] =  1.0*(_zPosHere[help] - z_sensor); 
                  derGL[((helphelp * 6) + 4)] =      0.0            ; 
                  derGL[((helphelp * 6) + 5)] = -1.0*(_xPosHere[help] - x_sensor);

                  derLC[0] = 0.0;
                  derLC[1] = 1.0;
                  derLC[2] = 0.0;
                  derLC[3] = _zPosHere[help] + _waferResidZ[help];
            
                  residual = _waferResidY[help];
                  
                  _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmay);
              
              
                  // shift in Z
                  derGL[((helphelp * 6) + 0)] =  0.0;
                  derGL[((helphelp * 6) + 1)] =  0.0;
                  derGL[((helphelp * 6) + 2)] = -1.0; 
                  // rotation in XY
                  derGL[((helphelp * 6) + 3)] = -1.0*(_yPosHere[help]-y_sensor); 
                  derGL[((helphelp * 6) + 4)] =  1.0*(_xPosHere[help]-x_sensor);
                  derGL[((helphelp * 6) + 5)] =  0.0;

                  derLC[0] = 0.0;
                  derLC[1] = 0.0;
                  derLC[2] = _xPosHere[help] + _waferResidX[help];
                  derLC[3] = _yPosHere[help] + _waferResidY[help];
            
                  residual = _waferResidZ[help];
                 
                  _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmaz);
                  _nMilleDataPoints++;

                } // end if plane is not excluded

              } // end loop over all planes

              // clean up


              delete [] derLC;
              delete [] derGL;
              delete [] label;

            }
        } else {

          streamlog_out ( ERROR2 ) << _alignMode << " is not a valid mode. Please choose 1,2 or 3." << endl;

        }

        _nGoodTracks++;

        // end local fit
        _mille->end();

        _nMilleTracks++;

        // Fill histograms for individual tracks
        // -------------------------------------

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

        string tempHistoName;

        if ( _histogramSwitch ) {
          if ( AIDA::IHistogram1D* chi2x_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_chi2XLocalname]) )
            chi2x_histo->fill(Chiquare[0]);
          else {
            streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _chi2XLocalname << endl;
            streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
            _histogramSwitch = false;
          }
        }

        if ( _histogramSwitch ) {
          if ( AIDA::IHistogram1D* chi2y_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_chi2YLocalname]) )
            chi2y_histo->fill(Chiquare[1]);
          else {
            streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _chi2YLocalname << endl;
            streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
            _histogramSwitch = false;
          }
        }


        // loop over all detector planes
        for(unsigned int iDetector = 0; iDetector < _nPlanes; iDetector++ ) {

          int sensorID = _orderedSensorID.at( iDetector );

          if ( 
              abs(_waferResidX[iDetector]) < 1e-06 &&  
              abs(_waferResidY[iDetector]) < 1e-06 &&  
              abs(_waferResidZ[iDetector]) < 1e-06 
             )   continue;


          if ( _histogramSwitch ) {
            tempHistoName = _residualXLocalname + "_d" + to_string( sensorID );
            if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
            {
                residx_histo->fill(_waferResidX[iDetector]);

                tempHistoName = _residualXvsYLocalname + "_d" + to_string( sensorID );
                AIDA::IProfile1D* residxvsY_histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMapProf1D[tempHistoName.c_str()]) ;
                residxvsY_histo->fill(_yPosHere[iDetector], _waferResidX[iDetector]);
 
                tempHistoName = _residualXvsXLocalname + "_d" + to_string( sensorID );
                AIDA::IProfile1D* residxvsX_histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMapProf1D[tempHistoName.c_str()]) ;
                residxvsX_histo->fill(_xPosHere[iDetector], _waferResidX[iDetector]);
            }
            else
            {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }

          if ( _histogramSwitch ) {
            tempHistoName = _residualYLocalname + "_d" + to_string( sensorID );
            if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
            {
              residy_histo->fill(_waferResidY[iDetector]);

              tempHistoName = _residualYvsYLocalname + "_d" + to_string( sensorID );
              AIDA::IProfile1D* residyvsY_histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMapProf1D[tempHistoName.c_str()]) ;
              residyvsY_histo->fill(_yPosHere[iDetector], _waferResidY[iDetector]);
 
              tempHistoName = _residualYvsXLocalname + "_d" + to_string( sensorID );
              AIDA::IProfile1D* residyvsX_histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMapProf1D[tempHistoName.c_str()]) ;
              residyvsX_histo->fill(_xPosHere[iDetector], _waferResidY[iDetector]);
            }
            else
            {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualYLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }
          if ( _histogramSwitch ) {
            tempHistoName = _residualZLocalname + "_d" + to_string( sensorID );
            if ( AIDA::IHistogram1D* residz_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
            {
              residz_histo->fill(_waferResidZ[iDetector]);

              tempHistoName = _residualZvsYLocalname + "_d" + to_string( sensorID );
              AIDA::IProfile1D* residzvsY_histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMapProf1D[tempHistoName.c_str()]) ;
              residzvsY_histo->fill(_yPosHere[iDetector], _waferResidZ[iDetector]);
 
              tempHistoName = _residualZvsXLocalname + "_d" + to_string( sensorID );
              AIDA::IProfile1D* residzvsX_histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMapProf1D[tempHistoName.c_str()]) ;
              residzvsX_histo->fill(_xPosHere[iDetector], _waferResidZ[iDetector]);
             }
            else 
            {
              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualZLocalname << endl;
              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
              _histogramSwitch = false;
            }
          }
        } // end loop over all detector planes


#endif

      } // end if apply track cuts

      // clean up
      delete [] _zPosHere;
      delete [] _yPosHere;
      delete [] _xPosHere;

    } // end loop over all track candidates

  } // end if only one track or no single track event

  streamlog_out ( MESSAGE1 ) << "Finished fitting tracks in event " << _iEvt << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  string tempHistoName;

  if ( _histogramSwitch ) {
    {
      stringstream ss;
      ss << _numberTracksLocalname << endl;
    }
    if ( AIDA::IHistogram1D* number_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[_numberTracksLocalname]) )
      number_histo->fill(_nGoodTracks);
    else {
      streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _numberTracksLocalname << endl;
      streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
      _histogramSwitch = false;
    }
  }

#endif

  // count events
  ++_iEvt;
  if ( isFirstEvent() ) _isFirstEvent = false;

}

TVector3 EUTelMille2::Line2Plane(int iplane, const TVector3& lpoint, const TVector3& lvector ) 
{
  TVector3 hitInPlane;
  TVector3 norm2Plane;
 
  if( _referenceHitVec == 0 )
  {
	int sensorID = _orderedSensorID[iplane];
  	hitInPlane.SetXYZ( geo::gGeometry().siPlaneXPosition(sensorID)*1000, geo::gGeometry().siPlaneYPosition(sensorID)*1000, geo::gGeometry().siPlaneZPosition(sensorID)*1000 );
	norm2Plane = geo::gGeometry().siPlaneNormal(sensorID);

//std::cout << "Retrived (offset): " << sensorID << ":" << hitInPlane(0) << ", " << hitInPlane(1) << ", " << hitInPlane(2) << std::endl;
//std::cout << "Retrived (angle): "<< sensorID  << ":" << norm2Plane(0) << ", " << norm2Plane(1) << ", " << norm2Plane(2) << std::endl;
  }
  else
  {
  	EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(iplane) ) ;
  	hitInPlane.SetXYZ( refhit->getXOffset()*1000, refhit->getYOffset()*1000, refhit->getZOffset()*1000);
	norm2Plane.SetXYZ( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
	} 

        TVector3 point( 1.,1.,1. );
          
        double linecoord_numenator   = norm2Plane.Dot(hitInPlane-lpoint);
        double linecoord_denumenator = norm2Plane.Dot(lvector);

        point = (linecoord_numenator/linecoord_denumenator)*lvector + lpoint;

  return point;
}

void EUTelMille2::end() {

  delete [] _telescopeResolY;
  delete [] _telescopeResolX;
  delete [] _telescopeResolZ;
  delete [] _yFitPos;
  delete [] _xFitPos;
  delete [] _waferResidY;
  delete [] _waferResidX;
  delete [] _waferResidZ;

  if(_alignMode == 3)
    {
      delete []  hitsarray2;
    }

  // close the output file
  delete _mille;

  // if write the pede steering file
  if (_generatePedeSteerfile) {

    streamlog_out ( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

    string tempHistoName;
    double *meanX = new double[_nPlanes];
    double *meanY = new double[_nPlanes];
    double *meanZ = new double[_nPlanes];

    // loop over all detector planes
    for(unsigned int iDetector = 0; iDetector < _nPlanes; iDetector++ ) {

      int sensorID = _orderedSensorID.at( iDetector );

      if ( _histogramSwitch ) {
        tempHistoName =  _residualXLocalname + "_d" + to_string( sensorID );
        if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          meanX[iDetector] = residx_histo->mean();
        else {
          streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
          streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
          _histogramSwitch = false;
        }
      }

      if ( _histogramSwitch ) {
        tempHistoName =  _residualYLocalname + "_d" + to_string( sensorID );
        if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          meanY[iDetector] = residy_histo->mean();
        else {
          streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualYLocalname << endl;
          streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
          _histogramSwitch = false;
        }
      }

      if ( _histogramSwitch ) {
        tempHistoName =  _residualZLocalname + "_d" + to_string( sensorID );
        if ( AIDA::IHistogram1D* residz_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
          meanZ[iDetector] = residz_histo->mean();
        else {
          streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualZLocalname << endl;
          streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
          _histogramSwitch = false;
        } 
      }
    } // end loop over all detector planes

    ofstream steerFile;
    steerFile.open(_pedeSteerfileName.c_str());

    if (steerFile.is_open()) {

      // find first and last excluded plane
      unsigned int firstnotexcl = _nPlanes;
      unsigned int lastnotexcl = 0;

      // loop over all planes
      for (unsigned int help = 0; help < _nPlanes; help++) 
      {

        int excluded = 0;

        // loop over all excluded planes
        for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
        {
          if (help == _excludePlanes[helphelp]) 
          {
//            excluded = 1;
          }
        } // end loop over all excluded planes

        if (excluded == 0 && firstnotexcl > help) 
        {
          firstnotexcl = help;
        }

        if (excluded == 0 && lastnotexcl < help) 
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
      for (unsigned int help = 0; help < _nPlanes; help++) {
 
        int excluded = 0; // flag for excluded planes

        // loop over all excluded planes
        for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {

          if (help == _excludePlanes[helphelp]) {
//            excluded = 1;
          }

        } // end loop over all excluded planes

        // if plane not excluded
        if (excluded == 0) {
          
          bool fixed = false;
          for(size_t i = 0;i< _FixedPlanes.size(); i++)
            {
              if(_FixedPlanes[i] == static_cast< int >(help))
                fixed = true;
            }
          
          if( fixed || (_FixedPlanes.empty() && (help == firstnotexcl || help == lastnotexcl) ) )
            {
		 if (_alignMode == 3) {
                steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
                steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
                steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
                steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
                steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
                steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
              }
              
            } else {
            
		 if (_alignMode == 3) {
              if (_usePedeUserStartValues == 0)
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
      for(StringVec::iterator it = _pedeSteerAddCmds.begin(); it != _pedeSteerAddCmds.end(); ++it) {
	// two backslashes will be interpreted as newline
	if (*it == "\\\\")
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

      streamlog_out ( MESSAGE5 ) << "File " << _pedeSteerfileName << " written." << endl;

    } else {

      streamlog_out ( ERROR2 ) << "Could not open steering file." << endl;

    }
    //cleaning up
    delete [] meanX;
    delete [] meanY;
    delete [] meanZ;
  } // end if write the pede steering file

  streamlog_out ( MESSAGE7 ) << "Number of data points used: " << _nMilleDataPoints << endl;
  streamlog_out ( MESSAGE7 ) << "Number of tracks used: " << _nMilleTracks << endl;

  // monitor the number of tracks in CDash when running tests
  CDashMeasurement meas_ntracks("ntracks",_nMilleTracks);  cout << meas_ntracks;  // output only if DO_TESTING is set

  // if running pede using the generated steering file
  if (_runPede == 1) {

    // check if steering file exists
    if (_generatePedeSteerfile == 1) {

      std::string command = "pede " + _pedeSteerfileName;

      streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str() << endl;

      bool encounteredError = false;
      
      // run pede and create a streambuf that reads its stdout and stderr
      redi::ipstream pede( command.c_str(), redi::pstreams::pstdout|redi::pstreams::pstderr ); 
      
      if (!pede.is_open()) {
	  streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl;
	  encounteredError = true;	  
      } else {

	// output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
	char buf[1024];
	std::streamsize n;
	std::stringstream pedeoutput; // store stdout to parse later
	std::stringstream pedeerrors;
	bool finished[2] = { false, false };
	while (!finished[0] || !finished[1])
	  {
	    if (!finished[0])
	      {
		while ((n = pede.err().readsome(buf, sizeof(buf))) > 0){
		  streamlog_out( ERROR5 ).write(buf, n).flush();
		  string error (buf, n);
		  pedeerrors << error;
		  encounteredError = true;
		}
		if (pede.eof())
		  {
		    finished[0] = true;
		    if (!finished[1])
		      pede.clear();
		  }
	      }

	    if (!finished[1])
	      {
		while ((n = pede.out().readsome(buf, sizeof(buf))) > 0){
		  streamlog_out( MESSAGE4 ).write(buf, n).flush();
		  string output (buf, n);
		  pedeoutput << output;
		}
		if (pede.eof())
		  {
		    finished[1] = true;
		    if (!finished[0])
		      pede.clear();
		  }
	      }
	  }

	// pede does not return exit codes on some errors (in V03-04-00)
	// check for some of those here by parsing the output
	{
	  const char * pch = strstr(pedeoutput.str().data(),"Too many rejects");
	  if (pch){
	    streamlog_out ( ERROR5 ) << "Pede stopped due to the large number of rejects. " << endl;
	    encounteredError = true;
	  }
	}
	
	{
	  const char* pch0 = strstr(pedeoutput.str().data(),"Sum(Chi^2)/Sum(Ndf) = ");
	  if (pch0 != 0){
	    streamlog_out ( DEBUG5 ) << " Parsing pede output for final chi2/ndf result.. " << endl;
	    // search for the equal sign after which the result for chi2/ndf is stated within the next 80 chars 
	    // (with offset of 22 chars since pch points to beginning of "Sum(..." string just found)
	    char* pch = (char*)((memchr (pch0+22, '=', 180)));
	    if (pch!=NULL){
	      char str[16];
	      // now copy the numbers after the equal sign
	      strncpy ( str, pch+1, 15 );
	      str[15] = '\0';   /* null character manually added */
	      // monitor the chi2/ndf in CDash when running tests
	      CDashMeasurement meas_chi2ndf("chi2_ndf",atof(str));  cout << meas_chi2ndf; // output only if DO_TESTING is set
	      streamlog_out ( MESSAGE6 ) << "Final Sum(Chi^2)/Sum(Ndf) = " << str << endl;
	    }	    
	  }
	}

        // wait for the pede execution to finish
        pede.close();

        // check the exit value of pede / react to previous errors
        if ( pede.rdbuf()->status() == 0 && !encounteredError) 
        {
          streamlog_out ( MESSAGE7 ) << "Pede successfully finished" << endl;
        } else {
          streamlog_out ( ERROR5 ) << "Problem during Pede execution, exit status: " << pede.rdbuf()->status() << ", error messages (repeated here): " << endl;
	  streamlog_out ( ERROR5 ) << pedeerrors.str() << endl;
	  // TODO: decide what to do now; exit? and if, how?
          streamlog_out ( ERROR5 ) << "Will exit now" << endl;
	  //exit(EXIT_FAILURE); // FIXME: can lead to (ROOT?) seg faults - points to corrupt memory? run valgrind...
	  return; // does fine for now
	}

        // reading back the millepede.res file and getting the
        // results.
        string millepedeResFileName = "millepede.res";

        streamlog_out ( MESSAGE6 ) << "Reading back the " << millepedeResFileName << endl
                                   << "Saving the alignment constant into " << _alignmentConstantLCIOFile << endl;

        // open the millepede ASCII output file
        ifstream millepede( millepedeResFileName.c_str() );


        // reopen the LCIO file this time in append mode
        LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

        try 
        {
          lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
        }
        catch ( IOException& e ) 
        {
          streamlog_out ( ERROR4 ) << e.what() << endl
                                   << "Sorry for quitting. " << endl;
          exit(-1);
        }


        // write an almost empty run header
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


        if ( millepede.bad() || !millepede.is_open() ) 
        {
          streamlog_out ( ERROR4 ) << "Error opening the " << millepedeResFileName << endl
                                   << "The alignment slcio file cannot be saved" << endl;
        }
        else 
        {
          vector<double > tokens;
          stringstream tokenizer;
          string line;

          // get the first line and throw it away since it is a
          // comment!
          getline( millepede, line );

          int counter = 0;

          while ( ! millepede.eof() ) {

            EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;

            bool goodLine = true;
            unsigned int numpars = 0;
              numpars = 6;

            for ( unsigned int iParam = 0 ; iParam < numpars ; ++iParam ) 
            {
              getline( millepede, line );

              if ( line.empty() ) {
                goodLine = false;
                continue;
              }

              tokens.clear();
              tokenizer.clear();
              tokenizer.str( line );

	      double buffer;
              // // check that all parts of the line are non zero
              while ( tokenizer >> buffer ) {
                tokens.push_back( buffer ) ;
              }

              if ( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) || (tokens.size() == 5) ) {
                goodLine = true;
              } else goodLine = false;

              bool isFixed = ( tokens.size() == 3 );

                 if ( iParam == 0 ) {
                    constant->setXOffset( tokens[1] / 1000. );
                    if ( ! isFixed ) constant->setXOffsetError( tokens[4] / 1000. ) ;                    
                  }
                  if ( iParam == 1 ) {
                    constant->setYOffset( tokens[1] / 1000. ) ;
                    if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000. ) ;
                  }
                  if ( iParam == 2 ) {
                    constant->setZOffset( tokens[1] / 1000. ) ;
                    if ( ! isFixed ) constant->setZOffsetError( tokens[4] / 1000. ) ;
                  }
                  if ( iParam == 3 ) {
                    constant->setAlpha( tokens[1]  ) ;
                    if ( ! isFixed ) constant->setAlphaError( tokens[4] ) ;
                  } 
                  if ( iParam == 4 ) {
                    constant->setBeta( tokens[1]  ) ;
                    if ( ! isFixed ) constant->setBetaError( tokens[4] ) ;
                  } 
                  if ( iParam == 5 ) {
                    constant->setGamma( tokens[1]  ) ;
                    if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
                  } 
            }

            // right place to add the constant to the collection
            if ( goodLine  ) {
//               constant->setSensorID( _orderedSensorID_wo_excluded.at( counter ) );
              constant->setSensorID( _orderedSensorID.at( counter ) );
              ++ counter;
              constantsCollection->push_back( constant );
              streamlog_out ( MESSAGE0 ) << (*constant) << endl;
            }
            else delete constant;
          }

        }



        event->addCollection( constantsCollection, _alignmentConstantCollectionName );
        lcWriter->writeEvent( event );
        delete event;

        lcWriter->close();

        millepede.close();

      }
    } else {

      streamlog_out ( ERROR2 ) << "Unable to run pede. No steering file has been generated." << endl;

    }


  } // end if running pede using the generated steering file

  streamlog_out ( MESSAGE2 ) << endl;
  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;
}

void EUTelMille2::bookHistos() {


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out ( MESSAGE2 ) << "Booking histograms..." << endl;

    const int    tracksNBin = 20  ;
    const double tracksMin  = -0.5;
    const double tracksMax  = 19.5;
    const int    Chi2NBin = 1000 ;
    const double Chi2Min  =   0.;
    const double Chi2Max  = 1000.;
    const int    NBin =   4000;
    const double Min  = -2000.;
    const double Max  =  2000.;

    AIDA::IHistogram1D * numberTracksLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_numberTracksLocalname,tracksNBin,tracksMin,tracksMax);
    if ( numberTracksLocal ) {
      numberTracksLocal->setTitle("Number of tracks after #chi^{2} cut");
      _aidaHistoMap.insert( make_pair( _numberTracksLocalname, numberTracksLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_numberTracksLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2XLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2XLocalname,Chi2NBin,Chi2Min,Chi2Max);
    if ( chi2XLocal ) {
      chi2XLocal->setTitle("Chi2 X");
      _aidaHistoMap.insert( make_pair( _chi2XLocalname, chi2XLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_chi2XLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    AIDA::IHistogram1D * chi2YLocal =
      AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2YLocalname,Chi2NBin,Chi2Min,Chi2Max);
    if ( chi2YLocal ) {
      chi2YLocal->setTitle("Chi2 Y");
      _aidaHistoMap.insert( make_pair( _chi2YLocalname, chi2YLocal ) );
    } else {
      streamlog_out ( ERROR2 ) << "Problem booking the " << (_chi2YLocalname) << endl;
      streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }

    string tempHistoName;
    string histoTitleXResid;
    string histoTitleYResid;
    string histoTitleZResid;

    for(unsigned int iDetector = 0; iDetector < _nPlanes; iDetector++ ){

      // this is the sensorID corresponding to this plane
      int sensorID = _orderedSensorID.at( iDetector );

      tempHistoName     =  _residualXLocalname + "_d" + to_string( sensorID );
      histoTitleXResid  =  "XResidual_d" + to_string( sensorID ) ;

      AIDA::IHistogram1D *  tempXHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempXHisto ) {
        tempXHisto->setTitle(histoTitleXResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempXHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      tempHistoName     =  _residualXvsXLocalname + "_d" + to_string( sensorID );
      histoTitleXResid  =  "XvsXResidual_d" + to_string( sensorID ) ;

      AIDA::IProfile1D *  tempX2dHisto =
        AIDAProcessor::histogramFactory(this)->createProfile1D(tempHistoName, 100, -10000., 10000.,  -10000., 10000. );
      if ( tempX2dHisto ) {
        tempX2dHisto->setTitle(histoTitleXResid);
        _aidaHistoMapProf1D.insert( make_pair( tempHistoName, tempX2dHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      tempHistoName     =  _residualXvsYLocalname + "_d" + to_string( sensorID );
      histoTitleXResid  =  "XvsYResidual_d" + to_string( sensorID ) ;

       tempX2dHisto =
        AIDAProcessor::histogramFactory(this)->createProfile1D(tempHistoName,  100, -10000., 10000.,  -10000., 10000.);
      if ( tempX2dHisto ) {
        tempX2dHisto->setTitle(histoTitleXResid);
        _aidaHistoMapProf1D.insert( make_pair( tempHistoName, tempX2dHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }


      tempHistoName     =  _residualYLocalname + "_d" + to_string( sensorID );
      histoTitleYResid  =  "YResidual_d" + to_string( sensorID ) ;

      AIDA::IHistogram1D *  tempYHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempYHisto ) {
        tempYHisto->setTitle(histoTitleYResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempYHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      tempHistoName     =  _residualYvsXLocalname + "_d" + to_string( sensorID );
      histoTitleYResid  =  "YvsXResidual_d" + to_string( sensorID ) ;

      AIDA::IProfile1D *  tempY2dHisto =
        AIDAProcessor::histogramFactory(this)->createProfile1D(tempHistoName, 100, -10000., 10000.,  -1000., 1000.);
      if ( tempY2dHisto ) {
        tempY2dHisto->setTitle(histoTitleYResid);
        _aidaHistoMapProf1D.insert( make_pair( tempHistoName, tempY2dHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      tempHistoName     =  _residualYvsYLocalname + "_d" + to_string( sensorID );
      histoTitleYResid  =  "YvsYResidual_d" + to_string( sensorID ) ;

        tempY2dHisto =
        AIDAProcessor::histogramFactory(this)->createProfile1D(tempHistoName, 100, -10000., 10000.,  -10000., 10000.);
      if ( tempY2dHisto ) {
        tempY2dHisto->setTitle(histoTitleYResid);
        _aidaHistoMapProf1D.insert( make_pair( tempHistoName, tempY2dHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      tempHistoName     =  _residualZLocalname + "_d" + to_string( sensorID );
      histoTitleZResid  =  "ZResidual_d" + to_string( sensorID ) ;

      AIDA::IHistogram1D *  tempZHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName,NBin, Min,Max);
      if ( tempZHisto ) {
        tempZHisto->setTitle(histoTitleZResid);
        _aidaHistoMap.insert( make_pair( tempHistoName, tempZHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }
 
 
      tempHistoName     =  _residualZvsYLocalname + "_d" + to_string( sensorID );
      histoTitleZResid  =  "ZvsYResidual_d" + to_string( sensorID ) ;

      AIDA::IProfile1D*  tempZ2dHisto =
        AIDAProcessor::histogramFactory(this)->createProfile1D(tempHistoName, 100, -10000., 10000.,  -10000., 10000.);
      if ( tempZ2dHisto ) {
        tempZ2dHisto->setTitle(histoTitleZResid);
        _aidaHistoMapProf1D.insert( make_pair( tempHistoName, tempZ2dHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }


      tempHistoName     =  _residualZvsXLocalname + "_d" + to_string( sensorID );
      histoTitleZResid  =  "ZvsXResidual_d" + to_string( sensorID ) ;

        tempZ2dHisto =
        AIDAProcessor::histogramFactory(this)->createProfile1D(tempHistoName, 100, -10000., 10000.,  -10000., 10000.);
      if ( tempZ2dHisto ) {
        tempZ2dHisto->setTitle(histoTitleZResid);
        _aidaHistoMapProf1D.insert( make_pair( tempHistoName, tempZ2dHisto ) );
      } else {
        streamlog_out ( ERROR2 ) << "Problem booking the " << (tempHistoName) << endl;
        streamlog_out ( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }
    }
  } catch (lcio::Exception& e ) {


#ifdef EUTEL_INTERACTIVE
    streamlog_out ( ERROR5 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
    string answer;
    while ( true ) {
      streamlog_out ( ERROR5 ) << "[q]/[c]" << endl;
      cin >> answer;
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" )
        _histogramSwitch = false;
      break;
    }
#else
    streamlog_out( WARNING6 ) << "No AIDAProcessor initialized. Continue without histogramming" << endl;

#endif


  }
#endif

}

#endif // USE_GEAR


