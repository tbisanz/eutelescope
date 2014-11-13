// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMULTILINEFIT2_H
#define EUTELMULTILINEFIT2_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelUtility.h"

//#include "TrackerHitImpl2.h"
#include "IMPL/TrackerHitImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "mille/Mille.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TMinuit.h>
#include <TSystem.h>
#include <TMath.h>
#include <TVector3.h>
class TMinuit;
#else
#error *** You need ROOT to compile this code.  *** 
#endif


namespace eutelescope {


  //! Straight line fit processor
  /*!
   *
   */

  class EUTelMille2 : public marlin::Processor {

  public:
    class hit2
    {
    public:
      hit2(){}

      hit2(double tx, double ty, double tz, double rx, double ry, double rz,int i)
      {
        x = tx;
        y = ty;
        z = tz;
        resolution_x = rx;
        resolution_y = ry;
        resolution_z = rz;
    
        planenumber = i;
      }
      double x;
      double y;
      double z;
      double resolution_x;
      double resolution_y;
      double resolution_z;
  
      int planenumber;
    };

   

    class trackfitter2
    {
    public:
      ~trackfitter2(){}
      trackfitter2(){}
      trackfitter2(hit2 *h, unsigned int num)
      {
        hitsarray2 = h;
        n = num;
      }

      double dot(const double *a, const double *b) const
      {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
      }
      double fit (double *x)
      {
        double chi2 = 0.0;
   
        const unsigned int n = 3;

        const double b0 = x[0];
        const double b1 = x[1];
        const double b2 = 0.0;
    
        const double alpha = x[2];
        const double beta =  x[3];

        const double c0 = TMath::Sin(beta);
        const double c1 = -1.0*TMath::Cos(beta) * TMath::Sin(alpha);
        const double c2 = TMath::Cos(alpha) * TMath::Cos(beta);
    
        double c[n] = {c0, c1, c2}; 
        for(size_t i = 0; i < n;i++)
          {
            const double p0 = hitsarray2[i].x;
            const double p1 = hitsarray2[i].y;
            const double p2 = hitsarray2[i].z;
        
            const double resol_x = hitsarray2[i].resolution_x;
            const double resol_y = hitsarray2[i].resolution_y;
            const double resol_z = hitsarray2[i].resolution_z;
        
            const double pmb[n] = {p0-b0, p1-b1, p2-b2}; //p - b
        
            const double coeff = dot(c, pmb);
            const double t[n] = {
              b0 + c0 * coeff - p0,
              b1 + c1 * coeff - p1,
              b2 + c2 * coeff - p2
            }; 
        
            //sum of distances divided by resolution^2
            chi2 += t[0]*t[0] / pow(resol_x,2) 
              + t[1]*t[1] / pow(resol_y,2) 
              + t[2]*t[2] / pow(resol_z,2);
          }
 
        return chi2;
      }
    private:
    hit2 *hitsarray2;
    unsigned int n;
    };

    //! Variables for hit parameters
    class HitsInPlane2 {
    public:
      HitsInPlane2(){
        measuredX = 0.0;
        measuredY = 0.0;
        measuredZ = 0.0;
      }
      HitsInPlane2(double x, double y, double z)
      {
        measuredX = x;
        measuredY = y;
        measuredZ = z;
      }
      bool operator<(const HitsInPlane2& b) const
      {
        return (measuredZ < b.measuredZ);
      }
      double measuredX;
      double measuredY;
      double measuredZ;
    };

    //recursive method which searches for track candidates - with omits!
    virtual void findtracks2(
                            int missinghits,
                            std::vector<IntVec > &indexarray, //resulting vector of hit indizes
                            IntVec vec, //for internal use
                            std::vector<std::vector<EUTelMille2::HitsInPlane2> > &_hitsArray, //contains all hits for each plane
                            unsigned int i, //plane number
                            int y //hit index number
                            );

    //! Returns a new instance of EUTelMille2
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMille2.
     */
    virtual Processor * newProcessor() {
      return new EUTelMille2;
    }

    //! Default constructor
    EUTelMille2 ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description.
     *
     *  The cluster center might be calculate using a standard linear
     *  charge center of gravity algortihm or applying a more
     *  sophisticated non linear eta function. This behaviour is
     *  regulated by the user from the steering file.
     *
     *  @throw UnknownDataTypeException if the cluster type is unknown
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);

    //! Track search method, fills class members _xPos and _trackResidX (Y,Z as well) 
    /** */
    void findMatchedHits( int &, Track* ) ;

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();

    virtual inline int getAllowedMissingHits(){return _allowedMissingHits;}
    virtual inline int getMimosa26ClusterChargeMin(){return _mimosa26ClusterChargeMin;}


  protected:


    //! Ordered sensor ID
    /*! Within the processor all the loops are done up to _nPlanes and
     *  according to their position along the Z axis (beam axis).
     *
     *  This vector is containing the sensorID sorted according to the
     *  same rule.
     */
    IntVec _orderedSensorID;
    IntVec _orderedSensorID_wo_excluded;


    //! reference HitCollection name 
    /*!
     */
    std::string      _referenceHitCollectionName;
    bool             _useReferenceHitCollection;
    LCCollectionVec* _referenceHitVec;    
 
    //! TrackerHit collection name
    /*! Input collection with hits.
     */
    StringVec _hitCollectionName;

    //! TRACK collection name
    /*! Output collection with fitted tracks.
     */
    std::string _trackCollectionName;

    //! Hot pixel collection name.
    /*! 
     * this collection is saved in a db file to be used at the clustering level
     */
    std::string _hotPixelCollectionName;

    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  first level key   sensor unique 
     *              value sensor map
     *  sensor map key    unique row number
     *             value  vector of column numbers.
     */
    
    std::map<std::string, bool > _hotPixelMap;

    //! Sensor ID vector
    IntVec _sensorIDVec;

    //! Sensor ID map (inverse sensorIDVec) 
    std::map< int, int > _sensorIDVecMap;
    //! Sensor ID vector, 
    /*! it's position along Z axis
     */ 
    IntVec _sensorIDVecZOrder;
    //! sensor ID to position along Z id
    /*!
     */
    std::map<int, int> _sensorIDtoZOrderMap;


    // parameters

    float _distanceMax;
    FloatVec _distanceMaxVec;
    std::vector<unsigned int > _excludePlanes; //only for internal usage
    IntVec _excludePlanes_sensorIDs; //this is going to be
                                                //set by the user.
    IntVec _FixedPlanes; //only for internal usage
    IntVec _FixedPlanes_sensorIDs; //this is going to be
    //set by the user.
    
    StringVec _pedeSteerAddCmds; // allows user-added commands in the pede steering file


    int _maxTrackCandidates;
    int _maxTrackCandidatesTotal;

    std::string _binaryFilename;

    float _telescopeResolution;
    bool _onlySingleHitEvents;
    bool _onlySingleTrackEvents;
    int _alignMode;
    bool _useResidualCuts;

    FloatVec _residualsXMin;
    FloatVec _residualsYMin;
    FloatVec _residualsXMax;
    FloatVec _residualsYMax;

    FloatVec _resolutionX;
    FloatVec _resolutionY;
    FloatVec _resolutionZ;

    IntVec _FixParameter;


    int _generatePedeSteerfile;
    std::string _pedeSteerfileName;
    bool _runPede;
    int _usePedeUserStartValues;
    FloatVec _pedeUserStartValuesX;
    FloatVec _pedeUserStartValuesY;
    FloatVec _pedeUserStartValuesZ;
    
    FloatVec _pedeUserStartValuesAlpha;
    FloatVec _pedeUserStartValuesBeta;
    FloatVec _pedeUserStartValuesGamma;

    int _inputMode;
    int _allowedMissingHits;
    int _mimosa26ClusterChargeMin;

    float _testModeSensorResolution;
    float _testModeXTrackSlope;
    float _testModeYTrackSlope;

    FloatVec _testModeSensorZPositions;

    FloatVec _testModeSensorXShifts;
    FloatVec _testModeSensorYShifts;
    FloatVec _testModeSensorGamma;
    FloatVec _testModeSensorAlpha;
    FloatVec _testModeSensorBeta;

    std::string _alignmentConstantLCIOFile;
    std::string _alignmentConstantCollectionName;

  private:

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    // Excluded planes
    int _nExcludePlanes;

    // Statistics
    int _nMilleDataPoints;
    int _nMilleTracks;

    // Mille
    Mille * _mille;

    //! Conversion ID map.
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the detectorID to the
     *  layerindex.
     */
    std::map< int, int > _conversionIdMap;

    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;
    size_t _nPlanes;

    std::vector<DoubleVec > _xPos;
    std::vector<DoubleVec > _yPos;
    std::vector<DoubleVec > _zPos;

    std::vector<DoubleVec > _trackResidX;
    std::vector<DoubleVec > _trackResidY;
    std::vector<DoubleVec > _trackResidZ;


//    std::vector<DoubleVec > _xPosTrack;
//    std::vector<DoubleVec > _yPosTrack;
//    std::vector<DoubleVec > _zPosTrack;

    double * _xPosHere;
    double * _yPosHere;
    double * _zPosHere;
    double * _waferResidX;
    double * _waferResidY;
    double * _waferResidZ;
    double * _xFitPos;
    double * _yFitPos;

    DoubleVec _siPlaneZPosition;

    //! Fill histogram switch
    /*! Only for debug reason
     */
    bool _histogramSwitch;
  };

  //! A global instance of the processor
  EUTelMille2 gEUTelMille2;

}
#endif
#endif

