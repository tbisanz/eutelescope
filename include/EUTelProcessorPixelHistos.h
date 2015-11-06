/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPROCESSORPIXELHISTOS_H
#define EUTELPROCESSORPIXELHISTOS_H

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelGenericSparsePixel.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>

// AIDA includes <.h>
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>

// system includes <>
#include <map>
#include <vector>

namespace eutelescope {

//! Processor to write out hot pixels 
/*! This processor is used to keep hot matrix out from the analysis
 *  procedure. It checks if pixels fired above a certain frequency
 *  and writes them into a TrackerData collection in an external 
 *  file if they did.
 *  You have to specify which plane IDs shall be processed (SensorIDVec)
 *  as well as which ones shall be excluded (ExcludedPlanes). This is
 *  necessary for histogramming and fast data processing.
 *
 *  @param NoOfEvents The amount of events to determine the firing frequency
 *
 *  @param SensorIDVec An integer vector containing the sensor IDs of the
 *  planes which are processed
 *
 *  @param MaxAllowedFiringFreq The firing frequency cut which shall be applied
 *
 *  @param HotPixelDBFile Name of the output file, currently appending to a file 
 *  does not work, use differnt files if multiple hot pixel collections are created
 *
 *  @param ExcludedPlanes Planes to be excluded from processing
 *
 *  @param HotPixelCollectionName The name of the collection in the output file
 */
class EUTelProcessorPixelHistos : public marlin::Processor {

public:
    //! Returns a new instance of EUTelProcessorPixelHistos
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelProcessorPixelHistos.
     */
    virtual Processor* newProcessor() {
      return new EUTelProcessorPixelHistos;
    }

    //! Default constructor
    EUTelProcessorPixelHistos();

    //! Default destructor
    ~EUTelProcessorPixelHistos() = default;

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and performs some asserts about
     *  the value of the provided parameters
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run LCRunHeader of the this current run
     *
     *  @throw InvalidParameterException if a paramter is wrongly set
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. If the current @c
     *  evt is flagged to be used for update, then the selected
     *  algorithm wrapper is called
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     *
     *  @throw InvalidParameterException if information in the cellID
     *  are inconsistence
     */
    virtual void processEvent(LCEvent * evt);

    //! Check call back
    /*! This method is called every event just after the processEvent
     *  one. For the time being it is just calling the pixel
     *  monitoring protected method
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    //virtual void check( LCEvent* event );

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished. Just printing a good bye message
     */
    virtual void end();


protected:
    //! book histogram method
    void bookHistos();
    std::vector<std::vector<AIDA::IHistogram1D*>> _arrSignal;
    std::vector<std::vector<AIDA::IHistogram1D*>> _arrTime;
    std::pair<int, int> mapPixels(int, int);

    int _maxX;
    int _minX;
    int _maxY;
    int _minY;
    
    //! Input collection name for ZS data
    /*! The input collection is the calibrated data one coming from
     *  the input data file. It is, usually, called
     *  "zsdata" and it is a collection of TrackerData
     */
    std::string _zsDataCollectionName;
    //! Number of events for update cycle
    int _noOfEvents;
   	
    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Sensor ID vector
    /*! Passed as a argument via the steering file, here you
     *  specify for which sensors hot pixels should be determined
     */
    int  _chosenSensorID;
};

//! A global instance of the processor
EUTelProcessorPixelHistos gEUTelProcessorPixelHistos;

}
#endif
