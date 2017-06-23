// Writen by Alexander Morton <alexander.morton@desy.de>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTelProcessorDUTGlobalEfficiency_h 
#define EUTelProcessorDUTGlobalEfficiency_h

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesLayerLayout.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <map>
#include <set>
#include <string>
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

namespace eutelescope {

  class EUTelProcessorDUTGlobalEfficiency : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelProcessorDUTGlobalEfficiency) 

  public:
    // Returns a new instance of EUTelProcessorDUTGlobalEfficiency
    virtual Processor *newProcessor() {
      return new EUTelProcessorDUTGlobalEfficiency;
    }

    EUTelProcessorDUTGlobalEfficiency(); // Default constructor
    // Called only at the begining of a job
    virtual void init();

    // Called every run
    virtual void processRunHeader(LCRunHeader *run);

    // Called every event.
    virtual void processEvent(LCEvent *event);

    // Called at the end of the job
    virtual void end();

  private:
    // Only names wit _(name) come from the steering file.
    // Collection names
	std::string _trackCollectionName;
    std::string _hitCollectionName;

	std::vector<int> _DUTIDs;
	std::map<int, std::vector<const double *>> _DUTHits;
	std::map<int, int> _DUTMatchedTrackes; 
	int _totalMatchedTracks;

//    bool _undoAlignment;

  }; // close class declaration

  //! A global instance of the processor
  EUTelProcessorDUTGlobalEfficiency gEUTelProcessorDUTGlobalEfficiency;

} // close eutelescope namespace scope

#endif // Gear check ifndef
#endif // EUTelLocaltoGlobalHitMaker ifndef
