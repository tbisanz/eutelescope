/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTelTestProcessor_H
#define EUTelTestProcessor_H 1

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"
#include "EUTelTripletGBLUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <map>
#include <set>
#include <string>
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#include <AIDA/IHistogram1D.h>

namespace eutelescope {

  class EUTelTestProcessor : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelTestProcessor)

  public:
    // Returns a new instance of EUTelTestProcessor
    virtual Processor *newProcessor() {
      return new EUTelTestProcessor;
    }

    EUTelTestProcessor(); // Default constructor
    // Called only at the begining of a job
    virtual void init();

    // Called every run
    virtual void processRunHeader(LCRunHeader *run);

    // Called every event.
    virtual void processEvent(LCEvent *event);

    // Called at the end of the job
    virtual void end();

    virtual void bookHistos();

    //unstaged change
  private:
    std::string _hitCollectionEUTInput;
    std::string _hitCollectionMCInput;
    //std::vector holding sensor IDs
    std::vector<int> _sensorIDVec;
    std::map<int, AIDA::IHistogram1D*> _HitDifHisto;
    std::map<int, AIDA::IHistogram1D*> _HitXDifHisto;
    std::map<int, AIDA::IHistogram1D*> _HitYDifHisto;
    std::map<int, AIDA::IHistogram1D*> _HitZDifHisto;
    std::map<int, AIDA::IHistogram1D*> _HitDifXCluster1Histo;
    std::map<int, AIDA::IHistogram1D*> _HitDifYCluster1Histo;
    std::map<int, AIDA::IHistogram1D*> _HitDifXCluster2Histo;
    std::map<int, AIDA::IHistogram1D*> _HitDifYCluster2Histo;
    std::map<int, AIDA::IHistogram1D*> _HitDifXCluster3Histo;
    std::map<int, AIDA::IHistogram1D*> _HitDifYCluster3Histo;
  }; // close class declaration
  //! A global instance of the processor
  EUTelTestProcessor gEUTelTestProcessor;

} // close eutelescope namespace scope

#endif // EUTelLocaltoGlobalHitMaker ifndef
