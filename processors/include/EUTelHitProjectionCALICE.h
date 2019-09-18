/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELHITPROJECTIONCALICE_H
#define EUTELHITPROJECTIONCALICE_H

// eutelescope includes ".h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// AIDA includes <.h>
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram1D.h>

// system includes <>
#include <string>

namespace eutelescope {

  class EUTelHitProjectionCALICE : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelHitProjectionCALICE)
    virtual void bookHistos();

  public:
    // Returns a new instance of EUTelHitProjectionCALICE
    virtual Processor *newProcessor() {
      return new EUTelHitProjectionCALICE;
    }

    // Default constructor
    EUTelHitProjectionCALICE();
    
    // Called only at the begining of a job
    virtual void init();

    // Called every run
    virtual void processRunHeader(LCRunHeader*){};

    // Called every event.
    virtual void processEvent(LCEvent *event);

    // Called at the end of the job
    virtual void end();

  private:
    // Collection names
    std::string _trackCollectionNameInput;
    std::string _caliceRawCollectionNameInput;
    std::map<int,std::map<int,AIDA::IHistogram2D*>> histoMap;

    AIDA::IHistogram2D* ineffHist = nullptr;
    AIDA::IHistogram2D* ADCHist = nullptr;
    AIDA::IHistogram2D* ADCCountHist = nullptr;
    AIDA::IHistogram2D* allTrackCountHist = nullptr;
    int singleTrackRejects = 0;
    int twoHitRejects = 0;
    int singleCaloHitRejects = 0;
  };

  //! A global instance of the processor
  EUTelHitProjectionCALICE gEUTelHitProjectionCALICE;

}
#endif
