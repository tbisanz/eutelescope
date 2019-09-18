/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCALICE_H
#define EUTELCALICE_H

// eutelescope includes ".h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// AIDA includes <.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram1D.h>

// system includes <>
#include <string>
#include <vector>

namespace eutelescope {

  class EUTelCALICE : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelCALICE)
    static std::pair<int,int> mapChipSiPM(int chip, int sipm);
    virtual void bookHistos();

  public:
    // Returns a new instance of EUTelCALICE
    virtual Processor *newProcessor() {
      return new EUTelCALICE;
    }

    // Default constructor
    EUTelCALICE();
    
    // Called only at the begining of a job
    virtual void init();

    // Called every run
    virtual void processRunHeader(LCRunHeader*){};

    // Called every event.
    virtual void processEvent(LCEvent *event);

    // Called at the end of the job
    virtual void end();

  private:
    std::string _CALICERawDataCollectionName;
    AIDA::IHistogram2D* hitsVsTracksAll;
    AIDA::IHistogram1D* nTracksHist;
    std::map<int,AIDA::IHistogram2D*> _hitsVsTracksHistos;
  };

  //! A global instance of the processor
  EUTelCALICE gEUTelCALICE;
}
#endif
