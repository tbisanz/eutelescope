/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPROCESSORTRANSFORMFROMGEAR_H
#define EUTELPROCESSORTRANSFORMFROMGEAR_H

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// system includes <>
#include <string>
#include <map>
#include <vector>

//EIGEN includes
#include <Eigen/Core>

namespace eutelescope {

struct GEAREntries
{
	double r1, r2, r3, r4;
	double offX, offY, offZ;
	double alpha, beta, gamma;
};

class EUTelProcessorTransformFromGEAR: public marlin::Processor
{
  public:
	virtual Processor* newProcessor() 
	{
		return new EUTelProcessorTransformFromGEAR;
	}

	//! Default constructor
	EUTelProcessorTransformFromGEAR();
	virtual void init();
	virtual void processRunHeader(LCRunHeader * run);
	virtual void processEvent(LCEvent * evt);
	virtual void end();

  private:
	std::map<int,GEAREntries> _GEAREntriesMap;
	std::map<int, Eigen::Matrix4d> _flipMatrix;
	std::map<int, Eigen::Matrix4d> _offsetMatrix;

  protected:
	std::string _inputHitCollectionName;
	std::string _outputHitCollectionName;
	
	int _iRun;
	int _iEvt;

	gear::SiPlanesParameters* _siPlanesParameters;
	gear::SiPlanesLayerLayout* _siPlanesLayerLayout;

	int _initialOutputCollectionSize;
};

//! A global instance of the processor
 EUTelProcessorTransformFromGEAR gEUTelProcessorTransformFromGEAR;
}
#endif
