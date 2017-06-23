/*
 *   This processor removes noisy pixels from a TrackerData collection
 *   been masked as noisy and written out into a noisy pixel database
 *
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorTimeCut.h"
#include "EUTELESCOPE.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "EUTelRunHeaderImpl.h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <LCIOTypes.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>

// system includes
#include <algorithm>
#include <memory>

namespace eutelescope {

  EUTelProcessorTimeCut::EUTelProcessorTimeCut()
      : Processor("EUTelProcessorTimeCut"), _inputCollectionName(""),
        _outputCollectionName(""), _sensorIDCut(0), _timeCut(0) {
    _description = "EUTelProcessorTimeCut removes noisy pixels "
                   "(TrackerData) from a collection. This processor requires a "
                   "noisy pixel collection.";

    registerInputCollection(LCIO::TRACKERDATA, "InputCollectionName",
                            "Input collection containing noisy raw data",
                            _inputCollectionName,
                            std::string("noisy_raw_data_collection"));
    registerOutputCollection(
        LCIO::TRACKERDATA, "OutputCollectionName",
        "Output collection where noisy pixels have been removed",
        _outputCollectionName, std::string("noisefree_raw_data_collection"));
 

		registerOptionalParameter("sensorIDCut", "TODO", _sensorIDCut, int(20));
		registerOptionalParameter("timeCut", "TODO", _timeCut, int(0));

 }

  void EUTelProcessorTimeCut::init() {
    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters();
  }

  void EUTelProcessorTimeCut::processRunHeader(LCRunHeader *rdr) {
    auto runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
    runHeader->addProcessor(type());
  }

  void EUTelProcessorTimeCut::processEvent(LCEvent *event) {

    // get the collection of interest from the event.
    LCCollectionVec *inputCollection = nullptr;

    try {
      inputCollection = dynamic_cast<LCCollectionVec *>(
          event->getCollection(_inputCollectionName));
    } catch (lcio::DataNotAvailableException &e) {
      return;
    }

    // now prepare output collection
    LCCollectionVec *outputCollection = nullptr;
    bool outputCollectionExists = false;
    size_t initialOutputCollectionSize = 0;

    try {
      outputCollection = dynamic_cast<LCCollectionVec *>(
          event->getCollection(_outputCollectionName));
      outputCollectionExists = true;
      initialOutputCollectionSize = outputCollection->size();
    } catch (lcio::DataNotAvailableException &e) {
      outputCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    }

    // prepare decoder for input data
    CellIDDecoder<TrackerDataImpl> inputDataDecoder(
        EUTELESCOPE::ZSDATADEFAULTENCODING);

    // read the encoding std::string from the input collection
    std::string encodingString =
        inputCollection->getParameters().getStringVal(LCIO::CellIDEncoding);
    outputCollection->parameters().setValue(LCIO::CellIDEncoding,
                                            encodingString);

    auto trackerData = std::make_unique<lcio::TrackerDataImpl>();

    for (size_t iEntry = 0; iEntry < inputCollection->size(); ++iEntry) {

      if (iEntry > 0) {
        outputCollection->push_back(trackerData.release());
        auto newTrackerData = std::make_unique<lcio::TrackerDataImpl>();
        trackerData = std::move(newTrackerData);
      }

      TrackerDataImpl *inputData = dynamic_cast<TrackerDataImpl *>(
          inputCollection->getElementAt(iEntry));

      int sensorID = inputDataDecoder(inputData)["sensorID"];
      //SparsePixelType pixelType = static_cast<SparsePixelType>(
      //    static_cast<int>(inputDataDecoder(inputData)["sparsePixelType"]));

      trackerData->setCellID0(inputData->getCellID0());
      trackerData->setCellID1(inputData->getCellID1());
      trackerData->setTime(inputData->getTime());



	// now prepare the EUTelescope interface to sparsified data.
    auto sparseData = std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(inputData);
    auto &pixelVec = sparseData->getPixels();

      // interface to sparsified data
      auto sparseOutputData = std::make_unique<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>>(trackerData.get());

      for (auto &pixel : pixelVec) {
        if ( !(sensorID == _sensorIDCut && pixel.getTime() != _timeCut) ){
          sparseOutputData->push_back(pixel);
        }
      }
    }
    outputCollection->push_back(trackerData.release());

    // add the collection if we created it and added elements
    if (!outputCollectionExists) {
      if (outputCollection->size() != initialOutputCollectionSize) {
        event->addCollection(outputCollection, _outputCollectionName);
      } else {
        delete outputCollection;
      }
    }
  }

  void EUTelProcessorTimeCut::end() {}
}
