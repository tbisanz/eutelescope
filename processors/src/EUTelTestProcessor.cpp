/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelTestProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelSparseClusterImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

#include "marlin/AIDAProcessor.h"
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;
using namespace marlin;

EUTelTestProcessor::EUTelTestProcessor():
Processor("EUTelTestProcessor"),
_hitCollectionEUTInput(),
_hitCollectionMCInput()	
{
		_description ="input collection from EUTelescope hitmaker and mc_hit";
		registerInputCollection(LCIO::TRACKERHIT, "hitCollectionEUTInput", "EUTelescope Input Collection", _hitCollectionEUTInput, std::string("hit"));
                registerInputCollection(LCIO::TRACKERHIT, "hitCollectionMCInput", "Monte Carlo Input Collection", _hitCollectionMCInput, std::string("mc_hit"));
}


void EUTelTestProcessor::init() {
		geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);

		_sensorIDVec = geo::gGeometry().sensorIDsVec();
		bookHistos();

		std::cout << "EUTelTestProcessor::init()" << std::endl;
}

void EUTelTestProcessor::processRunHeader(LCRunHeader* rdr)
{
	
}

void EUTelTestProcessor::processEvent(LCEvent* event)
{
		//Opens collection for input.
		LCCollection* inputCollection = nullptr;
		EUTelEventImpl *evt = static_cast<EUTelEventImpl *>(event);
	//	try {
				inputCollection = evt->getCollection(_hitCollectionEUTInput);
	//	} catch (DataNotAvailableException& e) {
	//			streamlog_out( WARNING2 ) << _hitCollectionNameInput << " collection not available" << std::endl;
	//			return;
	//	}
		LCCollection* inputCollection2 = nullptr;
		
				inputCollection2 = evt->getCollection(_hitCollectionMCInput);
		/*std::string encoding = inputCollection->getParameters().getStringVal( LCIO::CellIDEncoding );
		if(encoding.empty()) {
			encoding = EUTELESCOPE::HITENCODING;
		}*/
		lcio::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
		
		//preparing maps for input comparison

		std::vector<TrackerDataImpl*> registeredRawData;
		std::map<TrackerDataImpl*, TrackerHitImpl*> map1;
		std::map<TrackerDataImpl*, TrackerHitImpl*> map2;

		//filling maps

		for(int iHit = 0; iHit < inputCollection->getNumberOfElements(); ++iHit) {
                        TrackerHitImpl* inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));


                        auto objectVector1 = inputHit->getRawHits();
                        auto const & cluster_zs_data1 = static_cast<TrackerDataImpl*>(objectVector1[0]);
                        map1.insert(std::pair<TrackerDataImpl*, TrackerHitImpl*>(cluster_zs_data1, inputHit) );
		}


		
		for(int iHit = 0; iHit < inputCollection2->getNumberOfElements(); ++iHit) {
			TrackerHitImpl* inputHit2 = static_cast<TrackerHitImpl*>(inputCollection2->getElementAt(iHit));

			
			auto objectVector2 = inputHit2->getRawHits();
			auto const & cluster_zs_data2 = static_cast<TrackerDataImpl*>(objectVector2[0]);
			map2.insert(std::pair<TrackerDataImpl*, TrackerHitImpl*>(cluster_zs_data2, inputHit2));
				
		}

		//compare number of hits in the collections

		if(inputCollection->getNumberOfElements() >= inputCollection2->getNumberOfElements()){
			for(int iHit = 0; iHit < inputCollection->getNumberOfElements(); ++iHit) {
				
				TrackerHitImpl* inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));
				TrackerHitImpl* inputHit2 = static_cast<TrackerHitImpl*>(inputCollection2->getElementAt(iHit));
				int sensorID = hitDecoder(inputHit)["sensorID"];
				auto objectVector1 = inputHit->getRawHits();
                        	auto const & cluster_zs_data1 = static_cast<TrackerDataImpl*>(objectVector1[0]);
				if(map2.find(cluster_zs_data1) != map2.end()){
					//if hit in both collections compute difference in coordinates
				
					const double* inputPos = map1[cluster_zs_data1]->getPosition();
					const double* inputPos2 = map2[cluster_zs_data1]->getPosition();
					std::array<double, 3> inputArray1 = { inputPos[0], inputPos[1], inputPos[2] };
					std::array<double, 3> inputArray2 = { inputPos2[0], inputPos2[1], inputPos2[2] };

					int cluX = 0, cluY = 0;
			                float locX = 0, locY = 0;
					int clustersize = 0;	

        			        auto rawData1 = static_cast<TrackerDataImpl*>(inputHit->getRawHits()[0]);
       				        //if( inputHit->getType() == kEUTelSparseClusterImpl ){
                        			auto const & cluster1 = EUTelSparseClusterImpl<EUTelGenericSparsePixel>(rawData1);
                        			cluster1.getClusterSize(cluX, cluY);
                        			cluster1.getCenterOfGravity(locX, locY);
                        			clustersize = cluster1.size();
                			//}

					double difX = (inputArray1[0] - inputArray2[0])*1000;
                                        double difY = (inputArray1[1] - inputArray2[1])*1000;
                                        double difZ = (inputArray1[2] - inputArray2[2])*1000;
					double dif = sqrt(difX * difX+difY * difY+difZ*difZ);
					
						if(clustersize == 1){
							_HitDifXCluster1Histo.at(sensorID)->fill(difX);
							//_HitDifXCluster1Histo.at(sensorID)->SetFillColor(kRed);
							//_HitDifXCluster1Histo.at(sensorID)->SetMarkerStyle(21);
							//_HitDifXCluster1Histo.at(sensorID)->SetMarkerColor(kRed);

							_HitDifYCluster1Histo.at(sensorID)->fill(difY);
							//_HitDifYCluster1Histo.at(sensorID)->SetFillColor(kRed);
                                                        //_HitDifYCluster1Histo.at(sensorID)->SetMarkerStyle(21);
                                                        //_HitDifYCluster1Histo.at(sensorID)->SetMarkerColor(kRed);
						}
						if(clustersize == 2){
							if(cluX == 2){
								_HitDifXCluster2Histo.at(sensorID)->fill(difX);
								//_HitDifXCluster2Histo.at(sensorID)->SetFillColor(kBlue);
								//_HitDifXCluster2Histo.at(sensorID)->SetMarkerStyle(21);
								//_HitDifXCluster2Histo.at(sensorID)->SetMarkerColor(kBlue);
							}
							if(cluY == 2){
								_HitDifYCluster2Histo.at(sensorID)->fill(difY);
								//_HitDifYCluster2Histo.at(sensorID)->SetFillColor(kBlue);
                                                                //_HitDifYCluster2Histo.at(sensorID)->SetMarkerStyle(21);
                                                                //_HitDifYCluster2Histo.at(sensorID)->SetMarkerColor(kBlue);
							}
						}
						if(clustersize >= 3){
							if(cluX >= 2 && cluY >= 1){
								_HitDifXCluster3Histo.at(sensorID)->fill(difX);
								//_HitDifXCluster3Histo.at(sensorID)->SetFillColor(kGreen);
                                                                //_HitDifXCluster3Histo.at(sensorID)->SetMarkerStyle(21);
                                                                //_HitDifXCluster3Histo.at(sensorID)->SetMarkerColor(kGreen);
							}
							if(cluY >= 2 && cluX >=1){
								_HitDifYCluster3Histo.at(sensorID)->fill(difY);
								//_HitDifYCluster3Histo.at(sensorID)->SetFillColor(kGreen);
                                                                //_HitDifYCluster3Histo.at(sensorID)->SetMarkerStyle(21);
                                                                //_HitDifYCluster3Histo.at(sensorID)->SetMarkerColor(kGreen);
							}
						}	

					_HitDifHisto.at(sensorID)->fill(dif);

					_HitXDifHisto.at(sensorID)->fill(difX);
					//_HitXDifHisto.at(sensorID)->add(_HitDifXCluster1Histo.at(sensorID));
					//_HitXDifHisto.at(sensorID)->add(_HitDifXCluster2Histo.at(sensorID));
					//_HitXDifHisto.at(sensorID)->add(_HitDifXCluster3Histo.at(sensorID));

					_HitYDifHisto.at(sensorID)->fill(difY);
                                        //_HitYDifHisto.at(sensorID)->add(_HitDifYCluster1Histo.at(sensorID));
                                        //_HitYDifHisto.at(sensorID)->add(_HitDifYCluster2Histo.at(sensorID));
                                        //_HitYDifHisto.at(sensorID)->add(_HitDifYCluster3Histo.at(sensorID));
					

                                        _HitZDifHisto.at(sensorID)->fill(difZ);
				}
				else{
					streamlog_out(WARNING1) << "missmatch between number of hits in monte carlo and EUTelescope data sample" << std::endl;
				}
			}
		}
		if(inputCollection->getNumberOfElements() < inputCollection2->getNumberOfElements()){
			for(int iHit = 0; iHit < inputCollection2->getNumberOfElements(); ++iHit) {
                                TrackerHitImpl* inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));
                                TrackerHitImpl* inputHit2 = static_cast<TrackerHitImpl*>(inputCollection2->getElementAt(iHit));
                                int sensorID = hitDecoder(inputHit2)["sensorID"];
                                auto objectVector2 = inputHit2->getRawHits();
                                auto const & cluster_zs_data2 = static_cast<TrackerDataImpl*>(objectVector2[0]);
                                if(map1.find(cluster_zs_data2) != map1.end()){
					//if hit in both collections compute difference in coordinates

                                        const double* inputPos = map1[cluster_zs_data2]->getPosition();
                                        const double* inputPos2 = map2[cluster_zs_data2]->getPosition();
                                        std::array<double, 3> inputArray1 = { inputPos[0], inputPos[1], inputPos[2] };
                                        std::array<double, 3> inputArray2 = { inputPos2[0], inputPos2[1], inputPos2[2]};

                                        int cluX = 0, cluY = 0;
                                        float locX = 0, locY = 0;
					int clustersize = 0;

                                        auto rawData1 = static_cast<TrackerDataImpl*>(inputHit->getRawHits()[0]);
                                        //if( inputHit->getType() == kEUTelSparseClusterImpl ){
                                                auto const & cluster1 = EUTelSparseClusterImpl<EUTelGenericSparsePixel>(rawData1);
                                                cluster1.getClusterSize(cluX, cluY);
                                                cluster1.getCenterOfGravity(locX, locY);
                                                clustersize = cluster1.size();
                                        //}


                                        double difX = (inputArray1[0] - inputArray2[0])*1000;
                                        double difY = (inputArray1[1] - inputArray2[1])*1000;
                                        double difZ = (inputArray1[2] - inputArray2[2])*1000;
                                        double dif = sqrt(difX * difX+difY * difY+difZ*difZ);

                                                if(clustersize == 1){
                                                        _HitDifXCluster1Histo.at(sensorID)->fill(difX);
							//_HitDifXCluster1Histo.at(sensorID)->SetFillColor(kRed);
							//_HitDifXCluster1Histo.at(sensorID)->SetMarkerStyle(21);
							//_HitDifXCluster1Histo.at(sensorID)->SetMarkerColor(kRed);

                                                        _HitDifYCluster1Histo.at(sensorID)->fill(difY);
							//_HitDifYCluster1Histo.at(sensorID)->SetFillColor(kRed);
                                                        //_HitDifYCluster1Histo.at(sensorID)->SetMarkerStyle(21);
                                                        //_HitDifYCluster1Histo.at(sensorID)->SetMarkerColor(kRed);
                                                }
                                                if(clustersize == 2){
							if(cluX == 2){
                                                        	_HitDifXCluster2Histo.at(sensorID)->fill(difX);
								//_HitDifXCluster2Histo.at(sensorID)->SetFillColor(kBlue);
								//_HitDifXCluster2Histo.at(sensorID)->SetMarkerStyle(21);
								//_HitDifXCluster2Histo.at(sensorID)->SetMarkerColor(kBlue);
							}
							if(cluY == 2){
                                                        	_HitDifYCluster2Histo.at(sensorID)->fill(difY);
								//_HitDifYCluster2Histo.at(sensorID)->SetFillColor(kBlue);
                                                                //_HitDifYCluster2Histo.at(sensorID)->SetMarkerStyle(21);
                                                                //_HitDifYCluster2Histo.at(sensorID)->SetMarkerColor(kBlue);
							}
                                                }
                                                if(clustersize >= 3){
							if(cluX >= 2 && cluY >= 1){
                                                        	_HitDifXCluster3Histo.at(sensorID)->fill(difX);
								//_HitDifXCluster3Histo.at(sensorID)->SetFillColor(kGreen);
								//_HitDifXCluster3Histo.at(sensorID)->SetMarkerStyle(21);
								//_HitDifXCluster3Histo.at(sensorID)->SetMarkerColor(kGreen);
							}
							if(cluY >= 2 && cluX >= 1){
                                                       		_HitDifYCluster3Histo.at(sensorID)->fill(difY);
								//_HitDifYCluster3Histo.at(sensorID)->SetFillColor(kGreen);
                                                                //_HitDifYCluster3Histo.at(sensorID)->SetMarkerStyle(21);
                                                                //_HitDifYCluster3Histo.at(sensorID)->SetMarkerColor(kGreen);
							}
                                                }
                                                
                                        _HitDifHisto.at(sensorID)->fill(dif);

					_HitXDifHisto.at(sensorID)->fill(difX);
					//_HitXDifHisto.at(sensorID)->add(_HitDifXCluster1Histo.at(sensorID));
					//_HitXDifHisto.at(sensorID)->add(_HitDifXCluster2Histo.at(sensorID));
					//_HitXDifHisto.at(sensorID)->add(_HitDifXCluster3Histo.at(sensorID));
					
					_HitYDifHisto.at(sensorID)->fill(difY);
					//_HitYDifHisto.at(sensorID)->add(_HitDifYCluster1Histo.at(sensorID));
					//_HitYDifHisto.at(sensorID)->add(_HitDifYCluster2Histo.at(sensorID));
					//_HitYDifHisto.at(sensorID)->add(_HitDifYCluster3Histo.at(sensorID));

					_HitZDifHisto.at(sensorID)->fill(difZ);
                                }       	
				else{   	
					streamlog_out(WARNING1) << "missmatch between number of hits in monte carlo and EUTelescope data sample" << std::endl;
				}       	
                        }               	
                                                
		}                       	
}                                       	
                                        	
void EUTelTestProcessor::end()          	
{                                       	
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}                                       	
                                                
void EUTelTestProcessor::bookHistos() { 	
       	streamlog_out(MESSAGE1) << "Booking histograms " << std::endl;
                                                
 	 for (auto sensorID : _sensorIDVec) {
   		 auto basePath = "detector_" + to_string(sensorID);
		 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
  		 basePath.append("/");  	
                                        	
		 auto const countHistoName = "absolut hitposition difference" ;
                                        	
		 auto HitDifHisto =     	
       			 marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
           		 (basePath + countHistoName).c_str(), 100, 0, 100);
		 HitDifHisto->setTitle("absolut hitposition difference  ;difference in hit position [#mum];number of events");
                                        	
		 _HitDifHisto[sensorID] = HitDifHisto;
	}                               	
	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID);
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");  	
                                        	
                 auto const countHistoName = "hitposition X-difference";

                 auto HitXDifHisto =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 100, -20, 20);
		 HitXDifHisto->setTitle("hitposition X-difference  ;difference in hit position [#mum];number of events");

                 _HitXDifHisto[sensorID] = HitXDifHisto;
        }
	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID);
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "hitposition Y-difference";

                 auto HitYDifHisto =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 100, -20, 20);
		 HitYDifHisto->setTitle("hitposition Y-difference  ;difference in hit position [#mum];number of events");

                 _HitYDifHisto[sensorID] = HitYDifHisto;
        }
	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID);
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "hitposition Z-difference";

                 auto HitZDifHisto =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 100, -100, 100);
		 HitZDifHisto->setTitle("hitposition Z-difference  ;difference in hit position [#mum];number of events");

                 _HitZDifHisto[sensorID] = HitZDifHisto;
        }
	
	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID) + "/clustersize_1";
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "clustersize 1 X-difference";

                 auto HitDifXCluster1Histo =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 600, -150, 150);
                 HitDifXCluster1Histo->setTitle("clustersize 1 X-difference  ;difference in hit position [#mum];number of events");

                 _HitDifXCluster1Histo[sensorID] = HitDifXCluster1Histo;
        }
	
	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID) + "/clustersize_1";
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "clustersize 1 Y-difference";

                 auto HitDifYCluster1Histo =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 300, -100, 100);
                 HitDifYCluster1Histo->setTitle("clustersize 1 Y-difference  ;difference in hit position [#mum];number of events");

                 _HitDifYCluster1Histo[sensorID] = HitDifYCluster1Histo;
        }

	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID) + "/clustersize_2";
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "clustersize 2 X-difference";

                 auto HitDifXCluster2Histo =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 600, -150, 150);
                 HitDifXCluster2Histo->setTitle("clustersize 2 X-difference  ;difference in hit position [#mum];number of events");

                 _HitDifXCluster2Histo[sensorID] = HitDifXCluster2Histo;
        }

	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID) + "/clustersize_2";
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "clustersize 2 Y-difference";

                 auto HitDifYCluster2Histo =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 300, -100, 100);
                 HitDifYCluster2Histo->setTitle("clustersize 2 Y-difference  ;difference in hit position [#mum];number of events");

                 _HitDifYCluster2Histo[sensorID] = HitDifYCluster2Histo;
        }

	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID) + "/clustersize_>3";
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "clustersize >3 X-difference";

                 auto HitDifXCluster3Histo =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 100, -20, 20);
                 HitDifXCluster3Histo->setTitle("clustersize >3 X-difference  ;difference in hit position [#mum];number of events");

                 _HitDifXCluster3Histo[sensorID] = HitDifXCluster3Histo;
        }

	for (auto sensorID : _sensorIDVec) {
                 auto basePath = "detector_" + to_string(sensorID) + "/clustersize_>3";
                 marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
                 basePath.append("/");

                 auto const countHistoName = "clustersize >3 Y-difference";

                 auto HitDifYCluster3Histo =
                         marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(
                         (basePath + countHistoName).c_str(), 100, -20, 20);
                 HitDifYCluster3Histo->setTitle("clustersize >3 Y-difference  ;difference in hit position [#mu];number of events");

                 _HitDifYCluster3Histo[sensorID] = HitDifYCluster3Histo;
        }
}

	
