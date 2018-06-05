/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

#ifndef ALIBAVA_NAMESPACE_H
#define ALIBAVA_NAMESPACE_H

namespace alibava { }
#endif

#ifndef ALIBAVA_H
#define ALIBAVA_H

#ifdef USE_MARLIN
// streamlog include
#include "streamlog/streamlog.h"
#endif

// system includes <>
#include <iostream>
#include <string>
#include <vector>

namespace alibava
{

    // Global constants used in the Alibava package

    class ALIBAVA
    {

	public:

	    // Default destructor.
	    virtual ~ ALIBAVA ( )
	    {
		;
	    }

	public:

	    // System Parameters
	    // If you do change these, you might need to change things in AlibavaConverter processor

	    // Parameter to store/recall number of chips in the Alibava readout system
	    static const int NOOFCHIPS = 2;

	    // Parameter to store/recall number of channels in a chip in the Alibava readout system
	    static const int NOOFCHANNELS = 128;

	    // Parameter to store/recall length of a header in a chip in the Alibava readout system
	    static const int CHIPHEADERLENGTH = 16;

	    // Global Alibava Processor parameters
	    static const char * CHANNELSTOBEUSED;
	    static const char * SKIPMASKEDEVENTS;

	    // General Parameters
	    // Parameter key to store/recall the geometry identification number
	    static const char * GEOID;
	    static const char * NOTSET;

	    // Optional Parameters

	    // Data Header Parameters
	    // Parameter key to store/recall the header
	    static const char * HEADER;

	    // Parameter key to store/recall the header version number
	    static const char * HEADERVERSION;

	    // Parameter key to store/recall the number of events in the file
	    static const char * NOOFEVENT;
	
	    // Parameter key to store/recall the data type: Pedestal, Calibration, or BeamData
	    static const char * DATATYPE;

	    // Parameter to store/recall date of data
	    static const char * DATETIME;

	    // Parameter to store/recall pedestal values stored in Alibava run header
	    static const char * HEADERPEDESTAL;

	    // Parameter to store/recall noise values stored in Alibava run header
	    static const char * HEADERNOISE;

	    // Parameter to store/recall the applied processors
	    static const char * APPLIEDPROCESSOR;

	    // Parameter to store/recall if there is chip selection
	    static const char * ISTHERESELECTION;

	    // Parameter to store/recall which chip selected
	    static const char * SELECTEDCHIPNUM;

	    // Event Parameters
	    // Parameter key to store/recall the event type
	    static const char * EVENTTYPE;

	    // Parameter key to store/recall the "value" in Alibava event
	    static const char *   EVENTVALUE;

	    // Parameter key to store/recall event size
	    static const char *   EVENTSIZE;

	    // Parameter key to store/recall Event header code
	    static const char *   EVENTCODE;

	    // Parameter key to store/recall clock time (timestamp)
	    static const char *   EVENTCLOCK;

	    // Parameter key to store/recall TDC time (ns)
	    static const char *   EVENTTIME;

	    // Parameter key to store/recall Temperature in degrees
	    static const char *   EVENTTEMP;

	    // Parameter key to store/recall injected charge value (e) in Charge calibration run
	    static const char *   CALCHARGE;

	    // Parameter key to store/recall Delay value (ns) in Delay calibration run
	    static const char *   CALDELAY;

	    // Parameter key to store/recall the event mask
	    static const char * EVENTMASK;

	    // AlibavaData CellIDEncoding
	    // Parameter key to store/recall the CellIDEncoding string for all TrackerData before clustering
	    static const char * ALIBAVADATA_ENCODE;

	    // Parameter key to store/recall the part of CellIDEncode that returns chip number
	    static const char * ALIBAVADATA_ENCODE_CHIPNUM;

	    // AlibavaCluster CellIDEncoding
	    // Parameter key to store/recall the CellIDEncoding string for AlibavaCluster
	    static const char *   ALIBAVACLUSTER_ENCODE;

	    // Parameter key to store/recall the part of CellIDEncode that returns cluster ID
	    static const char *   ALIBAVACLUSTER_ENCODE_CLUSTERID;

	    // Parameter key to store/recall the part of CellIDEncode that returns if X axis is sensitive
	    static const char *   ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX;

	    // Parameter key to store/recall the part of CellIDEncode that returns if signal is negative
	    static const char *   ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE;

	    // Parameter key to store/recall the part of CellIDEncode that returns chip number
	    static const char *   ALIBAVACLUSTER_ENCODE_CHIPNUM;

	    // Parameter key to store/recall the part of CellIDEncode that returns seed channel number
	    static const char *   ALIBAVACLUSTER_ENCODE_SEED;

	    // Parameter key to store/recall the part of CellIDEncode that returns cluster size
	    static const char *   ALIBAVACLUSTER_ENCODE_CLUSTERSIZE;

    };

    /*
    enum DataType
    {
	kUNKNOWN = 0,
	kPedestal = 1,
	kChargeCal = 2,
	kDelayCal = 3,
	kBeamData = 4
    };
    */

    enum EventType
    {
	kNewFile = 0,
	kStartOfRun = 1,
	kDataBlock = 2,
	kCheckPoint = 3,
	kEndOfRun = 4
    };

} // end of alibava namespace

std::string trim_right ( const std::string &s );

std::string trim_left ( const std::string &s );

std::string trim_str ( const std::string &s );

std::string getSubStringUpToChar ( std::string search_string, const char* search_char, std::size_t start_position );

template < class T > std::string to_string ( const T& t )
{
    std::ostringstream ss;
    ss << t;
    return ss.str ( );
}

#endif
