#include "FEI4Staggered.h"

namespace eutelescope {
namespace geo {

FEI4Staggered::FEI4Staggered(): EUTelGenericPixGeoDescr(	20.00, 16.8, 0.025,		//size X, Y, Z
													0, 79, 0, 335,			//min max X,Y
													93.660734 )				//rad length					
{
	//Create the material for the sensor
	matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, _radLength, 45.753206 );
	Si = new TGeoMedium("FEI4Silicon",1, matSi);
	/* Make a box for the sensitive area
	Size is: x=2*400+78*250=20300 microns and y=336*50=16800 microns
	MakeBox takes the half of those values in mm as arguments */
	plane = _tGeoManager->MakeBox( "sensarea_fei4", Si, 10.0, 8.4, 0.0125 );
	
	//Prepare the repeating two row structure
	TGeoVolume* oneDoubleRow = _tGeoManager->MakeBox( "row", Si, 10.05, 0.05, 0.0125 );
	TGeoVolume* upperRow = _tGeoManager->MakeBox( "uRow", Si, 10.0, 0.025, 0.0125 );	
	TGeoVolume* lowerRow = _tGeoManager->MakeBox( "lRow", Si, 10.0, 0.025, 0.0125 );
	TGeoVolume* edgePixle = _tGeoManager->MakeBox( "ePx", Si, 0.375/2, 0.025, 0.0125 );
	
	//divisions
	upperRow->Divide("uPix", 1, 80, 0, 1, 0, "N");
	lowerRow->AddNode(edgePixle, 1, new TGeoTranslation(-10+0.375/2,0,0));
	lowerRow->AddNode(edgePixle, 2, new TGeoTranslation(+10-0.375/2,0,0));
	lowerRow->Divide("lPix", 1, 77, -10+0.375, 0.250, 0, "");
	
	//placement
	oneDoubleRow->AddNode(upperRow, 1, new TGeoTranslation(0,0.025,0));
	oneDoubleRow->AddNode(lowerRow, 1, new TGeoTranslation(0,-0.025,0));
			
	//Create plane which we divide in 168 rows (called largeRows);
	TGeoVolume* largeRows = plane->Divide("largerow", 2, 168, 0, 1, 0, "N");
	largeRows->AddNode(oneDoubleRow, 1);
}

FEI4Staggered::~FEI4Staggered()
{
	//delete matSi;
	//delete Si;
}

void  FEI4Staggered::createRootDescr(char const * planeVolume)
{
	//Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
	TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
	//Finaly add the sensitive area to the plane
	topplane->AddNode(plane, 1);
}

std::string FEI4Staggered::getPixName(int x , int y)
{
	char buffer [100];

	//since pixel 0|0 is located on the upper left corner we have to correct y by 335-y+1 
	//(one for the offset in TGeo which starts counting at 1)
	if (x == 0 )
	{
		snprintf( buffer, 100, "/sensarea_fei4_1/fei4edgeregion_1/fei4edgepixel_%d", 336-y);
	}

	else if ( x == 79 )
	{
		snprintf( buffer, 100, "/sensarea_fei4_1/fei4edgeregion_2/fei4edgepixel_%d", 336-y);
	}
	if(x > 0 && x < 79 )
	{
	
		snprintf( buffer, 100, "/sensarea_fei4_1/fei4centreregion_1/fei4centrerow_%d/fei4centrepixel_%d", 336-y, x);
	}
	//Return the full path
	return std::string( buffer ); 
}

//TODO: parse the path to a pixel number!
std::pair<int, int>  FEI4Staggered::getPixIndex(char const*){return std::make_pair(0,0); }

EUTelGenericPixGeoDescr* maker()
{
	FEI4Staggered* mPixGeoDescr = new FEI4Staggered();
	return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
}

} //namespace geo
} //namespace eutelescope

