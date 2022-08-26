//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include "VesselPlg.h"
#include "URPolygon.h"
#include "URPointXY.h"

void drawVessel( wchar_t *pwchShapeFileOut, const double valPos, const double valVeloc, const double valRastigenie )
{


	///
	// полигон корабля
	TURPolygon PlgVess(NumPartsVess, NumPointsVess,PartsVess,PointsVess);
	TURPointXY pntSdvig( 0., 0.);
	if (valVeloc < 0.)
	{
	  PlgVess = PlgVess.LinTransform( M_PI,  pntSdvig  , 1. );
	  PlgVess = PlgVess.SimOtragenieTransform();
    }
  pntSdvig.X =  valPos;
  TURPolygon PlgVess1 = PlgVess.LinTransform( 0.,  pntSdvig  , valRastigenie );

   PlgVess1.WriteSetSHPFiles(pwchShapeFileOut, &PlgVess1,1);


}

#pragma package(smart_init)
