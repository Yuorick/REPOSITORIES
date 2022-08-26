//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include "URPolygon.h"
#include "URPointXY.h"
#include "UABPlg.h"

void drawUAB( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig,const double valRastigenie)
{
	// полигон молнии
	TURPolygon PlgUAB(NumPartsUAB, NumPointsUAB,PartsUAB,PointsUAB);


  TURPolygon PlgUAB1 = PlgUAB.LinTransform( valAngRotation,  pntSdvig  , valRastigenie );

   PlgUAB1.WriteSetSHPFiles(pwcharrFileOut, &PlgUAB1,1);

}
#pragma package(smart_init)
