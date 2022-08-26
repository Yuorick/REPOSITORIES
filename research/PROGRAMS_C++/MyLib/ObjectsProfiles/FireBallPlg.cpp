//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include "URPolygon.h"
#include "URPointXY.h"
#include "FireBallPlg.h"

void drawFireBall( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig,const double valRastigenie)
{
	// полигон молнии
	TURPolygon PlgFB(NumPartsFB, NumPointsFB,PartsFB,PointsFB);


  TURPolygon PlgFB1 = PlgFB.LinTransform( valAngRotation,  pntSdvig  , valRastigenie );

   PlgFB1.WriteSetSHPFiles(pwcharrFileOut, &PlgFB1,1);

}
#pragma package(smart_init)
