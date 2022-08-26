//---------------------------------------------------------------------------


#pragma hdrstop

#include "JetPlg.h"

#include "URPolygon.h"
#include "URPointXY.h"


 void drawJet( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig,const double valRastigenie)
{

	TURPolygon PlgJet(NumPartsJet, NumPointsJet, PartsJet, PointsJet);


  TURPolygon PlgJet1 = PlgJet.LinTransform( valAngRotation,  pntSdvig  , valRastigenie );

   PlgJet1.WriteSetSHPFiles(pwcharrFileOut, &PlgJet1,1);

}

// перегруженная , с опцией зеркального отражения
void drawJet( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig
   ,const double valRastigenie, bool bOtragenie)
{

	TURPolygon PlgJet(NumPartsJet, NumPointsJet, PartsJet, PointsJet);
	if(bOtragenie)
	{
		PlgJet =PlgJet.SimOtragenieTransform();
	}



  TURPolygon PlgJet1 = PlgJet.LinTransform( valAngRotation,  pntSdvig  , valRastigenie );


   PlgJet1.WriteSetSHPFiles(pwcharrFileOut, &PlgJet1,1);

}


#pragma package(smart_init)
