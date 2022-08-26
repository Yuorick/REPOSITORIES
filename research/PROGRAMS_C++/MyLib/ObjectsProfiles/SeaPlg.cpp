//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include <math.h>
#include "VesselPlg.h"
#include "URPolygon.h"
#include "URPointXY.h"
#include "SeaPlg.h"
// ������� ������� ����������� � ���� ������
//
//
void drawSeaSin( wchar_t *pwchShapeFileOut, const double valAmpl
   ,const double valOmega, const double valPh0, const double valXMin
   , const double valXMax, const double valDeep, const int numPoints )
{


	TURPointXY  *pPoints  = new  TURPointXY  [numPoints];
	pPoints[numPoints - 2] =  TURPointXY (valXMax, valDeep);
	pPoints[numPoints - 1] =  TURPointXY (valXMin, valDeep);
	pPoints[0] =  pPoints[numPoints - 1];
	double step = (valXMax - valXMin)/ ((double)(numPoints -4));
	for (int i = 0; i < (numPoints-3); i++)
	{
	  double valx = valXMin + ((double)i)*  step;
	  pPoints[i + 1] =  TURPointXY(valx, valAmpl * sin( valOmega *  valx +valPh0));
	}
	// ������� ����
	 TURPolygon PlgSea(  numPoints,pPoints);


   PlgSea.WriteSetSHPFiles(pwchShapeFileOut, &PlgSea,1);
   delete []pPoints;


}

// ������� ������� ����������� � ���� ���������
//
//
void drawSeaPlane( wchar_t *pwchShapeFileOut, const double valXMin
   , const double valXMax, const double valDeep )
{


	TURPointXY  pPoints[5];

	pPoints[0] =  TURPointXY (valXMin, valDeep);
	pPoints[1] =  TURPointXY (valXMin, 0.);
	pPoints[2] =  TURPointXY (valXMax, 0.);
	pPoints[3] =  TURPointXY (valXMax, valDeep);
	pPoints[4] =  pPoints[0];

	// ������� ����
	TURPolygon PlgSea( 5, pPoints);

   PlgSea.WriteSetSHPFiles(pwchShapeFileOut, &PlgSea,1);



}

#pragma package(smart_init)
