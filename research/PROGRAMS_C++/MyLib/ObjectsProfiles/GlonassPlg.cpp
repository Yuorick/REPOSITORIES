//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include "GlonassPlg.h"
#include "URPolygon.h"
#include "URPointXY.h"
void drawGlonass( wchar_t *pwchShapeFileOut, double valRastigenie, double valXSdvig, double valYSdvig )
{
 

  TURPointXY pntSdvig( valXSdvig, valYSdvig);
  TURPolygon plgRez(NumPartsGlon, NumPointsGlon,PartsGlon,PointsGlon);
  TURPolygon plgRez1 =   plgRez.SdvigTransform(pntSdvig );
  plgRez1.WriteSetSHPFiles(pwchShapeFileOut, &plgRez1 , 1);

}
#pragma package(smart_init)
