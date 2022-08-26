//---------------------------------------------------------------------------

#ifndef ManulFunctionsH
#define ManulFunctionsH
#include "EtalonSign.h"
class  ETalonSign;
class TFar_2D;
class TManulFunctions
{

public:

static void __fastcall TManulFunctions::createSKZ_Graphs(wchar_t* wchOutFold, const double VAlTargH, const double VAlTargEPR
  , const TEtalonSign ETalonSign ,const double PowerPrd,const double KYPrd, const double VAlHAntenna, const double VAlZAhvatDist
 , TFar_2D Far_2D, const double VAlSigW, const double VAlSigEpsMMO, const double VAlMeasT
  , const double VAlVelocity0);
};
#endif
