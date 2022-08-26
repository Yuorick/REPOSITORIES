//---------------------------------------------------------------------------

#ifndef EmitterImageH
#define EmitterImageH
#include "URPolyLine.h"
#include "UrPointXY.h"
class TURPolyLine;
class TURPointXY;
class TEmitterImage
{
public:
 // отрезки
	TURPolyLine marrPlnSegm[3];
	TURPointXY  mpntCentre;

	__fastcall TEmitterImage();
	// Конструктор копирования
	TEmitterImage (const TEmitterImage &R2);
	// оператор присваивания
	TEmitterImage operator=(TEmitterImage  R2);
	__fastcall TEmitterImage(const double ValLength0 ,const double ValLength1,const double ValFi) ;

	 void createShapeFiles(wchar_t *wchFoldName) ;

 static	TEmitterImage   SdvigEmTransform( TEmitterImage ImageInp, TURPointXY pntSdvig );

};
#endif
