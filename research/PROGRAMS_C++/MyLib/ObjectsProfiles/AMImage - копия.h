//---------------------------------------------------------------------------

#ifndef AMImageH
#define AMImageH
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "EmitterImage.h"
#include "URPolygon.h"
class TURPolyLine;
class TURPointXY;
class TEmitterImage;
class TURPolygon ;
class TAMImage
{
public:


 // ����� ������
	TEmitterImage *mparrEmit ;
	int mquantEmit;
	TURPointXY mpntCentre;
	TURPolygon mplgHood;
	TURPolyLine mplnOut;

 // �������-�����
		__fastcall ~TAMImage() ;
 //	__fastcall TAMImage();
	// ����������� �����������
	TAMImage (const TAMImage &R2);
	// �������� ������������
	TAMImage operator=(TAMImage  R2);

	 __fastcall TAMImage( TEmitterImage Emit ,const int QuantEmit,const double DistEmit,
 const double SideWidth, const double BackWidth , const double LenghPlnOut );

   //	__fastcall TAMImage(const double ValLength0 ,const double ValLength1,const double ValFi) ;

   void  __fastcall createShapeFiles(wchar_t *wchFoldName);




};
#endif
