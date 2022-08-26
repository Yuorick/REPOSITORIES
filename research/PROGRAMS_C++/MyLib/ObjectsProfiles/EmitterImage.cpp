//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <string.h>
#include "EmitterImage.h"



__fastcall TEmitterImage::TEmitterImage()
{
 TURPointXY pnt[4];
 pnt[0] = TURPointXY( -1.-1./4.,0.);
 pnt[1] = TURPointXY( -1./4.,0.);
 pnt[2] = TURPointXY( 0.,1./4);
 pnt[3] = TURPointXY( 0.,-1./4);
 marrPlnSegm[0] =  TURPolyLine(pnt[0], pnt[1]);
 marrPlnSegm[1] =  TURPolyLine(pnt[1], pnt[2]);
 marrPlnSegm[2] =  TURPolyLine(pnt[1], pnt[3]);
 mpntCentre = TURPointXY(0.,0.);


}
// Конструктор копирования
 TEmitterImage::TEmitterImage (const TEmitterImage &R2)
 {
 memcpy(marrPlnSegm, R2.marrPlnSegm, 3 * sizeof (TURPolyLine));
 mpntCentre = R2.mpntCentre;

 }
// оператор присваивания
 TEmitterImage TEmitterImage::operator=(TEmitterImage  R2)
{
 memcpy(marrPlnSegm, R2.marrPlnSegm, 3 * sizeof (TURPolyLine));
 mpntCentre = R2.mpntCentre;

  return *this ;
}

// парам констр
 __fastcall TEmitterImage::TEmitterImage(const double ValLength0 ,const double ValLength1,const double ValFi)

 {
	mpntCentre = TURPointXY(0.,0.);
	TURPointXY pnt[4];
	pnt[0] = TURPointXY( -ValLength0 - ValLength1* cos(ValFi) ,0.);
	pnt[1] = TURPointXY( - ValLength1* cos(ValFi),0.);
	pnt[2] = TURPointXY( 0.,ValLength1* sin(ValFi));
	pnt[3] = TURPointXY( 0., -ValLength1* sin(ValFi));
	marrPlnSegm[0] =  TURPolyLine(pnt[0], pnt[1]);
	marrPlnSegm[1] =  TURPolyLine(pnt[1], pnt[2]);
	marrPlnSegm[2] =  TURPolyLine(pnt[1], pnt[3]);
 }


 void   TEmitterImage::createShapeFiles(wchar_t *wchFoldName)
{
   	wchar_t wchFileName[300]= {0};
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\Emitter.shp");
	TURPolyLine::WriteSetSHPFiles(wchFileName,marrPlnSegm, 3) ;
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\EmitPntCentre.shp");
	mpntCentre.WriteSetSHPFiles(wchFileName,&mpntCentre, 1) ;
}

// сдвиг эмиттера на вектор  pntSdvig
TEmitterImage   TEmitterImage::SdvigEmTransform( TEmitterImage ImageInp, TURPointXY pntSdvig )
{
  ImageInp.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\This");;
  TEmitterImage EmitRez ;
  EmitRez.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\EmitRez");


	   EmitRez.marrPlnSegm[0].Points[0].X =   ImageInp.marrPlnSegm[0].Points[0].X + pntSdvig.X;
	   ImageInp.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\This");
	   EmitRez.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\EmitRez");

	   EmitRez.marrPlnSegm[0].Points[0].Y =   ImageInp.marrPlnSegm[0].Points[0].Y + pntSdvig.Y;
	   ImageInp.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\This");
	   EmitRez.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\EmitRez");
  for (int i = 0; i < 3; i++)
  {
	 for (int j = 0; j < 2; j++)
	 {
	   ImageInp.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\This");
	   EmitRez.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\EmitRez");

	   EmitRez.marrPlnSegm[i].Points[j].X =   ImageInp.marrPlnSegm[i].Points[j].X + pntSdvig.X;
	   ImageInp.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\This");
	   EmitRez.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\EmitRez");

	   EmitRez.marrPlnSegm[i].Points[j].Y =   ImageInp.marrPlnSegm[i].Points[j].Y + pntSdvig.Y;
	   ImageInp.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\This");
	   EmitRez.createShapeFiles (   L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\EmitRez");

	 }
  }
  EmitRez.mpntCentre = TURPointXY(ImageInp.mpntCentre.X +  pntSdvig.X, ImageInp.mpntCentre.Y +  pntSdvig.Y);
  return EmitRez;
}

#pragma package(smart_init)
