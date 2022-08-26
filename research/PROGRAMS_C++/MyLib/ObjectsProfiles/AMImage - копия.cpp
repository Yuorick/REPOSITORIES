//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include "AMImage.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "EmitterImage.h"

 //---------------------------------------------------------------------------

__fastcall TAMImage::~TAMImage()
{
	if(mparrEmit) free(mparrEmit); ;
	mparrEmit = NULL ;

}
/*
__fastcall TAMImage::TAMImage()
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


}  */
// Конструктор копирования
 TAMImage::TAMImage (const TAMImage &R)
 {
	mquantEmit  = R.mquantEmit ;
	if (mparrEmit != NULL)
	{
	free(mparrEmit);
	mparrEmit = NULL;
	}
	if(R.mparrEmit != NULL)
	{
	mparrEmit = (TEmitterImage*)malloc( mquantEmit * sizeof(mquantEmit));
	if(mparrEmit == NULL) 	   {
	ShowMessage(L"Not memory for mparrEmit") ;
	Abort() ;
	}
	memcpy( mparrEmit,R.mparrEmit, R.mquantEmit  * sizeof(TEmitterImage));
	}

	mpntCentre = R.mpntCentre;
	mplgHood = R.mplgHood;;
	mplnOut = R.mplnOut;

 }
// оператор присваивания
 TAMImage TAMImage::operator=(TAMImage  R)
{
	mquantEmit  = R.mquantEmit ;
	if (mparrEmit != NULL)
	{
	free(mparrEmit);
	mparrEmit = NULL;
	}
	if(R.mparrEmit != NULL)
	{
	mparrEmit = (TEmitterImage*)malloc( mquantEmit * sizeof(mquantEmit));
	if(mparrEmit == NULL) 	   {
	ShowMessage(L"Not memory for mparrEmit") ;
	Abort() ;
	}
	memcpy( mparrEmit,R.mparrEmit, R.mquantEmit  * sizeof(TEmitterImage));
	}

	mpntCentre = R.mpntCentre;
	mplgHood = R.mplgHood;;
	mplnOut = R.mplnOut;
  return *this ;
}

// парам констр
// Emit - образец излучателя
// QuantEmit  - к-во излучателей
// DistEmit  - расстояние между издучателями
// SideWidth - толщина верхней и нижней частей кожуха
// BackWidth - толщина задней части кожуха
// LenghPlnOut  - длина линии выхода
 __fastcall TAMImage::TAMImage(TEmitterImage Emit ,const int QuantEmit,const double DistEmit,
 const double SideWidth, const double BackWidth , const double LenghPlnOut )

 {
	mpntCentre = TURPointXY(0.,0.);
	double valLen = ((double)(QuantEmit -1)) *  DistEmit/2.;
	mparrEmit =(TEmitterImage*) malloc( QuantEmit * sizeof( TEmitterImage));
	mquantEmit = QuantEmit;
	for (int i = 0; i < mquantEmit; i++)
	{
	  TURPointXY pntSdvig(0.,  valLen -  ((double)i) *DistEmit);
 //	  Emit.createShapeFiles( L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\FARPic2t") ;
	  mparrEmit[i] =  Emit.SdvigEmTransform(pntSdvig);
//	  mparrEmit[i].createShapeFiles( L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\FARPicti") ;
  //  mparrEmit[0].createShapeFiles( L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\FARPicti") ;
	  int iii = 0;
	}
//	mparrEmit[1].createShapeFiles( L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\FARPict1") ;
 //	mparrEmit[0].createShapeFiles( L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\FARPict0") ;
	mplgHood = TURPolygon ( 9) ;
	double valLenEmit = Emit.marrPlnSegm[0].Points[1].X- Emit.marrPlnSegm[0].Points[0].X ;
	mplgHood.Points[0] =  TURPointXY (Emit.marrPlnSegm[0].Points[1].X- 0.25 * valLenEmit,valLen +DistEmit/2. );
	mplgHood.Points[8] = mplgHood.Points[0] ;
	mplgHood.Points[3] = TURPointXY( mplgHood.Points[0].X ,- mplgHood.Points[0].Y);
	mplgHood.Points[1] =  TURPointXY (Emit.marrPlnSegm[0].Points[0].X,mplgHood.Points[0].Y);
	mplgHood.Points[2] = TURPointXY( mplgHood.Points[1].X, -mplgHood.Points[1].Y);
	mplgHood.Points[7] = TURPointXY (mplgHood.Points[0].X, mplgHood.Points[0].Y + SideWidth);
	mplgHood.Points[4] = TURPointXY (mplgHood.Points[7].X,-mplgHood.Points[7].Y);
	mplgHood.Points[6] = TURPointXY (mplgHood.Points[1].X -BackWidth   ,mplgHood.Points[7].Y);
	mplgHood.Points[5] = TURPointXY (mplgHood.Points[6].X,-mplgHood.Points[6].Y);
  //	TURPointXY::WriteSetSHPFiles(L"D:\\PROJECTS_C++\\НЛЦ_сборка3\\New\\FARPict\\poinys.shp",mplgHood.Points, 9) ;


	TURPointXY  pnt1(mplgHood.Points[5].X, 0.);
	TURPointXY  pnt2 (pnt1.X - LenghPlnOut, 0.);
	mplnOut = TURPolyLine(   pnt1,  pnt2) ;
 }


 void  __fastcall TAMImage::createShapeFiles(wchar_t *wchFoldName)
{
// 1.
   	wchar_t wchFileName[300]= {0};
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\Emitters.shp");
	TURPolyLine *ppln  = (TURPolyLine*)malloc(3 * mquantEmit * sizeof(TURPolyLine));
	for (int i =0; i < mquantEmit; i++)
	{
	memcpy(&ppln[3 * i], mparrEmit[i].marrPlnSegm, 3 * sizeof(TURPolyLine));
	}
	TURPolyLine::WriteSetSHPFiles(wchFileName,ppln, 3 * mquantEmit) ;
	///
// 2.
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\АМPntCentre.shp");
	mpntCentre.WriteSetSHPFiles(wchFileName,&mpntCentre, 1) ;
///
// 3.
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\HoodAMplg.shp");
	mplgHood.WriteSetSHPFiles(wchFileName,&mplgHood, 1) ;
	///
// 4.
    wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\OutLine.shp");
	mplnOut.WriteSetSHPFiles(wchFileName,&mplnOut, 1) ;

	free(ppln);
}


#pragma package(smart_init)
