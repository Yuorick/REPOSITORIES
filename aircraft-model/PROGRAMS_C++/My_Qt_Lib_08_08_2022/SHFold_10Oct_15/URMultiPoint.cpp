//---------------------------------------------------------------------------



#include <math.h>
#include <stdlib.h>
#include <wchar.h>
#include "URMultiPoint.h"
#include <stdio.h>
#include "URPolygon.h"

 TURMultiPoint::~TURMultiPoint()
{


	if(Points) //delete [] Points ;
	Points = NULL ;
  //	if(iPoints) delete iPoints ;
   //	iPoints = NULL ;
}
//---------------------------------------------------------------------------

 // оператор присваивания
 TURMultiPoint &TURMultiPoint::operator=(const TURMultiPoint &R)
 {
	NumPoints = R.NumPoints ;
	if (Points != NULL) Points = NULL ;//delete [] Points;
	if(R.Points != NULL)
	{
		Points = new TURPointXY[R.NumPoints];
		if(Points == NULL)
		{
        //ShowMessage(L"Not memory for Points") ;
		//Abort() ;
		}
		memcpy( Points,R.Points, R.NumPoints  * sizeof(TURPointXY));
	}
	 memcpy(Box,R.Box,4 * sizeof(double));
	return *this ;
 }

 // конструктор копирования
 TURMultiPoint::TURMultiPoint (const TURMultiPoint &R)
 {
	NumPoints = R.NumPoints ;
	if (Points != NULL) Points = NULL ;//delete [] Points;
	if(R.Points != NULL)
	{
		Points = new TURPointXY[R.NumPoints];
		if(Points == NULL)
		{
        //ShowMessage(L"Not memory for Points") ;
		//Abort() ;
		}
		memcpy( Points,R.Points, R.NumPoints  * sizeof(TURPointXY));
	}
	memcpy(Box,R.Box,4 * sizeof(double));

 }
//-----------------------------------------------------------------------------------
// парам констр
//--------------------------------------------------------------------------------------
 TURMultiPoint ::TURMultiPoint()
{
	NumPoints = 0 ;
	Points = NULL ;
}

 // парам констр
TURMultiPoint :: TURMultiPoint( TURPointXY *ppntArr, const int iNumPoints)

 {
	NumPoints = iNumPoints ;
	Points = NULL ;
	Points = new TURPointXY[iNumPoints];
	if(Points == NULL)
	{
    //ShowMessage(L"Not memory for Points") ;

	}
	memcpy(Points, ppntArr, iNumPoints * sizeof(TURPointXY));
 }


void TURMultiPoint::calcBoundBox()
{
	 double xMin = 1000000000 ;
	 double xMax = -1000000000 ;
	 double yMin =  1000000000 ;
	 double yMax =  - 1000000000 ;
	 for (int i =0; i < NumPoints; i++)
	 {
	   if (Points[i].X > xMax) xMax =  Points[i].X;
	   if (Points[i].X < xMin) xMin =  Points[i].X;
	   if (Points[i].Y > yMax) yMax =  Points[i].Y;
	   if (Points[i].Y < yMin) yMin =  Points[i].Y;

	 }
	Box[0] =  xMin;
	Box[1] =  yMin;
	Box[2] =  xMax ;
	Box[3] =  yMax;

}


void TURMultiPoint::createTargPointsArray(const int valTargCellSize, TURPointXY **ppTargPntArray, int *lenTargPntArray)
{
  *ppTargPntArray = (TURPointXY*)realloc( *ppTargPntArray, NumPoints* sizeof(TURPointXY));
  memcpy(*ppTargPntArray, Points, NumPoints * sizeof(TURPointXY));
  *lenTargPntArray = NumPoints;
}

	//http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
	// http://www.dbase.com/Knowledgebase/INT/db7_file_fmt.htm
	//http://www.clicketyclick.dk/databases/xbase/format/dbf.html#DBF_STRUCT

// изменение порядка следования байтов в 4 байтном слове (32 бита массив)
//  input : chstr - указатель на массив char[4]
// output: chstr - указатель на массив  char[4] c измененным порядком следования
// байтов  chstr[0] = chstr[3] ; chstr1] = chstr[2] ;  chstr[2] = chstr[1] ; chstr[3] = chstr[0] ;
void TURMultiPoint::ChangeByteOrder(int * pi0)
{
  char c;
  char * chstr = (char*)pi0;
  c = chstr[0];
  chstr[0] = chstr[3] ;
  chstr[3]  = c;
  c = chstr[1];
 chstr[1] = chstr[2] ;
  chstr[2]  = c;

}


void TURMultiPoint::WriteSetSHPFiles(wchar_t *wchFileName)
{

    int lenFileName = wcslen( wchFileName);

	if (!((wchFileName[lenFileName -1] == L'p')
	   && (wchFileName[lenFileName -2] == L'h')
	   && (wchFileName[lenFileName -3] == L's')))
	   {
         //ShowMessage(L" Error file name") ;
		 return;
	   }
	 wchar_t wchDBFFileName [200] ,wchSHXFileName [200] ;
	wcscpy(wchSHXFileName, wchFileName);
	wchSHXFileName[lenFileName -1] =  L'x' ;
	wcscpy(wchDBFFileName, wchFileName);
	wchDBFFileName [lenFileName -1] = L'f';
	wchDBFFileName [lenFileName -2] = L'b';
	wchDBFFileName [lenFileName -3] = L'd';
	WriteDBASEFile(wchDBFFileName) ;
	WriteMainFile(wchFileName) ;
	WriteIndexFile(wchSHXFileName) ;
}



void TURMultiPoint::WriteDBASEFile(wchar_t *wchFileName)
{
	 FILE  *fw ;
	 fw=_wfopen(wchFileName,L"wb");
    // if(!fw) ////ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;
	  //1. version number and Date of last update  - первые 0-3  байта
	  // char ch = '\x03';
	 int i0 = 201421059;
		  fwrite(&i0,sizeof(int),1 ,fw) ;
		// 2. Number of records  in data file   4-7 байты
		int numRecords = NumPoints ;
		   fwrite(&numRecords,sizeof(int),1 ,fw) ;
		// 3. bytes 8 -31
		int iarr[] = {458817,0,0,0,0,22272};
		  fwrite(iarr,sizeof(int),6 ,fw) ;
		// 4. описание полей bytes 32-63
		int iarr1[] = {25673
					   ,0
					   ,1308622848
					   ,0
					   ,6
					   ,0
					   ,0
					   ,0
					   };
		 fwrite(iarr1,sizeof(int),8 ,fw) ;
	  // 5.  Terminator (0Dh) byte 64
	  char ch = '\x0D';
		 fwrite(&ch,1,1 ,fw) ;
	   //
	 char arrch[7] = {'\x20','\x20','\x20','\x20','\x20','\x20','\x30'
	   };
	   for (int i =0 ; i < NumPoints; i++)  fwrite(arrch,1,7 ,fw) ;


	 fclose(fw);
}

void TURMultiPoint::WriteMainFile(wchar_t *wchFileName)
{
   const int SHAPE_TYPE = 1;
  FILE  *fw ;

	 fw=_wfopen(wchFileName,L"wb");
   //	 if(!fw) ////ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;

	 // Main File HEADER
	 //*************************************************************************************
	 //1. step.  byte 0 - 23
	 // its are 6 integer values, getting by copying real shp file
	 int iarr[6] = {170328064
					,0
					,0
					,0
					,0
					,0 };
		 fwrite(iarr,sizeof(int),6 ,fw) ;
	   ///

	   // 2. The value for file length is the total length of the file in 16-bit words (including the fifty
		//16-bit words that make up the header).   Byte4s 24-27
		int iFileLeng = (100 + 28 * NumPoints) /2;


		ChangeByteOrder(& iFileLeng);
		 fwrite(&iFileLeng,sizeof(int),1 ,fw) ;
		///

		// 3. version = 1000 shapetype =5  bytes 28-35
		   iarr[0] =  1000;
		   iarr[1] =  SHAPE_TYPE ;
		   fwrite(iarr,sizeof(int),2 ,fw) ;
		   ///

		// 4. Bounding box
		double Box[8] = {0};
		Box[0] = Points[0].X ; Box[1] = Points[0].Y ;
		Box[2] = Points[0].X ; Box[3] = Points[0].Y ;

		for (int i =0 ; i < NumPoints; i++)
		{

		   if (Points[i].X <  Box[0]) Box[0] =  Points[i].X;
		   if (Points[i].Y <  Box[1]) Box[1] =  Points[i].Y;
		   if (Points[i].X >  Box[2]) Box[2] =  Points[i].X;
		   if (Points[i].Y >  Box[3]) Box[3] =  Points[i].Y;
		}
		  fwrite(Box,sizeof(double),8 ,fw) ;
		///******************************************************************************************
		///******************************************************************************************

	   // RECORDS
	   for (int i = 0; i < NumPoints; i++)
	   {
		  // 1. Record Headers
		   // Byte 0 Record Number Record Number Integer Big
		   int iRecNum = i;
		   ChangeByteOrder(& iRecNum);
		   //Byte 4 Content Length Content Length Integer Big
		   int iContLeng = 20/2 ;
		   ChangeByteOrder(& iContLeng);
		   iarr[0] =  iRecNum;
		   iarr[1] =  iContLeng ;
			 fwrite(iarr,sizeof(int),2 ,fw) ;
		   //*****************************************************
		   //********** Content *******************************************
				int ishapetype =   SHAPE_TYPE ;
			  fwrite(&ishapetype ,sizeof(int),1 ,fw) ;
			  fwrite(&(Points[i].X),sizeof(double),1 ,fw) ;
			  fwrite(&(Points[i].Y),sizeof(double),1 ,fw) ;

	   }

  fclose(fw);
}

void TURMultiPoint::WriteIndexFile(wchar_t *wchFileName)
{
  FILE  *fw ;
	const int SHAPE_TYPE = 1;
	 fw=_wfopen(wchFileName,L"wb");
   //	 if(!fw) ////ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;
	  // Ibdex(SHX) File HEADER
	 //*************************************************************************************
	 //1. step.  byte 0 - 23
	 // its are 6 integer values, getting by copying real shp file
	 int iarr[6] = {170328064
					,0
					,0
					,0
					,0
					,0 };
		fwrite(iarr,sizeof(int),6 ,fw) ;
	   ///

	   // 2. The value for file length is the total length of the file in 16-bit words (including the fifty
		//16-bit words that make up the header).   Byte4s 24-27
		int iFileLeng = ( 100 + NumPoints *8 )/2 ;;
		ChangeByteOrder(& iFileLeng);
		fwrite(&iFileLeng,sizeof(int),1 ,fw) ;
		///

		// 3. version = 1000 shapetype =5  bytes 28-35
		   iarr[0] =  1000;
		   iarr[1] =  SHAPE_TYPE ;
		   fwrite(iarr,sizeof(int),2 ,fw) ;
		   ///

		// 4. Bounding box
		// 4. Bounding box
		double Box[8] = {0};
		Box[0] = Points[0].X ; Box[1] = Points[0].Y ;
		Box[2] = Points[0].X ; Box[3] = Points[0].Y ;

		for (int i =0 ; i < NumPoints; i++)
		{

		   if (Points[i].X <  Box[0]) Box[0] =  Points[i].X;
		   if (Points[i].Y <  Box[1]) Box[1] =  Points[i].Y;
		   if (Points[i].X >  Box[2]) Box[2] =  Points[i].X;
		   if (Points[i].Y >  Box[3]) Box[3] =  Points[i].Y;
		}
		  fwrite(Box,sizeof(double),8 ,fw) ;
		///******************************************************************************************
		///******************************************************************************************
	   // RECORDS
	   int offset = 100;
	   for (int i = 0; i < NumPoints; i++)
	   {
		  int ioffset = ( offset + i * 28 )/2;
		  ChangeByteOrder(& ioffset);
		   int iContLeng = 20/2 ;
		   ChangeByteOrder(& iContLeng);
			 fwrite(&ioffset,sizeof(int),1 ,fw) ;
			 fwrite(&iContLeng,sizeof(int),1 ,fw) ;

	   }

	   fclose(fw);

}

 /*******************************************************************************************************************************/

//void TURMultiPoint::createAimPointsArray(const int valAimCellSize, TURPointXY **ppAimPntArray, int *lenAimPntArray)
//{

//}
void TURMultiPoint::createUnatedPointsArray(TURPointXY **ppunatedPoints, int *quantUnitedPoints)
{
 *ppunatedPoints = (TURPointXY *)realloc(*ppunatedPoints, NumPoints * sizeof( TURPointXY));
 *quantUnitedPoints = NumPoints;
  memcpy(*ppunatedPoints, Points, NumPoints * sizeof( TURPointXY));

}

 /*******************************************************************************************************************************/
// линейное преобразование полигона  пергруженная
// INPUT:
// valAng - угол поворота
// pntSdvig - точка куда перемещается центр полигона
// valRastigenie - коэффициент растяжения
// OUTPUT:
// возвращает преобразованный полигон
void   TURMultiPoint::LinearTransformation(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie )
{
  double arrMtxPer[4] = {0.};
	arrMtxPer[0] = cos(valAng);
  arrMtxPer[1] = -sin(valAng);
  arrMtxPer[2] = -arrMtxPer[1];
	arrMtxPer[ 3] = arrMtxPer[0];

	for (int i =0; i < NumPoints; i++)
	{
		Points[i] =  Points[i].LinTransform( valAng ,  pntSdvig, valRastigenie );
	}

}

//------------------------------------------------------------------------------
// нахождение выпуклой оболочки

void   TURMultiPoint::ConvexShell(TURPolygon *pPolgConv)
{
	*pPolgConv =   TURPolygon::Conv(Points // массив точек, input
			, NumPoints // длина массива точек , input
				)  ;
}

double  TURMultiPoint::calcDiam(int *pnum0,  int *pnum1)
{
	return TURPointXY::calcDiam(Points, NumPoints, pnum0,  pnum1);
}



