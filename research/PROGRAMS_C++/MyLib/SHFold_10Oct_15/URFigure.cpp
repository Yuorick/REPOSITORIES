//---------------------------------------------------------------------------
  #pragma package(smart_init)
#include <vcl.h>

#pragma hdrstop
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "MatrixProccess.h"
#include <stdlib.h>
#include "URFigure.h"
#include "UrPointXY.h"
#include "URPolygon.h"
#include "URPolyLine.h"

//---------------------------------------------------------------------------



__fastcall TURFigure::TURFigure()
{
	/*RecNumber = 0 ;
	RecLength = 0 ;
	Type = ShapeType::NullShape ;
	FigureColor = (TColor)0 ;*/
}
//---------------------------------------------------------------------------

//__fastcall TURFigure::~TURFigure()
//{
//}

void TURFigure::createTargPointsArray(const int valTargCellSize, TURPointXY **ppTargPntArray, int *lenTargPntArray)
{

}

// формирование массива точек прицеливания
// INPUT:
// valTargCellSize - шаг ячейки (расстояния между точками)
// OUTPUT:
// *ppAimPntArray  - массив точек
// *lenAimPntArray  - его длина
void TURFigure::createAimPointsArray(const int valTargCellSize, TURPointXY **ppAimPntArray, int *lenAimPntArray)
{
int quantUnitedPoints = 1;
  TURPointXY *punatedPoints = (TURPointXY*)malloc(quantUnitedPoints * sizeof (TURPointXY));
  TURPointXY **ppunatedPoints = & punatedPoints;
  createUnatedPointsArray(ppunatedPoints,&quantUnitedPoints);
  ///
   TURPolygon PlgAim;
   TURPointXY arrPntTemp[3];

   switch(quantUnitedPoints)
   {
	   case 1:
	   arrPntTemp [0] = (*ppunatedPoints)[0];
	   arrPntTemp [1] = (*ppunatedPoints)[0];
	   arrPntTemp [2] = (*ppunatedPoints)[0];
	   PlgAim = TURPolygon(3, arrPntTemp);
	   break;

	   case 2:
	   arrPntTemp [0] = (*ppunatedPoints)[0];
	   arrPntTemp [1] = (*ppunatedPoints)[1];
	   arrPntTemp [2] = (*ppunatedPoints)[0];
	   PlgAim = TURPolygon(3, arrPntTemp);
	   break;

	   default:
	   PlgAim = 	TURPolygon::Conv(*ppunatedPoints // массив точек, input
			, quantUnitedPoints// длина массива точек , input
				);
	   break;


   }
	 free(punatedPoints);


	if (fabs(PlgAim.calcVectSq())< 400.)  //
	{               // это линия или точка
		int num0 = -1,  num1 = -1;
		double valDiam = TURPointXY::calcDiam(PlgAim.Points, PlgAim.NumPoints, &num0,  &num1);
		double valAngRotate = atan2(PlgAim.Points[num1].Y - PlgAim.Points[num0].Y, PlgAim.Points[num1].X - PlgAim.Points[num0].X);
		double  arrMtxPer[4] = {0.};
		arrMtxPer[0] = cos (valAngRotate);
		arrMtxPer[1] = -sin (valAngRotate);
		arrMtxPer[2] = -arrMtxPer[1];
		arrMtxPer[3] = arrMtxPer[0];
		TURPolygon PlgAimRotated = PlgAim.fncLinTransform(arrMtxPer )  ;
		PlgAimRotated.calcBoundBox() ;
		double temp = fabs(PlgAimRotated.Box[2] - PlgAimRotated.Box[0]);
		//( fabs(PlgAimRotated.Box[1]) > fabs(PlgAimRotated.Box[3]))? fabs(PlgAimRotated.Box[1]): fabs(PlgAimRotated.Box[3]);

		if (valDiam < 40.)
		{   // это точка
			*ppAimPntArray = (TURPointXY *)realloc(*ppAimPntArray, sizeof(TURPointXY));
			*lenAimPntArray = 1;
			(*ppAimPntArray)[0].X =  (PlgAim.Points[num1].X + PlgAim.Points[num0].X) /2.;
			(*ppAimPntArray)[0].Y =  (PlgAim.Points[num1].Y + PlgAim.Points[num0].Y) /2.  ;
			return;
		}
		else
		{ // это отрезок
			TURPolyLine  plnAim( PlgAim.Points[num0], PlgAim.Points[num1]) ;
			plnAim.createTargPointsArray(valTargCellSize, ppAimPntArray, lenAimPntArray);
			return;
		}

	}
	else
	{  // это полигон
		PlgAim.createTargPointsArray(valTargCellSize, ppAimPntArray, lenAimPntArray);
		return;
	}

}

void TURFigure::createUnatedPointsArray(TURPointXY **ppunatedPoints, int *quantUnitedPoints)
{

}

void TURFigure::find_Objects_Type_And_Quant(wchar_t *wchFileName, int *ipShapeType, int*ipQuant)
{
	 int lenFileName = wcslen( wchFileName);

	if (!((wchFileName[lenFileName -1] == L'p')
	   && (wchFileName[lenFileName -2] == L'h')
	   && (wchFileName[lenFileName -3] == L's')))
	   {
		 ShowMessage(L" Error file name") ;
		 return;
	   }
	wchar_t  wchSHXFileName[200] ;
	wcscpy(wchSHXFileName, wchFileName);
	wchSHXFileName[lenFileName -1] =  L'x' ;
	FILE  *fr0 ;

	 fr0=_wfopen(wchSHXFileName,L"rb");
	 if(!fr0) ShowMessage (L"TURPointXY::ReadSHPFile\nFile is not opened !") ;
	  int offset0 = 32;
	  int ishapetype = -1;

		fseek(fr0,offset0,SEEK_SET);
		fread(ipShapeType ,sizeof(int), 1,fr0) ;

	   offset0 = 24;
	  fseek(fr0,offset0,SEEK_SET);
	  int lenSHXFile = -1;
	  fread(&lenSHXFile ,sizeof(int), 1,fr0) ;
	  ChangeByteOrder( &lenSHXFile);
		*ipQuant =( 2 *  lenSHXFile - 100 )/ 8 ;
		 fclose(fr0);
}
//-----------------------------------------------------------------------------------
// Проверка того, что фигура является отрезком , точкой или полигоном
//
//
int  TURFigure::calcDimension(const double VAlTolerance)
{
	int quantUnitedPoints = 1;
	TURPointXY *punatedPoints = (TURPointXY*)malloc(quantUnitedPoints * sizeof (TURPointXY));
	TURPointXY **ppunatedPoints = & punatedPoints;
	createUnatedPointsArray(ppunatedPoints,&quantUnitedPoints);
	///
   TURPolygon PlgAim;
   TURPointXY arrPntTemp[3];

   switch(quantUnitedPoints)
   {
	   case 1:
	   arrPntTemp [0] = (*ppunatedPoints)[0];
	   arrPntTemp [1] = (*ppunatedPoints)[0];
	   arrPntTemp [2] = (*ppunatedPoints)[0];
	   PlgAim = TURPolygon(3, arrPntTemp);
	   break;

	   case 2:
	   arrPntTemp [0] = (*ppunatedPoints)[0];
	   arrPntTemp [1] = (*ppunatedPoints)[1];
	   arrPntTemp [2] = (*ppunatedPoints)[0];
	   PlgAim = TURPolygon(3, arrPntTemp);
	   break;

	   default:
	   PlgAim = 	TURPolygon::Conv(*ppunatedPoints // массив точек, input
			, quantUnitedPoints// длина массива точек , input
				);
	   break;


   }
	 free(punatedPoints);
	 if(fabs(PlgAim.calcVectSq())> VAlTolerance * VAlTolerance)
	 {
	 return 2;
	 }

	 int num0 = -1,  num1 = -1;
	 double valDiam = TURPointXY::calcDiam(PlgAim.Points, PlgAim.NumPoints, &num0,  &num1);
		if (valDiam < VAlTolerance * 2.)
		{
			return 0;
		}
	 /*	double valAngRotate = atan2(PlgAim.Points[num1].Y - PlgAim.Points[num0].Y, PlgAim.Points[num1].X - PlgAim.Points[num0].X);
		double  arrMtxPer[4] = {0.};
		arrMtxPer[0] = cos (valAngRotate);
		arrMtxPer[1] = -sin (valAngRotate);
		arrMtxPer[2] = -arrMtxPer[1];
		arrMtxPer[3] = arrMtxPer[0];
		TURPolygon PlgAimRotated = PlgAim.fncLinTransform(arrMtxPer )  ;
		PlgAimRotated.calcBoundBox() ;
		double temp = fabs(PlgAimRotated.Box[2] - PlgAimRotated.Box[0]);
		if (temp < VAlTolerance)
		{
			return 0;
		}  */

	return 1;

}

///
	//http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
	// http://www.dbase.com/Knowledgebase/INT/db7_file_fmt.htm
	//http://www.clicketyclick.dk/databases/xbase/format/dbf.html#DBF_STRUCT

// изменение порядка следования байтов в 4 байтном слове (32 бита массив)
//  input : chstr - указатель на массив char[4]
// output: chstr - указатель на массив  char[4] c измененным порядком следования
// байтов  chstr[0] = chstr[3] ; chstr1] = chstr[2] ;  chstr[2] = chstr[1] ; chstr[3] = chstr[0] ;
void TURFigure::ChangeByteOrder(int * pi0)
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

void TURFigure::calcBoundBox()
{
}


void TURFigure::WriteSetSHPFiles(wchar_t *wchFileName)
{

}

void   TURFigure::LinearTransformation(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie )
{

}

// нахождение выпуклой оболочки   полигона
void   TURFigure::ConvexShell(TURPolygon *pPolgConv)
{

}

// нахождение буфера выпуклой оболочки   массива точек
TURPolygon   TURFigure::Buffer(const double VAlBufferX,const double VAlBufferY)
{
	int quantUnitedPoints = 1;
	TURPointXY *punatedPoints = (TURPointXY*)malloc(quantUnitedPoints * sizeof (TURPointXY));
	TURPointXY **ppunatedPoints = & punatedPoints;
	createUnatedPointsArray(ppunatedPoints,&quantUnitedPoints);
	TURPolygon   plgConv0 = TURPolygon::Conv(*ppunatedPoints // массив точек, input
			, quantUnitedPoints // длина массива точек , input
				) ;
	free(punatedPoints);
	TURPointXY *pPntArr = new  TURPointXY[(plgConv0.NumPoints -1) *4];
	for (int i = 0; i < (plgConv0.NumPoints -1); i++)
	{
	 pPntArr[ i *4    ] =  TURPointXY(plgConv0.Points[i].X  -VAlBufferX, plgConv0.Points[i].Y  -VAlBufferY);
	 pPntArr[ i *4 + 1] =  TURPointXY(plgConv0.Points[i].X  -VAlBufferX, plgConv0.Points[i].Y  +VAlBufferY);
	 pPntArr[ i *4 + 2] =  TURPointXY(plgConv0.Points[i].X  +VAlBufferX, plgConv0.Points[i].Y  +VAlBufferY);
	 pPntArr[ i *4 + 3] =  TURPointXY(plgConv0.Points[i].X  +VAlBufferX, plgConv0.Points[i].Y  -VAlBufferY);
	}
	 int iii=0;
   //	 pPntArr[0].WriteSetSHPFiles(L"E:\\Ametist\\28-03-2018\\COAST\\Optimal_Strat_Rez\\pPntArr.shp", pPntArr,(plgConv0.NumPoints -1)* 4);
	TURPolygon   plgConv1 = TURPolygon::Conv(pPntArr // массив точек, input
			, (plgConv0.NumPoints -1) *4 // длина массива точек , input
				) ;

	delete []pPntArr;
	return plgConv1;

}

//-------------------------------------------------------------------------
// создание полигона буфера фигуры в соотвествии с корреляц матрицей
// вычисляются союбственные векторы и собственные числа заданной  корреляц матьрицы arrMtrxCorrSyst[4]
// точки фигуры переводяьтся в сиситему координаит осей эллипса корреляц матьрицы arrMtrxCorrSyst[4]
// - то есть в сиситему координат собственных векторов матрицы
// для каждой точки фигуры, перведенной в сиситему координат элиипса, строятся 4 точки, окружааающие ее
// на расстояниях sqrt(lambdax) * VAlCoeff по оси X и   sqrt(lambday) * VAlCoeff по оси Y
// далее, строится полигон-выпуклая оболочка этого массива точек и, наконец, переводится в исходную сиситему координат
TURPolygon   TURFigure::createBuffPolygonInAccordanceWithCorrMtrx(double *arrMtrxCorrSyst, const double VAlCoeff)
{
	TURPolygon plgConvGSK0;
	ConvexShell(&plgConvGSK0);
	 //	plgConvGSK0.WriteSetSHPFiles(L"E:\\Ametist\\28-03-2018\\COAST\\Optimal_Strat_Rez\\plgConvGSK0.shp", &plgConvGSK0, 1);
		///
		double arrF_Syst[4] = {0.} ,arrF_SystTr[4] = {0.} , arrMtrxLamb_Syst[4] = {0.}, arrPos00[2] = {0.};
		CalcProperVectors2(arrMtrxCorrSyst, arrF_Syst, arrMtrxLamb_Syst) ;
		MatrTransp(arrF_Syst, 2, 2, arrF_SystTr);
		TURPolygon plgConvGSK01=  plgConvGSK0.fncLinTransform(arrF_SystTr);
	 //	plgConvGSK01.WriteSetSHPFiles(L"E:\\Ametist\\28-03-2018\\COAST\\Optimal_Strat_Rez\\plgConvGSK01.shp", &plgConvGSK01, 1);
		///
		TURPolygon plgConvGSK02 = plgConvGSK01.Buffer(VAlCoeff * sqrt(arrMtrxLamb_Syst[0]), VAlCoeff * sqrt(arrMtrxLamb_Syst[3])) ;
		TURPolygon plgAim = plgConvGSK02.fncLinTransform(arrF_Syst);
		return plgAim;
}

/*
bool __fastcall TURFigure::LoadFromStream(TStream *Stream)
{
	if(4 != Stream->Read(&RecNumber, 4)) return false ;
	if(4 != Stream->Read(&RecLength, 4)) return false ;
	if(4 != Stream->Read(&Type, 4)) return false ;
	SwapInt(RecNumber) ;
	SwapInt(RecLength) ;

return true ;
}
//---------------------------------------------------------------------------

bool __fastcall TURFigure::SaveToStream(TStream *Stream)
{
	int RNumber = RecNumber ;
	int RLength = RecLength ;

	SwapInt(RNumber) ;
	SwapInt(RLength) ;

	if(4 != Stream->Write(&RNumber, 4)) return false ;
	if(4 != Stream->Write(&RLength, 4)) return false ;
	if(4 != Stream->Write(&Type, 4)) return false ;

return true ;
}
//---------------------------------------------------------------------------

void __fastcall TURFigure::SwapInt(int &Data)
{
	unsigned char b ;
	b = *((unsigned char *)&Data+0) ;
	*((unsigned char *)&Data+0) = *((unsigned char *)&Data+3) ;
	*((unsigned char *)&Data+3) = b ;

	b = *((unsigned char *)&Data+1) ;
	*((unsigned char *)&Data+1) = *((unsigned char *)&Data+2) ;
	*((unsigned char *)&Data+2) = b ;
}
//---------------------------------------------------------------------------

void __fastcall TURFigure::Draw(TCanvas *Canvas, bool bSelect)
{

}
//---------------------------------------------------------------------------

TURFigure* __fastcall TURFigure::PtInFigure(const TPointXY &pt)
{
return NULL ;
}
//---------------------------------------------------------------------------

TURFigure* __fastcall TURFigure::PtNearVisibleFigure(const TPointXY &pt)
{
return NULL ;
}
//---------------------------------------------------------------------------

double __fastcall TURFigure::GetArea()
{
return -1 ;
}
//---------------------------------------------------------------------------

double __fastcall TURFigure::GetLength()
{
return -1 ;
}
//---------------------------------------------------------------------------

int __fastcall TURFigure::GetSize()
{
return -1 ;
}
//---------------------------------------------------------------------------

void __fastcall TURFigure::SetOrgExt(const TPointXY &Offset, double kExt)
{
}
//---------------------------------------------------------------------------
 */

 #pragma package(smart_init)
