//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include <dir.h >
#include "ShipTarg.h"
#include <string.h>
#include <math.h>
#include "MatrixProccess.h"
 #include "YrWriteShapeFile.h"
 #include "HomingHead3D.h"


const int QUANT_COLS_BUFF_ShipTarg = 7;

 TShipTarg::~TShipTarg()
{
   if (mparrBuff != NULL)
   {
   free( mparrBuff) ;
   mparrBuff = NULL ;
   }
   if ( mpwcharrFoldReport != NULL)
   {
   free( mpwcharrFoldReport );
   mpwcharrFoldReport = NULL ;
   }
}
//---------------------------------------------------------------------------


TShipTarg::TShipTarg()
{
  // memset(marrVectSostGSK, 0, 9 * sizeof(double)) ;
	mquantParts = 1;
	long double arrt[10]={0};

	TPartData Part_temp(0,0,0,arrt);
	for (int i = 0; i < 10; i++)
	{
	  marrPartData[i] = Part_temp;
	}
	marrPartData[0] = TPartData() ;
	mTargData = TTargData();

  mpwcharrFoldReport = NULL ;
  mparrBuff   = NULL ;
  mQuantPntReport = 0;
  mLenMemoryAlloc = 0 ;

}

//---------------------------------------------------------------------------


// конструктор копирования
 TShipTarg ::TShipTarg (const TShipTarg &R)
 {
   //	memcpy(marrVectSostGSK, R.marrVectSostGSK, 9 * sizeof(double)) ;
	mquantParts = R.mquantParts;
	memcpy(marrPartData, R.marrPartData, 10 * sizeof(TPartData)) ;
	mT = R.mT;
	mTraject = R.mTraject;
	mTargData = R.mTargData ;
	// для отчета
	mLenMemoryAlloc = R.mLenMemoryAlloc ;
	mQuantPntReport = R.mQuantPntReport ;

	////
   //	if (mparrBuff != NULL)
 //	{
   //		free(mparrBuff);
   //		mparrBuff = NULL;

  //  }
	mparrBuff  = NULL;
	mpwcharrFoldReport = NULL ;
	if(R.mparrBuff  != NULL)
	{
		if((mparrBuff  = ( double*)malloc( QUANT_COLS_BUFF_ShipTarg  * R.mLenMemoryAlloc * sizeof( double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_ShipTarg  * R.mLenMemoryAlloc * sizeof( double));
		}
		else
		{
		 ShowMessage(L"Not memory for mparrBuff ") ;
		 Abort() ;
		}
	}

	//
	if(R.mpwcharrFoldReport != NULL)
	{
		if (mpwcharrFoldReport != NULL)
		{
		  free(mpwcharrFoldReport );
		  mpwcharrFoldReport = NULL;

		}
		if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		wcscpy(mpwcharrFoldReport, R.mpwcharrFoldReport);

		}
		else
		{
		 ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		 Abort() ;
		}
	}
 }
 // оператор присваивания
 TShipTarg TShipTarg::operator=(TShipTarg  R)
 {
   //	memcpy(marrVectSostGSK, R.marrVectSostGSK, 9 * sizeof(double)) ;
	mquantParts = R.mquantParts;
	memcpy(marrPartData, R.marrPartData, 10 * sizeof(TPartData)) ;
	mT = R.mT;
	mTraject = R.mTraject;
   	mTargData = R.mTargData ;

	// для отчета
	mLenMemoryAlloc = R.mLenMemoryAlloc ;
	mQuantPntReport = R.mQuantPntReport ;

	mparrBuff  = NULL;
	mpwcharrFoldReport = NULL ;
	if(R.mparrBuff  != NULL)
	{
		if((mparrBuff  = (  double*)malloc( QUANT_COLS_BUFF_ShipTarg  * R.mLenMemoryAlloc * sizeof( double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_ShipTarg  * R.mLenMemoryAlloc * sizeof( double));
		}
		else
		{
		 ShowMessage(L"Not memory for mparrBuff ") ;
		 Abort() ;
		}
	}

	//
	if(R.mpwcharrFoldReport != NULL)
	{
		if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		wcscpy(mpwcharrFoldReport, R.mpwcharrFoldReport);

		}
		else
		{
		 ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		 Abort() ;
		}
	}  //	memcpy(marrVectSostGSK, R.marrVectSostGSK, 9 * sizeof(double)) ;
	mquantParts = R.mquantParts;
	memcpy(marrPartData, R.marrPartData, 10 * sizeof(TPartData)) ;
	mT = R.mT;
	mTraject = R.mTraject;
	// для отчета
	mLenMemoryAlloc = R.mLenMemoryAlloc ;
	mQuantPntReport = R.mQuantPntReport ;

	mparrBuff  = NULL;
	mpwcharrFoldReport = NULL ;
	if(R.mparrBuff  != NULL)
	{
		if((mparrBuff  = ( double*)malloc( QUANT_COLS_BUFF_ShipTarg  * R.mLenMemoryAlloc * sizeof( double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_ShipTarg  * R.mLenMemoryAlloc * sizeof( double));
		}
		else
		{
		 ShowMessage(L"Not memory for mparrBuff ") ;
		 Abort() ;
		}
	}

	//
	if(R.mpwcharrFoldReport != NULL)
	{
		if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		wcscpy(mpwcharrFoldReport, R.mpwcharrFoldReport);

		}
		else
		{
		 ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		 Abort() ;
		}
	}

	return *this ;
 }

  // парам конструктор
  // Bearing -  угол пеленга цели в ГСК
  // TargCourse - угол курса цели в ГСК
  // TargZenitAng - угол между ветором скорости цели и напрвлением в зенит
  // V,R,H - скорость, дальность, высота
 TShipTarg::TShipTarg (const long double Bearing, const long double TargCourse
	, const long double TargZenitAng,  const long double V, const long double H
	,const long double R, const long double valT,const int quantParts
	,TPartData *arrPartData, wchar_t* pwcharrFoldReport)
 {
   long double valeps = asin(    H/R) ;


   long double  arrVectSostGSK[9] = {0} ;
   arrVectSostGSK[0] = R * cos(valeps) * sin(Bearing) ;
   arrVectSostGSK[1] = R * cos(valeps) * cos(Bearing) ;
   arrVectSostGSK[2] = H ;
   arrVectSostGSK[3] = V * sin(TargZenitAng) * sin(TargCourse) ;
   arrVectSostGSK[4] = V * sin(TargZenitAng) * cos(TargCourse) ;
   arrVectSostGSK[5] = V * cos(TargZenitAng) ;


	mquantParts = quantParts ;
	memcpy(marrPartData, arrPartData, quantParts * sizeof(TPartData)) ;
	mT = valT ;
	mTraject = TShipTraj(  arrPartData[0],arrVectSostGSK, valT) ;

  // дл яотчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 2000 ;
   if(pwcharrFoldReport != NULL)
   {
	   if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		wcscpy(mpwcharrFoldReport, pwcharrFoldReport);

		}
		else
		{
		 ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		 Abort() ;
		}

		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof( double))) != NULL)
		{

		memset( mparrBuff,  0,  QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof( double));
		}
		else
		{
		 ShowMessage(L"Not memory for mparrBuff ") ;
		 Abort() ;
		}

   }
	

 }
 /////
  // парам конструктор
  // Bearing -  угол пеленга цели в ГСК
  // TargCourse - угол курса цели в ГСК
  // TargZenitAng - угол между ветором скорости цели и напрвлением в зенит
  // V,R,H - скорость, дальность, высота
 TShipTarg::TShipTarg (const TInitTargData InitData

	,const int quantParts, TPartData *arrPartData, wchar_t* pwcharrFoldReport)
 {
   long double valeps = asin(InitData.mH /InitData.mR) ;


   long double  arrVectSostGSK[9] = {0} ;
   arrVectSostGSK[0] = InitData.mR * cos(valeps) * sin(InitData.mBearing  ) ;
   arrVectSostGSK[1] = InitData.mR * cos(valeps) * cos(InitData.mBearing ) ;
   arrVectSostGSK[2] = InitData.mH ;
   arrVectSostGSK[3] = InitData.mV * sin(InitData.mTargZenitAng) * sin(InitData.mTargCourse) ;
   arrVectSostGSK[4] = InitData.mV * sin(InitData.mTargZenitAng) * cos(InitData.mTargCourse) ;
   arrVectSostGSK[5] = InitData.mV * cos(InitData.mTargZenitAng) ;


	mquantParts = quantParts ;
	memcpy(marrPartData, arrPartData, quantParts * sizeof(TPartData)) ;
	mT = InitData.mT ;
	mTraject = TShipTraj(  arrPartData[0],arrVectSostGSK, mT) ;

  // дл яотчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;


	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 2000 ;
   if(pwcharrFoldReport != NULL)
   {
	   if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		wcscpy(mpwcharrFoldReport, pwcharrFoldReport);

		}
		else
		{
		 ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		 Abort() ;
		}

		if((mparrBuff  = ( double*)malloc( QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof( double))) != NULL)
		{

		memset( mparrBuff,  0,  QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof(double));
		}
		else
		{
		 ShowMessage(L"Not memory for mparrBuff ") ;
		 Abort() ;
		}

   }


 }
 /////
  // парам конструктор  3
  // Bearing -  угол пеленга цели в ГСК
  // TargCourse - угол курса цели в ГСК
  // TargZenitAng - угол между ветором скорости цели и напрвлением в зенит
  // V,R,H - скорость, дальность, высота
 TShipTarg::TShipTarg (const TInitTargData InitData

	,const int quantParts, TPartData *arrPartData,TTargData TargData, wchar_t* pwcharrFoldReport)
 {
   long double valeps = asin(InitData.mH /InitData.mR) ;


   long double  arrVectSostGSK[9] = {0} ;
   arrVectSostGSK[0] = InitData.mR * cos(valeps) * sin(InitData.mBearing  ) ;
   arrVectSostGSK[1] = InitData.mR * cos(valeps) * cos(InitData.mBearing ) ;
   arrVectSostGSK[2] = InitData.mH ;
   arrVectSostGSK[3] = InitData.mV * sin(InitData.mTargZenitAng) * sin(InitData.mTargCourse) ;
   arrVectSostGSK[4] = InitData.mV * sin(InitData.mTargZenitAng) * cos(InitData.mTargCourse) ;
   arrVectSostGSK[5] = InitData.mV * cos(InitData.mTargZenitAng) ;


	mquantParts = quantParts ;
	memcpy(marrPartData, arrPartData, quantParts * sizeof(TPartData)) ;
	mT = InitData.mT ;
	mTraject = TShipTraj(  arrPartData[0],arrVectSostGSK, mT) ;
	mTargData = TargData  ;
  // дл яотчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;


	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 2000 ;
   if(pwcharrFoldReport != NULL)
   {
	   if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		wcscpy(mpwcharrFoldReport, pwcharrFoldReport);

		}
		else
		{
		 ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		 Abort() ;
		}

		if((mparrBuff  = ( double*)malloc( QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof( double))) != NULL)
		{

		memset( mparrBuff,  0,  QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof(double));
		}
		else
		{
		 ShowMessage(L"Not memory for mparrBuff ") ;
		 Abort() ;
		}

   }


 }
 /////

 //
bool TShipTarg::recalcTrajPoint(const long double tNext)
{
	 int numNextPart = getNumCurrentPart(tNext);
	if(  getNumCurrentPart(mT)!= numNextPart)
	{
	   mTraject = TShipTraj (marrPartData[numNextPart],mTraject.marrVectSostGSK,mT);

	}

		bool breturn = mTraject.recalcTrajPoint(tNext ) ;
		mT =tNext;
		if(mpwcharrFoldReport)  updateReportData() ;
		return breturn ;



}

int TShipTarg:: getNumCurrentPart(const long double valt)
{
	long double tCurr = 0;
	int i =0 ;
	for (i =0 ; i < mquantParts; i++)
	{
	  tCurr += marrPartData[i].mTimePart ;
	  if ( valt <= (tCurr + 0.00001) ) return i;

	}
	return -1 ;
}

// вычисление типа маневра для момента времени valT
int TShipTarg::getNumMаneuvreType(const long double valt)
{
  int iNumPrt = getNumCurrentPart(valt);
  if (iNumPrt == -1) return -1 ;
  int iType =  marrPartData[iNumPrt].miTypePart;
  switch(iType)
  {
	  case 0:  // равном прямолин движение
	  return 0;
	  break ;
	  case 2: // движение по окружности
	  return 1 ;
	  break ;
	  case 4:  // самонаведение
	  return 2 ;
	  break ;
	  default:
	  return -1 ;
  }


}
long double TShipTarg:: getTimePartStarted(const long double valt)
{

	long double tCurr = 0;
	int i =0 ;
	for (i =0 ; i < mquantParts; i++)
	{
	  if ( valt <= (tCurr + marrPartData[i].mTimePart)  ) return tCurr;

	}
	return -1 ;
}
long double TShipTarg:: getTimePartWillFinish(const long double valt)
{
	long double tCurr = 0;
	int i =0 ;
	for (i =0 ; i < mquantParts; i++)
	{
	  if ( valt <= (tCurr + marrPartData[i].mTimePart)  ) return (tCurr + marrPartData[i].mTimePart  );

	}
	return -1 ;
}

 // занесение информации в массивы для отчета
 void TShipTarg::updateReportData()
 {
  try
  {
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = ( double *)realloc(mparrBuff, QUANT_COLS_BUFF_ShipTarg  * mLenMemoryAlloc * sizeof( double)) ;

	 }
	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF_ShipTarg ;

	   mparrBuff [num0]    =  mT ;

	   mparrBuff [num0 +1]  =  mTraject.marrVectSostGSK[0] ;

	   mparrBuff [num0 +2]  =  mTraject.marrVectSostGSK[1] ;

	   mparrBuff [num0 +3]  =  mTraject.marrVectSostGSK[2] ;

	   mparrBuff [num0 +4]  =  mTraject.marrVectSostGSK[3]  ;

	   mparrBuff [num0 +5]  =  mTraject.marrVectSostGSK[4] ;

	   mparrBuff [num0 +6]  =  mTraject.marrVectSostGSK[5]  ;



	   mQuantPntReport++ ;
	 }
	 catch(...)
	 {
	   //	 return;
     }

 }

// публикация  отчета
 void TShipTarg::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;
	 const int mColsBuff =  QUANT_COLS_BUFF_ShipTarg ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, mpwcharrFoldReport);
	 wcscat(wcharrPath, L"\\ShipTargReport");
		_wmkdir(wcharrPath);

		wcscat(wcharrPath, L"\\");



	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_ShipTarg * 30] = {0};

	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"TrgGSK_Real_X");
	 wcscpy( &wcharrFileNames[ 2 * 30], L"TrgGSK_Real_Y");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"TrgGSK_Real_Z");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"TrgGSK_Real_VX");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"TrgGSK_Real_VY");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"TrgGSK_Real_VZ");


	 for (int i = 1; i < QUANT_COLS_BUFF_ShipTarg; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_ShipTarg // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100. //  масштаб по оси X
								  ,1.  // масштаб по оси Y
								   )  ;
	 }

	 // траектория в плоскости OXY
	TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_ShipTarg // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,1  // номер переменной по оси X
								  ,2 // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,1.  // масштаб по оси Y
								   )  ;
	// траектория в плоскости OXZ
	TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_ShipTarg // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,1  // номер переменной по оси X
								  ,3 // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,1.  // масштаб по оси Y
								   )  ;
	// траектория в плоскости OYZ
	TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_ShipTarg // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,2  // номер переменной по оси X
								  ,3 // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,1.  // масштаб по оси Y
								   )  ;

 


	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;

 }



#pragma package(smart_init)
