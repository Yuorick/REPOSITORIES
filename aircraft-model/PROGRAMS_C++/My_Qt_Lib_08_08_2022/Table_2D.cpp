//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include <stdio.h>
#include <float.h>
#include "Table_2D.h"
#include "Table_1D.h"


__fastcall TTable_2D::~TTable_2D()
{
	if(mparrArg) //delete mparrArg ;
	mparrArg = NULL ;

	if(mpArrTable_1D) //delete []mpArrTable_1D ;
	mpArrTable_1D = NULL ;


}
//---------------------------------------------------------------------------

 // оператор присваивания
 TTable_2D TTable_2D::operator=(TTable_2D  R)
 {
	memcpy(Box, R.Box, 2 * sizeof(double));
	mNumCols = R.mNumCols ;
	if (mparrArg != NULL) mparrArg = NULL ;
	if(R.mparrArg != NULL)
	{
		mparrArg = new double[R.mNumCols];
		if(mparrArg == NULL)
		{
		//ShowMessage(L"Not memory for mparrArg") ;
		//Abort() ;
		}
		memcpy( mparrArg,R.mparrArg, R.mNumCols * sizeof(double));
	}

	if (mpArrTable_1D != NULL) mpArrTable_1D = NULL ;
	if(R.mpArrTable_1D != NULL)
	{
		mpArrTable_1D = new TTable_1D[R.mNumCols];
		if(mpArrTable_1D == NULL)
		{
		//ShowMessage(L"Not memory for mpArrTable_1D") ;
		//Abort() ;
		}
		memcpy( mpArrTable_1D,R.mpArrTable_1D, R.mNumCols * sizeof(TTable_1D));
	}
	return *this ;
 }

 // конструктор копирования
 TTable_2D::TTable_2D (const TTable_2D &R)
 {
	memcpy(Box, R.Box, 2 * sizeof(double));
	mNumCols = R.mNumCols ;
	if (mparrArg != NULL) mparrArg = NULL ;
	if(R.mparrArg != NULL)
	{
		mparrArg = new double[R.mNumCols];
		if(mparrArg == NULL)
		{
		//ShowMessage(L"Not memory for mparrArg") ;
		//Abort() ;
		}
		memcpy( mparrArg,R.mparrArg, R.mNumCols * sizeof(double));
	}

	if (mpArrTable_1D != NULL) mpArrTable_1D = NULL ;
	if(R.mpArrTable_1D != NULL)
	{
		mpArrTable_1D = new TTable_1D[R.mNumCols];
		if(mpArrTable_1D == NULL)
		{
		//ShowMessage(L"Not memory for mpArrTable_1D") ;
		//Abort() ;
		}
		memcpy( mpArrTable_1D,R.mpArrTable_1D, R.mNumCols * sizeof(TTable_1D));
	}
 }
//-----------------------------------------------------------------------------------
// парам констр
//--------------------------------------------------------------------------------------
 TTable_2D ::TTable_2D()
{
	mNumCols = 0;
	mparrArg = NULL;
	mpArrTable_1D = NULL;
}

 // парам констр
TTable_2D :: TTable_2D( double *parrArg, TTable_1D *pArrTable_1D, const int NumCols)

 {
	mNumCols = NumCols ;
	mparrArg = NULL ;
	mparrArg = new double[NumCols];
	if(mparrArg == NULL)
	{
	//ShowMessage(L"Not memory for mparrArg") ;
	//Abort() ;

	}
	memcpy(mparrArg, parrArg, NumCols * sizeof(double));


	mpArrTable_1D = NULL ;
	mpArrTable_1D = new TTable_1D[NumCols];
	if(mpArrTable_1D == NULL)
	{
	//ShowMessage(L"Not memory for mpArrTable_1D") ;
	//Abort() ;
	}
	memcpy(mpArrTable_1D, pArrTable_1D, NumCols * sizeof(TTable_1D));
 }


	// парам констр 2
 //	parrArgTab1 [ NumColsTab1  ] - массив аргументов таблицы 1 (УЗП при фиксированной дальности)
 //      для морской цели это целочисленный массив от 1 до к-ва попаданий при которогм вероятность поражения равна 1
 //	parrArgTab2 [ NumColsTab2  ] - массив аргументов таблицы 2 . Это дальности, прикоторых насчитаны УЗП. Высоты фиксвированы
 //      для морской цели это массив дальностей при которых насчитанны УЗП
 // parrVal[NumColsTab1 * NumColsTab2 ] - массив вероятностей
TTable_2D :: TTable_2D( double *parrArgTab1, const int NumColsTab1, double *parrArgTab2, const int NumColsTab2, double *parrVal )

 {
	 //1.
		mNumCols =  NumColsTab2;
		///

		// 2.
		mparrArg = new double[mNumCols];
		if(mparrArg == NULL)
		{
		//ShowMessage(L"Not memory for mparrArg") ;
		//Abort() ;

	 }
		memcpy(mparrArg, parrArgTab2,mNumCols * sizeof(double));
		///

		// 3.
		mpArrTable_1D = new TTable_1D[mNumCols];
		for (int i =0; i < mNumCols ; i++)
		{
			mpArrTable_1D[i] =	TTable_1D( parrArgTab1, &parrVal[i *   NumColsTab1 ], NumColsTab1) ;

		}

 }


 //-------------------------------------------------------------------------------
 double  TTable_2D ::calcValue(const double VAlRowArg, const double VAlColArg)
 {
	 calcBoundBox();
	 if (mNumCols  == 0)
	 {
	 return 0.;
	 }
	 if (VAlRowArg < Box[0])
   {
			if (mparrArg[mNumCols-1] > mparrArg[0])
			{
			return mpArrTable_1D[0].calcValue( VAlColArg);
			}
			else
			{
			return  mpArrTable_1D[mNumCols-1].calcValue( VAlColArg);
			}
	 }

	 if (VAlRowArg > Box[1])
   {
			if (mparrArg[mNumCols-1] > mparrArg[0])
			{
			 return mpArrTable_1D[mNumCols - 1].calcValue( VAlColArg);
			}
			else
			{
			return mpArrTable_1D[0].calcValue( VAlColArg);
			}
	 }
	 ///
	 double valReturn = -1.;
	 for (int i =0; i < mNumCols -1; i++)
	 {
			if ((VAlRowArg - mparrArg[i])*(VAlRowArg- mparrArg[i +1]) <= 0.)
			{

				if ((mparrArg[i] -  mparrArg[i + 1] )< 2.* DBL_MIN)
				{
				 valReturn =  (mpArrTable_1D[i].calcValue( VAlColArg) + mpArrTable_1D[i + 1].calcValue( VAlColArg)) /2.;
				 break;
				}

				double valTang  =  (mpArrTable_1D[i].calcValue( VAlColArg)-  mpArrTable_1D[i + 1].calcValue( VAlColArg) )
				/(mparrArg[i] -  mparrArg[i + 1] );
				valReturn =  valTang * (VAlRowArg -  mparrArg[i]) + mpArrTable_1D[i].calcValue( VAlColArg);
				break;

			}
	 }

	 return valReturn;
 }

 void TTable_2D::calcBoundBox()
{
	 double xMin = 1000000000 ;
	 double xMax = -1000000000 ;


	 for (int i =0; i < mNumCols; i++)
	 {
		 if (mparrArg[i] > xMax) xMax =  mparrArg[i];
		 if (mparrArg[i] < xMin) xMin =  mparrArg[i];

	 }
	Box[0] =  xMin;
	Box[1] =  xMax;
}
#pragma package(smart_init)
