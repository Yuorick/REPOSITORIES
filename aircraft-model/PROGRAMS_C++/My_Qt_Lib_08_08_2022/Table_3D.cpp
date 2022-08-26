//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include <stdio.h>
#include <float.h>
#include "Table_2D.h"
#include "Table_1D.h"
#include "Table_3D.h"



__fastcall TTable_3D::~TTable_3D()
{
	if(mparrArg) //delete mparrArg ;
	mparrArg = NULL ;

	if(mpArrTable_2D) //delete []mpArrTable_2D ;
	mpArrTable_2D = NULL ;


}
//---------------------------------------------------------------------------

 // оператор присваивания
 TTable_3D TTable_3D::operator=(TTable_3D  R)
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

	if (mpArrTable_2D != NULL) mpArrTable_2D = NULL ;
	if(R.mpArrTable_2D != NULL)
	{
		mpArrTable_2D = new TTable_2D[R.mNumCols];
		if(mpArrTable_2D == NULL)
		{
		//ShowMessage(L"Not memory for mpArrTable_2D") ;
		//Abort() ;
		}
		memcpy( mpArrTable_2D,R.mpArrTable_2D, R.mNumCols * sizeof(TTable_2D));
	}
	return *this ;
 }

 // конструктор копирования
 TTable_3D::TTable_3D (const TTable_3D &R)
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

	if (mpArrTable_2D != NULL) mpArrTable_2D = NULL ;
	if(R.mpArrTable_2D != NULL)
	{
		mpArrTable_2D = new TTable_2D[R.mNumCols];
		if(mpArrTable_2D == NULL)
		{
		//ShowMessage(L"Not memory for mpArrTable_2D") ;
		//Abort() ;
		}
		memcpy( mpArrTable_2D,R.mpArrTable_2D, R.mNumCols * sizeof(TTable_2D));
	}
 }
//-----------------------------------------------------------------------------------
// парам констр
//--------------------------------------------------------------------------------------
 TTable_3D::TTable_3D()
{
	mNumCols = 0;
	mparrArg = NULL;
	mpArrTable_2D = NULL;
}

 // парам констр
TTable_3D :: TTable_3D( double *parrArg, TTable_2D *pArrTable_2D, const int NumCols)

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


	mpArrTable_2D = NULL ;
	mpArrTable_2D = new TTable_2D[NumCols];
	if(mpArrTable_2D == NULL)
	{
	//ShowMessage(L"Not memory for mpArrTable_2D") ;
	//Abort() ;
	}
	memcpy(mpArrTable_2D, pArrTable_2D, NumCols * sizeof(TTable_2D));
 }

	// парам констр 2
 //	parrArgTab1 [ NumColsTab1  ] - массив аргументов таблицы 1 (УЗП при фиксированной высоте и дальности)
 //      иными словами, это массив векличин промахов, при которых насчтаны вероятности
 //	parrArgTab2 [ NumColsTab2  ] - массив аргументов таблицы 2 . Это дальности, прикоторых насчитаны УЗП. Высоты фиксвированы
 //	parrArgTab3[ NumColsTab3  ] - массив аргументов таблицы 3 . Это высоты.
 // parrVal[NumColsTab1 * NumColsTab2 * NumColsTab3] - масси в вероятностей
TTable_3D :: TTable_3D( double *parrArgTab1, const int NumColsTab1, double *parrArgTab2, const int NumColsTab2
 ,double *parrArgTab3, const int NumColsTab3, double *parrVal )

 {
	 //1.
		mNumCols =  NumColsTab3;
		///

		// 2.
		mparrArg = new double[mNumCols];
		if(mparrArg == NULL)
		{
		//ShowMessage(L"Not memory for mparrArg") ;
		//Abort() ;

	 }
		memcpy(mparrArg, parrArgTab3,mNumCols * sizeof(double));
		///

		// 3.
		mpArrTable_2D = new TTable_2D[mNumCols];
		for (int i =0; i < mNumCols ; i++)
		{
			 TTable_1D *pArrTable_1D = new TTable_1D[NumColsTab2];
				 for (int j =0; j < NumColsTab2; j++)
				 {
					pArrTable_1D[j] =	TTable_1D( parrArgTab1, &parrVal[i *  NumColsTab2 * NumColsTab1 + j * NumColsTab1], NumColsTab1) ;
				 }
			mpArrTable_2D[i] = TTable_2D( parrArgTab2, pArrTable_1D, NumColsTab2);
			 delete [] pArrTable_1D;
		}

 }

 //------------------------------------------------------
//VAlTab3Arg - высота
// VAlRowArg - дальность
// VAlColArg - ВЕЛИЧИНА ПРОМАХА
 double  TTable_3D::calcValue(const double VAlTab3Arg,const double VAlRowArg, const double VAlColArg)
 {
	 calcBoundBox();
	 if (mNumCols  == 0)
	 {
	 return 0.;
	 }
	 if (VAlTab3Arg < Box[0])
   {
			if (mparrArg[mNumCols-1] > mparrArg[0])
			{
			return mpArrTable_2D[0].calcValue(VAlRowArg, VAlColArg);
			}
			else
			{
			return  mpArrTable_2D[mNumCols-1].calcValue(VAlRowArg, VAlColArg);
			}
	 }

	 if (VAlTab3Arg > Box[1])
   {
			if (mparrArg[mNumCols-1] > mparrArg[0])
			{
			 return mpArrTable_2D[mNumCols - 1].calcValue(VAlRowArg, VAlColArg);
			}
			else
			{
			return mpArrTable_2D[0].calcValue(VAlRowArg, VAlColArg);
			}
	 }
	 ///
	 double valReturn = -1.;
	 for (int i =0; i < mNumCols -1; i++)
	 {
			if ((VAlTab3Arg - mparrArg[i])*(VAlTab3Arg- mparrArg[i +1]) <= 0.)
			{

				if ((mparrArg[i] -  mparrArg[i + 1] )< 2.* DBL_MIN)
				{
				 valReturn =  (mpArrTable_2D[i].calcValue(VAlRowArg, VAlColArg) + mpArrTable_2D[i + 1].calcValue(VAlRowArg, VAlColArg)) /2.;
				 break;
				}

				double valTang  =  (mpArrTable_2D[i].calcValue(VAlRowArg, VAlColArg)-  mpArrTable_2D[i + 1].calcValue(VAlRowArg, VAlColArg) )
				/(mparrArg[i] -  mparrArg[i + 1] );
				valReturn =  valTang * (VAlTab3Arg -  mparrArg[i]) + mpArrTable_2D[i].calcValue(VAlRowArg, VAlColArg);
				break;

			}
	 }

	 return valReturn;
 }


 //--------------------------------------------------------------------------------------------------
 void TTable_3D::calcBoundBox()
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
