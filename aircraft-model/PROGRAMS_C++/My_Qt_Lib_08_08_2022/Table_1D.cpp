#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "Table_1D.h"
#include "Gauss.h"

TTable_1D::~TTable_1D()
{
    if(mparrArg) delete []mparrArg ;
	mparrArg = NULL ;

    if(mparrVal) delete []mparrVal ;
	mparrVal = NULL ;


}
//---------------------------------------------------------------------------

 // оператор присваивания
 TTable_1D &TTable_1D::operator=(const TTable_1D  R)
 {
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

	if (mparrVal != NULL) mparrVal = NULL ;
	if(R.mparrVal != NULL)
	{
		mparrVal = new double[R.mNumCols];
		if(mparrVal == NULL)
		{
		//ShowMessage(L"Not memory for mparrVal") ;
		//Abort() ;
		}
		memcpy( mparrVal,R.mparrVal, R.mNumCols * sizeof(double));
	}
    memcpy(Box, R.Box, 4 * sizeof(double));
	return *this ;
 }

 // конструктор копирования
 TTable_1D::TTable_1D (const TTable_1D &R)
 {
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

	if (mparrVal != NULL) mparrVal = NULL ;
	if(R.mparrVal != NULL)
	{
		mparrVal = new double[R.mNumCols];
		if(mparrVal == NULL)
		{
		//ShowMessage(L"Not memory for mparrVal") ;
		//Abort() ;
		}
		memcpy( mparrVal,R.mparrVal, R.mNumCols * sizeof(double));
        memcpy(Box, R.Box, 4 * sizeof(double));

	}
 }
//-----------------------------------------------------------------------------------
// парам констр
//--------------------------------------------------------------------------------------
 TTable_1D ::TTable_1D()
{
	mNumCols = 0;
	mparrArg = NULL;
	mparrVal = NULL;
}

 //-------------------------------------------------------
  // парам констр
 TTable_1D :: TTable_1D( double *parrArg, double *parrVal, const int NumCols)
 {
     mNumCols = NumCols ;
     mparrArg = new double[NumCols];
     if(mparrArg == NULL)
     {
     //ShowMessage(L"Not memory for mparrArg") ;
     //Abort() ;

     }
     mparrVal = new double[NumCols];
     if(mparrVal == NULL)
     {
     //ShowMessage(L"Not memory for mparrVal") ;
     //Abort() ;
     }
     for (int i = 0; i < NumCols; ++i)
     {
         mparrArg[i] = parrArg[ i];
         mparrVal[i] = parrVal[ i];
     }
  }
 //-------------------------------------------------------
 // парам констр
TTable_1D :: TTable_1D( double *parrBuff, const int NumCols)
{
    mNumCols = NumCols ;
    mparrArg = new double[NumCols];
    if(mparrArg == NULL)
    {
    //ShowMessage(L"Not memory for mparrArg") ;
    //Abort() ;

    }
    mparrVal = new double[NumCols];
    if(mparrVal == NULL)
    {
    //ShowMessage(L"Not memory for mparrVal") ;
    //Abort() ;
    }
    for (int i = 0; i < NumCols; ++i)
    {
        mparrArg[i] = parrBuff[2 * i   ];
        mparrVal[i] = parrBuff[2 * i +1];
    }
 }
//-----------------------------------
 double  TTable_1D ::calcValue(const double VAlArg)
 {


     return calcLinearValueApprox(VAlArg);
 }

 void TTable_1D::calcBoundBox()
 {
      double xMin = 1000000000 ;
      double xMax = -1000000000 ;
      double yMin =  1000000000 ;
      double yMax =  - 1000000000 ;
      for (int i =0; i < mNumCols; i++)
      {
        if (mparrArg[i] > xMax) xMax =  mparrArg[i];
        if (mparrArg[i] < xMin) xMin =  mparrArg[i];
        if (mparrVal[i] > yMax) yMax =  mparrVal[i];
        if (mparrVal[i] < yMin) yMin =  mparrVal[i];

      }
     Box[0] =  xMin;
     Box[1] =  yMin;
     Box[2] =  xMax ;
     Box[3] =  yMax;

 }

// линейная аппроксимация
 double TTable_1D::calcLinearValueApprox(const double x)
 {
      calcBoundBox();
      if (mNumCols == 0)
      {
      return 0.;
      }
      if (x < Box[0])
    {
             if (mparrArg[mNumCols-1] > mparrArg[0])
             {
             return mparrVal[0];
             }
             else
             {
             return mparrVal[mNumCols-1];
             }
      }

    if (x > Box[2])
    {
             if (mparrArg[mNumCols-1] > mparrArg[0])
             {
              return mparrVal[mNumCols-1];
             }
             else
             {
             return mparrVal[0];
             }
      }
    ///
    double valReturn = -1.;                               // 11.01.2018
    for (int i =0; i < mNumCols -1; i++)
    {
             if ((x - mparrArg[i])*(x - mparrArg[i +1]) <= 0.)
             {

                 if ((fabs(mparrArg[i] -  mparrArg[i + 1] )< 2.* DBL_MIN))
                 {
                  valReturn =  (mparrVal[i]  +mparrVal[i +1]) /2.;  // 11.01.2018
                  break;                                             // 11.01.2018
                //	return (Points[i].Y  + Points[i +1].Y) /2.;     // 11.01.2018
                 }

                 double valTang  =  (mparrVal[i] -  mparrVal[i + 1] )/(mparrArg[i] -  mparrArg[ i + 1] );  // 11.01.2018
                 valReturn =  valTang * (x -  mparrArg[i]) + mparrVal[i];                                     // 11.01.2018
                 break;                                                                                    // 11.01.2018
                //	return   valTang * (x -  Points[i].X) + Points[i].Y;      // 11.01.2018
             }
      }

    return valReturn;       // 11.01.2018

 }
 //----------------------------------
 // вычисление номера оторезка в котором лежит x
 //если x <= min возвращает -1,
 // если x >= max возвращает -2
int TTable_1D::getSegmentNum(const double x)
  {
       calcBoundBox();
       if (mNumCols == 0)
       {
       return 0.;
       }
       if (x < Box[0])
       {
           return -1;
       }

     if (x > Box[2])
     {
             return -2;
      }
     ///
     for (int i =0; i < mNumCols -1; i++)
     {
              if ((x - mparrArg[i])*(x - mparrArg[i +1]) <= 0.)
              {
                  return i;
              }
       }

  return -3;
  }
//------------------------
double TTable_1D::calc_d_po_dx(const double x)
{
  int i = getSegmentNum( x);
  return (mparrVal[i + 1] - mparrVal[i])/(mparrArg[i + 1] - mparrArg[i]);
}
//------------------------------------------
void TTable_1D::multValue(const double coeff)
{
    for (int i =0; i <mNumCols; ++i )
    {
       mparrVal[i] *= coeff;
    }
}

//------------------------------------------
void TTable_1D::addValue(const double val)
{
    for (int i =0; i <mNumCols; ++i )
    {
       mparrVal[i] += val;
    }
}
//------------------------------------------
void TTable_1D::multRandValue(const double sig)
{
    for (int i =0; i <mNumCols; ++i )
    {
        double coeff =  getGauss( 1., sig );
       mparrVal[i] *= coeff;
    }
}
//-------------------------------------
TTable_1D TTable_1D::makeMiddleProfile()
{
    double parrArg[2] = {0.};
    double parrVal[2] = {0.};
    parrArg[0] = mparrArg[0];
    parrArg[1] = mparrArg[mNumCols -1];
    double valt = calcMiddleValue();
    parrVal[0] = valt;
    parrVal[1] = valt - 0.001;
    return TTable_1D(parrArg, parrVal, 2);
}
//-------------------------------
double TTable_1D::calcMiddleValue()
{
    double sum = 0.;
    for(int i =0; i< mNumCols -1; ++i)
    {
        sum += (mparrVal[i] + mparrVal[i +1]) * (mparrArg[i +1] - mparrArg[i ])/ 2.;
    }
    return sum / (mparrArg[mNumCols-1] - mparrArg[0]);
}



