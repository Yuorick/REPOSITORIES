#include "Hover.h"
#include  <string.h>
#include "MatrixProccess.h"

#include  <string.h>
#include  <math.h>
#include "BallanceCalc.h"
#include "MatrixProccess.h"



THover::THover():TPartHelicTraj()
{

    memset(marrC0, 0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
    memset(marrSteadySolution0, 0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
    memset(marrSteadyW, 0, 4  * sizeof(long double));

    memset(miarrNumsControlledVars, 0, (QUantCurNZSKVarsVS -1)  * sizeof(int));
    mY = 0.;
    mTHovering = 0.;
    //mPsi0 = 0.;

    mQuantControlledVars = QUantCurNZSKVarsVS -1;
    for (int i =0; i < (QUantCurNZSKVarsVS -1); ++i)
    {
        miarrNumsControlledVars[i] = i;
        miarrNumsControlledVarsTang[i] = -1;

    }

    mQuantControlledVarsTang = 6;
    miarrNumsControlledVarsTang[0] = 0;
    miarrNumsControlledVarsTang[1] = 1;
    miarrNumsControlledVarsTang[2] = 3;
    miarrNumsControlledVarsTang[3] = 4;
    miarrNumsControlledVarsTang[4] = 8;
    miarrNumsControlledVarsTang[5] = 10;




}

//---------------------------------------------------------------------------

// парам констр вращения
THover::THover( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                    , const long double VAlTHovering,const long double  TBegin)
  :TPartHelicTraj( Helic, Environment, TBegin, TBegin)
{

  memset(marrC0, 0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
  memset(marrSteadySolution0, 0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
  memset(marrSteadyW, 0, 4  * sizeof(long double));

  mY = VAlY;
  mTHovering = VAlTHovering;
 // mPsi0 = VAlPsi0;

  marrSteadySolution0[1] = VAlY;
  marrSteadySolution0[3] = 0.;

  long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
  long double arrXRez[6] = {0.};


  TBallanceCalc::calcBallParamsForSteadyLineMoving(Helic, Environment, VAlY, 0. ,0.
                                                ,arrX0, arrXRez );
   marrSteadySolution0[9] = arrXRez[0];
   marrSteadySolution0[10] = arrXRez[1];
   marrSteadySolution0[11] = 0. ;


   memcpy(marrSteadyW, &(arrXRez[2]), 4 * sizeof(long double));




   mQuantControlledVars = QUantCurNZSKVarsVS -1;
   for (int i =0; i < (QUantCurNZSKVarsVS -1); ++i)
   {
       miarrNumsControlledVars[i] = i;
       miarrNumsControlledVarsTang[i] = -1;

   }


   mQuantControlledVarsTang = 6;
   miarrNumsControlledVarsTang[0] = 0;
   miarrNumsControlledVarsTang[1] = 1;
   miarrNumsControlledVarsTang[2] = 3;
   miarrNumsControlledVarsTang[3] = 4;
   miarrNumsControlledVarsTang[4] = 8;
   miarrNumsControlledVarsTang[5] = 10;



}




 //---------------------------------------------------------------------------

 // оператор присваивания
 THover THover::operator=(THover  R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));

     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTimeCur = R.mTimeCur;

     memcpy(marrC0, R.marrC0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
     memcpy(marrSteadySolution0, R.marrSteadySolution0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSteadyW, R.marrSteadyW, 4  * sizeof(long double));
     mY = R.mY;
     mTHovering = R.mTHovering;
     //mPsi0 = R.mPsi0;

     mTBegin = R.mTBegin;
     mQuantControlledVars = R.mQuantControlledVars;
     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1)  * sizeof(int));
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1)  * sizeof(int));

     return *this;
 }

 // конструктор копирования
 THover::THover (const THover &R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));

     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTimeCur = R.mTimeCur;

     memcpy(marrC0, R.marrC0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
     memcpy(marrSteadySolution0, R.marrSteadySolution0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSteadyW, R.marrSteadyW, 4  * sizeof(long double));
     mY = R.mY;
     mTHovering = R.mTHovering;
     //mPsi0 = R.mPsi0;

     mTBegin = R.mTBegin;

     mQuantControlledVars = R.mQuantControlledVars;
     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1)  * sizeof(int));
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1)  * sizeof(int));

}


 //----------------------------------------------------------------------

 void THover::calc_W(long double *arrW)
 {

     long double  arrT0[QUantCurNZSKVarsVS -1] = {0.},  arrDelW[4] ={0.};

               MtrxMinusMatrx(marrPhaseVect, marrSteadySolution0,1, QUantCurNZSKVarsVS -1, arrT0);
               if (fabs(arrT0[11]) > 2.)
               {
               arrT0[11] = 2. * SIGNUM__(arrT0[11]);
               }

               if (fabs(arrT0[1]) > 10.)
               {
                 arrT0[1] = 10. * SIGNUM__(arrT0[1]);
               }
               if (fabs(arrT0[2]) > 40.)
               {
                 arrT0[2] = 40. * SIGNUM__(arrT0[2]);
               }
               if (fabs(arrT0[5]) > 5.)
               {
                 arrT0[5] = 5. * SIGNUM__(arrT0[5]);
               }

               if (fabs(arrT0[4]) > 10.)
               {
                 arrT0[4] = 10. * SIGNUM__(arrT0[4]);
               }



               MtrxMultMatrx(marrC0,4,  QUantCurNZSKVarsVS -1, arrT0,1, arrDelW) ;
               MtrxSumMatrx(marrSteadyW, arrDelW,1, 4, arrW) ;

               if (fabsl(arrW[1]) > mHelic.mFiMax)
               {
                   arrW[1]=  SIGNUM__(arrW[1])* mHelic.mFiMax;
               }

               if (fabsl(arrW[0]) > mHelic.mKappaTettaMax)
               {
                   arrW[0] = SIGNUM__(arrW[0])* mHelic.mKappaTettaMax;
               }

               if (fabsl(arrW[2]) > mHelic.mKappaTettaMax)
               {
                   arrW[2] = SIGNUM__(arrW[2])* mHelic.mKappaTettaMax;
               }

               if (fabsl(arrW[3]) > 5. * M_PI / 180.)
               {
                   arrW[3] = 5. * M_PI / 180.;
               }
 }

 //---------------------------------------------------------------------------------


//----------------------------------------
void THover::get_arrSteadyW(long double *arrW)
{
  memcpy(arrW, marrSteadyW, 4 * sizeof(long double));
}
//----------------------------------------
void THover::get_arrayOfControlledVars(int *piNum, int *iarr)
{
  *piNum = mQuantControlledVars;
    memcpy(iarr, miarrNumsControlledVars, mQuantControlledVars * sizeof(int));
}

//----------------------------------------
void THover::get_arrayOfControlledVarsTang(int *piNum, int *iarr)
{
  *piNum = mQuantControlledVarsTang;
    memcpy(iarr, miarrNumsControlledVarsTang, mQuantControlledVarsTang * sizeof(int));
}

//----------------------------------------
void THover::get_QuantOfControlledVarsTang(int *piNum)
{
  *piNum = mQuantControlledVarsTang;

}


//----------------------------------------
void THover::get_QuantOfControlledVars(int *piNum)
{
  *piNum = mQuantControlledVars;

}


//----------------------------------------
void THover::findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                    ,double *arrDataBuf ,const int maxQuant, int *pquantRows)
{
    long double *arrC = new long double[2* mQuantControlledVarsTang] ;
    long double *arrT0 = new long double[mQuantControlledVarsTang * mQuantControlledVarsTang] ;
    long double *arrT1 = new long double[mQuantControlledVarsTang * mQuantControlledVarsTang] ;



  memset(arrC, 0, 2* mQuantControlledVarsTang * sizeof(long double));
  memset(arrT0, 0, mQuantControlledVarsTang* mQuantControlledVarsTang * sizeof(long double));
  memset(arrT1, 0, mQuantControlledVarsTang *mQuantControlledVarsTang * sizeof(long double));

  bool bend = false;
  int quantRows = 0;

  const int QuantControls = 2;
  int iarrNumsControls[] = {0,1};
  for (int i6 = 0; i6 < 11; i6 ++) //++
  {
       arrC[5] = 0.05 + 0.1* ((long double) i6 );// Nu - Kappa

   for (int i5 =   0; i5 < 11; i5 ++) //++
   {
       arrC[4] =  0.05 + 0.1* ((long double) i5 );// OmZ - Kappa

    for (int i0 = 0; i0 < 10; ++i0)// ++
    {


        arrC[7] = -0.0005 - 0.01* ((long double) i0 );// Y - Fi

        for (int i1 = 0; i1< 10; ++i1)
        {

            arrC[9] =  -0.0005 - 0.01*((long double) i1 );// Vy - Fi
          for (int i2 = 0; i2 <  40; ++i2)
          {

                    arrC[2] =  -0.0005    -0.0025 *((long double) i2 );// Vx - Kappa



                  for (int i4 = 0; i4 < 60; ++i4)
                  {


                     arrC[0] = -0.0005    -0.005 *((long double) i4 ); // X - Kappa

                    MtrxMultMatrx(arr_dF_po_dW, mQuantControlledVarsTang, QuantControls, arrC, mQuantControlledVarsTang, arrT0) ;

                    MtrxSumMatrx(arr_dF_po_dx,arrT0 ,mQuantControlledVarsTang, mQuantControlledVarsTang, arrT1) ;
                  //  MtrxMinusMatrx(arr_dF_po_dx,arrT0 ,5, 5, arrT1) ;
                    if(TPartHelicTraj::IsStability_( valz1 , valz2, val_A,arrT1, mQuantControlledVarsTang))
                    {
                       // for (int j =0; j < 2 * mQuantControlledVarsTang; ++j)
                        //{
                        //   arrDataBuf [ 2 * mQuantControlledVarsTang * quantRows + j] =arrC[j];
                       // }
                        TPartHelicTraj::insertLocalGearMtrx_In_TotalGearMtrx(
                            mQuantControlledVarsTang, miarrNumsControlledVarsTang
                            ,QuantControls, iarrNumsControls
                            ,arrC, &(arrDataBuf[(QUantCurNZSKVarsVS -1) * 4 *quantRows]));
                        quantRows++;
                        if (quantRows == maxQuant)
                        {

                            bend = true;
                            break;
                        }
                    }
                    else
                    {
                      int jjj = 0;
                    }

                  }
                  if(bend)
                  {
                      break;
                  }
                }


          if(bend)
          {
              break;
          }
        }
        if(bend)
        {
            break;
        }
    }
    if(bend)
    {
        break;
    }

   }
   if(bend)
   {
       break;
   }
  }

  *pquantRows = quantRows;
  delete []arrC  ;
  delete []arrT0  ;
  delete []arrT1  ;
}

//--------------------------------------------------------------

void THover::findingCircle(long double valz1 ,  long double valz2, long double val_A
                  , long double *arrC00, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                 ,double *arrDataBuf ,const int maxQuant, int *pquantRows)
{

    bool bend = false;
    int quantRows = 0;
    int iarrNumsControls[4] = {0,1,2,3};

    long double *arrT00 = new long double[mQuantControlledVars * mQuantControlledVars];
    long double *arrT10 = new long double[mQuantControlledVars * mQuantControlledVars];
    memset(arrT00, 0, mQuantControlledVars * mQuantControlledVars * sizeof(long double));
    memset(arrT10, 0, mQuantControlledVars * mQuantControlledVars * sizeof(long double));

    //,0.000000		,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
   // ,-0.003680	,0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
  //  ,0.000000		,0.000000	,0.000000	,-0.101531	,-0.150000	,0.000000	,0.000000	,-0.500000	,0.000000	,0.000000
   // ,0.000000		,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000
    for (int i7 = 4; i7 < 30; ++i7)// 4
    {
          arrC00[43] = -0.001 - 0.005 * ((long double) i7 ); // это OmY - DelFi

      for (int i6 = 0; i6 < 20; ++i6)
      {

          arrC00[47] = -0.001- 0.005 * ((long double) i6 );// это Psi - DelFi

      for (int i0 = 0; i0 < 30; ++i0)//
      {
        arrC00[30] = -0.01 - 0.025 * ((long double) i0 ); // это OmX- Etta


          for (int i1 = 0; i1< 40; ++i1)
          {
            arrC00[33] = -0.01 - 0.025 * ((long double) i1 ); // это Gamma -Etta
          // arrC00[33] = -0.01 - 0.01 * ((long double) i1 ); // это Gamma -Etta

                     for (int i4 = 0; i4 < 20; ++i4)
                     {

                        arrC00[29] = - 0.001 -0.006* ( static_cast<long double>( i4 )) ;
                      // arrC00[29] = - 0.001 -0.03* ( static_cast<long double>( i4 )) ;   // это Vz - Etta

                       for (int i5 = 0; i5 < 20; ++i5)
                       {


                          arrC00[26] =   - 0.001 -0.006* ( static_cast<long double>( i5)) ;   // это Z -Etta
                         // arrC00[23] =-0.0001 - 0.0002 *((long double) i5 );
                          MtrxMultMatrx(arr_dF_po_dW, mQuantControlledVars, 4, arrC00, mQuantControlledVars, arrT00) ;

                          MtrxSumMatrx(arr_dF_po_dx,arrT00 ,mQuantControlledVars, mQuantControlledVars, arrT10) ;
                          if(TPartHelicTraj::IsStability_( valz1 ,  valz2, val_A,arrT10, mQuantControlledVars))
                          {

                              TPartHelicTraj::insertLocalGearMtrx_In_TotalGearMtrx(mQuantControlledVars, miarrNumsControlledVars
                                            ,4, iarrNumsControls,arrC00, &(arrDataBuf[quantRows * 4 *(QUantCurNZSKVarsVS -1)]));

                          quantRows++;
                          if (quantRows== maxQuant)
                          {

                          bend = true;
                          break;
                          }

                         }
                         else
                         {
                          int iii = 0;
                          }

                     }
                       if(bend)
                       {
                           break;
                       }
                  }



            if(bend)
            {
                break;
            }
          }
          if(bend)
          {
              break;
          }
      }
          if(bend)
          {
              break;
          }

      }
          if(bend)
          {
              break;
          }

   }
    *pquantRows = quantRows;
    delete []arrT00;
    delete []arrT10;

}


bool THover::IsEndOfPart()
{
    if (mTimeCur > (mTBegin + mTHovering ))
    {
        return true;
    }
    return false;
}




/*
THover::THover():TLineMove()
{

    mTHovering = 0.;
   // mPsi0 = 0.;

}

//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

 // парам констр равном прямолин движения
THover::THover( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                                , const long  double VAlVx, const long double VAlAngYaw
                , const long double VAlTHovering, const long double VAlTBegin)
    :TLineMove( Helic, Environment,  VAlY,  VAlVx,  VAlAngYaw, VAlTBegin)
{

    mTHovering = VAlTHovering;
  //  mPsi0 = 0.;

}



 //---------------------------------------------------------------------------

 // оператор присваивания
 THover THover::operator=(THover  R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTHovering = R.mTHovering;
  //   mPsi0 = R.mPsi0;
     mTBegin = R.mTBegin;
     mTimeCur =  R.mTimeCur;

     memcpy(marrC0, R.marrC0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
     memcpy(marrSteadySolution0, R.marrSteadySolution0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSteadyW, R.marrSteadyW, 4  * sizeof(long double));
     mY = R.mY;
     mVx = R.mVx;
     mAngYaw = R.mAngYaw;

     mQuantControlledVars = R.mQuantControlledVars;
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;

     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1) * sizeof(int));
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1) * sizeof(int));


     return *this;
 }

 // конструктор копирования
 THover::THover (const THover &R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTHovering = R.mTHovering;
   //  mPsi0 = R.mPsi0;
     mTBegin = R.mTBegin;
     mTimeCur =  R.mTimeCur;

     memcpy(marrC0, R.marrC0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
     memcpy(marrSteadySolution0, R.marrSteadySolution0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSteadyW, R.marrSteadyW, 4  * sizeof(long double));
     mY = R.mY;
     mVx = R.mVx;
     mAngYaw = R.mAngYaw;

     mQuantControlledVars = R.mQuantControlledVars;
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;

     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1) * sizeof(int));
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1) * sizeof(int));

}

 bool THover::IsEndOfPart()
 {
    if( (mTimeCur > (mTBegin + mTHovering )) && (Norm3(&(marrPhaseVect[3])) < 0.01) )
    {
        return true;
    }
    return false;
 }

*/
