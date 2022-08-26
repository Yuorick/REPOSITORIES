#include "LineMove.h"
#include  <string.h>
#include <QDebug>
#include "BallanceCalc.h"
#include "MatrixProccess.h"



TLineMove::TLineMove():TPartHelicTraj()
{

    memset(marrC0, 0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
    memset(marrSteadySolution0, 0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
    memset(marrSteadyW, 0, 4  * sizeof(long double));
    mQuantControlledVars = 11;
    memset(miarrNumsControlledVars, 0, (QUantCurNZSKVarsVS -1)  * sizeof(int));
    mY = 0.;
    mVx = 0.;
    mAngYaw = 0.;

    mQuantControlledVars = QUantCurNZSKVarsVS -2;
    for (int i =0; i < (QUantCurNZSKVarsVS -1); ++i)
    {
        miarrNumsControlledVars[i-1] = -1;
        miarrNumsControlledVarsTang[i-1] = -1;

    }

    for (int i =1; i < (QUantCurNZSKVarsVS -1); ++i)
    {
        miarrNumsControlledVars[i-1] = i;

    }

    mQuantControlledVarsTang = 5;

    miarrNumsControlledVarsTang[0] = 1;
    miarrNumsControlledVarsTang[1] = 3;
    miarrNumsControlledVarsTang[2] = 4;
    miarrNumsControlledVarsTang[3] = 8;
    miarrNumsControlledVarsTang[4] = 10;


}

//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

 // парам констр равном прямолин движения
TLineMove::TLineMove( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                                , const long  double VAlVx, const long double VAlAngYaw,const long double  TBegin)
    :TPartHelicTraj( Helic, Environment, TBegin, TBegin)
{

    memset(marrC0, 0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
    memset(marrSteadySolution0, 0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
    memset(marrSteadyW, 0, 4  * sizeof(long double));


    mY = VAlY;
    mVx = VAlVx;
    mAngYaw = VAlAngYaw;
    marrSteadySolution0[1] = VAlY;
    marrSteadySolution0[3] = VAlVx;    

    long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez[6] = {0.};


    TBallanceCalc::calcBallParamsForSteadyLineMoving (Helic, Environment, VAlY, VAlVx,VAlAngYaw
                                                          ,arrX0, arrXRez );


    marrSteadySolution0[11] = VAlAngYaw;


     marrSteadySolution0[9] = arrXRez[0];
     marrSteadySolution0[10] = arrXRez[1];
     memcpy(marrSteadyW, &(arrXRez[2]), 4 * sizeof(long double));

     mQuantControlledVars = QUantCurNZSKVarsVS -2;
     for (int i =0; i < (QUantCurNZSKVarsVS -1); ++i)
     {
         miarrNumsControlledVars[i] = -1;
         miarrNumsControlledVarsTang[i] = -1;

     }

     for (int i =1; i < (QUantCurNZSKVarsVS -1); ++i)
     {
         miarrNumsControlledVars[i-1] = i;

     }

     mQuantControlledVarsTang = 5;

     miarrNumsControlledVarsTang[0] = 1;
     miarrNumsControlledVarsTang[1] = 3;
     miarrNumsControlledVarsTang[2] = 4;
     miarrNumsControlledVarsTang[3] = 8;
     miarrNumsControlledVarsTang[4] = 10;



}



 //---------------------------------------------------------------------------

 // оператор присваивания
 TLineMove TLineMove::operator=(TLineMove  R)
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
     mVx = R.mVx;
     mAngYaw = R.mAngYaw;

     mQuantControlledVars = R.mQuantControlledVars;
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;

     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1) * sizeof(int));
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1) * sizeof(int));


     return *this;
 }

 // конструктор копирования
 TLineMove::TLineMove (const TLineMove &R)
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
     mVx = R.mVx;
     mAngYaw = R.mAngYaw;

     mQuantControlledVars = R.mQuantControlledVars;
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;

     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1) * sizeof(int));
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1) * sizeof(int));

}

 void TLineMove::calc_W(long double *arrW)
 {
     long double  arrT0[QUantCurNZSKVarsVS -1] = {0.},  arrDelW[4] ={0.};

               MtrxMinusMatrx(marrPhaseVect, marrSteadySolution0,1, QUantCurNZSKVarsVS -1, arrT0);
               if (fabs(arrT0[11]) > 1.)
               {
                 arrT0[11] =  SIGNUM__(arrT0[1]);
               }

               if (fabs(arrT0[1]) > 30.)
               {
                 arrT0[1] = 30. * SIGNUM__(arrT0[1]);
               }
               if (fabs(arrT0[2]) > 40.)
               {
                 arrT0[2] = 40. * SIGNUM__(arrT0[2]);
               }
               if (fabs(arrT0[5]) > 5.)
               {
                 arrT0[5] = 5. * SIGNUM__(arrT0[5]);
               }

               if (fabs(arrT0[4]) > 20.)
               {
                 arrT0[4] = 20. * SIGNUM__(arrT0[4]);
               }

               if (fabs(arrT0[11]) > 20.)
               {
                 arrT0[11] = 0.2 * SIGNUM__(arrT0[11]);
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
 //-------------------------------------------------------------------------

 void TLineMove::get_arrayOfControlledVars(int *piNum, int *iarr)
 {
   *piNum = mQuantControlledVars;
     memcpy(iarr, miarrNumsControlledVars, mQuantControlledVars * sizeof(int));
 }

 //-------------------------------------------------------------------------

 void TLineMove::get_arrayOfControlledVarsTang(int *piNum, int *iarr)
 {
   *piNum = mQuantControlledVarsTang;
     memcpy(iarr, miarrNumsControlledVarsTang, mQuantControlledVarsTang * sizeof(int));
 }

 //-------------------------------------------------------------------------

 void TLineMove::get_arrSteadyW(long double *arrW)
 {
   memcpy(arrW, marrSteadyW, 4 * sizeof(long double));
 }
 //-------------------------------------------------------------------------
  void TLineMove::get_QuantOfControlledVarsTang(int *piNum)
  {
    *piNum = mQuantControlledVarsTang;

  }

  //----------------------------------------
  //-------------------------------------------------------------------------
   void TLineMove::get_QuantOfControlledVars(int *piNum)
   {
     *piNum = mQuantControlledVars;

   }

   //----------------------------------------
 void TLineMove::findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
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

   for (int i6 =0; i6 < 8; i6 ++) //++
   {
         arrC[3] =  0.05 + 0.1* ((long double) i6 );// OmZ - Kappa

    for (int i5 =0; i5 < 8; i5 ++) //++
    {
        arrC[4] = 0.05 + 0.1* ((long double) i5 );// Nu - Kappa


     for (int i0 = 0; i0 < 20; ++i0)// ++
     {

         arrC[5] = -0.0005 - 0.005* ((long double) i0 );// Y - Fi

         for (int i1 = 0; i1< 20; ++i1)
         {
             arrC[7] =  -0.0005 - 0.005*((long double) i1 );// Vy - Fi

           for (int i2 = 0; i2 <  20; ++i2)
           {
               arrC[1] =  -0.0005    -0.0025 *((long double) i2 );// Vx - Kappa

                     MtrxMultMatrx(arr_dF_po_dW, mQuantControlledVarsTang, QuantControls, arrC, mQuantControlledVarsTang, arrT0) ;

                     MtrxSumMatrx(arr_dF_po_dx,arrT0 ,mQuantControlledVarsTang, mQuantControlledVarsTang, arrT1) ;

                     if(TPartHelicTraj::IsStability_( valz1 , valz2, val_A,arrT1, mQuantControlledVarsTang))
                     {
                         TPartHelicTraj::insertLocalGearMtrx_In_TotalGearMtrx(
                                    mQuantControlledVarsTang, miarrNumsControlledVarsTang
                                    ,QuantControls, iarrNumsControls,arrC
                                    , &(arrDataBuf[(QUantCurNZSKVarsVS -1) *4 * quantRows]));
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

   *pquantRows = quantRows;
   delete []arrC  ;
   delete []arrT0  ;
   delete []arrT1  ;
 }


 //---------------------------------------------------------------------------------------------
 void TLineMove::findingCircle(long double valz1 ,  long double valz2, long double val_A
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
     for (int i7 = 0; i7 < 10; ++i7)// 4
     {
           arrC00[39] = -0.0201 - 0.01 * ((long double) i7 ); // это OmY - DelFi
         // arrC00[39] = -0.01 - 0.01 * ((long double) i7 );
       for (int i6 =0; i6 < 20; ++i6)
       {
           arrC00[43] = -0.001- 0.005 * ((long double) i6 );// это Psi - DelFi
       for (int i0 = 0; i0 < 30; ++i0)// i0=8, 170,250, 0.1
       {
         arrC00[27] = -0.05 - 0.0025 * ((long double) i0 ); // это OmX- Etta
          // arrC00[27] = -0.2 - 0.05 * ((long double) i0 );
          //  arrC00[27] = -0.02 - 0.01 * ((long double) i0 );
           for (int i1 = 0; i1< 20; ++i1)
           {

            arrC00[30] = -0.01 - 0.01 * ((long double) i1 ); // это Gamma -Etta

                      for (int i4 =0; i4 < 10; ++i4)
                      {


                        arrC00[26] = - 0.001 -0.01* ( static_cast<long double>( i4 )) ;   // это Vz - Etta
                        //  arrC00[26] = -0.001 - 0.002 * ((long double)i4 + 1.);
                        for (int i5 =0; i5 < 20; ++i5)
                        {


                           arrC00[23] =   - 0.001 -0.01* ( static_cast<long double>( i5)) ;   // это Z -Etta
                          // arrC00[23] =-0.0001 - 0.0002 *((long double) i5 );
                           MtrxMultMatrx(arr_dF_po_dW, mQuantControlledVars, 4, arrC00, mQuantControlledVars, arrT00) ;

                           MtrxSumMatrx(arr_dF_po_dx,arrT00 ,mQuantControlledVars, mQuantControlledVars, arrT10) ;
                           if(TPartHelicTraj::IsStability_( valz1 ,  valz2, val_A,arrT10, mQuantControlledVars))
                           {
                               //for (int j=0; j < 44; ++j)
                               //{
                              //   arrDataBuf [ 44 * quantRows +j] = arrC00[j] ;
                              // }
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

 //------------------------------
 bool TLineMove::IsEndOfPart()
 {
     long double devV2 = (marrPhaseVect[3] - mVx)* (marrPhaseVect[3] - mVx) + marrPhaseVect[4] * marrPhaseVect[4]+marrPhaseVect[5] * marrPhaseVect[5];
    // long double devS2 = (marrPhaseVect[1] - mY) * (marrPhaseVect[1] - mY)  + marrPhaseVect[2] * marrPhaseVect[2];
     if (( marrPhaseVect[0] >= 0.) &&(devV2< 0.01) )
     {
         return true;
     }
     return false;
 }

