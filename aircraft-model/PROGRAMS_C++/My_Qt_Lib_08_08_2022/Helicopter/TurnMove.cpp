
#include "TurnMove.h"
#include  <string.h>
#include  <math.h>
#include "BallanceCalc.h"
#include "MatrixProccess.h"



TTurnMove::TTurnMove():TPartHelicTraj()
{

    memset(marrC0, 0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
    memset(marrSteadySolution0, 0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
    memset(marrSteadyW, 0, 4  * sizeof(long double));

    memset(miarrNumsControlledVars, 0, (QUantCurNZSKVarsVS -1)  * sizeof(int));
    mY = 0.;
    mVx = 0.;
    mAngYaw = 0.;
    mRadius = 0.;
    memset(marrRotCentre_CurNZSK, 0, 3 * sizeof(long double));

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

    mCourseEnd = M_PI/ 2.;




}

//---------------------------------------------------------------------------

// парам констр виража
TTurnMove::TTurnMove( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                              , const long  double VAlVx, const long double VAlAngYaw
                    , const long double VAlRadius,long double *arrRotCentre_CurNZSK
                    ,const long double  TBegin,const long double  CourseEnd)
  :TPartHelicTraj( Helic, Environment, TBegin, TBegin)
{
  mCourseEnd = CourseEnd;
  memset(marrC0, 0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
  memset(marrSteadySolution0, 0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
  memset(marrSteadyW, 0, 4  * sizeof(long double));

  mY = VAlY;
  mVx = VAlVx;
  mAngYaw = VAlAngYaw;
  mRadius  = VAlRadius;
  marrSteadySolution0[1] = VAlY;
  marrSteadySolution0[3] = VAlVx;

  long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
  long double arrXRez[6] = {0.};

  TBallanceCalc::calcBallParamsForTurn (Helic, Environment, VAlY, VAlVx,VAlAngYaw, VAlRadius
                                                        ,arrX0, arrXRez );

   marrSteadySolution0[9] = arrXRez[0];
   marrSteadySolution0[10] = arrXRez[1];
   marrSteadySolution0[11] = VAlAngYaw;
   TPartHelicTraj::fillUpOmega_StableTurn(VAlVx,fabsl(VAlRadius),arrXRez[0]
                                                 ,arrXRez[1], &(marrSteadySolution0[6]));
   marrSteadySolution0[9] =marrSteadySolution0[9] * SIGNUM__(VAlRadius);
   memcpy(marrSteadyW, &(arrXRez[2]), 4 * sizeof(long double));

   memcpy(marrRotCentre_CurNZSK, arrRotCentre_CurNZSK, 3 * sizeof(long double));


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

  mCurrentCourseIntegrated = 0.;

}




 //---------------------------------------------------------------------------

 // оператор присваивания
 TTurnMove TTurnMove::operator=(TTurnMove  R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTimeCur = R.mTimeCur;
     mCurrentCourseIntegrated= R.mCurrentCourseIntegrated;

     memcpy(marrC0, R.marrC0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
     memcpy(marrSteadySolution0, R.marrSteadySolution0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSteadyW, R.marrSteadyW, 4  * sizeof(long double));
     mY = R.mY;
     mVx = R.mVx;
     mAngYaw = R.mAngYaw;
     mRadius = R. mRadius;
     memcpy(marrRotCentre_CurNZSK, R.marrRotCentre_CurNZSK, 3 * sizeof(long double));
     mTBegin = R.mTBegin;


     mQuantControlledVars = R.mQuantControlledVars;
     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1)  * sizeof(int));
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1)  * sizeof(int));
     mCourseEnd = R.mCourseEnd;
     return *this;
 }

 // конструктор копирования
 TTurnMove::TTurnMove (const TTurnMove &R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTimeCur = R.mTimeCur;
     mCurrentCourseIntegrated= R.mCurrentCourseIntegrated;



     memcpy(marrC0, R.marrC0, (QUantCurNZSKVarsVS -1) * 4 * sizeof(long double));
     memcpy(marrSteadySolution0, R.marrSteadySolution0, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSteadyW, R.marrSteadyW, 4  * sizeof(long double));
     mY = R.mY;
     mVx = R.mVx;
     mAngYaw = R.mAngYaw;
     mRadius = R. mRadius;
     memcpy(marrRotCentre_CurNZSK, R.marrRotCentre_CurNZSK, 3 * sizeof(long double));
     mTBegin = R.mTBegin;

     mQuantControlledVars = R.mQuantControlledVars;
     memcpy(miarrNumsControlledVars, R.miarrNumsControlledVars, (QUantCurNZSKVarsVS -1)  * sizeof(int));
     mQuantControlledVarsTang = R.mQuantControlledVarsTang;
     memcpy(miarrNumsControlledVarsTang, R.miarrNumsControlledVarsTang, (QUantCurNZSKVarsVS -1)  * sizeof(int));
     mCourseEnd = R.mCourseEnd;


}


 //----------------------------------------------------------------------

 void TTurnMove::calc_W(long double *arrW)
 {

// положительное направления угла рыскания - против часовой стрелки при вращении относительно оси OY НЗСК !!!!!
     // угол поворота оси OX повернутой (скоростной) сиситемы координат относительно оси OX ТНЗСК
     long double val_Ang_OX_Rot = - mVx / mRadius* (mTimeCur - mTBegin);//888


     // матрица перехода из ТНЗСК в повернутую
     long double arrMtrxPer_CurNZSK_to_RotatedSK[9] = {0.};
     calcMtrxPer_CurNZSK_to_RotatedSK (val_Ang_OX_Rot,arrMtrxPer_CurNZSK_to_RotatedSK);

     // перевод фаового вектора из ТНЗСК в скоростную СК
     long double arrPhaseVect_RotatedSK[QUantCurNZSKVarsVS-1] = {0.};
     memcpy(arrPhaseVect_RotatedSK, marrPhaseVect, (QUantCurNZSKVarsVS-1) * sizeof(long double));
       // перевод вектора положения
     long double arrT[3] = {0.}, arrCurCentre[3] = {0.};
     arrCurCentre[1] = 0;
     arrCurCentre[0] =  mRadius * sinl(-val_Ang_OX_Rot);
     arrCurCentre[2] = mRadius - mRadius * cosl(-val_Ang_OX_Rot);
     MtrxMinusMatrx(marrPhaseVect, arrCurCentre, 3, 1, arrT);
     MtrxMultMatrx(arrMtrxPer_CurNZSK_to_RotatedSK,3, 3, arrT,1, arrPhaseVect_RotatedSK) ;
     ///

     // перевод вектора скорости
     MtrxMultMatrx(arrMtrxPer_CurNZSK_to_RotatedSK,3, 3, &(marrPhaseVect[3]),1, &(arrPhaseVect_RotatedSK[3])) ;
     ///
     arrPhaseVect_RotatedSK[11] = fnc_Minus_PI_Plus_PI(arrPhaseVect_RotatedSK[11] - val_Ang_OX_Rot);//8888



     long double  arrT0[QUantCurNZSKVarsVS -1] = {0.},  arrDelW[4] ={0.};

               MtrxMinusMatrx(arrPhaseVect_RotatedSK, marrSteadySolution0,1, QUantCurNZSKVarsVS -1, arrT0);
               arrT0[9] =  fnc_Minus_PI_Plus_PI(arrT0[9]);
               arrT0[10] =  fnc_Minus_PI_Plus_PI(arrT0[10]);
               arrT0[11] =  fnc_Minus_PI_Plus_PI(arrT0[11]);

               if (fabs(arrT0[2]) > 2.)
               {
                 arrT0[2] = 2. * SIGNUM__(arrT0[2]);
               }

               if (fabs(arrT0[5]) > 5.)
               {
                 arrT0[5] = 5. * SIGNUM__(arrT0[5]);
               }

               if (fabs(arrT0[4]) > 20.)
               {
                 arrT0[4] = 20. * SIGNUM__(arrT0[4]);
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
                  // arrW[3] = 5. * M_PI / 180.;
               }

 }


 //---------------------------------------------------------------------------------

 // матрица перехода из ТНЗСК в повернутую (скоростную)
void  TTurnMove::calcMtrxPer_CurNZSK_to_RotatedSK (const long double val_Ang_OX_Rot, long double *arrMtrxPer_CurNZSK_to_RotatedSK)
{
    memset (arrMtrxPer_CurNZSK_to_RotatedSK, 0., 9 * sizeof(long double));

 arrMtrxPer_CurNZSK_to_RotatedSK[4] = 1.;
 arrMtrxPer_CurNZSK_to_RotatedSK[0] = arrMtrxPer_CurNZSK_to_RotatedSK[8] = cosl(val_Ang_OX_Rot);
 arrMtrxPer_CurNZSK_to_RotatedSK[2] = -sinl(val_Ang_OX_Rot);  // 8888
 arrMtrxPer_CurNZSK_to_RotatedSK[6]  = -arrMtrxPer_CurNZSK_to_RotatedSK[2] ;// 8888
}
//----------------------------------------
void TTurnMove::get_arrSteadyW(long double *arrW)
{
  memcpy(arrW, marrSteadyW, 4 * sizeof(long double));
}
//----------------------------------------
void TTurnMove::get_arrayOfControlledVars(int *piNum, int *iarr)
{
  *piNum = mQuantControlledVars;
    memcpy(iarr, miarrNumsControlledVars, mQuantControlledVars * sizeof(int));
}

//----------------------------------------
void TTurnMove::get_arrayOfControlledVarsTang(int *piNum, int *iarr)
{
  *piNum = mQuantControlledVarsTang;
    memcpy(iarr, miarrNumsControlledVarsTang, mQuantControlledVarsTang * sizeof(int));
}

//----------------------------------------
void TTurnMove::get_QuantOfControlledVarsTang(int *piNum)
{
  *piNum = mQuantControlledVarsTang;

}


//----------------------------------------
void TTurnMove::get_QuantOfControlledVars(int *piNum)
{
  *piNum = mQuantControlledVars;

}


//----------------------------------------
void TTurnMove::findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
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
       arrC[5] = 0.05 + 0.1* ((long double) i6 );// Nu - Kappa

   for (int i5 =0; i5 < 8; i5 ++) //++
   {
       arrC[4] =  0.05 + 0.1* ((long double) i5 );// OmZ - Kappa

    for (int i0 = 0; i0 < 20; ++i0)// ++
    {


        arrC[7] = -0.0005 - 0.005* ((long double) i0 );// Y - Fi

        for (int i1 = 0; i1< 20; ++i1)
        {

            arrC[9] =  -0.0005 - 0.005*((long double) i1 );// Vy - Fi
          for (int i2 = 0; i2 <  20; ++i2)
          {

                    arrC[2] =  -0.0005    -0.0025 *((long double) i2 );// Vx - Kappa



                  for (int i4 = 0; i4 < 60; ++i4)
                  {


                     arrC[0] = -0.0005    -0.0025 *((long double) i4 ); // X - Kappa

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
// НЕ СДЕЛАНО!!!
void TTurnMove::findingCircle(long double valz1 ,  long double valz2, long double val_A
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
          arrC00[43] = -0.0201 - 0.01 * ((long double) i7 ); // это OmY - DelFi
        // arrC00[39] = -0.01 - 0.01 * ((long double) i7 );
      for (int i6 =0; i6 < 20; ++i6)
      {
          arrC00[47] = -0.001- 0.005 * ((long double) i6 );// это Psi - DelFi
      for (int i0 = 0; i0 < 20; ++i0)// i0=8, 170,250, 0.1
      {
        arrC00[30] = -0.05 - 0.025 * ((long double) i0 ); // это OmX- Etta
         // arrC00[27] = -0.2 - 0.05 * ((long double) i0 );
         //  arrC00[27] = -0.02 - 0.01 * ((long double) i0 );
          for (int i1 = 0; i1< 40; ++i1)
          {
           // arrC00[30] = -0.05 - 0.05 * ((long double) i1 );
           arrC00[33] = -0.01 - 0.05 * ((long double) i1 ); // это Gamma -Etta

                     for (int i4 =0; i4 < 10; ++i4)
                     {


                       arrC00[29] = - 0.001 -0.06* ( static_cast<long double>( i4 )) ;   // это Vz - Etta
                       //  arrC00[26] = -0.001 - 0.002 * ((long double)i4 + 1.);
                       for (int i5 =0; i5 < 10; ++i5)
                       {


                          arrC00[26] =   - 0.001 -0.06* ( static_cast<long double>( i5)) ;   // это Z -Etta
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


bool TTurnMove::IsEndOfPart()
{

 calcCurrentIntegtatedCourse();

 if ( (SIGNUM__(mRadius) *mCurrentCourseIntegrated) > (SIGNUM__(mRadius) * mCourseEnd ))
 {
     return true;
 }
 return false;
}

// вычисление текущего интегрированного курса
void TTurnMove::calcCurrentIntegtatedCourse()
{
    long double valCourseCur = atan2(marrPhaseVect[5], marrPhaseVect[3]);
    long double valCoursePrevious = fnc_Minus_PI_Plus_PI(mCurrentCourseIntegrated);
    long double valtemp = valCourseCur- valCoursePrevious;
    long double val_2PI_N = mCurrentCourseIntegrated - valCoursePrevious;

    if ( (SIGNUM(valCourseCur * valCoursePrevious) < 0.)&& (fabsl(valtemp) > 1.5 * M_PI))
    {
      mCurrentCourseIntegrated = val_2PI_N -SIGNUM__(valCourseCur) * 2. * M_PI + valCourseCur;
    }
    else
    {
      mCurrentCourseIntegrated = val_2PI_N  + valCourseCur;
    }


}


