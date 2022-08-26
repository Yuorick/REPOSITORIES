#include <dir.h>
#include "PeaceSins.h"

 #include <stdio.h>
 #include <math.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>

 #include <stdlib.h>
 #include <string.h>

 #include <float.h>
 #include "Environment.h"

 #include "YrWriteShapeFile.h"
 #include "Gauss.h"



 QPeaceSins ::QPeaceSins()
{
	mSig_Q      =  0.000582;
	mSig_Psi    = 0.000145;
	mSig_Tet    = 0.000145;
	mSig_dQdt   = 0.00116 ;
	mSig_dPsidt = 0.00116 ;
	mSig_dTetdt = 0.00116 ;
	mMaxSig_Q      =  0.000582;
	mMaxSig_Psi    = 0.000145;
	mMaxSig_Tet    = 0.000145;
	mMaxSig_dQdt   = 0.00116 ;
	mMaxSig_dPsidt = 0.00116 ;
	mMaxSig_dTetdt = 0.00116 ;
	mMaxSig_H = 0.1 ;
	mMaxSig_VH = 0.05 ;
	mK1         = 0.01 ;
	mSigV       = 0.2 * sqrt(2.) ;
   //	mSig_H       = 0.1 ;
     mTCur      = 0. ;
	mDelQ       = 0. ;
	mDelVQ      = 0. ;
	mDelPsi     = 0. ;
	mDelVPsi    = 0. ;
	mDelTet     = 0. ;
	mDelVTet    = 0. ;
	mDelH       = 0. ;
	mDelVH      = 0. ;
    mDelVVess   = 0. ;
    mTimeTemp = 0.;
    memset(marrSinsSyst, 0, 3 * sizeof(double));

}

 // конструктор копирования
  QPeaceSins ::QPeaceSins (const QPeaceSins &R)
 {
	mSig_Q  =  R.mSig_Q ;
	mSig_Psi = R.mSig_Psi;
	mSig_Tet = R.mSig_Tet ;
	mSig_dQdt = R.mSig_dQdt ;
	mSig_dPsidt = R.mSig_dPsidt;
	mSig_dTetdt = R.mSig_dTetdt;
	mMaxSig_Q  =  R.mMaxSig_Q ;
	mMaxSig_Psi = R.mMaxSig_Psi;
	mMaxSig_Tet = R.mMaxSig_Tet ;
	mMaxSig_dQdt = R.mMaxSig_dQdt ;
	mMaxSig_dPsidt = R.mMaxSig_dPsidt;
	mMaxSig_dTetdt = R.mMaxSig_dTetdt;
	mMaxSig_H = R.mMaxSig_H ;
	mMaxSig_VH = R.mMaxSig_VH ;

	mSig_H = R.mSig_H;
	mSig_VH = R.mSig_VH ;
	mK1 = R.mK1 ;
	mSigV = R.mSigV ;
     mTCur = R. mTCur ;
	mDelQ = R.mDelQ ;
	mDelVQ = R.mDelVQ ;
	mDelPsi = R.mDelPsi ;
	mDelVPsi = R.mDelVPsi ;
	mDelTet = R.mDelTet ;
	mDelVTet = R.mDelVTet ;
	mDelH  =R.mDelH  ;
	mDelVH = R.mDelVH ;
	mDelVVess = R.mDelVVess ;
	mEstQ = R.mEstQ;
	mEstVQ = R.mEstVQ;
	mEstPsi = R.mEstPsi;
	mEstVPsi = R.mEstVPsi;
	mEstTet = R.mEstTet;
	mEstVTet =R.mEstVTet ;
	mEstH = R.mEstH ;
	mEstVH = R.mEstVH;
	mEstVVess = R.mEstVVess ;
    mTimeTemp = R.mTimeTemp;
    memcpy(marrSinsSyst, R.marrSinsSyst, 3 * sizeof(double));

 } 

 // оператор присваивания
 QPeaceSins &QPeaceSins::operator=(const QPeaceSins  &R)
 {
     if(this == &R)
     {
         return *this;
     }

	mSig_Q  =  R.mSig_Q ;
	mSig_Psi = R.mSig_Psi;
	mSig_Tet = R.mSig_Tet ;
	mSig_dQdt = R.mSig_dQdt ;
	mSig_dPsidt = R.mSig_dPsidt;
	mSig_dTetdt = R.mSig_dTetdt;
	mMaxSig_Q  =  R.mMaxSig_Q ;
	mMaxSig_Psi = R.mMaxSig_Psi;
	mMaxSig_Tet = R.mMaxSig_Tet ;
	mMaxSig_dQdt = R.mMaxSig_dQdt ;
	mMaxSig_dPsidt = R.mMaxSig_dPsidt;
	mMaxSig_dTetdt = R.mMaxSig_dTetdt;
	mMaxSig_H = R.mMaxSig_H ;
	mMaxSig_VH = R.mMaxSig_VH ;

	mSig_H = R.mSig_H;
	mSig_VH = R.mSig_VH ;
	mK1 = R.mK1 ;
	mSigV = R.mSigV ;
     mTCur = R. mTCur ;
	mDelQ = R.mDelQ ;
	mDelVQ = R.mDelVQ ;
	mDelPsi = R.mDelPsi ;
	mDelVPsi = R.mDelVPsi ;
	mDelTet = R.mDelTet ;
	mDelVTet = R.mDelVTet ;
	mDelH  =R.mDelH  ;
	mDelVH = R.mDelVH ;
	mDelVVess = R.mDelVVess ;
	mEstQ = R.mEstQ;
	mEstVQ = R.mEstVQ;
	mEstPsi = R.mEstPsi;
	mEstVPsi = R.mEstVPsi;
	mEstTet = R.mEstTet;
	mEstVTet =R.mEstVTet ;
	mEstH = R.mEstH ;
	mEstVH = R.mEstVH;
    mEstVVess = R.mEstVVess ;
    mTimeTemp = R.mTimeTemp;
    memcpy(marrSinsSyst, R.marrSinsSyst, 3 * sizeof(double));

    return *this ;
 }



   // парам конструктор  1
 QPeaceSins::QPeaceSins (const TEnvironment Environment,const double MaxSig_Q
                         , const double MaxSig_Psi
                         , const double MaxSig_Tet
            ,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt
            ,const double MaxSig_H,const double MaxSig_VH,const double K1
            ,const double SigV,const double  VAlT, const double TimeTemp, double *arrSinsSyst )


 {
	const double valCoeffEnv = ((double)Environment.mBallWave)/9. ;

	mMaxSig_Q  = MaxSig_Q ;
	mMaxSig_Psi = MaxSig_Psi ;
	mMaxSig_Tet = MaxSig_Tet;
	mMaxSig_dQdt = MaxSig_dQdt ;
	mMaxSig_dPsidt = MaxSig_dPsidt ;
	mMaxSig_dTetdt = MaxSig_dTetdt ;
	mMaxSig_H = MaxSig_H ;
	mMaxSig_VH = MaxSig_VH ;
	mSig_Q  =  MaxSig_Q * valCoeffEnv;
	mSig_Psi = MaxSig_Psi * valCoeffEnv;
	mSig_Tet = MaxSig_Tet * valCoeffEnv;
	mSig_dQdt =MaxSig_dQdt * valCoeffEnv;
	mSig_dPsidt = MaxSig_dPsidt * valCoeffEnv;
	mSig_dTetdt = MaxSig_dTetdt * valCoeffEnv;
	mSig_H =  mMaxSig_H * valCoeffEnv;
	mSig_VH =  mMaxSig_VH * valCoeffEnv;
	mK1 = K1 ;
	mSigV = SigV ;
    mTCur = VAlT;
    mTimeTemp = TimeTemp;
    memcpy(marrSinsSyst, arrSinsSyst, 3 * sizeof(double) );

}


     // парам конструктор  2
 QPeaceSins::QPeaceSins (const TEnvironment Environment
                     ,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet
                     ,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt
                     ,const double MaxSig_H,const double MaxSig_VH,const double K1
                     ,const double SigV
                     ,const double  VAlT      // время
                     ,const double VAlTrueVVess   // ист скорость
                     ,const double  VAlTrueQ  // истинный угол курса
                     ,const double  VAlTruePsi  // истинный угол килевой качки
                     ,const double  VAlTrueTet  // истинный угол бортовой качки
                     ,const double  VAlTrueVQ  //истинный скорость изменения угла курса
                     ,const double  VAlTrueVPsi   //истинный скорость изменения угла килевой качки
                     ,const double  VAlTrueVTet   //истинный вкорость изменения угла бортовой качуки
                     ,const double VAlTrueH
                     ,const double VAlTrueVH
                     ,const double VAlT_Q
                     ,const double  VAlT_Psi
                     ,const double VAlT_Tet
                     , double *arrDelt
                     , const double TimeTemp
                     , double *arrSinsSyst )


 {
	const double valCoeffEnv = ((double)Environment.mBallWave)/9. ;

	mMaxSig_Q  = MaxSig_Q ;
	mMaxSig_Psi = MaxSig_Psi ;
	mMaxSig_Tet = MaxSig_Tet;
	mMaxSig_dQdt = MaxSig_dQdt ;
	mMaxSig_dPsidt = MaxSig_dPsidt ;
	mMaxSig_dTetdt = MaxSig_dTetdt ;
	mMaxSig_H = MaxSig_H ;
	mMaxSig_VH = MaxSig_VH ;
	mSig_Q  =  MaxSig_Q * valCoeffEnv;
	mSig_Psi = MaxSig_Psi * valCoeffEnv;
	mSig_Tet = MaxSig_Tet * valCoeffEnv;
	mSig_dQdt =MaxSig_dQdt * valCoeffEnv;
	mSig_dPsidt = MaxSig_dPsidt * valCoeffEnv;
	mSig_dTetdt = MaxSig_dTetdt * valCoeffEnv;
	mSig_H =  mMaxSig_H * valCoeffEnv;
	mSig_VH =  mMaxSig_VH * valCoeffEnv;
	mK1 = K1 ;
	mSigV = SigV ;
    mTimeTemp = TimeTemp;
    memcpy(marrSinsSyst, arrSinsSyst, 3 * sizeof(double) );

 fillValues_Delts_and_Ests(VAlT,  VAlTrueVVess,	 VAlTrueQ, VAlTruePsi,	 VAlTrueTet, VAlTrueVQ
	, VAlTrueVPsi,  VAlTrueVTet, VAlTrueH, VAlTrueVH, VAlT_Q,  VAlT_Psi, VAlT_Tet, arrDelt) ;


}
// моделирование в соответствии с лабораторной моделью
 void QPeaceSins::recalcPeaceSins_v0(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt)
{
     int n = (int)((valT - mTCur +1.E-20)/mTimeTemp);
     if ( n < 1)
     {
         return;
     }
     double valTCur = mTCur + ((double)n)* mTimeTemp;
     fillValues_Delts_and_Ests( valTCur,  VVess,	 Q, Psi,	 Tet, VQ
	, VPsi,  VTet, H, VH, T_Q,  T_Psi, T_Tet, arrDelt) ;



}

//формирование ошибок измерения и измерений СИНС
 void QPeaceSins::fillValues_Delts_and_Ests(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt)
{

    mDelQ = mSig_Q * sin(2.*M_PI/T_Q * valT - arrDelt[0]) ;
  // 11.ошибка определения скорости ищзменения угла курса:
  mDelVQ = 2.*M_PI/T_Q * mSig_Q * cos(2.*M_PI/T_Q * valT - arrDelt[0]);
  // 12.ошибка определения угла килевой качки :
  mDelPsi = mSig_Psi * sin(2.*M_PI/T_Psi * valT - arrDelt[1]);
  // 13.ошибка определения скорости ищзменения угла килевой качки:
  mDelVPsi = 2.*M_PI/T_Psi * mSig_Psi * cos(2.*M_PI/T_Psi * valT - arrDelt[1]);
  // 14.ошибка определения угла боротовой качки :
  mDelTet = mSig_Tet * sin(2.*M_PI/T_Tet * valT - arrDelt[2]);
  // 15. ошибка определения скорости ищзменения угла боротовой качки:
  mDelVTet = 2.*M_PI/T_Tet * mSig_Tet * cos(2.*M_PI/T_Tet * valT - arrDelt[2]);
  // 16.ошибка измерения высоты цетнтра корабля :
  mDelH = mSig_H * sin(2. * M_PI/ T_Tet * valT - arrDelt[3]) ;
	// 17.ошибка измерения скорости изменения  высоты цетнтра корабля :
  mDelVH = 2. * M_PI/ T_Tet * mSig_H * cos(2. * M_PI/ T_Tet *valT - arrDelt[3]) ;
  // 18.ошибка измерения скорости корабля:
  mDelVVess = VVess * mK1;

// ОЦЕНКИ (измерения) СИНС
   // QUANT_COLS_BUFF_PeaceSins.оценка ооценка угла курса :

   mEstQ = Q + mDelQ;  // ТЕСТ
  // 20.оценкая скорости ищзменения угла курса:
   mEstVQ = VQ + mDelVQ ;
  // 21.оценка угла килевой качки :
   mEstPsi = Psi +mDelPsi; // ТЕСТ
  // 22.оценка скорости ищзменения угла килевой качки:
   mEstVPsi = VPsi + mDelVPsi;
  // 23.оценка угла боротовой качки :
   mEstTet = Tet + mDelTet; // ТЕСТ
  // 24. оценка скорости ищзменения угла боротовой качки:
   mEstVTet = VTet + mDelVTet ;
  // 25.оценка высоты цетнтра корабля :
   mEstH = H + mDelH;
  // 26.оценка скорости изменения  высоты цетнтра корабля :
   mEstVH = VH + mDelVH;
  // 27.оценка скорости корабля:
   mEstVVess = VVess + mDelVVess;
   mTCur =  valT;
}
// Усовершенствованное моделирование  на основе системы 1-го порядка,
// позволяющее правильно имитировать ошибки оценивания скоростей изенения углов качек и курса
 void QPeaceSins::recalcPeaceSins_v1(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt)
{
     int n = (int)((valT - mTCur + 1.E-20)/mTimeTemp);
     if ( n < 1)
     {
         return;
     }
     double valTCur = mTCur + ((double)n)* mTimeTemp;
  const double h = valTCur -  mTCur ;
  // ошибка определения угла курса
  double valDelQ1 = callStepSys1(mDelQ, h, mSig_Q, mSig_dQdt) ;
  mDelVQ = valDelQ1 - mDelQ ;
  mDelQ =  valDelQ1 + marrSinsSyst [0];

  // 12.ошибки определения угла килевой качки :
  double valDelPsi1 = callStepSys1(mDelPsi, h, mSig_Psi, mSig_dPsidt) ;
  mDelVPsi = valDelPsi1 - mDelPsi ;
  mDelPsi =  valDelPsi1 + marrSinsSyst [1];

  // 14.ошибка определения угла боротовой качки :
  double valDelTet1 = callStepSys1(mDelTet, h, mSig_Tet, mSig_dTetdt) ;
  mDelVTet = valDelTet1 - mDelTet ;
  mDelTet =  valDelTet1+ marrSinsSyst [2];


  // 16.ошибка измерения высоты цетнтра корабля :
  mDelH = mSig_H * sin(2. * M_PI/ T_Tet * valT - arrDelt[3]) ;
  // 17.ошибка измерения скорости изменения  высоты цетнтра корабля :
  mDelVH = 2. * M_PI/ T_Tet * mSig_H * cos(2. * M_PI/ T_Tet *valT - arrDelt[3]) ;
  // 18.ошибка измерения скорости корабля:
  mDelVVess = VVess * mK1;

// ОЦЕНКИ (измерения) СИНС
   // QUANT_COLS_BUFF_PeaceSins.оценка ооценка угла курса :
   mEstQ = Q + mDelQ;
  // 20.оценкая скорости ищзменения угла курса:
   mEstVQ = VQ + mDelVQ ;
  // 21.оценка угла килевой качки :
   mEstPsi = Psi +mDelPsi;
  // 22.оценка скорости ищзменения угла килевой качки:
   mEstVPsi = VPsi + mDelVPsi;
  // 23.оценка угла боротовой качки :
   mEstTet = Tet + mDelTet;
  // 24. оценка скорости ищзменения угла боротовой качки:
   mEstVTet = VTet + mDelVTet ;
  // 25.оценка высоты цетнтра корабля :
   mEstH = H + mDelH;
  // 26.оценка скорости изменения  высоты цетнтра корабля :
   mEstVH = VH + mDelVH;
  // 27.оценка скорости корабля:
   mEstVVess = VVess + mDelVVess;

   // коррекция времени
    mTCur = valTCur;
}

 // расчет динамич системы 1-го порядка
 // INPUT:
 // valX- переменная. valX (t+h)= a * valX(t) + ksi
 //  h - шаг по времени
 //  valK - дисперсия valX
 //  valSigV*valSigV*h  - дисперсия разности valX[ t+h}- valX[ t+h]
 //
 // OUTPUT
 // valX (t + h)
 //
 //
 double QPeaceSins::callStepSys1(const double valX, const double h, const double SigX
    ,const double valSigV)
 {
   if (fabs(SigX ) < 0.00000001) return 0 ;

   const double sig2V =  h * valSigV * valSigV ;
   const double sig2X =  SigX * SigX ;
   const double a = 1. - sig2V / sig2X/ 2. ;
   if (( 1. - sig2V/sig2X/ 4.) > 0.)
   {
   const double sig2Ksi =  sig2V * ( 1. - sig2V/sig2X/ 4.) ;
     double ksi = getGauss(0, sqrt(sig2Ksi) );
   return a * valX + ksi ;
   }
   else
   {
     return getGauss(0, SigX );
   }
}
 //------------------------------
// формирование вектора оценок углов Эйлера в ПСК
void QPeaceSins::getEstArrEilers (const double VAlTCur,double *arrMu)
{
    double dt = VAlTCur - mTCur;
    arrMu[0] = mEstQ + dt *mEstVQ;
    arrMu[1] = mEstPsi + dt *mEstVPsi;
    arrMu[2] = mEstTet + dt *mEstVTet;
}
//------------------------------
// формирование вектора оценок угловой скорости вращения корабля в неподвижной ПСК
void QPeaceSins::getEstArr_dEilers_po_dt (const double VAlTCur,double *arr_dEilers_po_dt)
{
    arr_dEilers_po_dt[0] = mEstVPsi;
    arr_dEilers_po_dt[1] = mEstVTet;
    arr_dEilers_po_dt[2] = mEstVQ;
}

//---------------------------------------------------------
// формирование сводного вектора оценок угловой информации
// INPUT:
// VAlTCur - время на которое надо выдать информацию
// OUTPUT:
// arrCurSinsInfo[6] - сводный вектор
//arrCurSinsInfo[0], arrCurSinsInfo[1] arrCurSinsInfo[2] -mEstQ, mEstPsi, mEstTet
//arrCurSinsInfo[3], arrCurSinsInfo[4] arrCurSinsInfo[5] -mEstVQ, mEstVPsi, mEstVTet
void QPeaceSins::getCurVectInfo(const double VAlTCur,double *arrCurSinsInfo)
{
    getEstArrEilers (VAlTCur,arrCurSinsInfo);
    getEstArr_dEilers_po_dt (VAlTCur,&arrCurSinsInfo[3]);
}
//---------------------------------------------------

int QPeaceSins::createInputDataReport(wchar_t*FileName, const bool bHeader)

{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
    {
      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {
     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"***************************************\n");
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
 fprintf(fw,"***************************************\n");
 fprintf(fw,"      Инерциальная система \n");

fprintf(fw,"  максимальная СКО ошибки КУ (рад) = %6.5f\n",mMaxSig_Q);
fprintf(fw,"  максимальная СКО ошибки по углу КК (рад) = %6.5f\n",mMaxSig_Psi);
fprintf(fw,"  максимальная СКО СКО ошибки по углу БК (рад) %6.5f\n",mMaxSig_Tet);
fprintf(fw,"  максимальная СКО  скорости КУ (рад) = %6.5f\n",mMaxSig_dQdt);
fprintf(fw,"  максимальная СКО скорости угла КК (рад) = %6.5f\n",mMaxSig_dPsidt);
fprintf(fw,"  максимальная СКО  скорости угла БК (рад) = %6.5f\n",mMaxSig_dTetdt);
fprintf(fw,"  максимальная СКО высоты (м) = %4.2f\n",mMaxSig_H);
fprintf(fw,"  максимальная СКО скорости высоты (м) = %4.2f\n",mMaxSig_VH);
fprintf(fw,"  коэффициент ошибки определения скорости = %5.4f\n",mK1);

 fclose(fw);

}



