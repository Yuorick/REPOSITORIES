#include "HidroRLS.h"
#include <math.h>
#include <string.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "Gauss.h"
#include "TrueMeasParams.h"
#include "BigMeasure.h"



QHidroRLS::QHidroRLS()
{
    mSigT = 0.0001 ;
    mDiagWidth = 1. * M_PI/3000 ;
}


// конструктор копирования
 QHidroRLS ::QHidroRLS (const QHidroRLS &R)
 {
    mSigT = R.mSigT ;
	mDiagWidth = R.mDiagWidth ;
 }
 // оператор присваивания
 QHidroRLS &QHidroRLS::operator=(const QHidroRLS  &R)
 {
     if(this == &R)
     {
         return *this;
     }
     mSigT = R.mSigT ;
     mDiagWidth = R.mDiagWidth ;
     return *this ;
}


  // парам конструктор1
 QHidroRLS::QHidroRLS ( const double SigR,const double DiagWidth)
 {	
    mSigT = SigR;
	mDiagWidth = DiagWidth ;
 }


 // формирование измерения и ошибок измерения
 // имитация измерения антенны
 // INPUT:
 //trueMeasParams - структура истинных параметров, характеризующих положение цели
 // и обеспечивающих возможность наложения шумов и погрешностей
 //OUTPUT:
 //Meas - структура, в которой хранятся данные одного измерения
 // функция вычисляет следующие члены структуры Meas:
 // Meas.mTzapr - момент запросного сигнала
 // Meas.mTotv - момент ответного сигнала
 // Meas.mq - угловое измерение КУ
 // Meas.me- угловое измерение КУ

 void QHidroRLS::imitateMesure(const QTrueMeasParams trueMeasParams, QBigMeasure &Meas)
 {
    Meas.mTzaprZv = trueMeasParams.mTzapr ;
    Meas.mSig_t = calcSig_t(trueMeasParams);
    // ТЕСТ
    Meas.mTotvZv = trueMeasParams.mTotv + calcDeltaT(trueMeasParams);

    calc_qZv_and_eZv( trueMeasParams, Meas.mqzv, Meas.mezv
                      , Meas.mSig_q, Meas.mSig_e);

 }
 //----------------------------------------
 void QHidroRLS::imitateMesure(const QTrueMeasParams trueMeasParams, double *valTzaprZv
                               , double *valTotvZv, double *valSig_t, double *val_qzv
                               , double *val_ezv, double *valSig_q, double *valSig_e)
 {
    *valTzaprZv = trueMeasParams.mTzapr ;
    *valTotvZv = trueMeasParams.mTotv + calcDeltaT(trueMeasParams);
    *valSig_t = calcSig_t(trueMeasParams);
    calc_qZv_and_eZv( trueMeasParams, *val_qzv, *val_ezv
                      , *valSig_q, *valSig_e);

 }
 //-------------------------------------------------
 //
 double QHidroRLS::calcSig_t(const QTrueMeasParams trueMeasParams)
 {
    return mSigT;
 }
 // вычисление ошибки измерения времени прихода сигнала
 double QHidroRLS::calcDeltaT(const QTrueMeasParams trueMeasParams)
 {
     return getGauss(0, mSigT );
 }

 //-------------------------------------------------
 // вычисление  измерения КУ
 void QHidroRLS::calc_qZv_and_eZv(const QTrueMeasParams trueMeasParams, double &qZv, double &eZv
                                  , double &Sig_q, double &Sig_e)
 {
     qZv =  NODATA;
     eZv = NODATA;
     Sig_q =  NODATA;
     Sig_e =  NODATA;
 }

 //-----------------------------
 int QHidroRLS::createInputDataReport(wchar_t*FileName, const bool bHeader)
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

 fprintf(fw,"  тип устройства - ОЭК\n");
 fprintf(fw,"  характеристики ОЭК:\n");



 fprintf(fw,"  СКЗ ошибки по дальности(м) = %6.5f\n",mSigT);
 fprintf(fw,"  ширина диаграммы (рад) = %6.5f\n",mDiagWidth); 

  fclose(fw);

 }






