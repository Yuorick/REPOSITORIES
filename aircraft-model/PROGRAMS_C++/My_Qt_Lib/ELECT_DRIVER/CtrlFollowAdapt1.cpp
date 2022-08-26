#include "CtrlFollowAdapt1.h"
#include "MatrixProccess.h"
#include <wchar.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
//--------------------------------------------
QCtrlFollowAdapt1::QCtrlFollowAdapt1():QCtrlFollow()
{
 mFiltrW = QFiltr();
 mDispWTochka = 0.;
}

// конструктор копирования
 QCtrlFollowAdapt1::QCtrlFollowAdapt1(const QCtrlFollowAdapt1 &R):QCtrlFollow( R)
 {
   mFiltrW = R.mFiltrW;
   mDispWTochka = R.mDispWTochka;
 }

 // оператор присваивания
  QCtrlFollowAdapt1  &QCtrlFollowAdapt1::operator=( const QCtrlFollowAdapt1  &R)
 {
      if(this == &R)
      {
          return *this;
      }
     QCtrlFollow:: operator= (R);
     mFiltrW = R.mFiltrW;
     mDispWTochka = R.mDispWTochka;

     return *this ;
 }
//-----------------------------------

// парам конструктор
  // параметрический конструктор
QCtrlFollowAdapt1:: QCtrlFollowAdapt1( const QElectMotor ElectMotor ,  QLoad*  Load
                           , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                           ,const double VAlOmegaBegin , const double T0, const double TCur
                           , const double valh,TComp *pCmpArrLamb, const double *ARrFiltK_Begin
                                       ,const double DispWTochka)
    :QCtrlFollow (ElectMotor ,  Load, arrSpreadParams, MomOut, VAlTettaBegin
               ,VAlOmegaBegin , T0, TCur, valh,pCmpArrLamb)
  {

    double arrPhVectBegin[3] = {0.}, parrC[3] = {1.,0.,0.};
    arrPhVectBegin[0] = VAlTettaBegin;
    arrPhVectBegin[1] = VAlOmegaBegin;
    mFiltrW = QFiltr ( 3, 1,arrPhVectBegin
      , ARrFiltK_Begin,  TCur,  valh, true, parrC);
    mDispWTochka = DispWTochka;
}

//------------------------------------------------------------------
void QCtrlFollowAdapt1:: correctMomOut(const double VAlTettaZv
                          ,const double VAlDispTetta    , const double VAlW, const double VAlTcur)
{
  // экстраполяция фазового вектора на момент VAlTcur
    double valH = VAlTcur - mFiltrW.mTCur;
    // фундаментальная матрица (матрица перехода)
    double arrL[9] = {0.};
    arrL[0] = 1.;
    arrL[1] = valH;
    arrL[2] = valH * valH /2.;
    arrL[4] = 1.;
    arrL[5] = valH;
    arrL[8] = 1.;
    //
    double arrT0[3] = {0.};
    MtrxMultMatrx(arrL,3, 3,  mFiltrW.marrCurEst,1, arrT0);
    mFiltrW.marrCurEst[0] = arrT0[0] + VAlW *valH* valH /2.;
    mFiltrW.marrCurEst[1] = arrT0[1] + VAlW * valH ;
    mFiltrW.marrCurEst[2] = arrT0[2];

    double valTettaZv =VAlTettaZv;
    double valDispTetta = VAlDispTetta;
    double arrKww[9] = {0.};
            arrKww[0] = valH * valH * valH * valH * valH /5.;
            arrKww[1] = arrKww[3] = valH * valH * valH * valH /4.;
            arrKww[2] = arrKww[6] = arrKww[4] = valH * valH * valH /3.;
            arrKww[5] = arrKww[7] = valH * valH /2.;
            arrKww[8] = valH;
            MatrxMultScalar(arrKww, 3, 3, mDispWTochka,arrKww);
    mFiltrW.processMeasure_Type2(&valTettaZv,&valDispTetta,  VAlTcur
                                ,arrL,arrKww);
   // double valMomOut = mFiltrW.marrCurEst[2]/getInvSumJ0();
   // setParam(6, valMomOut);
}
//------------------------------------------------------------------
double QCtrlFollowAdapt1::estimateMomOut(const double VAlTettaZv
     ,const double VAlDispTetta    , const double VAlW, const double VAlTcur)
{
  // экстраполяция фазового вектора на момент VAlTcur
    double valH = VAlTcur - mFiltrW.mTCur;
    // фундаментальная матрица (матрица перехода)
    double arrL[9] = {0.};
    arrL[0] = 1.;
    arrL[1] = valH;
    arrL[2] = valH * valH /2.;
    arrL[4] = 1.;
    arrL[5] = valH;
    arrL[8] = 1.;
    //
    double arrT0[3] = {0.};
    MtrxMultMatrx(arrL,3, 3,  mFiltrW.marrCurEst,1, arrT0);
    mFiltrW.marrCurEst[0] = arrT0[0] + VAlW *valH* valH /2.;
    mFiltrW.marrCurEst[1] = arrT0[1] + VAlW * valH ;
    mFiltrW.marrCurEst[2] = arrT0[2];

    double valTettaZv =VAlTettaZv;
    double valDispTetta = VAlDispTetta;
    double arrKww[9] = {0.};
            arrKww[0] = valH * valH * valH * valH * valH /5.;
            arrKww[1] = arrKww[3] = valH * valH * valH * valH /4.;
            arrKww[2] = arrKww[6] = arrKww[4] = valH * valH * valH /3.;
            arrKww[5] = arrKww[7] = valH * valH /2.;
            arrKww[8] = valH;
            MatrxMultScalar(arrKww, 3, 3, mDispWTochka,arrKww);
    mFiltrW.processMeasure_Type2(&valTettaZv,&valDispTetta,  VAlTcur
                                ,arrL,arrKww);
    return mFiltrW.marrCurEst[2]/getInvSumJ0();

}
//----------------------------------------
double QCtrlFollowAdapt1::calcEstMomOut()
{
    double valMomOut = mFiltrW.marrCurEst[2]/getInvSumJ0();
    return valMomOut;
}
//-------------------------------------------------------
int QCtrlFollowAdapt1::createInputDataReport_CtrlFollow_Inherited(wchar_t*FileName, const bool bHeader)
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
fprintf(fw,"***************************************\n");
fprintf(fw,"  Подтип системы управления - \n");
fprintf(fw,"  следящая с адаптацией по внешнему моменту \n");
fprintf(fw,"  СКЗ априорной оценки скорости ускорения цели (м/(с*с*с)) =  %5.4\n",sqrt(mDispWTochka));
fclose(fw);
}

