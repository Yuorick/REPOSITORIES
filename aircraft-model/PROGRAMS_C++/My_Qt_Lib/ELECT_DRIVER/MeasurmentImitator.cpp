#include "MeasurmentImitator.h"

#include "Gauss.h"
#include <math.h>
#include <string.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>

extern bool BEZ_SHUMOV;



  QDriverMeasure:: QDriverMeasure()
  {
      memset(marrYzv, 0, 3 * sizeof(double));
      memset(marrKYzv, 0, 9 * sizeof(double));

  }

  //---------------------------------------------------------------------------

  // конструктор копирования
    QDriverMeasure :: QDriverMeasure (const  QDriverMeasure &R)
   {
        memcpy(marrYzv, R.marrYzv, 3 * sizeof(double));
        memcpy(marrKYzv, R.marrKYzv, 9 * sizeof(double));
    }

   // оператор присваивания
    QDriverMeasure  &QDriverMeasure::operator=( const QDriverMeasure  &R)
   {
        memcpy(marrYzv, R.marrYzv, 3 * sizeof(double));
        memcpy(marrKYzv, R.marrKYzv, 9 * sizeof(double));
       return *this ;
   }


    // парам конструктор
    QDriverMeasure:: QDriverMeasure (const double  *arrYzv,const double  *arrKYzv)
   {
        memcpy(marrYzv, arrYzv, 3 * sizeof(double));
        memcpy(marrKYzv, arrKYzv, 9 * sizeof(double));
   }


    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
  QMeasurmentImitator::QMeasurmentImitator()
  {
      // погрешность измерения тока Id, %
       mErrPerCent_Id = 0.;
      // погрешность измерения тока Iqu, %
       mErrPerCent_Iqu =0.;
      // СКЗ ошибки измерения кглового положения
       mSgmTetta =0.;
  }  
  //---------------------------------------------------------------------------

  // конструктор копирования
    QMeasurmentImitator :: QMeasurmentImitator (const  QMeasurmentImitator &R)
   {
           mErrPerCent_Id = R.mErrPerCent_Id;
           mErrPerCent_Iqu = R.mErrPerCent_Iqu ;
           mSgmTetta = R.mSgmTetta;
   }

   // оператор присваивания
    QMeasurmentImitator  &QMeasurmentImitator::operator=( const QMeasurmentImitator  &R)
   {
        mErrPerCent_Id = R.mErrPerCent_Id;
        mErrPerCent_Iqu = R.mErrPerCent_Iqu ;
        mSgmTetta = R.mSgmTetta;
       return *this ;
   }


    // парам конструктор
    QMeasurmentImitator:: QMeasurmentImitator (const double  ErrPerCent_Id , const double  ErrPerCent_Iqu,
        const double  SgmTetta)
   {
      mErrPerCent_Id = ErrPerCent_Id;
      mErrPerCent_Iqu  = ErrPerCent_Iqu ;
      mSgmTetta  = SgmTetta ;
   }



   // имитатор измерений
   // INPUT:
   // arrPfVect[4] - фазовый вектор сиситемы
   //OUTPUT:
   // arrMeasuredfVect[3] - вектор измерений - Id, Iqu, Tetta
   // arrKy[9] - еляционная матрица ошибок измерения
  void QMeasurmentImitator:: createMeasure(double *arrPfVect,  QDriverMeasure *measureCur)
  {
      measureCur->marrYzv[0] = arrPfVect[1] + getGauss(0., mErrPerCent_Id *fabs(arrPfVect[1])/100./3. );
      measureCur->marrYzv[1] =   arrPfVect[2] + getGauss(0., mErrPerCent_Iqu *fabs(arrPfVect[2])/100./3. );
      measureCur->marrYzv[2] = arrPfVect[3] + getGauss(0., mSgmTetta );
      memset(measureCur->marrKYzv, 0., 9 * sizeof(double));
      measureCur->marrKYzv[0] = (mErrPerCent_Id *fabs(arrPfVect[1])/100./3. +0.00001 )
              * (mErrPerCent_Id *fabs(arrPfVect[1])/100./3. +0.00001 );
      measureCur->marrKYzv[4] = (mErrPerCent_Iqu *fabs(arrPfVect[2])/100./3.  +0.00001)
              * (mErrPerCent_Iqu *fabs(arrPfVect[2])/100./3.  +0.00001);
      measureCur->marrKYzv[8] = mSgmTetta * mSgmTetta;
  }

  //--------------------------------------------
  int QMeasurmentImitator::createInputDataReport(wchar_t*FileName, const bool bHeader)
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
  fprintf(fw," Характеристики измерит. приборов приводов: \n");
  fprintf(fw," погрешность измерения тока Id(%) =  %5.4f\n",mErrPerCent_Id);
  fprintf(fw," погрешность измерения тока Iq(%) =  %5.4f\n",mErrPerCent_Iqu);
  fprintf(fw," СКЗ ошибки измерения углового положения (рад) =  %8.7f\n",mSgmTetta);
  fclose(fw);


  }



