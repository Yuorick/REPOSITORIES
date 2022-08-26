#include "Load.h"
#include <math.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>

#include "MatrixProccess.h"
#include "Segment.h"
#include "Environment.h"
#include "Plane.h"

extern const double VAL_TEST_COmx = 0.;
extern const double VAL_TEST_J = 0.1;
extern const double VAL_TEST_Cv = 0.01;
extern const double VAL_TEST_dOm_po_dt = 10000.;
extern const double VAL_TEST_Cx = 0.;
//---------------------------------------------------------------------------

   QLoad:: QLoad()
{
    //момент инерции нагрузки
    mJPayLoad = 0.;
    // коэффиц сопротивления воздуха профиля
    mCx = 0.;
     // коэффиц трения воздуха
    mCv = 0.;
    //
    mMax_dOm_po_dt = 0.;
}

//---------------------------------------------------------------------------

// конструктор копирования
     QLoad :: QLoad (const  QLoad &R)
 {
   mJPayLoad  = R.mJPayLoad ;
   mCx = R.mCx;
   mCv = R.mCv;
   mMax_dOm_po_dt = R.mMax_dOm_po_dt;

 }
     //-------------------------------------------------------------------

 // оператор присваивания
  QLoad  &QLoad::operator=( const QLoad  &R)
 {
      mJPayLoad  = R.mJPayLoad ;
      mCx = R.mCx;
      mCv = R.mCv;
      mMax_dOm_po_dt = R.mMax_dOm_po_dt;

     return *this ;
 }

  //-------------------------------------------------------------------

  // парам конструктор 1
 QLoad:: QLoad ( const double  JRotor, const double  COmx)
 {
    mJPayLoad  = JRotor;
    mCx = COmx;
    mCv = 0.;
    mMax_dOm_po_dt = 0.;

 }

//----------------------------------------------
   // парам конструктор 2
 QLoad:: QLoad ( const double  JRotor, const double  Cx, const double  Cv,const double  Max_dOm_po_dt)
  {
     mJPayLoad  = JRotor;
     mCx = Cx;
     mCv = Cv;
     mMax_dOm_po_dt = Max_dOm_po_dt;
  }



//-------------------------------------------------------------------
void  QLoad::setLength(const double  VAlLength)
{

}
//-------------------------------------------------------------------

double  QLoad::getLength()
{
}
//-------------------------------------------------------------------

double  QLoad::calcAirResistMom(  TEnvironment &Environment,const double ValTet, const double VAlOm)
{
     double temp = -mCv*VAlOm - mCx*VAlOm* VAlOm * sign_(VAlOm) ;
    return temp;
}

//-----------------------------------------------------------
double QLoad::calc_dAirResistMom_po_dOmega( TEnvironment &Environment,const double ValTet, const double VAlOm)
{
    return -2.* mCx * VAlOm - mCv;
}
//--------------------------------------------------
double QLoad::calc_dMa_po_Cx(const double VAlOm)
{
 return - VAlOm  * VAlOm  *  sign_ (VAlOm);
}
//-----------------------------------------------------
double QLoad::calc_dMa_po_Cv(const double VAlOm)
{
 return - VAlOm  ;
}
//-----------------------------------------------------
double QLoad::calc_d2Ma_po_dOm_po_Cv(const double VAlOm)
{
 return - 1.  ;
}
//--------------------------------------------------
double QLoad::calc_d2Ma_po_dOm_po_Cx(const double VAlOm)
{
 return -2.* VAlOm  ;
}
//---------------------------------
int QLoad::createInputDataReport(wchar_t*FileName, const bool bHeader)

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
fprintf(fw,"      Нагрузка\n");
fprintf(fw,"      тип - абстрактная\n");
fprintf(fw,"      Параметры   нагрузки\n");

fprintf(fw,"  момент инерции нагрузки(кг*м*м) = %5.2f\n",mJPayLoad);
fprintf(fw,"  коэффиц сопротивления воздуха = %5.2f\n",mCx);
fprintf(fw,"  коэффиц трения воздуха = %5.2f\n",mCv);

 fclose(fw);

}


