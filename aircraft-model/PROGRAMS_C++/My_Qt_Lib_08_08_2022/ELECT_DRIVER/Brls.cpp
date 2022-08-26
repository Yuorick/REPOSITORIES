#include <math.h>
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "Brls.h"
#include "Segment.h"
#include "Environment.h"
#include "Plane.h"
#include "MatrixProccess.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

extern const double VAL_BRLS_LENGTH = 2.54;
extern const double VAL_BRLS_HEIGHT = 0.16;
extern const double VAL_BRLS_MASS = 3.55;
extern const double VAL_BRLS_Cv = 0.;//0.001;
extern const double VAL_BRLS_dOm_po_dt = 2.* M_PI * 20./ 60./3.;
extern const double VAL_BRLS_Cx = 0.8;



QBrls::QBrls():QLoad()
{
 mLength = 0.;
 mbBarrier = false;
}

// конструктор копирования
  QBrls :: QBrls (const  QBrls &R):QLoad( R)
 {     
   mLength = R.mLength;
   mHeight = R.mHeight;
   mSgmBarrier = R.mSgmBarrier;
   mPlane = R.mPlane;
   mbBarrier = R.mbBarrier;
 }

 // оператор присваивания
  QBrls  &QBrls::operator=( const QBrls  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      QLoad:: operator= (R);

      mLength = R.mLength;
      mHeight = R.mHeight;
      mSgmBarrier = R.mSgmBarrier;
      mPlane = R.mPlane;
      mbBarrier = R.mbBarrier;

     return *this ;
 }

  // парам конструктор 1
  // создание БРЛС
 QBrls:: QBrls ( const double  Cv, const double   Cx, const double  VAlMass, const double  VAlLengthTotal
                 ,const double  VAlHeight
                 ,const double  Max_dOm_po_dt, const QSegment SgmBarrier, const TPlane Plane)
  :QLoad (0.,  Cx,   Cv, Max_dOm_po_dt)
 {
    mLength = VAlLengthTotal;
    mHeight = VAlHeight;
    mJPayLoad  = VAlMass * VAlLengthTotal*VAlLengthTotal/ 12.;

    mSgmBarrier = SgmBarrier;
    mPlane = Plane;
    mbBarrier = false;
 }

  void  QBrls::setLength(const double  VAlLength)
  {
      mLength = VAlLength;
  }

  double  QBrls::getLength()
  {
     return mLength ;
  }

  //----------------------------------------------
  // вычисление момента сопротивления воздуха по угловой скорости
  // INPUT:
  //Environment - внешняя среда - скорость и направления ветра в ГСК
  // Plane - плоскость вращения БРЛС в ГСК
  // SgmBarrier - сегмент барьера, препятствующего ветру в плоскости вращения БРЛС
  // ValTet - угловое положения БРЛС
  // VAlOm - угловая скорость вращения
  double  QBrls::calcAirResistMom( TEnvironment &Environment,const double ValTet, const double VAlOm)
  {

      double valMom0 =   calcAirResistFuncForBlade( Environment, ValTet,  VAlOm,calcIntegralAirMom);
      double valMom1 =   calcAirResistFuncForBlade( Environment, ValTet + M_PI,  VAlOm,calcIntegralAirMom);

      return valMom0 + valMom1-mCv*VAlOm;
  }
  //------------------------------------------------------------------------------
  //----------------------------------------------
  // вычисление производной момента сопротивления воздуха по угловой скорости
  // INPUT:
  //Environment - внешняя среда - скорость и направления ветра в ГСК
  // Plane - плоскость вращения БРЛС в ГСК
  // SgmBarrier - сегмент барьера, препятствующего ветру в плоскости вращения БРЛС
  // ValTet - угловое положения БРЛС
  // VAlOm - угловая скорость вращения
  double  QBrls::calc_dAirResistMom_po_dOmega( TEnvironment &Environment,const double ValTet, const double VAlOm)
  {
      double valMom0 =   calcAirResistFuncForBlade( Environment, ValTet,  VAlOm, calcIntegral_dAirMom_po_dOmega);
      double valMom1 =   calcAirResistFuncForBlade( Environment, ValTet + M_PI,  VAlOm, calcIntegral_dAirMom_po_dOmega);

      return valMom0 + valMom1 -mCv;
  }
  //------------------------------------------------------------------------------

  // вычисление коэффициента момента сопротивления воздуха вращающейся штанги без ветра
double  QBrls::calcCOmegaX ()
{
 return mCx * ATM_RoN0 * mHeight * mLength * mLength * mLength * mLength/ 64.;
}

//--------------------------------------------------------------------------------------

//вычисление момента сопротивления воздуха или производной момента для половины БРЛС (лопасти)
// INPUT:
//Environment - внешняя среда - скорость и направления ветра в ГСК
// Plane - плоскость вращения БРЛС в ГСК
// SgmBarrier - сегмент барьера, препятствующего ветру в плоскости вращения БРЛС
// ValTet - угловое положения БРЛС
// VAlOm - угловая скорость вращения
// pf - прототип функции или для момента или для его производной
double  QBrls::calcAirResistFuncForBlade(const TEnvironment Environment,const double ValTet, const double VAlOm
   ,double (*pf)(const double ,const double ,const double ,const double  ,const double,const double ))
{
    // 1. вектор скорости ветра в ГСК
    double arrWindV[3] = {0.};
    TEnvironment EnvironmentCur = Environment;
    EnvironmentCur.createVectWindV(arrWindV);
    ///

    // 2. перевод вектора скорости ветра из ГСК в СК
    double arrWindSKBrlsV[3] = {0.};
    MtrxTranspMultMatrx(mPlane.marrF,3, 3, arrWindV,1, arrWindSKBrlsV) ;
    ///

    // 3. перевод вектора скорости воздуха в систему координат лопасти
    double arrWindSKBrlsV0[2] = {0.};
    arrWindSKBrlsV0[0] = cos (ValTet) *arrWindSKBrlsV[0] + sin (ValTet) *arrWindSKBrlsV[1];
    arrWindSKBrlsV0[1] = -sin (ValTet) *arrWindSKBrlsV[0] + cos (ValTet) *arrWindSKBrlsV[1];
    ///

    // 4. вертикальная (нормальная к лопасти) составляющая скорости вектра
    double valW = arrWindSKBrlsV0[1];
    ///

    // скорость набегающего потока в точке лопасти с координатой x будет равна valW - VAlOm*x

    // определение незатененного участка лопасти
    // 5. перкеврд отрезка преграды в систему координат лопасти
    TURPointXY pntCentre(0., 0.);
    const TURPointXY pntSdvig(0.,0.);
    const double valRastigenie = 1.;
   TURPolyLine pln1 = mSgmBarrier.LinTransform(-ValTet, pntCentre //,??????????????????????????
                        , valRastigenie );
   QSegment SgmBarr0(pln1);
   ///

   // 6. сегмент лопасти
   TURPointXY pnt0(0.,0.);
   TURPointXY pnt1(mLength/ 2.,0.);
   QSegment BladeSgm(pnt0, pnt1);
   ///

   // нахождение открытых (незатенных) сегментов
   // их может быть 0,1 или 2
   QSegment arrOpenSegm[2],shadowSegm;
   int iquantOpenSgm = -1, iquantShadowSgm = -1;
   if (mbBarrier)
   {
   QSegment::calcOpen_and_Shadow_Segments(BladeSgm,  SgmBarr0,  arrWindSKBrlsV0
                       , arrOpenSegm, &iquantOpenSgm, &shadowSegm, &iquantShadowSgm);
   }
   else
   {
     iquantOpenSgm = 1;
     iquantShadowSgm = 0;
     arrOpenSegm[0].Points[0].X = 0.;
     arrOpenSegm[0].Points[0].Y = 0.;
     arrOpenSegm[0].Points[1].X = mLength/2.;
     arrOpenSegm[0].Points[1].Y = 0.;
   }

   //QSegment::calcOpen_and_Shadow_Segments(BladeSgm,  SgmBarr0,  arrWindSKBrlsV0
                 //      , arrOpenSegm, &iquantOpenSgm, &shadowSegm, &iquantShadowSgm);
   ///
   double valResistMomRez = 0.;
   for (int i =0; i < iquantOpenSgm; ++i)
   {
     valResistMomRez += calcAirResistFunc_With_Limits(VAlOm, valW
          , arrOpenSegm[i].Points[0].X,arrOpenSegm[i].Points[1].X,pf);
   }

   for (int i =0; i < iquantShadowSgm; ++i)
   {
     valResistMomRez += calcAirResistFunc_With_Limits(VAlOm, 0.
          , shadowSegm.Points[0].X,shadowSegm.Points[1].X,pf);
   }
   return valResistMomRez;
}

// вычисление момента сопротивления отрезка лопасти или его производной
// при условии постоянной скорости ветра на этом отрезке
// INPUT:
//VAlOm - угловая скорость вращзения БРЛС
//valW - нормальная составляющая скорости ветра к к лопасти в плоскости вращения
// ValLowLim - нижний предел интегрирования
// ValUpperLim - венрхний предел интегрирования
// pf - прототип функции
// ПОЯСНЕНИЯ:
// на заданном отрезке лопасти направления воздушного потока,
//  который обдувает лопасть может менять знак.
// Например, на дальнем краю лопасти, где линейная скорость элементарного участка лопасти больше
// направление воздушного потока может иметь знак "+"
// а на ближнем, знак "-"
// функция сопротивления воздуха на элементарном участке лопасти
// зависит от квадрата воздушной скорости
// поэтому, отрезок лопасти возможно надо разбить на 2 отрезка
// на каждом из которых направления воздушного потока постоянно
double  QBrls::calcAirResistFunc_With_Limits( const double VAlOm,  const double valW
                     ,  const double ValLowLim,  const double ValUpperLim
   ,double (*pf)(const double ,const double ,const double ,const double  ,const double,const double ))
{
    // 1. ннулвая угловая скорость
    if (fabs(VAlOm < 0.000001))
    {
        return sign__(valW ) * mCx*ATM_RoN0 *valW * valW/ 2.;
    }

    ///

    // 2. нахождение точки, в которой скоростной напор равен нулю
    // если эта точка лежит в пределах интегрирования, то отрезок интегрирования разбивается на 2.
    double valx0 = valW/ VAlOm;
    ///

    // 3. заготовка 2- отрезков интегрирования
    double arrLowLim[2] = {0.}, arrUpperLim[2] = {0.};
    arrLowLim[0] = ValLowLim;
    arrUpperLim[0] = ValUpperLim;
    int quantSegmItegr = 1;
    ///

    // 4. проверка, что точка valx0 лежит внутри и в этом случае переопределение
    if ((valx0 -ValLowLim )* (valx0 -ValUpperLim) < 0.)
    {
        quantSegmItegr = 2;
        arrUpperLim[0] = valx0;
        arrLowLim[1] = valx0;
        arrUpperLim[1] = ValUpperLim;
    }
    ///

    // 5. вычисление интегралов
    double valMomRez = 0.;
   // double valSignum =  sign__(valW - VAlOm *(valx0 - 1.));
    for (int i = 0; i < quantSegmItegr; ++i)
    {       
        double valMiddle = (arrLowLim[i]+ arrUpperLim[i])/ 2.;
        valMomRez += pf(mCx, VAlOm, valW, arrLowLim[i], arrUpperLim[i], mHeight) * sign__(valW - VAlOm *valMiddle);
    }

    return valMomRez;
}

//--------------------------------------
double  calcIntegralAirMom(const double Cx, const double VAlOm,const double  valW,const double  VAla
                                  ,const double  VAlb,const double  VAlh)
{
    double temp0= (VAlb * VAlb - VAla * VAla) * valW* valW /2.;
    double temp1= -(VAlb * VAlb* VAlb - VAla * VAla* VAla) * 2. * valW * VAlOm /3.;
    double temp2= (VAlb * VAlb* VAlb* VAlb - VAla * VAla*VAla* VAla) * VAlOm *  VAlOm / 4.;
    return  Cx * ATM_RoN0 / 2. *(temp0 - temp1 + temp2) *VAlh;
}
//-------------------------------------------------------
//--------------------------------------
double  calcIntegral_dAirMom_po_dOmega(const double Cx, const double VAlOm,const double  valW,const double  VAla
                                  ,const double  VAlb,const double  VAlh)
{
    double temp1= -(VAlb * VAlb* VAlb - VAla * VAla* VAla) * 2. * valW  /3.;
    double temp2= (VAlb * VAlb* VAlb* VAlb - VAla * VAla*VAla* VAla) * VAlOm  / 2.;
    return  Cx * ATM_RoN0 / 2. *(-temp1 + temp2)*VAlh;
}
//------------------------------------------------------------
double QBrls::calc_dMa_po_Cx(const double VAlOm)
{
    return  ATM_RoN0 * mHeight * mLength * mLength * mLength * mLength * VAlOm * VAlOm/ 128.;
}

//------------------------------------------------------------
double  QBrls::calc_d2Ma_po_dOm_po_Cx(const double VAlOm)
{

    const double VAlb = mLength/2.;
    double temp2= (VAlb * VAlb* VAlb* VAlb ) * VAlOm ;
    return   ATM_RoN0 / 2. * temp2 *mHeight;
}
//---------------------------------
int QBrls::createInputDataReport(wchar_t*FileName, const bool bHeader)

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
fprintf(fw,"      тип - БРЛС\n");
fprintf(fw,"      Параметры   нагрузки\n");

fprintf(fw,"  момент инерции нагрузки(кг*м*м) = %5.2f\n",mJPayLoad);
fprintf(fw,"  коэффиц сопротивления воздуха = %5.2f\n",mCx);
fprintf(fw,"  коэффиц трения воздуха = %5.2f\n",mCv);
fprintf(fw,"  полная длина штанги/диаметр вращения (м) = %4.2f\n",mLength);
fprintf(fw,"  высота штанги (м) = %4.2f\n",mCv);
fprintf(fw,"  коэффиц трения воздуха = %5.2f\n",mCv);
fprintf(fw,"  коэффиц трения воздуха = %5.2f\n",mCv);

 fclose(fw);

}



