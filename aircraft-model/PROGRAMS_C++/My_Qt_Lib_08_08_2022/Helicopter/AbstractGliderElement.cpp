#include "AbstractGliderElement.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include "MatrixProccess.h"


TAbstractGliderElement::TAbstractGliderElement()
{
    // площадь
    mS = 0.;

    // связанная система координат элемента планера в осях связанной системы координат вертолета

     mPlaneSvSK = TLongPlane();

    // координаты точки приложенияя аэродинамической силы в связанной системе координат элемента планера
     memset(marrAirP, 0, 3 * sizeof(long double));


    // Коэффициент лобового сопротивленияf
     mCx0 = 0.;

    // Коэффициент подъемной силы
     mCy0 = 0.;

    // Критичекий угол атаки, при котором коэффициент подъемной силы максималный
    mAlfCrit = 0.;

}

// Конструктор копирования
  TAbstractGliderElement::TAbstractGliderElement (const TAbstractGliderElement &R)
 {
      mS = R.mS;
      mPlaneSvSK   = R.mPlaneSvSK ;
      mCx0  = R.mCx0 ;
      mCy0 = R.mCy0;
      mAlfCrit = R.mAlfCrit;
      memcpy(marrAirP, R.marrAirP,  3 * sizeof(long double));
 }
 // оператор присваивания
  TAbstractGliderElement TAbstractGliderElement::operator=(TAbstractGliderElement  R)
 {
      mS = R.mS;
      mPlaneSvSK   = R.mPlaneSvSK ;
      mCx0  = R.mCx0 ;
      mCy0 = R.mCy0;
      mAlfCrit = R.mAlfCrit;
      memcpy(marrAirP, R.marrAirP,  3 * sizeof(long double));
      return *this ;
 }

TAbstractGliderElement:: TAbstractGliderElement(const long  double S
         , const TLongPlane PlaneSvSK, const long  double Cx0, const long  double Cy0, const long  double AlfCrit
                                                  ,long  double *arrAirP)
 {
     mS = S;
     mPlaneSvSK = PlaneSvSK;
     mCx0 = Cx0;
     mCy0 = Cy0;
     mAlfCrit = AlfCrit;
     memcpy(marrAirP, arrAirP,  3 * sizeof(long double));

}
 //-------------------------------------------------------
     // парам констр 2
     // arrInpDataGliders [13] - массив параметров
     //arrInpDataGliders  [0] - mS
     // arrInpDataGliders [1],[2],[3] - радиус вектор начала координат собственной СвСК элемента планера
     // arrInpDataGliders [4],[5],[6] - координаты точки приложенгия аэродинамической силы в собственной СвСК элем планера
     // arrInpDataGliders [7],[8],[9] - mCx0; mCy0, mAlfCrit
     // arrInpDataGliders [10],[11],[12] - эйлеровы углы задающие взаимное расположенгия осей СвСК и собственной СвСК элемента планнера
     //
TAbstractGliderElement::TAbstractGliderElement(long double *arrInpDataGliders)
{
    // площадь
     mS = arrInpDataGliders[0];

    // связанная система координат элемента планера в осях связанной системы координат вертолета

    TLongPlane PlaneSvSK;

    // координаты точки приложенияя аэродинамической силы в связанной системе координат элемента планера
    marrAirP[0] = arrInpDataGliders[4] ;
    marrAirP[1] = arrInpDataGliders[5] ;
    marrAirP[2] = arrInpDataGliders[6] ;

    // Коэффициент лобового сопротивления
   mCx0 = marrAirP[0] = arrInpDataGliders[7] ;


    // Коэффициент подъемной силы
   mCy0 = arrInpDataGliders[8];

    // Критичекий угол атаки, при котором коэффициент подъемной силы максималный
   mAlfCrit = arrInpDataGliders[9] * M_PI/ 180.;

   //
   PlaneSvSK.marrS0[0] =  arrInpDataGliders [1];
   PlaneSvSK.marrS0[1] =  arrInpDataGliders [2];
   PlaneSvSK.marrS0[2] =  arrInpDataGliders [3];
   ///

   //
   long  double matrPX[9] = {0.}, matrPY[9] = {0.}, matrPZ[9] = {0.};
   calcRotMatrX(arrInpDataGliders[10]/180. * M_PI, matrPX);
   calcRotMatrY(arrInpDataGliders[11]/180. * M_PI, matrPY);
   calcRotMatrZ(arrInpDataGliders[12]/180. * M_PI, matrPZ);
   MtrxMultMatrx_MultMatrx(matrPX,matrPY,matrPZ, 3, PlaneSvSK.marrF);

   mPlaneSvSK = PlaneSvSK;



}


//----------------------------------------------------------------------------------

 // вычисление матрицы дополнтельного разворота СвСК элемента планера
 // относительно СвСК вертолета
 // для элемента типа руля эта матрица вычисляется,
 // для простого элемента приравнивается единичной
 // INPUT:
 // VAlAlfa - угол поворота руля
 //OUTPUT:
 // arrMwave[9] - матрица 3х3
void TAbstractGliderElement::calMtrxMwave(const long double VAlAlfa,long  double *arrMwave)
{

}

//-----------------------------------------------------------------------------------
// функция вычисления суммарной аэродинамической силы
// и точки ее приложения в СвСК вертолета
// INPUT:
// arrUa[3] - вектор воздушной скорости вертолета в СвСК
// VAlFi_rv - угол поворота элемента планера
// VAlRo - плотность воздуха
// OUTPUT:
// arrAeroF[3] - вектор суммарной аэродинамической силы в СвСК
// arrAeroMom[3] - момент силы в СвСК
void TAbstractGliderElement::calAeroF_and_Mom(long double *arrUa, const long  double VAlRo
                                                   ,const long double VAlFi_rv,long  double *arrAeroF, long double *arrAeroMom)
{
    // модуль воздушной скорости
      long  double valModUa = Norm3(arrUa);
    // 0. проверка того, что воздушная скорость не равна нулю
       // если равна нулю, то аэродинамическая сила равна нулю
       if(valModUa < DBL_MIN )
       {
           memset(arrAeroF, 0, 3 * sizeof(long double));
           memset( arrAeroMom, 0, 3 * sizeof(long double));
           return;
       }
    // 1. Вычисление матрицы перехода из связанной сиситемы координат
    // элемента планера в СвСК вертолета
    long double arrMtrxTransf[9] ={0.}, arrMwave[9] ={0.};
    calMtrxMwave(VAlFi_rv, arrMwave);
    MtrxMultMatrx(mPlaneSvSK.marrF,3, 3, arrMwave,3, arrMtrxTransf);
    ///

    // 2. Вычисление вектора воздушной скорости элемента планера Ft*Ua
    long  double arrUaEl[3] = {0.};
    MtrxMultMatrx(arrUa,1, 3, arrMtrxTransf,3, arrUaEl);
    ///

    // 3.Вычисление угла атаки ЭП
     long  double temp = fabsl(arrUaEl[1]/ valModUa);

    long double valAlfEl = (temp > 0.001)? asinl(arrUaEl[1]/ sqrtl( arrUaEl[0] * arrUaEl[0] + arrUaEl[1] * arrUaEl[1])):0.;
    ///

    // 4. вычисление угла скольжения
   // double valBetEl = asin(arrUaEl[2]/ valModUa);
    ///

    // 5. вычисление площади миделя
   long  double valSMid = mS * sinl(fabsl(valAlfEl));// площадь проекции миделя на плоскость перпендикулярную вектору воздушной скорости
    ///

    // 6. вычисление модуля силы лобового сопротивления
    long double valFLob = mCx0 * VAlRo * valModUa  * valModUa/2. * valSMid;
    ///

    // 7. вычисление модуля подъемной силы
    long double valFPod = 0;
    valFPod  = fncPodFCoeff( valAlfEl) * VAlRo * valModUa  * valModUa/2. * valSMid;

    ///

    // 8.Вычисление вектора силы лобового сопротивления
    long double arrFLob[3] = {0.};
    MatrxMultScalar(arrUa, 1, 3, valFLob/(valModUa ),arrFLob);
    ///

    // 9. Вычисление вектора подъемной силы

    long  double arrFPod[3] = {0.};
     //double arr_ny[3] = {0.}, arrt[3] = {0.};
    // arr_ny[0] = arrMtrxTransf[1];
    // arr_ny[1] = arrMtrxTransf[4];
    // arr_ny[2] = arrMtrxTransf[7];
//
     //OuterProduct(arrUa , arr_ny, arrt) ;
     //MatrxMultScalar(arrt, 1, 3, valFPod/valModUa,arrFPod);
     arrFPod[0] = sinl(valAlfEl) * valFPod;
     arrFPod[1] = cosl(valAlfEl) * valFPod;
     ///

     // 10. Вычисление результирующего вектора аэродинамической силы
     MtrxSumMatrx(arrFPod, arrFLob,1, 3, arrAeroF) ;
     ///

     // 11. Функция вычисления радиус-вектора точки приложения
     // суммарной аэродинамической силы в СвСК вертолета
     long double arrPointAeroF[3] = {0.};
      MtrxMultMatrx(arrMtrxTransf,3, 3, mPlaneSvSK.marrS0, 1, arrPointAeroF) ;
      ///

      // 12. вычисление момента аэродинамической силы
      OuterProduct(arrPointAeroF , arrAeroF, arrAeroMom) ;

}

// функция коэффициента подъемной силы
long double TAbstractGliderElement::fncPodFCoeff(const long double valAlfEl)
{
    long double x =  valAlfEl/ mAlfCrit;
    if (x >=0.)
    {
        return mAlfCrit * fncTay(x);
    }
    else
    {
        return -mAlfCrit * fncTay(-x);
    }
}

long double TAbstractGliderElement::fncTay(const long  double x)
{
    if (x > 2.75)
    {
        return 0.;
    }
    if (x <= 0.75)
    {
        return x;
    }
    else
    {
        if ( x <= 1.25)
        {
            return -2. * x * x + 4. * x - 9./8.;
        }
        else
        {
            return 1. / 3. * (x - 2.75)* (x - 2.75);
        }
    }
}

// функция производной коэффициента подъемной силы по углу атаки
long double TAbstractGliderElement::fncDerivPodFCoeff(const long  double valAlfEl)
{
    double x =  valAlfEl/ mAlfCrit;
    if (x >=0.)
    {
        return mAlfCrit * fncDerivTay(x);
    }
    else
    {
        return mAlfCrit * fncDerivTay(-x);
    }
}

long double TAbstractGliderElement::fncDerivTay(const long double x)
{
    if (x > 2.75)
    {
        return 0.;
    }
    if (x <= 0.75)
    {
        return 1.;
    }
    else
    {
        if ( x <= 1.25)
        {
            return 4. * (1. -x);
        }
        else
        {
            return 2. / 3. * (x - 2.75);
        }
    }
}
//---------------------------------------------------------------------------
// вычисление матрицы ортов повернутой относительно оси X сиситемы координат
// на угол a против часовой стрелки
// Иными словами, если исходную правую сиситему координат
// повернуть относительно оси X на угол a, то это матрица столбцов координат ортов новой (развернутой)
// сиситемы координат в осях старой
void calcRotMatrX(const long  double a, long double*matrP)
{
  memset(matrP,0,9*sizeof(long  double));
  matrP[0] = 1;
  matrP[4] = cosl(a);
  matrP[5] = -sinl(a);
  matrP[7] = sinl(a);
  matrP[8] = cosl(a);

}

//---------------------------------------------------------------------------
// вычисление матрицы ортов повернутой относительно оси Y сиситемы координат
// на угол a против часовой стрелки
// Иными словами, если исходную правую сиситему координат
// повернуть относительно оси Y на угол a, то это матрица столбцов координат ортов новой (развернутой)
// сиситемы координат в осях старой
void  calcRotMatrY(const long  double a,long  double*matrP)
{
  memset(matrP,0,9*sizeof(long  double));
  matrP[4] = 1;
  matrP[0] = cosl(a);
  matrP[2] = sinl(a);
  matrP[6] = -sinl(a);
  matrP[8] = cosl(a);

}

//---------------------------------------------------------------------------
// вычисление матрицы ортов повернутой относительно оси Z сиситемы координат
// на угол a против часовой стрелки
// Иными словами, если исходную правую сиситему координат
// повернуть относительно оси Z на угол a, то это матрица столбцов координат ортов новой (развернутой)
// сиситемы координат в осях старой
void  calcRotMatrZ(const long  double a,long  double*matrP)
{
  memset(matrP,0,9*sizeof(long  double));
  matrP[8] = 1;
  matrP[0] = cosl(a);
  matrP[1] = -sinl(a);
  matrP[3] = sinl(a);
  matrP[4] = cosl(a);
}









