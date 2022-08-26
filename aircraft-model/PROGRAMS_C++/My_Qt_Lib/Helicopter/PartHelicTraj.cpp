#include "PartHelicTraj.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include "MatrixProccess.h"
#include "Comp.h"
#include "BallanceCalc.h"
#include <QDebug>







TPartHelicTraj::TPartHelicTraj()
{
    //  вертолет
    mHelic =THelic ();
    // атмосфера
    mEnvironment = TEnvironment ();

    //

    memset(marrPhaseVect, 0, QUantCurNZSKVarsVS * sizeof(long double));
    mTimeCur = 0.;

    mTBegin = 0.;

    memset(marrSvSK_Force, 0, 3  * sizeof(long double));


}

//---------------------------------------------------------------------------

 // парам констр
TPartHelicTraj::TPartHelicTraj( const THelic Helic,const TEnvironment  Environment
                                , const long double VAlTimeCur, const long double VAlTBegin)
{
   mHelic = Helic;
   mEnvironment =Environment;
   mTimeCur = VAlTimeCur;
   mTBegin = VAlTBegin;

   memset(marrPhaseVect, 0, QUantCurNZSKVarsVS * sizeof(long double));
   memset(marrSvSK_Force, 0, 3  * sizeof(long double));

}


//---------------------------------------------------------------------------

 // парам констр равном прямолин движения
TPartHelicTraj::TPartHelicTraj( const THelic Helic,const TEnvironment  Environment , const long double VAlY
                                , const long  double VAlVx, const long double VAlTimeCur, const long double VAlTBegin)
{
   mHelic = Helic;
   mEnvironment =Environment;
   mTimeCur = VAlTimeCur;
   mTBegin = VAlTBegin;

   memset(marrPhaseVect, 0, QUantCurNZSKVarsVS * sizeof(long double));

   memset(marrSvSK_Force, 0, 3  * sizeof(long double));
}



 //---------------------------------------------------------------------------

 // оператор присваивания
 TPartHelicTraj TPartHelicTraj::operator=(TPartHelicTraj  R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
     mTBegin = R.mTBegin;    // memcpy(marrSteadySolution, R.marrSteadySolution, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));
     mTimeCur = R.mTimeCur;

     return *this;
 }

 // конструктор копирования
 TPartHelicTraj::TPartHelicTraj (const TPartHelicTraj &R)
 {
     mHelic = R.mHelic;
     mEnvironment = R.mEnvironment;
     memcpy(marrPhaseVect, R.marrPhaseVect, QUantCurNZSKVarsVS * sizeof(long double));
     mTBegin = R.mTBegin;    // memcpy(marrSteadySolution, R.marrSteadySolution, (QUantCurNZSKVarsVS -1)  * sizeof(long double));
     memcpy(marrSvSK_Force, R.marrSvSK_Force, 3  * sizeof(long double));

     mTimeCur = R.mTimeCur;
}




//-------------------------------------------------------------------

// вычисление правой части системы дифференциальных уравнений твердого тела

void TPartHelicTraj::calcRightMemberOfDifEqSystem( const long double   valFi, const  long double   valKappa
                                                , const  long double   valEtta, const  long double   valDelFi
                                                , long double   * arrf, long double   *arrSvSK_Force, long double   *arrSvSK_Moment
                                               , long double   * arrSvSK_ShaftForce0, long double   *arrSvSK_ShaftMoment0
                                               , long double   *arrSvSK_AirForce0, long double   *arrSvSK_AirMoment0)
 {
    arrf[0] = marrPhaseVect[3];
    arrf[1] = marrPhaseVect[4];
    arrf[2] = marrPhaseVect[5];
    // вычисление плотности атмосферы
    marrPhaseVect[12] = mEnvironment.calcAirDensity(marrPhaseVect[1]);
   // long double   valAtmRo = mEnvironment.calcAirDensity(marrPhaseVect[1]);
    ///

    // вычисление вектора воздушной скорости в СвСК
       // вычисление вектора набегающего воздушного потока в  НЗСК (ось X на Север, Z на Восток, Y в Зенит)

    long double   arrNZSK_Va [3] ={0.}, arrWindV[3] = {0.};
        mEnvironment.createVectWindV(arrWindV);
        MtrxMinusMatrx(arrWindV, &(marrPhaseVect[3]),1, 3, arrNZSK_Va);
       // arrNZSK_Va [0] = mEnvironment.mWind_V * cosl(mEnvironment.mWind_Alf) - marrPhaseVect[3];
       // arrNZSK_Va [1] = mEnvironment.mWind_VertV - marrPhaseVect[4];
        //arrNZSK_Va [2] = mEnvironment.mWind_V * sinl(mEnvironment.mWind_Alf) - marrPhaseVect[5];
        ///

        // перевод вектора набегающего воздушного потока из НЗСК в СвСК
        long double   arrSvSK_Va [3] ={0.};
        long double   arrMtrxTransf_NZSK_SvSK[9] = {0.};
        calcMtrxTransf_CurNZSK_To_SvSK (arrMtrxTransf_NZSK_SvSK);
        MtrxMultMatrx(arrMtrxTransf_NZSK_SvSK,3, 3, arrNZSK_Va, 1, arrSvSK_Va);

    ///


  // вычисление   векторов равнодействующей аэродинамической силы и момента сил в СвСК
        // без учета силы тяжести (аэродинамика и винт учитываются только)

    mHelic.calcRezF_and_Moment_SvSK(&marrPhaseVect[6], valFi, valKappa, valEtta, valDelFi
            ,marrPhaseVect[12], arrSvSK_Va, arrSvSK_Force, arrSvSK_Moment, arrSvSK_ShaftForce0,arrSvSK_ShaftMoment0
            ,arrSvSK_AirForce0,arrSvSK_AirMoment0);
    ///

    // перевод равнодействующей силы в НЗСК
    long double   arrMtrxTransf[9] = {0.}, arrNZSK_Force[3] = {0.};
    calcMtrxTransf_SvSK_To_CurNZSK (arrMtrxTransf);
    MtrxMultMatrx(arrMtrxTransf,3, 3, arrSvSK_Force, 1, arrNZSK_Force);
    ///
    // учет силы тяжести
    arrNZSK_Force[1] -= G_ZEMLI * mHelic.mHelicMass; // 15.10.18 !!!!!!
    ///
    MatrxDivideScalar(arrNZSK_Force, 1, 3, mHelic.mHelicMass, &(arrf[3]));
   // arrBuffer[lenBuffer * QuantDebug_Buffer ] = ((double)lenBuffer)* 0.1;
   // arrBuffer[lenBuffer *  QuantDebug_Buffer  + 1] = arrf[4];



    ///

    // вычисление правой части 3-х уравнений угловых скоростей
    long double   arrT0[3] = {0.},  arrT1[3] ={0.},  arrT2[3] ={0.}, arrJ[9] = {0.}, arrJInv[9] = {0.};
    for (int j = 0; j < 9; ++j)
    {
      arrJ[j] = mHelic.marrJ[j];
    }
    MtrxMultMatrx(arrJ   ,3, 3, &(marrPhaseVect[6]), 1, arrT0) ; // J*omega
    OuterProduct(&(marrPhaseVect[6]) , arrT0, arrT1);
    MtrxMinusMatrx(arrSvSK_Moment, arrT1,1, 3,  arrT2);
    InverseMtrx3(arrJ,  arrJInv);
    MtrxMultMatrx(arrJInv   ,3, 3, arrT2, 1,&(arrf[6])) ;

    ///

    // вычисление правой части 3-х уравнений углов Эйлера
    long double   psi =   marrPhaseVect[11];
    long double   nu =    marrPhaseVect[10];
    long double   gamma = marrPhaseVect[9];

    calcRightPartForEilers( gamma,  nu, psi, &(marrPhaseVect[6]), &(arrf[9])); // ПРОВЕРИТЬ!!!!
    ///

    // вычисление раванодействующей силы в СвСК
    MtrxMultMatrx(arrMtrxTransf_NZSK_SvSK,3, 3, arrNZSK_Force, 1, arrSvSK_Force);
 }

//---------------------------------------------------

//---------------------------------------------------
void TPartHelicTraj::calcRightPartForEilers(const long double gamma, const  long double nu
              , const  long  double psi, long  double *arrOmega, long  double *arrRightPart)
{
     long double arrT3[9] = {0.};
    arrT3[0] = 1.;
    arrT3[1] = -tanl(nu) * cosl (gamma);
    arrT3[2] =  tanl(nu) * sinl (gamma);

    arrT3[3] = 0.;
    arrT3[4] = sinl (gamma);
    arrT3[5] = cosl (gamma);

    arrT3[6] =  0.;
    arrT3[7] =  cosl(gamma)/ cosl(nu);
    arrT3[8] = -sinl(gamma)/ cosl(nu);

    MtrxMultMatrx(arrT3 ,3, 3, arrOmega, 1, arrRightPart) ;
}
//---------------------------------------------------



//----------------------------------------------------------------------------------------
// вычисление вектора скорости набегающего воздушного потока в СвСК для сиситемы диф уравнений в НЗСК
// OUTPUT:
//arrUaSvSK[3] - вектор скорости набегающего воздушного потока в СвСК
void TPartHelicTraj::calcUaSvSK_Case_NZSK(long double *arrUaSvSK)
{
    // вычисление вектора набегающего воздушного потока  в  НЗСК (ось X на Север, Z на Восток, Y в Зенит)
    long double   arrNZSK_Va [3] ={0.}, arrWindV[3] = {0.};
        mEnvironment.createVectWindV(arrWindV);
        MtrxMinusMatrx(arrWindV, &(marrPhaseVect[3]),1, 3, arrNZSK_Va);
  //  long double arrNZSK_Va [3] ={0.};
   // arrNZSK_Va [0] = mEnvironment.mWind_V * cosl(mEnvironment.mWind_Alf) - marrPhaseVect[3];
   // arrNZSK_Va [1] = mEnvironment.mWind_VertV - marrPhaseVect[4];
  //  arrNZSK_Va [2] = mEnvironment.mWind_V * sinl(mEnvironment.mWind_Alf) - marrPhaseVect[5];
    ///
    // перевод вектора воздушной скорости из НЗСК в СвСК
    long double arrMtrxTransf_NZSK_SvSK[9] = {0.};
    calcMtrxTransf_CurNZSK_To_SvSK (arrMtrxTransf_NZSK_SvSK);
    MtrxMultMatrx(arrMtrxTransf_NZSK_SvSK,3, 3, arrNZSK_Va, 1,arrUaSvSK);
    ///
}





//--------------------------------------------------------------------------------


//-------------------------------------------------------------------------
// вычисление матрицы направляющих косинусов СвСК в осяз НЗСК
// или, что тоже самое, матрица перехода из СвСК в НЗСК
//OUTPUT:
//arrMtrxTransf[9]
void TPartHelicTraj::calcMtrxTransf_SvSK_To_CurNZSK (long double* arrMtrxTransf)
{
    long double psi =   marrPhaseVect[11];
    long double nu =    marrPhaseVect[10];
    long double gamma = marrPhaseVect[9];

    arrMtrxTransf[0] = cosl(psi) * cosl(nu);
    arrMtrxTransf[1] = sinl(psi) * sinl (gamma) - cosl(psi)* sinl(nu) * cosl(gamma);
    arrMtrxTransf[2] = sinl(psi) * cosl (gamma) + cosl(psi)* sinl(nu) * sinl(gamma);

    arrMtrxTransf[3] = sinl(nu);
    arrMtrxTransf[4] = cosl(nu) * cosl(gamma);
    arrMtrxTransf[5] = -cosl(nu) * sinl(gamma);

    arrMtrxTransf[6] = -sinl(psi)* cosl(nu);
    arrMtrxTransf[7] = cosl(psi) * sinl (gamma) + sinl (psi) * sinl (nu)* cosl(gamma);
    arrMtrxTransf[8] = cosl(psi) * cosl (gamma) - sinl(psi) * sinl (nu) * sinl(gamma);
}

//-------------------------------------------------------------------------

// вычисление матрицы направляющих косинусов СвСК в осяз НЗСК
// или, что тоже самое, матрица перехода из СвСК в НЗСК
//OUTPUT:
//arrMtrxTransf[9]
void TPartHelicTraj::calcMtrxTransf_SvSK_To_CurNZSK_ (const long double gamma
                     ,const long double nu,const long double psi, long double* arrMtrxTransf)
{
    arrMtrxTransf[0] = cosl(psi) * cosl(nu);
    arrMtrxTransf[1] = sinl(psi) * sinl (gamma) - cosl(psi)* sinl(nu) * cosl(gamma);
    arrMtrxTransf[2] = sinl(psi) * cosl (gamma) + cosl(psi)* sinl(nu) * sinl(gamma);

    arrMtrxTransf[3] = sinl(nu);
    arrMtrxTransf[4] = cosl(nu) * cosl(gamma);
    arrMtrxTransf[5] = -cosl(nu) * sinl(gamma);

    arrMtrxTransf[6] = -sinl(psi)* cosl(nu);
    arrMtrxTransf[7] = cosl(psi) * sinl (gamma) + sinl (psi) * sinl (nu)* cosl(gamma);
    arrMtrxTransf[8] = cosl(psi) * cosl (gamma) - sinl(psi) * sinl (nu) * sinl(gamma);
}

// вычисление матрицы направляющих косинусов  НЗСК в осях СвСК
// или, что тоже самое, матрица перехода из НЗСК  в СвСК
//OUTPUT:
//arrMtrxTransf[9]
void TPartHelicTraj::calcMtrxTransf_CurNZSK_To_SvSK (long double* arrMtrxTransf)
{
    long double arrT[9] = {0.};
    calcMtrxTransf_SvSK_To_CurNZSK (arrT) ;
    MatrTransp(arrT, 3, 3, arrMtrxTransf);
}

//-------------------------------------------------------------------------
// вычисление матрицы направляющих косинусов СвСК в осях Скоростной системы координат(СкорСК)
// или, что тоже самое, матрица перехода из СвСК в СкорСК в
//OUTPUT:
//arrMtrxTransf[9]
//Ось OY СкорСК лежит в плоскости симметрии ЛА и нерпендикулярна вектору скорости ЛА
// ось OX СкорСК направлена по вектору скорости
// Ось OZ СкорСК дополняет тройку векторов до правой
// ищемединичные  векторы скоростной сиситемы координат в осях СвСК
// и распрологаем их по строкам
bool TPartHelicTraj::calcMtrxTransf_SvSK_To_SkorSK (long double* arrMtrxTr)
{
    // 1.отсекаем случай когда скорость равна 0.
    // в этом случае скоростная система координат не определена
    if (Norm3( &(marrPhaseVect[3])) < 0.0001)
    {
        return false;
    }
    ///

    // 2. наъходим направляющие косинусы вектора скорости в осях СвСК
    // это будет вектор e1 СкорСК в осях СвСК
    // 2.1 для этого переводим вектор скорости из ТНЗСК в СвСК
    long double  arrMtrxPer_CurNZSK_to_SvSK[9] ={0.};
    calcMtrxTransf_CurNZSK_To_SvSK (arrMtrxPer_CurNZSK_to_SvSK);
    MtrxMultMatrx(arrMtrxPer_CurNZSK_to_SvSK, 3, 3, &(marrPhaseVect[3]),1, arrMtrxTr) ;
    //2.2 нормализуем вектор скорости в СвСК
    NormalizeVect3(arrMtrxTr);

    // 2.3 отсекаем случай когда вектор скорости параллелен вектору OY СвСК
    // в этом случае скоростная система координат не определена
    // это реализуется тогда, когда вторая координата вектора arrMtrxTr по модулю близка к 1
   // if (fabsl(arrMtrxTr[1]) > 0.9999)
   // {
    //    return false;
    //}



    // 3. формируем второй вектор OY
    // ищем вектор OY в осях СвСК

    // 3.1 отсекаем случай движения вверх (по оси OY СвСК)
    if (fabsl(arrMtrxTr[1]) > 0.999999)
    {
        arrMtrxTr[0] = 0.;
        arrMtrxTr[2] = 0.;
        arrMtrxTr[1] = SIGNUM(arrMtrxTr[1]);
      arrMtrxTr[3] = - SIGNUM(arrMtrxTr[1]) ;
      arrMtrxTr[8] = 1.;
      return true;
    }

    // 3.2 отсекаем случай движения в бок (по оси OZ СвСК)
    if (fabsl(arrMtrxTr[2]) > 0.999999)
    {
        arrMtrxTr[0] = 0.;
        arrMtrxTr[1] = 0.;
        arrMtrxTr[2] = SIGNUM(arrMtrxTr[2]);
        arrMtrxTr[4] = 1.;
        arrMtrxTr[6] = -arrMtrxTr[2];
        return true;

    }
    arrMtrxTr[3] = -arrMtrxTr[1];
    arrMtrxTr[4] = arrMtrxTr[0];
    NormalizeVect3(&(arrMtrxTr[3]));



    // 4. формируем третий вектор путем векторного произведения первого на второй
    OuterProduct(arrMtrxTr , &(arrMtrxTr[3]), &(arrMtrxTr[6])) ;
    return true;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// вычисление матрицы направляющих косинусов СвСК в осях Скоростной системы координат(СкорСК)
// или, что тоже самое, матрица перехода из СвСК в СкорСК в
//OUTPUT:
//arrMtrxTransf[9]
//Ось OY у СвСК и СкорСК совпадают
// ось OX СкорСК направлена по вектору скорости
// Ось OZ СкорСК дополняет тройку векторов до правой
bool TPartHelicTraj::calcMtrxTransf_SkorSK_To_SvSK (long double* arrMtrxTr)
{
    long double arrTemp[9] = {0.};
    if (!calcMtrxTransf_SvSK_To_SkorSK (arrTemp))
    {
        return false;
    }
    MatrTransp(arrTemp, 3, 3, arrMtrxTr);
    return true;

}

//-------------------------------------------------------------------------
// вычисление нормальной перегрузки
// В ГОСТе на динамику ЛА приведено неправильное определение нормальной перегрузки
// Здесь расчитывается правильно
// В теор механике вектор мгновенного ускорения тела раскладывается на составляющие по двум осям - оси мгновенной скорости
// и оси ценртростремительного (или центробежного) ускорения
long double TPartHelicTraj::calcPeregrNorm()
{
    // отсекаем случай когда СкорСК не определена и нормальная перегрузка
    // равна 0
    long double arrMtrxTr[9] = {0.};
    if (!calcMtrxTransf_SvSK_To_SkorSK ( arrMtrxTr))
    {
        return 0.;
    }
    ///

    // вычисление вектора равнодейсвующей силы в СкорСК
    long double arrSkorSK_Force[3] = {0.};
    MtrxMultMatrx(arrMtrxTr,3, 3, marrSvSK_Force,1, arrSkorSK_Force) ;
    ///

    return fabsl(arrSkorSK_Force[2]/ mHelic.mHelicMass/ G_ZEMLI);
}


//---------------------------------------------------------------
// вычисление числа Маха и вектора часных производных числа Маха по x
// INPUT:
//valTay - нормальная виртуальная  температура
// valVVozd - воздушная скорость

long double TPartHelicTraj::calcMach( long double valTay,long   double valVVozd)
{
 return valVVozd *sqrt(ATM_TAYN0/ valTay)/ ATM_AN0 ;

}
//-----------------------------------------------------------------------------------




   //-----------------------------------------------------------------------------------
    // перемещение вертолета в ссответсвии с полетным заданием
    // на время полета VAlFlyTime
    // INPUT:
    // VAlFlyTime - время полета
    // VAlStIntegr - шаг интегрирования
   // LEnarrC - длина массива передаточных чисел = 44
    //OUTPUT:
    // arrBuff - буфер (массив) для накопления информации
    // Если arrBuff==NULL то информация не накапливается
    // *pquantSteps - к-во шагов интегрирования

    // вектор выходной информации arrOutputVect[25]
   // составляющие воздушной скорости в связанной системе координат; 3
   // составляющие путевой скорости в связанной системе координат; 3
   // приборная скорость;                                          1
   // бароинерциальная вертикальная скорость;                      1
   // геометрическая высота;                                       1
   // абсолютная барометрическая высота;                           1
   // относительная барометрическая высота;                         1
   // геодезическая высота;                                        1
   // бароинерциальная высота;                                     1
   // угол сноса;                                                  1
   // угол атаки;                                                  1
   // угол скольжения;                                              1
   // скорость и направление ветра;                                2
   // крен;                                                        1
   // тангаж;                                                      1
   // курс истинный;                                                1
    //магнитное склонение;                                          1
   // температура наружного воздуха;                               1
   // текущие координаты местоположения объекта в географической системе координат.

   // задан текущзий фазовый вектор вертолета в *this
   // ьребуется обеспечить траекторию на постоянной высоте VAlY
   // со скоростью VAlVx и курсовым углом =0

    void TPartHelicTraj::Moving( const long double VAlFlyTime, const long double VAlStIntegr
                                 ,const long double VAlModelStep, double *arrBuff, int *pquantSteps)
    {
        // признак полета
        bool bFly = false;
        ///

        *pquantSteps = 0;




     int quantSt = (int)(VAlFlyTime/ VAlModelStep);
     long double arrSvSK_Force[3] = {0.},     arrSvSK_Moment[3] = {0.}
                ,      arrSvSK_ShaftForce0[3] = {0.},     arrSvSK_ShaftMoment0[3] = {0.}
                ,     arrSvSK_AirForce0[3] = {0.},     arrSvSK_AirMoment0[3] = {0.};

  //
    long double arrW[4] = {0.};


     for (int i = 0; i <quantSt; i++)
     {
           // заполнение буфера с вых информацией
           double *p = &(arrBuff[i * (QUantCurNZSKVarsVS + 4 +18)]);
           p[0]= mTimeCur;
           for (int j =0; j < (QUantCurNZSKVarsVS -1); ++j)
           {
            p[1 + j]  = marrPhaseVect[j];
           }
           p[3] = -p[3];
           //




           ///
           p[ QUantCurNZSKVarsVS    ] = arrW[1];
           p[ QUantCurNZSKVarsVS + 1] = arrW[0];
           p[ QUantCurNZSKVarsVS + 2] = arrW[2];
           p[ QUantCurNZSKVarsVS + 3] =arrW[3];

           for (int j =0; j <3; ++j)
           {
               p[ QUantCurNZSKVarsVS + 4 + j] = arrSvSK_Force[j];
               p[ QUantCurNZSKVarsVS + 7 + j] = arrSvSK_Moment[j];
               p[ QUantCurNZSKVarsVS + 10 + j] = arrSvSK_ShaftForce0[j];
               p[ QUantCurNZSKVarsVS + 13 + j] = arrSvSK_ShaftMoment0[j];
               p[ QUantCurNZSKVarsVS + 16 + j] = arrSvSK_AirForce0[j];
               p[ QUantCurNZSKVarsVS + 19 + j] = arrSvSK_AirMoment0[j];
           }

           (*pquantSteps)++;

           if ((fabsl(marrPhaseVect[9]) > 1.2) || (fabsl(marrPhaseVect[10]) > 1.2))
           {
               break;
           }
           const long double VAlTend = mTimeCur + VAlModelStep;
           int irez =replaceYourSelf( VAlTend,  VAlStIntegr);
          //int irez = makeStep(arrW[1], arrW[0], arrW[2],  arrW[3],  VAlStIntegr,  arrSvSK_Force,  arrSvSK_Moment
                      //         ,   arrSvSK_ShaftForce0,  arrSvSK_ShaftMoment0
                       //        ,  arrSvSK_AirForce0,  arrSvSK_AirMoment0);
         //  mTimeCur += VAlModelStep;
          // (*pquantSteps)++;

           if  (irez < 1)
           {
               marrPhaseVect [1] = 0.;
               marrPhaseVect [4] = 0.;
               memset(&marrPhaseVect[6], 0, (QUantCurNZSKVarsVS - 1) * sizeof(long double));
               if (bFly)
               {
                break; // вертолет приземлился
               }

           }
           else
           {
               bFly = true;
           }
     }

    }

   //-------------------------------------------------------------------------------------------
 // возвращает -1 если вертолет не может взлететь, или сел на землю
 // в противном случае возвращает 1
 int TPartHelicTraj::makeStep(const long double valFi, const long  double valKappa
                           , const long  double valEtta, const long  double valDelFi
                           , const long   double dt, long  double *arrSvSK_Force, long double *arrSvSK_Moment
                          ,long  double * arrSvSK_ShaftForce0,long  double *arrSvSK_ShaftMoment0
                          ,long  double *arrSvSK_AirForce0,long  double *arrSvSK_AirMoment0)
 {
  // 1. вычисление правой части системы дифференциальных уравнений твердого тела
   long  double *arrf = (long double *)malloc((QUantCurNZSKVarsVS-1) * sizeof(long double)) ;

     memset(arrf, 0,(QUantCurNZSKVarsVS-1) * sizeof(long double)) ;

     calcRightMemberOfDifEqSystem( valFi,  valKappa,  valEtta,valDelFi, arrf, arrSvSK_Force, arrSvSK_Moment
                                   ,   arrSvSK_ShaftForce0,  arrSvSK_ShaftMoment0
                                   ,  arrSvSK_AirForce0,  arrSvSK_AirMoment0) ;
     ///

  // 2. проверка того, что вертолет может взлететь или что он приземлился
    // if (marrPhaseVect[1] <= -10. * DBL_MIN )
    // {
   //      return -1;
   //  }

     // 2. умножение вектора arrf на dt
     long double *arrf1 = (long double *)malloc((QUantCurNZSKVarsVS-1) * sizeof(long double)) ;
     memset(arrf1, 0,(QUantCurNZSKVarsVS-1) * sizeof(long double)) ;
     MatrxMultScalar( arrf, 1, (QUantCurNZSKVarsVS-1), dt,arrf1);
     ///

     // 3. прибавлените к фазовому вектору
     MtrxSumMatrx(marrPhaseVect,arrf1,1, (QUantCurNZSKVarsVS-1), arrf) ;
     memcpy(marrPhaseVect, arrf,(QUantCurNZSKVarsVS-1) * sizeof(long double)) ;


     doAirDensity();
     free(arrf);
     free(arrf1);
     return 1;
 }

 //-------------------------------------------------------------------
//формирование матрицы частных производных 12х12
 // вектор функции правой части диф уроавнений
 // по фазовому вектору
 //INPUT:
 //arrW[4] - положенгие рулей, порядок - каппа, фи, етта, дельтаФи
// OUTPUT:
//arr_df_po_dS[12x12] - матрица частных производнных
 void TPartHelicTraj::create_df_po_dS(long  double *arrW, long  double *arr_df_po_dS)
 {
     long  double arr_df_po_dST[12*12] = {0.};

     // вектор приращений для расчета разностной     производной по переменным

     long double del[12] = {0.};
     del[0]  = del[1] = del[2] = 1.;
     del[3]  = del[4] = del[5] = 0.01;
     del[6]  = del[7] = del[8] =  0.0001;
     del[9]  = del[10] = del[11] = 0.0005;

     for (int i =0; i < 12; ++i)
     {
       calc_df_po_dSj(i, del[i],arrW, &(arr_df_po_dST[i * 12]) );
     }
     MatrTransp(arr_df_po_dST, 12, 12, arr_df_po_dS);
 }

 //--------------------------------------------------------------------
 //
 void TPartHelicTraj::calc_df_po_dSj(int j,long double del,long double *arrW, long double *arr_df_po_dSj )
 {
     long double    arrf0[12] = {0.};
     long double arrSvSK_Force[3] = {0.}, arrSvSK_Moment[3] = {0.}
                                             ,  arrSvSK_ShaftForce0[3] = {0.}, arrSvSK_ShaftMoment0[3] = {0.}
                                             ,arrSvSK_AirForce0[3] = {0.}, arrSvSK_AirMoment0[3] = {0.};
     calcRightMemberOfDifEqSystem( arrW[1], arrW[0]
                                         , arrW[2], arrW[3]
                                         , arrf0, arrSvSK_Force, arrSvSK_Moment
                                        ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                        , arrSvSK_AirForce0, arrSvSK_AirMoment0);

     TPartHelicTraj OperationCur = *this;
     OperationCur.marrPhaseVect[j] += del;
     OperationCur.marrPhaseVect[12] = OperationCur.mEnvironment.calcAirDensity(OperationCur.marrPhaseVect[1]);

     long double    arrf1[12] = {0.}, arrT[12] = {0.};
     OperationCur.calcRightMemberOfDifEqSystem( arrW[1], arrW[0]
                                          , arrW[2], arrW[3]
                                         , arrf1, arrSvSK_Force, arrSvSK_Moment
                                        ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                        , arrSvSK_AirForce0, arrSvSK_AirMoment0);
     MtrxMinusMatrx(arrf1, arrf0,1, 12, arrT);
     MatrxMultScalar(arrT, 1, 12, 1./del,arr_df_po_dSj);

 }

//----------------------------------------------------------
 //
 void TPartHelicTraj::create_df_po_dW(long double *arrW, long  double *arr_df_po_dW)
 {
     long  double arr_df_po_dWT[4*12] = {0.};

     // вектор приращений для расчета разностной     производной по переменным
     long double del[4] = {0.};
     del[0]  = del[1] = del[2] = 0.0001;
     del[3]  =  0.0001;

     for (int i =0; i < 4; ++i)
     {
       calc_df_po_dWj(i, del[i],arrW,&(arr_df_po_dWT[i * 12]) );
     }
     MatrTransp(arr_df_po_dWT, 4, 12, arr_df_po_dW);
 }

 //--------------------------------------------------------------------
 //
 void TPartHelicTraj::calc_df_po_dWj(int j,long double del,long double *arrW,long double *arr_df_po_dWj )
 {
     long double    arrf0[12] = {0.};
     long double arrSvSK_Force[3] = {0.}, arrSvSK_Moment[3] = {0.}
                                             ,  arrSvSK_ShaftForce0[3] = {0.}, arrSvSK_ShaftMoment0[3] = {0.}
                                             ,arrSvSK_AirForce0[3] = {0.}, arrSvSK_AirMoment0[3] = {0.};
     calcRightMemberOfDifEqSystem( arrW[1], arrW[0]
                                         , arrW[2], arrW[3]
                                         , arrf0, arrSvSK_Force, arrSvSK_Moment
                                        ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                        , arrSvSK_AirForce0, arrSvSK_AirMoment0);

     long double arrW1[4] = {0.};
     memcpy(arrW1, arrW, 4 * sizeof(long double));
     arrW1[j] += del;
     TPartHelicTraj OperationCur = *this;


     long double    arrf1[12] = {0.}, arrT[12] = {0.};
     calcRightMemberOfDifEqSystem( arrW1[1], arrW1[0]
                                         , arrW1[2], arrW1[3]
                                         , arrf1, arrSvSK_Force, arrSvSK_Moment
                                        ,  arrSvSK_ShaftForce0, arrSvSK_ShaftMoment0
                                        , arrSvSK_AirForce0, arrSvSK_AirMoment0);
     MtrxMinusMatrx(arrf1, arrf0,1, 12, arrT);
     MatrxMultScalar(arrT, 1, 12, 1./del,arr_df_po_dWj);
 }

//----------------------------------------------------------

 void TPartHelicTraj::doAirDensity()
 {
    marrPhaseVect[12]= mEnvironment.calcAirDensity(marrPhaseVect[1]);
 }



 bool TPartHelicTraj::replaceYourSelf(const long double VAlTend, const long double VAlStIntegr)
 {
     if (VAlTend > 50.)
     {
         int iii = 0;
     }
     int iStep = (VAlTend - mTimeCur)/VAlStIntegr;
     for (int i =0; i < iStep; ++i)
     {
        makeSingleStep(  VAlStIntegr) ;
     }


     makeSingleStep(  VAlTend - mTimeCur) ;
     return true;

 }


 void TPartHelicTraj::makeSingleStep( const long double VAlStIntegr)
 {


long double arrW[4] = {0.}, arrW_[4] = {0.};
calc_W(arrW);


 long double  arrSvSK_Moment[3] = {0.}
            ,      arrSvSK_ShaftForce0[3] = {0.},     arrSvSK_ShaftMoment0[3] = {0.}
            ,     arrSvSK_AirForce0[3] = {0.},     arrSvSK_AirMoment0[3] = {0.};


 int irez = makeStep(arrW[1], arrW[0], arrW[2],  arrW[3],  VAlStIntegr, marrSvSK_Force,  arrSvSK_Moment
         ,   arrSvSK_ShaftForce0,  arrSvSK_ShaftMoment0
         ,  arrSvSK_AirForce0,  arrSvSK_AirMoment0);

 mTimeCur += VAlStIntegr;
}

//--------------------------------------------------------------
 long double TPartHelicTraj::SIGNUM(long double x)
 {
    return (x >= 0.)?1.:-1.;
 }
//-------------------------------------------------------------------


 //------------------------------------------------

 //отображение лемнискаты на круг радиуса 1/2
 //z1, z2 - полюсы лемнискаты
 // val_A = a*a -стоит в правой части уравнения лемн6искаты |(z-z1)*(z-z2)| = a*a
 // уравнение леминискаты записыватся в виде |(z-z1)*(z-z2)|/(a*a) = 1
 // нам надо отобразить лемнискату на круг 1/2, стало быть,
 // |(z-z1)*(z-z2)|/(a*a) = 1/2
 // тогда конформное преобразование будет иметь вид
 // om = 1/A *z*z -(z1 + z2)/A * z + z1*z2/A
 // где A = a*a
 bool  TPartHelicTraj::IsStability(const long double z1 ,const long double  z2
                                  ,const long double val_A,long double*arrA,const  int dimArrA)
 {
     bool breturn = false;
    double *arrA0 = new double [dimArrA * dimArrA];
    for (int i = 0; i <dimArrA * dimArrA; ++i)
    {
        arrA0[i] = arrA[i] ;
    }

    double a = 1./ val_A;
    double b = -(z1 + z2)/ val_A;
    double c = z1 * z2/ val_A ;
    double *arrB = new double [dimArrA * dimArrA];
    double *arrT = new double [dimArrA * dimArrA];
    double *arrT0 = new double [dimArrA * dimArrA];
    caclMtrxPolyn2(  a,  b,c ,arrA0, dimArrA, arrB);
    memcpy(arrT, arrB, dimArrA * dimArrA * sizeof(double) );
    int i =0;
    double valSp =1000.;
    for ( i = 0; i < 60; ++i)
    {
      MtrxMultMatrx(arrB,dimArrA, dimArrA, arrT,dimArrA, arrT0) ;
      valSp =  fabs(Sp( arrT0, dimArrA));
      if (valSp <= FLT_MIN * 2.)
      {
          breturn = true;
          break;
      }

      memcpy(arrB, arrT0, dimArrA * dimArrA * sizeof(double) );
      memcpy(arrT, arrT0, dimArrA * dimArrA * sizeof(double) );
    }

    delete []arrA0;
    delete []arrB;
    delete []arrT;
    delete []arrT0;
    return breturn;
 }



 bool  TPartHelicTraj::IsStability_(const long double z1 ,const long double  z2
                                  ,const long double val_A,long double*arrA,const  int dimArrA)
 {
     bool breturn = false;
    float *arrA0 = new float [dimArrA * dimArrA];
    for (int i = 0; i <dimArrA * dimArrA; ++i)
    {
        arrA0[i] = arrA[i] ;
    }

    float a = 1./ val_A;
    float b = -(z1 + z2)/ val_A;
    float c = z1 * z2/ val_A;
    float *arrB = new float [dimArrA * dimArrA];
    float *arrT = new float [dimArrA * dimArrA];
    float *arrT0 = new float [dimArrA * dimArrA];
    caclMtrxPolyn2(  a,  b,c ,arrA0, dimArrA, arrB);
    memcpy(arrT, arrB, dimArrA * dimArrA * sizeof(float) );
    int i =0;
    float valSp =1000.;
    int iarr[3] = {0};
    for ( i = 0; i < 30; ++i)
    {
      MtrxMultMatrx(arrB,dimArrA, dimArrA, arrT,dimArrA, arrT0) ;
       // MtrxMultMatrxTransp(arrB,dimArrA, dimArrA, arrT,dimArrA, arrT0) ;
      valSp =  fabs(Sp( arrT0, dimArrA));


      if (valSp > sqrt(FLT_MAX))
      {
          break;
      }
      if (valSp < FLT_MIN )
      {
          int ind = iarr[0]+ iarr[1];
          if (ind ==2)
          {
          breturn = true;
          break;
          }
          else
           {
             if (iarr[0] == 1)
             {
              iarr[1] = 1.;
             }
             else
             {
              iarr[0] =1;
             }
            }
      }
      else
      {

       iarr[0] = iarr[1] = 0;

      }

      memcpy(arrB, arrT0, dimArrA * dimArrA * sizeof(float) );
      memcpy(arrT, arrT0, dimArrA * dimArrA * sizeof(float) );
    }

    delete []arrA0;
    delete []arrB;
    delete []arrT;
    delete []arrT0;
    return breturn;
 }

 // Требуетс    определить удовлтрятли корни хаактрстичессого урвнения матицы
 // arrA условиям: < Min< Re<Max и fabs(Im/Re) < MaxTang
 bool  TPartHelicTraj::IsRootsSuit_dim5(long double*arrA,const  int dimA, const  double VAlMinRe
                     ,const  double VAlMaxRe, const  double VAlMaxTang, TComp *cmparrRoots)
 {
     double arrA0[25] = {0.};
     for (int i =0;i < 25; ++i)
     {
       arrA0[i] = arrA[i];
     }
    if( !solvCharactEq_5(arrA0, cmparrRoots,VAlMinRe, VAlMaxRe))
    {
        return false;
    }

     for(int i =0; i < 5; ++i)
     {
         if(( cmparrRoots[i].m_Re - VAlMinRe) * ( cmparrRoots[i].m_Re - VAlMaxRe) > 0.)
         {
             return false;
         }
         if(fabs(cmparrRoots[i].m_Im/ cmparrRoots[i].m_Re)> VAlMaxTang)
         {
             return false;
         }
     }
     return true;
 }

 //---------------------------------------------------------

/*
 void TPartHelicTraj::fill_df_po_px_and_df_po_dW_LittleTask(const long double VAlY ,   const long double VAlVx
                                                           , long double *arr_dF_po_dx, long double *arr_dF_po_dW)
 {

     TPartHelicTraj OperationWork = *this;

     const long double VAlNu0 = -0.01;
     long double arrCntrRulsPos[4] = {0.}, arrCntrImpacts[4] = {0.};
     long double valNu = 10000.;
     OperationWork.calcBalanceParams_ForSteadyState_Move( VAlVx,   VAlY
                         ,  VAlNu0 ,arrCntrRulsPos
                         ,arrCntrImpacts,&valNu) ;


     memset(OperationWork.marrPhaseVect, 0, (QUantCurNZSKVarsVS -1) * sizeof(long double));
     OperationWork.marrPhaseVect[1] = VAlY;
     OperationWork.marrPhaseVect[3] = VAlVx;
     OperationWork.marrPhaseVect[10] = valNu;
     OperationWork.marrPhaseVect[12] = OperationWork.mEnvironment.calcAirDensity(OperationWork.marrPhaseVect[1]);

     // поиск передаточных чисел виртуального летчика

     long double valT_po_FI = 2. * OperationWork.mHelic.marrRotor[0].calc_Deriv_T_po_Fi(OperationWork.marrPhaseVect[12]);

     // заполнение матрицы dF_po_dx = A0
     memset(arr_dF_po_dx, 0, 25 * sizeof(long double));

     arr_dF_po_dx[7] = 1.;
     arr_dF_po_dx[19] = 1.;

     // массив номеров переменных в упрощенной задаче в сиситеме нумерации плоной задачи
     int iarrNums[5] = {3, 1, 4, 10, 8};

     // массив разностей для вычисления разностных производных аэродинм силы и момента
     long double arrDel[5] = {0.01, 10.,0.01,0.00005, 0.00001};

     long double arrFa[3] = {0.}, arrMa[3]= {0.}, arr_dFa_po_dxj[3]= {0.},arr_dMa_po_dxj[3]= {0.};

     for (int i = 0; i < 5; ++i)
     {
         OperationWork.calc_dFa_and_dMa_po_dxj( arrDel[i], iarrNums[i],arrFa ,arrMa, arr_dFa_po_dxj,arr_dMa_po_dxj);
         arr_dF_po_dx[i] = (arr_dFa_po_dxj[0] * cos(valNu)- arr_dFa_po_dxj[1]* sin(valNu))
                 / OperationWork.mHelic.mHelicMass;

         arr_dF_po_dx[10 + i] = (arr_dFa_po_dxj[0] * sin(valNu) + arr_dFa_po_dxj[1]* cos(valNu))
                 / OperationWork.mHelic.mHelicMass;

         arr_dF_po_dx[20 + i] = arr_dMa_po_dxj[2] / OperationWork.mHelic.marrJ[8] ;
     }
     arr_dF_po_dx[3] -= G_ZEMLI;
 ///
     // заполнение матрицы dF_po_dW = B0
   //  arrCntrRulsPos[4] - Fi, Kappa, Etta, Delfi
     long double valW1 = arrCntrRulsPos[1];// KAPPA
     long double valW2 = arrCntrRulsPos[0];// Fi
      memset(arr_dF_po_dW, 0, 10 * sizeof(long double));

     long double aggt = valNu - OperationWork.mHelic.marrRotor[0].mZaklinAng;
     long double valXcm = OperationWork.mHelic.marrRotor[0].mForceArmX;
     long double valYcm = (OperationWork.mHelic.marrRotor[0].mForceArmY + OperationWork.mHelic.marrRotor[1].mForceArmY) /2.;;


     arr_dF_po_dW[0] = valW2 * cos(aggt) / OperationWork.mHelic.mHelicMass;
     arr_dF_po_dW[1] = (valW1 * cos(aggt) - sin(aggt))/ OperationWork.mHelic.mHelicMass;

     arr_dF_po_dW[4] = valW2 * sin(aggt) / OperationWork.mHelic.mHelicMass;
     arr_dF_po_dW[5] = (valW1 * sin(aggt) + cos(aggt))/ OperationWork.mHelic.mHelicMass;

     arr_dF_po_dW[8] = -valYcm/OperationWork.mHelic.marrJ[8] * valW2;
     arr_dF_po_dW[9] = (valXcm - valYcm * valW1)/OperationWork.mHelic.marrJ[8] ;
     MatrxMultScalar(arr_dF_po_dW, 10, 1, valT_po_FI, arr_dF_po_dW);


 }
*/
/*
 void TPartHelicTraj::fill_df_po_px_and_df_po_dW_BigTask(const long double VAlY ,   const long double VAlVx
                                                        , long double *arr_dF_po_dx, long double *arr_dF_po_dW)
 {
     TPartHelicTraj OperationWork =  *this;
     memset(OperationWork.marrPhaseVect, 0, sizeof(double)* QUantCurNZSKVarsVS);
     OperationWork.marrPhaseVect[1] = VAlY;
     OperationWork.marrPhaseVect[3] = VAlVx;
     OperationWork.marrPhaseVect[12]= OperationWork.mEnvironment.calcAirDensity(VAlY);

     const long double VAlNu0 = -0.01;
     long double arrCntrRulsPos[4] = {0.}, arrCntrImpacts[4] = {0.};
     long double valNu = 10000.;
     OperationWork.calcBalanceParams_ForSteadyState_Move( VAlVx,   VAlY
                         ,  VAlNu0 ,arrCntrRulsPos
                         ,arrCntrImpacts,&valNu) ;

     memset(OperationWork.marrPhaseVect, 0, (QUantCurNZSKVarsVS -1) * sizeof(long double));
     OperationWork.marrPhaseVect[1] = VAlY;
     OperationWork.marrPhaseVect[3] = VAlVx;
     OperationWork.marrPhaseVect[10] = valNu;
     OperationWork.marrPhaseVect[12] = OperationWork.mEnvironment.calcAirDensity(OperationWork.marrPhaseVect[1]);
     //  arrCntrRulsPos[4] - Fi, Kappa, Etta, Delfi
       long double valW1 = arrCntrRulsPos[1];// KAPPA
       long double valW2 = arrCntrRulsPos[0];// Fi

       long double aggt = valNu - OperationWork.mHelic.marrRotor[0].mZaklinAng;
       long double valXcm = OperationWork.mHelic.marrRotor[0].mForceArmX;
       long double valYcm = (OperationWork.mHelic.marrRotor[0].mForceArmY + OperationWork.mHelic.marrRotor[1].mForceArmY) /2.;

        long double valT_po_FI = 2. * OperationWork.mHelic.marrRotor[0].calc_Deriv_T_po_Fi(OperationWork.marrPhaseVect[12]);


     // ФОРМИРОВАНИЕ МАТРИЦ ДЛЯ ПОЛЛНОЙ ЗАДАЧИ
       // вектор опорного управления
     long double arrW[4] = {0.};
     arrW[0] = arrCntrRulsPos[1];  // каппа
     arrW[1] = arrCntrRulsPos[0];  // фи
      // матрица df_po_dS
     long double arr_df_po_dS[12 * 12] = {0.};
     OperationWork.create_df_po_dS(arrW, arr_df_po_dS);
     // выделение матрицы 11х11

     for (int i = 0; i < 11; ++i)
          for (int j =0; j< 11;++j)
     {
        arr_dF_po_dx[ i * 11 + j] = arr_df_po_dS[ (i + 1) * 12 + j + 1] ;
     }
     // матрица df_po_dW
     long double arr_df_po_dW[12*4] = {0.};
     OperationWork.create_df_po_dW(arrW, arr_df_po_dW);

     // выделение матрицы 11х4

     memcpy(arr_dF_po_dW, &(arr_df_po_dW[4]), 11 * 4 * sizeof(long double));


 }
*/

 //-----------------------------------------------------------------------------------
  // вираж безобратных связей
  // на время полета VAlFlyTime
  // INPUT:
  // VAlFlyTime - время полета
  // VAlStIntegr - шаг интегрирования
  //OUTPUT:
  // arrBuff - буфер (массив) для накопления информации
  // Если arrBuff==NULL то информация не накапливается
  // *pquantSteps - к-во шагов интегрирования

  // вектор выходной информации arrOutputVect[25]
 // составляющие воздушной скорости в связанной системе координат; 3
 // составляющие путевой скорости в связанной системе координат; 3
 // приборная скорость;                                          1
 // бароинерциальная вертикальная скорость;                      1
 // геометрическая высота;                                       1
 // абсолютная барометрическая высота;                           1
 // относительная барометрическая высота;                         1
 // геодезическая высота;                                        1
 // бароинерциальная высота;                                     1
 // угол сноса;                                                  1
 // угол атаки;                                                  1
 // угол скольжения;                                              1
 // скорость и направление ветра;                                2
 // крен;                                                        1
 // тангаж;                                                      1
 // курс истинный;                                                1
  //магнитное склонение;                                          1
 // температура наружного воздуха;                               1
 // текущие координаты местоположения объекта в географической системе координат. 2

  bool TPartHelicTraj::TryHorizontalTurn(const long double VAlR,const long double VAlPsi,const long double VAlFlyTime, const long double VAlStIntegr
                           , double *arrBuff, int *pquantSteps)
  {
      // признак полета
      bool bFly = false;
      ///

      *pquantSteps = 0;

   // управляющие воздействия
      // 1.Общий шаг
      long double valFi = 0. ;
      // 2. Ход ручки ППУ по тангажу
      long double valKappa = 0.;
      // 3. Ход ручки ППУ по крену
      long double valEtta = 0.;
      // 4. Ход педалей
      long double valDelFi = 0.;
      long double arrCntrRulsPos[4] = {0.}, arrCntrImpacts[4] = {0.};
      long double valNu = 10000.;
     // calcBalanceParams_ForSteadyState_Move(marrPhaseVect[3], marrPhaseVect[1]
      //                    , 0.01 , arrCntrRulsPos
       //                   ,arrCntrImpacts,&valNu );


      long double arrX0[6] = {0.,0.,0., 0.3,0.,0.};
      long double arrXRez [6] = {0.};
      if(!TBallanceCalc::calcBallParamsForTurn(mHelic,mEnvironment
                                                            ,marrPhaseVect[1],marrPhaseVect[3],VAlPsi
                                                            , VAlR,arrX0, arrXRez ))
      {
              return false;
     }
      marrPhaseVect[6] = VAlR/marrPhaseVect[3] * sinl(arrXRez[1]);
      marrPhaseVect[7] = VAlR/marrPhaseVect[3] * cosl(arrXRez[1]) * cosl(arrXRez[1]);
      marrPhaseVect[8] = -VAlR/marrPhaseVect[3] * cosl(arrXRez[1]) * sinl(arrXRez[1]);

      marrPhaseVect[9] = arrXRez[0];
      marrPhaseVect[10] = arrXRez[1];
      marrPhaseVect[11] = VAlPsi;
      valFi = arrXRez[3];
      valKappa = arrXRez[2];
      valEtta = arrXRez[4];
      valDelFi = arrXRez[5];

   int quantSt = (int)(VAlFlyTime/ VAlStIntegr);
   long double arrSvSK_Force[3] = {0.},     arrSvSK_Moment[3] = {0.}
              ,      arrSvSK_ShaftForce0[3] = {0.},     arrSvSK_ShaftMoment0[3] = {0.}
              ,     arrSvSK_AirForce0[3] = {0.},     arrSvSK_AirMoment0[3] = {0.};



   for (int i = 0; i <quantSt; i++) {


         // заполнение буфера с вых информацией
         double *p = &(arrBuff[i * 17]);
         p[0]= mTimeCur;
         for (int j =0; j < (QUantCurNZSKVarsVS -1); ++j)
         {
          p[1 + j]  = marrPhaseVect[j];
         }
         p[ QUantCurNZSKVarsVS    ] = valFi;
         p[ QUantCurNZSKVarsVS + 1] = valKappa;
         p[ QUantCurNZSKVarsVS + 2] = valEtta ;
         p[ QUantCurNZSKVarsVS + 3] = valDelFi;

        /* for (int j =0; j <3; ++j)
         {
             p[ QUantCurNZSKVarsVS + 4 + j] = arrSvSK_Force[j];
             p[ QUantCurNZSKVarsVS + 7 + j] = arrSvSK_Moment[j];
             p[ QUantCurNZSKVarsVS + 10 + j] = arrSvSK_ShaftForce0[j];
             p[ QUantCurNZSKVarsVS + 13 + j] = arrSvSK_ShaftMoment0[j];
             p[ QUantCurNZSKVarsVS + 16 + j] = arrSvSK_AirForce0[j];
             p[ QUantCurNZSKVarsVS + 19 + j] = arrSvSK_AirMoment0[j];
         }
*/

         int irez = makeStep(valFi, valKappa, valEtta,  valDelFi,  VAlStIntegr,  arrSvSK_Force,  arrSvSK_Moment
                             ,   arrSvSK_ShaftForce0,  arrSvSK_ShaftMoment0
                             ,  arrSvSK_AirForce0,  arrSvSK_AirMoment0);
         mTimeCur += VAlStIntegr;
        // (*pquantSteps)++;

         if  (irez < 1)
         {
             marrPhaseVect [1] = 0.;
             marrPhaseVect [4] = 0.;
             memset(&marrPhaseVect[6], 0, (QUantCurNZSKVarsVS - 1) * sizeof(long double));
             if (bFly)
             {
              break; // вертолет приземлился
             }
             if((marrPhaseVect[4] <= 0.) && (i > 1))
             {
              int iii = 0;
              break;
             }
         }
         else
         {
             bFly = true;
         }
        ( *pquantSteps)++;
   }
return true;
  }


  //-------------------------------------------------------------------------------------------
// заполнение матриц большой линеаризованной задачи для виража
  // INPUT:
  // VAlRad - радиус виража
  // arrW [] - вектор положения рулей, - Kappa,Fi,  Etta, Delfi
  //OUTPUT:
  // arr_dF_po_dx[100] - матрица частных производных по фазовому вектору
  // arr_dF_po_dW [10 * 4] - матрица частных производных по управлениям
  // в фазовом векторе хранятся остальные параметры - высота, скорость, угол рыскания Psi
  void TPartHelicTraj::fill_df_po_px_and_df_po_dW_Turn(const long double VAlRad, long double *arrW
                                                         , long double *arr_dF_po_dx, long double *arr_dF_po_dW)
  {
      doAirDensity();
       // матрица df_po_dS
      long double arr_df_po_dS[12 * 12] = {0.};
      create_df_po_dS(arrW, arr_df_po_dS);
      // выделение матрицы 11x11
      // номера пременных в фазовом векторе:
      // 1 - высота
      // 2 - Z
      // 3 - Vx
      // 4 - Vy
      // 5 - Vz
      // 6 - Omx
      // 7 - Omy
      // 8 - Omz
      // 9 - Gamma (крен)
      // 10 - NU (тангаж)
      // 11 - Psi (рыскание)

      // массив номеров фазового вектора
      int arrNUms[11 ] = {1, 2,3, 4, 5, 6, 7, 8, 9, 10, 11};

      for (int i = 0; i < 11; ++i)
      {
          int numi = arrNUms[i] ;
      for (int j =0; j< 11;++j)
      {
          int numj = arrNUms[j] ;
         arr_dF_po_dx[ i * 11 + j] = arr_df_po_dS[ numi * 12 + numj] ;
      }
      }
      // матрица df_po_dW
      long double arr_df_po_dW[12*4] = {0.};
      create_df_po_dW(arrW, arr_df_po_dW);

      // выделение матрицы 10х4

      for (int i = 0; i < 11; ++i)
      {
       int numi = arrNUms[i] ;
       memcpy(&( arr_dF_po_dW[ i * 4]), &(arr_df_po_dW[ numi * 4]), 4 * sizeof(long double));

      }

  }


  //-----------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------
// заполнение матриц малой  линеаризованной задачи для виража
  // по каналу тангажа
  // INPUT:
  // VAlRad - радиус виража
  // arrW [] - вектор положения рулей, - Kappa,Fi,  Etta, Delfi
  //OUTPUT:
  // arr_dF_po_dx[25] - матрица частных производных по фазовому вектору
  // arr_dF_po_dW [5 * 2] - матрица частных производных по управлениям
  // в фазовом векторе хранятся остальные параметры - высота, скорость, угол рыскания Psi
  void TPartHelicTraj::fill_df_po_px_and_df_po_dW_Turn_Little(const long double VAlRad, long double *arrW
                                                         , long double *arr_dF_po_dx, long double *arr_dF_po_dW)
  {
      doAirDensity();
       // матрица df_po_dS
      long double arr_df_po_dS[12 * 12] = {0.};
      create_df_po_dS(arrW, arr_df_po_dS);
      // выделение матрицы 5x5
      // номера пременных в фазовом векторе:
      // 3 - Vx
      // 1 - высота Y
      // 4 - Vy
      // 10 - NU (тангаж)
      // 8 - Omz



      // массив номеров фазового вектора
      int iarrNums[5] = {3, 1, 4, 10, 8};

      for (int i = 0; i < 5; ++i)
      {
          int numi = iarrNums[i] ;
      for (int j =0; j< 5;++j)
      {
          int numj = iarrNums[j] ;
         arr_dF_po_dx[ i * 5 + j] = arr_df_po_dS[ numi * 12 + numj] ;
      }
      }
      // матрица df_po_dW
      long double arr_df_po_dW[12*4] = {0.};
      create_df_po_dW(arrW, arr_df_po_dW);

      // выделение матрицы 5х2

      for (int i = 0; i < 5; ++i)
      {
       int numi = iarrNums[i] ;
       memcpy(&( arr_dF_po_dW[ i * 2]), &(arr_df_po_dW[ numi * 4]), 2 * sizeof(long double));

      }

  }


  //-----------------------------------------------------------------------------------
  // вычисление вектора угловых скоростей вращения в СвСК
  // для стационарного виража с радиусом VAlRad
  // INPUT:
  // VAlV - модуль скорости
  // VAlRad - радиус
  // VAlGamma - угол крена
  // VAlNu - угол тангажа
  //OUTPUT:
  // arrOmega[3] - вектор угловой скорости в СвСК
  void TPartHelicTraj::fillUpOmega_StableTurn(const long double VAlV,const long double VAlRad,const long double VAlGamma
                                              ,const long double VAlNu, long double *arrOmega)
  {

    long double val_dPsi_po_dt = VAlV / VAlRad;
    arrOmega[0] = val_dPsi_po_dt * sinl(VAlNu);
    arrOmega[1] = val_dPsi_po_dt * cosl(VAlNu) * cosl(VAlGamma);
    arrOmega[2] = -val_dPsi_po_dt * cosl(VAlNu) * sinl(VAlGamma);

  }

  void TPartHelicTraj::calc_W(long double *arrW)

  {
 ;
  }


  void TPartHelicTraj::fill_df_po_px_and_df_po_dW_TangCanale(long double *arr_dF_po_dx, long double *arr_dF_po_dW)
  {


          doAirDensity();
           // матрица df_po_dS
          long double arr_df_po_dS[(QUantCurNZSKVarsVS -1) * (QUantCurNZSKVarsVS -1)] = {0.};
          long double arrW[4] = {0.};
          get_arrSteadyW(arrW);
          create_df_po_dS(arrW, arr_df_po_dS);
          // выделение матрицы 5x5
          // номера пременных в фазовом векторе:
          // 3 - Vx
          // 1 - высота Y
          // 4 - Vy
          // 10 - NU (тангаж)
          // 8 - Omz


          int iarrNums[QUantCurNZSKVarsVS -1] = {0};
          int iNumVars = -1;
          get_arrayOfControlledVarsTang(&iNumVars , iarrNums);


          for (int i = 0; i < iNumVars; ++i)
          {
              int numi = iarrNums[i] ;
          for (int j =0; j< iNumVars;++j)
          {
              int numj = iarrNums[j] ;
             arr_dF_po_dx[ i * iNumVars + j] = arr_df_po_dS[ numi * (QUantCurNZSKVarsVS -1) + numj] ;
          }
          }
          // матрица df_po_dW
          long double arr_df_po_dW[(QUantCurNZSKVarsVS -1)*4] = {0.};
          create_df_po_dW(arrW, arr_df_po_dW);

          // выделение матрицы iNumVars х 2

          for (int i = 0; i < iNumVars; ++i)
          {
           int numi = iarrNums[i] ;
           memcpy(&( arr_dF_po_dW[ i * 2]), &(arr_df_po_dW[ numi * 4]), 2 * sizeof(long double));

          }

  }


  void TPartHelicTraj::get_arrSteadyW(long double *arrW)
  {

  }

  //-------------------------------------------------------------------------------------------
// заполнение матриц большой линеаризованной задачи для виража
  // INPUT:
  // VAlRad - радиус виража
  // arrW [] - вектор положения рулей, - Kappa,Fi,  Etta, Delfi
  //OUTPUT:
  // arr_dF_po_dx[100] - матрица частных производных по фазовому вектору
  // arr_dF_po_dW [10 * 4] - матрица частных производных по управлениям
  // в фазовом векторе хранятся остальные параметры - высота, скорость, угол рыскания Psi
  void TPartHelicTraj::fill_df_po_px_and_df_po_dW( long double *arr_dF_po_dx, long double *arr_dF_po_dW)
  {

      doAirDensity();
       // матрица df_po_dS
      long double arr_df_po_dS[(QUantCurNZSKVarsVS -1) * (QUantCurNZSKVarsVS -1)] = {0.};
      long double arrW[4] = {0.};
      get_arrSteadyW(arrW);
      create_df_po_dS(arrW, arr_df_po_dS);
      // выделение матрицы 5x5
      // номера пременных в фазовом векторе:
      // 3 - Vx
      // 1 - высота Y
      // 4 - Vy
      // 10 - NU (тангаж)
      // 8 - Omz


      int iarrNums[QUantCurNZSKVarsVS -1] = {0};
      int iNumVars = -1;
      get_arrayOfControlledVars(&iNumVars , iarrNums);


      for (int i = 0; i < iNumVars; ++i)
      {
          int numi = iarrNums[i] ;
      for (int j =0; j< iNumVars;++j)
      {
          int numj = iarrNums[j] ;
         arr_dF_po_dx[ i * iNumVars + j] = arr_df_po_dS[ numi * (QUantCurNZSKVarsVS -1) + numj] ;
      }
      }
      // матрица df_po_dW
      long double arr_df_po_dW[(QUantCurNZSKVarsVS -1)*4] = {0.};
      create_df_po_dW(arrW, arr_df_po_dW);

      // выделение матрицы iNumVars х 4

      for (int i = 0; i < iNumVars; ++i)
      {
       int numi = iarrNums[i] ;
       memcpy(&( arr_dF_po_dW[ i * 4]), &(arr_df_po_dW[ numi * 4]), 4 * sizeof(long double));

      }
  }

//---------------------------------------------------------------------------------
  void TPartHelicTraj::get_arrayOfControlledVars(int *piNum, int *iarr)
  {

  }


  void TPartHelicTraj::get_arrayOfControlledVarsTang(int *piNum, int *iarr)
  {

  }
  //----------------------------------------------------------------------------

    void TPartHelicTraj::get_QuantOfControlledVarsTang(int *piNumr)
    {


    }

    //----------------------------------------------------------------------------

    void TPartHelicTraj::get_QuantOfControlledVars(int *piNumr)
    {

    }

 //----------------------------------------------------------------------------
  void TPartHelicTraj::selectGearsTang(long double valz1 ,  long double valz2, long double val_A
                                       ,double *arrDataBuf ,const int maxQuant, int *quantRows)
  {
      int quantControlledVars = -1;
      get_QuantOfControlledVarsTang(&quantControlledVars);
      long double *arr_dF_po_dx = new long double[quantControlledVars * quantControlledVars] ;
      long double *arr_dF_po_dW = new long double[quantControlledVars * 2] ;
      memset(arr_dF_po_dx, 0, quantControlledVars * quantControlledVars * sizeof(long double));
      memset(arr_dF_po_dW, 0, 2* quantControlledVars * sizeof(long double));

        fill_df_po_px_and_df_po_dW_TangCanale(arr_dF_po_dx, arr_dF_po_dW);
        //


      // поиск передаточных чисел виртуального летчика



   memset(arrDataBuf, 0, maxQuant * (QUantCurNZSKVarsVS -1) *4 * sizeof( double));

   findingCircleTang(valz1 , valz2,  val_A,arr_dF_po_dx, arr_dF_po_dW ,arrDataBuf ,maxQuant, quantRows);

    delete []arr_dF_po_dx;
    delete []arr_dF_po_dW;
  }


  //-----------------------------------------------------------------------------------
  void TPartHelicTraj::findingCircleTang(long double valz1 ,  long double valz2, long double val_A, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                                         ,double *arrDataBuf ,const int maxQuant, int *pquantRows)
  {

  }
  //вычисление положения рулей для упрощенной модели винта
  //
  //INPUT:
  //arrCntrImpacts[4] - силы H,T, S и момент винтов
  //OUTPUT:
  // - Fi, Kappa, Etta, DelFi
  void TPartHelicTraj::calcCntrRulsPos_New(long double *arrCntrImpacts,long  double *pvalFi
                    ,long  double *pvalKappa ,long  double *pvalEtta,long  double *pvalDelFi)
  {
      // вычисление плотности атмосферы
     //long  double valAtmRo = mEnvironment.calcAirDensity(marrPhaseVect[1]);

     long  double arrUaSvSK[3] ={0.};
     // calcUaSvSK_Case_SvSK(arrUaSvSK);
     calcUaSvSK_Case_NZSK(arrUaSvSK);
     long  double valSquare =  M_PI * mHelic.marrRotor[0].mBlade.mBladeR * mHelic.marrRotor[0].mBlade.mBladeR;
     long  double temp = mHelic.marrRotor[0].mCt *marrPhaseVect[12] * valSquare
                   * (mHelic.marrRotor[0].mOmega * mHelic.marrRotor[0].mBlade.mBladeR )* (mHelic.marrRotor[0].mOmega * mHelic.marrRotor[0].mBlade.mBladeR );
       *pvalFi =  arrCntrImpacts[1]/ temp;
       if ((*pvalFi) > mHelic.mFiMax)
       {
          // *pvalFi = mHelic.mFiMax;
       }

       // 8. ВЫЧИСЛЕНИЕ УГЛА ПЕДАЛЕЙ valDelFi
      long  double Va = sqrtl( arrUaSvSK[0] * arrUaSvSK[0] +  arrUaSvSK[2] * arrUaSvSK[2]);
      long  double valQp = mHelic.marrRotor[0].calcQreact( marrPhaseVect[12],  Va);

      *pvalDelFi = arrCntrImpacts[3]/ valQp/ 2.;

      *pvalKappa = arrCntrImpacts[0]/ arrCntrImpacts[1];
      *pvalEtta  = arrCntrImpacts[2]/ arrCntrImpacts[1];
  }
  //----------------------------------------------------------------
void TPartHelicTraj::insertLocalGearMtrx_In_TotalGearMtrx(const int QuantControlledVars, int *iarrNumsControlledVars
                                                              ,const int QuantControls, int *iarrNumsControls
                                                              ,long double *arrLocalGearMtrx, double *arrRezTotalGearMtrx)
{
  memset (arrRezTotalGearMtrx, 0, 4 *(QUantCurNZSKVarsVS-1)  * sizeof(double) );
  for (int i =0; i < QuantControls; ++i)
  {
      int numTotalRow = iarrNumsControls[i];
      for(int j =0; j <QuantControlledVars; ++j )
      {
         int numTotalCol = iarrNumsControlledVars[j];
         arrRezTotalGearMtrx[numTotalRow * (QUantCurNZSKVarsVS-1) + numTotalCol] = arrLocalGearMtrx[ i *QuantControlledVars + j];
      }
  }
}


//----------------------------------------------------------------
void TPartHelicTraj::extractLocalGearMtrx_From_TotalGearMtrx(double *arrInpTotalGearMtrx,const int QuantControlledVars, int *iarrNumsControlledVars
                                                            ,const int QuantControls, int *iarrNumsControls
                                                            ,long double *arrOutLocalGearMtrx )
{

for (int i =0; i < QuantControls; ++i)
{
    int numTotalRow = iarrNumsControls[i];
    for(int j =0; j <QuantControlledVars; ++j )
    {
       int numTotalCol = iarrNumsControlledVars[j];
       arrOutLocalGearMtrx[ i *QuantControlledVars + j] = arrInpTotalGearMtrx[numTotalRow * (QUantCurNZSKVarsVS-1) + numTotalCol];
    }
}
}
//--------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------
 void TPartHelicTraj::selectGears(long double valz1 ,  long double valz2, long double val_A
                                 ,  double *arrC, double *arrDataBuf ,const int maxQuant, int *quantRows)
 {
     int quantControlledVars = -1;
     int iarrNumsControlledVars[QUantCurNZSKVarsVS -1] = {0};
     get_arrayOfControlledVars(&quantControlledVars, iarrNumsControlledVars);

     long double *arr_dF_po_dx = new long double[quantControlledVars * quantControlledVars] ;
     long double *arr_dF_po_dW = new long double[quantControlledVars * 4] ;
     memset(arr_dF_po_dx, 0, quantControlledVars * quantControlledVars * sizeof(long double));
     memset(arr_dF_po_dW, 0, 4* quantControlledVars * sizeof(long double));

       fill_df_po_px_and_df_po_dW(arr_dF_po_dx, arr_dF_po_dW);
       //
       int iarrNumsControls[4] = {0,1,2,3};
       const int QuantControls =4;
       long double arrC00[(QUantCurNZSKVarsVS -1) *4] = {0.};
       extractLocalGearMtrx_From_TotalGearMtrx(arrC, quantControlledVars, iarrNumsControlledVars
                                                                              , QuantControls, iarrNumsControls
                                                                              ,arrC00);



     // поиск передаточных чисел виртуального летчика



  memset(arrDataBuf, 0, maxQuant * (QUantCurNZSKVarsVS -1) * 4 *sizeof( double));

  findingCircle(valz1 , valz2,  val_A, arrC00, arr_dF_po_dx, arr_dF_po_dW ,arrDataBuf ,maxQuant, quantRows);

   delete []arr_dF_po_dx;
   delete []arr_dF_po_dW;
 }


 //-----------------------------------------------------------------------------------

 void TPartHelicTraj::findingCircle(long double valz1 ,  long double valz2, long double val_A
                   , long double *arrC, long double *arr_dF_po_dx, long double *arr_dF_po_dW
                  ,double *arrDataBuf ,const int maxQuant, int *pquantRows)
 {
 }
/*
  //формирование буфера для построения графика функции стоящей
  // в правой части уравнения относительно ню
  // это относится к нахождению параметров стабилизированнонго
  //горизонтального движения с постоянной скоростью на заданной высоте
  // INPUT:
  //VAlVx - скорость по X
  //VAlY - высота
  //VAlRangNu
  //QUantPoints - к-во точек
  //OUTPUT:
  // arrBuff1[QUantPoints *2] - массив с выходной информацией
  void TPartHelicTraj::calcGraph_FGr_ForSteadyState_Move(const long double VAlVx, const long double  VAlY
                                    ,const long double   VAlRangNu, const int QUantPoints, double *arrBuff)
  {
      TPartHelicTraj PartHelicTrajCur = *this;
      memset(PartHelicTrajCur.marrPhaseVect, 0, QUantCurNZSKVarsVS * sizeof(long double));
      PartHelicTrajCur.marrPhaseVect[1] = VAlY;
      PartHelicTrajCur.marrPhaseVect[3] = VAlVx;
      PartHelicTrajCur.doAirDensity();

      double valStep = 2. *VAlRangNu/ ((double)QUantPoints - 1.);
      long double valFa = 0.,valMa = 0.;
      long double arrAirForce[3] = {0.}, arrAirMom[3] = {0.};
      for (int i = 0; i < QUantPoints; ++i)
      {
        arrBuff[ i * 4    ] = -VAlRangNu + ((double) i) *  valStep;
        arrBuff[ i * 4 + 1] = PartHelicTrajCur.fncFGr_HorizBalance(arrBuff[ i * 4], &valFa, &valMa, arrAirForce, arrAirMom);
        arrBuff[ i * 4 + 2] = valFa;
        arrBuff[ i * 4 + 3] = valMa;

      }
  }
*/
/*
  // вычисление функции вертикальногшо балланса
  // для уппрощенной модели винта+-
  //INPUT:
  //VAlNu - угол тангажа
  //OUTPUT:
  //*pvalFRotX - сила винта по оси X БСК
  //*pvalFRotY - сила винта по оси Y БСК (сила тяги T)
  long double TPartHelicTraj::fncFGr_HorizBalance(const long double VAlNu
                           ,long  double *pvalFRotX,long  double *pvalFRotY
                           ,long  double *arrSvSK_AirForce,long  double * arrSvSK_AirMoment)
  {

      //
      TPartHelicTraj OperationCur = *this;
      OperationCur.marrPhaseVect[10] = VAlNu;
      // вычисление плотности атмосферы
      marrPhaseVect[QUantCurNZSKVarsVS -1] = OperationCur.mEnvironment.calcAirDensity(marrPhaseVect[1]);
      ///

      //вычисление вектора воздушной скорости в NZSK
      long double arrUaSvSK[3]= {0.};
      OperationCur.calcUaSvSK_Case_NZSK(arrUaSvSK);
      ///

      // вычисление векторов силы аэродинамич сопротивления и момента этой силы в СвСК

      OperationCur.mHelic.calcRezAirF_and_Moment_SvSK(&(OperationCur.marrPhaseVect[6]), 0.
                           , OperationCur.marrPhaseVect[QUantCurNZSKVarsVS -1], arrUaSvSK, arrSvSK_AirForce, arrSvSK_AirMoment);
      ///
      long double alfZ = mHelic.marrRotor[0].mZaklinAng;
      *pvalFRotX = -cos(alfZ)* arrSvSK_AirForce[0] + sin(alfZ)* arrSvSK_AirForce[1]
              + sin(VAlNu -alfZ )* OperationCur.mHelic.mHelicMass * G_ZEMLI;

      *pvalFRotY = -sin(alfZ)* arrSvSK_AirForce[0] - cos(alfZ)* arrSvSK_AirForce[1]
              + cos(VAlNu -alfZ )* OperationCur.mHelic.mHelicMass * G_ZEMLI;

      long double valXcm = mHelic.marrRotor[0].mForceArmX;
      long double valYcm = (mHelic.marrRotor[0].mForceArmY + mHelic.marrRotor[1].mForceArmY) /2.;

      return -arrSvSK_AirMoment[2] - (*pvalFRotY) * valXcm  + (*pvalFRotX) * valYcm ;


  }
*/
 /*
  // вычисление разностной производной функции вертикальногшо балланса
  long double TPartHelicTraj::fncDerivF_HorizBalance(const long double VAlNu, const long  double VAl_dx)
  {
      long double valFRotX, valFRotY;
      long double arrASirForce[3] = {0.}, arrASirMom[3] = {0.};
      return (fncFGr_HorizBalance(VAlNu + VAl_dx, &valFRotX, &valFRotY, arrASirForce, arrASirMom)
              - fncFGr_HorizBalance(VAlNu, &valFRotX, &valFRotY,  arrASirForce, arrASirMom))/VAl_dx ;
  }
*/
 /*
  // нахождение параметров балансировки при горизонтальном полете  с заданной скоростью
  // INPUT:
  //VAlVx - горизонт.  скорость
  // VAlY- высота
  //OUTPUT:
  // arrCntrRulsPos [4] - вектор положения управляющих рулей (общ шаг, каппа, этта, педали)
  // arrCntrImpacts[4] - вектор управляющих воздействий (продольная сила, сила тяги, бовокая сила, крутящий момент
  //*pvalNu - угол тангажа
  bool TPartHelicTraj::calcBalanceParams_ForSteadyState_Move(const long double VAlVx, const long double  VAlY
                      , const long double VAlNu0 ,long double *arrCntrRulsPos
                      ,long double *arrCntrImpacts,long double *pvalNu)
  {
      memset(arrCntrRulsPos, 0, 4 * sizeof(long double));
      memset(arrCntrImpacts, 0, 4 * sizeof(long double));
     TPartHelicTraj PartHelicTrajWork = *this;
     memset(PartHelicTrajWork.marrPhaseVect, 0, sizeof(long double)* QUantCurNZSKVarsVS);
     PartHelicTrajWork.marrPhaseVect[1] = VAlY;
     PartHelicTrajWork.marrPhaseVect[3] = VAlVx;
     PartHelicTrajWork.marrPhaseVect[12] = PartHelicTrajWork.mEnvironment.calcAirDensity(VAlY);
     PartHelicTrajWork.marrPhaseVect[10] = VAlNu0;
     PartHelicTrajWork.doAirDensity();

     // решение уравненитя вертикальной балансировки по углу ню
     bool breturn = false;
     long double valFRotX, valFRotY;
     long double arrAirForce[3] = {0.}, arrAirMom[3] = {0.};
     int i = 0;
     for ( i = 0; i < 50;++i)
     {
         if (i ==4)
         {
             int iii=0;
         }
         long double valDel = - PartHelicTrajWork.fncFGr_HorizBalance(PartHelicTrajWork.marrPhaseVect[10], &valFRotX, &valFRotY,arrAirForce , arrAirMom)
                 / PartHelicTrajWork.fncDerivF_HorizBalance(PartHelicTrajWork.marrPhaseVect[10], 0.0001);
         PartHelicTrajWork.marrPhaseVect[10] += valDel;
         if (fabs(valDel) < 0.000001)
         {
              breturn = true;
              break;
         }
     }
     if (!breturn)
     {
      return false;
     }
     *pvalNu = PartHelicTrajWork.marrPhaseVect[10];
     arrCntrImpacts[0] = valFRotX;
     arrCntrImpacts[1] = valFRotY;
     PartHelicTrajWork.calcCntrRulsPos_New(arrCntrImpacts,&(arrCntrRulsPos[0])
                       ,&(arrCntrRulsPos[1]) ,&(arrCntrRulsPos[2]),&(arrCntrRulsPos[3]));
  // контроль
     long double valFi = arrCntrRulsPos[0];
     long double valKappa = arrCntrRulsPos[1];
     long double valEtta = arrCntrRulsPos[2];
     long double valDelFi = arrCntrRulsPos[3];
     memset(PartHelicTrajWork.marrPhaseVect, 0, sizeof(long double) * 13);
     PartHelicTrajWork.marrPhaseVect[1] = VAlY;
     PartHelicTrajWork.marrPhaseVect[3] = VAlVx;
     PartHelicTrajWork.marrPhaseVect[10] = *pvalNu;

     PartHelicTrajWork.marrPhaseVect[12] = PartHelicTrajWork.mEnvironment.calcAirDensity(VAlY);
     long double arrUaSvSK[3] ={0.};
     PartHelicTrajWork.calcUaSvSK_Case_NZSK(arrUaSvSK);

     long double arrSvSK_Force[3] = {0.}, arrSvSK_Moment[3] = {0.}
              , arrSvSK_ShaftForce0[3] = {0.},arrSvSK_ShaftMoment0[3] = {0.}
                ,arrSvSK_AirForce0[3] = {0.},arrSvSK_AirMoment0[3] = {0.};
     PartHelicTrajWork.mHelic.calcRezF_and_Moment_SvSK(&(PartHelicTrajWork.marrPhaseVect[6]), valFi, valKappa, valEtta, valDelFi
             ,PartHelicTrajWork.marrPhaseVect[12], arrUaSvSK, arrSvSK_Force, arrSvSK_Moment, arrSvSK_ShaftForce0,arrSvSK_ShaftMoment0
             ,arrSvSK_AirForce0,arrSvSK_AirMoment0);

     ///
     return true;
  }

*/
/*
  // вычисление вектора частных производных по переменной с номером j
  // вектор функций Fa и Ma. разностным методом
  // INPUT:
  //Val_dxj - разность по аргументу
  // j - номер переменной по которой вычисляется производная
  // arrFa[3] - вектор аэродинамической силы
  //arrMa[3] - вектор аэродинамического момента
  //
  // OUTPUT:
  //arr_dFa_po_dxj[3] - производная вектора силы по переменной с номером j
  //arr_dFa_po_dxj[3] - производная вектора момента по переменной с номером j
  void TPartHelicTraj::calc_dFa_and_dMa_po_dxj(const long double  Val_dxj, const int j,long double  *arrFa
                  ,long double  *arrMa,long double  * arr_dFa_po_dxj,long double  *arr_dMa_po_dxj)
  {
          TPartHelicTraj PartHelicTrajTemp = *this;
          PartHelicTrajTemp.marrPhaseVect[j] += Val_dxj;


         long double  arrUaSvSK[3] = {0.};
         // PartHelicTrajTemp.calcUaSvSK_Case_SvSK(arrUaSvSK);
          PartHelicTrajTemp.calcUaSvSK_Case_NZSK(arrUaSvSK);
          // вычисление плотности атмосферы
         PartHelicTrajTemp.marrPhaseVect[12] = PartHelicTrajTemp.mEnvironment.calcAirDensity(PartHelicTrajTemp.marrPhaseVect[1]);
          ///

         long double   arrSvSK_AirForce[3] = {0.}, arrSvSK_AirMoment[3] = {0.}, arrT0[3] = {0.};
          PartHelicTrajTemp.mHelic.calcRezAirF_and_Moment_SvSK(&(PartHelicTrajTemp.marrPhaseVect[6]), 0.
                            , PartHelicTrajTemp.marrPhaseVect[12], arrUaSvSK, arrSvSK_AirForce, arrSvSK_AirMoment);
          MtrxMinusMatrx(arrSvSK_AirForce, arrFa,3, 1, arrT0);
          MatrxDivideScalar(arrT0, 3, 1, Val_dxj, arr_dFa_po_dxj);

          MtrxMinusMatrx(arrSvSK_AirMoment, arrMa,3, 1, arrT0);
          MatrxDivideScalar(arrT0, 3, 1, Val_dxj, arr_dMa_po_dxj);

  }

*/
//-----------------------------------------
// проверка того, что часть траектории закончилась
  bool TPartHelicTraj::IsEndOfPart()
  {

  }

  //---------------------------------------------------------------------------------
   double  TPartHelicTraj::fncMod2Pi__( double a_fVal )
   {
       double  fAngle = fmod( a_fVal, 2. * M_PI );
       if ( fAngle<0 )
           fAngle += 2. * M_PI;
       return  fAngle;
   }

   //---------------------------------------------------------------------------------
   double  TPartHelicTraj::fnc_Minus_PI_Plus_PI( double a_fVal )
   {

       double  fAngle = fmod( a_fVal, 2. * M_PI );
       if ( fabs(fAngle) > M_PI )
       {
           fAngle -= 2. * M_PI * SIGNUM__(fAngle);
       }

       return  fAngle;
   }

   //---------------------------------------------------------------------------------






