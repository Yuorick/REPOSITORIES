#include <math.h>
#include <string.h>
#include <float.h>
#include "Rotor.h"
#include "MatrixProccess.h"
#include "Environment.h"


  TRotor::TRotor()
{
// лопасть винта
	 mBlade = TBlade();
	// к-во лопастей
	 mQuantBlades = 0;
	 //
     mBasePLane  = TLongPlane ();
     //
     mOmega = 0.;
     // коэффициент тяги винта CT
     mCt = 0.;;
    // коэффициент компенсатора взмаха, константа,   0.6 > k > 0.4
     mKWave = 0.;
     mZaklinAng = 0.;

   //  mForceArmX = 0.;
    // mForceArmY = 0.;
     mbFullTypeoffEq = false;
}
// Конструктор копирования
  TRotor::TRotor (const TRotor &R)
 {
	 mBlade = R.mBlade;
	 mQuantBlades   = R.mQuantBlades ;
	 mBasePLane  = R.mBasePLane ;

     mOmega = R.mOmega;
     mCt = R.mCt;
     mKWave = R.mKWave;
     mZaklinAng = R.mZaklinAng;
   //  mForceArmX = R.mForceArmX;
   //  mForceArmY = R.mForceArmY;
     mbFullTypeoffEq = R.mbFullTypeoffEq;

 }
 // оператор присваивания
  TRotor TRotor::operator=(TRotor  R)
 {
	 mBlade = R.mBlade;
	 mQuantBlades   = R.mQuantBlades ;
	 mBasePLane  = R.mBasePLane ;

     mOmega = R.mOmega;
     mCt = R.mCt;
     mKWave = R.mKWave;
     mZaklinAng = R.mZaklinAng;
    // mForceArmX = R.mForceArmX;
    // mForceArmY = R.mForceArmY;
     mbFullTypeoffEq = R.mbFullTypeoffEq;

	 return *this ;
 }

/*
 // парам констр
 // Blade - допасть
 // QuantBlades - к-во лопастей
 // VAlZaklinAng - угол заклинения винта
 // arrCoord0[3] - координаты точки осноования оси вращения винта в СвСК
 // L - расстояние олт точки лснования винта до центра вращения
 // FiMax - максимально допустимый общий шаг НВ
 TRotor::TRotor(const TBlade Blade, const int QuantBlades, const long double VAlZaklinAng
                , const long double VAlCoordX, const long double VAlCoordY,const long  double L, const long double Omega)
 {
     mBlade = Blade;
     mQuantBlades = QuantBlades;

     long double arrS0[3] = {0.};// координаты центра вращения в СвСК
     long double arrF[9] = {0.}; // матрица ортов БСК в осях СвСК
     arrS0[0] = VAlCoordX + L * sinl(VAlZaklinAng);
     arrS0[1] = VAlCoordY + L * cosl(VAlZaklinAng);
     arrF[0] =  cosl(VAlZaklinAng);
     arrF[1] =  sinl(VAlZaklinAng);
     arrF[3] = -sinl(VAlZaklinAng);
     arrF[4] =  cos(VAlZaklinAng);
     arrF[8] =  1.;
     mBasePLane = TLongPlane(  arrS0,  arrF);

     mOmega = Omega;
     // коэффициент тяги винта CT
     mCt = 0.;
    // коэффициент компенсатора взмаха, константа,   0.6 > k > 0.4
     mKWave = 0.;
     mZaklinAng = VAlZaklinAng;
     calcForceArms(&mForceArmX, &mForceArmY) ;
     mbFullTypeoffEq = false;

 }
 */
 //-------------------------------------------------------


  // парам констр
  // Blade - допасть
  // QuantBlades - к-во лопастей
  // VAlZaklinAng - угол заклинения винта
  // arrCoord0[3] - координаты точки осноования оси вращения винта в СвСК
  // L - расстояние олт точки лснования винта до центра вращения
  // FiMax - максимально допустимый общий шаг НВ
  TRotor::TRotor(const TBlade Blade, const int QuantBlades, const long double VAlZaklinAng
                 , const long double VAlCoordX, const long double VAlCoordY,const long  double L
                 , const long double Omega, const bool bFullType)
  {
      mBlade = Blade;
      mQuantBlades = QuantBlades;

      long double arrS0[3] = {0.};// координаты центра вращения в СвСК
      long double arrF[9] = {0.}; // матрица ортов БСК в осях СвСК
      arrS0[0] = VAlCoordX + L * sinl(VAlZaklinAng);
      arrS0[1] = VAlCoordY + L * cosl(VAlZaklinAng);
      arrF[0] =  cosl(VAlZaklinAng);
      arrF[1] =  sinl(VAlZaklinAng);
      arrF[3] = -sinl(VAlZaklinAng);
      arrF[4] =  cos(VAlZaklinAng);
      arrF[8] =  1.;
      mBasePLane = TLongPlane(  arrS0,  arrF);

      mOmega = Omega;
      // коэффициент тяги винта CT
      mCt = 0.;
     // коэффициент компенсатора взмаха, константа,   0.6 > k > 0.4
      mKWave = 0.;
      mZaklinAng = VAlZaklinAng;
    //  calcForceArms(&mForceArmX, &mForceArmY) ;
      mbFullTypeoffEq = bFullType;

  }
  //-------------------------------------------------------


  // парам констр
  // Blade - допасть
  // QuantBlades - к-во лопастей
  // VAlZaklinAng - угол заклинения винта

  // L - расстояние от центра масс вертолета  до центра вращения винта

  TRotor::TRotor(const TBlade Blade, const int QuantBlades, const long double VAlZaklinAng
                 ,const long  double L, const long double Omega, const bool bFullType)
  {
      mBlade = Blade;
      mQuantBlades = QuantBlades;

      long double arrS0[3] = {0.};// координаты центра вращения в СвСК
      long double arrF[9] = {0.}; // матрица ортов БСК в осях СвСК
      arrS0[0] =  L * sinl(VAlZaklinAng);
      arrS0[1] =  L * cosl(VAlZaklinAng);
      arrF[0] =  cosl(VAlZaklinAng);
      arrF[1] =  sinl(VAlZaklinAng);
      arrF[3] = -sinl(VAlZaklinAng);
      arrF[4] =  cos(VAlZaklinAng);
      arrF[8] =  1.;
      mBasePLane = TLongPlane(  arrS0,  arrF);

      mOmega = Omega;
      // коэффициент тяги винта CT
      mCt = 0.;
     // коэффициент компенсатора взмаха, константа,   0.6 > k > 0.4
      mKWave = 0.;
      mZaklinAng = VAlZaklinAng;
    //  calcForceArms(&mForceArmX, &mForceArmY) ;
      mbFullTypeoffEq = bFullType;

  }
  //-------------------------------------------------------


 // Вычисление результирующей силы и момента создаваемых винтом в СвСК
 //INPUT:
 // VAlFi, VAl_DnKappa, VAl_DnEtta - общий шаг, угол отклонения автомата перекроса по тангажу * Dn
                                  // , угол отклонения автомата перекроса по крену * Dn
 // arrUa[3] - вектор воздушной скорости вертолета в СвСК
 // VAlRo - плотность воздуха
 // OUTPUT:
 // arrRezF[3] - результирующий вектор сил
 // arrRezMom[3] - результирующий вектор момента сил
void TRotor::calcRezF_and_RezMom(const long double VAlFi, const  long double VAl_DnKappa
      ,const long  double VAl_DnEtta ,long double *arrSvSKUa, const long  double VAlRo
       ,long  double *arrOMegaSvSK, long double *arrRezSvSK_F, long  double *arrRezSvSK_Mom)
{
    if(mbFullTypeoffEq)
    {
        calcRezF_and_RezMom_FullType(VAlFi,  VAl_DnKappa
              , VAl_DnEtta ,arrSvSKUa, VAlRo
               ,arrOMegaSvSK, arrRezSvSK_F, arrRezSvSK_Mom);
    }
    else
    {
        calcRezF_and_RezMom_SimpleType(VAlFi,  VAl_DnKappa
              , VAl_DnEtta ,arrSvSKUa, VAlRo
               ,arrOMegaSvSK, arrRezSvSK_F, arrRezSvSK_Mom);
    }


}

// вычисление угла атаки конструктивной плоскости вращения винта
// arrSvSKUa[3] - воздушная скорость веротолета в СвСК
long double TRotor::calcAttackAng(long double *arrSvSKUa)
{
     long double valUa = Norm3( arrSvSKUa) ;
     if(valUa < 0.000001)
     {
         return 0.;
     }
    // пересчет вектора воздушной скорости в БСК
  long double arrUaBSK[3] = {0.};
  MtrxTranspMultMatrx(mBasePLane.marrF,3, 3, arrSvSKUa, 1, arrUaBSK) ;
  long double temp = arrUaBSK[1]/ valUa;
  if(fabs(temp) >= 1.)
  {
      temp = (arrUaBSK[1] > 0.)?0.9999999999:-0.9999999999;
  }
  return asinl (temp ) ;
}


long double TRotor::calcQreact(const long double VAlRo, const long double Va)
{

    double temp0 = (mBlade.mBladeR * mBlade.mBladeR * mBlade.mBladeR * mBlade.mBladeR  - mBlade.mRadHorizHsarnir* mBlade.mRadHorizHsarnir* mBlade.mRadHorizHsarnir* mBlade.mRadHorizHsarnir);
    double temp1 = (mBlade.mBladeR * mBlade.mBladeR  -  mBlade.mRadHorizHsarnir* mBlade.mRadHorizHsarnir);
    return 3./8. * mBlade.mCX0* mBlade.mBlade_b * VAlRo
                   *(temp0*mOmega*mOmega + temp1  *Va *Va) ;
}






//------------------------------------------------

// вычисление плеча сил винта в БСК
//
void  TRotor::calcForceArms(long double *pForceArmX,long  double *pForceArmY)
{
    long double arrS0Wave[3] = {0.};
    MtrxTranspMultMatrx(mBasePLane.marrF,3, 3, mBasePLane.marrS0, 1, arrS0Wave) ;
    *pForceArmX = arrS0Wave[0];
    *pForceArmY = arrS0Wave[1];
}

//--------------------------------------------------
// вычисление константы при горизонтальных шарнирах моменета
long double TRotor::calcHorSharConst()
{
    return ((double) mQuantBlades)*mBlade.mSg * mBlade.mRadHorizHsarnir * mOmega * mOmega/ 2.;
}

//------------------------------------------------
void TRotor::calCoef_a1_and_b1(const long double VAlFi, const long double VAlRo
              ,long double * arrSvSKUa, long double *pval_a1, long double *pval_b1)
{

    // 1. вычисление силы тяги
      long  double valSquare =  M_PI * mBlade.mBladeR* mBlade.mBladeR;
      long  double valTiaga = 1./2. * mCt * VAlRo * valSquare
               * (mOmega * mBlade.mBladeR )* (mOmega * mBlade.mBladeR )* VAlFi;
       ///
//2. индуктивная скорость всасывания

  long  double val_v = sqrt(valTiaga/ 2./VAlRo / valSquare);

   //3. угол атаки конструктивной плоскости вращения
   const long double valAttackAng =   calcAttackAng(arrSvSKUa);
   ///

   // 4. вычисление безразмерных коэффициентов лямбда и мю
   const long double valUa = Norm3( arrSvSKUa) ;
   long double valOmegaR =  mOmega * mBlade.mBladeR;
   long double valLamb = (valUa * sinl(valAttackAng) - val_v)/ valOmegaR;
   long double valMu = valUa * cosl(valAttackAng) / valOmegaR;
   ///

   // 5. вычисление массвой характ лопасти
   long double valGamma_lamb = mBlade.mCyalfa * VAlRo * mBlade.mBlade_b
           *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR
           /mBlade.mIg / 2.;
   // 6. вычисление коэффициентов махового движения без учета компенсатора взмаха
   long double a0 = valGamma_lamb * (valLamb / 3. + VAlFi/ 4.)- mBlade.mSg * G_ZEMLI/mBlade.mIg/ (mOmega*mOmega);

   long double a10 = 2. * valMu * (valLamb + 4./3.*VAlFi);

   long double b10 = 4./3.* valMu *a0;

   // 7. вычисление коэффициентов махового движения с учетом компенсатора взмаха
   *pval_a1 = (a10 + mKWave * b10)/ (1. + mKWave * mKWave);
   *pval_b1 = (b10 - mKWave * a10)/ (1. + mKWave * mKWave)* Sign(mOmega);
}

long double TRotor::Sign(long double a)
{
   return (fabsl(a) < DBL_MIN)?-1.:0.;
}

//---------------------------------------------------------------------------
// вычисление производной силы тяги винта по общему шагу
// VAlRo - плотность воздуха
long double TRotor::calc_Deriv_T_po_Fi(const long double VAlRo)
{
     // радиус ометаемой площади несущего винта
     long double valR = mBlade.mBladeR;
     ///

     long double valSquare =  M_PI * valR * valR;
     return  mCt * VAlRo * valSquare
                 * (mOmega * valR )* (mOmega * valR )/ 2.;
}

//--------------------------------------------------------------
void TRotor::calcRezF_and_RezMom_FullType(const long double VAlFi, const  long double VAl_DnKappa
      ,const long  double VAl_DnEtta ,long double *arrSvSKUa, const long  double VAlRo
       ,long  double *arrOMegaSvSK, long double *arrRezSvSK_F, long  double *arrRezSvSK_Mom)
{
    // 1. вычисление силы тяги
      long  double valSquare =  M_PI * mBlade.mBladeR* mBlade.mBladeR;
      long  double valTiaga = 1./2. * mCt * VAlRo * valSquare
               * (mOmega * mBlade.mBladeR )* (mOmega * mBlade.mBladeR )* VAlFi;
       ///

    //2. индуктивная скорость всасывания

      long  double val_v = sqrt(valTiaga/ 2./VAlRo / valSquare);

       //3. угол атаки конструктивной плоскости вращения
       const long double valAttackAng =   calcAttackAng(arrSvSKUa);
       ///

       // 4. вычисление безразмерных коэффициентов лямбда и мю
       const long double valUa = Norm3( arrSvSKUa) ;
       long double valOmegaR =  mOmega * mBlade.mBladeR;
       long double valLamb = (valUa * sinl(valAttackAng) - val_v)/ valOmegaR;
       long double valMu = valUa * cosl(valAttackAng) / valOmegaR;
       ///

       // 5. вычисление массвой характ лопасти
       long double valGamma_lamb = mBlade.mCyalfa * VAlRo * mBlade.mBlade_b
               *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR
               /mBlade.mIg / 2.;
       // 6. вычисление коэффициентов махового движения без учета компенсатора взмаха
       long double a0 = valGamma_lamb * (valLamb / 3. + VAlFi/ 4.)- mBlade.mSg * G_ZEMLI/mBlade.mIg/ (mOmega*mOmega);

       long double a10 = 2. * valMu * (valLamb + 4./3.*VAlFi);

       long double b10 = 4./3.* valMu *a0;

       //long double a10 = 0., b10 = 0.;
       calCoef_a1_and_b1( VAlFi,  VAlRo
                     ,arrSvSKUa, &a10, &b10) ;


       // 7. вычисление коэффициентов махового движения с учетом компенсатора взмаха
       long double a1 = (a10 + mKWave * b10)/ (1. + mKWave * mKWave);
       long double b1 = (b10 - mKWave * a10)/ (1. + mKWave * mKWave);
       ///

       //8. вычисление вращательных производных
           //  вычисление массвой характ лопасти
           // long double valGamma_lamb = mBlade.mCyalfa * VAlRo * mBlade.mBlade_b
            //   *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR
          //    /mBlade.mIg / 2.;
           ///
       long double a1Omegaz = -(mKWave *valGamma_lamb + 8.)/valGamma_lamb/mOmega/(1. + mKWave * mKWave);
       long double a1Omegax =  (8. * mKWave  - valGamma_lamb)/valGamma_lamb/mOmega/(1. + mKWave * mKWave);
       long double b1Omegax = -a1Omegaz;
       long double b1Omegaz = -a1Omegax;
       ///

       // 9. пересчет вектора мгновенной угловой скорости вертолета из СвСК в БСК
       long double arrOMegaBSK[3] = {0.};
       MtrxTranspMultMatrx(mBasePLane.marrF,3, 3, arrOMegaSvSK, 1, arrOMegaBSK) ;
       ///

       // 10. расчет добавок в углы смещения конуса вращения, вызванных кориолисовыми ускорекниями
       double delta_a1 =  a1Omegax *  arrOMegaBSK[0] + a1Omegaz*  arrOMegaBSK[2];
       double delta_b1 = -a1Omegaz *  arrOMegaBSK[0] - a1Omegax*  arrOMegaBSK[2];
       ///

       // 11. вычисление углов нормали конуса вращения лопастей
      long  double aRez = a1 + delta_a1 + VAl_DnKappa; // НЕ ЗАБЫТЬ!!!!!!!
      long  double bRez = b1 + delta_b1 + VAl_DnEtta;
      // long double aRez =  VAl_DnKappa; //
      // long double bRez =  VAl_DnEtta;
       ///

       // 12. вычисление результирующей силы в БСК
      long  double arrRezBSK_F[3] = {0.};
       arrRezBSK_F[0] = valTiaga * aRez;
       arrRezBSK_F[1] = valTiaga;
       arrRezBSK_F[2] = valTiaga * bRez;
       ///

       /////////  ВЫЧИСЛЕНИЕ МОМЕНТОВ //////////////////////////////////////////////////////////////////
    // ТАК БЫЛО
       /*
       // 13. ПЕРЕСЧЕТ ВЕКТОРА начала координат БСК из СвСК в сиситему координат с осями параллельными БСК и
       // началом совпадающим с центорм СвСК (центром масс вертолета)
       long double arrS0Wave[3] = {0.};
       MtrxTranspMultMatrx(mBasePLane.marrF,3, 3, mBasePLane.marrS0, 1, arrS0Wave) ;

       ///




       // 14. вычисление момента сил силы arrRezBSK_F
       long double arrRezF_Mom[3] = {0.};
        OuterProduct(arrS0Wave , arrRezBSK_F, arrRezF_Mom) ;
        */

       // по другому
        // 13. ПЕРЕСЧЕТ ВЕКТОРА силы в СвСК

        // началом совпадающим с центорм СвСК (центром масс вертолета)
        long double arrF_SvSK[3] = {0.};
        MtrxMultMatrx(mBasePLane.marrF,3, 3, arrRezBSK_F, 1, arrF_SvSK) ;

        ///




        // 14. вычисление момента сил силы arrRezBSK_F
        long double arrRezF_Mom_SvSK[3] = {0.};
         OuterProduct(mBasePLane.marrS0 , arrF_SvSK, arrRezF_Mom_SvSK) ;
        ///

        // 15.вычисление	управляющих моментов, создаваемых центробежной силой вращения лопастей из-за разносов горизонтальных шарниров.
       long  double arrHorShar_Mom[3] = {0.};
        arrHorShar_Mom[0] = bRez;//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
        arrHorShar_Mom[2] = aRez;//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
        long double valScal = calcHorSharConst();

        MatrxMultScalar(arrHorShar_Mom, 1, 3,  valScal, arrHorShar_Mom);//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
        ///
        // так было
        /*
        // 16. вычисление суммарного управляющего момента
        long double arrSumControl_Mom[3] = {0.};
        MtrxSumMatrx(arrHorShar_Mom, arrRezF_Mom,1, 3, arrSumControl_Mom) ;
        */
        ///
        // так стало
        long double arrSumControl_Mom[3] = {0.};
        memcpy(arrSumControl_Mom,arrHorShar_Mom, 3 * sizeof(long double) );
        // 17.вычисление демпфирующего момента
      //  long double valYae = arrS0Wave[1]*(1.
        //      + ((long double)mQuantBlades) *  mBlade.mRadHorizHsarnir * mBlade.mSg * mOmega * mOmega/ 2. / arrS0Wave[1]/valTiaga );
        long  double arrDemp_Mom[3] = {0.};
        arrDemp_Mom[0] = arrSumControl_Mom[0] * b1Omegax * arrOMegaBSK[0];//!!!!!!!!!
        arrDemp_Mom[2] = arrSumControl_Mom[0] * a1Omegaz * arrOMegaBSK[2];//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
        ///

        // 18.вычисление инерционных моментов
        long  double arrK_Mom[3] = {0.};
         arrK_Mom[0] = -((double) mQuantBlades)*mBlade.mSg * mBlade.mRadHorizHsarnir * mOmega * arrOMegaBSK[2];
         arrK_Mom[2] = -((double) mQuantBlades)*mBlade.mSg * mBlade.mRadHorizHsarnir * mOmega * arrOMegaBSK[0];
         ///

         // 19. вычисление реактивного момента
       long  double arrY_Mom[3] = {0.};
           // 19.1 вычисление модуля проекции воздушной скорости на плоскость OXZ БСК
        // const long double valUa = Norm3( arrSvSKUa) ;

        // const long double valAttackAng =   calcAttackAng(arrSvSKUa);// угол атаки конструктивной плоскости вращения
         ///
         long   double Va = valUa * cos (valAttackAng);
          arrY_Mom[1] = 3./8. * mBlade.mCX0* mBlade.mBlade_b * VAlRo * VAlFi
                  *(
                    (mBlade.mPofile_d1* mBlade.mPofile_d1* mBlade.mPofile_d1* mBlade.mPofile_d1 - mBlade.mPofile_d0* mBlade.mPofile_d0* mBlade.mPofile_d0* mBlade.mPofile_d0)
                       *mOmega*mOmega
                    + ( mBlade.mPofile_d1* mBlade.mPofile_d1 -  mBlade.mPofile_d0* mBlade.mPofile_d0) *Va *Va
                     ) * ( mOmega/ fabs(mOmega));
           arrY_Mom[1] = calcQreact( VAlRo,  Va) * ( mOmega/ fabsl(mOmega))*  VAlFi;


         // 19. суммирование суммарного управляющего момента, демпфирующего и инерционного
        long double arrRezBSK_Mom[3] = {0.}, arrTemp0[3] = {0.}, arrTemp1[3] = {0.};
         MtrxSumMatrx(arrSumControl_Mom, arrDemp_Mom,1, 3, arrTemp0) ;
         MtrxSumMatrx(arrTemp0, arrK_Mom,1, 3, arrTemp1) ;
         MtrxSumMatrx(arrTemp1, arrY_Mom,1, 3,arrRezBSK_Mom) ;

         // персчет векторов силы и момента в СвСК

         MtrxMultMatrx(mBasePLane.marrF,3, 3, arrRezBSK_F, 1, arrRezSvSK_F) ;
         MtrxMultMatrx(mBasePLane.marrF,3, 3, arrRezBSK_Mom, 1, arrRezSvSK_Mom) ;//
         // так стало
         MtrxSumMatrx(arrRezSvSK_Mom, arrRezF_Mom_SvSK,1, 3,arrTemp1) ;
         memcpy(arrRezSvSK_Mom, arrTemp1, 3 * sizeof(long double)) ;
}

//-----------------------------------------------------------------------

// Вычисление результирующей силы и момента создаваемых винтом в СвСК
// для упрощенной модели винта
//INPUT:
// VAlFi, VAl_DnKappa, VAl_DnEtta - общий шаг, угол отклонения автомата перекроса по тангажу * Dn
                                 // , угол отклонения автомата перекроса по крену * Dn
// arrUa[3] - вектор воздушной скорости вертолета в СвСК
// VAlRo - плотность воздуха
// OUTPUT:
// arrRezF[3] - результирующий вектор сил
// arrRezMom[3] - результирующий вектор момента сил
void TRotor::calcRezF_and_RezMom_SimpleType(const long double VAlFi, const  long double VAl_DnKappa
     ,const long  double VAl_DnEtta ,long double *arrSvSKUa, const long  double VAlRo
      ,long  double *arrOMegaSvSK, long double *arrRezSvSK_F, long  double *arrRezSvSK_Mom)
{

// 1. вычисление силы тяги
  long  double valSquare =  M_PI * mBlade.mBladeR* mBlade.mBladeR;
  long  double valTiaga = 1./2. * mCt * VAlRo * valSquare
           * (mOmega * mBlade.mBladeR )* (mOmega * mBlade.mBladeR )* VAlFi;
   ///
 /*
//2. индуктивная скорость всасывания

  long  double val_v = sqrt(valTiaga/ 2./VAlRo / valSquare);

   //3. угол атаки конструктивной плоскости вращения
   const long double valAttackAng =   calcAttackAng(arrSvSKUa);
   ///

   // 4. вычисление безразмерных коэффициентов лямбда и мю
   const long double valUa = Norm3( arrSvSKUa) ;
   long double valOmegaR =  mOmega * mBlade.mBladeR;
   long double valLamb = (valUa * sinl(valAttackAng) - val_v)/ valOmegaR;
   long double valMu = valUa * cosl(valAttackAng) / valOmegaR;
   ///

   // 5. вычисление массвой характ лопасти
   long double valGamma_lamb = mBlade.mCyalfa * VAlRo * mBlade.mBlade_b
           *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR
           /mBlade.mIg / 2.;
   // 6. вычисление коэффициентов махового движения без учета компенсатора взмаха
   long double a0 = valGamma_lamb * (valLamb / 3. + VAlFi/ 4.)- mBlade.mSg * G_ZEMLI/mBlade.mIg/ (mOmega*mOmega);

   long double a10 = 2. * valMu * (valLamb + 4./3.*VAlFi);

   long double b10 = 4./3.* valMu *a0;
   */
   long double a10 = 0., b10 = 0.;
   calCoef_a1_and_b1( VAlFi,  VAlRo
                 ,arrSvSKUa, &a10, &b10) ;


   // 7. вычисление коэффициентов махового движения с учетом компенсатора взмаха
   long double a1 = (a10 + mKWave * b10)/ (1. + mKWave * mKWave);
   long double b1 = (b10 - mKWave * a10)/ (1. + mKWave * mKWave);
   ///

   //8. вычисление вращательных производных
       //  вычисление массвой характ лопасти
        long double valGamma_lamb = mBlade.mCyalfa * VAlRo * mBlade.mBlade_b
           *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR *mBlade.mBladeR
           /mBlade.mIg / 2.;
       ///
   long double a1Omegaz = -(mKWave *valGamma_lamb + 8.)/valGamma_lamb/mOmega/(1. + mKWave * mKWave);
   long double a1Omegax =  (8. * mKWave  - valGamma_lamb)/valGamma_lamb/mOmega/(1. + mKWave * mKWave);
   long double b1Omegax = -a1Omegaz;
   long double b1Omegaz = -a1Omegax;
   ///

   // 9. пересчет вектора мгновенной угловой скорости вертолета из СвСК в БСК
   long double arrOMegaBSK[3] = {0.};
   MtrxTranspMultMatrx(mBasePLane.marrF,3, 3, arrOMegaSvSK, 1, arrOMegaBSK) ;
   ///

   // 10. расчет добавок в углы смещения конуса вращения, вызванных кориолисовыми ускорекниями
   double delta_a1 =  a1Omegax *  arrOMegaBSK[0] + a1Omegaz*  arrOMegaBSK[2];
   double delta_b1 = -a1Omegaz *  arrOMegaBSK[0] - a1Omegax*  arrOMegaBSK[2];
   ///

   // 11. вычисление углов нормали конуса вращения лопастей
  //long  double aRez = a1 + delta_a1 + VAl_DnKappa; // НЕ ЗАБЫТЬ!!!!!!!
  //long  double bRez = b1 + delta_b1 + VAl_DnEtta;
   long double aRez =  VAl_DnKappa; //
   long double bRez =  VAl_DnEtta;
   ///

   // 12. вычисление результирующей силы в БСК
  long  double arrRezBSK_F[3] = {0.};
   arrRezBSK_F[0] = valTiaga * aRez;
   arrRezBSK_F[1] = valTiaga;
   arrRezBSK_F[2] = valTiaga * bRez;
   ///

   /////////  ВЫЧИСЛЕНИЕ МОМЕНТОВ //////////////////////////////////////////////////////////////////

   // 13. ПЕРЕСЧЕТ ВЕКТОРА начала координат БСК из СвСК в сиситему координат с осями параллельными БСК и
   // началом совпадающим с центорм СвСК (центром масс вертолета)
   long double arrS0Wave[3] = {0.};
   MtrxTranspMultMatrx(mBasePLane.marrF,3, 3, mBasePLane.marrS0, 1, arrS0Wave) ;
   //MatrxMultScalar(arrS0Wave, 1, 3, -1.,arrS0Wave);
   ///




   // 14. вычисление момента сил силы arrRezBSK_F
   long double arrRezF_Mom[3] = {0.};
    OuterProduct(arrS0Wave , arrRezBSK_F, arrRezF_Mom) ;
    ///

    // 15.вычисление	управляющих моментов, создаваемых центробежной силой вращения лопастей из-за разносов горизонтальных шарниров.
   long  double arrHorShar_Mom[3] = {0.};
   // arrHorShar_Mom[0] = bRez;//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
   // arrHorShar_Mom[2] = aRez;//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
    long double valScal = calcHorSharConst();

    MatrxMultScalar(arrHorShar_Mom, 1, 3,  valScal, arrHorShar_Mom);//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
    ///

    // 16. вычисление суммарного управляющего момента
    long double arrSumControl_Mom[3] = {0.};
    MtrxSumMatrx(arrHorShar_Mom, arrRezF_Mom,1, 3, arrSumControl_Mom) ;
    ///

    // 17.вычисление демпфирующего момента
    long double valYae = arrS0Wave[1]*(1.
          + ((long double)mQuantBlades) *  mBlade.mRadHorizHsarnir * mBlade.mSg * mOmega * mOmega/ 2. / arrS0Wave[1]/valTiaga );
   long  double arrDemp_Mom[3] = {0.};
 //   arrDemp_Mom[0] = arrSumControl_Mom[0] * b1Omegax * arrOMegaBSK[0];//!!!!!!!!!
 //   arrDemp_Mom[2] = arrSumControl_Mom[0] * a1Omegaz * arrOMegaBSK[2];//!!!!!!!!! НЕ ЗАБЫТЬ!!!!!!!!!!
    ///

    // 18.вычисление инерционных моментов
    long  double arrK_Mom[3] = {0.};
   //  arrK_Mom[0] = -((double) mQuantBlades)*mBlade.mSg * mBlade.mRadHorizHsarnir * mOmega * arrOMegaBSK[2];
   //  arrK_Mom[2] = -((double) mQuantBlades)*mBlade.mSg * mBlade.mRadHorizHsarnir * mOmega * arrOMegaBSK[0];
     ///

     // 19. вычисление реактивного момента
   long  double arrY_Mom[3] = {0.};
       // 19.1 вычисление модуля проекции воздушной скорости на плоскость OXZ БСК
     const long double valUa = Norm3( arrSvSKUa) ;

     const long double valAttackAng =   calcAttackAng(arrSvSKUa);// угол атаки конструктивной плоскости вращения
     ///
     long   double Va = valUa * cos (valAttackAng);
       //arrY_Mom[1] = 3./8. * mBlade.mCX0* mBlade.mBlade_b * VAlRo * VAlFi
             // *(
             //   (mBlade.mPofile_d1* mBlade.mPofile_d1* mBlade.mPofile_d1* mBlade.mPofile_d1 - mBlade.mPofile_d0* mBlade.mPofile_d0* mBlade.mPofile_d0* mBlade.mPofile_d0)
              //     *mOmega*mOmega
              //  + ( mBlade.mPofile_d1* mBlade.mPofile_d1 -  mBlade.mPofile_d0* mBlade.mPofile_d0) *Va *Va
               //  ) * ( mOmega/ fabs(mOmega));
       arrY_Mom[1] = calcQreact( VAlRo,  Va) * ( mOmega/ fabsl(mOmega))*  VAlFi;


     // 19. суммирование суммарного управляющего момента, демпфирующего и инерционного
    long double arrRezBSK_Mom[3] = {0.}, arrTemp0[3] = {0.}, arrTemp1[3] = {0.};
     MtrxSumMatrx(arrSumControl_Mom, arrDemp_Mom,1, 3, arrTemp0) ;
     MtrxSumMatrx(arrTemp0, arrK_Mom,1, 3, arrTemp1) ;
     MtrxSumMatrx(arrTemp1, arrY_Mom,1, 3,arrRezBSK_Mom) ;

     // персчет векторов силы и момента в СвСК

     MtrxMultMatrx(mBasePLane.marrF,3, 3, arrRezBSK_F, 1, arrRezSvSK_F) ;
     MtrxMultMatrx(mBasePLane.marrF,3, 3, arrRezBSK_Mom, 1, arrRezSvSK_Mom) ;
}

