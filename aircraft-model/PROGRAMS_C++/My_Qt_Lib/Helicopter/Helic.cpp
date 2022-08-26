//---------------------------------------------------------------------------

#include <math.h>
#include <string.h>

#include "Helic.h"
#include "Blade.h"
#include "AbstractGliderElement.h"
#include "MatrixProccess.h"
#include "HelicConstants.h"

extern  long double constArrGearSets[];
extern  long double constArrHelicPlanersData[8 * 13];
extern  long double constArrHelicPlanersData_[8 * 13];
extern const int NUmRowsArrHelicPlanersData;
extern const int NUmColsArrHelicPlanersData;
extern const long double constBladeR ; // радиус ометаемой винтом площади
extern const long double constRadHorizHsarnir ; //расстояние от центра гориз шарнира до оси вращения вала винта marrDblSpinBoxBlade[0]
extern const long double constPofile_d0 ; // высота профиля у оси вала винта marrDblSpinBoxBlade[1]
extern const long double constPofile_d1 ; // высота профиля на конце marrDblSpinBoxBlade[2]
extern const long double constBlade_b ; // хорда     marrDblSpinBoxBlade[3]
extern const long double constBladeM ; // масса лопасти
extern const long double constBladeCX0 ;
extern const int constQuantBlades ;
extern const long double constFiMax ;
extern const long double constRotorOmega  ;
//extern const long double constShaftAxeX ;
extern const long double constShaftAxeX_ ;
extern const long double constShaftAxeY_;
//extern const long double constShaftAxeY;
//extern const long double constZaclinAng ;
extern const long double constZaclinAng_ ;
//extern const long double constShaftUpperL ;
//extern const long double constShaftLowerL ;
extern const long double constShaftUpperL_ ;
extern const long double constShaftLowerL_ ;
extern const long double constHelicMass;
//extern  long double ARrINertMtrx[9];
extern  long double ARrINertMtrx_[9];
extern  long double ARrINertMtrxPrived[9];
extern const  long double constAlfaZaklNew;
extern const  long double constXdistLow;
extern const  long double  constDeltaXUp_Low;





 THelic::THelic()
{

  /*   memset(marrRotor, 0, 2 * sizeof(TRotor));
	// частота вращениеия винта
    // memset(marrRotorOmega, 0, 2 * sizeof(double));
	// масса вертолетв
	 mHelicMass = 0.;

   //  memset(marrComGlEl, 0, 7 * sizeof(TCommonGliderElement));
     memset(marrJ, 0, 9 * sizeof(long double));
     mDn = 0.;
     marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
     mFiMax = 0.;
     mKappaTettaMax = 0.;
     mRuleAngMax = 0.;
     mbFullType = false;*/

}
// Конструктор копирования
 THelic::THelic (const THelic &R)
 {
     memcpy(marrRotor, R.marrRotor, 2 * sizeof(TRotor));
    // memcpy(marrRotorOmega, R.marrRotorOmega, 2 * sizeof(double));
     memcpy(marrComGlEl, R.marrComGlEl, 7 * sizeof(TCommonGliderElement));
     memcpy(marrJ, R.marrJ, 9 * sizeof(long double));
     mRuleGlEl = R.mRuleGlEl;
     mHelicMass  = R.mHelicMass ;
     mDn = R.mDn;
     mFiMax = R.mFiMax;
     mKappaTettaMax = R.mKappaTettaMax;
     mRuleAngMax = R.mRuleAngMax;
     mbFullType = R.mbFullType;
 }
 // оператор присваивания
  THelic THelic::operator=(THelic  R)
 {
      memcpy(marrRotor, R.marrRotor, 2 * sizeof(TRotor));

      memcpy(marrComGlEl, R.marrComGlEl, 7 * sizeof(TCommonGliderElement));
      memcpy(marrJ, R.marrJ, 9 * sizeof(long double));
      mRuleGlEl = R.mRuleGlEl;
      mHelicMass  = R.mHelicMass ;
      mDn = R.mDn;
      mFiMax = R.mFiMax;
      mKappaTettaMax = R.mKappaTettaMax;
      mRuleAngMax = R.mRuleAngMax;
      mbFullType = R.mbFullType;

     return *this ;
 }

  // парам констр 1
  THelic::THelic(TRotor *arrRotor, TCommonGliderElement *arrComGlEl
                 ,TRuleGliderElement RuleGlEl,const long double HelicMass
                 , const long double FiMax, const long  double KappaTettaMax)
  {
      memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
      memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
      mRuleGlEl = RuleGlEl;
      mHelicMass  = HelicMass ;
      memcpy(marrJ, ARrINertMtrx_, 9 * sizeof(long double));

      // доопределяем константы винтов из РЛЭ
      const long double TKel0 = 25.75 + 273.15;
      const long double HeliH = 4500.;
      mDn = 1.;//1.65;
      mFiMax = FiMax;
      marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

      marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
      ///

     // доопределяем константы лопастей на основе РЛЭ
     marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
     ///
     mKappaTettaMax = KappaTettaMax;
     mRuleAngMax = M_PI / 6.;
     mbFullType = false;
  }

 //---------------------------------------------------------------------------

  // парам констр 1 упрощенная аэродинамика и винт
  THelic::THelic(TRotor *arrRotor, TCommonGliderElement *arrComGlEl
                 ,TRuleGliderElement RuleGlEl,const long double HelicMass
                 , const long double FiMax, const long  double KappaTettaMax,const long  double RuleAngMax )
  {
      memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
      memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
      mRuleGlEl = RuleGlEl;
      mHelicMass  = HelicMass ;
      memcpy(marrJ, ARrINertMtrx_, 9 * sizeof(long double));

      // доопределяем константы винтов из РЛЭ
      const long double TKel0 = 25.75 + 273.15;
      const long double HeliH = 4500.;
      mDn = 1.;//1.65;
      mFiMax = FiMax;
      marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

      marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
      ///

     // доопределяем константы лопастей на основе РЛЭ
     marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
     ///
     mKappaTettaMax = KappaTettaMax;
     mRuleAngMax = RuleAngMax;
     mbFullType = false;
  }
/*
  //--------------------------------------------------------------------------------
  // парам конст создания Ка-52
  THelic::THelic(const bool bFullType, const bool bTypeofRotor)
  {
      mbFullType = bFullType;
      long double valZaclinAng = 0.;
      long double valShaftAxeX = 0.;
      if (!bFullType)
      { // упрощенная модель
       memcpy(marrJ, ARrINertMtrx_, 9 * sizeof(long double));
       valZaclinAng = constZaclinAng_;
       valShaftAxeX = constShaftAxeX_;
      }
      else
      {
          memcpy(marrJ, ARrINertMtrx, 9 * sizeof(long double));
          valZaclinAng = constZaclinAng;
          valShaftAxeX = constShaftAxeX;
      }
      // создание  лопасти винта4
           TBlade Blade (constBladeR, constRadHorizHsarnir, constPofile_d0, constPofile_d1, constBladeM, constBlade_b, constBladeCX0);

           // 0. Ввод файла с аэродинамическими характеристиками
         long double *arrInpDataGliders = new long double [NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData];
         memcpy(arrInpDataGliders, constArrHelicPlanersData_, NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData * sizeof(long double));
         ///
         double arrRotorOmega[2] ={0.};
         arrRotorOmega[0] = constRotorOmega * 2. * M_PI;
         arrRotorOmega[1] = -arrRotorOmega[0];


         TRotor arrRotor[2];
         arrRotor[0] = TRotor(Blade, constQuantBlades, valZaclinAng, valShaftAxeX, constShaftAxeY_,constShaftUpperL_,arrRotorOmega[0], bTypeofRotor);
         arrRotor[1] = TRotor(Blade, constQuantBlades, valZaclinAng, valShaftAxeX, constShaftAxeY_,constShaftLowerL_,  arrRotorOmega[1], bTypeofRotor);

         ///

         // 4. Вертолет
         // создание массива аэродинамическихэлементов планера
         TCommonGliderElement arrComGlEl[7];
         for (int i = 0; i < 7; i++)
         {
         arrComGlEl [i] = TCommonGliderElement(&(arrInpDataGliders [i* 13]));
         }

         TRuleGliderElement RuleGlEl  (&(arrInpDataGliders [7* 13]));
         ///


        // mHelic =  THelic(arrRotor, arrComGlEl
          //         ,RuleGlEl, constHelicMass, constFiMax, constFiMax);
         /// вертолет создан
         ///


         memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
         memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
         mRuleGlEl = RuleGlEl;
         mHelicMass = constHelicMass;


         // доопределяем константы винтов из РЛЭ
         const long double TKel0 = 25.75 + 273.15;
         const long double HeliH = 4500.;
         mDn = 1.;//1.65;
         mFiMax = constFiMax;
         marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

         marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
         ///

        // доопределяем константы лопастей на основе РЛЭ
        marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
        ///
        mKappaTettaMax = constFiMax;
        mRuleAngMax = M_PI / 6.;
        delete []arrInpDataGliders;
  }
*/

  // парам конст создания Ка-52
  THelic::THelic(const bool bFullType, const bool bTypeofRotor, const long double VAlHelicMass)
  {
      mbFullType = bFullType;
      mHelicMass = VAlHelicMass;
      long double valZaclinAng = 0.;
      // создание  лопасти винта4
       TBlade Blade (constBladeR, constRadHorizHsarnir, constPofile_d0, constPofile_d1, constBladeM, constBlade_b, constBladeCX0);
    // long double valShaftAxeX = 0.;
    // long double valShaftAxeY = 0.;
    // long double valShaftUpperL = 0.;
     //long double valShaftLowerL = 0.;
     long double arrT[9] = {0.};
     long double *arrInpDataGliders = new long double [NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData];
     TRotor arrRotor[2];
      if (!bFullType)
      { // упрощенная модель
          memcpy(arrT, ARrINertMtrx_, 9 * sizeof(long double));
          MatrxMultScalar(arrT, 3, 3, mHelicMass/constHelicMass ,marrJ);

       valZaclinAng = constZaclinAng_;

      memcpy(arrInpDataGliders, constArrHelicPlanersData_, NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData * sizeof(long double));
      arrRotor[0] = TRotor(Blade, constQuantBlades, valZaclinAng, constShaftAxeX_, constShaftAxeY_,constShaftUpperL_ ,constRotorOmega * 2. * M_PI, bTypeofRotor);
      arrRotor[1] = TRotor(Blade, constQuantBlades, valZaclinAng, constShaftAxeX_, constShaftAxeY_,constShaftLowerL_,  -constRotorOmega * 2. * M_PI, bTypeofRotor);


      }
      else
      {


          memcpy(arrT, ARrINertMtrxPrived, 9 * sizeof(long double));
          MatrxMultScalar(arrT, 3, 3, mHelicMass,marrJ);

         // valZaclinAng = constAlfaZaklNew;


        memcpy(arrInpDataGliders, constArrHelicPlanersData, NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData * sizeof(long double));
        arrRotor[0] =  TRotor(Blade, constQuantBlades, constAlfaZaklNew
                          ,constXdistLow + constDeltaXUp_Low, constRotorOmega * 2. * M_PI, bTypeofRotor);

        arrRotor[1] =  TRotor(Blade, constQuantBlades, constAlfaZaklNew
                          ,constXdistLow , -constRotorOmega * 2. * M_PI, bTypeofRotor);
      }






         ///

         // 4. Вертолет
         // создание массива аэродинамическихэлементов планера
         TCommonGliderElement arrComGlEl[7];
         for (int i = 0; i < 7; i++)
         {
         arrComGlEl [i] = TCommonGliderElement(&(arrInpDataGliders [i* 13]));
         }

         TRuleGliderElement RuleGlEl  (&(arrInpDataGliders [7* 13]));
         ///


         memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
         memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
         mRuleGlEl = RuleGlEl;



         // доопределяем константы винтов из РЛЭ
         const long double TKel0 = 25.75 + 273.15;
         const long double HeliH = 4500.;
         mDn = 1.;//1.65;
         mFiMax = constFiMax;
         marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

         marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
         ///

        // доопределяем константы лопастей на основе РЛЭ
        marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
        ///
        mKappaTettaMax = constFiMax;
        mRuleAngMax = M_PI / 6.;
        delete []arrInpDataGliders;
  }


  // парам конст создания Ка-52
  THelic::THelic(const bool bFullType, const bool bTypeofRotor, const long double VAlAlfaZakl, const long double VAlHelicMass)
  {
      mbFullType = bFullType;
      mHelicMass = VAlHelicMass;

      // создание  лопасти винта
       TBlade Blade (constBladeR, constRadHorizHsarnir, constPofile_d0, constPofile_d1, constBladeM, constBlade_b, constBladeCX0);

     long double arrT[9] = {0.};
     long double *arrInpDataGliders = new long double [NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData];
     TRotor arrRotor[2];
      if (!bFullType)
      { // упрощенная модель
          memcpy(arrT, ARrINertMtrx_, 9 * sizeof(long double));
          MatrxMultScalar(arrT, 3, 3, mHelicMass/constHelicMass ,marrJ);



      memcpy(arrInpDataGliders, constArrHelicPlanersData_, NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData * sizeof(long double));
      arrRotor[0] = TRotor(Blade, constQuantBlades, VAlAlfaZakl, constShaftAxeX_, constShaftAxeY_,constShaftUpperL_ ,constRotorOmega * 2. * M_PI, bTypeofRotor);
      arrRotor[1] = TRotor(Blade, constQuantBlades, VAlAlfaZakl, constShaftAxeX_, constShaftAxeY_,constShaftLowerL_,  -constRotorOmega * 2. * M_PI, bTypeofRotor);


      }
      else
      {


          memcpy(arrT, ARrINertMtrxPrived, 9 * sizeof(long double));
          MatrxMultScalar(arrT, 3, 3, mHelicMass,marrJ);

         // valZaclinAng = constAlfaZaklNew;


        memcpy(arrInpDataGliders, constArrHelicPlanersData, NUmRowsArrHelicPlanersData * NUmColsArrHelicPlanersData * sizeof(long double));
        arrRotor[0] =  TRotor(Blade, constQuantBlades, VAlAlfaZakl
                          ,constXdistLow + constDeltaXUp_Low, constRotorOmega * 2. * M_PI, bTypeofRotor);

        arrRotor[1] =  TRotor(Blade, constQuantBlades, VAlAlfaZakl
                          ,constXdistLow , -constRotorOmega * 2. * M_PI, bTypeofRotor);
      }






         ///

         // 4. Вертолет
         // создание массива аэродинамическихэлементов планера
         TCommonGliderElement arrComGlEl[7];
         for (int i = 0; i < 7; i++)
         {
         arrComGlEl [i] = TCommonGliderElement(&(arrInpDataGliders [i* 13]));
         }

         TRuleGliderElement RuleGlEl  (&(arrInpDataGliders [7* 13]));
         ///


         memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
         memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
         mRuleGlEl = RuleGlEl;



         // доопределяем константы винтов из РЛЭ
         const long double TKel0 = 25.75 + 273.15;
         const long double HeliH = 4500.;
         mDn = 1.;//1.65;
         mFiMax = constFiMax;
         marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

         marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
         ///

        // доопределяем константы лопастей на основе РЛЭ
        marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
        ///
        mKappaTettaMax = constFiMax;
        mRuleAngMax = M_PI / 6.;
        delete []arrInpDataGliders;
  }

/*
  // парам конст создания Ка-52
  THelic::THelic( const bool bTypeofRotor, long double *arrInpDataGliders)
  {
      mbFullType = true;
      mHelicMass = constHelicMass;
     long double valZaclinAng = 0.;
     long double valShaftAxeX = 0.;
     long double valShaftAxeY = 0.;
     long double valShaftUpperL = 0.;
     long double valShaftLowerL = 0.;

     long double arrT[9] = {0.};
     memcpy(arrT, ARrINertMtrxPrived, 9 * sizeof(long double));
     MatrxMultScalar(arrT, 3, 3, mHelicMass,marrJ);

     valZaclinAng = constZaclinAng;
     valShaftAxeX = constShaftAxeX;
     valShaftAxeY = constShaftAxeY;
     valShaftUpperL = constShaftUpperL;
     valShaftLowerL = constShaftLowerL;
      // создание  лопасти винта4
       TBlade Blade (constBladeR, constRadHorizHsarnir, constPofile_d0, constPofile_d1, constBladeM, constBlade_b, constBladeCX0);





       ///
      TRotor arrRotor[2];
         arrRotor[0] = TRotor(Blade, constQuantBlades, valZaclinAng, valShaftAxeX, valShaftAxeY,valShaftUpperL,constRotorOmega * 2. * M_PI, bTypeofRotor);
         arrRotor[1] = TRotor(Blade, constQuantBlades, valZaclinAng, valShaftAxeX, valShaftAxeY,valShaftLowerL,  -constRotorOmega * 2. * M_PI, bTypeofRotor);

         ///

         // 4. Вертолет
         // создание массива аэродинамическихэлементов планера
         TCommonGliderElement arrComGlEl[7];
         for (int i = 0; i < 7; i++)
         {
         arrComGlEl [i] = TCommonGliderElement(&(arrInpDataGliders [i* 13]));
         }

         TRuleGliderElement RuleGlEl  (&(arrInpDataGliders [7* 13]));
         ///


         memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
         memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
         mRuleGlEl = RuleGlEl;



         // доопределяем константы винтов из РЛЭ
         const long double TKel0 = 25.75 + 273.15;
         const long double HeliH = 4500.;
         mDn = 1.;//1.65;
         mFiMax = constFiMax;
         marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

         marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
         ///

        // доопределяем константы лопастей на основе РЛЭ
        marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
        ///
        mKappaTettaMax = constFiMax;
        mRuleAngMax = M_PI / 6.;

  }
*/
 //---------------------------------------------------------------------------

  // парам конст создания Ка-52
  // arrInpDataGliders - массив описывающий элеиенты планера относительно некоторого центра тяжести
  // arrCoordCentreMass - координаты нового центра тяжести в изначальной СвСК (
  // VAlAlfaZakl - угол заклинения. В новой СвСК цетры вращения винтов лежат на прямой имеющей в новой СвСК угол осью OY VAlAlfaZakl
  //VAlRotorUpDist - расстояние от центра масс до верхнего винта
  // VAlRotorLowDist - - расстояние от центра масс до нижнего винта
//  VAlCentrvka - это центровка. после того, как все сделано, начало координат каждого элемента планера и винтов смещается пол оси X
  // в точку VAlCentrvka
  //
  //сначала создается весь планер. затем смещается в arrCoordCentreMass
  //затем приделывается 2 винта
  //затем делается центровка
  THelic::THelic( const bool bTypeofRotor, const long double VAlMass, long double *arrINertMtrxPrived, long double *arrInpDataGliders
                  , long double* arrCoordCentreMass,const long double VAlAlfaZakl, const long double VAlRotorUpDist
                  , const long double VAlRotorLowDist,const long double VAlCentrvka)
  {
      mbFullType = true;

     MatrxMultScalar(arrINertMtrxPrived, 3, 3, VAlMass ,marrJ);

      // создание  лопасти винта
       TBlade Blade (constBladeR, constRadHorizHsarnir, constPofile_d0, constPofile_d1, constBladeM, constBlade_b, constBladeCX0);

       ///
      TRotor arrRotor[2];

      arrRotor[0] = TRotor(Blade, constQuantBlades, VAlAlfaZakl,VAlRotorUpDist,constRotorOmega * 2. * M_PI, bTypeofRotor);
      arrRotor[1] = TRotor(Blade, constQuantBlades, VAlAlfaZakl,VAlRotorLowDist,-constRotorOmega * 2. * M_PI, bTypeofRotor);
      arrRotor[0].mBasePLane.marrS0[0] -= VAlCentrvka;
      arrRotor[1].mBasePLane.marrS0[0] -= VAlCentrvka;

         ///

         // 4. Вертолет
         // создание массива аэродинамическихэлементов планера исходных
         TCommonGliderElement arrComGlEl[7];
         for (int i = 0; i < 7; i++)
         {
         arrComGlEl [i] = TCommonGliderElement(&(arrInpDataGliders [i* 13]));

         arrComGlEl[i].mPlaneSvSK.marrS0[0] -= arrCoordCentreMass[0] ;
         arrComGlEl[i].mPlaneSvSK.marrS0[1] -= arrCoordCentreMass[1];
         arrComGlEl[i].mPlaneSvSK.marrS0[2] -= arrCoordCentreMass[2];
         arrComGlEl[i].mPlaneSvSK.marrS0[0] -=VAlCentrvka;
         }

         TRuleGliderElement RuleGlEl  (&(arrInpDataGliders [7* 13]));
         RuleGlEl.mPlaneSvSK.marrS0[0] -= arrCoordCentreMass[0];
         RuleGlEl.mPlaneSvSK.marrS0[1] -= arrCoordCentreMass[1];
         RuleGlEl.mPlaneSvSK.marrS0[2] -= arrCoordCentreMass[2];
         RuleGlEl.mPlaneSvSK.marrS0[0] -=VAlCentrvka;

         ///


         memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
         memcpy(marrComGlEl, arrComGlEl, 7 * sizeof(TCommonGliderElement));
         mRuleGlEl = RuleGlEl;
         mHelicMass = VAlMass;


         // доопределяем константы винтов из РЛЭ
         const long double TKel0 = 25.75 + 273.15;
         const long double HeliH = 4500.;
         mDn = 1.;//1.65;
         mFiMax = constFiMax;
         marrRotor[0].mCt = marrRotor[1].mCt = calcCt( TKel0, HeliH);

         marrRotor[0].mKWave = marrRotor[1].mKWave = 0.5;
         ///

        // доопределяем константы лопастей на основе РЛЭ
        marrRotor[0].mBlade.mCyalfa = marrRotor[1].mBlade.mCyalfa = calc_C_y_alfa( TKel0, HeliH);
        ///
        mKappaTettaMax = constFiMax;
        mRuleAngMax = M_PI / 6.;

  }

 //---------------------------------------------------------------------------
// вычисление константы Ct по данным РЛЭ
// задана точка еа высоте висения (максимальной), масса вертолета, температура у поверхнсти земли
// требуется накйти константу Ct
    // TKel0 - температура на поверхности земли, град Кельвина
    // HeliH -  высота  висения, м
    // BladeOmega  - частота вращениеия винта
long double THelic::calcCt(const long double TKel0,const long  double HeliH)
{
      long double valTKel = TKel0  - 0.0065 * HeliH;
      long double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
      long double valRo =   valPAm * 0.028964/8.31/  valTKel;

      // радиус ометаемой площади несущего винта
      long double valR = marrRotor[0].mBlade.mBladeR;
      ///

      long double valCt =  9.81 * mHelicMass/
      (valRo   *M_PI * valR* valR *valR *valR
            *marrRotor[0].mOmega * marrRotor[0].mOmega)/ mFiMax;
      return valCt;
}

//--------------------------------------------------------------
//-------------------------------------------------------------
	// вычисление коэффициента подъемной силы лопасти  по данным РЛЭ
	// Задана точка висения при заданных данных - масса вертолета
	// , температура забортного воздуха и высота висения
	// На основании этих данных можно вычислить плотность воздуха, подъемную силу
	// и найти коэффициент Cyalfa
    // TKel0 - температура на поверхности земли, град Кельвина
	// HeliH -  высота  висения, м
	// BladeOmega  - частота вращениеия винта
	// максимальный общий шаг задан 15 град !!!!
long double THelic::calc_C_y_alfa(const long double TKel0,const long double HeliH)
{
	double valTKel = TKel0  - 0.0065 * HeliH;   
	double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
	double valRo =   valPAm * 0.028964/8.31/  valTKel;

	double valTemp = marrRotor[0].mBlade.mBladeR *  marrRotor[0].mBlade.mBladeR *  marrRotor[0].mBlade.mBladeR
	- marrRotor[0].mBlade.mRadHorizHsarnir * marrRotor[0].mBlade.mRadHorizHsarnir * marrRotor[0].mBlade.mRadHorizHsarnir;
	double valC_y_alfa = 6. * mHelicMass * 9.81
	/ (((double)(marrRotor[0].mQuantBlades)) *  valRo
    * marrRotor[0].mBlade.mBlade_b *mFiMax
    * marrRotor[0].mOmega *  marrRotor[0].mOmega * valTemp );
   return valC_y_alfa ;
}

//-------------------------------------------------------------
	// вычисление массовой характеристики лопасти  по данным РЛЭ
	// Задана точка висения при заданных данных - масса вертолета
	// , температура забортного воздуха и высота висения
	// На основании этих данных можно вычислить плотность воздуха, подъемную силу
	// и найти коэффициент Cyalfa и вычислить Gamma_L
    // TKel0 - температура на поверхности земли, град Кельвина
	// HeliH -  высота  висения, м
long double THelic::calc_Gamma_L(const long double TKel0,const long  double HeliH)
{
    double valC_y_alfa  = calc_C_y_alfa( TKel0,HeliH) ;
    double valTKel = TKel0  - 0.0065 * HeliH;
	double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
	double valRo =   valPAm * 0.028964/8.31/  valTKel;
	return  valC_y_alfa *  valRo * marrRotor[0].mBlade.mBlade_b
	*marrRotor[0].mBlade.mBladeR *  marrRotor[0].mBlade.mBladeR *  marrRotor[0].mBlade.mBladeR*  marrRotor[0].mBlade.mBladeR
	/marrRotor[0].mBlade.calcInertiaMoment()/2.;
}
/*
//-------------------------------------------------------------
// вычисление равнодействующей силы и момента силы
//INPUT:
// VAlAtmR - плотность воздуха
// arrSvSK_Va[3] - вектор воздушной скорости в СвСК
//OUTPUT:
// arrSvSK_Force[3] - результирующая сила
// arrSvSK_Moment[3] - результирующий момент
// arrSvSK_ShaftForce0[3] - результирующая сила винта
// arrSvSK_ShaftMoment0[3] - результирующий момент винта
// arrSvSK_AirForce0[3] - результирующая аэродинамическая сила
// arrSvSK_AirMoment0[3] - результирующий аэродинамический момент
//
void THelic::calcRezF_and_Moment_SvSK(double *arrOmegaSvSK, const double valFi, const  double   valKappa
                                      , const  double   valEtta, const  double   valDelFi
                                      ,const double VAlAtmRo, double *arrSvSK_Ua
                                      , double *arrSvSK_Force, double *arrSvSK_Moment
                                      , double * arrSvSK_ShaftForce0, double *arrSvSK_ShaftMoment0
                                      , double *arrSvSK_AirForce0, double *arrSvSK_AirMoment0)
{
    memset(arrSvSK_Force, 0, 3 * sizeof(double));
    memset(arrSvSK_Moment, 0, 3 * sizeof(double));

    double arrSvSK_ShaftForce[3] = {0.}, arrSvSK_ShaftMoment[3] = {0.};
    calcRezShaftF_and_Moment_SvSK(arrOmegaSvSK, valFi,   valKappa,   valEtta,    valDelFi
                      , VAlAtmRo, arrSvSK_Ua, arrSvSK_ShaftForce, arrSvSK_ShaftMoment);

    double arrSvSK_AirForce[3] = {0.}, arrSvSK_AirMoment[3] = {0.};
    calcRezAirF_and_Moment_SvSK(arrOmegaSvSK, valDelFi
                      , VAlAtmRo, arrSvSK_Ua, arrSvSK_AirForce, arrSvSK_AirMoment);

    MtrxSumMatrx(arrSvSK_ShaftForce, arrSvSK_AirForce,1, 3, arrSvSK_Force) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
    MtrxSumMatrx(arrSvSK_ShaftMoment, arrSvSK_AirMoment,1, 3, arrSvSK_Moment) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
    if (arrSvSK_ShaftForce0)
    {
        memcpy(arrSvSK_ShaftForce0, arrSvSK_ShaftForce, 3 * sizeof(double));
        memcpy(arrSvSK_ShaftMoment0, arrSvSK_ShaftMoment, 3 * sizeof(double));
        memcpy(arrSvSK_AirForce0, arrSvSK_AirForce, 3 * sizeof(double));
        memcpy(arrSvSK_AirMoment0, arrSvSK_AirMoment, 3 * sizeof(double));
    }

}

//-------------------------------------------------------------
// вычисление равнодействующей аэродинамич силы и момента аэродинамич силы
//INPUT:
// VAlAtmR - плотность воздуха
// arrSvSK_Va[3] - вектор воздушной скорости в СвСК
//OUTPUT:
// arrSvSK_Force[3] - результирующая сила
// arrSvSK_Moment[3] - результирующий момент

void THelic::calcRezAirF_and_Moment_SvSK(double *arrOmegaSvSK, const  double   valDelFi
                                      ,const double VAlAtmRo, double *arrSvSK_Ua
                                      , double *arrSvSK_AirForce, double *arrSvSK_AirMoment)
{
    memset (arrSvSK_AirForce, 0, 3 * sizeof(double));
    memset (arrSvSK_AirMoment, 0, 3 * sizeof(double));
   // вычисление массива аэродинамических сил и моментов этих сил для аэродинамических элементов планера
    double arrF[3] = {0.}, arrArmMoment[3] = {0.}, arrT[3] ={0.};

    TAbstractGliderElement *pAbstr[8];
    for (int i = 0; i < 7; i++)
    {
      pAbstr[i] = &(marrComGlEl[i]);
    }
    pAbstr[7] = &mRuleGlEl;
    for (int i = 0; i < 8; i++)
    {
     // pAbstr[i]->calAeroF_and_Mom(arrSvSK_Ua, VAlAtmRo,   valDelFi, arrF, arrArmMoment);// !!! НЕ ЗАБЫТЬ!!!
      pAbstr[i]->calAeroF_and_Mom(arrSvSK_Ua, VAlAtmRo,   0., arrF, arrArmMoment);


      MtrxSumMatrx(arrSvSK_AirForce, arrF,1, 3, arrT) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ АЭРОДИН СИЛЫ !!!
      memcpy(arrSvSK_AirForce, arrT, 3 * sizeof(double));


      MtrxSumMatrx(arrArmMoment, arrSvSK_AirMoment,1, 3, arrT) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
      memcpy(arrSvSK_AirMoment, arrT, 3 * sizeof(double));// ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
    }

}

//-------------------------------------------------------------
// вычисление равнодействующей силы винтов и момента силы винтов
//INPUT:
// VAlAtmR - плотность воздуха
// arrSvSK_Va[3] - вектор воздушной скорости в СвСК
//OUTPUT:
// arrSvSK_ShaftForce[3] - результирующая сила
// arrSvSK_ShaftMoment[3] - результирующий момент
void THelic::calcRezShaftF_and_Moment_SvSK(double *arrOmegaSvSK, const double valFi, const  double   valKappa
                                      , const  double   valEtta, const  double   valDelFi
                                      ,const double VAlAtmRo, double *arrSvSK_Ua
                                      , double *arrSvSK_ShaftForce, double *arrSvSK_ShaftMoment)
{
    memset (arrSvSK_ShaftForce, 0, 3 * sizeof(double));
    memset (arrSvSK_ShaftMoment, 0, 3 * sizeof(double));
    double arrF[3] = {0.}, arrArmMoment[3] = {0.}, arrT[3] ={0.};

    //  вычисление  аэродинамических сил и моментов этих сил для 2-х винтов
    for (int i = 0; i< 2; i++)
    {
        marrRotor[i].calcRezF_and_RezMom( valFi + valDelFi - 2.*((double)i) * valDelFi,    valKappa
                ,   valEtta,   arrSvSK_Ua,VAlAtmRo, arrOmegaSvSK,arrF, arrArmMoment);

        MtrxSumMatrx(arrSvSK_ShaftForce, arrF,1, 3, arrT) ;
        memcpy(arrSvSK_ShaftForce, arrT, 3 * sizeof(double));


        MtrxSumMatrx(arrArmMoment, arrSvSK_ShaftMoment,1, 3, arrT) ;  // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
        memcpy(arrSvSK_ShaftMoment, arrT, 3 * sizeof(double));
    }

}


//----------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
*/
//-------------------------------------------------------------
// вычисление равнодействующей силы и момента силы
//INPUT:
// VAlAtmR - плотность воздуха
// arrSvSK_Va[3] - вектор воздушной скорости в СвСК
//OUTPUT:
// arrSvSK_Force[3] - результирующая сила
// arrSvSK_Moment[3] - результирующий момент
// arrSvSK_ShaftForce0[3] - результирующая сила винта
// arrSvSK_ShaftMoment0[3] - результирующий момент винта
// arrSvSK_AirForce0[3] - результирующая аэродинамическая сила
// arrSvSK_AirMoment0[3] - результирующий аэродинамический момент
//
void THelic::calcRezF_and_Moment_SvSK(long double *arrOmegaSvSK, const long double valFi, const   long double   valKappa
                                      , const   long double   valEtta, const   long double   valDelFi
                                      ,const long double VAlAtmRo,long  double *arrSvSK_Ua
                                      , long double *arrSvSK_Force,long  double *arrSvSK_Moment
                                      , long double * arrSvSK_ShaftForce0,long  double *arrSvSK_ShaftMoment0
                                      , long double *arrSvSK_AirForce0,long double *arrSvSK_AirMoment0)
{
    memset(arrSvSK_Force, 0, 3 * sizeof(long double));
    memset(arrSvSK_Moment, 0, 3 * sizeof(long double));

    long double arrSvSK_ShaftForce[3] = {0.}, arrSvSK_ShaftMoment[3] = {0.};
    calcRezShaftF_and_Moment_SvSK(arrOmegaSvSK, valFi,   valKappa,   valEtta,    valDelFi
                      , VAlAtmRo, arrSvSK_Ua, arrSvSK_ShaftForce, arrSvSK_ShaftMoment);

    long double arrSvSK_AirForce[3] = {0.}, arrSvSK_AirMoment[3] = {0.};
    calcRezAirF_and_Moment_SvSK(arrOmegaSvSK, valDelFi
                      , VAlAtmRo, arrSvSK_Ua, arrSvSK_AirForce, arrSvSK_AirMoment);

    MtrxSumMatrx(arrSvSK_ShaftForce, arrSvSK_AirForce,1, 3, arrSvSK_Force) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
    MtrxSumMatrx(arrSvSK_ShaftMoment, arrSvSK_AirMoment,1, 3, arrSvSK_Moment) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
    if (arrSvSK_ShaftForce0)
    {
        memcpy(arrSvSK_ShaftForce0, arrSvSK_ShaftForce, 3 * sizeof(long double));
        memcpy(arrSvSK_ShaftMoment0, arrSvSK_ShaftMoment, 3 * sizeof(long double));
        memcpy(arrSvSK_AirForce0, arrSvSK_AirForce, 3 * sizeof(long double));
        memcpy(arrSvSK_AirMoment0, arrSvSK_AirMoment, 3 * sizeof(long double));
    }

}

//-------------------------------------------------------------
// вычисление равнодействующей аэродинамич силы и момента аэродинамич силы
//INPUT:
// VAlAtmR - плотность воздуха
// arrSvSK_Va[3] - вектор воздушной скорости в СвСК
//OUTPUT:
// arrSvSK_Force[3] - результирующая сила
// arrSvSK_Moment[3] - результирующий момент

void THelic::calcRezAirF_and_Moment_SvSK(long double *arrOmegaSvSK, const   long double   valDelFi
                                      ,const long double VAlAtmRo, long double *arrSvSK_Ua
                                      ,long  double *arrSvSK_AirForce, long double *arrSvSK_AirMoment)
{
    memset (arrSvSK_AirForce, 0, 3 * sizeof(long double));
    memset (arrSvSK_AirMoment, 0, 3 * sizeof(long double));
   // вычисление массива аэродинамических сил и моментов этих сил для аэродинамических элементов планера
    long double arrF[3] = {0.}, arrArmMoment[3] = {0.}, arrT[3] ={0.};

    TAbstractGliderElement *pAbstr[8];
    for (int i = 0; i < 7; i++)
    {
      pAbstr[i] = &(marrComGlEl[i]);
    }
    pAbstr[7] = &mRuleGlEl;

    long double valRuleAng = valDelFi * mRuleAngMax / mFiMax; // поворот руля направления

    for (int i = 0; i < 8; i++)
    {
       pAbstr[i]->calAeroF_and_Mom(arrSvSK_Ua, VAlAtmRo,   valRuleAng, arrF, arrArmMoment);//
   //   pAbstr[i]->calAeroF_and_Mom(arrSvSK_Ua, VAlAtmRo,   0., arrF, arrArmMoment);


      MtrxSumMatrx(arrSvSK_AirForce, arrF,1, 3, arrT) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ АЭРОДИН СИЛЫ !!!
      memcpy(arrSvSK_AirForce, arrT, 3 * sizeof(long double));


      MtrxSumMatrx(arrArmMoment, arrSvSK_AirMoment,1, 3, arrT) ; // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
      memcpy(arrSvSK_AirMoment, arrT, 3 * sizeof(long double));// ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
    }

}

//-------------------------------------------------------------
// вычисление равнодействующей силы винтов и момента силы винтов
//INPUT:
// VAlAtmR - плотность воздуха
// arrSvSK_Va[3] - вектор воздушной скорости в СвСК
//OUTPUT:
// arrSvSK_ShaftForce[3] - результирующая сила
// arrSvSK_ShaftMoment[3] - результирующий момент
void THelic::calcRezShaftF_and_Moment_SvSK(long double *arrOmegaSvSK, const long  double valFi, const   long double   valKappa
                                      , const   long double   valEtta, const   long double   valDelFi
                                      ,const double VAlAtmRo, long double *arrSvSK_Ua
                                      , long double *arrSvSK_ShaftForce,long  double *arrSvSK_ShaftMoment)
{
    memset (arrSvSK_ShaftForce, 0, 3 * sizeof(long double));
    memset (arrSvSK_ShaftMoment, 0, 3 * sizeof(long double));
    long double arrF[3] = {0.}, arrArmMoment[3] = {0.}, arrT[3] ={0.};

    //  вычисление  аэродинамических сил и моментов этих сил для 2-х винтов
    for (int i = 0; i< 2; i++)
    {
        marrRotor[i].calcRezF_and_RezMom( valFi + valDelFi - 2.*((long double)i) * valDelFi,    valKappa
                ,   valEtta,   arrSvSK_Ua,VAlAtmRo, arrOmegaSvSK,arrF, arrArmMoment);

        MtrxSumMatrx(arrSvSK_ShaftForce, arrF,1, 3, arrT) ;
        memcpy(arrSvSK_ShaftForce, arrT, 3 * sizeof(long double));


        MtrxSumMatrx(arrArmMoment, arrSvSK_ShaftMoment,1, 3, arrT) ;  // ДЛЯ/ ОТЛАДКИ УБРАНЫ МОМЕНТЫ!!!
        memcpy(arrSvSK_ShaftMoment, arrT, 3 * sizeof(long double));
    }

}

// изменение координат центра масс
// в текущей СвСК заданы координаты нового ценотра масс
// требуется пересчитать все векторы начал СвСКЭП в новую СвСК(с центром в новой точке).
// INPUT:
// VAlX0, VAlY0, VAlZ0 - координаты нового ыентра
void THelic::changeCentreMass(const long double VAlX0,const long double VAlY0,const long double VAlZ0 )
{
    long double arrCM_SvSK[3] = {0.}, parrRez[3] = {0.};
    arrCM_SvSK[0] = VAlX0;
    arrCM_SvSK[1] = VAlY0;
    arrCM_SvSK[2] = VAlZ0;
    MtrxMinusMatrx(mRuleGlEl.mPlaneSvSK.marrS0, arrCM_SvSK,1, 3, parrRez);
    memcpy(mRuleGlEl.mPlaneSvSK.marrS0, parrRez, 3 * sizeof(long double));
    for (int j = 0; j < 7; ++j)
    {
        MtrxMinusMatrx(marrComGlEl[j].mPlaneSvSK.marrS0, arrCM_SvSK,1, 3, parrRez);
        memcpy(marrComGlEl[j].mPlaneSvSK.marrS0, parrRez, 3 * sizeof(long double));

    }

    MtrxMinusMatrx(marrRotor[0].mBasePLane.marrS0, arrCM_SvSK,1, 3, parrRez);
    memcpy(marrRotor[0].mBasePLane.marrS0, parrRez, 3 * sizeof(long double));


    MtrxMinusMatrx(marrRotor[1].mBasePLane.marrS0, arrCM_SvSK,1, 3, parrRez);
    memcpy(marrRotor[1].mBasePLane.marrS0, parrRez, 3 * sizeof(long double));
}

