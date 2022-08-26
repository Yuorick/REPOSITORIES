#include "ElectMotor.h"
#include "math.h"
#include "Equations.h"
#include <wchar.h>
#include <stdio.h>
#include <time.h>
#include "MatrixProccess.h"
#include "Comp.h"
#include <string.h>
#include "Gauss.h"

// коэффициент лобового сопротивления пластины
//#define VAL_CX 0.3
// плотность атмосферы
#define RO_ATM 1.22
// к-во пар обмоток
extern const double CONST_ZP = 8.;
extern const double CONST_DBM_185_10_03_2_OmegaMax = 33.5;
extern const double CONST_DBM_185_10_03_2_UMax = 60.;
extern const double CONST_DBM_185_10_03_2_AMax = 220.;
//extern const double CONST_MOM_DRY_FRICTION = 0.;// 0.01;
//extern const double CONST_GARMON_COEFF = 0.9;
extern const double CONST_J_ROTOR = 0.015;

//---------------------------------------------------------------------------

   QElectMotor:: QElectMotor()
{
    // индуктивность катушки статора
    mInductL = 0.;
    // сопротивление катушки статора
    mResist = 0.;

    // индуктивный поток
    mPsi_f = 0.;
    // момент инерции ротора
    mJ0 = 0.;

    // к-во пар
    mZp = 8;
    //   
    mUMax = 0.;
    // предельно допустимая амплитуда тока
    mAMax = 0.;
   // предельно допустимая частота вращения
   mOmegaMax = 0.;

   mMomDryFriction = 0.;   

   mGrmnPh0 = 0.;

   mInvL = 0.;

  mMaxAmpGrmn = 0.;
  mGrmnPh0 = 0.;
  mMaxDispRezid = 0.;
  mPracticalMaxDispRezid = 0.;


}

//---------------------------------------------------------------------------

// конструктор копирования
     QElectMotor :: QElectMotor (const  QElectMotor &R)
 {
         mInductL = R.mInductL;
         mResist = R.mResist ;
         mPsi_f = R.mPsi_f;
         mJ0 = R.mJ0;         
         mZp = R.mZp;
         mUMax = R.mUMax;
         mAMax = R.mAMax;
         mOmegaMax = R.mOmegaMax;
         mMomDryFriction = R.mMomDryFriction;         
         mGrmnPh0 = R.mGrmnPh0;
         mMaxAmpGrmn = R.mMaxAmpGrmn;
         mMaxDispRezid = R.mMaxDispRezid;
         mInvL = R.mInvL;
         mPracticalMaxDispRezid = R.mPracticalMaxDispRezid;

 }

 // оператор присваивания
  QElectMotor  &QElectMotor::operator=( const QElectMotor  &R)
 {
      mInductL = R.mInductL;
      mResist = R.mResist ;
      mPsi_f = R.mPsi_f;
      mJ0 = R.mJ0;
      mZp = R.mZp;
      mUMax = R.mUMax;
      mAMax = R.mAMax;
      mOmegaMax = R.mOmegaMax;
      mMomDryFriction = R.mMomDryFriction;
      mGrmnPh0 = R.mGrmnPh0;
      mMaxAmpGrmn = R.mMaxAmpGrmn;
      mMaxDispRezid = R.mMaxDispRezid;
      mInvL = R.mInvL;
      mPracticalMaxDispRezid = R.mPracticalMaxDispRezid;

     return *this ;
 }

//--------------------------------------------------------------------
  // парам конструктор
 /*    QElectMotor:: QElectMotor (const double  InductL , const double  Resist ,
      const double  Psi_f , const double   JRotor , const double MaxAmpGrmn
      , const double  GrmnPh0,const double MaxDispMomRezid, const double PracticalMaxDispRezid)
 {
    mInductL = InductL;
    mResist  = Resist ;
    mPsi_f  = Psi_f ;
    mJ0 = JRotor;
    mMaxAmpGrmn = MaxAmpGrmn;
    mZp = CONST_ZP;
    mOmegaMax = CONST_DBM_185_10_03_2_OmegaMax ;
    mUMax =  CONST_DBM_185_10_03_2_UMax;
    mAMax = CONST_DBM_185_10_03_2_AMax ;
    mMomDryFriction = CONST_MOM_DRY_FRICTION;
    mGrmnPh0 = GrmnPh0;
    mInvL = 1./ InductL ;
    mMaxDispRezid = MaxDispMomRezid;
    mPracticalMaxDispRezid = PracticalMaxDispRezid;

 }*/
     // парам конструктор
QElectMotor:: QElectMotor (const double  InductL , const double  Resist ,
         const double  Psi_f , const double   JRotor , const double MaxAmpGrmn
         , const double  GrmnPh0,const double MaxDispMomRezid, const double PracticalMaxDispRezid, const double MomDryFriction)
    {
       mInductL = InductL;
       mResist  = Resist ;
       mPsi_f  = Psi_f ;
       mJ0 = JRotor;
       mMaxAmpGrmn = MaxAmpGrmn;
       mZp = CONST_ZP;
       mOmegaMax = CONST_DBM_185_10_03_2_OmegaMax ;
       mUMax =  CONST_DBM_185_10_03_2_UMax;
       mAMax = CONST_DBM_185_10_03_2_AMax ;
       mMomDryFriction = MomDryFriction;
       mGrmnPh0 = GrmnPh0;
       mInvL = 1./ InductL ;
       mMaxDispRezid = MaxDispMomRezid;
       mPracticalMaxDispRezid = PracticalMaxDispRezid;

    }

//-----------------------------------------
//Вычисление Psi_f на основании "Основных технических характеристик двигателя"
//VAlZp = 8 - к-во пар
// VAlCoeff = 0,8-0,87 - коэффициент момента Н*М/А
// показывает сколько момента можно получить с 1 ампера
// Psi_f = 2/3 * VAlCoeff/ VAlZp =  0.0198- 0.01815
 double    QElectMotor:: calc_Psi_f (const double  VAlZp, const double   VAlCoeff )
 {
     return 2./3. * VAlCoeff/ VAlZp;
 }

 //----------------------------------------
 double    QElectMotor:: calcSumMomNoise (const double VAlTetta, const double VAlOmega)
 {
     double valMomGarmon = calc_HarmMom(VAlTetta,  VAlOmega);
     double valDispNoise =  mPracticalMaxDispRezid  * VAlOmega * VAlOmega /mOmegaMax * mOmegaMax;
     return -mMomDryFriction *SIGNUM_(VAlOmega) +valMomGarmon +  getGauss(0., sqrt(valDispNoise));
 }
//----------------------------
 double QElectMotor::SIGNUM_(const double a)
 {
   if (fabs(a) < 0.0000000001)
   {
       return 0.;
   }
   return (a > 0.)?1.:-1.;
 }
 //-----------------------------------------------------
 //----------------------------------------
 double    QElectMotor:: calcDisp_MomNoise ( const double VAlOmega)
 {
     double temp = calc_AmpHarmMom ( 0.,  VAlOmega);
     double temp1 = calcDisp_FluctMomNoise ( VAlOmega);
     return temp * temp/ 2.  +temp1;
             //( mMaxAmpGrmn * mMaxAmpGrmn /2. + mMaxDispRezid) * VAlOmega  * VAlOmega /(mOmegaMax * mOmegaMax);

 }
 //----------------------------------------
 double    QElectMotor:: calcDisp_FluctMomNoise ( const double VAlOmega)
 {
     return  mMaxDispRezid * VAlOmega  * VAlOmega /(mOmegaMax * mOmegaMax);

 }

 //----------------------------------------
 double    QElectMotor:: getMaxHarmonAmp ()
 {
     return   mMaxAmpGrmn;

 }
 //----------------------------------------
 double    QElectMotor:: getHarmonPh ()
 {
     return   mGrmnPh0 ;

 }
 //----------------------------------------
 double    QElectMotor:: getLInv ()
 {
     return   mInvL ;

 }
 //----------------------------------------
 double    QElectMotor:: getJ0 ()
 {
     return   mJ0 ;

 }
 //----------------------------------------
 //----------------------------------------
 double    QElectMotor:: getL ()
 {
     return   mInductL ;

 }
 //----------------------------------------
 double    QElectMotor:: getMomDryFriction ()
 {
     return   mMomDryFriction ;

 }
 void   QElectMotor:: setLInv (const double Linv)
 {
    mInductL =  1./Linv ;
    mInvL = Linv;

 }


 //----------------------------------------
 void   QElectMotor:: setMaxAmpGrmn (const double MaxAmpGrmn)
 {
        mMaxAmpGrmn = MaxAmpGrmn;

 }
 //----------------------------------------
 void   QElectMotor:: setMaxDispMomRezid (const double MaxDispMomRezid)
 {
        mMaxDispRezid =  MaxDispMomRezid;

 }
 //----------------------------------------
 void    QElectMotor:: setHarmonPh (const double GrmnPh0)
 {
       mGrmnPh0 = GrmnPh0;

 }
 //---------------------------------------
 // производгная нармониечкого момента по Om
 double QElectMotor:: calc_dHarmMom_po_dOm (const double VAlTetta)
 {
     return mMaxAmpGrmn /mOmegaMax * cos (mZp * VAlTetta + mGrmnPh0);
 }
 //---------------------------------------
 // производгная нармониечкого момента по Tetta
 double QElectMotor:: calc_dHarmMom_po_dTetta (const double VAlTetta, const double VAlOm)
 {
     double temp = calc_AmpHarmMom ( VAlTetta,  VAlOm);
     return -temp * mZp * sin (mZp * VAlTetta + mGrmnPh0);
 }
 //---------------------------------------
 // производгная нармонического момента по амплитуде
 double QElectMotor:: calc_dHarmMom_po_dAmp (const double VAlTetta, const double VAlOm)
 {
     return VAlOm/mOmegaMax * cos (mZp * VAlTetta + mGrmnPh0);
 }

 //---------------------------------------
 // производгная нармонического момента по фазе
 double QElectMotor:: calc_dHarmMom_po_dPh (const double VAlTetta, const double VAlOm)
 {
     double temp = calc_AmpHarmMom ( VAlTetta,  VAlOm);
     return -temp * sin (mZp * VAlTetta + mGrmnPh0);
 }
 //---------------------------------------
 // вычисление гармонического момента
 double QElectMotor:: calc_HarmMom (const double VAlTetta, const double VAlOm)
 {
     double temp = calc_AmpHarmMom ( VAlTetta,  VAlOm);
     return temp * cos (mZp * VAlTetta + mGrmnPh0);
 }

 // ----------------------------------------------------------
 // вычисление дисперсии случайной составляющей остаточного момента
  double QElectMotor:: calc_MaxDispMomRezid (const double MaxAmpGrmn)
  {
     return MaxAmpGrmn * MaxAmpGrmn *(1. - 0.9 * 0.9/2.)/ 9.;// = 0.26 * 0.26 *A *A
  }
  //---------------------------------------
  // вычисление гармонического момента
  double QElectMotor:: calc_AmpHarmMom (const double VAlTetta, const double VAlOm)
  {
      return mMaxAmpGrmn /mOmegaMax * VAlOm  ;
  }
  //--------------------------------------------
  int QElectMotor::createInputDataReport(wchar_t*FileName, const bool bHeader)
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
  fprintf(fw,"  Электромотор\n");


  fprintf(fw,"  индуктивность катушки статора(мГ) = %5.4f\n",mInductL * 1000.);
  fprintf(fw,"  момент инерции ротора (кг*м*м) = %6.4f\n",mJ0);
  fprintf(fw,"  макс. амп. гармон. сост. остаточного момента = %6.4f\n",mMaxAmpGrmn);
  fprintf(fw,"  нач. фаза гармонической составляющей = %6.4f\n",mGrmnPh0);
  fprintf(fw,"  макс. дисп. шума остаточного момента теор. = %8.6f\n",mMaxDispRezid);
  fprintf(fw,"  макс. дисп. шума остаточного момента практ. = %6.4f\n",mPracticalMaxDispRezid);
  fprintf(fw,"  момент сухого трения = %6.4f\n",mMomDryFriction);
  fprintf(fw,"  сопротивление катушки статора (Ом) = %6.4f\n",mResist);
  fprintf(fw,"  индуктивный поток = %6.4f\n",mPsi_f);
  fprintf(fw,"  к-во пар = %3.1f\n",(double)mZp);
  fprintf(fw,"  предельно допустимое напряжение (В) = %6.1f\n",mUMax);
  fprintf(fw,"  предельно допустимая амплитуда тока (А) = %6.4f\n",mAMax);
  fprintf(fw,"  предельно допустимая частота вращения (рад/с) = %6.4f\n",mOmegaMax);

  fclose(fw);
  }



