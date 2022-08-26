#include "Blade.h"


 TBlade::TBlade()
{
  // радиус винта
	 mBladeR  = 0.;
	// расстояние от оси вращения до  втулки ГШ
	 mRadHorizHsarnir =0.;
	// высота вертикального сечения  лопасти ближе к центру вращения
	 mPofile_d0 =0.;
	// высота вертикального сечения  лопасти дальше от центра вращения
	 mPofile_d1=0.;
	// масса лопасти
	 mBladeM=0.;
	// хорда лопасти
	 mBlade_b=0.;
     // коэффициент поджъемной силы
     mCyalfa=0.;
     // коэффициент лобового сопротивления
     mCX0 = 0.;;

     // статический момент
     mSg=0.;

     // момент инерции
      mIg=0.;

}
// Конструктор копирования
 TBlade::TBlade (const TBlade &R)
 {
     // радиус винта
      mBladeR  = R.mBladeR ;
     // расстояние от оси вращения до  втулки ГШ
      mRadHorizHsarnir = R.mRadHorizHsarnir;
     // высота вертикального сечения  лопасти ближе к центру вращения
      mPofile_d0 = R.mPofile_d0;
     // высота вертикального сечения  лопасти дальше от центра вращения
      mPofile_d1 = R.mPofile_d1;
     // масса лопасти
      mBladeM = R.mBladeM;
     // хорда лопасти
      mBlade_b = R.mBlade_b;
      // коэффициент подъемной силы
      mCyalfa = R.mCyalfa;
      // статический момент
      mSg = R.mSg;
      // момент инерции
       mIg = R. mIg;
       // коэффициент лобового сопротивления
       mCX0  = R.mCX0;
   }
 // оператор присваивания
  TBlade TBlade::operator=(TBlade  R)
 {
      // радиус винта
       mBladeR  = R.mBladeR ;
      // расстояние от оси вращения до  втулки ГШ
       mRadHorizHsarnir = R.mRadHorizHsarnir;
      // высота вертикального сечения  лопасти ближе к центру вращения
       mPofile_d0 = R.mPofile_d0;
      // высота вертикального сечения  лопасти дальше от центра вращения
       mPofile_d1 = R.mPofile_d1;
      // масса лопасти
       mBladeM = R.mBladeM;
      // хорда лопасти
       mBlade_b = R.mBlade_b;
       // коэффициент подъемной силы
       mCyalfa = R.mCyalfa;
       // статический момент
       mSg = R.mSg;
       // момент инерции
        mIg = R. mIg;
        // коэффициент лобового сопротивления
        mCX0  = R.mCX0;
     return *this ;
 }

 // парам констр 1
 TBlade::TBlade(const long double   BladeR,const long double   RadHorizHsarnir
   ,const long double   Pofile_d0,const long double   Pofile_d1,const long double   BladeM
                ,const long double   Blade_b, const long double  CX0)
 {
   mBladeR = BladeR;
   mRadHorizHsarnir = RadHorizHsarnir;
   mPofile_d0 = Pofile_d0 ;
   mPofile_d1 = Pofile_d1;
   mBladeM = BladeM;
   mBlade_b = Blade_b;
   mCyalfa = 0.; // будет определена в Helic.cpp
   mSg = calc_X_StatMoment_Sg() ;
   mIg = calcInertiaMomentHorSharnir();
   mCX0 = CX0;

 }
//-------------------------------------------------------------------------------
// вычисление расстояния от точки приложения аэродинамической силы
// сопротивления до центра втулки  оси НВ
long double  TBlade::calc_rQ()
{
   long double  x0 =  mRadHorizHsarnir;
   long double  x1 = mBladeR;
   long double  valK = (mPofile_d1  - mPofile_d0)/(x1 - x0);
   long double  temp =  mPofile_d0 - valK * x0;
   long double  temp2 = (x1 * x1 - x0 * x0) / 2.;
   long double  temp3 = (x1 * x1* x1 - x0 * x0* x0) / 3.;
   long double  temp4 = (x1 * x1* x1 * x1 - x0 * x0 * x0 * x0) / 4.;
   long double  val_rQ = (temp * temp3 + valK * temp4 )/ ( temp * temp2 + valK * temp3);
   return val_rQ;
}

//-------------------------------------------------------------
// вычисление расстояния от центра масс лопасти  до ближайшей точке лопасти к центру вращения
long double  TBlade::calc_X_BladeCentreMass()
{
    const long double  VAl_L = mBladeR - mRadHorizHsarnir;// VAl_L - длина лопасти
	return VAl_L * (mPofile_d0 + 2. * mPofile_d1)/ (mPofile_d0 +  mPofile_d1) /3.;
}

//-------------------------------------------------------------
// вычисление статического момента относитедбьно центра вращения
long double  TBlade::calc_X_StatMoment_Sg()
{
    long double  valCentreMass = calc_X_BladeCentreMass();
	return   mBladeM * ( valCentreMass + mRadHorizHsarnir);
}

//-------------------------------------------------------------
// вычисление момента инерции  относительно центра тяжести
long double  TBlade::calcInertiaMoment0()
{
    const long double  VAlL = mBladeR - mRadHorizHsarnir;// VAl_L - длина лопасти
	return   mBladeM * VAlL * VAlL *( mPofile_d0 * mPofile_d0 + mPofile_d0 * mPofile_d1 * 4.
	 + mPofile_d1* mPofile_d1)
	   / (mPofile_d0 + mPofile_d1) / (mPofile_d0 + mPofile_d1)/18.;
}

//-------------------------------------------------------------
// вычисление момента инерции  относительно центра втулки горизонтального шарнира
long double  TBlade::calcInertiaMomentHorSharnir()
{
    const long double  VAl_L = mBladeR - mRadHorizHsarnir;// VAl_L - длина лопасти
	return   mBladeM* VAl_L * VAl_L *( mPofile_d0  +  3. *mPofile_d1)
	   / (mPofile_d0 + mPofile_d1) /6.;
}

//-------------------------------------------------------------
// вычисление момента инерции  относительно центра вращения лопасти
long double   TBlade::calcInertiaMoment()
{
    long double  valXCentre = calc_X_BladeCentreMass() ;
	return   calcInertiaMoment0() + (valXCentre  +  mRadHorizHsarnir)
	 * (valXCentre  +  mRadHorizHsarnir) * mBladeM;
}



