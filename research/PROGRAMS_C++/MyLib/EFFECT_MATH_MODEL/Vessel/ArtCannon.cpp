//---------------------------------------------------------------------------


#pragma hdrstop

#include "ArtCannon.h"


TArtCannon::TArtCannon()
{
 menumCannonType = A_192M;
 mRateOfFire = 0.;
 mAngGroupedFire =0.;
 mSigRabT = 0.;
 mRabT=0.;
}



// конструктор копирования
 TArtCannon ::TArtCannon (const TArtCannon &R)
 {
	menumCannonType = R.menumCannonType;
	mRateOfFire  =  R.mRateOfFire ;
	mAngGroupedFire = R.mAngGroupedFire;
	mSigRabT = R.mSigRabT;
	mRabT=R.mRabT;

 }
 // оператор присваивания
 TArtCannon &TArtCannon::operator=(const TArtCannon  &R)
 {
  menumCannonType = R.menumCannonType;
	mRateOfFire  =  R.mRateOfFire ;
	mAngGroupedFire = R.mAngGroupedFire;
	mSigRabT = R.mSigRabT;
	mRabT=R.mRabT;

	return *this ;
 }

  // парам конструктор1
 TArtCannon::TArtCannon (const double RateOfFire, const double AngGroupedFire
  ,const double SigmaDelayT, const double RabT)
 {
	mRateOfFire  = RateOfFire  ;
	mAngGroupedFire = AngGroupedFire ;
	mSigRabT = SigmaDelayT;
	mRabT = RabT;

 }

 // парам конструктор1
 TArtCannon::TArtCannon (enumCannonType EnumCannonType, const double SigmaDelayT
 , const double RabT)
 {
	menumCannonType = EnumCannonType;
	mSigRabT = SigmaDelayT;
	mRabT = RabT;
	switch( EnumCannonType)
	{
		case A_192M:
		mRateOfFire = 2.;

		mAngGroupedFire = 0.001;
		break;
		case A_190_01:
		 mRateOfFire = 60. / 75.;

		 mAngGroupedFire = 0.001;
		break;
		case AK_630:
		mRateOfFire = 60. / 4980.;

		mAngGroupedFire = 0.0003;
		break;
		case AK_630_2M:
		mRateOfFire = 60. / 4980./ 2.;

		mAngGroupedFire = 0.003;
		break;

		case AK_176:
		mRateOfFire = 60./125. ;

		mAngGroupedFire = 0.001;
		break;

		default:
		break;

	}

 }
#pragma package(smart_init)
