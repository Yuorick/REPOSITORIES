//---------------------------------------------------------------------------


//---------------------------------------------------------------------------


#pragma hdrstop

#include "Bomb.h"
#include <string.h>
#include <math.h>
#include "Atmosphere.h"

 // коэффиц лобового сопротивления при М = 150/ 331
  const long double Cx00 =  0.2;
   // число аха для  mCx00
  const long double Maсh00 =150./ 331.;
 // коэффиц лобового сопротивления при М = 1.
  const long double Cx01 =  0.4;
   // число аха для  mCx00
  const long double  Maсh01  = 1.2;

  // коэффиц подъемной силы изодированного фюзеляжа (оживало + цилиндр)
  const long double CyFuz =  0.035;
  // коэфиициент подъемной силы Cу_корм = 0,8 * ( (D1/D)* (D1/D) -1) *2 * cos(tet1)  ** cos(tet1)
  // D = 0.45, D1 = 0.65,  tet1 = 30 град. Отнесен к единице площади миделя
  const long double CyKorma =  1.303;
   // табличное значение Cyкрыла при sqrt(1 -M*М) = 0.
  const long double Cy10 = 0.027;
 // коэффиц лобового сопротивления при sqrt(1 -M*М) = Lamb_k
  const long double Cy11 =  0.018;



//---------------------------------------------------------------------------
TBomb::TBomb()
{
 // площадь иделя фюзеляжа
	 mMidSq = M_PI * 0.58 * 0.58 / 4.;//M_PI * 0.45 * 0.45 / 4.;
// площадь боковой поверзности фюзеляжа ( M_PI * D * L)
	 mSurfFuzSq = M_PI * sqrt( 4. * mMidSq / M_PI)* 4.63 ;
// масса
	mMass = 2200.;
// макс угол атаки:
	mMaxAt = 20. * M_PI / 180. ;
 // коэфф удлинения крыла L*L/S= 1.5*1.5/ (1.5*0.75)
	mLamb_k = 2.25/4.5;
 // плошщадь крвыльев S = 3. * 1.5
	mSq_k = 4.5;  //
 // признак наличия ГСН
   //	mbHomeHead =true ;
 // ГСН
	 mHomingHead = THomingHead();
 // ЭМБ
	 mEbomb = TEbomb();
 //  пост времени системы стабилизации по углу атаки
	 mTStabSyst = 0.5;

}

//---------------------------------------------------------------------------


// конструктор копирования
 TBomb ::TBomb (const TBomb &R)
 {
	 mMidSq = R.mMidSq ;
	 mSurfFuzSq = R.mSurfFuzSq;
	 mMass = R.mMass;
	 mMaxAt = R.mMaxAt ;
	 mSq_k = R.mSq_k ;
	 mLamb_k = R.mLamb_k ;
	 mHomingHead = R.mHomingHead ;
	 mEbomb = R.mEbomb ;
	 mTStabSyst = R.mTStabSyst ;

 }

 // оператор присваивания
 TBomb TBomb::operator=(TBomb  R)
 {
	 mMidSq = R.mMidSq ;
	 mSurfFuzSq = R.mSurfFuzSq;
	 mMass = R.mMass;
	 mMaxAt = R.mMaxAt ;
	 mSq_k = R.mSq_k ;
	 mLamb_k = R.mLamb_k ;
	 mHomingHead = R.mHomingHead ;
	 mEbomb = R.mEbomb ;
	 mTStabSyst = R.mTStabSyst ;

	return *this ;
 }
// параметрическийц конструктор 1
 TBomb::TBomb(
	 long double MidSq // площадь иделя фюзеляжа
	,long double SurfFuzSq // площадь боковой поверзности фюзеляжа
	,long double Mass// масса
	,long double MaxAt // макс угол атаки
	,long double Lamb_k // уцдлинение крыльев
	,long double Sq_k // площадб крыльев, м кв
	,bool bHomeHead // признак наличия ГСН
	,THomingHead HomingHead // ГСН
	,TEbomb Ebomb// ЭМБ
	)
{
	mMidSq = MidSq;
	mSurfFuzSq = SurfFuzSq ;
	mMass = Mass ;
	mMaxAt = MaxAt ;
	mLamb_k = Lamb_k ;
	mSq_k  = Sq_k  ;
	mHomingHead = HomingHead ;
	mEbomb  = Ebomb ;

}
// параметрическийц конструктор 2
 TBomb::TBomb(THomingHead HomingHead )
{
	*this =  TBomb();
	mHomingHead = HomingHead ;

}

// параметрическийц конструктор 3
 TBomb::TBomb(THomingHead HomingHead,TEbomb Ebomb )
{
	*this =  TBomb();
	mHomingHead = HomingHead ;
	mEbomb = Ebomb ;
}

// параметрическийц конструктор 4
 TBomb::TBomb(THomingHead HomingHead,TEbomb Ebomb, const long double TStabSyst )
{
	*this =  TBomb();
	mHomingHead = HomingHead ;
	mEbomb = Ebomb ;
	mTStabSyst = TStabSyst;
}



// расчет коэфф лобового сопрот
 long double TBomb::fncCx( long double valMach)
{
	if (valMach >= Maсh01 )
	 return  Cx01;

	if (valMach >= Maсh00 ) return  Cx00 + (Cx01 - Cx00)/(Maсh01 - Maсh00) * (valMach - Maсh00);
	return  Cx00;
}


// расчет производной коэфф лобового сопрот  по числу Маха
 long double TBomb::fnc_dCx_po_dMach( long double valMach)
{
	if (valMach >= Maсh01 ) return  0.;

	 return  2. *(Cx01 - Cx00)/(Maсh01 - Maсh00) ;

}

// расчет коэфф подьемной силы
 long double TBomb::fncCy( long double valMach)
{
 /*  try
   {
	long double vala = 0.;
	if (valMach >= 1. )
	{
	  vala =  sqrt ( valMach * valMach - 1. ) ;
	}
	else
	{
	 vala = sqrt ( 1. -  valMach * valMach) ;
    }

	return 1.42 *((CyFuz * mSurfFuzSq/ mMidSq + CyKorma +  mSq_k/ mMidSq * mLamb_k * (Cy10 + (Cy11 - Cy10) * vala)));
   //	return 1. *((CyFuz * mSurfFuzSq/ mMidSq + CyKorma +  mSq_k/ mMidSq * mLamb_k * (Cy10 + (Cy11 - Cy10) * vala)));
	}
	catch(...)
	{
		int iii = 0;
	} */
	return 4.5;

}

// расчет производной коэфф подьемной силы  по числу Маха
 long double TBomb::fnc_dCy_po_dMach( long double valMach)
{
 /* try
  {
	if (valMach >= 1. ) return  0.;
	double vala = sqrt ( 1. -  valMach * valMach) ;

	return 1.42*( -mSq_k/ mMidSq * mLamb_k * (Cy11 - Cy10) *valMach /sqrt ( 1. -  valMach * valMach)) ;
   //return 1.*( -mSq_k/ mMidSq * mLamb_k * (Cy11 - Cy10) *valMach /sqrt ( 1. -  valMach * valMach)) ;
  }
  catch(...)
  {
	  int ii = 0;
  }  */
  return 0.;

}


//---------------------------------------------------------------------------
#pragma package(smart_init)
