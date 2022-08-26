//---------------------------------------------------------------------------

#ifndef RotorH
#define RotorH
#include "LongPlane.h"
//#include "Plane.h"
#include "Blade.h"

// класс описывает винт вертолета
class TBlade;
class TLongPlane;
class TRotor
{
public:
	// лопасть винта
	TBlade mBlade;
	// к-во лопастей
	int mQuantBlades;
    // максимально допустимый общий шаг НВ
    //double mFiMax;
	// базовачя система координат привязанная к оси вращения вала  винта
    // ось Y направлена по оси вращения вала
     TLongPlane mBasePLane;

     // частота вращения
     long double mOmega;

     // коэффициент тяги винта CT
     long double mCt;
    // коэффициент компенсатора взмаха, константа,   0.6 > k > 0.4
    long double mKWave;
    // угол заклинения винта
    long double mZaklinAng;

    // плечо силы T от носительнооси Z БСК
   // long double mForceArmX;
    // плечо сил H и S от носительнооси оси X БСК
   // long double mForceArmY;

    // тип уравнений описывающих винт
    // если mbFullTypeoffEq == false, то учитываются только лишь 3 силы - H,T,S И их моменты
    // если mbFullTypeoffEq == true, то используется аорогная модель
    bool mbFullTypeoffEq;



  TRotor() ;
// Конструктор копирования
  TRotor (const TRotor &R2) ;

 // оператор присваивания
 TRotor   operator=(TRotor  R2) ;

 // TRotor(const TBlade Blade, const int QuantBlades, const long double VAlZaklinAng
 //              , const long double VAlCoordX, const long double VAlCoordY,const long double L,  const long double Omega);

  TRotor(const TBlade Blade, const int QuantBlades, const long double VAlZaklinAng
                  , const long double VAlCoordX, const long double VAlCoordY,const long  double L
                  , const long double Omega, const bool bFullType);

  TRotor(const TBlade Blade, const int QuantBlades, const long double VAlZaklinAng
                   ,const long  double L, const long double Omega, const bool bFullType);

  void calcRezF_and_RezMom(const long double VAlFi, const  long double VAl_DnKappa
                           ,const  long double VAl_DnEtta,long  double *arrSvSKUa, const long double VAlRo
                            ,long  double *arrOMegaSvSK,long double *arrRezSvSK_F,long  double *arrRezSvSK_Mom);

  long double calcAttackAng(long double *arrSvSKUa);

 long  double calcQreact(const long  double VAlRo, const long double Va);

  void  calcForceArms(long double *pForceArmX,long  double *pForceArmY);

  long double calcHorSharConst();

  void calCoef_a1_and_b1(const long double VAlFi, const long double VAlRo
                ,long double * arrSvSKUa, long double *pval_a1, long double *pval_b1);

  static long double Sign(long double a);

 long double calc_Deriv_T_po_Fi(const long double VAlRo);

 void calcRezF_and_RezMom_FullType(const long double VAlFi, const  long double VAl_DnKappa
       ,const long  double VAl_DnEtta ,long double *arrSvSKUa, const long  double VAlRo
        ,long  double *arrOMegaSvSK, long double *arrRezSvSK_F, long  double *arrRezSvSK_Mom);

 void calcRezF_and_RezMom_SimpleType(const long double VAlFi, const  long double VAl_DnKappa
      ,const long  double VAl_DnEtta ,long double *arrSvSKUa, const long  double VAlRo
       ,long  double *arrOMegaSvSK, long double *arrRezSvSK_F, long  double *arrRezSvSK_Mom);




};
#endif
