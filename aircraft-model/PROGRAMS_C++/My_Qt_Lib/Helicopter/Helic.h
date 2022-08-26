//---------------------------------------------------------------------------

#ifndef HelicH
#define HelicH
#include "Rotor.h"
#include "CommonGliderElement.h"
#include "RuleGliderElement.h"

// класс описывает вертолет
// у вертолета может быть 2 винта
// в соосной схеме винт с номером 0 описывает верхний винт, а с номером 1 нижний
// в схеме МИ-8 винт с номером 0 описыает несущий винт, а с номером 1 рулевой
class TRotor;
class TCommonGliderElement;
class TRuleGliderElement;
class THelic
{
public:
    // тип вертолета
    // если mbFullType == true, то это точная модель
    // если mbFullType == false, то это упррощенная модель
    // она отличается от полной модели тем, что тензор инерции диагональный,
    // угол заклинения равен нулю, основание минта находится прямо над центром масс
    bool mbFullType;
	//  винты вертолета
	TRotor marrRotor[2];	

	// масса вертолетв
    long double mHelicMass;

    // аэродинамические элементы планера вертолета простого типа
    // marrComGlEl[0] - горизонтальная проекция фюзелфжа
    // marrComGlEl[1] - Передняя проекция фюзеляжа
    // marrComGlEl[2] - Вертикальная проекция фюзеляжа
    // marrComGlEl[3] - Правый стабилизатор
    // marrComGlEl[4] - Левый стабилизатор
    // marrComGlEl[5] - Правое крыло
    // marrComGlEl[6] - Левое крыло
    TCommonGliderElement marrComGlEl[7];
       // Воздушный руль
    TRuleGliderElement mRuleGlEl;

    // матрица момент инерции в СвСК
    long double marrJ[9];

    // коэффициент передачи циклического шага к углам
    long double mDn;

    // максимально допустимый общий шаг НВ
    long double mFiMax;
    // максимально допустимый циклический шаг НВ
    long double mKappaTettaMax;

    //максимально допустимый угол поворота руля направления
    long double mRuleAngMax;




  THelic() ;
// Конструктор копирования
  THelic (const THelic &R2) ;

 // оператор присваивания
 THelic  operator=(THelic  R2) ;

  // парам констр
 THelic(TRotor *arrRotor, TCommonGliderElement *arrComGlEl
                  ,TRuleGliderElement RuleGlEl,const long double HelicMass
                  , const long double FiMax, const long  double KappaTettaMax);

 THelic(TRotor *arrRotor, TCommonGliderElement *arrComGlEl
                  ,TRuleGliderElement RuleGlEl,const long double HelicMass
                  , const long double FiMax, const long  double KappaTettaMax,const long  double RuleAngMax );

 THelic( const bool bTypeofRotor, const long double VAlMass, long double *arrINertMtrxPrived, long double *arrInpDataGliders
                   , long double* arrCoordCentreMass,const long double VAlAlfaZakl, const long double VAlRotorUpDist
                   , const long double VAlRotorLowDist,const long double VAlCentrvka);



 THelic(const bool bFullType, const bool bTypeofRotor, const long double VAlHelicMass);

//THelic( const bool bTypeofRotor, long double *arrInpDataGliders);

 THelic(const bool bFullType, const bool bTypeofRotor, const long double VAlAlfaZakl, const long double VAlHelicMass);

long double calcCt(const long double TKel0,const long  double HeliH);


long double calc_C_y_alfa(const long double TKel0,const long double HeliH) ;

	// HeliH -  высота  висения, м
long double calc_Gamma_L(const long double TKel0,const long  double HeliH);
/*
void calcRezF_and_Moment_SvSK(double *arrOmegaSvSK, const double valFi, const  double   valKappa
                              , const  double   valEtta, const  double   valDelFi
                              ,const double VAlAtmRo, double *arrSvSK_Ua
                              , double *arrSvSK_Force, double *arrSvSK_Moment
                              , double * arrSvSK_ShaftForce0, double *arrSvSK_ShaftMoment0
                              , double *arrSvSK_AirForce0, double *arrSvSK_AirMoment0);

void calcRezAirF_and_Moment_SvSK(double *arrOmegaSvSK, const  double valDelFi
                                      ,const double VAlAtmRo, double *arrSvSK_Ua
                                      , double *arrSvSK_AirForce, double *arrSvSK_AirMoment);

void calcRezShaftF_and_Moment_SvSK(double *arrOmegaSvSK, const double valObSh, const  double valXB
                                      , const  double valXK, const  double valXH
                                      ,const double VAlAtmRo, double *arrSvSK_Ua
                                      , double *arrSvSK_ShaftForce, double *arrSvSK_ShaftMoment);
*/
void calcRezF_and_Moment_SvSK(long double *arrOmegaSvSK, const long double valFi, const   long double   valKappa
                                      , const   long double   valEtta, const   long double   valDelFi
                                      ,const long double VAlAtmRo,long  double *arrSvSK_Ua
                                      , long double *arrSvSK_Force,long  double *arrSvSK_Moment
                                      , long double * arrSvSK_ShaftForce0,long  double *arrSvSK_ShaftMoment0
                                      , long double *arrSvSK_AirForce0,long double *arrSvSK_AirMoment0);

void calcRezAirF_and_Moment_SvSK(long double *arrOmegaSvSK, const   long double   valDelFi
                                      ,const long double VAlAtmRo, long double *arrSvSK_Ua
                                      ,long  double *arrSvSK_AirForce, long double *arrSvSK_AirMoment);

void calcRezShaftF_and_Moment_SvSK(long double *arrOmegaSvSK, const long  double valFi, const   long double   valKappa
                                      , const   long double   valEtta, const   long double   valDelFi
                                      ,const double VAlAtmRo, long double *arrSvSK_Ua
                                      , long double *arrSvSK_ShaftForce,long  double *arrSvSK_ShaftMoment);

void changeCentreMass(const long double VAlX0,const long double VAlY0,const long double VAlZ0 );






};
#endif
