//---------------------------------------------------------------------------

#ifndef TargBearing0H
#define TargBearing0H

// хранение данных по исходному пеленгу цели
class TTargBearing0
{
public:
	//  угол пеленга цели в ГСК
	double mBearing ;
	// угол курса цели в ГСК
	double mTargCourse ;
	// угол между ветором скорости цели и горизонтом
	double mTargHorizAng ;
	//- скорость
	double mV;
	// дальность,
	double mR ;
	// высота
	double mH ;
	//
  //	double mMy ;
	//
	// ЭПР
   //	double mTargEPR;
   //	enumTargetType mTargType;

	 __fastcall ~TTargBearing0() ;
	// конструктор по умолчанию
	TTargBearing0 () ;
	// конструктор копирования
	TTargBearing0  (const TTargBearing0  &R) ;
	// оператор присваивания
	TTargBearing0  operator=(TTargBearing0   R2) ;

	// парам конструктор1
	TTargBearing0 (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
	 ,const double  R,const double H);//, const double My , enumTargetType TargType );

  //	 TTargBearing0 (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
  //	 ,const double  R,const double H, const double My, const double TargEPR , enumTargetType TargType );

	 void raschet_nach_coord (double *pX) ;



}  ;
#endif
