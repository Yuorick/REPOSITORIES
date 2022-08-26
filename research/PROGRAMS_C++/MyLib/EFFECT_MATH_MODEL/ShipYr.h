//---------------------------------------------------------------------------

#ifndef ShipYrH
#define ShipYrH
class TShipYr
{
public:


	double marrParral [3] ;//  вектор параллакса

	// парамеитры движения
	double mTShip ; // время привязки траекторной информации

	double mQ; // угол курса
	double mPsi; // угол килевой качки
	double mTet; // угол бортовой качки
	double mVQ; //скорость изменения угла курса
	double mVPsi ; // скорость изменения угла килевой качки
	double mVTet ; // вкорость изменения угла бортовой качуки


	// оценка вектора состояния корабля в ГСК
	double marrEstVectSost[9] ;

  /*
	 // ОТЧЕТ
	 // К-ВО ТОЧЕК В БУФУЕРЕ
	int mQuantPntReport ;
	 // ПАРАМЕТР ЗАРЕЗЕРВИРОВАННОЙ ПАМАЯТИ
	int mLenMemoryAlloc ;
	// БУФЕР ПАМЯТИ
	double *mparrBuff    ;
	//  ПУТЬ К ПАПКЕ С ОТЧЕТОМ
	wchar_t *mpwcharrFoldReport ;

	*/

	 __fastcall ~TShipYr() ;
	// конструктор по умолчанию
	TShipYr () ;
	// конструктор копирования
	TShipYr  (const TShipYr  &R) ;
	// оператор присваивания
	TShipYr  operator=(TShipYr   R2) ;

  // парам конструктор1
 TShipYr ( double *arrPar,const  double TShip
		 ,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet
		 , const double VShip, const double ZShip, const double ZVShip ) ;

 // парам конструктор 2
 TShipYr ( double *arrParral);

 void recalcVS(const  double valT
		 ,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet
		 , const double VShip, const double ZShip, const double ZVShip );

 void LinExtrap(const  double valT  , TShipYr &ShipYrExtr ) ;

	/*
	// парам конструктор1
	TShipYr (const TSins Sins, const TMeasurer Measurer, const TDriverMech Driver
		 ,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
		 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
		 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess,const  double TVess
		 , double *arrVectSost, double *arrEstVectSost,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet, double *arrDelt, wchar_t *pwcharrFoldReport);

     // парам конструктор3
   TShipYr (const double Bearing, const double TargCourse
	, const double TargZenitAng,  const double V, const double H ,
	const double R,const double valT, wchar_t *pwcharrFoldReport);

 // парам конструктор 4
 TShipYr (const TSins Sins, const TMeasurer Measurer, const TDriverMech Driver , const TEnvironment Environment
		 ,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
		 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
		 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
		 ,double *arrDelt ,const TInitTargData InitTargData, wchar_t *pwcharrFoldReport);


	 void calcAngles(const double valT);
	 void recalcVess(const double valT);
	 static void recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps) ;

	 void updateReportData() ;
	 void WriteReport() ;

	 void Move(const double valT,const double valStep) ;

	 void VSProlong(const double valTExtr, double *arrVSShipYrExtr)  ;

	 void  GetZamer_IN_KGSK (double *pTZam, double *arrZam_KGSK, double *arrKZam_KGSK ) ;

	  void  calcSummarizedCorMtrx_ErrMes_In_GSK (double *arrCorrMtrx_GSK );

	  void  GetZamer_IN_GSK (double *pTZam, double *arrZam_GSK, double *arrKZam_GSK );


   */



}  ;
#endif
