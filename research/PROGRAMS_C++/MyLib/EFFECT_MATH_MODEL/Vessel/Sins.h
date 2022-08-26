//---------------------------------------------------------------------------

#ifndef SinsH
#define SinsH
class TEnvironment ;
class TSins

{
public:
     // КОНСТАНТЫ
	// ХАРАКТЕРИСТИКИ -  ДАННЫЕ В ТЗ
//  1.ограничение сверху на скз (точность)  по углу курса
	 double mMaxSig_Q;
//  2.ограничение сверху на скз (точность)  определения угла кормовой качки
	 double mMaxSig_Psi;
//  3.ограничение сверху на скз (точность)  определения угла бортовой качки
	double mMaxSig_Tet;
//  4.ограничение сверху на скз (точность)  определения скорости изменеия угла курса
	double mMaxSig_dQdt;
 // 5.ограничение сверху на скз (точность)  определения скорости изменеия угла кормовой качки
	double mMaxSig_dPsidt;
//  6.ограничение сверху на скз (точность) определения скорости изменеия угла бортовой качки
	double mMaxSig_dTetdt;
//  7.ограничение сверху на скз (точность) определения высоты  центра тяжести корабля
	double mMaxSig_H;
//  8.ограничение сверху на скз (точность) определения скорости изменения высоты  центра тяжести корабля
	double mMaxSig_VH;
	 // РЕАЛИЗОВАВШИЕСЯ ХАРАКТЕРИСТИКИ
//  9.ТОчность по углу курса
	double mSig_Q;
//  10.точность определения угла кормовой качки
	double mSig_Psi;
//  11.точность определения угла бортовой качки
	double mSig_Tet;
//  12.точность определения скорости изменеия угла курса
	double mSig_dQdt;
 // 13.точность определения скорости изменеия угла кормовой качки
	double mSig_dPsidt;
//  14.точность определения скорости изменеия угла бортовой качки
	double mSig_dTetdt;
//  15. точность определения высоты корабля
	double mSig_H ;
//  16. точность определения высоты корабля
	double mSig_VH ;
	///
 // 17. коэффициент ошибки определения скорости
	double mK1 ;
//  18. точность определения скорости
	double mSigV ;

///

// ДИНАМИЧЕСКАЯ ИНФОРМАЦИЯ
 // 19.Текущее время
  double mTSins ;
 //
  // 20.ошибка определения угла курса :
  double mDelQ;
  // 21.ошибка определения скорости ищзменения угла курса:
  double mDelVQ;
  // 22.ошибка определения угла килевой качки :
  double mDelPsi;
  // 23.ошибка определения скорости ищзменения угла килевой качки:
  double mDelVPsi;
  // 24.ошибка определения угла боротовой качки :
  double mDelTet;
  // 25. ошибка определения скорости ищзменения угла боротовой качки:
  double mDelVTet ;
  // 26.ошибка измерения высоты цетнтра корабля :
  double mDelH ;
  // 27.ошибка измерения скорости изменения  высоты цетнтра корабля :
  double mDelVH;
  // 28.ошибка измерения скорости корабля:
  double mDelVVess ;

// ОЦЕНКИ СИНС
   // 29.оценка ооценка угла курса :
  double mEstQ;
  // 30.оценкая скорости ищзменения угла курса:
  double mEstVQ;
  // 31.оценка угла килевой качки :
  double mEstPsi;
  // 32.оценка скорости ищзменения угла килевой качки:
  double mEstVPsi;
  // 33.оценка угла боротовой качки :
  double mEstTet;
  // 34. оценка скорости ищзменения угла боротовой качки:
  double mEstVTet ;
  // 35.оценка высоты цетнтра корабля :
  double mEstH ;
  // 36.оценка скорости изменения  высоты цетнтра корабля :
  double mEstVH;
  // 37.оценка скорости корабля:
  double mEstVVess ;

 // ОТЧЕТ ПО СИННС

	 // 38.максимальное возможное к-во измерений в франенте
	int mQuantPntReport ;

	 //39. текущее количество измерений ы фрагменте
	int mLenMemoryAlloc ;

	//40. массив измерений,  длина mQuantY *3  в ГСК
	double *mparrBuff    ;

	// 41. массив привязок измерений по времени, длина mQuantY
	wchar_t *mpwcharrFoldReport ;





		__fastcall ~TSins() ;


		__fastcall	 TSins () ;

		// конструктор копирования
		__fastcall  TSins  (const TSins  &R) ;
		// оператор присваивания
		TSins  &operator=(const TSins   &R2) ;
		// парам констр  1
		__fastcall TSins( const int QuantPntReport, const int LenMemoryAlloc
   , 	wchar_t *pwcharrFoldReport, double *parrBuff ) ;

   TSins (const TEnvironment Environment
		 ,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet
		 ,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt
		 ,const double MaxSig_H,const double MaxSig_VH,const double K1
			,const double SigV , wchar_t* pwcharrFoldReport);

	 // парам конструктор  3
 TSins (const TEnvironment Environment
	,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet
	,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt
	,const double MaxSig_H,const double MaxSig_VH,const double K1
	,const double SigV
	,const double  VAlT      // время
	,const double VAlTrueVVess   // ист скорость
	,const double  VAlTrueQ  // истинный угол курса
	,const double  VAlTruePsi  // истинный угол килевой качки
	,const double  VAlTrueTet  // истинный угол бортовой качки
	,const double  VAlTrueVQ  //истинный скорость изменения угла курса
	,const double  VAlTrueVPsi   //истинный скорость изменения угла килевой качки
	,const double  VAlTrueVTet   //истинный вкорость изменения угла бортовой качуки
	,const double VAlTrueH,const double VAlTrueVH,const double VAlT_Q
	,const double  VAlT_Psi,const double VAlT_Tet, double *arrDelt
	,wchar_t* pwcharrFoldReport);


  	void recalcSins_v0(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt);

   double callStepSys1(const double valX, const double h, const double valK
   ,const double valSigV) ;

   void recalcSins_v1(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt);

	 void WriteReport() ;
	 void WriteReport(wchar_t *pwcharrPath);

	 void fillValues_Delts_and_Ests(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt);

	  void updateReportData() ;




};
#endif
