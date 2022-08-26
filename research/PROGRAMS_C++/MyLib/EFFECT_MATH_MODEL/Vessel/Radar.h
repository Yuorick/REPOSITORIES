//---------------------------------------------------------------------------

#ifndef RadarH
#define RadarH
//---------------------------------------------------------------------------
class TRadar
{
public:
// ТОчность по углу курса
	double mSigV;
//  точность определения угла места
	double mSigU;
// точность определения дальности (шаг дискретизации)
	double mSigR;
// тактовая частота

	double mh ; // интервал между измерениями- нужен для стационарнфых расчетов

	double mRzv; //  замер дальности
	double mVzv;  // замер по углу V (курс)
	double mUzv ; // замер по углу U (угол места)
	double mT ; // время привязки замера
	double mDelR ;// ошибка измерения дальности
	double mDelV ; // ошибка ищзмерения угла V
	double mDelU ; // ошибка измерения угда U


	// ОТЧЕТ

	 // К-ВО ТОЧЕК В БУФУЕРЕ
	int mQuantPntReport ;

	 // ПАРАМЕТР ЗАРЕЗЕРВИРОВАННОЙ ПАМАЯТИ
	int mLenMemoryAlloc ;

	// БУФЕР ПАМЯТИ
	double *mparrBuff    ;

	//  ПУТЬ К ПАПКЕ С ОТЧЕТОМ
	wchar_t *mpwcharrFoldReport ;

	 __fastcall ~TRadar() ;
	// конструктор по умолчанию
	TRadar () ;
	// конструктор копирования
	TRadar  (const TRadar  &R) ;
	// оператор присваивания
	TRadar  operator=(TRadar   R2) ;

	// парам конструктор1
	 TRadar::TRadar (const double SigV, const double SigU, const double SigR ,const double T,const double h, wchar_t *pwcharrFoldReport);

	// парам конструктор2
	TRadar (const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU , wchar_t *pwcharrFoldReport);


	  void createMeasure (const double valR,const double valV, const double valU, const double valT) ;

   void updateReportData() ;
   void WriteReport() ;
   void WriteReport(wchar_t *pwcharrPath) ;




}  ;
#endif
