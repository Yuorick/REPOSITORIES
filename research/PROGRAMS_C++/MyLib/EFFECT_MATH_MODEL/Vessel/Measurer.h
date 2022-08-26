//---------------------------------------------------------------------------

#ifndef MeasurerH
#define MeasurerH
#include "OEChannel.h"
 #include "Radar.h"
//---------------------------------------------------------------------------
//class TRadar;
//class TOEChannel;
class TMeasurer
{
public:

	// ТОчность по углу курса
	double mSigV;
//  точность определения угла места
	double mSigU;
// точность определения дальности (шаг дискретизации)
	double mSigR;

	double mRzv; //  замер дальности
	double mVzv;  // замер по углу V (курс)
	double mUzv ; // замер по углу U (угол места)
	double mT ; // время привязки замера
	double mDelR ;// ошибка измерения дальности
	double mDelV ; // ошибка ищзмерения угла V
	double mDelU ; // ошибка измерения угда U
// private:
 TRadar mRadar;
 TOEChannel mOEChannel;
 // признак включения ОЭК в состав измерительных средств:
	bool mbOEChannel;
 // е6сли в формировании замера ОЭК не принимал участие, то   mAlf = 1
 // если принимал, то mAlf = 0.5
 // mAlf нужна для обнаружения маневра типа кабрирование
  double mAlf;

	// ОТЧЕТ

	 // К-ВО ТОЧЕК В БУФУЕРЕ
	int mQuantPntReport ;

	 // ПАРАМЕТР ЗАРЕЗЕРВИРОВАННОЙ ПАМАЯТИ
	int mLenMemoryAlloc ;

	// БУФЕР ПАМЯТИ
	double *mparrBuff    ;

	//  ПУТЬ К ПАПКЕ С ОТЧЕТОМ
	wchar_t *mpwcharrFoldReport ;

	 __fastcall ~TMeasurer() ;
	// конструктор по умолчанию
	TMeasurer () ;
	// конструктор копирования
	TMeasurer  (const TMeasurer  &R) ;
	// оператор присваивания
	TMeasurer  operator=(TMeasurer   R2) ;


	TMeasurer (const double SigV, const double SigU, const double SigR ,const double T,const double h, wchar_t *pwcharrFoldReport);
	// парам конструктор
	TMeasurer (const double SigV, const double SigU, const double SigR ,const double T,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU
	,const bool bOEChannel , const TRadar Radar, const TOEChannel OEChannel, wchar_t *pwcharrFoldReport);

	 // парам конструк   2
	TMeasurer(const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU , wchar_t *pwcharrFoldReport) ;
     // парам конструктор 4
	  TMeasurer ( const TRadar Radar, const TOEChannel OEChannel, wchar_t *pwcharrFoldReport) ;


	  void createMeasure (const double valR,const double valV, const double valU, const double valT) ;

   void updateReportData() ;
   void WriteReport() ;
   void takeRadarData() ;
   void uniteDatas() ;
   static void uniteTwoMeasures(const double valMeasure0, const double valSig0, const double valDel0
  ,const double valMeasure1, const double valSig1,  const double valDel1
  , double *pMeasureRez,  double *pSigRez,  double *pDelRez) ;

   void WriteReport(wchar_t *pwcharrPath)  ;





}  ;
#endif
