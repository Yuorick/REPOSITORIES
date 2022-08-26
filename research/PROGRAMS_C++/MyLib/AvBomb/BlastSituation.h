//---------------------------------------------------------------------------

#ifndef BlastSituationH
#define BlastSituationH
class TBlastSituation {
public:
	//угол визирования на участке подлета к цели
	double mAlfVis;
	// радиус поражения м
	double mRDefeat;
	// угол при вершине конуса отражателя
	double mFiDefeat;
	// СКЗ ошибки по дальности
	double mSigR;
	// СКЗ ошибки по углу  визирования
	double mSigAlf;
	// длина опасной части цели
	double mLDanger;



	// конструктор по умолчанию
	TBlastSituation();
	// конструктор копирования
	TBlastSituation(const TBlastSituation &R);

	// оператор присваивания
	TBlastSituation operator = (TBlastSituation R2);

	// парам конструктор
	TBlastSituation(const double AlfVis, const double RDefeat,const double FiDefeat
	 ,const double LDanger ,const double SigR
	 ,const double SigAlf
	 );

	 bool fncCreateGraph_Probab_FROM_FiDiagr(wchar_t  *pwchFileName
	 , const double valAngStep,const double valAngMin,const double valAngMax) ;


	 void findOptBlastTreshold_and_createProbabGraph(double &valDTresh
	,double &valProbTresh, wchar_t *pwchFileName);

	double calcProb(const double valDCur) ;
	double fncUslovnP(const double vald) ;
	static double fncNormRaspr( const double a, const double b);
	void createGrapUslovnP_from_D( wchar_t *pwchFileName) ;






};

#endif
