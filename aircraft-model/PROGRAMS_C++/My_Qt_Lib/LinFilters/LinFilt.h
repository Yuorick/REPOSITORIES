//---------------------------------------------------------------------------

#ifndef LinFiltH
#define LinFiltH
class TMeasure;
class TComp;
class TLinFilt
{
public:
	// размерность вектора состояния объекта
	int mDimX;

	// размерность вектора наблюдений
	int mDimY;
	// Динамическая информация
//  1.Время привязки оценок вектора состояния
	double mTf;
//  2.Оценка вектора состояния угла  на момент mTf
	double *mparrEstX;

 // 3.корреляционная матрица ошибок оценивания
	double *mparrMtrxK;

	TLinFilt();
	~TLinFilt();
	TLinFilt(const int DimX, const int DimY, const double ValT,   double *parrEstX, double *parrMtrxK);
	virtual void createMtrxA(const TMeasure Measure,const double Val_h, double *pMtrxA) = 0;//{};

	virtual void createMtrxCorrU(const TMeasure Measure,const double Val_h,const double Val_Tau2, double *pMtrxCorrU) = 0;//{};


	virtual void createMtrxC(const TMeasure Meas,const double Val_h, double *pMtrxC) = 0;//{};




	virtual	int fncStep( TMeasure Measure, const double Val_Tau2);

	void fncStabSolutionDynamic(const double Val_h, const double Val_Tau2
	 ,const TMeasure Measure);

	 virtual  void fncAnalyticSolution(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT){};
};
#endif
