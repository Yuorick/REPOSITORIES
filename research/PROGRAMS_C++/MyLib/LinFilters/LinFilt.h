//---------------------------------------------------------------------------

#ifndef LinFiltH
#define LinFiltH
class TMeasure;
class TComp;
class TLinFilt
{
public:
	// ����������� ������� ��������� �������
	int mDimX;

	// ����������� ������� ����������
	int mDimY;
	// ������������ ����������
//  1.����� �������� ������ ������� ���������
	double mTf;
//  2.������ ������� ��������� ����  �� ������ mTf
	double *mparrEstX;

 // 3.�������������� ������� ������ ����������
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
