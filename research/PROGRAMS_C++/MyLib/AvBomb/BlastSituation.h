//---------------------------------------------------------------------------

#ifndef BlastSituationH
#define BlastSituationH
class TBlastSituation {
public:
	//���� ����������� �� ������� ������� � ����
	double mAlfVis;
	// ������ ��������� �
	double mRDefeat;
	// ���� ��� ������� ������ ����������
	double mFiDefeat;
	// ��� ������ �� ���������
	double mSigR;
	// ��� ������ �� ����  �����������
	double mSigAlf;
	// ����� ������� ����� ����
	double mLDanger;



	// ����������� �� ���������
	TBlastSituation();
	// ����������� �����������
	TBlastSituation(const TBlastSituation &R);

	// �������� ������������
	TBlastSituation operator = (TBlastSituation R2);

	// ����� �����������
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
