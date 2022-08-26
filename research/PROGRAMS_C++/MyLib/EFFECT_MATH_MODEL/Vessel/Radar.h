//---------------------------------------------------------------------------

#ifndef RadarH
#define RadarH
//---------------------------------------------------------------------------
class TRadar
{
public:
// �������� �� ���� �����
	double mSigV;
//  �������� ����������� ���� �����
	double mSigU;
// �������� ����������� ��������� (��� �������������)
	double mSigR;
// �������� �������

	double mh ; // �������� ����� �����������- ����� ��� ������������� ��������

	double mRzv; //  ����� ���������
	double mVzv;  // ����� �� ���� V (����)
	double mUzv ; // ����� �� ���� U (���� �����)
	double mT ; // ����� �������� ������
	double mDelR ;// ������ ��������� ���������
	double mDelV ; // ������ ���������� ���� V
	double mDelU ; // ������ ��������� ���� U


	// �����

	 // �-�� ����� � �������
	int mQuantPntReport ;

	 // �������� ����������������� �������
	int mLenMemoryAlloc ;

	// ����� ������
	double *mparrBuff    ;

	//  ���� � ����� � �������
	wchar_t *mpwcharrFoldReport ;

	 __fastcall ~TRadar() ;
	// ����������� �� ���������
	TRadar () ;
	// ����������� �����������
	TRadar  (const TRadar  &R) ;
	// �������� ������������
	TRadar  operator=(TRadar   R2) ;

	// ����� �����������1
	 TRadar::TRadar (const double SigV, const double SigU, const double SigR ,const double T,const double h, wchar_t *pwcharrFoldReport);

	// ����� �����������2
	TRadar (const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU , wchar_t *pwcharrFoldReport);


	  void createMeasure (const double valR,const double valV, const double valU, const double valT) ;

   void updateReportData() ;
   void WriteReport() ;
   void WriteReport(wchar_t *pwcharrPath) ;




}  ;
#endif
