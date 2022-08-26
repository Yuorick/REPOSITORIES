//---------------------------------------------------------------------------

#ifndef OEChannelH
#define OEChannelH
//---------------------------------------------------------------------------
class TOEChannel
{
public:
// �������� �� ���� �����
	double mSigV;
//  �������� ����������� ���� �����
	double mSigU;
// �������� ����������� ��������� (��� �������������)
	double mSigR;
// �������� �������

	double mh ; // �������� ����� �����������
 // ������ ��������� (���)
	double mDiagWidth ;


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

	 __fastcall ~TOEChannel() ;
	// ����������� �� ���������
	TOEChannel () ;
	// ����������� �����������
	TOEChannel  (const TOEChannel  &R) ;
	// �������� ������������
	TOEChannel  operator=(TOEChannel   R2) ;

	// ����� �����������1
	 TOEChannel::TOEChannel (const double SigV, const double SigU, const double SigR ,const double T,const double h
 , const double DiagWidth, wchar_t *pwcharrFoldReport);

	// ����� �����������2
	TOEChannel (const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU
	, const double DiagWidth  , wchar_t *pwcharrFoldReport);


	  void createMeasure (const double valR,const double valV, const double valU, const double valT) ;

   void updateReportData() ;
   void WriteReport() ;
   void WriteReport(wchar_t *pwcharrPath)  ;




}  ;
#endif
