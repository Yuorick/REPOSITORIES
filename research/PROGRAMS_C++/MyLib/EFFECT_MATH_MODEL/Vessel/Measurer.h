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

	// �������� �� ���� �����
	double mSigV;
//  �������� ����������� ���� �����
	double mSigU;
// �������� ����������� ��������� (��� �������������)
	double mSigR;

	double mRzv; //  ����� ���������
	double mVzv;  // ����� �� ���� V (����)
	double mUzv ; // ����� �� ���� U (���� �����)
	double mT ; // ����� �������� ������
	double mDelR ;// ������ ��������� ���������
	double mDelV ; // ������ ���������� ���� V
	double mDelU ; // ������ ��������� ���� U
// private:
 TRadar mRadar;
 TOEChannel mOEChannel;
 // ������� ��������� ��� � ������ ������������� �������:
	bool mbOEChannel;
 // �6��� � ������������ ������ ��� �� �������� �������, ��   mAlf = 1
 // ���� ��������, �� mAlf = 0.5
 // mAlf ����� ��� ����������� ������� ���� ������������
  double mAlf;

	// �����

	 // �-�� ����� � �������
	int mQuantPntReport ;

	 // �������� ����������������� �������
	int mLenMemoryAlloc ;

	// ����� ������
	double *mparrBuff    ;

	//  ���� � ����� � �������
	wchar_t *mpwcharrFoldReport ;

	 __fastcall ~TMeasurer() ;
	// ����������� �� ���������
	TMeasurer () ;
	// ����������� �����������
	TMeasurer  (const TMeasurer  &R) ;
	// �������� ������������
	TMeasurer  operator=(TMeasurer   R2) ;


	TMeasurer (const double SigV, const double SigU, const double SigR ,const double T,const double h, wchar_t *pwcharrFoldReport);
	// ����� �����������
	TMeasurer (const double SigV, const double SigU, const double SigR ,const double T,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU
	,const bool bOEChannel , const TRadar Radar, const TOEChannel OEChannel, wchar_t *pwcharrFoldReport);

	 // ����� ��������   2
	TMeasurer(const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU , wchar_t *pwcharrFoldReport) ;
     // ����� ����������� 4
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
