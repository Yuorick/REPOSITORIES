//---------------------------------------------------------------------------

#ifndef DriverMechH
#define DriverMechH
//---------------------------------------------------------------------------
class TDriverMech
{
public:
// �������� �� ���� �����
	double mSigBet;
//  �������� ����������� ���� �����
	double mSigEps;
// �������� ��������� ���� r�����
	double mSigDrBet;
// ��������  ��������� ���� �����
	double mSigDrEps;

// ������������ ����������
//protected:
  double mTDr ; // ������� ����� �������� ����������
  double mEstEps ;// ������(���������) ���� Eps
  double mEstBet ; // ������(���������) ���� Bet
  double mRealEps ;// �������� ���� Eps
  double mRealBet ; // �������� ���� Bet
  double mDelEps ;// ������ ��  Eps
  double mDelBet ; // ������ �� Bet
public:
  // �����

	 // �-�� ����� � �������
	int mQuantPntReport ;

	 //36. �������� ����������������� �������
	int mLenMemoryAlloc ;

	//37. ����� ������
	double *mparrBuff    ;

	// 38. ���� � ����� � �������
	wchar_t *mpwcharrFoldReport ;


	__fastcall ~TDriverMech() ;
	// ����������� �� ���������
	TDriverMech () ;
	// ����������� �����������
  	TDriverMech  (const TDriverMech  &R) ;
	// �������� ������������
	TDriverMech  &operator=(const TDriverMech   &R2) ;
	// ����� ������ 1
	TDriverMech( const double SigBet, const double SigEps, const double SigDrBet
		, const double SigDrEps, wchar_t *pwcharrFoldReport);

	 TDriverMech ( const double DriverSigBet // �������� ��������� ���� Bet �������
								,const double DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,const double DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,const double DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,const double TDr  // ������� ����� �������� ����������
								,const double VAlEstEps // ������(���������) ���� Eps
								,const double VAlEstBet  // ������(���������) ���� Bet
								,const double VAlRealEps // �������� ���� Eps
								,const double VAlRealBet  // �������� ���� Bet
								,const double VAlDelEps // ������ ��  Eps
								,const double VAlDelBet  // ������ �� Bet
								,wchar_t *pwcharrFoldReport);

	  void init (const double TDr  // ������� ����� �������� ����������
				,const double VAlEstEps // ������(���������) ���� Eps
				,const double VAlEstBet  // ������(���������) ���� Bet
				,const double VAlRealEps // �������� ���� Eps
				,const double VAlRealBet  // �������� ���� Bet
				,const double VAlDelEps // ������ ��  Eps
				,const double VAlDelBet  // ������ �� Bet
				) ;


	void recalcDriver(const double valT,const double valEps0,const double valBet0 ) ;

   void updateReportData() ;
   void WriteReport() ;
	void WriteReport(wchar_t *pwcharrPath);




}  ;

#endif
