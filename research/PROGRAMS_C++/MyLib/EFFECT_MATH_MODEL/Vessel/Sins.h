//---------------------------------------------------------------------------

#ifndef SinsH
#define SinsH
class TEnvironment ;
class TSins

{
public:
     // ���������
	// �������������� -  ������ � ��
//  1.����������� ������ �� ��� (��������)  �� ���� �����
	 double mMaxSig_Q;
//  2.����������� ������ �� ��� (��������)  ����������� ���� �������� �����
	 double mMaxSig_Psi;
//  3.����������� ������ �� ��� (��������)  ����������� ���� �������� �����
	double mMaxSig_Tet;
//  4.����������� ������ �� ��� (��������)  ����������� �������� �������� ���� �����
	double mMaxSig_dQdt;
 // 5.����������� ������ �� ��� (��������)  ����������� �������� �������� ���� �������� �����
	double mMaxSig_dPsidt;
//  6.����������� ������ �� ��� (��������) ����������� �������� �������� ���� �������� �����
	double mMaxSig_dTetdt;
//  7.����������� ������ �� ��� (��������) ����������� ������  ������ ������� �������
	double mMaxSig_H;
//  8.����������� ������ �� ��� (��������) ����������� �������� ��������� ������  ������ ������� �������
	double mMaxSig_VH;
	 // ��������������� ��������������
//  9.�������� �� ���� �����
	double mSig_Q;
//  10.�������� ����������� ���� �������� �����
	double mSig_Psi;
//  11.�������� ����������� ���� �������� �����
	double mSig_Tet;
//  12.�������� ����������� �������� �������� ���� �����
	double mSig_dQdt;
 // 13.�������� ����������� �������� �������� ���� �������� �����
	double mSig_dPsidt;
//  14.�������� ����������� �������� �������� ���� �������� �����
	double mSig_dTetdt;
//  15. �������� ����������� ������ �������
	double mSig_H ;
//  16. �������� ����������� ������ �������
	double mSig_VH ;
	///
 // 17. ����������� ������ ����������� ��������
	double mK1 ;
//  18. �������� ����������� ��������
	double mSigV ;

///

// ������������ ����������
 // 19.������� �����
  double mTSins ;
 //
  // 20.������ ����������� ���� ����� :
  double mDelQ;
  // 21.������ ����������� �������� ���������� ���� �����:
  double mDelVQ;
  // 22.������ ����������� ���� ������� ����� :
  double mDelPsi;
  // 23.������ ����������� �������� ���������� ���� ������� �����:
  double mDelVPsi;
  // 24.������ ����������� ���� ��������� ����� :
  double mDelTet;
  // 25. ������ ����������� �������� ���������� ���� ��������� �����:
  double mDelVTet ;
  // 26.������ ��������� ������ ������� ������� :
  double mDelH ;
  // 27.������ ��������� �������� ���������  ������ ������� ������� :
  double mDelVH;
  // 28.������ ��������� �������� �������:
  double mDelVVess ;

// ������ ����
   // 29.������ ������� ���� ����� :
  double mEstQ;
  // 30.������� �������� ���������� ���� �����:
  double mEstVQ;
  // 31.������ ���� ������� ����� :
  double mEstPsi;
  // 32.������ �������� ���������� ���� ������� �����:
  double mEstVPsi;
  // 33.������ ���� ��������� ����� :
  double mEstTet;
  // 34. ������ �������� ���������� ���� ��������� �����:
  double mEstVTet ;
  // 35.������ ������ ������� ������� :
  double mEstH ;
  // 36.������ �������� ���������  ������ ������� ������� :
  double mEstVH;
  // 37.������ �������� �������:
  double mEstVVess ;

 // ����� �� �����

	 // 38.������������ ��������� �-�� ��������� � ��������
	int mQuantPntReport ;

	 //39. ������� ���������� ��������� � ���������
	int mLenMemoryAlloc ;

	//40. ������ ���������,  ����� mQuantY *3  � ���
	double *mparrBuff    ;

	// 41. ������ �������� ��������� �� �������, ����� mQuantY
	wchar_t *mpwcharrFoldReport ;





		__fastcall ~TSins() ;


		__fastcall	 TSins () ;

		// ����������� �����������
		__fastcall  TSins  (const TSins  &R) ;
		// �������� ������������
		TSins  &operator=(const TSins   &R2) ;
		// ����� ������  1
		__fastcall TSins( const int QuantPntReport, const int LenMemoryAlloc
   , 	wchar_t *pwcharrFoldReport, double *parrBuff ) ;

   TSins (const TEnvironment Environment
		 ,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet
		 ,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt
		 ,const double MaxSig_H,const double MaxSig_VH,const double K1
			,const double SigV , wchar_t* pwcharrFoldReport);

	 // ����� �����������  3
 TSins (const TEnvironment Environment
	,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet
	,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt
	,const double MaxSig_H,const double MaxSig_VH,const double K1
	,const double SigV
	,const double  VAlT      // �����
	,const double VAlTrueVVess   // ��� ��������
	,const double  VAlTrueQ  // �������� ���� �����
	,const double  VAlTruePsi  // �������� ���� ������� �����
	,const double  VAlTrueTet  // �������� ���� �������� �����
	,const double  VAlTrueVQ  //�������� �������� ��������� ���� �����
	,const double  VAlTrueVPsi   //�������� �������� ��������� ���� ������� �����
	,const double  VAlTrueVTet   //�������� �������� ��������� ���� �������� ������
	,const double VAlTrueH,const double VAlTrueVH,const double VAlT_Q
	,const double  VAlT_Psi,const double VAlT_Tet, double *arrDelt
	,wchar_t* pwcharrFoldReport);


  	void recalcSins_v0(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt);

   double callStepSys1(const double valX, const double h, const double valK
   ,const double valSigV) ;

   void recalcSins_v1(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt);

	 void WriteReport() ;
	 void WriteReport(wchar_t *pwcharrPath);

	 void fillValues_Delts_and_Ests(const double valT, const double VVess,	const double Q
	,const double Psi,	const double Tet,const double VQ
	,const double VPsi, const double VTet,const double H,const double VH,const double T_Q
	,const double  T_Psi,const double T_Tet, double *arrDelt);

	  void updateReportData() ;




};
#endif
