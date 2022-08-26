//---------------------------------------------------------------------------

#ifndef ShipYrH
#define ShipYrH
class TShipYr
{
public:


	double marrParral [3] ;//  ������ ����������

	// ���������� ��������
	double mTShip ; // ����� �������� ����������� ����������

	double mQ; // ���� �����
	double mPsi; // ���� ������� �����
	double mTet; // ���� �������� �����
	double mVQ; //�������� ��������� ���� �����
	double mVPsi ; // �������� ��������� ���� ������� �����
	double mVTet ; // �������� ��������� ���� �������� ������


	// ������ ������� ��������� ������� � ���
	double marrEstVectSost[9] ;

  /*
	 // �����
	 // �-�� ����� � �������
	int mQuantPntReport ;
	 // �������� ����������������� �������
	int mLenMemoryAlloc ;
	// ����� ������
	double *mparrBuff    ;
	//  ���� � ����� � �������
	wchar_t *mpwcharrFoldReport ;

	*/

	 __fastcall ~TShipYr() ;
	// ����������� �� ���������
	TShipYr () ;
	// ����������� �����������
	TShipYr  (const TShipYr  &R) ;
	// �������� ������������
	TShipYr  operator=(TShipYr   R2) ;

  // ����� �����������1
 TShipYr ( double *arrPar,const  double TShip
		 ,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet
		 , const double VShip, const double ZShip, const double ZVShip ) ;

 // ����� ����������� 2
 TShipYr ( double *arrParral);

 void recalcVS(const  double valT
		 ,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet
		 , const double VShip, const double ZShip, const double ZVShip );

 void LinExtrap(const  double valT  , TShipYr &ShipYrExtr ) ;

	/*
	// ����� �����������1
	TShipYr (const TSins Sins, const TMeasurer Measurer, const TDriverMech Driver
		 ,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
		 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
		 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess,const  double TVess
		 , double *arrVectSost, double *arrEstVectSost,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet, double *arrDelt, wchar_t *pwcharrFoldReport);

     // ����� �����������3
   TShipYr (const double Bearing, const double TargCourse
	, const double TargZenitAng,  const double V, const double H ,
	const double R,const double valT, wchar_t *pwcharrFoldReport);

 // ����� ����������� 4
 TShipYr (const TSins Sins, const TMeasurer Measurer, const TDriverMech Driver , const TEnvironment Environment
		 ,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
		 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
		 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
		 ,double *arrDelt ,const TInitTargData InitTargData, wchar_t *pwcharrFoldReport);


	 void calcAngles(const double valT);
	 void recalcVess(const double valT);
	 static void recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps) ;

	 void updateReportData() ;
	 void WriteReport() ;

	 void Move(const double valT,const double valStep) ;

	 void VSProlong(const double valTExtr, double *arrVSShipYrExtr)  ;

	 void  GetZamer_IN_KGSK (double *pTZam, double *arrZam_KGSK, double *arrKZam_KGSK ) ;

	  void  calcSummarizedCorMtrx_ErrMes_In_GSK (double *arrCorrMtrx_GSK );

	  void  GetZamer_IN_GSK (double *pTZam, double *arrZam_GSK, double *arrKZam_GSK );


   */



}  ;
#endif
