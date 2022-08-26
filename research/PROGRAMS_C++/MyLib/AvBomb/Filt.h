//---------------------------------------------------------------------------

#ifndef FiltH
#define FiltH

class TShipTarg;
class TBombTraj;
//---------------------------------------------------------------------------
class TFilt
{
public:



 // ������������ ����������
//  1.����� �������� ������ ������� ���������
	long double mTf;
//  2.������ ������� ��������� ���� ������������ ����  � ��� �� ������ mTf
	long double marrEstAngVis[2];

 // 3.�������������� ������� ������ ���������� � ���
	long double marrK_AngVis [4];

  //  4.����� ������� ��������� ��������� ����  � ��� �� ������ mTf
	long double marrEstR[2];

 // 5. �������������� ������� ������ ���������� ��������� � ���
	long double marrK_R [4];

// 6. ������ ������� ��������� ����  � ��� �� ��� X �� ������ mTf
	 long double marrEstTargX[2];

// 7.�������������� ������� ������ ���������� ������� ��������� ����  � ��� �� ��� X �� ������ mTf
	long double marrK_TargX [4];

 // 6. ������� ����, ��� ������ ������ �������������
	bool mbInit ;






	// ����������� �� ���������
	TFilt () ;
	// ����������� �����������
	TFilt  (const TFilt  &R) ;
	// �������� ������������
	TFilt  operator=(TFilt   R2) ;

	TFilt(TShipTarg ShipTarg, TBombTraj BombTraj);

	void OneStepGolubev(long double valFiVisZv,long double valRZv,long double valXZv,long double valTZv,long double sigmAngVis_w
	   ,long double sigmAngVis_bmo 	,long double sigmAngVis_mmo, long double sigmR_w,long double sigmR_bmo ,long double sigmR_mmo
	   , long double sigmX_w,long double sigmX_bmo ,long double sigmX_mmo);

	void 	fncFiltStep(long double valYZv,long double valTZv,long double sigm_w,long double sigm_bmo
,long double sigm_mmo,long double *arrXEst, long double *arrK );

    void OneStepFltrKlm_Real(long double valYZv,long double valTZv,long double sigm_w,long double sigm_bmo
,long double sigm_mmo,long double *arrXEst, long double *arrK );






}  ;
#endif
