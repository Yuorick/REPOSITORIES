//---------------------------------------------------------------------------

#ifndef BombTrajH
#define BombTrajH
#include "Bomb.h"
#include "Filt.h"
#include "CntlFuncPar.h"
// ����� ��������� �������������� �������
class TBomb;
class TFilt ;
class TBombTraj
{
public:
 // ��������� �������
	long double mTStart;  // ��� ������
	long double mTet0;  // ��� ���� ��������
	long double mV0 ;  // ��� ��������
	long double mAltit ; // ��� ������
	long double mX0 ; // ���. ��������� �� ��� X
	long double mTCur ;  // ������ ��������� �������� �������
	long double mStepInt ; // ���� �������������� �� ������� ��� ��������� ��������

	TBomb mBomb ;  // ���
	// �������  ������  � ����������� ��:
// 1. marrStrSK_VS [0]- x
// 2. marrStrSK_VS [1]-  y
// 3. marrStrSK_VS [2]-  ������� ��������� ��������� ��
// 4. marrStrSK_VS [3]-  ������� �������� V
// 5. marrStrSK_VS [4]-  ���� �����

	long double marrStrSK_VS [5] ;

// ������  ���������� � ������ ����� ������ ���������,  marrWSigSq[0] = marrWSigSq[1] = 0.;  0< marrWSigSq[0] < 1
//    DelW[i] = marrCfW[i] * W, i = 2,3,4. �� ����, ��� ����������� ������������ � ����� ��������� ��� ������ ��������� ���������
	 long double marrCfW [5] ;

	// ����������� �� ���������
	TBombTraj () ;
	// ����������� �����������
	TBombTraj  (const TBombTraj  &R) ;

	// �������� ������������
	TBombTraj  operator=(TBombTraj   R2) ;
	// ����� �����������
	TBombTraj (	long double Tet0  // ��� ���� ��������
	,long double V0  // ��� ��������
	,long double Altit  // ��� ������
	,long double X0  // ���. ��������� �� ��� X
	,long double StepInt  // ���� �������������� �� ������� ��� ��������� ��������
	 ,long double TStart
	);

	// ����� ����������� 2
   TBombTraj (
	long double Tet0  // ��� ���� ��������
	,long double V0  // ��� ��������
	,long double Altit  // ��� ������
	 ,long double X0  // ���. ��������� �� ��� X
	,long double StepInt  // ���� �������������� �� ������� ��� ��������� ��������
	 ,long double TStart
	,TBomb Bomb
	);

		 // ����� ����������� 3
   TBombTraj (
	long double Tet0  // ��� ���� ��������
	,long double V0  // ��� ��������
	,long double Altit  // ��� ������
	 ,long double X0  // ���. ��������� �� ��� X
	,long double StepInt  // ���� �������������� �� ������� ��� ��������� ��������
	 ,long double TStart
	,TBomb Bomb
	,long double *arrCfW
	);

	 void fncFillNachalnieUsloviaVS() ;
	 void fncCalc_q(long double valMach,long double  &val_q);
	 void fncCalc_F(const long double valU, long double *arrF);

	 void fncCalcMach(long double valTay, long double valDerivTay ,long double &valMach) ;

	 void fncMovePhasVector(const long double valU, const double valTNext );



	 void fncCalcMach_and_GradMach(long double valTay, long double valDerivTay
	 ,long double &valMach, long double *arrGradMach);

	 void fncCalc_q_and_Grad_q(long double  &val_q, long double *arrGrad_q);

	 void fncCalc_Cx_and_Grad_Cx(long double valMach,long double  *arrGradMach
		,long double  &val_Cx, long double *arrGrad_Cx);

	 void fncCalc_Cy_and_Grad_Cy(long double valMach,long double  *arrGradMach
		,long double  &val_Cy, long double *arrGrad_Cy);

	 void fncCalc_F_and_dF_po_dx_and_dF_po_dU(const long double valU, long double *arrF,long double *arr_dF_po_dx,long double *arr_dF_po_dU);

	 void fncGradF0(long double *arr_dF0_po_dx);
	 void fncGradF1(long double *arr_dF1_po_dx);
	 void fncGradF2(long double valF2,long double valTay, long double valDerivTay , long double *arr_dF2_po_dx)  ;
	 void fncGradF3(long double val_q , long double *arrGrad_q
	, long double val_Cx, long double* arrGrad_Cx, long double *arr_dF3_po_dx);

	void fncGradF4(long double val_q , long double *arrGrad_q
	, long double val_Cy, long double* arrGrad_Cy, long double valU, long double *arr_dF4_po_dx);

	void fncCalc_dF_po_dU(long double val_Cy, long double val_q, long double *arr_dF_po_dU);

	void fncEilerStep(const long double valU, const long double valStepInt);

}  ;
#endif
