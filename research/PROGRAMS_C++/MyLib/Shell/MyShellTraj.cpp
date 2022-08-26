//---------------------------------------------------------------------------


#pragma hdrstop

#include "MyShellTraj.h"

#include <string.h>
#include <math.h>
#include "Atmosphere.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "URPolygon.h"
#include "URPointXY.h"
#include "Environment.h"
#include "Constants.h"



extern	const  double HT[5][3];


	// ����������� ��
  // ������� ���� ��� ������������� �� ��� OX ����������� �� ������ ������� ������� !!!
		 //TShellBody  mShellBodyBody
//---------------------------------------------------------------------------
TMyShellTraj::TMyShellTraj()
{
	mTStart = 0.;  // ��� ������
	mTet0 = 17./180.* M_PI;  // ��� ���� ��������  17 ����


	mAlfDir =  M_PI / 2.; //������������� ���� ��������� ��������
	mPsi0 = 0.; //��� ������ (���� ����) � ����������� �� mPsi0

	mAltit =0.; // ��� ������
	mLatitude = M_PI/3.;  // ������
	mTCur =0. ;  // ������ ��������� �������� �������
  //	StepInt = 0.0001 ;
	mShellBody = TShellBody () ;  // ������
	// ������ ����� �� ������ ���������� �������� �������  - ��� �������� ������� �����������
 //	 double arrDelta [10] = { 10.
 //							   ,10.
 //							   ,10.
 //							   ,0.001
 //							   ,10.
	//						   ,0.01
	//						   ,1.
	//						   ,0.01
	 //						   ,0.001
		//					   ,0.01
		//					   } ;
 //	memcpy(marrDelta, arrDelta, 10 * sizeof( double)) ;

 //	marrDispWindParamsScatters[0] = 0.01;
 //	marrDispWindParamsScatters[2] = 0.01;
	marrDispWindParamsScatters[1] = 0.005 * 0.005;
 //
	// �������  ������  � ����������� ��:
//  marrStrSK_VS [0]- x
//  marrStrSK_VS [1]-  y
//  marrStrSK_VS [2]-  z
//  marrStrSK_VS [3]-  ���� ���
//  marrStrSK_VS [4]-  ������� ��������� �����
//  marrStrSK_VS [5]-  ��������� ��������� ��
//  marrStrSK_VS [6]-  ������� �������� V
//  marrStrSK_VS [7]-  ���� �����
	fncFillNachalnieUsloviaVS();


}

//---------------------------------------------------------------------------


// ����������� �����������
 TMyShellTraj ::TMyShellTraj (const TMyShellTraj &R)
 {
	mTStart = R.mTStart;  // ��� ������
	mTet0 = R.mTet0;  // ��� ���� ��������  17 ����
	mAlfDir = R.mAlfDir; //��� ������ ��������  ������
	mPsi0 = R.mPsi0 ;  //��� ������ (���� ����) � ����������� �� mPsi0

  //	StepInt = R.StepInt ;


	mAltit = R.mAltit; // ��� ������
	mLatitude = R.mLatitude;  // ������
	mTCur = R.mTCur ;  // ������ ��������� �������� �������

	mShellBody = R.mShellBody ;  // ������

	memcpy(marrStrSK_VS, R.marrStrSK_VS, 8 * sizeof( double));
 //	memcpy(marrStrSK_Jac, R.marrStrSK_Jac, 64 * sizeof( double));
 //	memcpy(marrStrSK_JacCoeffFom, R.marrStrSK_JacCoeffFom, 8 * sizeof( double));
 //	memcpy(marrStrSK_JacMass, R.marrStrSK_JacMass, 8 * sizeof( double));
//	memcpy(marrDelta, R.marrDelta, 10 * sizeof( double));
//	memcpy(marrDispWindParamsScatters, R.marrDispWindParamsScatters , 3 * sizeof(double));

 }

 // �������� ������������
 TMyShellTraj &TMyShellTraj::operator=(const TMyShellTraj  &R)
 {
	mTStart = R.mTStart;  // ��� ������
	mTet0 = R.mTet0;  // ��� ���� ��������  17 ����
	mAlfDir = R.mAlfDir; //��� ������ ��������  ������
	mPsi0 = R.mPsi0 ;  //��� ������ (���� ����) � ����������� �� mPsi0




	mAltit = R.mAltit; // ��� ������
	mLatitude = R.mLatitude;  // ������
	mTCur = R.mTCur ;  // ������ ��������� �������� �������

	mShellBody = R.mShellBody ;  // ������

	memcpy(marrStrSK_VS, R.marrStrSK_VS, 8 * sizeof( double));
 //	memcpy(marrStrSK_Jac, R.marrStrSK_Jac, 64 * sizeof( double));
 //	memcpy(marrStrSK_JacCoeffFom, R.marrStrSK_JacCoeffFom, 8 * sizeof( double));
 //	memcpy(marrStrSK_JacMass, R.marrStrSK_JacMass, 8 * sizeof( double));
 //	memcpy(marrDelta, R.marrDelta, 10 * sizeof( double));
 //	memcpy(marrDispWindParamsScatters, R.marrDispWindParamsScatters , 3 * sizeof(double));

	return *this ;
 }


  // ����� �����������
 TMyShellTraj::TMyShellTraj (const  double TStart, const  double Tet0
			  , const  double AlfDir
			  ,const  double Psi0 ,const  double  Altit,const  double Latitude
			  ,const  double  TCur)
 {
	 mTStart = TStart;  // ��� ������
	 mTet0 = Tet0;  // ��� ���� ��������
   //	 double mBet0 ; //��� ������ ��������
	 mAlfDir = AlfDir ; //������������� ���� ��������� ��������
	 mPsi0 = Psi0 ; // ��� ������ (���� ����)
 
	 mAltit = Altit ; // ��� ������
	 mLatitude = Latitude ;  // ������
	 mTCur = TCur ;  // ������ ��������� �������� �������

	 mShellBody = TShellBody() ;  // ������


	 fncFillNachalnieUsloviaVS();

 }

  // ����� ����������� 2
 TMyShellTraj::TMyShellTraj (const TShellBody ShellBody, const  double TStart, const  double Tet0
			 , const  double AlfDir
			  ,const  double Psi0 ,const  double  Altit,const  double Latitude
			  ,const  double  TCur)
 {
	 mTStart = TStart;  // ��� ������
	 mTet0 = Tet0;  // ��� ���� ��������
   //	 double mBet0 ; //��� ������ ��������
	 mAlfDir = AlfDir ; //������������� ���� ��������� ��������
	 mPsi0 = Psi0 ; // ��� ������ (���� ����)




	 mAltit = Altit ; // ��� ������
	 mLatitude = Latitude ;  // ������
	 mTCur = TCur ;  // ������ ��������� �������� �������

	 mShellBody = ShellBody;  // ������

	 fncFillNachalnieUsloviaVS();

 }
 ///


  // ����� ����������� 3
 TMyShellTraj::TMyShellTraj (const TShellBody ShellBody, const  double TStart, const  double Tet0
			 , const  double AlfDir,const  double Psi0 ,const  double  Altit,const  double Latitude
			  ,const  double  TCur, double *arrDispWindParamsScatters)
 {
	 mTStart = TStart;  // ��� ������
	 mTet0 = Tet0;  // ��� ���� ��������
   //	 double mBet0 ; //��� ������ ��������
	 mAlfDir = AlfDir ; //������������� ���� ��������� ��������
	 mPsi0 = Psi0 ; // ��� ������ (���� ����)




	 mAltit = Altit ; // ��� ������
	 mLatitude = Latitude ;  // ������
	 mTCur = TCur ;  // ������ ��������� �������� �������

	 mShellBody = ShellBody;  // ������
	 // ������ ����� �� ������ ���������� �������� �������  - ��� �������� ������� �����������
	// double arrDelta [10] = { 10.
	//						   ,10.
	////						   ,10.
	 //						   ,0.001
	//						   ,10.
	 //						   ,0.01
	 //						   ,1.
	 //						   ,0.01
		//					   ,0.001
		 //					   ,0.01
		//					   } ;
 //	memcpy(marrDelta, arrDelta, 10 * sizeof( double)) ;
 //	memcpy(marrDispWindParamsScatters, arrDispWindParamsScatters, 3 * sizeof( double)) ;

	 fncFillNachalnieUsloviaVS();

 }
 ///


  // ����� ����������� 4
  // Eps0 - ���� ����� �������� � ����
  // Bet0-  ����. ����  �������� � ����
  // arrVesselVElocity [3] - �������� ������� (���)
  // ShellBody -  �������
 TMyShellTraj::TMyShellTraj (double *arrVesselVelocity
 , const TShellBody ShellBody,  const  double Eps0,  const  double Bet0 )


 {
	 mTStart = 0.;  // ��� ������
	 mTCur = 0.;
	 // ������ ��������� �������� � ���
	  double arrGSKV0[3] = {0.};
	  arrGSKV0[0] =  arrVesselVelocity[0] + ShellBody.mV0 * cos(Eps0) * sin( Bet0);
	  arrGSKV0[1] =  arrVesselVelocity[1] + ShellBody.mV0 * cos(Eps0) * cos( Bet0);
	  arrGSKV0[2] =  arrVesselVelocity[2] + ShellBody.mV0 * sin(Eps0) ;
	 ///
	 //
	 mShellBody = ShellBody;  // ������
	 mShellBody.mV0 =  Norm3( arrGSKV0) ;

	 ///

	 mTet0 = asin( arrGSKV0[2]/ mShellBody.mV0);
	 mAlfDir = atan2( arrGSKV0[0], arrGSKV0[1]);
	 mPsi0 = 0. ; // ��� ������ (���� ����)

	 mAltit = 0. ; // ��� ������
	 mLatitude = M_PI /3. ;  // ������

	 fncFillNachalnieUsloviaVS();
 }
//



void TMyShellTraj::fncFillNachalnieUsloviaVS()
{
  marrStrSK_VS[0] = 0.;
  marrStrSK_VS[1] = mAltit ;
  marrStrSK_VS[2] = 0.;
  marrStrSK_VS[3] =  mPsi0 ;
  marrStrSK_VS[4] = mShellBody.mOmega0 ;
 // marrStrSK_VS[5] = ATM_PN0 ;//
  marrStrSK_VS[5] = 1 ;//
  marrStrSK_VS[6] = mShellBody.mV0 ;
  marrStrSK_VS[7] = mTet0 ;
 // memset(marrStrSK_Jac, 0, 64 * sizeof( double));

 // memset(marrStrSK_JacCoeffFom, 0, 8 * sizeof( double));
//  memset(marrStrSK_JacMass, 0, 8 * sizeof( double));

	 mCoefCx = mShellBody.mplnCix.LinearValueApprox(mTet0 / M_PI * 180.) ;
	 mCoefCy = mShellBody.mplnCiy.LinearValueApprox(mTet0 / M_PI * 180.) ;
	 mCoefCz = mShellBody.mplnCiz.LinearValueApprox(mTet0 / M_PI * 180.) ;


}



// ���������� ������ ������� ������ ����� ������� ��� ��������� -arrF
// ������� ������� ����������� ������ ����� �� x  - mtr_dF_po_dx
// ������� ������� ����������� ������ ����� �� ��������� ��������   - mtr_dF_po_dz
void TMyShellTraj::fncCalc_F_and_H_and_HI( double *arrF,  double *mtr_dF_po_dx
	,  double *mtr_dF_po_dz, double *mtr_dF_po_di)
{
	/* memset ( mtr_dF_po_dz, 0, 64 * sizeof( double)) ;
   memset ( mtr_dF_po_di, 0, 16 * sizeof( double)) ;
  // ������� ������
   const  double valSm = M_PI * mShellBody.mDm * mShellBody.mDm / 4. ;
  // ���������� ���������� ����������� �����������  � �� ���������
   double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///

  // ���������� ����� ���� � ������� ��� ������� �����������
   double valMach = 0.,arrGradMach[8] ={0.} ;
  fncCalcMach_and_GradMach(valTay,valDerivTay, valMach, arrGradMach);
  ///

  // ���������� ������ q � ������� �� ���������
   double val_q = 0.,arrGrad_q[8] ={0.} ;
  fncCalc_q_and_Grad_q(valMach,arrGradMach, val_q, arrGrad_q);
  ///

  // ������� ����������� ��������
   double valVpr = marrStrSK_VS[6] * cos(marrStrSK_VS [7]) * (R_ZEMLI + mAltit)/
	   (R_ZEMLI + marrStrSK_VS [1])  ;
	   ///

  // ���������� �������������� ������� Knm � ������� �� ����������� �� x
  // �� ������� 1-�� � 6-�� ���������� ( �������� ���������� � 0)
	 double valKnm = 0., arrGradKnm[8] = {0.} ;
	fncCalcKnm_and_Grad_Knm( valMach, arrGradMach
			   ,valKnm, arrGradKnm ) ;
	///

 // ���������� ��������������� ������ iz(z2) � �� �����������
	 double val_iz = 0., val_Deriv_iz = 0. ;
	mShellBody.fnkIz0(mTet0,val_iz, val_Deriv_iz) ;
	///

// ���������� ��������������� ������ ix(z2) � �� �����������
	 double val_ix = 0., val_Deriv_ix = 0. ;
	mShellBody.fnkIx0(mTet0,val_ix, val_Deriv_ix) ;
	///

// ���������� �������������� ������� mxwx  � �� ����������� �� M
	 double valMxOmegax = 0., arrGradMxOmegax[8] = {0.} ;
	fncCalcMxOmegax_and_Grad_MxOmegax( valMach, arrGradMach
			   ,valMxOmegax, arrGradMxOmegax ) ;
	///

// ���������� �������������� ������� CxEtal  � �� ����������� �� M
	 double valCxEtal = 0., arrGradCxEtal [8] = {0.} ;
	fncCalcCxEtal_and_Grad_CxEtal(valMach, arrGradMach
			   ,valCxEtal, arrGradCxEtal ) ;
	///

// ���������� ���������� ���������� ������������ ������ ����� ��� ����� �������� �����
    double valDeltaTettaTochka =0.          //�������  �������� ���� ������� ������ DeltaTettaTochka
			  ,valDerivDeltaTettaTochkaPoPsi= 0. // ����������� DeltaTettaTochka �� Psi (= marrStrSK_VS [3])
			  ,valDeltaPsiTochka =0.           //�������  �������� ���� ���� DeltaPsiTochka
			  ,valDerivDeltaPsiTochkaPoPsi= 0. // ����������� DeltaPsiTochka �� Psi (= marrStrSK_VS [3])
			  ,valDerivDeltaPsiTochkaPoTetta= 0.;  // ����������� DeltaPsiTochka �� Tetta (=  marrStrSK_VS [7])
   fncCalcDeltaTettaTochka_and_DerivPoPsi(valDeltaTettaTochka, valDerivDeltaTettaTochkaPoPsi) ;
   fncCalcDeltaPsiTochka_and_DerivPoPsi_and_DerivPoTetta(valDeltaPsiTochka
		   ,valDerivDeltaPsiTochkaPoPsi, valDerivDeltaPsiTochkaPoTetta);
	 ///

   arrF[0] =  valVpr * cos(marrStrSK_VS [3]) ;
   arrF[1] =  marrStrSK_VS [6] * sin(marrStrSK_VS [7]) ;
   arrF[2] =  -valVpr * sin(marrStrSK_VS [3]) ;

  /// mtr_dF_po_di[3 * 2 + 1] = mShellBody.mvalIx0 * mShellBody.mL * val_iz *valKnm * marrStrSK_VS [4] *
   //		(-G_ZEMLI/ marrStrSK_VS [6]/ marrStrSK_VS [6] + 1./(R_ZEMLI + marrStrSK_VS [1]))/
   //		(mShellBody.mMass * mShellBody.mDm * mShellBody.mh_gob);
  // arrF[3] =  mtr_dF_po_di[3 * 2 + 1] + valDeltaPsiTochka;
  // arrF[3] =  (mShellBody.mvalIx0 * mShellBody.mL / (mShellBody.mMass * mShellBody.mDm * mShellBody.mh_gob))
	//	  * val_iz *valKnm * marrStrSK_VS [4]
	 //	  *(-G_ZEMLI/ marrStrSK_VS [6]/ marrStrSK_VS [6] + 1./(R_ZEMLI + marrStrSK_VS [1]))+ valDeltaPsiTochka;
  // mtr_dF_po_di[3 * 2 + 1] = ( arrF[3] - valDeltaPsiTochka)/ val_iz ;
   // ��������� ���������� ��� ��������� �����
  // const  double val_k0 = -valSm * mShellBody.mL * mShellBody.mL/ mShellBody.mvalIx0 / ATM_AN0 / 2.;
  // arrF[4] = val_k0 *valMxOmegax *val_q * marrStrSK_VS [4] / valMach ;
	 const  double val_qsbm = 0.474 * mShellBody.mDm* mShellBody.mDm/ mShellBody.mMass * marrStrSK_VS [6] * marrStrSK_VS [6] * marrStrSK_VS[5];
	 const  double val_mmx = (valMxOmegax * mShellBody.mDm /2./ mShellBody.mV0);
	 const  double val_mn2  = val_mmx* val_qsbm  * mShellBody.mMass * mShellBody.mL  / mShellBody.mvalIx0 ;
	arrF[4] = -val_mn2 * marrStrSK_VS [4] ;
   ///

   //arrF[5] = -G_ZEMLI *  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) / valTay /ATM_R_UNIVER;
   arrF[5] = -  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) / valTay*(G_ZEMLI /ATM_R_UNIVER + valDerivTay);//HT[0][2]);

  mtr_dF_po_di[6 * 2 ] =  -  valCxEtal * val_q * valSm / mShellBody.mMass ;

    double valTemp = val_ix * valCxEtal * val_q * marrStrSK_VS [6]*marrStrSK_VS [6] ;
   arrF[6] = -G_ZEMLI * sin(marrStrSK_VS [7])   -  val_ix * valCxEtal * val_q * marrStrSK_VS [6]*marrStrSK_VS [6] ;
   arrF[7] = -G_ZEMLI * cos(marrStrSK_VS [7])/ marrStrSK_VS [6]  +
			marrStrSK_VS [6] * cos(marrStrSK_VS [7]) / (R_ZEMLI + marrStrSK_VS [1]) +
			valDeltaTettaTochka;

	arrF[3] =  (mShellBody.mvalIx0 * mShellBody.mL / (mShellBody.mMass * mShellBody.mDm * mShellBody.mh_gob))
		  * val_iz *valKnm * marrStrSK_VS [4]* arrF[7]/marrStrSK_VS [6]/ cos(marrStrSK_VS [7])+ valDeltaPsiTochka; ;
		//  *(-G_ZEMLI/ marrStrSK_VS [6]/ marrStrSK_VS [6] + 1./(R_ZEMLI + marrStrSK_VS [1]))+ valDeltaPsiTochka;
  // mtr_dF_po_di[3 * 2 + 1] = ( arrF[3] - valDeltaPsiTochka)/ val_iz ;

 */

   /*
// �������  ������� ������� ����������� ������ ����� �� ��������� ��������   - mtr_dF_po_dz
   memset(mtr_dF_po_dx, 0, 64 * sizeof( double))  ;

    double arrGradF3_po_z [8] = {0} , arrGradF6_po_z [8] = {0.} ;
   // arr_dFTransp - ��������������� ������ ��� �������� ���������� ������� Fi. ��������� �������� ���������.
   // �������, ��� ��������� �������  mtr_dF_po_dx  ������ mtr_dF_po_dx ���� ���������������
    double	arr_dFTransp[64]= {0.} ;
   fncCalcGradF0(  arr_dFTransp    );       // ������ ��������� ������� F0
   fncCalcGradF1( &arr_dFTransp[ 8]); // ������ ��������� ������� F1
   fncCalcGradF2( &arr_dFTransp[16]); // ������ ��������� ������� F2
   fncCalcGradF7( valDerivDeltaTettaTochkaPoPsi ,&arr_dFTransp[7 * 8]) ;// ������ ��������� ������� F5

   fncCalcGradF3(&arr_dFTransp[7 * 8],arrF[7], valDeltaPsiTochka      // ������ ��������� ������� F3
		   ,valDerivDeltaPsiTochkaPoPsi, valDerivDeltaPsiTochkaPoTetta
		   ,val_iz, val_Deriv_iz
		   ,valKnm, arrGradKnm
		   ,&arr_dFTransp[ 3 * 8],arrGradF3_po_z) ;

   fncCalcGradF4(valMach, arrGradMach    // ������ ��������� ������� F4
				 ,val_q, arrGrad_q
				 ,valMxOmegax, arrGradMxOmegax
				 ,&arr_dFTransp[ 4 * 8]) ;

   fncCalcGradF5( valTay, valDerivTay, &arr_dFTransp[ 5 * 8]) ;



   fncCalcGradF6(val_q, arrGrad_q
							   ,valCxEtal, arrGradCxEtal
							   , val_ix,  val_Deriv_ix
							   ,&arr_dFTransp[ 6 * 8],arrGradF6_po_z); // ������ ��������� ������� F6
  ///


 //MatrTransp(arr_dFTransp, 8, 8, mtr_dF_po_dx);
 memcpy( mtr_dF_po_dx, arr_dFTransp, 64 * sizeof(double)) ;
 // ��� ������� !!! �� ������
 mtr_dF_po_dz[ 3 * 8 +7 ] = arrGradF3_po_z[7] ;
 mtr_dF_po_dz[ 6 * 8 +7 ] = arrGradF6_po_z[7] ;
  */


}

void TMyShellTraj::fncCalc_F(TEnvironment Environment, double *arrF)
{
  // ���������� ���������� ����������� �����������  � �� ���������
   double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///
   // ���������� ������� �������� ����� � ���
   double arrWindV_SSK[3] ={0.};
   calcVectWindV(Environment, arrWindV_SSK);
   ///

   // ���������� ��������� ����������
   //const double VAlVozdV = sqrt(marrStrSK_VS[6] * marrStrSK_VS[6] - 2. * arrWindV_SSK[0] * marrStrSK_VS[6]
  //	+ arrWindV_SSK[1]  * arrWindV_SSK[1]  + Environment.mWind_VertV * Environment.mWind_VertV);

	 const double VAlVozdV = sqrt(marrStrSK_VS[6] * marrStrSK_VS[6] - 2. * arrWindV_SSK[0] * marrStrSK_VS[6]
	+ arrWindV_SSK[1]  * arrWindV_SSK[1]  + arrWindV_SSK[2]  * arrWindV_SSK[2]);
	///

  // ���������� ����� ���� � ������� ��� ������� �����������
   double valMach = calcMach(  valTay,  VAlVozdV);

  ///


	///
	double valCx = -1., val_Deriv_Cx = 0. ;
	mShellBody.fnkCxEtal (valMach,  valCx, val_Deriv_Cx) ;


  /*
  // ��� ��� � �����
  // ���������� ������������� deltaCx,deltaCy, deltaCz
  double valSinEps2 = - arrWindV_SSK[1] /VAlVozdV;
  double valCosEps2 = sqrt(1. - valSinEps2 * valSinEps2);
 // double valSinEps1 = arrWindV_SSK[2] /VAlVozdV /valCosEps2;
  double valSinEps1 = arrWindV_SSK[2] /VAlVozdV *valCosEps2;
  double valCosEps1 = sqrt(1. - valSinEps1 * valSinEps1);


  double deltaCx  = valCx * (valCosEps1 * valCosEps2 -1. );
  double deltaCy =  valCx * valSinEps2 ;
 // double deltaCz =  valCx * valSinEps1 * valCosEps2;
 double valCz =   mLearnShellBody.fnkCz (valMach);
 double deltaCz =  valCz * valSinEps1 * valCosEps2;  */
  ///
   ///
 //

 // ���������� ������� �������� ����� � ���
   double arrWindV_SSK1[3] ={0.},arrWindV_SSK2[3] ={0.}, arrWindV_GSK[3] = {0.};
   arrWindV_GSK[0] =  -Environment.mWind_V * sin(Environment.mWind_Alf);
   arrWindV_GSK[1] =  -Environment.mWind_V * cos(Environment.mWind_Alf);
	arrWindV_GSK[2] =  -Environment.mWind_VertV;
   calcVectWindV(Environment, arrWindV_SSK1);
   transform_xyzGSK_To_xyzSSK( 3, arrWindV_GSK, arrWindV_SSK2);

   double arrAirV_SSK[3] = {0.};
   arrAirV_SSK[0] = marrStrSK_VS[6] -  arrWindV_SSK2[0];
   arrAirV_SSK[1] = -  arrWindV_SSK2[1];
   arrAirV_SSK[2] = -  arrWindV_SSK2[2];

   double temp00 = arrAirV_SSK[2]/ Norm3( arrAirV_SSK);
   if (fabs(temp00) > 0.99999999)
   {
	 temp00 =  0.99999999 * SIGNUM( temp00);
   }
   double epsW = asin( temp00 );

   temp00 = arrAirV_SSK[1]/ sqrt(arrAirV_SSK[1] *  arrAirV_SSK[1] + arrAirV_SSK[0] *  arrAirV_SSK[0] );
   if (fabs(temp00) > 0.99999999)
   {
	 temp00 =  0.99999999 * SIGNUM( temp00);
   }
   double betW = asin( temp00);


   temp00 = arrAirV_SSK[2]/ sqrt(arrAirV_SSK[2] *  arrAirV_SSK[2] + arrAirV_SSK[0] *  arrAirV_SSK[0] );
   if (fabs(temp00) > 0.99999999)
   {
	 temp00 =  0.99999999 * SIGNUM( temp00);
   }
   double gamW = asin(temp00);


	double deltaCx  = valCx * (cos(epsW) * cos(betW) -1. );
	double deltaCy =  valCx * sin(betW) ;
 // double deltaCz =  valCx * sin(epsW);
	double valCz =   mShellBody.fnkCz (valMach);
 double deltaCz =  valCz *   sin(gamW);
 ////////

  // ���������� ������ q

	 double val_q = calc_q(VAlVozdV, valTay) ;
  ///

  // ������� ����������� ��������
   double valVpr = marrStrSK_VS[6] * cos(marrStrSK_VS [7]) * (R_ZEMLI + mAltit)/
	   (R_ZEMLI + marrStrSK_VS [1])  ;
	   ///

  // ���������� �������������� ������� Knm � ������� �� ����������� �� x
  // �� ������� 1-�� � 6-�� ���������� ( �������� ���������� � 0)
   //	 double valKnm = 0., valDerivKnm = 0. ;

	//  mShellBody.fnkKnm(valMach, valKnm, valDerivKnm);
	///

 // ���������� ��������������� ������ iz(z2) � �� �����������
   //	 double val_iz = 0., val_Deriv_iz = 0. ;
 //	mShellBody.fnkIz0(mTet0,val_iz, val_Deriv_iz) ;
	///
	// ���������� ��������������� ������ iz(z2) � �� �����������
	 double valKnm = 0., val_Deriv_Knm = 0. ;
	 mShellBody.fnkKnm (valMach, valKnm, val_Deriv_Knm )  ;


// ���������� �������������� ������� mxwx  � �� ����������� �� M

	double valMxOmegax = 0., valDerivMxOmegax = 0;
	 mShellBody.fnkMxOmegax(valMach,valMxOmegax, valDerivMxOmegax ) ;
	///

  // ������� ������
	const double VAlSmid = mShellBody.mDm * mShellBody.mDm / 4. * M_PI;

// ���������� ���������� ���������� ������������ ������ ����� ��� ����� �������� �����
    double valDeltaTettaTochka =0.          //�������  �������� ���� ������� ������ DeltaTettaTochka
			  ,valDerivDeltaTettaTochkaPoPsi= 0. // ����������� DeltaTettaTochka �� Psi (= marrStrSK_VS [3])
			  ,valDeltaPsiTochka =0.           //�������  �������� ���� ���� DeltaPsiTochka
			  ,valDerivDeltaPsiTochkaPoPsi= 0. // ����������� DeltaPsiTochka �� Psi (= marrStrSK_VS [3])
			  ,valDerivDeltaPsiTochkaPoTetta= 0.;  // ����������� DeltaPsiTochka �� Tetta (=  marrStrSK_VS [7])
   fncCalcDeltaTettaTochka_and_DerivPoPsi(valDeltaTettaTochka, valDerivDeltaTettaTochkaPoPsi) ;
   fncCalcDeltaPsiTochka_and_DerivPoPsi_and_DerivPoTetta(valDeltaPsiTochka
		   ,valDerivDeltaPsiTochkaPoPsi, valDerivDeltaPsiTochkaPoTetta);
	 ///

   arrF[0] =  valVpr * cos(marrStrSK_VS [3]) ;
   arrF[1] =  marrStrSK_VS [6] * sin(marrStrSK_VS [7]) ;
   arrF[2] =  -valVpr * sin(marrStrSK_VS [3]) ;

	// const  double val_qsbm = 0.474 * mShellBody.mDm* mShellBody.mDm/ mShellBody.mMass * marrStrSK_VS [6] * marrStrSK_VS [6] * marrStrSK_VS[5];
 //	 const  double val_mmx = (valMxOmegax * mShellBody.mDm /2./ mShellBody.mV0);
	// const  double val_mn2  = val_mmx* val_qsbm  * mShellBody.mMass * mShellBody.mL  / mShellBody.mvalIx0 ;
	arrF[4] = -valMxOmegax *mShellBody.mL *mShellBody.mL/ mShellBody.mvalIx0 * val_q
			 * VAlSmid /valMach /ATM_AN0 * marrStrSK_VS [4];
	// -val_mn2 * marrStrSK_VS [4] ;
	 ///

   arrF[5] = -  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) / valTay*(G_ZEMLI /ATM_R_UNIVER + valDerivTay);//HT[0][2]);



  //  double valTemp = val_ix * valCxEtal * val_q * marrStrSK_VS [6]*marrStrSK_VS [6] ;
   arrF[6] = -G_ZEMLI * sin(marrStrSK_VS [7])   -   (valCx + deltaCx) * VAlSmid/ mShellBody.mMass * val_q ;
   arrF[7] = -G_ZEMLI * cos(marrStrSK_VS [7])/ marrStrSK_VS [6]  +
			marrStrSK_VS [6] * cos(marrStrSK_VS [7]) / (R_ZEMLI + marrStrSK_VS [1])
			-  deltaCy * VAlSmid/ mShellBody.mMass * val_q / marrStrSK_VS [6]
			+ valDeltaTettaTochka;

	arrF[3] =  (mShellBody.mvalIx0 * mShellBody.mL / (mShellBody.mMass * mShellBody.mDm * mShellBody.mh_gob))
		  * valKnm * marrStrSK_VS [4]* arrF[7]/marrStrSK_VS [6]/ cos(marrStrSK_VS [7])+ valDeltaPsiTochka
		  - deltaCz * VAlSmid/ mShellBody.mMass * val_q / cos(marrStrSK_VS [7])/marrStrSK_VS [6] ;




}




// ���������� ������� �������� ����� � ���

void TMyShellTraj::calcVectWindV(TEnvironment Environment, double *arrWindV_SSK)
{
 arrWindV_SSK[0] = -Environment.mWind_V * cos(Environment.mWind_Alf - (mAlfDir - marrStrSK_VS [3]))
	 * cos(marrStrSK_VS [7] ) + Environment.mWind_VertV * sin (marrStrSK_VS [7]);
 arrWindV_SSK[1] =  Environment.mWind_V * cos(Environment.mWind_Alf - (mAlfDir - marrStrSK_VS [3]))
	 * sin(marrStrSK_VS [7] ) + Environment.mWind_VertV * cos (marrStrSK_VS [7]);
 arrWindV_SSK[2] = -Environment.mWind_V * sin(Environment.mWind_Alf - (mAlfDir - marrStrSK_VS [3])) ;
}

//-=-------------------------------------------------------
double TMyShellTraj::SIGNUM( double temp00)
{

	if(temp00 > 0)
	{
	return 1.;

	}
	else
	{
	return -1.;
	}
}
// ���������� ����� ���� � ������� ������ ����������� ����� ���� �� x
// INPUT:
//valTay - ���������� �����������  �����������
// valDerivTay - ����������� �� ������ ���������� �����������  �����������
// OUTPUT
// valMach - ����� ����
// arrGradMach - ������ ��������� ����� ���� �� ������� ����������
// ���� ���������� !!!!!
void TMyShellTraj::fncCalcMach_and_GradMach(const double valTay, const double valDerivTay
	 , double &valMach,  double *arrGradMach)
{
  /* memset( arrGradMach, 0, 8 * sizeof( double));
   arrGradMach [6] = sqrt(ATM_TAYN0/ valTay)/ ATM_AN0 ;
   valMach = marrStrSK_VS[6] * arrGradMach [6] ;
   arrGradMach [1] = - 0.5 * valMach * valDerivTay/ valTay ; */

}

// ���������� ����� ���� � ������� ������ ����������� ����� ���� �� x
// INPUT:
//valTay - ���������� �����������  �����������
// valVVozd - ��������� ��������

double TMyShellTraj::calcMach( double valTay,  double valVVozd)
{
 return valVVozd *sqrt(ATM_TAYN0/ valTay)/ ATM_AN0 ;

}

// ���������� Knm(M(x7,x2))�  ������� ������� ����������� ������� Knm(M(x7,x2)) �� x
// INPUT:
//valMach - ����� ����
// arrGradMach - �������� ����� ���� �� x
// OUTPUT
// valKnm - ������� Knm
// arrGradKnm - ������ ��������� ������� Knm �� ������� ����������
void TMyShellTraj::fncCalcKnm_and_Grad_Knm(const  double valMach,  double *arrGradMach
			   , double &valKnm,  double *arrGradKnm )
{
  /* memset( arrGradKnm, 0, 8 * sizeof( double));
	double  valDerivKnm = 0;
   mShellBody.fnkKnm(valMach, valKnm, valDerivKnm);
   arrGradKnm [1] = valDerivKnm * arrGradMach [1] ;
   arrGradKnm [6] = valDerivKnm * arrGradMach [6] ;
  */
}




 double valMxOmegax = 0., valDerivMxOmegax = 0. ;

// ���������� MxOmegax(M(x7,x2))�  ������� ������� ����������� ������� MxOmegax(M(x7,x2)) �� x
// INPUT:
//valMach - ����� ����
// arrGradMach - �������� ����� ���� �� x
// OUTPUT
// valMxOmegax - ������� valMxOmegax
// arrGradMxOmegax - ������ ��������� ������� MxOmegax �� ������� ����������
void TMyShellTraj::fncCalcMxOmegax_and_Grad_MxOmegax(const  double valMach,  double *arrGradMach
			   , double &valMxOmegax,  double *arrGradMxOmegax )
{
  /* memset( arrGradMxOmegax, 0, 8 * sizeof( double));
	double  valDerivMxOmegax = 0;
   mShellBody.fnkMxOmegax(valMach,valMxOmegax, valDerivMxOmegax ) ;
   arrGradMxOmegax [1] = valDerivMxOmegax * arrGradMach [1] ;
   arrGradMxOmegax [6] = valDerivMxOmegax * arrGradMach [6] ;
   */

}
void TMyShellTraj::fncCalcCxEtal_and_Grad_CxEtal(const  double valMach,  double *arrGradMach
			   , double &valCxEtal,  double *arrGradCxEtal )
{
   memset( arrGradCxEtal, 0, 8 * sizeof( double));
	double  valDerivCxEtal = 0;
   mShellBody.fnkCxEtal(valMach,valCxEtal, valDerivCxEtal ) ;
   arrGradCxEtal [1] = valDerivCxEtal * arrGradMach [1] ;
   arrGradCxEtal [6] = valDerivCxEtal * arrGradMach [6] ;

}



// ����������� ����������� ������ � ������� ������� ����������� �������� ������ �� x
// INPUT:
// valMach - ����� ����
// arrGradMach[8] - ������ ��������� ����� ���� �� ������� ����������
// OUTPUT :
// val_q - ���������� �����
// arrGrad_q[8]  -  ������ ��������� c���������� ������  �� ������� ����������
// ���������� !!!!!
void TMyShellTraj::fncCalc_q_and_Grad_q( double valMach, double  *arrGradMach
		, double  &val_q,  double *arrGrad_q)
{

/*   double temp = ATM_RoN0 * ATM_AN0 *ATM_AN0 * valMach;
 //val_q =  temp * marrStrSK_VS [5] * valMach / 2. ;
  memset( arrGrad_q, 0, 8 * sizeof( double));
  arrGrad_q [1] = temp *marrStrSK_VS [5] *arrGradMach [1] ;
  arrGrad_q [5] = temp * valMach / 2. ;
  arrGrad_q [6] =  temp *marrStrSK_VS [5] *arrGradMach [6] ;
 val_q =  mShellBody.mDm* mShellBody.mDm * 0.474 / mShellBody.mMass * marrStrSK_VS [5] ;
*/
}



// ����������� ����������� ������ � ������� ������� ����������� �������� ������ �� x
// INPUT:
// valMach - ����� ����
// arrGradMach[8] - ������ ��������� ����� ���� �� ������� ����������
// OUTPUT :
// val_q - ���������� �����
// arrGrad_q[8]  -  ������ ��������� c���������� ������  �� ������� ����������
// ���������� !!!!!
double TMyShellTraj::calc_q(const  double VAlVVozd, const  double VAlTay)
{
return ATM_RoN0 / 2. * ATM_TAYN0 /VAlTay * marrStrSK_VS [5] *  VAlVVozd * VAlVVozd;
}

// ������ ���������� �� ��������� ����� ��� ����� �������� �������� �����
 void TMyShellTraj::fncCalcDeltaTettaTochka_and_DerivPoPsi( double &valDeltaTettaTochka
		   ,  double &valDerivDeltaTettaTochkaPoPsi)
 {
	  double temp = 2.* OMEGA_ZEMLI * cos (mLatitude);
	 valDeltaTettaTochka  =  temp * sin(mAlfDir - marrStrSK_VS [3]) ;
	 valDerivDeltaTettaTochkaPoPsi = - temp * cos(mAlfDir - marrStrSK_VS [3]) ;
	   //��� ������� !!!
   //	  valDeltaTettaTochka = 0. ;
   //	 valDerivDeltaTettaTochkaPoPsi = 0. ;
 }

 // ������ ���������� �� ��������� ����� ��� ����� �������� �������� �����
 void TMyShellTraj::fncCalcDeltaPsiTochka_and_DerivPoPsi_and_DerivPoTetta( double &valDeltaPsiTochka
		   ,  double &valDerivDeltaPsiTochkaPoPsi,  double &valDerivDeltaPsiTochkaPoTetta)
 {
	 valDeltaPsiTochka =  -2.* OMEGA_ZEMLI *
	 ( sin (mLatitude) - cos (mLatitude) * cos(mAlfDir - marrStrSK_VS [3])* tan(marrStrSK_VS [7])) ;
	 valDerivDeltaPsiTochkaPoPsi = 2.* OMEGA_ZEMLI * cos (mLatitude)*sin(mAlfDir - marrStrSK_VS [3])*tan(marrStrSK_VS [7]) ;
	 valDerivDeltaPsiTochkaPoTetta = 2.* OMEGA_ZEMLI * cos (mLatitude)*cos(mAlfDir - marrStrSK_VS [3])
	   /cos(marrStrSK_VS [7])/cos(marrStrSK_VS [7]) ;
	  //��� ������� !!!
	 // valDeltaPsiTochka = 0. ;
	 //valDerivDeltaPsiTochkaPoPsi = 0. ;
	// valDerivDeltaPsiTochkaPoTetta = 0. ;
 }

 // ���������� ������� ������� ����������� (���������) ������� F0
 // OUTPUT:
 // arrGradF0 [8] - ������ �����������
 void TMyShellTraj:: fncCalcGradF0(  double *arrGradF0)
 {
	memset( arrGradF0,0, 8 * sizeof( double)) ;
	 double valTemp0 = (mAltit + R_ZEMLI) /(marrStrSK_VS [1] + R_ZEMLI); // ���������� ����������
	 double valTemp1 = valTemp0 * cos(marrStrSK_VS [7]); // ���������� ����������
	 double valTemp2 = cos(marrStrSK_VS [3] ); // ���������� ����������
	arrGradF0[1] = - marrStrSK_VS [6]* valTemp1 * valTemp2 /(marrStrSK_VS [1] + R_ZEMLI);
	arrGradF0[3] =  - sin( marrStrSK_VS [3] )* marrStrSK_VS [6] *valTemp1 ;
	arrGradF0[6] = valTemp2 * valTemp1 ;
	arrGradF0[7] = - marrStrSK_VS [6]* sin(marrStrSK_VS [7])* valTemp0 * valTemp2;
 }

  // ���������� ������� ������� ����������� (���������) ������� F1
 // OUTPUT:
 // arrGradF1 [8] - ������ �����������
 void TMyShellTraj::fncCalcGradF1(  double *arrGradF1)
 {
	memset( arrGradF1,0, 8 * sizeof( double)) ;
	arrGradF1[6] = sin(marrStrSK_VS [7]) ;
	arrGradF1[7] =  marrStrSK_VS [6]* cos(marrStrSK_VS [7]);
 }

 // ���������� ������� ������� ����������� (���������) ������� F2
 // INPUT:
 // valDerivDeltaTettaTochkaPoPsi - ����������� ����������� �����(�������� �����) �� ���
 //
 // OUTPUT:
 // arrGradF2 [8] - ������ �����������
 void TMyShellTraj::fncCalcGradF2(  double *arrGradF2)
 {
	memset( arrGradF2,0, 8 * sizeof( double)) ;
	 double valTemp0 = (mAltit + R_ZEMLI) /(marrStrSK_VS [1] + R_ZEMLI); // ���������� ����������
	 double valTemp1 = valTemp0 * cos(marrStrSK_VS [7]); // ���������� ����������
	 double valTemp2 = sin(marrStrSK_VS [3] ); // ���������� ����������
	arrGradF2[1] =  marrStrSK_VS [6]* valTemp1 * valTemp2 /(marrStrSK_VS [1] + R_ZEMLI);
	arrGradF2[3] =  - cos( marrStrSK_VS [3] )* marrStrSK_VS [6] *valTemp1 ;
	arrGradF2[6] =  - valTemp2 * valTemp1 ;
	arrGradF2[7] =  marrStrSK_VS [6]* sin(marrStrSK_VS [7])* valTemp0 * valTemp2;
 }

  // ���������� ������� ������� ����������� (���������) ������� F7(x1,x4,x7,x8)- ���� TETTA
 // OUTPUT:
 // arrGradF3 [8] - ������ �����������
 void TMyShellTraj::fncCalcGradF7(const  double valDerivDeltaTettaTochkaPoPsi , double *arrGradF7)
 {
	memset( arrGradF7, 0, 8 * sizeof( double)) ;
	 double valTemp0 = 1. /(marrStrSK_VS [1] + R_ZEMLI); // ���������� ����������
	 double valTemp1 = valTemp0 * cos(marrStrSK_VS [7]); // ���������� ����������

	arrGradF7[1] =  -marrStrSK_VS [6]* valTemp1 * valTemp0 ;
	arrGradF7[6] =  G_ZEMLI * cos(marrStrSK_VS [7])/marrStrSK_VS [6]/marrStrSK_VS [6] + valTemp1 ;
	arrGradF7[7] =  sin(marrStrSK_VS [7])*( G_ZEMLI /marrStrSK_VS [6] - marrStrSK_VS [6] *valTemp0 ) ;
	arrGradF7[3] = valDerivDeltaTettaTochkaPoPsi ;
 }
  // ���������� ������� ������� ����������� (���������) ������� F3  - ���
  // INPUT:
  // arrGradF7 -  �������� ������� F7
  // valF7 - �������� ������� F7
  // valDeltaPsiTochka -  ����������� ������ ���� ������� ����������
  // valDerivDeltaPsiTochkaPoPsi - ������� ����������� ����������� ������ ���� ������� ���������� �� Psi ( �� marrStrSK_VS [3]) )
  // valDerivDeltaPsiTochkaPoTetta - ������� ����������� ����������� ������ ���� ������� ���������� �� �etta ( �� marrStrSK_VS [7]) )
 // OUTPUT:
 // arrGradF3 [8] - ������ ����������� �� x
 // arrGradF3_po_z - ������ ����� ����������� �� ��������� �������� , �� ������� ������� �� z2 (x8)
 void TMyShellTraj::fncCalcGradF3( double *arrGradF7,const  double valF7,const  double valDeltaPsiTochka
		   ,const  double valDerivDeltaPsiTochkaPoPsi, const  double valDerivDeltaPsiTochkaPoTetta
		   ,const  double val_iz, const  double val_Deriv_iz
		   ,const  double valKnm,  double *arrGradKnm
		   ,  double *arrGradF3, double *arrGradF3_po_z)
 {
	memset( arrGradF3,0, 8 * sizeof( double)) ;
	memset( arrGradF3_po_z,0, 8 * sizeof( double)) ;
	// ��������� ���������� ���������
	const  double val_k0 = mShellBody.mvalIx0 * mShellBody.mL /(mShellBody.mMass * mShellBody.mDm * mShellBody.mh_gob);
	///
	const  double  val_cosx8 = cos( marrStrSK_VS [7]) ;  // ��������� �����
	const  double  val_x7cosx8 = val_cosx8 *  marrStrSK_VS [6] ; // ��������� �����

	arrGradF3 [1] = val_k0 * val_iz / val_x7cosx8 * ( arrGradKnm[1] * valF7 + valKnm  * arrGradF7 [1] )
			* marrStrSK_VS [4] ;

	arrGradF3 [3] = val_k0 * val_iz / val_x7cosx8 * valKnm  * arrGradF7 [3] +  valDerivDeltaPsiTochkaPoPsi
					 * marrStrSK_VS [4] ;

	arrGradF3 [4] = val_k0 * val_iz * valKnm * valF7 /val_x7cosx8 ;

	arrGradF3 [6] = val_k0 * val_iz / val_cosx8 * marrStrSK_VS [4]  * (
	 ( arrGradKnm[6] * valF7 + valKnm  * arrGradF7 [6] ) * marrStrSK_VS [6] - valKnm * valF7
												  ) / ( marrStrSK_VS [6] *marrStrSK_VS [6]) ;

	arrGradF3 [7] = val_k0 * val_iz / marrStrSK_VS [6] * ( arrGradF7[7] * val_cosx8 +  valF7 * sin(marrStrSK_VS [7] ) )
				   * marrStrSK_VS [4] 	/ (val_cosx8 * val_cosx8 ) +  valDerivDeltaPsiTochkaPoTetta ;

	arrGradF3_po_z [7] = val_k0 * valKnm *valF7/ val_x7cosx8 * val_Deriv_iz  * marrStrSK_VS [4] ;

 }
// ���������� ������� ������� ����������� (���������) ������� F4  - �����(������� �������� �������� )
  // INPUT:
  // valMach, arrGradMach  - ����� ���� � ��� ��������
  // val_q, arrGrad_q  - ���������� ����� � ��� ��������
  // valMxOmegax, arrGrad_MxOmegax -  ����� ������������� ������� � ��� ��������
 // OUTPUT:
 // arrGradF4 [8] - ������ ����������� �� x
 void TMyShellTraj::fncCalcGradF4(const  double valMach,  double *arrGradMach
							   ,const  double val_q,  double *arrGrad_q
							   ,const  double valMxOmegax,  double *arrGrad_MxOmegax
							   ,  double *arrGradF4)


 {
	memset( arrGradF4,0, 8 * sizeof( double)) ;
	const  double valSm = M_PI * mShellBody.mDm * mShellBody.mDm / 4. ;
	const  double val_k0  = -valSm *  mShellBody.mL* mShellBody.mL / ATM_AN0 / 2. ;
	arrGradF4 [1] =  val_k0 * marrStrSK_VS [4] * ( (arrGrad_MxOmegax [1] * val_q  + valMxOmegax * arrGrad_q[1] ) * valMach
				- valMxOmegax * val_q * arrGradMach [1]) / (valMach * valMach) ;
	arrGradF4 [4] =  val_k0  * valMxOmegax * val_q / valMach ;
	arrGradF4 [5] =  val_k0  * valMxOmegax * marrStrSK_VS [4] * arrGrad_q [5] / valMach ;
	arrGradF4 [6] =   val_k0 * marrStrSK_VS [4] * ( (arrGrad_MxOmegax [6] * val_q  + valMxOmegax * arrGrad_q[6] ) * valMach
				- valMxOmegax * val_q * arrGradMach [6]) / (valMach * valMach) ;

 }


 // ���������� ������� ������� ����������� (���������) ������� F5  - pi (������� ��������� �������)
  // INPUT:
  // valTay, valDerivTay  - ������� ����������� ����������� � �� ����������� �� ������

 // OUTPUT:
 // arrGradF5 [8] - ������ ����������� �� x
void TMyShellTraj::fncCalcGradF5(const  double valTay, const  double valDerivTay,  double *arrGradF5)
 {
	memset( arrGradF5,0, 8 * sizeof( double)) ;
	arrGradF5 [1] =  G_ZEMLI  *  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) *valDerivTay
		 / valTay / valTay /ATM_R_UNIVER ;
	arrGradF5 [5] =  -G_ZEMLI  *   marrStrSK_VS [6] * sin(marrStrSK_VS [7])/ valTay /ATM_R_UNIVER;
	arrGradF5 [6] =  -G_ZEMLI * marrStrSK_VS [5] *  sin(marrStrSK_VS [7]) / valTay /ATM_R_UNIVER;
	arrGradF5 [7] =  -G_ZEMLI * marrStrSK_VS [5] * marrStrSK_VS [6] * cos(marrStrSK_VS [7]) / valTay /ATM_R_UNIVER;
 }



 // ���������� ������� ������� ����������� (���������) ������� F6  - Vk
  // INPUT:
  // val_q, arrGrad_q  - ���������� ����� � ��� ��������
  // valCxEtal, arrGradCxEtal -  ����� ���������� ���� ������������� ������� ���������
  // val_ix, val_Deriv_ix - ������������ ����� � �� ����������� �� ����� ����
 // OUTPUT:
 // arrGradF6 [8] - ������ ����������� �� x
 // arrGradF6_po_z - ������ ����� ����������� �� ��������� �������� , �� ������� ������� �� z2 (x8)
 void TMyShellTraj::fncCalcGradF6(const  double val_q,  double *arrGrad_q
							   ,const  double valCxEtal,  double *arrGradCxEtal
							   ,const  double val_ix, const  double val_Deriv_ix
							   , double *arrGradF6, double *arrGradF6_po_z)
 {
	memset( arrGradF6,0, 8 * sizeof( double)) ;
	memset( arrGradF6_po_z,0, 8 * sizeof( double)) ;
	const  double valSm = M_PI * mShellBody.mDm * mShellBody.mDm / 4. ;
	// ��������� ���������� ���������
   const  double val_k0 = - val_ix *valSm / mShellBody.mMass ;

	///


	arrGradF6 [1] = val_k0 * ( arrGradCxEtal[1] * val_q + valCxEtal * arrGrad_q [1]) ;



	arrGradF6 [6] = val_k0 * ( arrGradCxEtal[6] * val_q + valCxEtal * arrGrad_q [6]) ;


	arrGradF6 [7] = -G_ZEMLI ;

	arrGradF6_po_z [7] = -valSm / mShellBody.mMass * valCxEtal * val_q * val_Deriv_ix ;

 }

 // ��� ����� ������ �� ���� valStepInt
 void TMyShellTraj::fncEilerStep(TEnvironment Environment,const  double valStepInt)
 {
    double arrF[8] ={0.} ;
   // ��� �������������� �������� �������
   fncCalc_F(Environment,arrF) ;
	double arrT0[64] = {0.},arrT1[64] = {0.};
   MatrxMultScalar(arrF, 8, 1, valStepInt,arrT0); // f * dt
   MtrxSumMatrx(arrT0, marrStrSK_VS,8, 1, arrT1) ;
   memcpy(marrStrSK_VS, arrT1, 8 * sizeof( double)) ;   ///

  mTCur += valStepInt ;


 }

// ������������� ��������� �� �������  valTNext
void TMyShellTraj::fncMovePhasVector(TEnvironment Environment,const double VAlStepInt, const double valTNext )
{
  int iCirc = (valTNext - mTCur) / VAlStepInt ;
   double valTTemp = mTCur + (( double)iCirc) * VAlStepInt ;
  for (int i = 0; i < iCirc; i++)
  {
	  fncEilerStep(Environment,VAlStepInt);
  }

  fncEilerStep(Environment,valTNext - valTTemp);

}


// ������������� ��������� �� �������  �������
// OUTPUT:
// valDHoriz - ��������������� ����� ����� �������
void TMyShellTraj::fncMoveShell_TO_ZeroAlt(TEnvironment Environment,const double VAlStepInt, double &valDHoriz )
{
  int iCirc = 1000. / VAlStepInt ;

  int i = 0 ;
  for ( i = 0; i < iCirc; i++)
  {
	  fncEilerStep(Environment,VAlStepInt);
	  if (marrStrSK_VS [1] < 0.) break ;


  }
   valDHoriz = sqrtl(marrStrSK_VS [0]*marrStrSK_VS [0]+ marrStrSK_VS [2]*marrStrSK_VS [2]) ;


	 //fncCalcMtrxPartialDeriv(Environment,VAlStepInt );
	// fncCalcVectPartialDeriv_CoeffForm(Environment,VAlStepInt, 0.001 ) ;
	 //fncCalcVectPartialDeriv_Mass(Environment,VAlStepInt, 0.01 * mShellBody.mMass ) ;

}


// ���������� ������� ��������� � ������� ����������� ��
// OUTPUT:
// arrVS_PrStSK [6]   - ������ ��������� � ������� ����������� ��
void TMyShellTraj::fncCalcVS_v_PrStSK(   double *arrVS_PrStSK  )
{
	arrVS_PrStSK[0] = marrStrSK_VS [0] ;
	arrVS_PrStSK[1] = marrStrSK_VS [1] ;
	arrVS_PrStSK[2] = marrStrSK_VS [2] ;
	arrVS_PrStSK[3] = marrStrSK_VS [6] * cos(marrStrSK_VS [7])* cos(marrStrSK_VS [3]);
	arrVS_PrStSK[4] = marrStrSK_VS [6] * sin(marrStrSK_VS [7]) ;
	arrVS_PrStSK[5] = - marrStrSK_VS [6] * cos(marrStrSK_VS [7])* sin(marrStrSK_VS [3]);
}



 // ���������� ������� ������� ����������� ��� �������� ��� ������� ��
 // ����������� �� � ������������� ����������� ��
 // OUTPUT:
 // arrJac_PrStSK [6*8]
void TMyShellTraj::fncCalcJacobi_PrStSK(  double *arrJac_PrStSK  )
{
	memset ( arrJac_PrStSK, 0, 48 * sizeof( double)) ;
  for (int i = 0; i < 3; i++) arrJac_PrStSK[ i * 8 + i ] = 1. ;
	arrJac_PrStSK[ 3 * 8 + 3 ] = - marrStrSK_VS [6] * cos(marrStrSK_VS [7] ) * sin (marrStrSK_VS [3]);
	arrJac_PrStSK[ 3 * 8 + 6 ] = cos(marrStrSK_VS [7] )  * cos(marrStrSK_VS [3] ) ;
	arrJac_PrStSK[ 3 * 8 + 7 ] = - marrStrSK_VS [6] * sin(marrStrSK_VS [7] ) * cos (marrStrSK_VS [3]);

	arrJac_PrStSK[ 4 * 8 + 6 ] =  sin(marrStrSK_VS [7] );
	arrJac_PrStSK[ 4 * 8 + 7 ] =  marrStrSK_VS [6] * cos(marrStrSK_VS [7] ) ;

	arrJac_PrStSK[ 5 * 8 + 3 ] =  - marrStrSK_VS [6] * cos(marrStrSK_VS [7] ) * cos (marrStrSK_VS [3]);
	arrJac_PrStSK[ 5 * 8 + 6 ] =  - cos(marrStrSK_VS [7] )  * sin(marrStrSK_VS [3] ) ;
	arrJac_PrStSK[ 5 * 8 + 7 ] =    marrStrSK_VS [6] * sin(marrStrSK_VS [7] ) * sin (marrStrSK_VS [3]);

}


// ������������� ��������� �� �������  �������
// OUTPUT:
// valDHoriz - ��������������� ����� ����� �������
void TMyShellTraj::fncMoveClass_TO_ZeroAlt_AND_ShowGraphs(TEnvironment Environment,const double VAlStepInt, wchar_t *wcharrPath1,  double &valDHoriz )
{
  const int QUANT_COLS = 9 , QUANT_POINTS_MAX = 3000;
  const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ������������ ����� ����� ����������
   wchar_t 	wcharrPath[400] = {0};;


   if (wcharrPath)
   {
   wcscpy(wcharrPath, wcharrPath1 );
   if(wcharrPath[wcslen(wcharrPath) -1] != L'\\')
   {
	   wcharrPath[wcslen(wcharrPath) ] = L'\\';
	   wcharrPath[wcslen(wcharrPath) +1] = 0;

   }


	parrBuff = new double [QUANT_COLS * QUANT_POINTS_MAX] ;
	memset (parrBuff, 0, QUANT_COLS * QUANT_POINTS_MAX * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ QUANT_COLS * lenName] ;
	memset (pwcharrFileNames, 0, QUANT_COLS * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"Z");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"Psi");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"Omega");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"PI");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"V");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Tetta");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;
	pscaleY[3] = 0.1;
	pscaleY[4] = 1000.;
	pscaleY[5] = 1.;
	pscaleY[6] = 1.;
	pscaleY[7] = 1.;
	pscaleY[8] = 1000.;
  }

  int iCirc = 1000. / VAlStepInt ;
 //  double valTTemp = mTCur + (( double)iCirc) * VAlStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  for ( i = 0; i < iCirc; i++)
  {
	  fncEilerStep(Environment, VAlStepInt);
	  if (marrStrSK_VS [1] < 0.) break ;

	  if (((double)mTCur > (valTOut + DEL_T - 1E-15 ))&& ( iNupPointsOut < QUANT_POINTS_MAX ))
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)marrStrSK_VS [0];
	   p[2] = (double)marrStrSK_VS [1];
	   p[3] = (double)marrStrSK_VS [2];
	   p[4] = (double)marrStrSK_VS [3];
	   p[5] = (double)marrStrSK_VS [4];
	   p[6] = (double)marrStrSK_VS [5];
	   p[7] = (double)marrStrSK_VS [6];
	   p[8] = (double)marrStrSK_VS [7];
		}

	  iNupPointsOut++;

	  }



  }
   valDHoriz = sqrtl(marrStrSK_VS [0]*marrStrSK_VS [0]+ marrStrSK_VS [2]*marrStrSK_VS [2]) ;

	 if (wcharrPath)
	 {
	 for (int j = 1; j < QUANT_COLS -1; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,iNupPointsOut //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,j  // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� X
								  ,pscaleY[j]  // ������� �� ��� Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,iNupPointsOut //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,2  // ����� ���������� �� ��� Y
								  ,1.//  ������� �� ��� X
								  ,1.  // ������� �� ��� Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,iNupPointsOut //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,3  // ����� ���������� �� ��� Y
								  ,1.//  ������� �� ��� X
								  ,1.  // ������� �� ��� Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete []parrBuff ;
  delete []pwcharrFileNames ;
  delete []pscaleY ;

   }

	// fncCalcMtrxPartialDeriv( Environment, VAlStepInt );
	// fncCalcVectPartialDeriv_CoeffForm(Environment, VAlStepInt, marrDelta[8] ) ;
	// fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt,marrDelta[9] * mShellBody.mMass ) ;

}





// ������������� ��������� �� �������  VAlFixedT

void TMyShellTraj::fncMoveClass_TO_FixedTime_AND_ShowGraphs(TEnvironment Environment
,const double VAlStepInt, const double VAlFixedT ,wchar_t *wcharrPath1  )
{
	const int QUANT_COLS = 9 , QUANT_POINTS  = (VAlFixedT - mTCur) / VAlStepInt ;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ������������ ����� ����� ����������
	wchar_t 	wcharrPath[400] = {0};

   if (wcharrPath1)
   {
   wcscpy(wcharrPath, wcharrPath1 );
   if(wcharrPath[wcslen(wcharrPath) -1] != L'\\')
   {
	   wcharrPath[wcslen(wcharrPath) ] = L'\\';
	   wcharrPath[wcslen(wcharrPath) +1] = 0;

   }


	parrBuff = new double [QUANT_COLS * QUANT_POINTS] ;
	memset (parrBuff, 0, QUANT_COLS * QUANT_POINTS * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ QUANT_COLS * lenName] ;
	memset (pwcharrFileNames, 0, QUANT_COLS * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"Z");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"Psi");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"Omega");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"PI");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"V");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Tetta");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;
	pscaleY[3] = 0.1;
	pscaleY[4] = 1000.;
	pscaleY[5] = 1.;
	pscaleY[6] = 1.;
	pscaleY[7] = 1.;
	pscaleY[8] = 1000.;


	 }



  int i = 0 ;

  for ( i = 0; i < QUANT_POINTS; i++)
  {
	  fncEilerStep(Environment, VAlStepInt);

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[i *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)marrStrSK_VS [0];
	   p[2] = (double)marrStrSK_VS [1];
	   p[3] = (double)marrStrSK_VS [2];
	   p[4] = (double)marrStrSK_VS [3];
	   p[5] = (double)marrStrSK_VS [4];
	   p[6] = (double)marrStrSK_VS [5];
	   p[7] = (double)marrStrSK_VS [6];
	   p[8] = (double)marrStrSK_VS [7];
		}

  }

	 if (wcharrPath1)
	 {
	 for (int j = 1; j < QUANT_COLS -1; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,j  // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� X
								  ,pscaleY[j]  // ������� �� ��� Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,2  // ����� ���������� �� ��� Y
								  ,1.//  ������� �� ��� X
								  ,1.  // ������� �� ��� Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,3  // ����� ���������� �� ��� Y
								  ,1.//  ������� �� ��� X
								  ,1.  // ������� �� ��� Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete []parrBuff ;
  delete []pwcharrFileNames ;
  delete []pscaleY ;

   }

// fncCalcMtrxPartialDeriv( Environment, VAlStepInt );
//  fncCalcVectPartialDeriv_CoeffForm(Environment, VAlStepInt, marrDelta[8] ) ;
 //  fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt,marrDelta[9] * mShellBody.mMass ) ;

}


// ���������� ������� ������� ���������� �������� ������� ��
// ���������� ������� � ������� i
void TMyShellTraj::fncCalcVectPartialDeriv(TEnvironment Environment
  ,const double VAlStepInt, const int iVarNum,  double valDelta,  double *arrVectPartialDeriv )
{
  TMyShellTraj shTrTemp = *this ;
  shTrTemp.mTCur = (*this).mTStart  ;
  shTrTemp.fncFillNachalnieUsloviaVS();
  shTrTemp.marrStrSK_VS [iVarNum] += valDelta;
  shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
  MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartialDeriv);
}

// ���������� ������� ������� ���������� �������� ������� ��
// ����� ������������ ����� ix
void TMyShellTraj::fncCalcVectPartialDeriv_Coef_Cx(TEnvironment Environment,const double VAlStepInt
	,const  double valDelta, double *arrVectPartDerivCoeff_Cx)
{
	TMyShellTraj shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.mShellBody.mplnCx.stretchDiagrAlongXY(1., ( 1. + valDelta)) ;
  shTrTemp.mCoefCx = shTrTemp.mShellBody.mplnCx.LinearValueApprox(mTet0 / M_PI * 180.) ;
  shTrTemp.mCoefCy = shTrTemp.mShellBody.mplnCiy.LinearValueApprox(mTet0 / M_PI * 180.) ;
  shTrTemp.mCoefCz = shTrTemp.mShellBody.mplnCz.LinearValueApprox(mTet0 / M_PI * 180.) ;
 // for (int i = 0; i < 31; i++)
 // {
	//	shTrTemp.mShellBody.marr_ix[i][1] =  mShellBody.marr_ix[i][1] * ( 1. + valDelta);
 //}

  shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartDerivCoeff_Cx);
}

//-------------------------------------------------------------------------------------
// ���������� ������� ������� ���������� �������� ������� ��
// ����� Cz
void TMyShellTraj::fncCalcVectPartialDeriv_Coef_Cz(TEnvironment Environment,const double VAlStepInt
	,const  double valDelta, double *arrVectPartDerivCoeff_Cz)
{
	TMyShellTraj shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.mShellBody.mplnCz.stretchDiagrAlongXY(1., ( 1. + valDelta)) ;
  shTrTemp.mCoefCx = shTrTemp.mShellBody.mplnCx.LinearValueApprox(mTet0 / M_PI * 180.) ;
  shTrTemp.mCoefCy = shTrTemp.mShellBody.mplnCiy.LinearValueApprox(mTet0 / M_PI * 180.) ;
  shTrTemp.mCoefCz = shTrTemp.mShellBody.mplnCz.LinearValueApprox(mTet0 / M_PI * 180.) ;
 // for (int i = 0; i < 31; i++)
 // {
	//	shTrTemp.mShellBody.marr_ix[i][1] =  mShellBody.marr_ix[i][1] * ( 1. + valDelta);
 //}

  shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartDerivCoeff_Cz);
}

// ���������� ������� ������� ���������� �������� ������� ��
// �����
void TMyShellTraj::fncCalcVectPartialDeriv_Mass(TEnvironment Environment
	,const double VAlStepInt,const  double valDelta, double *arrVectPartDerivMass)
{
  TMyShellTraj shTrTemp = *this ;
  shTrTemp.mTCur = (*this).mTStart  ;
  shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.mShellBody.mMass = mShellBody.mMass + valDelta;
	//shTrTemp.mShellBody.mV0 = (shTrTemp.mShellBody.mMass - valDelta )/ shTrTemp.mShellBody.mMass * shTrTemp.mShellBody.mV0;
	//shTrTemp.mShellBody.mOmega0 = (shTrTemp.mShellBody.mMass - valDelta )/ shTrTemp.mShellBody.mMass * shTrTemp.mShellBody.mOmega0;
	shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartDerivMass);
}
/*
// ������ ������� ������� ����������� �� ��� ��������
void TMyShellTraj::fncCalcMtrxPartialDeriv(TEnvironment Environment, const double VAlStepInt )
{
 // arrStrSK_JacT - ��������������� ������ ��� �������� ����������. ��������� �������� ���������.
	 // �������, ��� ��������� �������  marrStrSK_Jac  ������ arrStrSK_JacT ���� ���������������
	double	arrStrSK_JacT[64]= {0.} ;

	for (int i = 0; i < 8; i++)
	{
	 fncCalcVectPartialDeriv(Environment,VAlStepInt, i, marrDelta[i], &arrStrSK_JacT[i*8] ) ;
	}

	MatrTransp(arrStrSK_JacT, 8, 8, marrStrSK_Jac);

 }
 */

// ������ ������� ������� ����������� �� ���������� �����
// arrMtrx_Wind_PartialDeriv[8*3]  - ������� ������� ����������� �� ����� ��������, ������������, ������ ��������
void TMyShellTraj::fncCalcMtrxTransp_Wind_PartialDeriv(TEnvironment Environment
, const double VAlStepInt, double *arrMtrxTransp_Wind_PartialDeriv )
{
	// arrStrSK_JacT - ��������������� ������ ��� �������� ����������. ��������� �������� ���������.
	// �������, ��� ��������� �������  marrStrSK_Jac  ������ arrStrSK_JacT ���� ���������������

	memset(arrMtrxTransp_Wind_PartialDeriv, 0, 24 * sizeof(double));

	TEnvironment Environment1 =  Environment;
	Environment1.mWind_V =  Environment.mWind_V + 0.1 ;
	TMyShellTraj shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.fncMovePhasVector(Environment1,VAlStepInt, mTCur );
	double arrT0[8] = {0.} ;
	MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1.0 /0.1 ,arrMtrxTransp_Wind_PartialDeriv);

	Environment1 =  Environment;
	Environment1.mWind_Alf =  Environment1.mWind_Alf + 0.005;
	shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.fncMovePhasVector(Environment1,VAlStepInt, mTCur );
	MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /0.005 ,&arrMtrxTransp_Wind_PartialDeriv[8]);

	Environment1 =  Environment;
	Environment1.mWind_VertV =  Environment.mWind_VertV + 0.1 ;
	shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.fncMovePhasVector(Environment1,VAlStepInt, mTCur );
	MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /0.1  ,&arrMtrxTransp_Wind_PartialDeriv[16]);

 }

 // --------------------------------------------------------------------------------------------------
 // �������� ������� �� ��� � ���
 // ���� LEnArrVS = 3 �� ��������������� ������ ���������
 // ���� LEnArrVS = 6 �� ��������������� ������ ���������  � ������ ��������
 //
void TMyShellTraj::transform_xyzGSK_To_xyzSSK( const int LEnArrVS, double *arrGSKInp, double *arrSSKOut)
{
	arrSSKOut[1] = arrGSKInp[2];
	arrSSKOut[0] = sin(mAlfDir) * arrGSKInp[0] + cos(mAlfDir) * arrGSKInp[1];
	arrSSKOut[2] = cos(mAlfDir) * arrGSKInp[0] - sin(mAlfDir) * arrGSKInp[1];
	if(LEnArrVS == 6)
	{
	arrSSKOut[4] = arrGSKInp[5];
	arrSSKOut[3] = sin(mAlfDir) * arrGSKInp[3] + cos(mAlfDir) * arrGSKInp[4];
	arrSSKOut[5] = cos(mAlfDir) * arrGSKInp[3] - sin(mAlfDir) * arrGSKInp[4];

		}
}

// ------------------------------------------------------------------------------
 // �������� ������� �� ������������� ��� �  ���
 // ���� LEnArrVS = 3 �� ��������������� ������ ���������
 // ���� LEnArrVS = 6 �� ��������������� ������ ���������  � ������ ��������

 void TMyShellTraj::transform_xyzSSK_To_xyzGSK( const int LEnArrVS, double *arrSSKInp, double *arrGSKOut)
{
	arrGSKOut[2] = arrSSKInp[1];
	arrGSKOut[0] = sin(mAlfDir) * arrSSKInp[0] + cos(mAlfDir) * arrSSKInp[2];
	arrGSKOut[1] = cos(mAlfDir) * arrSSKInp[0] - sin(mAlfDir) * arrSSKInp[2];
	if(LEnArrVS == 6)
	{
	arrGSKOut[5] = arrSSKInp[4];
	arrGSKOut[3] = sin(mAlfDir) * arrSSKInp[3] + cos(mAlfDir) * arrSSKInp[5];
	arrGSKOut[4] = cos(mAlfDir) * arrSSKInp[3] - sin(mAlfDir) * arrSSKInp[5];

    }
}
// ------------------------------------------------------------------------------
// ������������ ������� �������� �� ������� ����������� �� � ���
// ���� LEnArrVS = 3 �� arrMtrxTransformOut ������������� ������ ���������
 // ���� LEnArrVS  = 6 �� arrMtrxTransformOut  ������������� ������ ���������  � ������ ��������
 void TMyShellTraj::createMtrxTransform_xyzSSK_To_xyzGSK( const int LEnArrVS,  double *arrMtrxTransformOut)
 {
	 if(LEnArrVS == 3)
	 {
		 memset(arrMtrxTransformOut,0,  9 * sizeof(double));
		 arrMtrxTransformOut [0] = sin(mAlfDir);
		 arrMtrxTransformOut [2] = cos(mAlfDir);
		 arrMtrxTransformOut [3] = cos(mAlfDir);
		 arrMtrxTransformOut [5] = -sin(mAlfDir);
		 arrMtrxTransformOut [7] = 1.;
		 return;
	 }

	 if(LEnArrVS == 6)
	 {
		 memset(arrMtrxTransformOut, 0,  36 * sizeof(double));
		 arrMtrxTransformOut [0] = sin(mAlfDir);
		 arrMtrxTransformOut [2] = cos(mAlfDir);
		 arrMtrxTransformOut [6] = cos(mAlfDir);
		 arrMtrxTransformOut [8] = -sin(mAlfDir);
		 arrMtrxTransformOut [13] = 1.;

		 arrMtrxTransformOut [21] = sin(mAlfDir);
		 arrMtrxTransformOut [23] = cos(mAlfDir);
		 arrMtrxTransformOut [27] = cos(mAlfDir);
		 arrMtrxTransformOut [29] = -sin(mAlfDir);
		 arrMtrxTransformOut [34] = 1.;
		 return;
	 }
 }


// ------------------------------------------------------------------------------

// ���������� ����� ������������ ���������� ����� ������������ ������� � ����
// INPUT:
// arrTargVS_SSK0[6] - ������ ��������� ���� � ��� �� ������ ������ �������� �������
// VAl_dtInt - ��� ��������������
// ����������:
//  ����������� ����������
double TMyShellTraj::calcPointMissMinimum(TEnvironment Environment, double *arrTargVS_SSK0, const double VAl_dtInt)
{

	double valF0prev = 10000000.0; // �������- ��-� �����. ����� �������� � ����� �� ������ T
	double valF0; // �������- ��-� �����. ����� �������� � ����� �� ������ T+dt

	double valMaxT = 600.;
	int iMqxQuantIter = valMaxT/ VAl_dtInt;
	for (int i =0; i < iMqxQuantIter; i++)
	{
		fncEilerStep( Environment, VAl_dtInt) ;
		// ���������� ����������� �� ������ valTCur ��������� ����
		double arrTargPosExtr[3] ={0.}, arrT[3] ={0.};
		MatrxMultScalar(&arrTargVS_SSK0[3], 3, 1, mTCur,arrT);
		MtrxSumMatrx(arrTargVS_SSK0, arrT,3, 1, arrTargPosExtr) ;
		///

		double arrDelt[3] ={0.};
		MtrxMinusMatrx(arrTargPosExtr, marrStrSK_VS ,3, 1, arrDelt);
		valF0 = Norm3( arrDelt) ;
		///

		if (valF0 > valF0prev)
		{
			break;
		}
		 valF0prev = valF0;
	}

	return  valF0;
}


// ���������� ������� ������� ����������� �������� ������� � ����������� �������� ���������
// �� ���������� �������� - ��������� ��������, ����� �����, ���� �����
void TMyShellTraj::calcJacobian_8x10 (TEnvironment Environment,const double VAlStepInt
  , double *arrStrSK_Jacobian)
{
	 // LEN_ARR_SCATTERS - ����� ������� ���������� �� �������� ��������� ��������
		// ��������� �� ������� ��������� ��������
// 0. marrStrSK_VS [3]-  ���� ��� (�������)
// 1. marrStrSK_VS [5]-  ������� ��������� ��������� ��
// 2. marrStrSK_VS [7]-  ���� �����
// 3. marrStrSK_VS [6]-  ��������
// 4. �����
// 5.����� ����� Cx
// 6. ���� ����� �� ��� Z
// 7. ������ �������������� �������� �����
// 8. ����������� ��������������� �����
// 9. ������ ������������� �����

		double	arrStrSK_JacT[LEN_ARR_SCATTERS * 8]= {0.};

		double arrDelta[7] = {
								 0.001
								 ,0.01
								 ,0.001
								 ,1.
								 ,0.01
								 ,0.01
								 ,0.01
								 } ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 3, arrDelta[0], arrStrSK_JacT ) ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 5, arrDelta[1], &arrStrSK_JacT[8] ) ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 7, arrDelta[2], &arrStrSK_JacT[2 *8] ) ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 6, arrDelta[3], &arrStrSK_JacT[3 *8] ) ;
	fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt, arrDelta[4]* mShellBody.mMass, &arrStrSK_JacT[4 * 8] )  ;
	fncCalcVectPartialDeriv_Coef_Cx( Environment,VAlStepInt
	, arrDelta[5],  &arrStrSK_JacT[5 * 8] )  ;
	fncCalcVectPartialDeriv_Coef_Cz(Environment, VAlStepInt
	,  arrDelta[4],  &arrStrSK_JacT[6 * 8] )  ;
	fncCalcMtrxTransp_Wind_PartialDeriv( Environment
	,VAlStepInt, &arrStrSK_JacT[7* 8]) ;

	MatrTransp(arrStrSK_JacT, LEN_ARR_SCATTERS, 8, arrStrSK_Jacobian);

}


// ���������� �������������� ������� �������� ������� ��������� �������
// � ���
// INPUT:
// Environment -������� �����
// VAlStepInt - ��� ��������������
// VAlTFlight - �������� �����
// arrMtrxShellDisp  - ������������ ������� ���������� ��������� ����������
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - ������� ������� ����������� �������� �������
// ������� � ������������� �������� ��������� �� ����������
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - �������������� ������� � ���
// arrShellScatteringsCorMtrxPos_SSK [3*3] - �������������� ������� ������ ��������� ��������� � ���
// arrShellVS_GSK[6] - ������ ��������� ������� � ���
void TMyShellTraj::calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (TEnvironment Environment,const double VAlStepInt
	 ,const double VAlTFlight,double* arrMtrxShellDisp,  double *arrStrSK_Jacobian
	,double* arrShellScatteringsCorMtarx_GSK, double *arrShellVS_GSK, double* arrShellScatteringsCorMtrxPos_SSK)
{
	 fncMovePhasVector(Environment, VAlStepInt, VAlTFlight);
	 calcJacobian_8x10 ( Environment, VAlStepInt,  arrStrSK_Jacobian);

	 // ���������� ������ �������� ������� ��������� ������� � ����������� �� - arrCorMtrxTrajCK [8*8]
	 double arrT0[LEN_ARR_SCATTERS * 8] = {0.}, arrCorMtrxTrajCK[64] = {0.};
		MtrxMultMatrx(arrStrSK_Jacobian,8, LEN_ARR_SCATTERS, arrMtrxShellDisp,LEN_ARR_SCATTERS, arrT0) ;
		MtrxMultMatrxTransp(arrT0, 8, LEN_ARR_SCATTERS, arrStrSK_Jacobian,8, arrCorMtrxTrajCK) ;
		///

		// ���������� ������� ������� ����������� ������� ������� ��������� ������� ��
		// ����������� �� � �������. ��
		double arrJac_PrStSK [ 6*8] = {0.};
		fncCalcJacobi_PrStSK(  arrJac_PrStSK  );
		///

		// ���������� ��������. ������� ������ � ������� ����������� ��  - arrCorMtrxPrStCK
		double arrTemp0 [ 6*8] = {0.}, arrCorMtrxPrStCK [6*6] = {0.};
		MtrxMultMatrx(arrJac_PrStSK,6, 8, arrCorMtrxTrajCK ,8, arrTemp0) ;
		MtrxMultMatrxTransp(arrTemp0, 6, 8, arrJac_PrStSK,6, arrCorMtrxPrStCK) ;
		for (int i =0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				arrShellScatteringsCorMtrxPos_SSK[ i * 3 + j] =   arrCorMtrxPrStCK [ i * 6 + j];
			}
		}
		///

		// ���������� ������� �������� �� ������� ����������� �� � ���
		double arrMtrxTransf_From_StPrSK_To_GSK [36];
		createMtrxTransform_xyzSSK_To_xyzGSK( 6, arrMtrxTransf_From_StPrSK_To_GSK);
		double arrTemp1 [ 36] = {0.};
		MtrxMultMatrx(arrMtrxTransf_From_StPrSK_To_GSK,6, 6, arrCorMtrxPrStCK ,6, arrTemp1) ;
		MtrxMultMatrxTransp(arrTemp1, 6, 6, arrMtrxTransf_From_StPrSK_To_GSK,6, arrShellScatteringsCorMtarx_GSK) ;
		///

		// ���������� ������� ��������� ������� � ���
		double arrVS_PrStSK [6] = {0.};
		fncCalcVS_v_PrStSK(arrVS_PrStSK) ;
		///
		MtrxMultMatrx(arrMtrxTransf_From_StPrSK_To_GSK,6, 6, arrVS_PrStSK ,1, arrShellVS_GSK) ;

}




// ����������  ������� ��������� �������   � ���
// � ������� ����� (������� ����������� ������� ��������� � �������� �����)
// �� 10 ����������� ����������
// INPUT:
// Environment -������� �����
// VAlStepInt - ��� ��������������
// VAlTFlight - �������� �����
// arrMtrxShellDisp  - ������������ ������� ���������� ��������� ����������
//OUTPUT:
// arrGSK_Jacobian [ 6 * LEN_ARR_SCATTERS ] - ������� ������� ����������� �������� �������
// ������� � ��� �� ����������
// arrShellVS_GSK[6] - ������ ��������� ������� � ���
void TMyShellTraj::calc_VS_GSK_And_Jacobian_6x10_GSK (TEnvironment Environment
	 ,const double  VAlStepInt,const double VAlTFlight,  double *arrGSK_Jacobian
	, double *arrShellVS_GSK)
{

	double arrStrSK_Jacobian[ 8 *  LEN_ARR_SCATTERS] = {0.};
	fncMovePhasVector(Environment, VAlStepInt, VAlTFlight);
	calcJacobian_8x10 ( Environment, VAlStepInt,  arrStrSK_Jacobian);

	// ���������� ������� ������� ����������� ������� ������� ��������� ������� ��
	// ����������� �� � �������. ��
	double arrJac_PrStSK [ 6*8] = {0.};
	fncCalcJacobi_PrStSK(  arrJac_PrStSK  );
	///



	// ���������� ��������. ������� ������ � ������� ����������� ��  - arrCorMtrxPrStCK
	double arrTemp [ 6 * LEN_ARR_SCATTERS ]  = {0.};

	MtrxMultMatrx(arrJac_PrStSK,6, 8, arrStrSK_Jacobian ,LEN_ARR_SCATTERS, arrTemp) ;
	///

	// ���������� ������� �������� �� ������� ����������� �� � ���
	double arrMtrxTransf_From_StPrSK_To_GSK [36];
	createMtrxTransform_xyzSSK_To_xyzGSK( 6, arrMtrxTransf_From_StPrSK_To_GSK);
	MtrxMultMatrx(arrMtrxTransf_From_StPrSK_To_GSK,6, 6, arrTemp ,LEN_ARR_SCATTERS, arrGSK_Jacobian) ;


	// ���������� ������� ��������� ������� � ���
	double arrVS_PrStSK [6] = {0.};
	fncCalcVS_v_PrStSK(arrVS_PrStSK) ;
	///
	MtrxMultMatrx(arrMtrxTransf_From_StPrSK_To_GSK,6, 6, arrVS_PrStSK ,1, arrShellVS_GSK) ;

}



double  TMyShellTraj::fncPartialDerivD_TochkiPadenia_po_Tetta0(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);
	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TMyShellTraj ShellTraj1(arrVesselVelocity, mShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TMyShellTraj ShellTraj2(arrVesselVelocity, mShellBody, VAlTetta0 + 0.001, M_PI/2.);
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) / 0.001;
	 return valPartialDeriv;
}



double  TMyShellTraj::fncPartialDerivD_TochkiPadenia_po_Mass(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TMyShellTraj ShellTraj1(arrVesselVelocity, mShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TMyShellTraj ShellTraj2(arrVesselVelocity, mShellBody, VAlTetta0 , M_PI/2.);
	 double delm = ShellTraj1.mShellBody.mMass  * 0.01;
	 ShellTraj2.mShellBody.mMass += delm ;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) / delm ;
	 return valPartialDeriv;
}


double  TMyShellTraj::fncPartialDerivD_TochkiPadenia_po_V0(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TMyShellTraj ShellTraj1(arrVesselVelocity, mShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TMyShellTraj ShellTraj2(arrVesselVelocity, mShellBody, VAlTetta0 , M_PI/2.);

	 ShellTraj2. marrStrSK_VS[6] += 1.;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) ;
	 return valPartialDeriv;
}


double  TMyShellTraj::fncPartialDerivD_TochkiPadenia_po_Cx(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TMyShellTraj ShellTraj1(arrVesselVelocity, mShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TMyShellTraj ShellTraj2(arrVesselVelocity, mShellBody, VAlTetta0 , M_PI/2.);

	ShellTraj2.mShellBody.mplnCx = ShellTraj1.mShellBody.mplnCx.MultScalar(1.01 ) ;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	double  valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) /0.01;

	 return valPartialDeriv;
}

double  TMyShellTraj::fncPartialDerivZ_TochkiPadenia_po_Psi0(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TMyShellTraj ShellTraj1(arrVesselVelocity, mShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TMyShellTraj ShellTraj2(arrVesselVelocity, mShellBody, VAlTetta0 , M_PI/2.);
	 ShellTraj2. marrStrSK_VS[3] += 0.001;

	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	double  valPartialDeriv = (ShellTraj2.marrStrSK_VS[2] -  ShellTraj1.marrStrSK_VS[2] ) /0.001;

	 return valPartialDeriv;
}







#pragma package(smart_init)


