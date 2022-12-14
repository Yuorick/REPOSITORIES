//---------------------------------------------------------------------------


#pragma hdrstop

#include "LearnShellTraj.h"

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


	// ??????????? ??
  // ??????? ???? ??? ????????????? ?? ??? OX ??????????? ?? ?????? ??????? ??????? !!!
		 //TShellBody  mShellBodyBody
//---------------------------------------------------------------------------
TLearnShellTraj::TLearnShellTraj()
{
	mTStart = 0.;  // ??? ??????
	mTet0 = 17./180.* M_PI;  // ??? ???? ????????  17 ????


	mAlfDir =  M_PI / 2.; //????????????? ???? ????????? ????????
	mPsi0 = 0.; //??? ?????? (???? ????) ? ??????????? ?? mPsi0

	mAltit =0.; // ??? ??????
	mLatitude = M_PI/3.;  // ??????
	mTCur =0. ;  // ?????? ????????? ???????? ???????
  //	StepInt = 0.0001 ;
	mLearnShellBody = TLearnShellBody () ;  // ??????
	// ?????? ????? ?? ?????? ?????????? ???????? ???????  - ??? ???????? ??????? ???????????
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
	// ???????  ??????  ? ??????????? ??:
//  marrStrSK_VS [0]- x
//  marrStrSK_VS [1]-  y
//  marrStrSK_VS [2]-  z
//  marrStrSK_VS [3]-  ???? ???
//  marrStrSK_VS [4]-  ??????? ????????? ?????
//  marrStrSK_VS [5]-  ????????? ????????? ??
//  marrStrSK_VS [6]-  ??????? ???????? V
//  marrStrSK_VS [7]-  ???? ?????
	fncFillNachalnieUsloviaVS();


}

//---------------------------------------------------------------------------


// ??????????? ???????????
 TLearnShellTraj ::TLearnShellTraj (const TLearnShellTraj &R)
 {
	mTStart = R.mTStart;  // ??? ??????
	mTet0 = R.mTet0;  // ??? ???? ????????  17 ????
	mAlfDir = R.mAlfDir; //??? ?????? ????????  ??????
	mPsi0 = R.mPsi0 ;  //??? ?????? (???? ????) ? ??????????? ?? mPsi0

  //	StepInt = R.StepInt ;


	mAltit = R.mAltit; // ??? ??????
	mLatitude = R.mLatitude;  // ??????
	mTCur = R.mTCur ;  // ?????? ????????? ???????? ???????

	mLearnShellBody = R.mLearnShellBody ;  // ??????

	memcpy(marrStrSK_VS, R.marrStrSK_VS, 8 * sizeof( double));
 //	memcpy(marrStrSK_Jac, R.marrStrSK_Jac, 64 * sizeof( double));
 //	memcpy(marrStrSK_JacCoeffFom, R.marrStrSK_JacCoeffFom, 8 * sizeof( double));
 //	memcpy(marrStrSK_JacMass, R.marrStrSK_JacMass, 8 * sizeof( double));
//	memcpy(marrDelta, R.marrDelta, 10 * sizeof( double));
//	memcpy(marrDispWindParamsScatters, R.marrDispWindParamsScatters , 3 * sizeof(double));

 }

 // ???????? ????????????
 TLearnShellTraj TLearnShellTraj::operator=(TLearnShellTraj  R)
 {
	mTStart = R.mTStart;  // ??? ??????
	mTet0 = R.mTet0;  // ??? ???? ????????  17 ????
	mAlfDir = R.mAlfDir; //??? ?????? ????????  ??????
	mPsi0 = R.mPsi0 ;  //??? ?????? (???? ????) ? ??????????? ?? mPsi0




	mAltit = R.mAltit; // ??? ??????
	mLatitude = R.mLatitude;  // ??????
	mTCur = R.mTCur ;  // ?????? ????????? ???????? ???????

	mLearnShellBody = R.mLearnShellBody ; // ??????

	memcpy(marrStrSK_VS, R.marrStrSK_VS, 8 * sizeof( double));
 //	memcpy(marrStrSK_Jac, R.marrStrSK_Jac, 64 * sizeof( double));
 //	memcpy(marrStrSK_JacCoeffFom, R.marrStrSK_JacCoeffFom, 8 * sizeof( double));
 //	memcpy(marrStrSK_JacMass, R.marrStrSK_JacMass, 8 * sizeof( double));
 //	memcpy(marrDelta, R.marrDelta, 10 * sizeof( double));
 //	memcpy(marrDispWindParamsScatters, R.marrDispWindParamsScatters , 3 * sizeof(double));

	return *this ;
 }


  // ????? ???????????
 TLearnShellTraj::TLearnShellTraj (const  double TStart, const  double Tet0
			  , const  double AlfDir
			  ,const  double Psi0 ,const  double  Altit,const  double Latitude
			  ,const  double  TCur)
 {
	 mTStart = TStart;  // ??? ??????
	 mTet0 = Tet0;  // ??? ???? ????????
   //	 double mBet0 ; //??? ?????? ????????
	 mAlfDir = AlfDir ; //????????????? ???? ????????? ????????
	 mPsi0 = Psi0 ; // ??? ?????? (???? ????)
 
	 mAltit = Altit ; // ??? ??????
	 mLatitude = Latitude ;  // ??????
	 mTCur = TCur ;  // ?????? ????????? ???????? ???????

	 mLearnShellBody = TLearnShellBody() ;  // ??????


	 fncFillNachalnieUsloviaVS();

 }

  // ????? ??????????? 2
 TLearnShellTraj::TLearnShellTraj (const TLearnShellBody LearnShellBody, const  double TStart, const  double Tet0
			 , const  double AlfDir
			  ,const  double Psi0 ,const  double  Altit,const  double Latitude
			  ,const  double  TCur)
 {
	 mTStart = TStart;  // ??? ??????
	 mTet0 = Tet0;  // ??? ???? ????????
   //	 double mBet0 ; //??? ?????? ????????
	 mAlfDir = AlfDir ; //????????????? ???? ????????? ????????
	 mPsi0 = Psi0 ; // ??? ?????? (???? ????)




	 mAltit = Altit ; // ??? ??????
	 mLatitude = Latitude ;  // ??????
	 mTCur = TCur ;  // ?????? ????????? ???????? ???????

	 mLearnShellBody = LearnShellBody;  // ??????

	 fncFillNachalnieUsloviaVS();

 }
 ///


  // ????? ??????????? 3
 TLearnShellTraj::TLearnShellTraj (const TLearnShellBody LearnShellBody, const  double TStart, const  double Tet0
			 , const  double AlfDir,const  double Psi0 ,const  double  Altit,const  double Latitude
			  ,const  double  TCur, double *arrDispWindParamsScatters)
 {
	 mTStart = TStart;  // ??? ??????
	 mTet0 = Tet0;  // ??? ???? ????????
   //	 double mBet0 ; //??? ?????? ????????
	 mAlfDir = AlfDir ; //????????????? ???? ????????? ????????
	 mPsi0 = Psi0 ; // ??? ?????? (???? ????)




	 mAltit = Altit ; // ??? ??????
	 mLatitude = Latitude ;  // ??????
	 mTCur = TCur ;  // ?????? ????????? ???????? ???????

	 mLearnShellBody = LearnShellBody;  // ??????

	 fncFillNachalnieUsloviaVS();

 }
 ///


  // ????? ??????????? 4
  // Eps0 - ???? ????? ???????? ? ????
  // Bet0-  ????. ????  ???????? ? ????
  // arrVesselVElocity [3] - ???????? ??????? (????)
  // LearnShellBody -  ???????
 TLearnShellTraj::TLearnShellTraj (double *arrVesselVelocity, const TLearnShellBody LearnShellBody,  const  double Eps0,  const  double Bet0 )


 {
	 mTStart = 0.;  // ??? ??????
	 mTCur = 0.;
	 // ?????? ????????? ???????? ? ???
	  double arrGSKV0[3] = {0.};
	  arrGSKV0[0] =  arrVesselVelocity[0] + LearnShellBody.mV0 * cos(Eps0) * sin( Bet0);
	  arrGSKV0[1] =  arrVesselVelocity[1] + LearnShellBody.mV0 * cos(Eps0) * cos( Bet0);
	  arrGSKV0[2] =  arrVesselVelocity[2] + LearnShellBody.mV0 * sin(Eps0) ;
	 ///
	 //
	 mLearnShellBody = LearnShellBody;  // ??????
	 mLearnShellBody.mV0 =  Norm3( arrGSKV0) ;

	 ///

	 mTet0 = arcSin( arrGSKV0[2]/ mLearnShellBody.mV0);
	 mAlfDir = atan2( arrGSKV0[0], arrGSKV0[1]);
	 mPsi0 = 0. ; // ??? ?????? (???? ????)

	 mAltit = 0. ; // ??? ??????
	 mLatitude = M_PI /3. ;  // ??????

	 fncFillNachalnieUsloviaVS();
 }
//



void TLearnShellTraj::fncFillNachalnieUsloviaVS()
{
  marrStrSK_VS[0] = 0.;
  marrStrSK_VS[1] = mAltit ;
  marrStrSK_VS[2] = 0.;
  marrStrSK_VS[3] =  mPsi0 ;
  marrStrSK_VS[4] = mLearnShellBody.mOmega0 ;
 // marrStrSK_VS[5] = ATM_PN0 ;//
  marrStrSK_VS[5] = 1 ;//
  marrStrSK_VS[6] = mLearnShellBody.mV0 ;
  marrStrSK_VS[7] = mTet0 ;


}



// ?????????? ?????? ??????? ?????? ????? ??????? ??? ????????? -arrF
// ??????? ??????? ??????????? ?????? ????? ?? x  - mtr_dF_po_dx
// ??????? ??????? ??????????? ?????? ????? ?? ????????? ????????   - mtr_dF_po_dz
void TLearnShellTraj::fncCalc_F_and_H_and_HI( double *arrF,  double *mtr_dF_po_dx
	,  double *mtr_dF_po_dz, double *mtr_dF_po_di)
{
	/* memset ( mtr_dF_po_dz, 0, 64 * sizeof( double)) ;
   memset ( mtr_dF_po_di, 0, 16 * sizeof( double)) ;
  // ??????? ??????
   const  double valSm = M_PI * mLearnShellBody.mDm * mLearnShellBody.mDm / 4. ;
  // ?????????? ?????????? ??????????? ???????????  ? ?? ?????????
   double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///

  // ?????????? ????? ???? ? ??????? ??? ??????? ???????????
   double valMach = 0.,arrGradMach[8] ={0.} ;
  fncCalcMach_and_GradMach(valTay,valDerivTay, valMach, arrGradMach);
  ///

  // ?????????? ?????? q ? ??????? ?? ?????????
   double val_q = 0.,arrGrad_q[8] ={0.} ;
  fncCalc_q_and_Grad_q(valMach,arrGradMach, val_q, arrGrad_q);
  ///

  // ??????? ??????????? ????????
   double valVpr = marrStrSK_VS[6] * cos(marrStrSK_VS [7]) * (R_ZEMLI + mAltit)/
	   (R_ZEMLI + marrStrSK_VS [1])  ;
	   ///

  // ?????????? ?????????????? ??????? Knm ? ??????? ?? ??????????? ?? x
  // ?? ??????? 1-?? ? 6-?? ?????????? ( ???????? ?????????? ? 0)
	 double valKnm = 0., arrGradKnm[8] = {0.} ;
	fncCalcKnm_and_Grad_Knm( valMach, arrGradMach
			   ,valKnm, arrGradKnm ) ;
	///

 // ?????????? ??????????????? ?????? iz(z2) ? ?? ???????????
	 double val_iz = 0., val_Deriv_iz = 0. ;
	mLearnShellBody.fnkIz0(mTet0,val_iz, val_Deriv_iz) ;
	///

// ?????????? ??????????????? ?????? ix(z2) ? ?? ???????????
	 double val_ix = 0., val_Deriv_ix = 0. ;
	mLearnShellBody.fnkIx0(mTet0,val_ix, val_Deriv_ix) ;
	///

// ?????????? ?????????????? ??????? mxwx  ? ?? ??????????? ?? M
	 double valMxOmegax = 0., arrGradMxOmegax[8] = {0.} ;
	fncCalcMxOmegax_and_Grad_MxOmegax( valMach, arrGradMach
			   ,valMxOmegax, arrGradMxOmegax ) ;
	///

// ?????????? ?????????????? ??????? CxEtal  ? ?? ??????????? ?? M
	 double valCxEtal = 0., arrGradCxEtal [8] = {0.} ;
	fncCalcCxEtal_and_Grad_CxEtal(valMach, arrGradMach
			   ,valCxEtal, arrGradCxEtal ) ;
	///

// ?????????? ?????????? ?????????? ???????????? ?????? ????? ??? ????? ???????? ?????
    double valDeltaTettaTochka =0.          //???????  ???????? ???? ??????? ?????? DeltaTettaTochka
			  ,valDerivDeltaTettaTochkaPoPsi= 0. // ??????????? DeltaTettaTochka ?? Psi (= marrStrSK_VS [3])
			  ,valDeltaPsiTochka =0.           //???????  ???????? ???? ???? DeltaPsiTochka
			  ,valDerivDeltaPsiTochkaPoPsi= 0. // ??????????? DeltaPsiTochka ?? Psi (= marrStrSK_VS [3])
			  ,valDerivDeltaPsiTochkaPoTetta= 0.;  // ??????????? DeltaPsiTochka ?? Tetta (=  marrStrSK_VS [7])
   fncCalcDeltaTettaTochka_and_DerivPoPsi(valDeltaTettaTochka, valDerivDeltaTettaTochkaPoPsi) ;
   fncCalcDeltaPsiTochka_and_DerivPoPsi_and_DerivPoTetta(valDeltaPsiTochka
		   ,valDerivDeltaPsiTochkaPoPsi, valDerivDeltaPsiTochkaPoTetta);
	 ///

   arrF[0] =  valVpr * cos(marrStrSK_VS [3]) ;
   arrF[1] =  marrStrSK_VS [6] * sin(marrStrSK_VS [7]) ;
   arrF[2] =  -valVpr * sin(marrStrSK_VS [3]) ;

  /// mtr_dF_po_di[3 * 2 + 1] = mLearnShellBody.mvalIx0 * mLearnShellBody.mL * val_iz *valKnm * marrStrSK_VS [4] *
   //		(-G_ZEMLI/ marrStrSK_VS [6]/ marrStrSK_VS [6] + 1./(R_ZEMLI + marrStrSK_VS [1]))/
   //		(mLearnShellBody.mMass * mLearnShellBody.mDm * mLearnShellBody.mh_gob);
  // arrF[3] =  mtr_dF_po_di[3 * 2 + 1] + valDeltaPsiTochka;
  // arrF[3] =  (mLearnShellBody.mvalIx0 * mLearnShellBody.mL / (mLearnShellBody.mMass * mLearnShellBody.mDm * mLearnShellBody.mh_gob))
	//	  * val_iz *valKnm * marrStrSK_VS [4]
	 //	  *(-G_ZEMLI/ marrStrSK_VS [6]/ marrStrSK_VS [6] + 1./(R_ZEMLI + marrStrSK_VS [1]))+ valDeltaPsiTochka;
  // mtr_dF_po_di[3 * 2 + 1] = ( arrF[3] - valDeltaPsiTochka)/ val_iz ;
   // ????????? ?????????? ??? ????????? ?????
  // const  double val_k0 = -valSm * mLearnShellBody.mL * mLearnShellBody.mL/ mLearnShellBody.mvalIx0 / ATM_AN0 / 2.;
  // arrF[4] = val_k0 *valMxOmegax *val_q * marrStrSK_VS [4] / valMach ;
	 const  double val_qsbm = 0.474 * mLearnShellBody.mDm* mLearnShellBody.mDm/ mLearnShellBody.mMass * marrStrSK_VS [6] * marrStrSK_VS [6] * marrStrSK_VS[5];
	 const  double val_mmx = (valMxOmegax * mLearnShellBody.mDm /2./ mLearnShellBody.mV0);
	 const  double val_mn2  = val_mmx* val_qsbm  * mLearnShellBody.mMass * mLearnShellBody.mL  / mLearnShellBody.mvalIx0 ;
	arrF[4] = -val_mn2 * marrStrSK_VS [4] ;
   ///

   //arrF[5] = -G_ZEMLI *  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) / valTay /ATM_R_UNIVER;
   arrF[5] = -  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) / valTay*(G_ZEMLI /ATM_R_UNIVER + valDerivTay);//HT[0][2]);

  mtr_dF_po_di[6 * 2 ] =  -  valCxEtal * val_q * valSm / mLearnShellBody.mMass ;

    double valTemp = val_ix * valCxEtal * val_q * marrStrSK_VS [6]*marrStrSK_VS [6] ;
   arrF[6] = -G_ZEMLI * sin(marrStrSK_VS [7])   -  val_ix * valCxEtal * val_q * marrStrSK_VS [6]*marrStrSK_VS [6] ;
   arrF[7] = -G_ZEMLI * cos(marrStrSK_VS [7])/ marrStrSK_VS [6]  +
			marrStrSK_VS [6] * cos(marrStrSK_VS [7]) / (R_ZEMLI + marrStrSK_VS [1]) +
			valDeltaTettaTochka;

	arrF[3] =  (mLearnShellBody.mvalIx0 * mLearnShellBody.mL / (mLearnShellBody.mMass * mLearnShellBody.mDm * mLearnShellBody.mh_gob))
		  * val_iz *valKnm * marrStrSK_VS [4]* arrF[7]/marrStrSK_VS [6]/ cos(marrStrSK_VS [7])+ valDeltaPsiTochka; ;
		//  *(-G_ZEMLI/ marrStrSK_VS [6]/ marrStrSK_VS [6] + 1./(R_ZEMLI + marrStrSK_VS [1]))+ valDeltaPsiTochka;
  // mtr_dF_po_di[3 * 2 + 1] = ( arrF[3] - valDeltaPsiTochka)/ val_iz ;

 */

   /*
// ???????  ??????? ??????? ??????????? ?????? ????? ?? ????????? ????????   - mtr_dF_po_dz
   memset(mtr_dF_po_dx, 0, 64 * sizeof( double))  ;

    double arrGradF3_po_z [8] = {0} , arrGradF6_po_z [8] = {0.} ;
   // arr_dFTransp - ??????????????? ?????? ??? ???????? ?????????? ??????? Fi. ????????? ???????? ?????????.
   // ???????, ??? ????????? ???????  mtr_dF_po_dx  ?????? mtr_dF_po_dx ???? ???????????????
    double	arr_dFTransp[64]= {0.} ;
   fncCalcGradF0(  arr_dFTransp    );       // ?????? ????????? ??????? F0
   fncCalcGradF1( &arr_dFTransp[ 8]); // ?????? ????????? ??????? F1
   fncCalcGradF2( &arr_dFTransp[16]); // ?????? ????????? ??????? F2
   fncCalcGradF7( valDerivDeltaTettaTochkaPoPsi ,&arr_dFTransp[7 * 8]) ;// ?????? ????????? ??????? F5

   fncCalcGradF3(&arr_dFTransp[7 * 8],arrF[7], valDeltaPsiTochka      // ?????? ????????? ??????? F3
		   ,valDerivDeltaPsiTochkaPoPsi, valDerivDeltaPsiTochkaPoTetta
		   ,val_iz, val_Deriv_iz
		   ,valKnm, arrGradKnm
		   ,&arr_dFTransp[ 3 * 8],arrGradF3_po_z) ;

   fncCalcGradF4(valMach, arrGradMach    // ?????? ????????? ??????? F4
				 ,val_q, arrGrad_q
				 ,valMxOmegax, arrGradMxOmegax
				 ,&arr_dFTransp[ 4 * 8]) ;

   fncCalcGradF5( valTay, valDerivTay, &arr_dFTransp[ 5 * 8]) ;



   fncCalcGradF6(val_q, arrGrad_q
							   ,valCxEtal, arrGradCxEtal
							   , val_ix,  val_Deriv_ix
							   ,&arr_dFTransp[ 6 * 8],arrGradF6_po_z); // ?????? ????????? ??????? F6
  ///


 //MatrTransp(arr_dFTransp, 8, 8, mtr_dF_po_dx);
 memcpy( mtr_dF_po_dx, arr_dFTransp, 64 * sizeof(double)) ;
 // ??? ??????? !!! ?? ??????
 mtr_dF_po_dz[ 3 * 8 +7 ] = arrGradF3_po_z[7] ;
 mtr_dF_po_dz[ 6 * 8 +7 ] = arrGradF6_po_z[7] ;
  */


}

void TLearnShellTraj::fncCalc_F(TEnvironment Environment, double *arrF)
{
  // ?????????? ?????????? ??????????? ???????????  ? ?? ?????????
   double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///



   // ?????????? ??????? ???????? ????? ? ???
   double arrWindV_SSK[3] ={0.};
   calcVectWindV(Environment, arrWindV_SSK);
   ///

   // ?????????? ????????? ??????????
   //const double VAlVozdV = sqrt(marrStrSK_VS[6] * marrStrSK_VS[6] - 2. * arrWindV_SSK[0] * marrStrSK_VS[6]
  //	+ arrWindV_SSK[1]  * arrWindV_SSK[1]  + Environment.mWind_VertV * Environment.mWind_VertV);

	 const double VAlVozdV = sqrt(marrStrSK_VS[6] * marrStrSK_VS[6] - 2. * arrWindV_SSK[0] * marrStrSK_VS[6]
	+ arrWindV_SSK[1]  * arrWindV_SSK[1]  + arrWindV_SSK[2]  * arrWindV_SSK[2]);
	///

  // ?????????? ????? ???? ? ??????? ??? ??????? ???????????
   double valMach = calcMach(  valTay,  VAlVozdV);

  ///


	///
	double valCx = -1., val_Deriv_Cx = 0. ;
	mLearnShellBody.fnkCxEtal (valMach,  valCx, val_Deriv_Cx) ;


  /*
  // ??? ??? ? ?????
  // ?????????? ????????????? deltaCx,deltaCy, deltaCz
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

 // ?????????? ??????? ???????? ????? ? ???
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
   double epsW = arcSin( arrAirV_SSK[2]/ Norm3( arrAirV_SSK) );
   double betW = arcSin( arrAirV_SSK[1]/ sqrt(arrAirV_SSK[1] *  arrAirV_SSK[1] + arrAirV_SSK[0] *  arrAirV_SSK[0] ));
   double gamW = arcSin(arrAirV_SSK[2]/ sqrt(arrAirV_SSK[2] *  arrAirV_SSK[2] + arrAirV_SSK[0] *  arrAirV_SSK[0] ));
	double deltaCx  = valCx * (cos(epsW) * cos(betW) -1. );
  double deltaCy =  valCx * sin(betW) ;
 // double deltaCz =  valCx * sin(epsW);
	double valCz =   mLearnShellBody.fnkCz (valMach);
 double deltaCz =  valCz *   sin(gamW);
 ////////

  // ?????????? ?????? q

	 double val_q = calc_q(VAlVozdV, valTay) ;
  ///

  // ??????? ??????????? ????????
   double valVpr = marrStrSK_VS[6] * cos(marrStrSK_VS [7]) * (R_ZEMLI + mAltit)/
	   (R_ZEMLI + marrStrSK_VS [1])  ;
	   ///

  // ?????????? ?????????????? ??????? Knm ? ??????? ?? ??????????? ?? x
  // ?? ??????? 1-?? ? 6-?? ?????????? ( ???????? ?????????? ? 0)
   //	 double valKnm = 0., valDerivKnm = 0. ;

	//  mLearnShellBody.fnkKnm(valMach, valKnm, valDerivKnm);
	///

 // ?????????? ??????????????? ?????? iz(z2) ? ?? ???????????
   //	 double val_iz = 0., val_Deriv_iz = 0. ;
 //	mLearnShellBody.fnkIz0(mTet0,val_iz, val_Deriv_iz) ;
	///
	// ?????????? ??????????????? ?????? iz(z2) ? ?? ???????????
	 double valKnm = 0., val_Deriv_Knm = 0. ;
	 mLearnShellBody.fnkKnm (valMach, valKnm, val_Deriv_Knm )  ;


// ?????????? ?????????????? ??????? mxwx  ? ?? ??????????? ?? M

	double valMxOmegax = 0., valDerivMxOmegax = 0;
	 mLearnShellBody.fnkMxOmegax(valMach,valMxOmegax, valDerivMxOmegax ) ;
	///

  // ??????? ??????
	const double VAlSmid = mLearnShellBody.mDm * mLearnShellBody.mDm / 4. * M_PI;

// ?????????? ?????????? ?????????? ???????????? ?????? ????? ??? ????? ???????? ?????
    double valDeltaTettaTochka =0.          //???????  ???????? ???? ??????? ?????? DeltaTettaTochka
			  ,valDerivDeltaTettaTochkaPoPsi= 0. // ??????????? DeltaTettaTochka ?? Psi (= marrStrSK_VS [3])
			  ,valDeltaPsiTochka =0.           //???????  ???????? ???? ???? DeltaPsiTochka
			  ,valDerivDeltaPsiTochkaPoPsi= 0. // ??????????? DeltaPsiTochka ?? Psi (= marrStrSK_VS [3])
			  ,valDerivDeltaPsiTochkaPoTetta= 0.;  // ??????????? DeltaPsiTochka ?? Tetta (=  marrStrSK_VS [7])
   fncCalcDeltaTettaTochka_and_DerivPoPsi(valDeltaTettaTochka, valDerivDeltaTettaTochkaPoPsi) ;
   fncCalcDeltaPsiTochka_and_DerivPoPsi_and_DerivPoTetta(valDeltaPsiTochka
		   ,valDerivDeltaPsiTochkaPoPsi, valDerivDeltaPsiTochkaPoTetta);
	 ///

   arrF[0] =  valVpr * cos(marrStrSK_VS [3]) ;
   arrF[1] =  marrStrSK_VS [6] * sin(marrStrSK_VS [7]) ;
   arrF[2] =  -valVpr * sin(marrStrSK_VS [3]) ;

	// const  double val_qsbm = 0.474 * mLearnShellBody.mDm* mLearnShellBody.mDm/ mLearnShellBody.mMass * marrStrSK_VS [6] * marrStrSK_VS [6] * marrStrSK_VS[5];
 //	 const  double val_mmx = (valMxOmegax * mLearnShellBody.mDm /2./ mLearnShellBody.mV0);
	// const  double val_mn2  = val_mmx* val_qsbm  * mLearnShellBody.mMass * mLearnShellBody.mL  / mLearnShellBody.mvalIx0 ;
	arrF[4] = -valMxOmegax *mLearnShellBody.mL *mLearnShellBody.mL/ mLearnShellBody.mvalIx0 * val_q
			 * VAlSmid /valMach /ATM_AN0 * marrStrSK_VS [4];
	// -val_mn2 * marrStrSK_VS [4] ;
	 ///

   arrF[5] = -  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) / valTay*(G_ZEMLI /ATM_R_UNIVER + valDerivTay);//HT[0][2]);



  //  double valTemp = val_ix * valCxEtal * val_q * marrStrSK_VS [6]*marrStrSK_VS [6] ;
   arrF[6] = -G_ZEMLI * sin(marrStrSK_VS [7])   -   (valCx + deltaCx) * VAlSmid/ mLearnShellBody.mMass * val_q ;
   arrF[7] = -G_ZEMLI * cos(marrStrSK_VS [7])/ marrStrSK_VS [6]  +
			marrStrSK_VS [6] * cos(marrStrSK_VS [7]) / (R_ZEMLI + marrStrSK_VS [1])
			-  deltaCy * VAlSmid/ mLearnShellBody.mMass * val_q / marrStrSK_VS [6]
			+ valDeltaTettaTochka;

	arrF[3] =  (mLearnShellBody.mvalIx0 * mLearnShellBody.mL / (mLearnShellBody.mMass * mLearnShellBody.mDm * mLearnShellBody.mh_gob))
		  * valKnm * marrStrSK_VS [4]* arrF[7]/marrStrSK_VS [6]/ cos(marrStrSK_VS [7])+ valDeltaPsiTochka
		  - deltaCz * VAlSmid/ mLearnShellBody.mMass * val_q / cos(marrStrSK_VS [7])/marrStrSK_VS [6] ;



}

// ?????????? ??????? ???????? ????? ? ???

void TLearnShellTraj::calcVectWindV(TEnvironment Environment, double *arrWindV_SSK)
{
 arrWindV_SSK[0] = -Environment.mWind_V * cos(Environment.mWind_Alf - (mAlfDir - marrStrSK_VS [3]))
	 * cos(marrStrSK_VS [7] ) + Environment.mWind_VertV * sin (marrStrSK_VS [7]);
 arrWindV_SSK[1] =  Environment.mWind_V * cos(Environment.mWind_Alf - (mAlfDir - marrStrSK_VS [3]))
	 * sin(marrStrSK_VS [7] ) + Environment.mWind_VertV * cos (marrStrSK_VS [7]);
 arrWindV_SSK[2] = -Environment.mWind_V * sin(Environment.mWind_Alf - (mAlfDir - marrStrSK_VS [3])) ;
}

// ?????????? ????? ???? ? ??????? ?????? ??????????? ????? ???? ?? x
// INPUT:
//valTay - ?????????? ???????????  ???????????
// valDerivTay - ??????????? ?? ?????? ?????????? ???????????  ???????????
// OUTPUT
// valMach - ????? ????
// arrGradMach - ?????? ????????? ????? ???? ?? ??????? ??????????
// ???? ?????????? !!!!!
void TLearnShellTraj::fncCalcMach_and_GradMach(const double valTay, const double valDerivTay
	 , double &valMach,  double *arrGradMach)
{
  /* memset( arrGradMach, 0, 8 * sizeof( double));
   arrGradMach [6] = sqrt(ATM_TAYN0/ valTay)/ ATM_AN0 ;
   valMach = marrStrSK_VS[6] * arrGradMach [6] ;
   arrGradMach [1] = - 0.5 * valMach * valDerivTay/ valTay ; */

}

// ?????????? ????? ???? ? ??????? ?????? ??????????? ????? ???? ?? x
// INPUT:
//valTay - ?????????? ???????????  ???????????
// valVVozd - ????????? ????????

double TLearnShellTraj::calcMach( double valTay,  double valVVozd)
{
 return valVVozd *sqrt(ATM_TAYN0/ valTay)/ ATM_AN0 ;

}

// ?????????? Knm(M(x7,x2))?  ??????? ??????? ??????????? ??????? Knm(M(x7,x2)) ?? x
// INPUT:
//valMach - ????? ????
// arrGradMach - ???????? ????? ???? ?? x
// OUTPUT
// valKnm - ??????? Knm
// arrGradKnm - ?????? ????????? ??????? Knm ?? ??????? ??????????
void TLearnShellTraj::fncCalcKnm_and_Grad_Knm(const  double valMach,  double *arrGradMach
			   , double &valKnm,  double *arrGradKnm )
{
  /* memset( arrGradKnm, 0, 8 * sizeof( double));
	double  valDerivKnm = 0;
   mLearnShellBody.fnkKnm(valMach, valKnm, valDerivKnm);
   arrGradKnm [1] = valDerivKnm * arrGradMach [1] ;
   arrGradKnm [6] = valDerivKnm * arrGradMach [6] ;
  */
}






// ?????????? MxOmegax(M(x7,x2))?  ??????? ??????? ??????????? ??????? MxOmegax(M(x7,x2)) ?? x
// INPUT:
//valMach - ????? ????
// arrGradMach - ???????? ????? ???? ?? x
// OUTPUT
// valMxOmegax - ??????? valMxOmegax
// arrGradMxOmegax - ?????? ????????? ??????? MxOmegax ?? ??????? ??????????
void TLearnShellTraj::fncCalcMxOmegax_and_Grad_MxOmegax(const  double valMach,  double *arrGradMach
			   , double &valMxOmegax,  double *arrGradMxOmegax )
{
  /* memset( arrGradMxOmegax, 0, 8 * sizeof( double));
	double  valDerivMxOmegax = 0;
   mLearnShellBody.fnkMxOmegax(valMach,valMxOmegax, valDerivMxOmegax ) ;
   arrGradMxOmegax [1] = valDerivMxOmegax * arrGradMach [1] ;
   arrGradMxOmegax [6] = valDerivMxOmegax * arrGradMach [6] ;
   */

}
void TLearnShellTraj::fncCalcCxEtal_and_Grad_CxEtal(const  double valMach,  double *arrGradMach
			   , double &valCxEtal,  double *arrGradCxEtal )
{
   memset( arrGradCxEtal, 0, 8 * sizeof( double));
	double  valDerivCxEtal = 0;
   mLearnShellBody.fnkCxEtal(valMach,valCxEtal, valDerivCxEtal ) ;
   arrGradCxEtal [1] = valDerivCxEtal * arrGradMach [1] ;
   arrGradCxEtal [6] = valDerivCxEtal * arrGradMach [6] ;

}



// ??????????? ??????????? ?????? ? ??????? ??????? ??????????? ???????? ?????? ?? x
// INPUT:
// valMach - ????? ????
// arrGradMach[8] - ?????? ????????? ????? ???? ?? ??????? ??????????
// OUTPUT :
// val_q - ?????????? ?????
// arrGrad_q[8]  -  ?????? ????????? c?????????? ??????  ?? ??????? ??????????
// ?????????? !!!!!
void TLearnShellTraj::fncCalc_q_and_Grad_q( double valMach, double  *arrGradMach
		, double  &val_q,  double *arrGrad_q)
{

/*   double temp = ATM_RoN0 * ATM_AN0 *ATM_AN0 * valMach;
 //val_q =  temp * marrStrSK_VS [5] * valMach / 2. ;
  memset( arrGrad_q, 0, 8 * sizeof( double));
  arrGrad_q [1] = temp *marrStrSK_VS [5] *arrGradMach [1] ;
  arrGrad_q [5] = temp * valMach / 2. ;
  arrGrad_q [6] =  temp *marrStrSK_VS [5] *arrGradMach [6] ;
 val_q =  mLearnShellBody.mDm* mLearnShellBody.mDm * 0.474 / mLearnShellBody.mMass * marrStrSK_VS [5] ;
*/
}



// ??????????? ??????????? ?????? ? ??????? ??????? ??????????? ???????? ?????? ?? x
// INPUT:
// valMach - ????? ????
// arrGradMach[8] - ?????? ????????? ????? ???? ?? ??????? ??????????
// OUTPUT :
// val_q - ?????????? ?????
// arrGrad_q[8]  -  ?????? ????????? c?????????? ??????  ?? ??????? ??????????
// ?????????? !!!!!
double TLearnShellTraj::calc_q(const  double VAlVVozd, const  double VAlTay)
{
return ATM_RoN0 / 2. * ATM_TAYN0 /VAlTay * marrStrSK_VS [5] *  VAlVVozd * VAlVVozd;
}

// ?????? ?????????? ?? ????????? ????? ??? ????? ???????? ???????? ?????
 void TLearnShellTraj::fncCalcDeltaTettaTochka_and_DerivPoPsi( double &valDeltaTettaTochka
		   ,  double &valDerivDeltaTettaTochkaPoPsi)
 {
	  double temp = 2.* OMEGA_ZEMLI * cos (mLatitude);
	 valDeltaTettaTochka  =  temp * sin(mAlfDir - marrStrSK_VS [3]) ;
	 valDerivDeltaTettaTochkaPoPsi = - temp * cos(mAlfDir - marrStrSK_VS [3]) ;
	   //??? ??????? !!!
   //	  valDeltaTettaTochka = 0. ;
   //	 valDerivDeltaTettaTochkaPoPsi = 0. ;
 }

 // ?????? ?????????? ?? ????????? ????? ??? ????? ???????? ???????? ?????
 void TLearnShellTraj::fncCalcDeltaPsiTochka_and_DerivPoPsi_and_DerivPoTetta( double &valDeltaPsiTochka
		   ,  double &valDerivDeltaPsiTochkaPoPsi,  double &valDerivDeltaPsiTochkaPoTetta)
 {
	 valDeltaPsiTochka =  -2.* OMEGA_ZEMLI *
	 ( sin (mLatitude) - cos (mLatitude) * cos(mAlfDir - marrStrSK_VS [3])* tan(marrStrSK_VS [7])) ;
	 valDerivDeltaPsiTochkaPoPsi = 2.* OMEGA_ZEMLI * cos (mLatitude)*sin(mAlfDir - marrStrSK_VS [3])*tan(marrStrSK_VS [7]) ;
	 valDerivDeltaPsiTochkaPoTetta = 2.* OMEGA_ZEMLI * cos (mLatitude)*cos(mAlfDir - marrStrSK_VS [3])
	   /cos(marrStrSK_VS [7])/cos(marrStrSK_VS [7]) ;
	  //??? ??????? !!!
	 // valDeltaPsiTochka = 0. ;
	 //valDerivDeltaPsiTochkaPoPsi = 0. ;
	// valDerivDeltaPsiTochkaPoTetta = 0. ;
 }

 // ?????????? ??????? ??????? ??????????? (?????????) ??????? F0
 // OUTPUT:
 // arrGradF0 [8] - ?????? ???????????
 void TLearnShellTraj:: fncCalcGradF0(  double *arrGradF0)
 {
	memset( arrGradF0,0, 8 * sizeof( double)) ;
	 double valTemp0 = (mAltit + R_ZEMLI) /(marrStrSK_VS [1] + R_ZEMLI); // ?????????? ??????????
	 double valTemp1 = valTemp0 * cos(marrStrSK_VS [7]); // ?????????? ??????????
	 double valTemp2 = cos(marrStrSK_VS [3] ); // ?????????? ??????????
	arrGradF0[1] = - marrStrSK_VS [6]* valTemp1 * valTemp2 /(marrStrSK_VS [1] + R_ZEMLI);
	arrGradF0[3] =  - sin( marrStrSK_VS [3] )* marrStrSK_VS [6] *valTemp1 ;
	arrGradF0[6] = valTemp2 * valTemp1 ;
	arrGradF0[7] = - marrStrSK_VS [6]* sin(marrStrSK_VS [7])* valTemp0 * valTemp2;
 }

  // ?????????? ??????? ??????? ??????????? (?????????) ??????? F1
 // OUTPUT:
 // arrGradF1 [8] - ?????? ???????????
 void TLearnShellTraj::fncCalcGradF1(  double *arrGradF1)
 {
	memset( arrGradF1,0, 8 * sizeof( double)) ;
	arrGradF1[6] = sin(marrStrSK_VS [7]) ;
	arrGradF1[7] =  marrStrSK_VS [6]* cos(marrStrSK_VS [7]);
 }

 // ?????????? ??????? ??????? ??????????? (?????????) ??????? F2
 // INPUT:
 // valDerivDeltaTettaTochkaPoPsi - ??????????? ??????????? ?????(???????? ?????) ?? ???
 //
 // OUTPUT:
 // arrGradF2 [8] - ?????? ???????????
 void TLearnShellTraj::fncCalcGradF2(  double *arrGradF2)
 {
	memset( arrGradF2,0, 8 * sizeof( double)) ;
	 double valTemp0 = (mAltit + R_ZEMLI) /(marrStrSK_VS [1] + R_ZEMLI); // ?????????? ??????????
	 double valTemp1 = valTemp0 * cos(marrStrSK_VS [7]); // ?????????? ??????????
	 double valTemp2 = sin(marrStrSK_VS [3] ); // ?????????? ??????????
	arrGradF2[1] =  marrStrSK_VS [6]* valTemp1 * valTemp2 /(marrStrSK_VS [1] + R_ZEMLI);
	arrGradF2[3] =  - cos( marrStrSK_VS [3] )* marrStrSK_VS [6] *valTemp1 ;
	arrGradF2[6] =  - valTemp2 * valTemp1 ;
	arrGradF2[7] =  marrStrSK_VS [6]* sin(marrStrSK_VS [7])* valTemp0 * valTemp2;
 }

  // ?????????? ??????? ??????? ??????????? (?????????) ??????? F7(x1,x4,x7,x8)- ???? TETTA
 // OUTPUT:
 // arrGradF3 [8] - ?????? ???????????
 void TLearnShellTraj::fncCalcGradF7(const  double valDerivDeltaTettaTochkaPoPsi , double *arrGradF7)
 {
	memset( arrGradF7, 0, 8 * sizeof( double)) ;
	 double valTemp0 = 1. /(marrStrSK_VS [1] + R_ZEMLI); // ?????????? ??????????
	 double valTemp1 = valTemp0 * cos(marrStrSK_VS [7]); // ?????????? ??????????

	arrGradF7[1] =  -marrStrSK_VS [6]* valTemp1 * valTemp0 ;
	arrGradF7[6] =  G_ZEMLI * cos(marrStrSK_VS [7])/marrStrSK_VS [6]/marrStrSK_VS [6] + valTemp1 ;
	arrGradF7[7] =  sin(marrStrSK_VS [7])*( G_ZEMLI /marrStrSK_VS [6] - marrStrSK_VS [6] *valTemp0 ) ;
	arrGradF7[3] = valDerivDeltaTettaTochkaPoPsi ;
 }
  // ?????????? ??????? ??????? ??????????? (?????????) ??????? F3  - ???
  // INPUT:
  // arrGradF7 -  ???????? ??????? F7
  // valF7 - ???????? ??????? F7
  // valDeltaPsiTochka -  ??????????? ?????? ???? ??????? ??????????
  // valDerivDeltaPsiTochkaPoPsi - ??????? ??????????? ??????????? ?????? ???? ??????? ?????????? ?? Psi ( ?? marrStrSK_VS [3]) )
  // valDerivDeltaPsiTochkaPoTetta - ??????? ??????????? ??????????? ?????? ???? ??????? ?????????? ?? ?etta ( ?? marrStrSK_VS [7]) )
 // OUTPUT:
 // arrGradF3 [8] - ?????? ??????????? ?? x
 // arrGradF3_po_z - ?????? ????? ??????????? ?? ????????? ???????? , ?? ??????? ??????? ?? z2 (x8)
 void TLearnShellTraj::fncCalcGradF3( double *arrGradF7,const  double valF7,const  double valDeltaPsiTochka
		   ,const  double valDerivDeltaPsiTochkaPoPsi, const  double valDerivDeltaPsiTochkaPoTetta
		   ,const  double val_iz, const  double val_Deriv_iz
		   ,const  double valKnm,  double *arrGradKnm
		   ,  double *arrGradF3, double *arrGradF3_po_z)
 {
	memset( arrGradF3,0, 8 * sizeof( double)) ;
	memset( arrGradF3_po_z,0, 8 * sizeof( double)) ;
	// ????????? ?????????? ?????????
	const  double val_k0 = mLearnShellBody.mvalIx0 * mLearnShellBody.mL /(mLearnShellBody.mMass * mLearnShellBody.mDm * mLearnShellBody.mh_gob);
	///
	const  double  val_cosx8 = cos( marrStrSK_VS [7]) ;  // ????????? ?????
	const  double  val_x7cosx8 = val_cosx8 *  marrStrSK_VS [6] ; // ????????? ?????

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
// ?????????? ??????? ??????? ??????????? (?????????) ??????? F4  - ?????(??????? ???????? ???????? )
  // INPUT:
  // valMach, arrGradMach  - ????? ???? ? ??? ????????
  // val_q, arrGrad_q  - ?????????? ????? ? ??? ????????
  // valMxOmegax, arrGrad_MxOmegax -  ????? ????????????? ??????? ? ??? ????????
 // OUTPUT:
 // arrGradF4 [8] - ?????? ??????????? ?? x
 void TLearnShellTraj::fncCalcGradF4(const  double valMach,  double *arrGradMach
							   ,const  double val_q,  double *arrGrad_q
							   ,const  double valMxOmegax,  double *arrGrad_MxOmegax
							   ,  double *arrGradF4)


 {
	memset( arrGradF4,0, 8 * sizeof( double)) ;
	const  double valSm = M_PI * mLearnShellBody.mDm * mLearnShellBody.mDm / 4. ;
	const  double val_k0  = -valSm *  mLearnShellBody.mL* mLearnShellBody.mL / ATM_AN0 / 2. ;
	arrGradF4 [1] =  val_k0 * marrStrSK_VS [4] * ( (arrGrad_MxOmegax [1] * val_q  + valMxOmegax * arrGrad_q[1] ) * valMach
				- valMxOmegax * val_q * arrGradMach [1]) / (valMach * valMach) ;
	arrGradF4 [4] =  val_k0  * valMxOmegax * val_q / valMach ;
	arrGradF4 [5] =  val_k0  * valMxOmegax * marrStrSK_VS [4] * arrGrad_q [5] / valMach ;
	arrGradF4 [6] =   val_k0 * marrStrSK_VS [4] * ( (arrGrad_MxOmegax [6] * val_q  + valMxOmegax * arrGrad_q[6] ) * valMach
				- valMxOmegax * val_q * arrGradMach [6]) / (valMach * valMach) ;

 }


 // ?????????? ??????? ??????? ??????????? (?????????) ??????? F5  - pi (??????? ????????? ???????)
  // INPUT:
  // valTay, valDerivTay  - ??????? ??????????? ??????????? ? ?? ??????????? ?? ??????

 // OUTPUT:
 // arrGradF5 [8] - ?????? ??????????? ?? x
void TLearnShellTraj::fncCalcGradF5(const  double valTay, const  double valDerivTay,  double *arrGradF5)
 {
	memset( arrGradF5,0, 8 * sizeof( double)) ;
	arrGradF5 [1] =  G_ZEMLI  *  marrStrSK_VS [5] * marrStrSK_VS [6] * sin(marrStrSK_VS [7]) *valDerivTay
		 / valTay / valTay /ATM_R_UNIVER ;
	arrGradF5 [5] =  -G_ZEMLI  *   marrStrSK_VS [6] * sin(marrStrSK_VS [7])/ valTay /ATM_R_UNIVER;
	arrGradF5 [6] =  -G_ZEMLI * marrStrSK_VS [5] *  sin(marrStrSK_VS [7]) / valTay /ATM_R_UNIVER;
	arrGradF5 [7] =  -G_ZEMLI * marrStrSK_VS [5] * marrStrSK_VS [6] * cos(marrStrSK_VS [7]) / valTay /ATM_R_UNIVER;
 }



 // ?????????? ??????? ??????? ??????????? (?????????) ??????? F6  - Vk
  // INPUT:
  // val_q, arrGrad_q  - ?????????? ????? ? ??? ????????
  // valCxEtal, arrGradCxEtal -  ????? ?????????? ???? ????????????? ??????? ?????????
  // val_ix, val_Deriv_ix - ???????????? ????? ? ?? ??????????? ?? ????? ????
 // OUTPUT:
 // arrGradF6 [8] - ?????? ??????????? ?? x
 // arrGradF6_po_z - ?????? ????? ??????????? ?? ????????? ???????? , ?? ??????? ??????? ?? z2 (x8)
 void TLearnShellTraj::fncCalcGradF6(const  double val_q,  double *arrGrad_q
							   ,const  double valCxEtal,  double *arrGradCxEtal
							   ,const  double val_ix, const  double val_Deriv_ix
							   , double *arrGradF6, double *arrGradF6_po_z)
 {
	memset( arrGradF6,0, 8 * sizeof( double)) ;
	memset( arrGradF6_po_z,0, 8 * sizeof( double)) ;
	const  double valSm = M_PI * mLearnShellBody.mDm * mLearnShellBody.mDm / 4. ;
	// ????????? ?????????? ?????????
   const  double val_k0 = - val_ix *valSm / mLearnShellBody.mMass ;

	///


	arrGradF6 [1] = val_k0 * ( arrGradCxEtal[1] * val_q + valCxEtal * arrGrad_q [1]) ;



	arrGradF6 [6] = val_k0 * ( arrGradCxEtal[6] * val_q + valCxEtal * arrGrad_q [6]) ;


	arrGradF6 [7] = -G_ZEMLI ;

	arrGradF6_po_z [7] = -valSm / mLearnShellBody.mMass * valCxEtal * val_q * val_Deriv_ix ;

 }

 // ??? ????? ?????? ?? ???? valStepInt
 void TLearnShellTraj::fncEilerStep(TEnvironment Environment,const  double valStepInt)
 {
    double arrF[8] ={0.} ;
   // ??? ?????????????? ???????? ???????
   fncCalc_F(Environment,arrF) ;
	double arrT0[64] = {0.},arrT1[64] = {0.};
   MatrxMultScalar(arrF, 8, 1, valStepInt,arrT0); // f * dt
   MtrxSumMatrx(arrT0, marrStrSK_VS,8, 1, arrT1) ;
   memcpy(marrStrSK_VS, arrT1, 8 * sizeof( double)) ;   ///

  mTCur += valStepInt ;


 }

// ????????????? ????????? ?? ???????  valTNext
void TLearnShellTraj::fncMovePhasVector(TEnvironment Environment,const double VAlStepInt, const double valTNext )
{
  int iCirc = (valTNext - mTCur) / VAlStepInt ;
   double valTTemp = mTCur + (( double)iCirc) * VAlStepInt ;
  for (int i = 0; i < iCirc; i++)
  {
	  fncEilerStep(Environment,VAlStepInt);
  }

  fncEilerStep(Environment,valTNext - valTTemp);

}

 /*
// ????????????? ????????? ?? ???????  ???????
// OUTPUT:
// valDHoriz - ??????????????? ????? ????? ???????
void TLearnShellTraj::fncMoveShell_TO_ZeroAlt(TEnvironment Environment,const double VAlStepInt
, double &valDHoriz )
{
  int iCirc = 1000. / VAlStepInt ;

  int i = 0 ;
  for ( i = 0; i < iCirc; i++)
  {
	  fncEilerStep(Environment,VAlStepInt);
	  if (marrStrSK_VS [1] < 0.) break ;


  }
   valDHoriz = sqrt(marrStrSK_VS [0]*marrStrSK_VS [0]+ marrStrSK_VS [2]*marrStrSK_VS [2]) ;


}
*/
// ????????????? ????????? ?? ???????  ??????? ?? ????? ? ?????? ???????? ?????
// OUTPUT:
// valGeoDist - ????????????? ??????????????? ????? ????? ???????
// valGeoTetta - ????????????? ???? ??????? ?????????? ? ????? ???????
void TLearnShellTraj::fncMoveShell_TO_ZeroAlt(TEnvironment Environment,const double VAlStepInt
, double &valGeoDist )
{
  int iCirc = 1000. / VAlStepInt ;

  int i = 0 ;
  for ( i = 0; i < iCirc; i++)
  {
	  fncEilerStep(Environment,VAlStepInt);

	  // ?????????? ?? ?????? ?????
	  double val_r = sqrt((marrStrSK_VS [1] + R_ZEMLI)* (marrStrSK_VS [1] + R_ZEMLI)
						  + marrStrSK_VS [0]* marrStrSK_VS [0]);
	  ///

	  if (val_r <= R_ZEMLI)
	  {
		break;
	  }

  }
	  double valGeoAlf = atan(marrStrSK_VS [0]/ (marrStrSK_VS [1] + R_ZEMLI));
	  valGeoDist =  R_ZEMLI * valGeoAlf;


}


// ?????????? ??????? ????????? ? ??????? ??????????? ??
// OUTPUT:
// arrVS_PrStSK [6]   - ?????? ????????? ? ??????? ??????????? ??
void TLearnShellTraj::fncCalcVS_v_PrStSK(   double *arrVS_PrStSK  )
{
	arrVS_PrStSK[0] = marrStrSK_VS [0] ;
	arrVS_PrStSK[1] = marrStrSK_VS [1] ;
	arrVS_PrStSK[2] = marrStrSK_VS [2] ;
	arrVS_PrStSK[3] = marrStrSK_VS [6] * cos(marrStrSK_VS [7])* cos(marrStrSK_VS [3]);
	arrVS_PrStSK[4] = marrStrSK_VS [6] * sin(marrStrSK_VS [7]) ;
	arrVS_PrStSK[5] = - marrStrSK_VS [6] * cos(marrStrSK_VS [7])* sin(marrStrSK_VS [3]);
}


 // ?????????? ??????? ??????? ??????????? ??? ???????? ??? ??????? ??
 // ??????????? ?? ? ????????????? ??????????? ??
 // OUTPUT:
 // arrJac_PrStSK [6*8]
void TLearnShellTraj::fncCalcJacobi_PrStSK(  double *arrJac_PrStSK  )
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


// ????????????? ????????? ?? ???????  ???????
// OUTPUT:
// valDHoriz - ??????????????? ????? ????? ???????
void TLearnShellTraj::fncMoveClass_TO_ZeroAlt_AND_ShowGraphs(TEnvironment Environment
,const double VAlStepInt, wchar_t *wcharrPath1,  double &valGeoDist )
{
  const int QUANT_COLS = 9 , QUANT_POINTS_MAX = 300000;
  const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ???????????? ????? ????? ??????????
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

	  /*
		 double val_r = sqrt((marrStrSK_VS [1] + R_ZEMLI)* (marrStrSK_VS [1] + R_ZEMLI)
						  + marrStrSK_VS [0]* marrStrSK_VS [0]);
	  ///

	  if (val_r <= R_ZEMLI)
	  {
		break;
	  }

  }
	  double valGeoAlf = atan(marrStrSK_VS [0]/ (marrStrSK_VS [1] + R_ZEMLI));
	  valGeoDist =  R_ZEMLI * valGeoAlf;

	  */
	  double val_r = sqrt((marrStrSK_VS [1] + R_ZEMLI)* (marrStrSK_VS [1] + R_ZEMLI)
						  + marrStrSK_VS [0]* marrStrSK_VS [0]);
	  ///

	  if (val_r <= R_ZEMLI)
	  {
		break;
	  }

	  //if (((double)mTCur > (valTOut + DEL_T - 1E-15 ))&& ( iNupPointsOut < QUANT_POINTS_MAX ))
	  if ( iNupPointsOut < QUANT_POINTS_MAX )
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
   double valGeoAlf = atan(marrStrSK_VS [0]/ (marrStrSK_VS [1] + R_ZEMLI));
	  valGeoDist =  R_ZEMLI * valGeoAlf;

	 if (wcharrPath)
	 {
	 for (int j = 1; j < QUANT_COLS -1; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,0  // ????? ?????????? ?? ??? X
								  ,j  // ????? ?????????? ?? ??? Y
								  ,100 //  ??????? ?? ??? X
								  ,pscaleY[j]  // ??????? ?? ??? Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,2  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,3  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete parrBuff ;
  delete pwcharrFileNames ;
  delete pscaleY ;

   }

	// fncCalcMtrxPartialDeriv( Environment, VAlStepInt );
	// fncCalcVectPartialDeriv_CoeffForm(Environment, VAlStepInt, marrDelta[8] ) ;
	// fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt,marrDelta[9] * mLearnShellBody.mMass ) ;

}



// ????????????? ????????? ?? ???????  ???????
// OUTPUT:
// valDHoriz - ??????????????? ????? ????? ???????
void TLearnShellTraj::fncMoveClass_TO_ZeroAlt_AND_ShowGraphs(TEnvironment Environment
,const double VAlStepInt, wchar_t *wcharrPath1,  double &valGeoDist,  double &valH )
{
  const int QUANT_COLS = 9 , QUANT_POINTS_MAX = 300000;
  const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ???????????? ????? ????? ??????????
   wchar_t 	wcharrPath[400] = {0};;
   valH = -1.;

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


	  double val_r = sqrt((marrStrSK_VS [1] + R_ZEMLI)* (marrStrSK_VS [1] + R_ZEMLI)
						  + marrStrSK_VS [0]* marrStrSK_VS [0]);
	  ///

	  if (val_r <= R_ZEMLI)
	  {
		break;
	  }


	  if ( iNupPointsOut < QUANT_POINTS_MAX )
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

	   double val_HCur = sqrt((marrStrSK_VS [1] + R_ZEMLI) * (marrStrSK_VS [1] + R_ZEMLI)
		 + marrStrSK_VS [0] * marrStrSK_VS [0]) - R_ZEMLI;
	   if (val_HCur > valH )
	   {
		valH = val_HCur;
	   }
	   else
	   {
	   int uuu=0;
	   }


	  iNupPointsOut++;

	  }
	 }


  }
   double valGeoAlf = atan(marrStrSK_VS [0]/ (marrStrSK_VS [1] + R_ZEMLI));
   valGeoDist =  R_ZEMLI * valGeoAlf;

	 if (wcharrPath)
	 {
	 for (int j = 1; j < QUANT_COLS ; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,0  // ????? ?????????? ?? ??? X
								  ,j  // ????? ?????????? ?? ??? Y
								  ,100 //  ??????? ?? ??? X
								  ,pscaleY[j]  // ??????? ?? ??? Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,2  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,3  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete parrBuff ;
  delete pwcharrFileNames ;
  delete pscaleY ;

   }

	// fncCalcMtrxPartialDeriv( Environment, VAlStepInt );
	// fncCalcVectPartialDeriv_CoeffForm(Environment, VAlStepInt, marrDelta[8] ) ;
	// fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt,marrDelta[9] * mLearnShellBody.mMass ) ;

}


// ???????????? !!!
// ????????????? ????????? ?? ???????  ???????
// VAlScaleTime - ?????????? ?? ??? ???????
// VAlScalePos - ?????????? ?? ???? ????????? xyz
// OUTPUT:
// valDHoriz - ??????????????? ????? ????? ???????
void TLearnShellTraj::fncMoveClass_TO_ZeroAlt_AND_ShowGraphs(TEnvironment Environment
,const double VAlStepInt, wchar_t *wcharrPath1,const double VAlScaleTime
, const double VAlScalePos,  double &valGeoDist,double &valH )
{
  const int QUANT_COLS = 9 , QUANT_POINTS_MAX = 300000;
  const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ???????????? ????? ????? ??????????
   wchar_t 	wcharrPath[400] = {0};;
   valH = -1.;

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
	pscaleY[1] = VAlScalePos;
	pscaleY[2] = VAlScalePos;
	pscaleY[3] = VAlScalePos;
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


	  double val_r = sqrt((marrStrSK_VS [1] + R_ZEMLI)* (marrStrSK_VS [1] + R_ZEMLI)
						  + marrStrSK_VS [0]* marrStrSK_VS [0]);
	  ///

	  if (val_r <= R_ZEMLI)
	  {
		break;
	  }


	  if ( iNupPointsOut < QUANT_POINTS_MAX )
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

	   double val_HCur = sqrt((marrStrSK_VS [1] + R_ZEMLI) * (marrStrSK_VS [1] + R_ZEMLI)
		 + marrStrSK_VS [0] * marrStrSK_VS [0]) - R_ZEMLI;
	   if (val_HCur > valH )
	   {
		valH = val_HCur;
	   }
	   else
	   {
	   int uuu=0;
	   }


	  iNupPointsOut++;

	  }
	 }


  }
   double valGeoAlf = atan(marrStrSK_VS [0]/ (marrStrSK_VS [1] + R_ZEMLI));
   valGeoDist =  R_ZEMLI * valGeoAlf;

	 if (wcharrPath)
	 {
	 for (int j = 1; j < QUANT_COLS ; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,0  // ????? ?????????? ?? ??? X
								  ,j  // ????? ?????????? ?? ??? Y
								  ,VAlScaleTime//  ??????? ?? ??? X
								  ,pscaleY[j]  // ??????? ?? ??? Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,2  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,iNupPointsOut //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,3  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete parrBuff ;
  delete pwcharrFileNames ;
  delete pscaleY ;

   }

	// fncCalcMtrxPartialDeriv( Environment, VAlStepInt );
	// fncCalcVectPartialDeriv_CoeffForm(Environment, VAlStepInt, marrDelta[8] ) ;
	// fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt,marrDelta[9] * mLearnShellBody.mMass ) ;

}





// ????????????? ????????? ?? ???????  VAlFixedT

void TLearnShellTraj::fncMoveClass_TO_FixedTime_AND_ShowGraphs(TEnvironment Environment
,const double VAlStepInt, const double VAlFixedT ,wchar_t *wcharrPath1  )
{
	const int QUANT_COLS = 9 , QUANT_POINTS  = (VAlFixedT - mTCur) / VAlStepInt ;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ???????????? ????? ????? ??????????
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


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,QUANT_POINTS //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,0  // ????? ?????????? ?? ??? X
								  ,j  // ????? ?????????? ?? ??? Y
								  ,100 //  ??????? ?? ??? X
								  ,pscaleY[j]  // ??????? ?? ??? Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,QUANT_POINTS //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,2  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ???? ? ?????
								  , parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,QUANT_COLS // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,QUANT_POINTS //  - ?-?? ?????
								  ,pwcharrFileNames //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,1  // ????? ?????????? ?? ??? X
								  ,3  // ????? ?????????? ?? ??? Y
								  ,1.//  ??????? ?? ??? X
								  ,1.  // ??????? ?? ??? Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete parrBuff ;
  delete pwcharrFileNames ;
  delete pscaleY ;

   }

// fncCalcMtrxPartialDeriv( Environment, VAlStepInt );
//  fncCalcVectPartialDeriv_CoeffForm(Environment, VAlStepInt, marrDelta[8] ) ;
 //  fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt,marrDelta[9] * mLearnShellBody.mMass ) ;

}


// ?????????? ??????? ??????? ?????????? ???????? ??????? ??
// ?????????? ??????? ? ??????? i
void TLearnShellTraj::fncCalcVectPartialDeriv(TEnvironment Environment
  ,const double VAlStepInt, const int iVarNum,  double valDelta,  double *arrVectPartialDeriv )
{
  TLearnShellTraj shTrTemp = *this ;
  shTrTemp.mTCur = (*this).mTStart  ;
  shTrTemp.fncFillNachalnieUsloviaVS();
  shTrTemp.marrStrSK_VS [iVarNum] += valDelta;
  shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
  MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartialDeriv);
}

// ?????????? ??????? ??????? ?????????? ???????? ??????? ??
// ????? ???????????? ????? ix
void TLearnShellTraj::fncCalcVectPartialDeriv_Coef_Cx(TEnvironment Environment,const double VAlStepInt
	,const  double valDelta, double *arrVectPartDerivCoeff_Cx)
{
	TLearnShellTraj shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.mLearnShellBody.mplnCx.stretchDiagrAlongXY(1., ( 1. + valDelta)) ;

  shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartDerivCoeff_Cx);
}

//-------------------------------------------------------------------------------------
// ?????????? ??????? ??????? ?????????? ???????? ??????? ??
// ????? Cz
void TLearnShellTraj::fncCalcVectPartialDeriv_Coef_Cz(TEnvironment Environment,const double VAlStepInt
	,const  double valDelta, double *arrVectPartDerivCoeff_Cz)
{
	TLearnShellTraj shTrTemp = *this ;
	shTrTemp.mTCur = (*this).mTStart  ;
	shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.mLearnShellBody.mplnKnm.stretchDiagrAlongXY(1., ( 1. + valDelta)) ;


  shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartDerivCoeff_Cz);
}

// ?????????? ??????? ??????? ?????????? ???????? ??????? ??
// ?????
void TLearnShellTraj::fncCalcVectPartialDeriv_Mass(TEnvironment Environment
	,const double VAlStepInt,const  double valDelta, double *arrVectPartDerivMass)
{
  TLearnShellTraj shTrTemp = *this ;
  shTrTemp.mTCur = (*this).mTStart  ;
  shTrTemp.fncFillNachalnieUsloviaVS();
	shTrTemp.mLearnShellBody.mMass = mLearnShellBody.mMass + valDelta;
	//shTrTemp.mLearnShellBody.mV0 = (shTrTemp.mLearnShellBody.mMass - valDelta )/ shTrTemp.mLearnShellBody.mMass * shTrTemp.mLearnShellBody.mV0;
	//shTrTemp.mLearnShellBody.mOmega0 = (shTrTemp.mLearnShellBody.mMass - valDelta )/ shTrTemp.mLearnShellBody.mMass * shTrTemp.mLearnShellBody.mOmega0;
	shTrTemp.fncMovePhasVector(Environment,VAlStepInt, mTCur );
   double arrT0[8] = {0.} ;
  MtrxMinusMatrx( shTrTemp.marrStrSK_VS, marrStrSK_VS,8, 1, arrT0);
	MatrxMultScalar(arrT0, 8, 1, 1. /valDelta ,arrVectPartDerivMass);
}
/*
// ?????? ??????? ??????? ??????????? ?? ??? ????????
void TLearnShellTraj::fncCalcMtrxPartialDeriv(TEnvironment Environment, const double VAlStepInt )
{
 // arrStrSK_JacT - ??????????????? ?????? ??? ???????? ??????????. ????????? ???????? ?????????.
	 // ???????, ??? ????????? ???????  marrStrSK_Jac  ?????? arrStrSK_JacT ???? ???????????????
	double	arrStrSK_JacT[64]= {0.} ;

	for (int i = 0; i < 8; i++)
	{
	 fncCalcVectPartialDeriv(Environment,VAlStepInt, i, marrDelta[i], &arrStrSK_JacT[i*8] ) ;
	}

	MatrTransp(arrStrSK_JacT, 8, 8, marrStrSK_Jac);

 }
 */

// ?????? ??????? ??????? ??????????? ?? ?????????? ?????
// arrMtrx_Wind_PartialDeriv[8*3]  - ??????? ??????? ??????????? ?? ????? ????????, ????????????, ?????? ????????
void TLearnShellTraj::fncCalcMtrxTransp_Wind_PartialDeriv(TEnvironment Environment
, const double VAlStepInt, double *arrMtrxTransp_Wind_PartialDeriv )
{
	// arrStrSK_JacT - ??????????????? ?????? ??? ???????? ??????????. ????????? ???????? ?????????.
	// ???????, ??? ????????? ???????  marrStrSK_Jac  ?????? arrStrSK_JacT ???? ???????????????

	memset(arrMtrxTransp_Wind_PartialDeriv, 0, 24 * sizeof(double));

	TEnvironment Environment1 =  Environment;
	Environment1.mWind_V =  Environment.mWind_V + 0.1 ;
	TLearnShellTraj shTrTemp = *this ;
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
 // ???????? ??????? ?? ??? ? ???
 // ???? LEnArrVS = 3 ?? ??????????????? ?????? ?????????
 // ???? LEnArrVS = 6 ?? ??????????????? ?????? ?????????  ? ?????? ????????
 //
void TLearnShellTraj::transform_xyzGSK_To_xyzSSK( const int LEnArrVS, double *arrGSKInp, double *arrSSKOut)
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
 // ???????? ??????? ?? ????????????? ??? ?  ???
 // ???? LEnArrVS = 3 ?? ??????????????? ?????? ?????????
 // ???? LEnArrVS = 6 ?? ??????????????? ?????? ?????????  ? ?????? ????????

 void TLearnShellTraj::transform_xyzSSK_To_xyzGSK( const int LEnArrVS, double *arrSSKInp, double *arrGSKOut)
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
// ???????????? ??????? ???????? ?? ??????? ??????????? ?? ? ???
// ???? LEnArrVS = 3 ?? arrMtrxTransformOut ????????????? ?????? ?????????
 // ???? LEnArrVS  = 6 ?? arrMtrxTransformOut  ????????????? ?????? ?????????  ? ?????? ????????
 void TLearnShellTraj::createMtrxTransform_xyzSSK_To_xyzGSK( const int LEnArrVS,  double *arrMtrxTransformOut)
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

// ?????????? ????? ???????????? ?????????? ????? ???????????? ??????? ? ????
// INPUT:
// arrTargVS_SSK0[6] - ?????? ????????? ???? ? ??? ?? ?????? ?????? ???????? ???????
// VAl_dtInt - ??? ??????????????
// ??????????:
//  ??????????? ??????????
double TLearnShellTraj::calcPointMissMinimum(TEnvironment Environment, double *arrTargVS_SSK0, const double VAl_dtInt)
{

	double valF0prev = 10000000.0; // ???????- ??-? ?????. ????? ???????? ? ????? ?? ?????? T
	double valF0; // ???????- ??-? ?????. ????? ???????? ? ????? ?? ?????? T+dt

	double valMaxT = 600.;
	int iMqxQuantIter = valMaxT/ VAl_dtInt;
	for (int i =0; i < iMqxQuantIter; i++)
	{
		fncEilerStep( Environment, VAl_dtInt) ;
		// ?????????? ??????????? ?? ?????? valTCur ????????? ????
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


// ?????????? ??????? ??????? ??????????? ???????? ??????? ? ??????????? ???????? ?????????
// ?? ?????????? ???????? - ????????? ????????, ????? ?????, ???? ?????
void TLearnShellTraj::calcJacobian_8x10 (TEnvironment Environment,const double VAlStepInt
  , double *arrStrSK_Jacobian)
{
	 // LEN_ARR_SCATTERS - ????? ??????? ?????????? ?? ???????? ????????? ????????
		// ????????? ?? ??????? ????????? ????????
// 0. marrStrSK_VS [3]-  ???? ??? (???????)
// 1. marrStrSK_VS [5]-  ??????? ????????? ????????? ??
// 2. marrStrSK_VS [7]-  ???? ?????
// 3. marrStrSK_VS [6]-  ????????
// 4. ?????
// 5.????? ????? Cx
// 6. ???? ????? ?? ??? Z
// 7. ?????? ?????????????? ???????? ?????
// 8. ??????????? ??????????????? ?????
// 9. ?????? ????????????? ?????

		double	arrStrSK_JacT[LEN_ARR_SCATTERS * 8]= {0.};

		double arrDelta[] = {
								 0.001
								 ,0.01
								 ,0.001
								 ,1.
								 ,0.01
								 ,0.01
								 } ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 3, arrDelta[0], arrStrSK_JacT ) ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 5, arrDelta[1], &arrStrSK_JacT[8] ) ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 7, arrDelta[2], &arrStrSK_JacT[2 *8] ) ;
	fncCalcVectPartialDeriv(Environment,VAlStepInt, 6, arrDelta[3], &arrStrSK_JacT[3 *8] ) ;
	fncCalcVectPartialDeriv_Mass( Environment, VAlStepInt, arrDelta[4]* mLearnShellBody.mMass, &arrStrSK_JacT[4 * 8] )  ;
	fncCalcVectPartialDeriv_Coef_Cx( Environment,VAlStepInt
	, arrDelta[5],  &arrStrSK_JacT[5 * 8] )  ;
	fncCalcVectPartialDeriv_Coef_Cz(Environment, VAlStepInt
	,  arrDelta[4],  &arrStrSK_JacT[6 * 8] )  ;
	fncCalcMtrxTransp_Wind_PartialDeriv( Environment
	,VAlStepInt, &arrStrSK_JacT[7* 8]) ;

	MatrTransp(arrStrSK_JacT, LEN_ARR_SCATTERS, 8, arrStrSK_Jacobian);

}


// ?????????? ?????????????? ??????? ???????? ??????? ????????? ???????
// ? ???
// INPUT:
// Environment -??????? ?????
// VAlStepInt - ??? ??????????????
// VAlTFlight - ???????? ?????
// arrMtrxShellDisp  - ???????????? ??????? ?????????? ????????? ??????????
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - ??????? ??????? ??????????? ???????? ???????
// ??????? ? ????????????? ???????? ????????? ?? ??????????
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - ?????????????? ??????? ? ???
// arrShellScatteringsCorMtrxPos_SSK [3*3] - ?????????????? ??????? ?????? ????????? ????????? ? ???
// arrShellVS_GSK[6] - ?????? ????????? ??????? ? ???
void TLearnShellTraj::calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (TEnvironment Environment,const double VAlStepInt
	 ,const double VAlTFlight,double* arrMtrxShellDisp,  double *arrStrSK_Jacobian
	,double* arrShellScatteringsCorMtarx_GSK, double *arrShellVS_GSK, double* arrShellScatteringsCorMtrxPos_SSK)
{
	 fncMovePhasVector(Environment, VAlStepInt, VAlTFlight);
	 calcJacobian_8x10 ( Environment, VAlStepInt,  arrStrSK_Jacobian);

	 // ?????????? ?????? ???????? ??????? ????????? ??????? ? ??????????? ?? - arrCorMtrxTrajCK [8*8]
	 double arrT0[LEN_ARR_SCATTERS * 8] = {0.}, arrCorMtrxTrajCK[64] = {0.};
		MtrxMultMatrx(arrStrSK_Jacobian,8, LEN_ARR_SCATTERS, arrMtrxShellDisp,LEN_ARR_SCATTERS, arrT0) ;
		MtrxMultMatrxTransp(arrT0, 8, LEN_ARR_SCATTERS, arrStrSK_Jacobian,8, arrCorMtrxTrajCK) ;
		///

		// ?????????? ??????? ??????? ??????????? ??????? ??????? ????????? ??????? ??
		// ??????????? ?? ? ???????. ??
		double arrJac_PrStSK [ 6*8] = {0.};
		fncCalcJacobi_PrStSK(  arrJac_PrStSK  );
		///

		// ?????????? ????????. ??????? ?????? ? ??????? ??????????? ??  - arrCorMtrxPrStCK
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

		// ?????????? ??????? ???????? ?? ??????? ??????????? ?? ? ???
		double arrMtrxTransf_From_StPrSK_To_GSK [36];
		createMtrxTransform_xyzSSK_To_xyzGSK( 6, arrMtrxTransf_From_StPrSK_To_GSK);
		double arrTemp1 [ 36] = {0.};
		MtrxMultMatrx(arrMtrxTransf_From_StPrSK_To_GSK,6, 6, arrCorMtrxPrStCK ,6, arrTemp1) ;
		MtrxMultMatrxTransp(arrTemp1, 6, 6, arrMtrxTransf_From_StPrSK_To_GSK,6, arrShellScatteringsCorMtarx_GSK) ;
		///

		// ?????????? ??????? ????????? ??????? ? ???
		double arrVS_PrStSK [6] = {0.};
		fncCalcVS_v_PrStSK(arrVS_PrStSK) ;
		///
		MtrxMultMatrx(arrMtrxTransf_From_StPrSK_To_GSK,6, 6, arrVS_PrStSK ,1, arrShellVS_GSK) ;

}

//-------------------------------------------------------------------------------------


double TLearnShellTraj::findOptimalTab_Cx_76_(double *arrTab, int numRows, int numCols
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter)
{

	TLearnShellTraj ShellTrajCur = *this;
	//double valy = -1.;
	const int J0 = 3;

	double valFGr0 = ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;

	if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}




	for (int i =0;i < NUmIter; i++)
	{


		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ numRows);

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[23] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < 23; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnCx.Points[J0 +j].Y += 0.0001;
		double valy1 = ShellTrajCur1.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
		 arrGrad[j] = (valy1 - valFGr0)/ 0.0001;
	   }
	   ///

	   // ????? ????
	   double step = 0.005;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {
	   if (n ==39)
	   {
		 int iii =0;
	   }
		  step = step * 0.66;
		 double arrT0[23] = {0.};
		 MatrxMultScalar(arrGrad, 1, 23, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 for (int k =0; k < 23; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnCx.Points[J0 +k].Y  -=   arrT0[ k];
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}

	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = mLearnShellBody.mplnKnm;
	*pplnMxOmx = mLearnShellBody.mplnMxOmx;

	 return valFGr0;
}

//-------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------


double TLearnShellTraj::findOptimalTab_MxOmegax_76_(double *arrTab, int numRows, int numCols
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp, double *parrObjFnc, const int NUmIter)
{

	TLearnShellTraj ShellTrajCur = *this;
	//double valy = -1.;
	const int J0 = 3;

	double valFGr0 = ShellTrajCur.calcFGr_ForKnm_76(arrTab,  numRows,  numCols) ;
	if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}

	for (int i =0;i < NUmIter; i++)
	{
		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ numRows);

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[23] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < 23; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +j].Y += 0.01;
		double valy1 =ShellTrajCur1.calcFGr_ForKnm_76(arrTab,  numRows,  numCols) ;
		 arrGrad[j] = (valy1 - valFGr0)/ 0.01;
	   }
	   NormalizeVect(arrGrad, 23)  ;
	   ///

	   // ????? ????
	   double step = 0.01;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {
	   if (n ==39)
	   {
		 int iii =0;
	   }
		  step = step * 0.66;
		 double arrT0[23] = {0.};
		 MatrxMultScalar(arrGrad, 1, 23, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 bool bcontinue = false;
		 for (int k =0; k < 23; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +k].Y  -=   arrT0[ k];
		   if (ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +k].Y <=0.)
		   {
			 bcontinue = true;
			 break;
		   }
		 }

		 if (bcontinue)
		 {
           continue;
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForKnm_76(arrTab,  numRows,  numCols) ;
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;
	*pplnCx = mLearnShellBody.mplnCx;
	*pplnKnm = mLearnShellBody.mplnKnm;

	 return valFGr0;
}

//-------------------------------------------------------------------------------------


double TLearnShellTraj::findOptimalTab_Knm_76_(double *arrTab, int numRows, int numCols
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp, double *parrObjFnc, const int NUmIter)
{

	TLearnShellTraj ShellTrajCur = *this;
	//double valy = -1.;
	const int J0 = 3;

	double valFGr0 = ShellTrajCur.calcFGr_ForKnm_76(arrTab,  numRows,  numCols) ;
	if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}

	for (int i =0;i < NUmIter; i++)
	{
	   if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ numRows);

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[23] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < 23; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +j].Y += 0.01;
		double valy1 =ShellTrajCur1.calcFGr_ForKnm_76(arrTab,  numRows,  numCols) ;
		 arrGrad[j] = (valy1 - valFGr0)/ 0.01;
	   }
	   ///

	   // ????? ????
	   double step = 0.1;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {
	   if (n ==39)
	   {
		 int iii =0;
	   }
		  step = step * 0.66;
		 double arrT0[23] = {0.};
		 MatrxMultScalar(arrGrad, 1, 23, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 for (int k =0; k < 23; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +k].Y  -=   arrT0[ k];
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForKnm_76(arrTab,  numRows,  numCols) ;
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnCx = mLearnShellBody.mplnCx;
	*pplnMxOmx = mLearnShellBody.mplnMxOmx;
	 return valFGr0;
}

//-------------------------------------------------------------------------------------//-------------------------------------------------------------------------------------


double TLearnShellTraj::findOptimalTab_Cz_76_(double *arrTab, int numRows, int numCols
   , TURPolyLine *pplnCz, double *arrDisp)
{

	TLearnShellTraj ShellTrajCur = *this;
	//double valy = -1.;
	const int J0 = 3;

	double valFGr0 = ShellTrajCur.calcFGr_ForCz_76(arrTab,  numRows,  numCols) ;
	for (int i =0;i < 5; i++)
	{
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[23] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < 23; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnCz.Points[J0 +j].Y += 0.001;
		double valy1 =ShellTrajCur1.calcFGr_ForCz_76(arrTab,  numRows,  numCols) ;
		 arrGrad[j] = (valy1 - valFGr0)/ 0.001;
	   }
	   ///

	   // ????? ????
	   double step = 0.01;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 20; n++)
	   {
	   if (n ==39)
	   {
		 int iii =0;
	   }
		  step = step * 0.66;
		 double arrT0[23] = {0.};
		 MatrxMultScalar(arrGrad, 1, 23, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 for (int k =0; k < 23; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnCz.Points[J0 +k].Y  -=   arrT0[ k];
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForCz_76(arrTab,  numRows,  numCols) ;
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnCz = ShellTrajCur.mLearnShellBody.mplnCz;
	 return valFGr0;
}

//-------------------------------------------------------------------------------------
double TLearnShellTraj::calcFGr_ForCx_76(double *arrTab, int numRows, int numCols)
{
	TEnvironment Environment(0., M_PI/2., 0.);
	double sum = 0.;
	for (int i = 0; i < numRows; i++)
	{
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);


		double valDHoriz = -1.;
		ShellTrajCur.fncMoveShell_TO_ZeroAlt(Environment,0.001, valDHoriz );
		double delD = ShellTrajCur.marrStrSK_VS[0] -  arrTab [i *numCols +2];
		//double delVk = ShellTrajCur.marrStrSK_VS[6] -  arrTab [i *numCols +3];
		//double delTetk = ( -ShellTrajCur.marrStrSK_VS[7] -  arrTab [i *numCols +4]) * 1000.;

		//double delH = ShellTrajCur.marrStrSK_VS[1];
		sum = sum + delD* delD ;//+ delVk * delVk + delTetk * delTetk  + delH*delH;
	}
	return sum;
}

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
double TLearnShellTraj::calcFGr_ForCz_76(double *arrTab, int numRows, int numCols)
{
	TEnvironment Environment(0., M_PI/2., 0.);
	TEnvironment Environment1(10., 0., 0.);
	double sum = 0.;
	for (int i = 0; i < numRows; i++)
	{
        // ?????????? ???????? ?????
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		double valDHoriz = -1.;
		ShellTrajCur.fncMoveShell_TO_ZeroAlt(Environment,0.001, valDHoriz );

		 TLearnShellTraj ShellTrajCur1(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		 ShellTrajCur1.fncMoveShell_TO_ZeroAlt(Environment1,0.001, valDHoriz );

		 double valSnos =  fabs(ShellTrajCur1.marrStrSK_VS[2] -ShellTrajCur.marrStrSK_VS[2]);
		 ///
		double delD = valSnos -  arrTab [i *numCols +6];
		//double delVk = ShellTrajCur.marrStrSK_VS[6] -  arrTab [i *numCols +3];
		//double delTetk = ( -ShellTrajCur.marrStrSK_VS[7] -  arrTab [i *numCols +4]) * 1000.;

		//double delH = ShellTrajCur.marrStrSK_VS[1];
		sum = sum + delD* delD ;//+ delVk * delVk + delTetk * delTetk  + delH*delH;
	}
	return sum;
}

//-------------------------------------------------------------------------------------

double TLearnShellTraj::calcFGr_ForKnm_76(double *arrTab, int numRows, int numCols)
{
	TEnvironment Environment(0., M_PI/2., 0.);
	double sum = 0.;
	for (int i = 0; i < numRows; i++)
	{
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		double valDHoriz = -1.;
		ShellTrajCur.fncMoveShell_TO_ZeroAlt(Environment,0.001, valDHoriz );

		double delZ = ShellTrajCur.marrStrSK_VS[2] - arrTab[ i *numCols +5] ;
		sum = sum + delZ* delZ;
	}
	return sum;
}

//-------------------------------------------------------------------------------------
void TLearnShellTraj::calcSoglasCoeffsCix_76(double *arrTab, int  numRows,int  numCols, TURPolyLine *pplnCix)
{
   TEnvironment Environment(0., M_PI/2., 0.);



  (*pplnCix).Points[0].X =  0.;
  (*pplnCix).Points[0].Y =  1.;
  (*pplnCix).Points[numRows + 1].X =  M_PI /2.;
  (*pplnCix).Points[numRows + 1].Y =  1.;

  for (int i =0; i < numRows; i++)
  {
  (*pplnCix).Points[1 + i].X =  arrTab[ i *numCols +1] * 180. / M_PI;

	TLearnShellTraj ShellTraj0 = *this;
	double valDHoriz = -1.;

   ShellTraj0.fncMoveShell_TO_ZeroAlt(Environment,0.001, valDHoriz );
   double deltaf = ShellTraj0.marrStrSK_VS[0] - arrTab[ i *numCols +2];


   TLearnShellTraj ShellTrajCur = *this;
   ShellTrajCur.mLearnShellBody.mplnCx =   mLearnShellBody.mplnCx.MultScalar(1.01) ;
   ShellTrajCur.fncMoveShell_TO_ZeroAlt(Environment,0.001, valDHoriz );

   double valk = (ShellTrajCur.marrStrSK_VS[0] - ShellTraj0.marrStrSK_VS[0]) / 0.01;
   if (fabs(valk) < 0.0000001 )
   {
	 (*pplnCix).Points[1 + i].Y = 1.;
	 continue;
   }
   (*pplnCix).Points[1 + i].Y = 1. - deltaf /valk ;

   // ????????
   //	ShellTrajCur.mLearnShellBody.mplnCx =   mLearnShellBody.mplnCx.MultScalar(0. ) ;
  // ShellTrajCur.fncMovePhasVector(Environment,0.001, arrTab [i *numCols] );
  // double deltaf1 = ShellTrajCur.marrStrSK_VS[0] - arrTab[ i *6 +2];

	int ii=0;



  }

}


void TLearnShellTraj::estimateDisp_76(double *arrTab, const int numRows, const int numCols
	,const double  VAlSigTechTet, const double  VAlSigTechPsi
		, double  *pvalDispV, double  *pvalDispCx, double  *pvalDispM)
{
	TEnvironment Environment(0., M_PI/2., 0.);

	double arrVesselVelocity[3] = {0.};

	double valDHoriz = -1.;

	double *arrA   = new double [numRows  * 2];
	double *arrAT   = new double [numRows  * 2];
	double *arrb = new double [numRows ];
	 const double VAlStepInt =0.01;

	for (int i =0 ; i < numRows; i++)
	{
	   // ???????? ???????????
	   TLearnShellTraj ShellTraj0(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);

		ShellTraj0.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );



		///

		// ???????? ?? ?????
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] + 0.001, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj1.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) / 0.001;
	  ///

	  // ???????? ?? M
		  TLearnShellTraj ShellTraj4(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj4.mLearnShellBody.mMass += 0.0011;
	 ShellTraj4.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	 double   valPartialDeriv1 = (ShellTraj4.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) ;

		 ///

		 // ??????????  arrb
		 double temp0 = valPartialDeriv * VAlSigTechTet ;

		 arrb [i] = arrTab[ i *numCols +12] * arrTab[ i *numCols +12] -  temp0* temp0 - valPartialDeriv1*valPartialDeriv1;

		 if  (arrb [i] < 0.) arrb [i] = 0.;

		 // ???????? ?? V0
		   TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj2.marrStrSK_VS[6] += 5;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) / 5.;

		 arrA [ 2 * i ] =  valPartialDeriv * valPartialDeriv ;
		// arrA [3 * numRows + 3 * i ] =  arrVectPartialDeriv1[2] * arrVectPartialDeriv1[2];
		 ///




		  // ???????? ?? Cx
		  TLearnShellTraj ShellTraj3(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj3.mLearnShellBody.mplnCx = ShellTraj0.mLearnShellBody.mplnCx.MultScalar(1.01 ) ;
	 ShellTraj3.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj3.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) ;

		 arrA [ 2 * i +1] =  valPartialDeriv * valPartialDeriv ;
	   //	 arrA [3 * numRows + 3 * i + 1] =  arrVectPartialDeriv1[2] * arrVectPartialDeriv1[2]* 0.01* 0.01;
		 ///


	}

	// ????? ???? ?????????
		 double arrATA[4] = {0.}, arrATA_Inv[4] = {0.};
		 double arrATb[2] = {0.}, arrRez [2] = {0.};
		 MatrTransp(arrA, numRows , 2, arrAT);
		 MtrxMultMatrx(arrAT,2, numRows , arrA,2, arrATA) ;

		  MtrxMultMatrx(arrAT,2, numRows , arrb,1, arrATb) ;
		 InverseMtrx2(arrATA, arrATA_Inv);
		 MtrxMultMatrx(arrATA_Inv,2,2, arrATb,1, arrRez) ;
	  *pvalDispV = arrRez[0];
	  *pvalDispCx= arrRez[1];
	  *pvalDispM = 0.0011* mLearnShellBody.mMass   * 0.0011* mLearnShellBody.mMass  ;


	delete [] arrA;
	delete [] arrb;
}



void TLearnShellTraj::estimateDisp_76(double *arrTab, const int numRows, const int numCols
	,const double  VAlSigTechTet, const double  VAlSigTechPsi
		, double  *pvalDispV, double  *pvalDispCx, double  *pvalDispM
		, double  *pvalDispMxOmx)
{
	TEnvironment Environment(0., M_PI/2., 0.);

	double arrVesselVelocity[3] = {0.};

	double valDHoriz = -1.;

	double *arrA   = new double [numRows  * 3];
	double *arrAT   = new double [numRows  * 3];
	double *arrb = new double [numRows ];
	 const double VAlStepInt =0.01;

	for (int i =0 ; i < numRows; i++)
	{
	   // ???????? ???????????
	   TLearnShellTraj ShellTraj0(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);

		ShellTraj0.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );



		///

		// ???????? ?? ?????
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] + 0.001, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj1.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) / 0.001;
	  ///

	  // ???????? ?? M
		  TLearnShellTraj ShellTraj4(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj4.mLearnShellBody.mMass += 0.0011;
	 ShellTraj4.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	 double   valPartialDeriv1 = (ShellTraj4.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) ;

		 ///

		 // ??????????  arrb
		 double temp0 = valPartialDeriv * VAlSigTechTet ;

		 arrb [i] = arrTab[ i *numCols +12] * arrTab[ i *numCols +12] -  temp0* temp0 - valPartialDeriv1*valPartialDeriv1;

		 if  (arrb [i] < 0.) arrb [i] = 0.;

		 // ???????? ?? V0
		   TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj2.marrStrSK_VS[6] += 5;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) / 5.;

		 arrA [ 3 * i ] =  valPartialDeriv * valPartialDeriv ;
		// arrA [3 * numRows + 3 * i ] =  arrVectPartialDeriv1[2] * arrVectPartialDeriv1[2];
		 ///




		  // ???????? ?? Cx
		  TLearnShellTraj ShellTraj3(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj3.mLearnShellBody.mplnCx = ShellTraj0.mLearnShellBody.mplnCx.MultScalar(1.01 ) ;
	 ShellTraj3.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj3.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) ;

		 arrA [ 3 * i +1] =  valPartialDeriv * valPartialDeriv ;
	   //	 arrA [3 * numRows + 3 * i + 1] =  arrVectPartialDeriv1[2] * arrVectPartialDeriv1[2]* 0.01* 0.01;
		 ///
	   // ???????? ?? MxOmx
		  TLearnShellTraj ShellTraj5(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj5.mLearnShellBody.mplnKnm = ShellTraj0.mLearnShellBody.mplnKnm.MultScalar(1.01 ) ;
	 ShellTraj5.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj5.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] )/0.01 ;

		 arrA [ 3 * i +2] =  valPartialDeriv * valPartialDeriv ;

	}

	// ????? ???? ?????????
		 double arrATA[9] = {0.}, arrATA_Inv[9] = {0.};
		 double arrATb[3] = {0.}, arrRez [3] = {0.};
		 MatrTransp(arrA, numRows , 3, arrAT);
		 MtrxMultMatrx(arrAT,3, numRows , arrA,2, arrATA) ;

		  MtrxMultMatrx(arrAT,3, numRows , arrb,1, arrATb) ;
		 InverseMtrx3(arrATA, arrATA_Inv);
		 MtrxMultMatrx(arrATA_Inv,3,3, arrATb,1, arrRez) ;
	  *pvalDispV     = arrRez[0];
	  *pvalDispCx    = arrRez[1];
	  *pvalDispMxOmx = arrRez[2];
	  *pvalDispM = 0.0011* mLearnShellBody.mMass   * 0.0011* mLearnShellBody.mMass  ;


	delete [] arrA;
	delete [] arrb;
}

/*
void TLearnShellTraj::estimateDisp_76(double *arrTab, const int numRows, const int numCols
	,const double  VAlSigTechTet, const double  VAlSigTechPsi
		, double  *pvalDispV, double  *pvalDispCx, double  *pvalDispM)
{
	TEnvironment Environment(0., M_PI/2., 0.);

	double arrVesselVelocity[3] = {0.};

	double valDHoriz = -1.;

	double *arrA   = new double [numRows  * 2];
	double *arrAT   = new double [numRows  * 2];
	double *arrb = new double [numRows ];
	 const double VAlStepInt =0.01;

	for (int i =0 ; i < numRows; i++)
	{
	   // ???????? ???????????
	   TLearnShellTraj ShellTraj0(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);

		ShellTraj0.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );



		///

		// ???????? ?? ?????
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] + 0.001, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj1.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) / 0.001;
	  ///

	  // ???????? ?? M
		  TLearnShellTraj ShellTraj4(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj4.mLearnShellBody.mMass += 0.0011;
	 ShellTraj4.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	 double   valPartialDeriv1 = (ShellTraj4.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) ;

		 ///

		 // ??????????  arrb
		 double temp0 = valPartialDeriv * VAlSigTechTet ;

		 arrb [i] = arrTab[ i *numCols +12] * arrTab[ i *numCols +12] -  temp0* temp0 - valPartialDeriv1*valPartialDeriv1;

		 if  (arrb [i] < 0.) arrb [i] = 0.;

		 // ???????? ?? V0
		   TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj2.marrStrSK_VS[6] += 5;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) / 5.;

		 arrA [ 2 * i ] =  valPartialDeriv * valPartialDeriv ;
		// arrA [3 * numRows + 3 * i ] =  arrVectPartialDeriv1[2] * arrVectPartialDeriv1[2];
		 ///




		  // ???????? ?? Cx
		  TLearnShellTraj ShellTraj3(arrVesselVelocity, mLearnShellBody, arrTab [i *numCols +1] , M_PI/2.);
		   ShellTraj3.mLearnShellBody.mplnCx = ShellTraj0.mLearnShellBody.mplnCx.MultScalar(1.01 ) ;
	 ShellTraj3.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	   valPartialDeriv = (ShellTraj3.marrStrSK_VS[0] -  ShellTraj0.marrStrSK_VS[0] ) ;

		 arrA [ 2 * i +1] =  valPartialDeriv * valPartialDeriv ;
	   //	 arrA [3 * numRows + 3 * i + 1] =  arrVectPartialDeriv1[2] * arrVectPartialDeriv1[2]* 0.01* 0.01;
		 ///


	}

	// ????? ???? ?????????
		 double arrATA[4] = {0.}, arrATA_Inv[4] = {0.};
		 double arrATb[2] = {0.}, arrRez [2] = {0.};
		 MatrTransp(arrA, numRows , 2, arrAT);
		 MtrxMultMatrx(arrAT,2, numRows , arrA,2, arrATA) ;

		  MtrxMultMatrx(arrAT,2, numRows , arrb,1, arrATb) ;
		 InverseMtrx2(arrATA, arrATA_Inv);
		 MtrxMultMatrx(arrATA_Inv,2,2, arrATb,1, arrRez) ;
	  *pvalDispV = arrRez[0];
	  *pvalDispCx= arrRez[1];
	  *pvalDispM = 0.0011* mLearnShellBody.mMass   * 0.0011* mLearnShellBody.mMass  ;


	delete [] arrA;
	delete [] arrb;
}
*/

double  TLearnShellTraj::fncPartialDerivD_TochkiPadenia_po_Tetta0(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);
	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, VAlTetta0 + 0.001, M_PI/2.);
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) / 0.001;
	 return valPartialDeriv;
}



double  TLearnShellTraj::fncPartialDerivD_TochkiPadenia_po_Mass(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, VAlTetta0 , M_PI/2.);
	 double delm = ShellTraj1.mLearnShellBody.mMass  * 0.01;
	 ShellTraj2.mLearnShellBody.mMass += delm ;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) / delm ;
	 return valPartialDeriv;
}


double  TLearnShellTraj::fncPartialDerivD_TochkiPadenia_po_V0(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, VAlTetta0 , M_PI/2.);

	 ShellTraj2. marrStrSK_VS[6] += 1.;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	  double valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) ;
	 return valPartialDeriv;
}


double  TLearnShellTraj::fncPartialDerivD_TochkiPadenia_po_Cx(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, VAlTetta0 , M_PI/2.);

	ShellTraj2.mLearnShellBody.mplnCx = ShellTraj1.mLearnShellBody.mplnCx.MultScalar(1.01 ) ;
	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	double  valPartialDeriv = (ShellTraj2.marrStrSK_VS[0] -  ShellTraj1.marrStrSK_VS[0] ) /0.01;

	 return valPartialDeriv;
}

double  TLearnShellTraj::fncPartialDerivZ_TochkiPadenia_po_Psi0(const double VAlTetta0)
{
	 TEnvironment Environment(0.,0.,0.);

	 double arrVesselVelocity[3] = {0.}, valDHoriz = -1., VAlStepInt = 0.001;
	 TLearnShellTraj ShellTraj1(arrVesselVelocity, mLearnShellBody, VAlTetta0, M_PI/2.);
	 ShellTraj1.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );

	 TLearnShellTraj ShellTraj2(arrVesselVelocity, mLearnShellBody, VAlTetta0 , M_PI/2.);
	 ShellTraj2. marrStrSK_VS[3] += 0.001;

	 ShellTraj2.fncMoveShell_TO_ZeroAlt(Environment,VAlStepInt, valDHoriz );
	double  valPartialDeriv = (ShellTraj2.marrStrSK_VS[2] -  ShellTraj1.marrStrSK_VS[2] ) /0.001;

	 return valPartialDeriv;
}


double TLearnShellTraj::arcSin(const double x)
{
double y = x;
if (fabs(y) > 0.99999999)
{
y = 0.99999999 * SIGNUM0(y);
}
return asin(y);
}

double TLearnShellTraj::SIGNUM0(const double y)
{
 return (y >=0.)?1.:-1.;
}

//---------------------------------------------------------
// 0. marrStrSK_VS [3]-  ???? ??? (???????)
// 1. marrStrSK_VS [5]-  ??????? ????????? ????????? ??
// 2. marrStrSK_VS [7]-  ???? ?????
// 3. marrStrSK_VS [6]-  ????????
// 4. ?????
// 5.????? ????? Cx
// 6. ???? ????? ?? ??? Z
// 7. ?????? ?????????????? ???????? ?????
// 8. ??????????? ??????????????? ?????
// 9. ?????? ????????????? ?????
void TLearnShellTraj::createPictures_Barrier_No1 (wchar_t *pwchGraphDir
,TEnvironment Environment,const double VAlStepInt
	 ,const double VAlTFlight,const double VAlScale
	 ,const double VAlSigPsi,const double VAlSigPi,const double VAlSigTetta
	 ,const double VAlSigV0,const double VAlSigMass,const double VAlSigCx
	 ,const double VAlSigWindV,const double VAlSigWindAlf,const double VAlSigTime
	 )
{
  // 1. ?????????? ?????????? ?? ????? ???????
  TLearnShellTraj LearnShellTraj = *this;
	 double valDHoriz = -1., val_H = 0.;
	 const double VAlScaleTime = 1.;
	 const double VAlScalePos = 1.;

	 LearnShellTraj.fncMoveClass_TO_ZeroAlt_AND_ShowGraphs
		(Environment, VAlStepInt, pwchGraphDir, VAlScaleTime,VAlScalePos , valDHoriz, val_H);
		///


// ?????????? ???? ????????? (???????????) ?? ??????? ?????????? ????????
// ? ???????? ?????? ???????
double arrMtrxShellDisp [LEN_ARR_SCATTERS*LEN_ARR_SCATTERS] = {0.};
arrMtrxShellDisp[0] = VAlSigPsi* VAlSigPsi; // Psi
arrMtrxShellDisp[LEN_ARR_SCATTERS   + 1] = VAlSigPi * VAlSigPi; // Pi
arrMtrxShellDisp[2*LEN_ARR_SCATTERS + 2] = VAlSigTetta * VAlSigTetta; // Tetta
arrMtrxShellDisp[3*LEN_ARR_SCATTERS + 3] = VAlSigV0 * VAlSigV0 ;
arrMtrxShellDisp[4*LEN_ARR_SCATTERS + 4] = VAlSigMass * VAlSigMass;
arrMtrxShellDisp[5*LEN_ARR_SCATTERS + 5] = VAlSigCx * VAlSigCx ;
arrMtrxShellDisp[7*LEN_ARR_SCATTERS + 7] = VAlSigWindV * VAlSigWindV;
arrMtrxShellDisp[8*LEN_ARR_SCATTERS + 8] = VAlSigWindAlf * VAlSigWindAlf;


///

double arrStrSK_Jacobian[8 *LEN_ARR_SCATTERS] ={0.};
double arrShellScatteringsCorMtarx_GSK[36] = {0.};
double arrShellVS_GSK[6] = {0.};
double arrShellScatteringsCorMtrxPos_SSK[9] = {0.};



	  TLearnShellTraj LearnShellTrajCur = *this;
		 LearnShellTrajCur.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK(Environment, VAlStepInt
	 ,VAlTFlight,  arrMtrxShellDisp,
	 arrStrSK_Jacobian, arrShellScatteringsCorMtarx_GSK,
	  arrShellVS_GSK,  arrShellScatteringsCorMtrxPos_SSK);
	  double p[100] = {0.};
	  p[0] = VAlTFlight;

	  p[1] = fabs(  arrStrSK_Jacobian[0]) *VAlSigPsi;// dX/dPsi
	  p[2] = fabs(  arrStrSK_Jacobian[2])* VAlSigTetta; // dX/dTet
	  p[3] = fabs(  arrStrSK_Jacobian[3])* VAlSigV0 ; // dX/dV
	  p[4] = fabs(  arrStrSK_Jacobian[4])* VAlSigMass; // dX/dm
	  p[5] = fabs(  arrStrSK_Jacobian[5])* VAlSigCx ; // dX/dCx
	  p[6] = fabs(  arrStrSK_Jacobian[1])* VAlSigPi; // dX/dPi
	  p[7] = fabs(  arrStrSK_Jacobian[7])* VAlSigWindV;// dx/dVwind
	  p[8] = fabs(  arrStrSK_Jacobian[8])* VAlSigWindAlf;  // dx/dAlfWind
	  p[9] = fabs(   LearnShellTrajCur.marrStrSK_VS[6]
		  * cos(LearnShellTrajCur.marrStrSK_VS[7])) *VAlSigTime ; // Vx

	  p[10]  = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS    ])*VAlSigPsi;// dY/dPsi
	  p[11]  = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS + 2])* VAlSigTetta; // dY/dTet
	  p[12]  = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS + 3])* VAlSigV0 ; // dY/dV
	  p[13]  = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS + 4])* VAlSigMass; // dY/dm
	  p[14] = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS + 5])* VAlSigCx ; // dY/dCx
	  p[15] = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS  + 1])* VAlSigPi; // dY/dPi
	  p[16] = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS  + 7])* VAlSigWindV; // dY/dVwind
	  p[17] = fabs(  arrStrSK_Jacobian[LEN_ARR_SCATTERS  + 8])* VAlSigWindAlf;  // dY/dAlfWind
	  p[18] = fabs(    LearnShellTrajCur.marrStrSK_VS[6]
		   * sin(LearnShellTrajCur.marrStrSK_VS[7]))*VAlSigTime ;   // Vy

	  p[19]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS    ])*VAlSigPsi;// dZ/dPsi
	  p[20]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS + 2])* VAlSigTetta; // dZ/dTet
	  p[21]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS + 3])* VAlSigV0 ;// dZ/dV
	  p[22]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS + 4])* VAlSigMass; // dZ/dm
	  p[23]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS + 5])* VAlSigCx ; // dZ/dCx
	  p[24]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS + 1])* VAlSigPi; // dz/dPi
	  p[25]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS  + 7])* VAlSigWindV; // dz/dVwind
	  p[26]  = fabs(  arrStrSK_Jacobian[2 *LEN_ARR_SCATTERS  + 8])* VAlSigWindAlf;  // dz/dAlfWind


	 TURPolyLine Pln_DelXY_po_Psi0   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[1], p[10],  VAlScale) ;
	TURPolyLine Pln_DelXZ_po_Psi0   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[2]
	, p[1], p[19], VAlScale) ;

	TURPolyLine Pln_DelXY_po_Tetta0   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[2], p[11],  VAlScale) ;
	TURPolyLine Pln_DelXZ_po_Tetta0   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[2]
	, p[2], p[20],  VAlScale) ;

	TURPolyLine Pln_DelXY_po_V0   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[3], p[12],  VAlScale) ;

	TURPolyLine Pln_DelXYCx   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[5], p[14], VAlScale) ;

	TURPolyLine Pln_DelXY_m   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[4], p[13],  VAlScale) ;

	TURPolyLine Pln_DelXY_pi   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[6], p[15],  VAlScale) ;

	TURPolyLine Pln_DelXY_Vwind   =  TURPolyLine( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[7], p[16],  VAlScale) ;

	TURPolyLine Pln_DelXZ_Vwind   =  TURPolyLine( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[2]
	, p[7], p[25],  VAlScale) ;

	TURPolyLine Pln_DelXY_Alfwind   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[8], p[17],  VAlScale) ;

	TURPolyLine Pln_DelXZ_Alfwind  =  TURPolyLine( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[2]
	, p[8], p[26],  VAlScale) ;

	TURPolyLine Pln_DelXY_t   ( LearnShellTrajCur.marrStrSK_VS[0], LearnShellTrajCur.marrStrSK_VS[1]
	, p[9], p[18],  VAlScale) ;

  wchar_t  wchFileName4[400] = {0};
  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_Psi0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_po_Psi0, 1) ;

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXZ_po_Psi0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXZ_po_Psi0, 1) ;

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_Tetta0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_po_Tetta0, 1) ;

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXZ_po_Tetta0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXZ_po_Tetta0, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_V0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_po_V0, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_Cx.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXYCx, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_m.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_m, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_pi.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_pi, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_Vwind.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_Vwind, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXZ_po_Vwind.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXZ_Vwind, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_Alfwind.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_Alfwind, 1) ;

  wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXZ_po_Alfwind.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXZ_Alfwind, 1) ;

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Pln_DelXY_po_t.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName4,&Pln_DelXY_t, 1) ;
  ///


  // ?????????? ??????? ????????? ??? arrShellScatteringsCorMtrxPos_SSK ? ???????????? ?????????
  double arrElK_XY0[4] = {0.};
  arrElK_XY0[0] = arrShellScatteringsCorMtrxPos_SSK [0];
  arrElK_XY0[1] = arrShellScatteringsCorMtrxPos_SSK [1];
  arrElK_XY0[2] = arrElK_XY0[1];
  arrElK_XY0[3] = arrShellScatteringsCorMtrxPos_SSK [4];
  double valLambX = -1.,valLambY = -1.;
  calcSobstChislaPoOsiam(arrElK_XY0,&valLambX,&valLambY);
  double t0 = sqrt(valLambX);
  double t1 = sqrt(valLambY);


  TURPolygon EllipsXY0 =   TURPolygon::createEllips(arrElK_XY0, 1, 1001);
  TURPointXY pntSdvigXY(LearnShellTrajCur.marrStrSK_VS[0]  , LearnShellTrajCur.marrStrSK_VS[1]);
  EllipsXY0 =  EllipsXY0.LinTransform(0. ,  pntSdvigXY, VAlScale );

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\EllipsXY0.shp");
  TURPolygon::WriteSetSHPFiles(wchFileName4,&EllipsXY0, 1) ;

  double arrElK_XZ0[4] = {0.};
  arrElK_XZ0[0] = arrShellScatteringsCorMtrxPos_SSK [0];
  arrElK_XZ0[1] = arrShellScatteringsCorMtrxPos_SSK [2];
  arrElK_XZ0[2] = arrElK_XZ0[1];
  arrElK_XZ0[3] = arrShellScatteringsCorMtrxPos_SSK [8];

  TURPolygon EllipsXZ0 =  TURPolygon::createEllips(arrElK_XZ0, 1, 1001);
  TURPointXY pntSdvigXZ(LearnShellTrajCur.marrStrSK_VS[0]  , LearnShellTrajCur.marrStrSK_VS[2]);
  EllipsXZ0 =  EllipsXZ0.LinTransform(0. ,  pntSdvigXZ, VAlScale );

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\EllipsXZ0.shp");
  TURPolygon::WriteSetSHPFiles(wchFileName4,&EllipsXZ0, 1) ;
  ///

  // ?????????? ??????? ??? ????? ???????? ?? V0
  double arrElK_withoutV0[4] = {0.};
  arrElK_withoutV0[0] = arrShellScatteringsCorMtrxPos_SSK [0]- p[3] * p[3] ;
  arrElK_withoutV0[1] = arrShellScatteringsCorMtrxPos_SSK [1]-  p[3] * p[12];
  arrElK_withoutV0[2] = arrElK_withoutV0[1];
  arrElK_withoutV0[3] = arrShellScatteringsCorMtrxPos_SSK [4] -  p[12] * p[12];
  calcSobstChislaPoOsiam(arrElK_withoutV0,&valLambX,&valLambY);
  double t2 = sqrt(valLambX);
  double t3 = sqrt(valLambY);

  TURPolygon EllipswithoutV0 =   TURPolygon::createEllips(arrElK_withoutV0, 1, 1001);

  EllipswithoutV0 =  EllipswithoutV0.LinTransform(0. ,  pntSdvigXY, VAlScale );

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\EllipsXY_WithoutV0.shp");
  TURPolygon::WriteSetSHPFiles(wchFileName4,&EllipswithoutV0, 1) ;
 ///

 // ?????????? ??????? ? ?????? ???????? ???????
 double arrElK_WithDeltaT[4] = {0.};
  arrElK_WithDeltaT[0] = arrShellScatteringsCorMtrxPos_SSK [0]+ p[9] * p[9] ;
  arrElK_WithDeltaT[1] = arrShellScatteringsCorMtrxPos_SSK [1]+ p[9] * p[18];
  arrElK_WithDeltaT[2] = arrElK_WithDeltaT[1];
  arrElK_WithDeltaT[3] = arrShellScatteringsCorMtrxPos_SSK [4] +  p[18] * p[18];

  TURPolygon EllipsWithDeltaT =   TURPolygon::createEllips(arrElK_WithDeltaT, 1, 1001);
  EllipsWithDeltaT =  EllipsWithDeltaT.LinTransform(0. ,  pntSdvigXY, VAlScale );

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\EllipsWithDeltaT.shp");
  TURPolygon::WriteSetSHPFiles(wchFileName4,&EllipsWithDeltaT, 1) ;

   double arrElK_With_dV0_With_dt_XZ[4] = {0.};
  arrElK_With_dV0_With_dt_XZ[0] = arrElK_WithDeltaT [0];
  arrElK_With_dV0_With_dt_XZ[1] = arrShellScatteringsCorMtrxPos_SSK [2];
  arrElK_With_dV0_With_dt_XZ[2] = arrElK_With_dV0_With_dt_XZ[1];
  arrElK_With_dV0_With_dt_XZ[3] = arrShellScatteringsCorMtrxPos_SSK [8];

  TURPolygon EllipsXZ_With_dV0_With_dt =  TURPolygon::createEllips(arrElK_With_dV0_With_dt_XZ, 1, 1001);

  EllipsXZ_With_dV0_With_dt =  EllipsXZ_With_dV0_With_dt.LinTransform(0. ,  pntSdvigXZ, VAlScale );

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\EllipsXZ_With_dV0_With_dt.shp");
  TURPolygon::WriteSetSHPFiles(wchFileName4,&EllipsXZ_With_dV0_With_dt, 1) ;

  ///
  // ?????????? ??????? ? ?????? ???????? ??????? ?? ??? ????? ???????? ????????? ????????
 double arrElK_WithoutdV0_but_WithDeltaT[4] = {0.};
  arrElK_WithoutdV0_but_WithDeltaT[0] = arrShellScatteringsCorMtrxPos_SSK [0]+ p[9] * p[9] - p[3] * p[3] ;;
  arrElK_WithoutdV0_but_WithDeltaT[1] = arrShellScatteringsCorMtrxPos_SSK [1]+ p[9] * p[18]-  p[3] * p[12];
  arrElK_WithoutdV0_but_WithDeltaT[2] = arrElK_WithoutdV0_but_WithDeltaT[1];
  arrElK_WithoutdV0_but_WithDeltaT[3] = arrShellScatteringsCorMtrxPos_SSK [4] +  p[18] * p[18]-  p[12] * p[12];

  TURPolygon Ellips_Without_dV0_WithDeltaT =   TURPolygon::createEllips(arrElK_WithoutdV0_but_WithDeltaT, 1, 1001);
  Ellips_Without_dV0_WithDeltaT =  Ellips_Without_dV0_WithDeltaT.LinTransform(0. ,  pntSdvigXY, VAlScale );

   wcscpy(wchFileName4, pwchGraphDir );
  wcscat(wchFileName4, L"\\Ellips_Without_dV0_WithDeltaT.shp");
  TURPolygon::WriteSetSHPFiles(wchFileName4,&Ellips_Without_dV0_WithDeltaT, 1) ;

  ///

 // ?????????? ???????? ???????? ??????????? ? ??????????? ?? ????????? ???????

 double valT = LearnShellTraj.mTCur;
 double valSt =  0.1;
 int NumRows = valT/valSt;

 int NumCols = 13;
 double *parrBuff = new double [NumRows * NumCols] ;
 memset (parrBuff, 0, NumRows * NumCols * sizeof(double)) ;

	  const int lenName =50 ;// ???????????? ????? ????? ??????????

	wchar_t *pwcharrFileNames1 = new wchar_t [ NumCols * lenName] ;
	memset (pwcharrFileNames1, 0, NumCols * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames1[ 0 * lenName], L"t");

	wcscpy( &pwcharrFileNames1[ 1 * lenName], L"AxeX_EllipsXY0");
	wcscpy( &pwcharrFileNames1[ 2* lenName],  L"AxeY_EllipsXY0");
	wcscpy( &pwcharrFileNames1[ 3 * lenName], L"AxeX_EllipswithoutV0");
	wcscpy( &pwcharrFileNames1[ 4 * lenName], L"AxeY_EllipswithoutV0");
	wcscpy( &pwcharrFileNames1[ 5 * lenName], L"AxeX_EllipsWithDeltaT");
	wcscpy( &pwcharrFileNames1[ 6 * lenName], L"AxeY_EllipsWithDeltaT");
	wcscpy( &pwcharrFileNames1[ 7 * lenName], L"AxeX_Ell_Not_dV0_WithDeltaT");
	wcscpy( &pwcharrFileNames1[ 8 * lenName], L"AxeY_Ell_Not_dV0_WithDeltaT");
	wcscpy( &pwcharrFileNames1[ 9 * lenName], L"GaneX");
	wcscpy( &pwcharrFileNames1[ 10 * lenName], L"GaneY");
	wcscpy( &pwcharrFileNames1[ 11 * lenName], L"AxeZ_EllipsXZ0");
	wcscpy( &pwcharrFileNames1[ 12* lenName], L"d");


	double arrMtrxShellDisp_1 [LEN_ARR_SCATTERS*LEN_ARR_SCATTERS] = {0.};
	memcpy (arrMtrxShellDisp_1, arrMtrxShellDisp, LEN_ARR_SCATTERS *LEN_ARR_SCATTERS * sizeof(double));
	arrMtrxShellDisp_1[3*LEN_ARR_SCATTERS + 3] = 0.;// ?? ??? ?????????
for (int i =0; i < NumRows; i++)
{
  double *pp = &parrBuff[i * NumCols];
  double valTCur = (1. + ((double)i)) *  valSt;

  pp[0] = valTCur;

   TLearnShellTraj LearnShellTrajCur = *this;

   LearnShellTrajCur.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK(Environment, VAlStepInt
	 ,valTCur,  arrMtrxShellDisp,
	 arrStrSK_Jacobian, arrShellScatteringsCorMtarx_GSK,
	  arrShellVS_GSK,  arrShellScatteringsCorMtrxPos_SSK);
  // 1 ??? ??????? ? ?????? ???????? ?? V0 ??? ????? ???????? ???????
	  double arrElK_XY0[4] = {0.};
  arrElK_XY0[0] = arrShellScatteringsCorMtrxPos_SSK [0];
  arrElK_XY0[1] = arrShellScatteringsCorMtrxPos_SSK [1];
  arrElK_XY0[2] = arrElK_XY0[1];
  arrElK_XY0[3] = arrShellScatteringsCorMtrxPos_SSK [4];
  double valLambX = -1.,valLambY = -1.;
  calcSobstChislaPoOsiam(arrElK_XY0,&valLambX,&valLambY);
  pp[1] = sqrt(valLambX);
  pp[2] = sqrt(valLambY);
  ///

  // 2 ?????????? ??????? ?  ?????? ???????? ?? V0  ? ? ?????? ???????? ???????
 double arrElK_With_dV0_With_dt[4] = {0.};
 double val_dX_po_dt =    LearnShellTrajCur.marrStrSK_VS[6]
		  * cos(LearnShellTrajCur.marrStrSK_VS[7]) *VAlSigTime ; // Vx
 double val_dY_po_dt =     LearnShellTrajCur.marrStrSK_VS[6]
		   * sin(LearnShellTrajCur.marrStrSK_VS[7])*VAlSigTime ;   // Vy
  arrElK_With_dV0_With_dt[0] = arrShellScatteringsCorMtrxPos_SSK [0]+ val_dX_po_dt * val_dX_po_dt ;
  arrElK_With_dV0_With_dt[1] = arrShellScatteringsCorMtrxPos_SSK [1]+ val_dX_po_dt * val_dY_po_dt;
  arrElK_With_dV0_With_dt[2] = arrElK_With_dV0_With_dt[1];
  arrElK_With_dV0_With_dt[3] = arrShellScatteringsCorMtrxPos_SSK [4] +  val_dY_po_dt * val_dY_po_dt;
   calcSobstChislaPoOsiam(arrElK_With_dV0_With_dt,&valLambX,&valLambY);
   pp[5] = sqrt(valLambX);
  pp[6] = sqrt(valLambY);
//

   // 3 ?????????? ??????? ??? ????? ???????? ?? V0 ? ??? ???????? ???????

   TLearnShellTraj LearnShellTrajCur_1 = *this;
  LearnShellTrajCur_1.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK(Environment, VAlStepInt
	 ,valTCur,  arrMtrxShellDisp_1,
	 arrStrSK_Jacobian, arrShellScatteringsCorMtarx_GSK,
	  arrShellVS_GSK,  arrShellScatteringsCorMtrxPos_SSK);
  // ??? ??????? ? ?????? ???????? ?? V0 ??? ????? ???????? ???????
	  double arrElK_Not_dV0_Not_dt[4] = {0.};
  arrElK_Not_dV0_Not_dt[0] = arrShellScatteringsCorMtrxPos_SSK [0];
  arrElK_Not_dV0_Not_dt[1] = arrShellScatteringsCorMtrxPos_SSK [1];
  arrElK_Not_dV0_Not_dt[2] = arrElK_Not_dV0_Not_dt[1];
  arrElK_Not_dV0_Not_dt[3] = arrShellScatteringsCorMtrxPos_SSK [4];

  calcSobstChislaPoOsiam(arrElK_Not_dV0_Not_dt,&valLambX,&valLambY);

  pp[3] = sqrt(valLambX);
  pp[4] = sqrt(valLambY);

  // 4 ?????????? ??????? ??? ????? ???????? ?? V0  ? ? ?????? ???????? ???????
 double arrElK_Not_dV0_With_dt[4] = {0.};
// double val_dX_po_dt =    LearnShellTrajCur.marrStrSK_VS[6]
 //		  * cos(LearnShellTrajCur.marrStrSK_VS[7]) *VAlSigTime ; // Vx
 //double val_dY_po_dt =     LearnShellTrajCur.marrStrSK_VS[6]
 //		   * sin(LearnShellTrajCur.marrStrSK_VS[7])*VAlSigTime ;   // Vy
  arrElK_Not_dV0_With_dt[0] = arrShellScatteringsCorMtrxPos_SSK [0]+ val_dX_po_dt * val_dX_po_dt ;
  arrElK_Not_dV0_With_dt[1] = arrShellScatteringsCorMtrxPos_SSK [1]+ val_dX_po_dt * val_dY_po_dt;
  arrElK_Not_dV0_With_dt[2] = arrElK_Not_dV0_With_dt[1];
  arrElK_Not_dV0_With_dt[3] = arrShellScatteringsCorMtrxPos_SSK [4] +  val_dY_po_dt * val_dY_po_dt;
   calcSobstChislaPoOsiam(arrElK_Not_dV0_With_dt,&valLambX,&valLambY);

  pp[7] = sqrt(valLambX);
  pp[8] = sqrt(valLambY);
  ///
  pp[9] = (pp[5] - pp[7])/pp[5] * 100;
  pp[10] = (pp[6] - pp[8])/pp[6] * 100;

  pp[11] = sqrt(arrShellScatteringsCorMtrxPos_SSK [8]);
  pp[12]=  Norm3( arrShellVS_GSK) ;


  }

  for (int j = 1; j < (NumCols-1) ; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(pwchGraphDir  // ???? ? ?????
								  ,parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,NumCols // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,NumRows //  - ?-?? ?????
								  ,pwcharrFileNames1 //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,0  // ????? ?????????? ?? ??? X
								  ,j  // ????? ?????????? ?? ??? Y
								  ,1//  ??????? ?? ??? X
								  ,1  // ??????? ?? ??? Y
								   );
	/// ???????? ?? ????????? ?????????
	TYrWriteShapeFile::WriteOneReport(pwchGraphDir  // ???? ? ?????
								  ,parrBuff // ?????? ? ??????????? - ??????? nBuffRows x nBuffCols
								  ,NumCols // - ?-?? ?????????? ? ??????? ????????? ?????????? ? ??????
								  ,NumRows //  - ?-?? ?????
								  ,pwcharrFileNames1 //??????? ? ???????? ?????????? - ??????? nBuffCols x lenName
								  ,lenName // ???????????? ????? ????? ??????????
								  ,NumCols-1  // ????? ?????????? ?? ??? X
								  ,j  // ????? ?????????? ?? ??? Y
								  ,0.01//  ??????? ?? ??? X
								  ,1  // ??????? ?? ??? Y
								   );
	 }

  delete []parrBuff;
  delete []pwcharrFileNames1;


}




		 ///


 void TLearnShellTraj::calcSobstChislaPoOsiam(double *arrElK_XY0
 ,double *pvalLambX ,double *pvalLambY)
 {

 //   arrF - ??????? ??????????? ????????  ?????? ???????
	double arrF[4] = {0.} , arrMtrxLamb[4] = {0.};
	CalcProperVectors2(arrElK_XY0, arrF , arrMtrxLamb) ;

	 ///

	 // ?????? ???? ??????????? ? ???????? ?????????? ??????????? ?????

	   double det = arrF[0] * arrF[3]- arrF[1] * arrF[2] ;
	   if (det < 0.)
	   {
		double temp = arrMtrxLamb[3];
		arrMtrxLamb[3] = arrMtrxLamb[0];
		arrMtrxLamb[0] = temp;
	   }
	   *pvalLambX =  arrMtrxLamb[0];
	   *pvalLambY =  arrMtrxLamb[3];
 }


 //-----------------------------------------------------------------
double TLearnShellTraj::findOptimalTab76_No2_Cx(double *arrTab, int numRows, int numCols
   ,double valDistLow,double valDistUp,double valHLow,double valHUp
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter, int *pnumPoints)
{

	TLearnShellTraj ShellTrajCur = *this;
	//double valy = -1.;
	const int J0 = 0;

	double valFGr0 = ShellTrajCur.calcFGr_ForCx_Tab76_No2(arrTab,  numRows,  numCols
	,valDistLow, valDistUp, valHLow, valHUp, pnumPoints) ;
    if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}


	  for (int i =0;i < NUmIter; i++){


		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ (*pnumPoints));

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[100] = {0.};
	  int NumWorkingPoints = 27;
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < NumWorkingPoints; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnCx.Points[J0 +j].Y += 0.0001;
		double valy1 = ShellTrajCur1.calcFGr_ForCx_Tab76_No2(arrTab,  numRows,  numCols
	,valDistLow, valDistUp, valHLow, valHUp, pnumPoints) ;

		 arrGrad[j] = (valy1 - valFGr0)/ 0.0001;
	   }
	   ///

	   // ????? ????
	   double step = 0.002;
	   double valF = 1000000000000., valF1 = 10000000000000.;
	   int n =0;
	   for ( n =0; n < 40; n++)
	   {

		 step = step * 0.66;
		 double arrT0[100] = {0.};
		 MatrxMultScalar(arrGrad, 1, NumWorkingPoints, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 for (int k =0; k < NumWorkingPoints; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnCx.Points[J0 +k].Y  -=   arrT0[ k];
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForCx_Tab76_No2(arrTab,  numRows,  numCols
	,valDistLow, valDistUp, valHLow, valHUp, pnumPoints) ;
			 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}

	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;


	 return valFGr0;
}

//-------------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------------
double TLearnShellTraj::calcFGr_ForCx_Tab76_No2(double *arrTab, int numRows, int numCols
					   ,double valDistLow,double valDistUp,double valHLow,double valHUp, int *pnumPoints)
{
	TEnvironment Environment(0., M_PI/2., 0.);
	double sum = 0.;
	*pnumPoints = 0;
	for (int i = 0; i < numRows; i++)
	{

	  if(

	  ((arrTab [i *numCols +5] - valDistLow)* (arrTab [i *numCols +5] - valDistUp) >0.)
	  ||
	  ((arrTab [i *numCols +3] - valHLow)* (arrTab [i *numCols +3] - valHUp) >0.)
	  )
	  {
	  continue;
	  }
	   (*pnumPoints)++;
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody
		, arrTab [i *numCols +1] , M_PI/2.);


		double valDHoriz = -1.;
		ShellTrajCur.fncMovePhasVector( Environment, 0.001,	arrTab [i *numCols]);
		double arrDelS[2] = {0.};
		MtrxMinusMatrx(ShellTrajCur.marrStrSK_VS, &(arrTab [i *numCols +2]),1, 2, arrDelS);
		double delD = ScalProduct(arrDelS , arrDelS, 2) ;

		sum = sum + delD ;//+ delVk * delVk + delTetk * delTetk  + delH*delH;
	}
	return sum;
}

//------------------------------------------------------------
//-----------------------------------------------------------------
double TLearnShellTraj::findOptimalTab76_No2_MxOmegax(double *arrTab, int numRows, int numCols
   ,double valDistLow,double valDistUp,double valHLow,double valHUp
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter, int *pnumPoints)
{
	TLearnShellTraj ShellTrajCur = *this;
	const int J0 = 0;
	double valFGr0 = ShellTrajCur.calcFGr_ForKnm_Tab76_No2(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
	if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}

	int NumWorkingPoints = 27;
	for (int i =0;i < NUmIter; i++)
	{
		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ (*pnumPoints));

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[100] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < NumWorkingPoints; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +j].Y += 0.01;
		double valy1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No2(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 arrGrad[j] = (valy1 - valFGr0)/ 0.01;
	   }
	   NormalizeVect(arrGrad, NumWorkingPoints)  ;
	   ///

	   // ????? ????
	   double step = 0.002;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {
	   if (n ==39)
	   {
		 int iii =0;
	   }
		  step = step * 0.66;
		 double arrT0[100] = {0.};
		 MatrxMultScalar(arrGrad, 1, NumWorkingPoints, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 bool bcontinue = false;
		 for (int k =0; k < NumWorkingPoints; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +k].Y  -=   arrT0[ k];
		   if (ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +k].Y <=0.)
		   {
			 bcontinue = true;
			 break;
		   }
		 }

		 if (bcontinue)
		 {
           continue;
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No2(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;

	 return valFGr0;

}

//-------------------------------------------------------------------------------------

 //-------------------------------------------------------------------------------------
double TLearnShellTraj::calcFGr_ForKnm_Tab76_No2(double *arrTab, int numRows, int numCols
					   ,double valDistLow,double valDistUp,double valHLow,double valHUp, int *pnumPoints)
{
	TEnvironment Environment(0., M_PI/2., 0.);
	double sum = 0.;
	*pnumPoints = 0;
	for (int i = 0; i < numRows; i++)
	{

	  if(

	  ((arrTab [i *numCols +5] - valDistLow)* (arrTab [i *numCols +5] - valDistUp) >0.)
	  ||
	  ((arrTab [i *numCols +3] - valHLow)* (arrTab [i *numCols +3] - valHUp) >0.)
	  )
	  {
	  continue;
	  }
	   (*pnumPoints)++;
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody
		, arrTab [i *numCols +1] , M_PI/2.);


		double valDHoriz = -1.;
		ShellTrajCur.fncMovePhasVector( Environment, 0.001,	arrTab [i *numCols]);
		double delZ = ShellTrajCur.marrStrSK_VS[2] - arrTab[ i *numCols + 4] ;
		sum = sum + delZ* delZ;

	}
	return sum;
}
//-----------------------------------------------------------------

double TLearnShellTraj::findOptimalTab76_No2_Knm(double *arrTab, int numRows, int numCols
   ,double valDistLow,double valDistUp,double valHLow,double valHUp
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter, int *pnumPoints)
{
	TLearnShellTraj ShellTrajCur = *this;
	const int J0 = 0;
	double valFGr0 = ShellTrajCur.calcFGr_ForKnm_Tab76_No2(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);

	if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}


	int NumWorkingPoints = 27;
	for (int i =0;i < NUmIter; i++)
	{
		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ (*pnumPoints));

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[100] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < NumWorkingPoints; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +j].Y += 0.01;
		double valy1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No2(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 arrGrad[j] = (valy1 - valFGr0)/ 0.01;
	   }
	   NormalizeVect(arrGrad, NumWorkingPoints)  ;
	   ///

	   // ????? ????
	   double step = 0.005;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {

		  step = step * 0.66;
		 double arrT0[100] = {0.};
		 MatrxMultScalar(arrGrad, 1, NumWorkingPoints, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 bool bcontinue = false;
		 for (int k =0; k < NumWorkingPoints; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +k].Y  -=   arrT0[ k];
		   if (ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +k].Y <=0.)
		   {
		   //	 bcontinue = true;
			// break;
			int iii =0;
		   }
		 }

		// if (bcontinue)
		// {
		 //  continue;
		 //}
		 valF1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No2(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;

	 return valFGr0;

}

//----------------------------------------------
//----------------------------------------------
//----------------------------------------------
double TLearnShellTraj::findOptimalTab76_No3_Cx(double *arrTab, int numRows, int numCols
   ,double valDistLow,double valDistUp,double valHLow,double valHUp
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter, int *pnumPoints)
{

	TLearnShellTraj ShellTrajCur = *this;
	//double valy = -1.;
	const int J0 = 0;

	double valFGr0 = ShellTrajCur.calcFGr_ForCx_Tab76_No3(arrTab,  numRows,  numCols
	,valDistLow, valDistUp, valHLow, valHUp, pnumPoints) ;
    if (0 == NUmIter)
	{
	  *pplnCx = mLearnShellBody.mplnCx;
	  *pplnKnm = mLearnShellBody.mplnKnm;
	  *pplnMxOmx = mLearnShellBody.mplnMxOmx;
	  return valFGr0;
	}

	for (int i =0;i < NUmIter; i++){


		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ (*pnumPoints));

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[100] = {0.};
	  int NumWorkingPoints = 27;
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < NumWorkingPoints; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnCx.Points[J0 +j].Y += 0.0001;
		double valy1 = ShellTrajCur1.calcFGr_ForCx_Tab76_No3(arrTab,  numRows,  numCols
	,valDistLow, valDistUp, valHLow, valHUp, pnumPoints) ;

		 arrGrad[j] = (valy1 - valFGr0)/ 0.0001;
	   }
	   ///

	   // ????? ????
	   double step = 0.002;
	   double valF = 1000000000000., valF1 = 10000000000000.;
	   int n =0;
	   for ( n =0; n < 40; n++)
	   {

		 step = step * 0.66;
		 double arrT0[100] = {0.};
		 MatrxMultScalar(arrGrad, 1, NumWorkingPoints, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 for (int k =0; k < NumWorkingPoints; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnCx.Points[J0 +k].Y  -=   arrT0[ k];
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForCx_Tab76_No3(arrTab,  numRows,  numCols
	,valDistLow, valDistUp, valHLow, valHUp, pnumPoints) ;
			 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;


	 return valFGr0;
}

//-------------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------------
double TLearnShellTraj::calcFGr_ForCx_Tab76_No3(double *arrTab, int numRows, int numCols
					   ,double valDistLow,double valDistUp,double valHLow,double valHUp, int *pnumPoints)
{
	TEnvironment Environment(0., M_PI/2., 0.);
	double sum = 0.;
	*pnumPoints = 0;
	for (int i = 0; i < numRows; i++)
	{

	  if(

	  ((arrTab [i *numCols +4] - valDistLow)* (arrTab [i *numCols +4] - valDistUp) >0.)
	  ||
	  ((arrTab [i *numCols +3] - valHLow)* (arrTab [i *numCols +3] - valHUp) >0.)
	  )
	  {
	  continue;
	  }
	   (*pnumPoints)++;
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody
		, arrTab [i *numCols +1] , M_PI/2.);


		double valDHoriz = -1.;
		ShellTrajCur.fncMovePhasVector( Environment, 0.001,	arrTab [i *numCols]);
		double *ps = ShellTrajCur.marrStrSK_VS;
		double temp = sqrt(ps[0] * ps[0] + ps[2] * ps[2]);
		double delD = (ps[1]- arrTab [i *numCols +3]) * (ps[1]- arrTab [i *numCols +3])
		+ (temp - arrTab [i *numCols +2] )* (temp - arrTab [i *numCols +2] );
		sum = sum + delD ;
	}
	return sum;
}

//------------------------------------------------------------
//-----------------------------------------------------------------
double TLearnShellTraj::findOptimalTab76_No3_MxOmegax(double *arrTab, int numRows, int numCols
   ,double valDistLow,double valDistUp,double valHLow,double valHUp
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter, int *pnumPoints)
{
 /*	TLearnShellTraj ShellTrajCur = *this;
	const int J0 = 0;
	double valFGr0 = ShellTrajCur.calcFGr_ForKnm_Tab76_No3(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
	int NumWorkingPoints = 27;
	for (int i =0;i < NUmIter; i++)
	{
		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ (*pnumPoints));

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[100] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < NumWorkingPoints; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +j].Y += 0.01;
		double valy1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No3(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 arrGrad[j] = (valy1 - valFGr0)/ 0.01;
	   }
	   NormalizeVect(arrGrad, NumWorkingPoints)  ;
	   ///

	   // ????? ????
	   double step = 0.005;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {
	   if (n ==39)
	   {
		 int iii =0;
	   }
		  step = step * 0.66;
		 double arrT0[100] = {0.};
		 MatrxMultScalar(arrGrad, 1, NumWorkingPoints, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 bool bcontinue = false;
		 for (int k =0; k < NumWorkingPoints; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +k].Y  -=   arrT0[ k];
		   if (ShellTrajCur1.mLearnShellBody.mplnMxOmx.Points[J0 +k].Y <=0.)
		   {
			 bcontinue = true;
			 break;
		   }
		 }

		 if (bcontinue)
		 {
           continue;
		 }
		 valF1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No3(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;

	 return valFGr0;
	*/
	return 0.;
}

//-------------------------------------------------------------------------------------

 //-------------------------------------------------------------------------------------
double TLearnShellTraj::calcFGr_ForKnm_Tab76_No3(double *arrTab, int numRows, int numCols
					   ,double valDistLow,double valDistUp,double valHLow,double valHUp, int *pnumPoints)
{
 /*	TEnvironment Environment(0., M_PI/2., 0.);
	double sum = 0.;
	*pnumPoints = 0;
	for (int i = 0; i < numRows; i++)
	{

	  if(

	  ((arrTab [i *numCols +5] - valDistLow)* (arrTab [i *numCols +5] - valDistUp) >0.)
	  ||
	  ((arrTab [i *numCols +3] - valHLow)* (arrTab [i *numCols +3] - valHUp) >0.)
	  )
	  {
	  continue;
	  }
	   (*pnumPoints)++;
		double arrVesselVelocity[3] = {0.};
		TLearnShellTraj ShellTrajCur(arrVesselVelocity, mLearnShellBody
		, arrTab [i *numCols +1] , M_PI/2.);


		double valDHoriz = -1.;
		ShellTrajCur.fncMovePhasVector( Environment, 0.001,	arrTab [i *numCols]);
		double delZ = ShellTrajCur.marrStrSK_VS[2] - arrTab[ i *numCols + 4] ;
		sum = sum + delZ* delZ;

	}
	return sum; */
		return 0.;

}
//-----------------------------------------------------------------

double TLearnShellTraj::findOptimalTab76_No3_Knm(double *arrTab, int numRows, int numCols
   ,double valDistLow,double valDistUp,double valHLow,double valHUp
   , TURPolyLine *pplnCx, TURPolyLine *pplnKnm, TURPolyLine *pplnMxOmx, double *arrDisp
   , double *parrObjFnc, const int NUmIter, int *pnumPoints)
{
	/*
	TLearnShellTraj ShellTrajCur = *this;
	const int J0 = 0;
	double valFGr0 = ShellTrajCur.calcFGr_ForKnm_Tab76_No3(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
	int NumWorkingPoints = 27;
	for (int i =0;i < NUmIter; i++)
	{
		if (parrObjFnc)
		{
		parrObjFnc[i] =  sqrt(valFGr0/ (*pnumPoints));

		}
	   bool bend = true;
	  // ?????????? ?????????
	  double arrGrad[100] = {0.};
	 // valy =ShellTrajCur.calcFGr_ForCx_76(arrTab,  numRows,  numCols) ;
	  for (int j = 0; j < NumWorkingPoints; j++)
	  {
		TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +j].Y += 0.01;
		double valy1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No3(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 arrGrad[j] = (valy1 - valFGr0)/ 0.01;
	   }
	   NormalizeVect(arrGrad, NumWorkingPoints)  ;
	   ///

	   // ????? ????
	   double step = 0.005;
	   double valF = 1000000000000., valF1 = 10000000000000.;;
	   for (int n =0; n < 40; n++)
	   {

		  step = step * 0.66;
		 double arrT0[100] = {0.};
		 MatrxMultScalar(arrGrad, 1, NumWorkingPoints, step,arrT0);
		 TLearnShellTraj ShellTrajCur1 = ShellTrajCur;
		 bool bcontinue = false;
		 for (int k =0; k < NumWorkingPoints; k++)
		 {
		   ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +k].Y  -=   arrT0[ k];
		   if (ShellTrajCur1.mLearnShellBody.mplnKnm.Points[J0 +k].Y <=0.)
		   {
		   //	 bcontinue = true;
			// break;
			int iii =0;
		   }
		 }

		// if (bcontinue)
		// {
		 //  continue;
		 //}
		 valF1 =ShellTrajCur1.calcFGr_ForKnm_Tab76_No3(arrTab, numRows,  numCols
					 ,valDistLow, valDistUp, valHLow, valHUp, pnumPoints);
		 if (valF1 < valFGr0)
		 {
			ShellTrajCur = ShellTrajCur1;
			valFGr0 =  valF1;
			bend = false;
			break;
		 }
	   }
	   if (bend)
	   {
		break;
	   }

	}
	*pplnCx = ShellTrajCur.mLearnShellBody.mplnCx;
	*pplnKnm = ShellTrajCur.mLearnShellBody.mplnKnm;
	*pplnMxOmx = ShellTrajCur.mLearnShellBody.mplnMxOmx;

	 return valFGr0;
	*/
		return 0.;

}



#pragma package(smart_init)


