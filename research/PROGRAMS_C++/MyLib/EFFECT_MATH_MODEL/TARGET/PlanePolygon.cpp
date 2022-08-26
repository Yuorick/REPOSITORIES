//---------------------------------------------------------------------------


#pragma hdrstop
#include  <string.h>
#include "PlanePolygon.h"
#include "URPointXY.h"
#include "MatrixProccess.h"
#include "CalcCorMatrx.h"


//--------------------------------------------------------------------------------------
 TPlanePolygon ::TPlanePolygon()
{
 menumTypeTargBody = PLANE_BODY;
 mPlane = TPlane ();
 mPolygon = TURPolygon();

}

 // �������� ������������
 TPlanePolygon &TPlanePolygon::operator=( const TPlanePolygon  &R)
 {
	mPlane = R.mPlane ;
	mPolygon = R.mPolygon ;
	menumTypeTargBody = R.menumTypeTargBody;
	return *this ;
 }

 // ����������� �����������
 TPlanePolygon::TPlanePolygon (const TPlanePolygon &R)
 {
	mPlane = R.mPlane ;
	mPolygon = R.mPolygon ;
	menumTypeTargBody = R.menumTypeTargBody;

 }

 // ����� ������
TPlanePolygon :: TPlanePolygon( const TPlane Plane, const TURPolygon Polygon, const enumTypeTargBody EnumTypeTargBody)
{
	 mPlane = Plane ;
	 TURPolygon PolygonTemp  =Polygon;
	 double valS = PolygonTemp.calcVectSq();
	 if(valS < 0.)
	 {
	   PolygonTemp.flip();
     }
	 mPolygon = PolygonTemp ;
	 menumTypeTargBody = EnumTypeTargBody;
}


 // ����� ������  2
TPlanePolygon :: TPlanePolygon( double* arrS0, double* arrF, const TURPolygon Polygon, const enumTypeTargBody EnumTypeTargBody)

 {
	 memcpy(mPlane.marrS0, arrS0, 3 * sizeof(double));
	 memcpy(mPlane.marrF, arrF, 9 * sizeof(double));
	 mPolygon = Polygon ;
	 menumTypeTargBody = EnumTypeTargBody;
 }



 // ����������� ���������� ����� � �������� ��������
 // INPUT:
 // arrPos - ��������� ����� �����
 // arrVelo - ������ �������� (��������������� ������ ������������� ���������)
 // ���������� true ���� ���� ������������ � false  � ���� ������
 //
bool TPlanePolygon::isLineIntersectPlanePolygon(double *arrPos, double *arrVelo)
{
	// � ����������� �� ���� �������� �������� - PLANE_BODY ��� AXIALLY_SYMMETRICAL_BODY
	// ��������� �������� ���� �������� ����� ��� ��� ���� ( � ������ PLANE_BODY)
	// ��� ���������� � �������� ( � ������ AXIALLY_SYMMETRICAL_BODY)
	// ������� ��������� ���������� ��������� ������� ��������  arrVelo
	// ���� � ����� ���������� ������������ ������� ������� ��������
	// ��� ������  AXIALLY_SYMMETRICAL_BODY
	double arrPosWorking [3] ={0.}, arrVeloWorking [3] ={0.};
	memcpy(arrPosWorking, arrPos, 3 * sizeof(double));
	memcpy(arrVeloWorking, arrVelo, 3 * sizeof(double));
	if (menumTypeTargBody == AXIALLY_SYMMETRICAL_BODY)
	{
	 formMatrxE(3,   mPlane.marrF);

	 double arrMtrxPer[9] = {0.}, arr_v[2] = {0.};
	 arr_v[0] = arrVelo[1];
	 arr_v[1] = arrVelo[2];
	 double val_v = NormVect2(arr_v);
	 if (val_v < 0.00000001)
	 {
			 return false;
	 }
	 arr_v[0] = arr_v[0] / val_v ;
	 arr_v[1] = arr_v[1] / val_v ;
	 arrMtrxPer[0] = 1.;
	 arrMtrxPer[4] =  arr_v[1];
	 arrMtrxPer[5] = -arr_v[0];
	 arrMtrxPer[7] =  arr_v[0];
	 arrMtrxPer[8] =  arr_v[1];
	 double arrPosTemp[3] = {0.};
	 MtrxMinusMatrx(arrPos,   mPlane.marrS0,3,1,  arrPosTemp);

	 MtrxMultMatrx(arrMtrxPer,3, 3, arrPosTemp,1, arrPosWorking) ;
	 MtrxMultMatrx(arrMtrxPer,3, 3, arrVelo,1, arrVeloWorking) ;
	   mPlane.marrS0[0] = 0.;
	   mPlane.marrS0[1] = 0.;
	   mPlane.marrS0[2] = 0.;
	}

	TURPointXY pntIntersect;
	if( !mPlane.findIntersectingPoint_with_Line(arrPosWorking, arrVeloWorking, &pntIntersect))
	{
		return false;
	}

	return mPolygon.PtInPlygon( pntIntersect) ;

}

//---------------------------------------------------------------------------
 // �������� �������� ���� ��� �������� ��������  �� �������������� ���������
 // INPUT:
 // arrTargV[3] -  ������ �������� ���� � ���
 // arrVectMissV [3] - ������ ������������� ��������
 //
 //
TURPolygon TPlanePolygon::createShadowPlg_For_PlaneBody(double *arrTargV,  double *arrVectMissV)
{
	TURPolygon plgRez = mPolygon;
	double  arrMtrxPer_SkSK_GSK[9] = {0.}, arrL[9] = {0.}, arrTemp0[9] ={0.}, arrMtrxPer[9] = {0.};
	calcMatrxPer_from_SSK_To_DecartPrSK(arrTargV, arrMtrxPer_SkSK_GSK);
	arrL[0] = 1.;
	arrL[2] = - arrVectMissV[0] / arrVectMissV[2];
	arrL[4] = 1.;
	arrL[5] = - arrVectMissV[1] / arrVectMissV[2];

	MtrxMultMatrx(arrL ,3, 3, arrMtrxPer_SkSK_GSK,3, arrTemp0) ;
	MtrxMultMatrx(arrTemp0 ,3, 3, mPlane.marrF,3, arrMtrxPer) ;

	double arrS0[3] = {0.};
	MtrxMultMatrx(arrTemp0 ,3, 3, mPlane.marrS0,1, arrS0) ;

	for (int i =0; i < mPolygon.NumPoints; i++)
	{
	 double arrS[3] = {0.}, arrTemp1[3] = {0.}, arrTemp2[3] = {0.};
	 arrS[0] =  mPolygon.Points[i].X;
	 arrS[1] =  mPolygon.Points[i].Y;


	 MtrxMultMatrx(arrMtrxPer ,3, 3, arrS,1, arrTemp1) ;
	 MtrxSumMatrx(arrTemp1, arrS0,1, 3, arrTemp2) ;
	 plgRez.Points[i].X =  arrTemp2[0];
	 plgRez.Points[i].Y =  arrTemp2[1];

	}

	if (plgRez.calcVectSq()< 0.)
	{
	plgRez.flip();
	}
	return plgRez;
}

TURPolygon TPlanePolygon::createShadowPlg(double *arrTargV,  double *arrVectMissV)
{
	if ( menumTypeTargBody == PLANE_BODY)
	{
		return createShadowPlg_For_PlaneBody(arrTargV, arrVectMissV);
	}
	TPlanePolygon PlanePolygonCur = *this;
	PlanePolygonCur.menumTypeTargBody =  PLANE_BODY;
	// ������������ ��������� � ������� ����� ������� ������
 /*	double arrF[9] = {0.};

 	memcpy(arrF, mPlane.marrF, 9 * sizeof(double));
	double e1 = mPlane.marrF[0];
	double e2 = mPlane.marrF[3];
	double e3 = mPlane.marrF[6];
	arrF [1] = -e1 * e2;
	arrF [2] =  e2;
	arrF [4] = -e2 * e3;
	arrF [5] = -e1 * ( 1. + e2 * e3 - e3 * e3);
	arrF [7] = 1. - e3 * e3;
	arrF [8] = e1 * e2 * (e2 - e3);
	double temp = calcDet3( arrF);
	memcpy(PlanePolygonCur.mPlane.marrF,arrF , 9 * sizeof(double));  */

	memset(PlanePolygonCur.mPlane.marrF, 0, 9 * sizeof(double));
	PlanePolygonCur.mPlane.marrF [0] = 1.;
	PlanePolygonCur.mPlane.marrF [4] = 1.;
	PlanePolygonCur.mPlane.marrF [4] = 1.;
	return  PlanePolygonCur.createShadowPlg_For_PlaneBody(arrTargV, arrVectMissV);


}


 TURPolygon TPlanePolygon::createProjectionOfPolygon_To_CartinPlane(double *arrTargV, double *arrVectMissV)
 {


	// ���� ����� ��������� ��������
	double valAng = calcAngBetweenVect (arrTargV, arrVectMissV,3);


	// ���� ������� �������� �����������, �� �������� ������������������ ���� �� ��������� ���������
	// ����� ����� �������� � ���� ����� - ������ ��������� . ��� ������������ ������  fabs(cos(valAng)) < 0.999999999
	if((menumTypeTargBody == AXIALLY_SYMMETRICAL_BODY) &&(fabs(cos(valAng)) < 0.999999999))
	{

	 double arrBasis[9] ={0.}, arrBasis0[9] ={0.};
	 createOrthogBasis_dim3 (arrTargV,arrVectMissV,  arrBasis);
	 // �������� ��������� �����
	 for (int i = 0; i < 3; i++)
	 {
		double temp = arrBasis[i * 3 +1];
		arrBasis[i * 3 +1] = arrBasis[i * 3 + 2] ;
		arrBasis[i * 3 + 2] = - temp;
	 }

	 double  parrMtrxPer_GSK_SkSK[9] = {0.};
	 calcMatrxPer_from_DecartPrSK_To_SSK(arrTargV, parrMtrxPer_GSK_SkSK) ;
	 MtrxMultMatrx(parrMtrxPer_GSK_SkSK,3, 3, arrBasis,3, arrBasis0) ;
	 memcpy(mPlane.marrF,arrBasis0, 9 * sizeof(double));
	}
 return	createProjectionOf_PlaneBodyPolygon_To_CartinPlane(arrTargV, arrVectMissV) ;


 }

 //------------------------------------------------------------------------------------------
TURPolygon TPlanePolygon::createProjectionOf_PlaneBodyPolygon_To_CartinPlane(double *arrTargV, double *arrVectMissV)
{
	TURPolygon plgRez = mPolygon;
   //	int inumpoints = mPolygon.NumPoints;
   //	int inumparts =  mPolygon.NumParts;
   //	TURPolygon plgRez(  inumpoints);

	double  arrMtrxPer_SkSK_GSK[9] = {0.}, arrL[9] = {0.};
	// ������� �������� �� ���� ���� � ���
	calcMatrxPer_from_SSK_To_DecartPrSK(arrTargV, arrMtrxPer_SkSK_GSK);
	 // ������� �������� �� �� ��������� � ���
	MtrxMultMatrx(arrMtrxPer_SkSK_GSK,3, 3, mPlane.marrF,3, arrL) ;


	// ������� �������� �� ��� � �� ��������� ���������
	 double  arrMtrxPer_GSK_CartinSK[9] = {0.}; //
	 calcMatrxPer_from_DecartPrSK_To_SSK(arrVectMissV,  arrMtrxPer_GSK_CartinSK) ;
	 // ��������� ������� ��������
	 double arrMtrxPer_Rez[9] = {0.};
	 MtrxMultMatrx(arrMtrxPer_GSK_CartinSK, 3, 3, arrL, 3, arrMtrxPer_Rez) ;


	for (int   i =0; i < mPolygon.NumPoints; i++)
	{
	 double arrS[3] = {0.}, arrTemp1[3] = {0.};
	 arrS[0] =  mPolygon.Points[i].X;
	 arrS[1] =  mPolygon.Points[i].Y;


	 MtrxMultMatrx(arrMtrxPer_Rez ,3, 3, arrS,1, arrTemp1) ;

	 plgRez.Points[i].X =  arrTemp1[1];
	 plgRez.Points[i].Y =  arrTemp1[2];

	}

	if (plgRez.calcVectSq()< 0.)
	{
	plgRez.flip();
	}
  return plgRez;
}

//------------------------------------------------------------------------------------------
// ���������� ���������� �� ����� �� �������� ��������
double TPlanePolygon::calcDistanse(double *arrPosOtnSSK)
{
 // 1. ������ ���������� ����� �� ������� �������� ��������� ���� � ������� ��������� ���������
 double arrPosOtnSKP[3] = {0.};
 mPlane.transform_xyzSSK_to_xyzSKP(arrPosOtnSSK, arrPosOtnSKP) ;
 ///

 double valPerpendiculareLength = 0.; // ����� ���������������, ���������� �� ����� �� ���������
 double valProjectionLength = 0.; // ����������� ���������� �� ��������� �������������� �� ��������
								  // �������� � ���������
 double arrPerpendic[2] = {0.}; // ��������� ��������������



 switch(menumTypeTargBody)
 {
   case PLANE_BODY:
   valPerpendiculareLength = arrPosOtnSKP[2];
   arrPerpendic[0] = arrPosOtnSKP[0];
   arrPerpendic[1] = arrPosOtnSKP[1];
   break;

   case AXIALLY_SYMMETRICAL_BODY:
   arrPerpendic[0] = arrPosOtnSKP[0];
   arrPerpendic[1] = sqrt(arrPosOtnSKP[1] * arrPosOtnSKP[1] + arrPosOtnSKP[2] * arrPosOtnSKP[2]) ;

   break;

   default:
   break;
 };
 TURPointXY pnt(arrPerpendic[0], arrPerpendic[1]);
 valProjectionLength = mPolygon.calcDistFromPoint(pnt);

 return sqrt(valProjectionLength  * valProjectionLength  + valPerpendiculareLength * valPerpendiculareLength);
}

#pragma package(smart_init)
