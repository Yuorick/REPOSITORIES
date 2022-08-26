//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include  <math.h>
#include "Plane.h"
#include "MatrixProccess.h"
#include "UrPointXY.h"

//--------------------------------------------------------------------------------------
 TPlane ::TPlane()
{
	double arrT [9] = {0.};
	memcpy(marrS0, arrT, 3 * sizeof(double));
	memcpy(marrF, arrT, 9 * sizeof(double));
	
}

 // �������� ������������
 TPlane &TPlane::operator=(const TPlane  &R)
 {
	memcpy(marrS0, R.marrS0, 3 * sizeof(double));
	memcpy(marrF, R.marrF, 9 * sizeof(double));
	return *this ;
 }

 // ����������� �����������
 TPlane::TPlane (const TPlane &R)
 {
	memcpy(marrS0, R.marrS0, 3 * sizeof(double));
	memcpy(marrF, R.marrF, 9 * sizeof(double));
 }


 // ����� ������
TPlane ::TPlane( double* arrS0, double* arrF)
{
	memcpy(marrS0, arrS0, 3 * sizeof(double));
	memcpy(marrF, arrF, 9 * sizeof(double));

}

bool TPlane::findIntersectingPoint_with_Line(double *arrPosWorking, double *arrVeloWorking, TURPointXY *ppntIntersect)
{
	// 1. ���������� ����� ��������� ������ � ��������� � �������� �������� ���������
	double arrN[3] = {0.} ; // ������� ������� � ���������
	fillNormalVect(arrN);
	///
	double arrPointINtersect_PrSK[3] = {0.}; // ��� ����� ������ ����� ����������� ��������� ���� ��
	double arrDelta [3] ={0.};
	MtrxMinusMatrx(arrPosWorking, marrS0,3, 1, arrDelta);
	double val_t0 = ScalProduct(arrN ,  arrDelta, 3) ;
	double val_t = ScalProduct(arrN , arrVeloWorking, 3) ;
	if (fabs(val_t0) < 0.00000001)
	{  // ����� arrPosWorking ����� �� ��������� ��� ������ �������� ���������� ��������� � ����������� ���
		memcpy(arrPointINtersect_PrSK, arrPosWorking, 3 * sizeof(double));
		return false;
	}
	else
	{
		if ( fabs(val_t )/ Norm3(arrVeloWorking) < 0.0000001)
		{ // ����� �� ��������� �� �����, ������ ����������� ���������
		return false;
		}
		double val_par =  -val_t0 / val_t ;
		double arrt0[3] = {0.};
		MatrxMultScalar(arrVeloWorking, 1, 3, val_par,arrt0);
		MtrxSumMatrx(arrPosWorking, arrt0, 3, 1, arrPointINtersect_PrSK) ;
	}

	double arrPointINtersect_SKP [3] ={0.};
	transform_xyzSSK_to_xyzSKP(arrPointINtersect_PrSK, arrPointINtersect_SKP);
	(*ppntIntersect).X = arrPointINtersect_SKP[0] ;
	(*ppntIntersect).Y = arrPointINtersect_SKP[1 ];
	if (fabs( arrPointINtersect_SKP[2]) > 0.00000001)
	{
		ShowMessage(L" ERROR TPlane::findIntersectingPoint_with_Line");
	}
	return true;
}
//-----------------------------------------------
// ������������ ������� ������� - 3-�� ������� �������  marrF
//
void TPlane::fillNormalVect(double *arrN)
{
 arrN [0]	=  marrF [ 2];
 arrN [1]	=  marrF [ 5];
 arrN [2]	=  marrF [ 8];
}

//----------------------------------------------
// �������������� ������� �� �������� ������� �������� ��������� � �������� ��������� ���������
//INPUT:
// arrPoint_PrSK[3]
// OUTPUT:
// arrPoint_SKP[3]
void TPlane::transform_xyzSSK_to_xyzSKP(double *arrPoint_PrSK, double *arrPoint_SKP)
{
	double arrDelta [3] ={0.};
	MtrxMinusMatrx(arrPoint_PrSK, marrS0,3, 1, arrDelta);
	MtrxTranspMultMatrx(marrF,3, 3, arrDelta, 1, arrPoint_SKP) ;
}
//----------------------------------------------
// �������������� ������� �� �������� ��������� ��������� � �������� ������� �������� ���������
//INPUT:
// arrPoint_PrSK[3]
// OUTPUT:
// arrPoint_SKP[3]
void TPlane::transform_xyzSKP_to_xyzSSK(double *arrPoint_SKP, double *arrPoint_PrSK)
{
	double arrDelta [3] ={0.};
	MtrxMultMatrx(marrF,3, 3, arrPoint_SKP, 1, arrDelta) ;
	MtrxSumMatrx(arrDelta, marrS0,3, 1, arrPoint_PrSK);
}
#pragma package(smart_init)
