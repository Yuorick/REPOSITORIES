//---------------------------------------------------------------------------


#pragma hdrstop

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TruncatedConeCircle.h"
#include "Plane.h"


  TTruncatedConeCircle ::TTruncatedConeCircle()
{

	mR = 0. ;
	mr = 0. ;
	mH = 0.;
   //	mM = 0.;


}
TTruncatedConeCircle ::TTruncatedConeCircle(const TPlane Plane, const double R, const double r,const double H, const double M):TSimpleBody_3D( Plane, M)
{
 mR = R;
 mr = r;
 mH = H;
// mM = M;

}

// ����������� �����������
 TTruncatedConeCircle ::TTruncatedConeCircle (const TTruncatedConeCircle &R)
 {
	mR  = R.mR;
	mr = R.mr ;
	mH = R.mH;
	mM = R.mM   ;
	mPlane = R.mPlane;
 }
 // �������� ������������
 TTruncatedConeCircle TTruncatedConeCircle::operator=(TTruncatedConeCircle  R)
 {
	mR  = R.mR;
	mr = R.mr ;
	mH = R.mH;
	mM = R.mM   ;
	mPlane = R.mPlane;
	return *this ;
  }


// ���������� ������
double TTruncatedConeCircle::calcCapacity()
{
 return  M_PI * mH * (mR * mR + mR * mr + mr * mr)/3.;
}

// ���������� ��������� ������ �������
// ������ ���������� �� ����������� ��������, ��� 46
void TTruncatedConeCircle::calcCentreOfGravity(double *arrCentreGrav)
{
  arrCentreGrav [0] = mH * (mR * mR + 2. * mR * mr + 3. * mr * mr)/ (mR * mR + mR * mr + mr * mr) / 4.;
  arrCentreGrav [1] = 0.;
  arrCentreGrav [2] = 0.;
}

// ���������� ������� �������� �������
// ��� ������������ �������
// ������������ ��� �������� ��������� �� �������
// ������������ ���� X � Z ���������� ������� ��������
//  ��������� ����� ������������� �� ������������ (��������).
// �� ���� ����������� ��������� ����� ������.
// ����������� ���������� ������ ������� ����������, �������� � ������ �������.
// ������������ ������� ������� �������� � ������ �������.
// �� ������� �������� ��������� ������ ������� ��������� ������ ������������ ����  X � Z
void TTruncatedConeCircle::calcInertiaMtrx(double *arrInertMtrx)
{
  // 1. ���������� ������� ������� ������������ ��� X
  // ������ ���������� �� ����������� ��������, ��� 46
  memset(arrInertMtrx, 0, 9 * sizeof(double));
  double valMass = 1.;
  arrInertMtrx[0] = 3. * valMass * (mR*mR*mR*mR*mR - mr*mr*mr*mr*mr)/ (mR*mR*mR - mr*mr*mr) / 10.;
  ////


  double arrCentreGravLittle [3] = {0.}  ,arrCentreGravBig [3] = {0.};
  // 2. ���������� ����� ������ � �������� �������
   double valHLittle = mH * mr / (mR - mr);
	double valHBig =  valHLittle + mH;
	///


   //  3. ����������  ��������� ������ ���� ������ ������
   calcCentreOfGravity_For_ConeCircle(valHLittle, arrCentreGravLittle) ;
   arrCentreGravLittle [0] += mH;
   ///

   // 4. ����������  ��������� ������ ���� �������� ������
   calcCentreOfGravity_For_ConeCircle(valHBig, arrCentreGravBig) ;
   ///


   //       5.   ���������� ����� ������ ������
  double valMLittle = valMass / calcCapacity() * calcCapacity_For_ConeCircle(mr, valHLittle);
   ///


   //    6. ���������� ����� , ��������  ������
  double valMBig = valMass + valMLittle;
  ///



   // 7. ����������  ��������� ������ ���� ����������  ������
   double arrCentreGrav[3] = {0.};
   calcCentreOfGravity(arrCentreGrav)  ;
   ///

   // 8.���������� ����� �� ������ � �������� �������
   double valDelRLittle =  arrCentreGravLittle[0] - arrCentreGravBig[0];
   ///

   // 9.���������� ����� �� ����������  � �������� �������
   double valDelR =  arrCentreGrav[0] - arrCentreGravBig[0];
   ///

   // ������� ��������

   arrInertMtrx[4] = calcJyy_For_ConeCircle(mR, valHBig,valMBig) - calcJyy_For_ConeCircle(mr, valHLittle,valMLittle)
				  -valMLittle * valDelRLittle * valDelRLittle - valMass * valDelR * valDelR;

   arrInertMtrx[8] = arrInertMtrx[4];


}

// ���������� ������� ������� ������� ��������� ������ ������������ ��� ���������
// ������ ���������� �� ����������� ��������, ��� 237
double TTruncatedConeCircle::calcJxx_For_ConeCircle(const double VAlR, const double VAlM) // +
{
return 3. * VAlM * VAlR * VAlR /10.;

}

// ���������� ������� ������� ������� ��������� ������ ������������ ���
// ���������������� ��� ���������  � ����������� ����� ����� ����
// ������ ���������� �� ����������� ��������, ��� 237
double TTruncatedConeCircle::calcJyy_For_ConeCircle(const double VAlR, const double VAlH, const double VAlM)  // +
{
return 3. * VAlM * (VAlR * VAlR   +  VAlH * VAlH / 4.)/20.;

}
// ���������� ������  ������� ��������� ������
double TTruncatedConeCircle::calcCapacity_For_ConeCircle(const double VAlR, const double VAlH)  // +
{
  return  M_PI * VAlR *VAlR * VAlH / 3.;
}

// ���������� ������� ������ ���� ������� ��������� ������ ������������ ���
// ���������������� ��� ���������  � ����������� ����� ����� ����
// ������ ���������� �� ����������� ��������, ��� 47
void TTruncatedConeCircle::calcCentreOfGravity_For_ConeCircle(const double VAlH
, double *arrCentreGrav)   // +
{
  arrCentreGrav [0] = VAlH / 4. ;
  arrCentreGrav [1] = 0.;
  arrCentreGrav [2] = 0.;
}
// ���������� ������� �������� �������  ��� ������� ��������� ������  ������������ ������ �������
void TTruncatedConeCircle::calcInertiaMtrx_For_ConeCircle(const double VAlR, const double VAlH, const double VAlM  // +
,double *arrInertMtrx)
{
   memset(arrInertMtrx, 0, 9 * sizeof(double));
   arrInertMtrx [0] = calcJxx_For_ConeCircle(VAlR,VAlM);
   arrInertMtrx [4] = calcJyy_For_ConeCircle( VAlR,VAlH, VAlM);
   arrInertMtrx [8] = arrInertMtrx [4] ;

}
#pragma package(smart_init)