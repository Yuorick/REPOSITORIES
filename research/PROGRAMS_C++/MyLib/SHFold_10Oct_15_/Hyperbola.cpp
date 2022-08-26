//---------------------------------------------------------------------------


#pragma hdrstop

#include "Hyperbola.h"

#include <vcl.h>
#include "ArcEllipse.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include <stdio.h>
 #include <math.h>
 #include "URPolyLineZ.h"
 #include "Equations.h"
 #include "URPolygon.h"
 #include "Comp.h"



//---------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
 THyperbola ::THyperbola()
{
	ma = 4.;
	mb = 3. ;

}

 // ����� ������
THyperbola :: THyperbola( const double a, const double b)

 {

   ma = a;
   mb = b ;


 }


 // �������� ������������
 THyperbola THyperbola::operator=(THyperbola  R)
 {
	ma  = R.ma ;
	mb  = R.mb;
	return *this ;
 }

 // ����������� �����������
 THyperbola::THyperbola (const THyperbola &R)
 {
   ma  = R.ma ;
   mb  = R.mb;

 }

double  THyperbola::fncFocus ()
{
	return sqrt( ma * ma + mb * mb);
}


void THyperbola::ShowMe(wchar_t *FileName, const int iNumPoints0, const double VAlDiap)
{
   TURPolyLine line( *this,  iNumPoints0,  VAlDiap) ;


   line.WriteSetSHPFiles(FileName,&line, 1 ) ;
}

TURPointXY  THyperbola::ProjectPointOnHyperbola (TURPointXY  pntInp)
{

  double xZv = pntInp.X / ma;
  double yZv = pntInp.Y / mb;

  // ���������� �������� �����  pntTemp �� ��������� ���������

  const double a =  xZv* xZv - yZv * yZv;
  const double b = 2. *  (xZv* xZv -+ yZv * yZv);
  const double c  = a -1.;
  TComp cmpx1,cmpx2;
  // ������� ����������� ��������� a*x*x + b*x + c =0
// ����������
// 0 - 2 �������������� ��������� �����
// 1 - 2 �������������� ������� �����
// 2- ������� �� ������� ���� ���� ������� ������ a!= 0, c=0
// 3 - 2 ���������� ����������� �����
// 4 -  1 �������������� ������ (a =0)
// 5  - �������������� a=b=0, c!=0
// 6 - ��������� )a=b=c=0
   SolvEq2( a, b, c,cmpx1,cmpx2);
  double q1 = cmpx1.m_Re;
  double q2 = cmpx2.m_Re;
  double q = q1;
  if (fncTemp( q2,  xZv, yZv, pntInp) < fncTemp( q1,  xZv, yZv, pntInp))
  {
	q = q2;
  }

  TURPointXY  pntRez( xZv / ( 1. -q) * ma, yZv / ( 1. +q) * mb);
  return  pntRez;


}

double THyperbola::fncTemp(double q, double xZv, double yZv, TURPointXY  pntInp)
{
 double x =  xZv / ( 1. -q) * ma;
 double y =  yZv / ( 1. +q) * mb;
 return (x - pntInp.X)*(x - pntInp.X) + (y - pntInp.Y)*(y - pntInp.Y);

}


// ���������� ����� ������� ������ ����������� �� ��������� ������  pntInp
// � ���������
// pntInp -  �������� ����� �� ���������
// OUTPUT:
// pntOut0,pntOut1 -   ���� �������
// ��������� �-�� ����� �������
//
int THyperbola::findTangencyPoints(TURPointXY  pntInp, TURPointXY  &pntOut0, TURPointXY  &pntOut1)
{

  // �������� � ������ "��������� ���������"
  TURPointXY  pntInp0( pntInp.X/ ma, pntInp.Y/ mb);
  ///
  if ((pntInp0.X * pntInp0.X - pntInp0.Y * pntInp0.Y -1.) > 0.)
  {
	return 0;
  }
  // ���������� ����� ������� ��������� ��������� x*x - y*y =1
  double vala = pntInp0.Y *  pntInp0.Y  -   pntInp0.X * pntInp0.X;
  double valb = 2. * pntInp0.X;
  double valc = -(1. +  pntInp0.Y *  pntInp0.Y);
  TComp x1,x2;
  // ������� ����������� ��������� a*x*x + b*x + c =0
// ����������
// 0 - 2 �������������� ��������� �����
// 1 - 2 �������������� ������� �����
// 2- ������� �� ������� ���� ���� ������� ������ a!= 0, c=0
// 3 - 2 ���������� ����������� �����
// 4 -  1 �������������� ������ (a =0)
// 5  - �������������� a=b=0, c!=0
// 6 - ��������� )a=b=c=0
int irez =  SolvEq2(vala, valb, valc,x1,x2) ;
switch (irez)
{
	case 0:
	case 2:
	 pntOut0.X =  x1.m_Re;
	 pntOut0.Y = sqrt( pntOut0.X * pntOut0.X  -1) ;
	 if(!chekRoot (pntInp0, pntOut0))
	 {
	   pntOut0.Y = -pntOut0.Y;
	 }
	 pntOut0.Y = pntOut0.Y * mb;
	 pntOut0.X = pntOut0.X * ma;


	 pntOut1.X =  x2.m_Re;
	 pntOut1.Y = sqrt( pntOut1.X * pntOut1.X  -1);
	 if(!chekRoot (pntInp0, pntOut1))
	 {
	   pntOut1.Y = -pntOut1.Y;
	 }
	 pntOut1.Y = pntOut1.Y * mb;
	 pntOut1.X = pntOut1.X * ma;
	 return 2;
   //	break;
	case 1:
     pntOut0.X =  x1.m_Re;
	 pntOut0.Y = sqrt( pntOut0.X * pntOut0.X  -1) * mb;
	 pntOut0.X = pntOut0.X * ma;
	 return 1;
	default:
	return 0;
   //	break;
}
}


//------------------------------------------------
// �������� ����� ����������� ���������
bool chekRoot (TURPointXY  pntInp0, TURPointXY pntOut)
{
	double b = ( pntOut.X -   pntInp0.X)/ ( pntOut.Y -   pntInp0.Y);
	if (fabs(pntOut.X - pntOut.Y * b) > 0.0000001)
	{
	  return false;
	}
	return true;
}

#pragma package(smart_init)