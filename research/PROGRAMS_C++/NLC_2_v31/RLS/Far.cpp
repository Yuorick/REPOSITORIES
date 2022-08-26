//---------------------------------------------------------------------------


#pragma hdrstop
#include "Far.h"

 #include "Faceta.h"
//---------------------------------------------------------------------------

#include <vcl.h>


#include <math.h>
#include "Comp.h"
#include <stdio.h>
#include <stdlib.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"

TFar::TFar()
{
 // ���������� ������� (�����)
  m_N = 4;
 // ���������� ����� ��������
   m_D = 32.8;
 // ����� �����
   mLambda = 3.15;;
 // ������ �����
  mFaceta = TFaceta();

}
// ����������� �����������
TFar::TFar (const TFar &R2)
 {
 // ���������� ������� (�����)
  m_N = R2.m_N;
 // ���������� ����� ��������
  m_D = R2.m_D;
 // ����� �����
  mLambda = R2.mLambda;

  mFaceta = R2.mFaceta;

 }

// ����� ������
 __fastcall TFar::TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta)
 {
	 m_N = N ;
	 m_D= D;
	 mLambda = Lambda ;
	 mFaceta = Faceta;

 }

 // ����� ������
 __fastcall TFar::TFar(const int N)
 {

 // ���������� ������� (�����)
  m_N = N;
 // ���������� ����� ��������
   m_D = 32.8;
 // ����� �����
   mLambda = 3.15;;
 // ������ �����
  mFaceta = TFaceta();

 }
// {

// }

 // �������� ������������
 TFar TFar::operator=(TFar  R2)
{
 // ���������� ������� (�����)
  m_N = R2.m_N;
 // ���������� ����� ��������
  m_D = R2.m_D;
 // ����� �����
  mLambda = R2.mLambda;

  mFaceta = R2.mFaceta;

  return *this ;
}



#pragma package(smart_init)
