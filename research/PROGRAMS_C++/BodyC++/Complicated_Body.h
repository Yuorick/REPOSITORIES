//---------------------------------------------------------------------------
 // ����� ��������� ������� �������� ����, ��������� �� ��������� ������� ���
 // ������� ���� ��� ��� ������ �������� �������� ��� ������ �������� ��������� ������
 // ��������� ������� ��� ����� ���� ��������� ��� �������������
 // MAX_POSSIBLE_QUANT_OF_SIMPLE_BODY - ��� ����������� ��������� �-�� ������� ���
 // marrCircleCylinder - ������ ���������
 // mQuantCircleCylinder - �-�� ���������
 // marrTruncatedConeCircle - ��������� ��������� �������
 // mQuantTruncatedConeCircle - ����������  ��������� �������
 // ������ ������� ���� ����� ���� ��������� ������� ���������
 // ������� ���� ����� ���� ��������� �������� ���������
 // �������� ������������ ��������� ������ ��������� ������� �������� ����
 // � �������� ���� ��������� � ����� ������ TSimpleBody_3D ��� ��������� mPlane
 // ��� ������� ��������� ����� ��������� �������� ��������� �������� ����
 // � ������ ��������� ������ ���� �������� ���������
 // ����� ������������ ��� ������ ������� �������� ������� �������� ����
 // � ������ ��������� ������ ������� ������ �������� ���� �
 // � ��������� �������� ��������� �������� ����
 //

#ifndef Complicated_BodyH
#define Complicated_BodyH
#include "SimpleBody_3D.h"
#include "CircleCylinder.h"
#include "TruncatedConeCircle.h"


#define MAX_POSSIBLE_QUANT_OF_SIMPLE_BODY 20
class TCircleCylinder;
class TTruncatedConeCircle;
class TComplicated_Body
{
public:

	  TCircleCylinder marrCircleCylinder[MAX_POSSIBLE_QUANT_OF_SIMPLE_BODY];
	  int mQuantCircleCylinder;
	  TTruncatedConeCircle marrTruncatedConeCircle[MAX_POSSIBLE_QUANT_OF_SIMPLE_BODY];
	  int mQuantTruncatedConeCircle;

	  // ���������� ������ ���� � ��������� ��������� ������� ���������
	  double marrSdvig[3];

	__fastcall TComplicated_Body() ;

__fastcall TComplicated_Body(const int  QuantCircleCylinder
  , TCircleCylinder *arrCircleCylinder, const int  QuantTruncatedConeCircle
  , TTruncatedConeCircle *arrTruncatedConeCircle);

  __fastcall TComplicated_Body(const double VAlM,const double VAlH1
  , const double VAlR1, const double VAlH2, const double VAlH3
  , const double VAlR3, const double VAlH4, const double VAlX5
  , const double VAlR5, const double VAlH5) ;

	TComplicated_Body  (const TComplicated_Body &R);

	TComplicated_Body operator=(TComplicated_Body  R);

	__fastcall calcCapacity();

	void calcCentreOfGravity(double *arrCentreGrav);

	double calcMass();

	void setupPointersArray(TSimpleBody_3D **SimpleBody_3D );

	void doCentreUp();

	void calcInertiaMtrx(double *arrMtrxInertia);



};


#endif
