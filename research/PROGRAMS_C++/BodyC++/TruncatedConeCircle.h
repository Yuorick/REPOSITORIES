//---------------------------------------------------------------------------

#ifndef TruncatedConeCircleH
#define TruncatedConeCircleH
#include "SimpleBody_3D.h"
// ����� ���������  ���������   �������� ����� � ���������� R � r. ���������� ����� ����������� H
// ��������� � �������� R ����������� � ������������ ��������� OXZ
// ������ ��������� ����������� � ������ ���������.
// ��������� � �������� r ����������� � ������������ ��������� ��������� �� ��� X �� H
// �������� ����������� ��������� ������ �����:
// OX - � �����
// OY - �����
// OZ - �� ���
// ����� ������ ����� 1
class TTruncatedConeCircle : public TSimpleBody_3D
{
public:

	double mR;
	double mr;
	double mH;
   //	double mM; // �����



	//__fastcall ~TTruncatedConeCircle() ;


	 TTruncatedConeCircle () ;

	 // ����������� �����������
	  TTruncatedConeCircle  (const TTruncatedConeCircle  &R) ;

	  // �������� ������������
	  TTruncatedConeCircle  operator=(TTruncatedConeCircle   R2) ;

	  TTruncatedConeCircle(const TPlane Plane, const double R, const double r,const double H,const double M);

	  virtual double calcCapacity() ;

	  virtual void calcCentreOfGravity(double *arrCentreGrav)  ;

	  virtual  void calcInertiaMtrx(double *arrInertMtrx) ;

	  double calcJyy_For_ConeCircle(const double VAlR, const double VAlH, const double VAlM) ;

	  double calcJxx_For_ConeCircle(const double VAlR, const double VAlM)  ;


	  double calcCapacity_For_ConeCircle(const double VAlR, const double VAlH);

	  void calcCentreOfGravity_For_ConeCircle(const double VAlH, double *arrCentreGrav);

	  void calcInertiaMtrx_For_ConeCircle(const double VAlR, const double VAlH, const double VAlM
	 ,double *arrInertMtrx);



};

	  #endif