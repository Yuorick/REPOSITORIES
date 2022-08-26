//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "ConstantsCalculation.h"

//---------------------------------------------------------------------------
// ���������� ��������� Ct
// ������ ����� �� ������ ������� (������������), ����� ���������, ����������� � ���������� �����
// ��������� ������ ��������� Ct
	// TKel0 - ����������� �� ����������� �����, ���� ��������
	// HeliH -  ������  �������, �
	// BladeOmega  - ������� ���������� �����
	// BladeR - ������ �����
	// HelicMas - ����� ���������
double calcCt(const double TKel0,const double HeliH, const double BladeOmega
	, const double BladeR,const	double HelicMass)
{
	  double valTKel = TKel0  - 0.0065 * 5500;
	  double valTCel = valTKel- 273.15;
	  double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
	  double valRo =   valPAm * 0.028964/8.31/  valTKel;
	  double valCt = 2.* 9.81 * HelicMass/  (valRo   *M_PI * BladeR * BladeR *BladeR *BladeR
			*BladeOmega *BladeOmega);
	  return valCt;
}
//--------------------------------------------------------------
// ���������� ���������� �� ����� ���������� ���������������� ����
// ������������� �� ������ ������  ��� ��
// VAlR - ������ ��������� �������
// VAlRadHorizHsarnir - ������ ���������������� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
double calc_rQ(const double VAlR, const double VAlRadHorizHsarnir
	, const double VAlPofile_d0,const double VAlPofile_d1 )
{
   double x0 =  VAlRadHorizHsarnir;
   double x1 = VAlR;
   double valK = (VAlPofile_d1  - VAlPofile_d0)/(x1 - x0);
   double temp =  VAlPofile_d0 - valK * x0;
   double temp2 = (x1 * x1 - x0 * x0) / 2.;
   double temp3 = (x1 * x1* x1 - x0 * x0* x0) / 3.;
   double temp4 = (x1 * x1* x1 * x1 - x0 * x0 * x0 * x0) / 4.;
   double val_rQ = (temp * temp3 + valK * temp4 )/ ( temp * temp2 + valK * temp3);
   return val_rQ;
}
//-------------------------------------------------------------
// ���������� ���������� �� ������ ���� �������  �� ������ ����� ������� � ������ ��������
// VAl_L - ����� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
double calc_X_BladeCentreMass(const double VAl_L, const double VAlPofile_d0,const double VAlPofile_d1 )
{
	return VAl_L * (VAlPofile_d0 + 2. * VAlPofile_d1)/ (VAlPofile_d0 +  VAlPofile_d1) /3.;
}
//-------------------------------------------------------------
// ���������� ������������ ������� ������������� ������ ��������
// VAlR - ������ ��������� �������
// VAlRadHorizHsarnir - ������ ���������������� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
// VAlM - ����� �������
double calc_X_StatMoment_Sg(const double VAlR, const double VAlRadHorizHsarnir
	, const double VAlPofile_d0,const double VAlPofile_d1 , const double VAlM )
{
	double valCentreMass = calc_X_BladeCentreMass(VAlR - VAlRadHorizHsarnir,  VAlPofile_d0, VAlPofile_d1 );
	return   VAlM * ( valCentreMass + VAlRadHorizHsarnir);
}

//-------------------------------------------------------------
// ���������� ������� �������  ������������� ������ �������
// VAl_L - ����� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
// VAlM - ����� �������
double calcInertiaMoment0(const double VAlL, const double VAlPofile_d0,const double VAlPofile_d1,
  const double VAlM )
{

	return   VAlM * VAlL * VAlL *( VAlPofile_d0 * VAlPofile_d0 + VAlPofile_d0 * VAlPofile_d1 * 4.
	 + VAlPofile_d1* VAlPofile_d1)
	   / (VAlPofile_d0 + VAlPofile_d1) / (VAlPofile_d0 + VAlPofile_d1)/18.;
}
//-------------------------------------------------------------
// ���������� ������� �������  ������������� ������ ������ ��������������� �������
// VAl_L - ����� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
// VAlM - ����� �������
double calcInertiaMomentHorSharnir(const double VAlL, const double VAlPofile_d0,const double VAlPofile_d1,
  const double VAlM )
{

	return   VAlM * VAlL * VAlL *( VAlPofile_d0  +  3. *VAlPofile_d1)
	   / (VAlPofile_d0 + VAlPofile_d1) /6.;
}
//-------------------------------------------------------------
// ���������� ������� �������  ������������� ������ �������� �������
// VAlR - ������ ��������� �������
// VAlRadHorizHsarnir - ������ ���������������� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
// VAlM - ����� �������
double calcInertiaMoment(const double VAlR, const double VAlRadHorizHsarnir
	, const double VAlPofile_d0,const double VAlPofile_d1 , const double VAlM  )
{
	double valXCentre = calc_X_BladeCentreMass(VAlR - VAlRadHorizHsarnir,  VAlPofile_d0, VAlPofile_d1 ) ;
	return   calcInertiaMoment0(VAlR - VAlRadHorizHsarnir
	  ,  VAlPofile_d0, VAlPofile_d1,VAlM ) + (valXCentre  +  VAlRadHorizHsarnir)* (valXCentre  +  VAlRadHorizHsarnir) * VAlM;
}

//-------------------------------------------------------------
// ���������� ������������ ��������� ���� �������
// VAlR - ������ ��������� �������
// VAlRadHorizHsarnir - ������ ���������������� �������
// VAlPofile_d0 - ������ ������� ������������� ������� �������  � ����� ��������� � ������ ��������
// VAlPofile_d1 - ������ ������� ������������� ������� �������  � ������� ����� �� ������ ��������
// VAlM - ����� �������
double calc_C_y_alfa(const double TKel0,const double HeliH, const double BladeOmega
	, const double BladeR , const double VAlRadHorizHsarnir, const double VAlBlade_b
	, const int QUantBlades, const	double HelicMass)
{
	double valTKel = TKel0  - 0.0065 * 5500;
	double valTCel = valTKel- 273.15;
	double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
	double valRo =   valPAm * 0.028964/8.31/  valTKel;
	double valFiMax = 15./180. * M_PI; // ������������ ��� 15 ����
	double valBladeS =  VAlBlade_b * (BladeR - VAlRadHorizHsarnir);
	double valC_y_alfa = 2. * HelicMass * 9.81
	/ (((double)QUantBlades) * valFiMax * valRo * valBladeS *BladeR *BladeR
			*BladeOmega *BladeOmega );
   return valC_y_alfa ;
}

#pragma package(smart_init)