//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <string.h>

#include "Helic.h"


__fastcall  THelic::THelic()
{

	 memset(marrRotor, 0, 2 * sizeof(TRotor));
	// ������� ���������� �����
	 memset(marrRotorOmega, 0, 2 * sizeof(double));
	// ����� ���������
	 mHelicMass = 0.;

}
// ����������� �����������
__fastcall  THelic::THelic (const THelic &R)
 {
	 memcpy(marrRotor, R.marrRotor, 2 * sizeof(TRotor));
	 memcpy(marrRotorOmega, R.marrRotorOmega, 2 * sizeof(double));
	 mHelicMass  = R.mHelicMass ;

 }
 // �������� ������������
  THelic THelic::operator=(THelic  R)
 {
	 memcpy(marrRotor, R.marrRotor, 2 * sizeof(TRotor));
	 memcpy(marrRotorOmega, R.marrRotorOmega, 2 * sizeof(double));
	 mHelicMass  = R.mHelicMass ;
	 return *this ;
 }

 // ����� ������ 1
 __fastcall THelic::THelic(TRotor *arrRotor, double *arrRotorOmega
   ,const double HelicMass)
 {
	 memcpy(marrRotor,arrRotor, 2 * sizeof(TRotor));
	 memcpy(marrRotorOmega, arrRotorOmega, 2 * sizeof(double));
	 mHelicMass  = HelicMass ;
 }
 //-------------------------------------------------------
 //---------------------------------------------------------------------------
// ���������� ��������� Ct �� ������ ���
// ������ ����� �� ������ ������� (������������), ����� ���������, ����������� � ���������� �����
// ��������� ������ ��������� Ct
	// TKel0 - ����������� �� ����������� �����, ���� ��������
	// HeliH -  ������  �������, �
	// BladeOmega  - ������� ���������� �����
double THelic::calcCt(const double TKel0,const double HeliH)
{
	  double valTKel = TKel0  - 0.0065 * 5500;
	  double valTCel = valTKel- 273.15;
	  double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
	  double valRo =   valPAm * 0.028964/8.31/  valTKel;

	  // ������ ��������� ������� �������� �����
	  double valR = marrRotor[0].mBlade.mBladeR;
	  ///

	  double valCt = 2.* 9.81 * mHelicMass/
	  (valRo   *M_PI * valR* valR *valR *valR
			*marrRotorOmega[0] * marrRotorOmega[0]);
	  return valCt;
}
//--------------------------------------------------------------
//-------------------------------------------------------------
	// ���������� ������������ ��������� ���� �������  �� ������ ���
	// ������ ����� ������� ��� �������� ������ - ����� ���������
	// , ����������� ���������� ������� � ������ �������
	// �� ��������� ���� ������ ����� ��������� ��������� �������, ��������� ����
	// � ����� ����������� Cyalfa
	/// TKel0 - ����������� �� ����������� �����, ���� ��������
	// HeliH -  ������  �������, �
	// BladeOmega  - ������� ���������� �����
	// ������������ ����� ��� ����� 15 ���� !!!!
double THelic::calc_C_y_alfa(const double TKel0,const double HeliH)
{
	double valTKel = TKel0  - 0.0065 * 5500;
	double valTCel = valTKel- 273.15;
	double valPAm = 101325. * exp(log(valTKel/ TKel0) * 9.81* 0.028964/ 8.31/ 0.0065);
	double valRo =   valPAm * 0.028964/8.31/  valTKel;
	double valFiMax = 15./180. * M_PI; // ������������ ��� 15 ����   !!!!!
	double valBladeS =  marrRotor[0].mBlade.mBlade_b * (marrRotor[0].mBlade.mBladeR - marrRotor[0].mBlade.mRadHorizHsarnir);
	double valC_y_alfa = 2. * mHelicMass * 9.81
	/ (((double)(marrRotor[0].mQuantBlades)) * valFiMax * valRo
	* valBladeS *marrRotor[0].mBlade.mBladeR * marrRotor[0].mBlade.mBladeR
			* marrRotorOmega[0] * marrRotorOmega[0] );
   return valC_y_alfa ;
}
#pragma package(smart_init)
