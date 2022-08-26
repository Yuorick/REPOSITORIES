//---------------------------------------------------------------------------

#ifndef LearnShellBodyH
#define LearnShellBodyH
#include "URPolyLine.h"
#include "Detonator.h"
#include "Constants.h"



class TURPolyLine;
class TDetonator;



// ����� ��������� �������������� �������
class TLearnShellBody
{
public:
// ��� �������(������)
	enumShellType  mEnumShellType  ;
 // ����������� ����� �������, �
	 double mL;
// ������� ������, �
	 double mDm ;
// ���������� �� ��������� �� �� ��
	 double mLc;
 // ����� ��
	 double mLg ;
// ��������� h_��� ��� ������� ��������� �������  mh_gob =  Lc + 0.57*Lg -0.16*Dm
	 double mh_gob ;
// �����
	 double mMass ;
// ������ �������
	 double mvalIx0 ;

//   �������� ������������� ����� ������� � ����������� �� ����� ����
	 TURPolyLine mplnCx;

//   �������� ������������� ����� ������� � ����������� �� ����� ����  �� ���������
	 TURPolyLine mplnKnm;

//   �������� ������������� ����� ������� � ����������� �� ����� ����  �� ������� ��������
	 TURPolyLine mplnMxOmx;

	 //   �������� ������������� ��������� ���� ����������� �� ����� ����
	 TURPolyLine mplnCz;

	 // ��������� �������� ���������� �������� �������  ��������
	 double mDispOmega0;

	 // ��������� �������� ���������� �������� ������� ��������
	 double mDispV0;

	 // ��������� ��������������  �������� Cx (1 + delta)
	 double mDispCx;

	 // ��������� ��������������  �������� Cz (1 + delta)
	 double mDispCz;

	 // ���������  �������� �����
	 double mDispMass0;






	 double mV0;   // ��� ��������

	 double mOmega0 ;  // ��� ������� ��������
//



	// ����������� �� ���������
	TLearnShellBody () ;
	// ����������� �����������
	TLearnShellBody  (const TLearnShellBody  &R) ;

	// �������� ������������
	TLearnShellBody  operator=(TLearnShellBody   R2) ;
   // ���������������� �����������
	 TLearnShellBody( enumShellType EnumShellType) ;
	 // ���������������� �����������
	TLearnShellBody( enumShellType EnumShellType, const TDetonator Detonator);

	TLearnShellBody( enumShellType EnumShellType, const TURPolyLine PLN_CxEtal0
	 , const TURPolyLine PLN_Knm0 ,const TURPolyLine PLN_MxOmegax0);


	 void fnkCxEtal (const  double valM,  double &val_CxEtal, double &val_Grad_CxEtal) ;
	 void fnkMxOmegax (const  double valM,  double &val_MxOmegax, double &val_Grad_MxOmegax) ;
	 void fnkKnm (const  double valM,  double &val_Knm, double &val_Grad_Knm)  ;
	 double calcBallisticCoeff();
     double fnkCz (const  double valM);

	 void fillShellVozmDispMatr_ShootingEarth (const double VAlSigTechAngAU
, const double VAlSigPiAtm, double *arrMtrxShellDisp);



}  ;

#endif