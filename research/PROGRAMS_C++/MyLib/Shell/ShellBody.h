//---------------------------------------------------------------------------

#ifndef ShellBodyH
#define ShellBodyH
#include "URPolyLine.h"
#include "Detonator.h"
#include "Constants.h"




class TURPolyLine;
class TDetonator;



// ����� ��������� �������������� �������
class TShellBody
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
 // ����� ��
	 double mExplMass;
// ���������
TDetonator mDetonator;
//   �������� ������������� ������������ ����� �������� ����������� �� �� ��� �������� �� ��� X
	 TURPolyLine mplnCix;
//   �������� ������������� ������������ ����� �������� ����������� �� �� ��� �������� �� ��� Y
	 TURPolyLine mplnCiy;
//   �������� ������������� ������������ ����� �������� ����������� �� �� ��� �������� �� ��� Z
	 TURPolyLine mplnCiz;


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
	TShellBody () ;
	// ����������� �����������
	TShellBody  (const TShellBody  &R) ;

	// �������� ������������
	TShellBody  &operator=(const TShellBody   &R2) ;
   // ���������������� �����������
	 TShellBody( enumShellType EnumShellType) ;
	 // ���������������� �����������
	TShellBody( enumShellType EnumShellType, const TDetonator Detonator);


	  void fnkCxEtal (const  double valM,  double &val_CxEtal, double &val_Grad_CxEtal) ;
	  void fnkMxOmegax (const  double valM,  double &val_MxOmegax, double &val_Grad_MxOmegax) ;
	  void fnkKnm (const  double valM,  double &val_Knm, double &val_Grad_Knm)  ;
	 double calcBallisticCoeff();
	 void fnkIx0 (const  double valTet0,  double &val_ix, double &val_Grad_ix);
	 void fnkIz0 (const  double valTet0,  double &val_iz, double &val_Grad_iz) ;

	 void fillShellVozmDispMatr_ShootingEarth (const double VAlSigTechAngAU
, const double VAlSigPiAtm, double *arrMtrxShellDisp);
	 double fnkCz (const  double valM);
}  ;


#endif