//---------------------------------------------------------------------------

#ifndef TargBearing0H
#define TargBearing0H

// �������� ������ �� ��������� ������� ����
class TTargBearing0
{
public:
	//  ���� ������� ���� � ���
	double mBearing ;
	// ���� ����� ���� � ���
	double mTargCourse ;
	// ���� ����� ������� �������� ���� � ����������
	double mTargHorizAng ;
	//- ��������
	double mV;
	// ���������,
	double mR ;
	// ������
	double mH ;
	//
  //	double mMy ;
	//
	// ���
   //	double mTargEPR;
   //	enumTargetType mTargType;

	 __fastcall ~TTargBearing0() ;
	// ����������� �� ���������
	TTargBearing0 () ;
	// ����������� �����������
	TTargBearing0  (const TTargBearing0  &R) ;
	// �������� ������������
	TTargBearing0  operator=(TTargBearing0   R2) ;

	// ����� �����������1
	TTargBearing0 (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
	 ,const double  R,const double H);//, const double My , enumTargetType TargType );

  //	 TTargBearing0 (const double Bearing, const double TargCourse, const double TargZenitAng ,const double V
  //	 ,const double  R,const double H, const double My, const double TargEPR , enumTargetType TargType );

	 void raschet_nach_coord (double *pX) ;



}  ;
#endif
