//---------------------------------------------------------------------------

#ifndef BladeH
#define BladeH
// ����� ��������� ������� ����� ���������
class TBlade
{
public:

	// ������ ��������� �������
	double mBladeR ;
	// ���������� �� ��� �������� ��  ������ ��
	double mRadHorizHsarnir;
	// ������ ������������� �������  ������� ����� � ������ ��������
	double mPofile_d0;
	// ������ ������������� �������  ������� ������ �� ������ ��������
	double mPofile_d1;
	// ����� �������
	double mBladeM;

	// ����� �������
	double mBlade_b;



	__fastcall  TBlade() ;
	// ����������� �����������
	__fastcall  TBlade (const TBlade &R2) ;

	// �������� ������������
	TBlade   operator=(TBlade  R2) ;

	// ����� ������
	 __fastcall TBlade(const double  BladeR,const double  RadHorizHsarnir
   ,const double  Pofile_d0,const double  Pofile_d1,const double  BladeM  ,const double  Blade_b);

   double calc_rQ();

   double calc_X_BladeCentreMass()  ;

   double calc_X_StatMoment_Sg() ;

   double calcInertiaMoment0();

   double calcInertiaMomentHorSharnir() ;

   double  calcInertiaMoment() ;
};
#endif