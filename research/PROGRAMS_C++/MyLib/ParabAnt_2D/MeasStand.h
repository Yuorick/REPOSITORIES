//---------------------------------------------------------------------------

#ifndef MeasStandH
#define MeasStandH
#include "MeasStand.h"
#include "Comp.h"
class TComp;
class TMeasStand
{
// ������ ��������� ������ �������� �� ������ �������
// mRadAxeAnt - ���� ������� ����������� ��� ������� ������������ ���������
// ����������� ��� �������� ����� �������� ����, ������������� ����� 3 � 4 ��������
// mpcmparrMeas[3] - ������ ��������
public:
 // ���� ������� ��� �������
 double	 mRadAxeAnt ;
 // ����� �������� � ���������
 int mCountDgrs;
 // ������ ��������� �������� � �����
TComp  *mpcmparrMeas ;



 __fastcall~TMeasStand();
 TMeasStand() ;


// ����������� �����������
 TMeasStand(const TMeasStand &R) ;
 TMeasStand operator=(TMeasStand  R2) ;
 // ����� ������
__fastcall TMeasStand(const double RadAxeAnt, TComp *pcmparrMeas);
__fastcall TMeasStand(const double RadAxeAnt,  int CountDgrs);



};
#endif
