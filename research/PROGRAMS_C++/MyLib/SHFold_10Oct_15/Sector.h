//---------------------------------------------------------------------------

#ifndef SectorH
#define SectorH
#include "UrPointXY.h"
class TURPointXY ;
class TSector
{
public:

 TURPointXY mPntCentre; // �����
 double mR ;  //  ������
 double mFi0 ; // ������� ����  mFi0 > mFi1 , ����������� �� ������� �������
 double mFi1 ; //  �����������


   TSector();
 // ����������� �����������
 TSector (const TSector &R) ;
 // �������� ������������
 TSector &operator=(const TSector  &R);
 // �������� ������
  TSector(const TURPointXY PntCentre, const double R,  const double Fi0,  const double Fi1 );
 	// ����� ������ 2

  void ShowMe(wchar_t *FileName);


};

#endif
