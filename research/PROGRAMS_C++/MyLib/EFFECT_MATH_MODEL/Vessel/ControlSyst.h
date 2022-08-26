//---------------------------------------------------------------------------

#ifndef ControlSystH
#define ControlSystH

#include "Table_1D.h"

class TTable_1D;
class TControlSyst
{
public:
	//���� ����������
	double mFiltT;
	// �������� ����
	double mSinsDelayT;

	// ���� ������� ������ ���
	double mRzvT;

     // ������������� �������
	 TTable_1D mTblPereletSimulated;



	// ����������� �� ���������
	TControlSyst () ;
	// ����������� �����������
	TControlSyst  (const TControlSyst  &R) ;
	// �������� ������������
	TControlSyst  &operator=(const TControlSyst   &R2) ;

 //	TControlSyst (const double FiltT, const double SinsDelayT)  ;

	TControlSyst (const double FiltT, const double SinsDelayT, const double RzvT)  ;

    TControlSyst (const double FiltT, const double SinsDelayT
  , const double RzvT, const TTable_1D TblPereletSimulated);




}  ;
#endif
