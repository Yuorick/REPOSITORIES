//---------------------------------------------------------------------------

#ifndef LineMoveH
#define LineMoveH
class TLineMove
{
public:
// ���� �������� ����������
	double mT0;
//  ������ ��������� � ������ mT
	double marrVS[9];


	 __fastcall ~TLineMove() ;
	// ����������� �� ���������
	TLineMove () ;
	// ����������� �����������
	TLineMove  (const TLineMove  &R) ;
	// �������� ������������
	TLineMove  operator=(TLineMove   R2) ;

	// ����� �����������1
	TLineMove (const double T,  double *arrVS) ;

	void ExtrapolateVS (const double tExtr, double *arrVSTargExtr_GSK);



}  ;
#endif
