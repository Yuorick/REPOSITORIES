//---------------------------------------------------------------------------

#ifndef ZamerH
#define ZamerH
// ��������� ������
class TZamer
{
public:

	 // ������ ���������
	double marrMeas[3];
	// �������������� ������� ������
	double marrCorr[9];
	// ����� ������
	double mT;





 __fastcall  TZamer() ;
// ����������� �����������
__fastcall  TZamer (const TZamer &R2) ;
 // ����� ������
 __fastcall TZamer ( double *arrMeas,double *arrCorr, const double VAlT) ;

 // �������� ������������
 TZamer   &operator=(const TZamer  &R2) ;




};
#endif
