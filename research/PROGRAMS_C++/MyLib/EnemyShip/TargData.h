//---------------------------------------------------------------------------

#ifndef TargDataH
#define TargDataH
class TTargData
{
public:
	 double mLDanger ;// ����� ������� ����� ����
	 double mLTarg  ; // ����� ����

	// ����������� �� ���������
	TTargData () ;
	// ����������� �����������
	TTargData  (const TTargData  &R) ;
	// �������� ������������
	TTargData  operator=(TTargData   R2) ;

	// ����� �����������
	 TTargData (const double LDanger, const double LTarg)  ;

}  ;
#endif
