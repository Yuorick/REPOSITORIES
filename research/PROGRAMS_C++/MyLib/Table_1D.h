//---------------------------------------------------------------------------

#ifndef Table_1DH
#define Table_1DH
// ����� ��������� ��������� �������
// ������� ������������ �� ���� ������ (�������) �������� ��������� ������� � ������
// �������� ������� �������� � �������  mparrVal
// �������� ���������� � ������� mparrArg
// mNumCols - ����� ������ � �������� ������ �������

//
class TURPolyLine;
class TTable_1D
{
public:

	int mNumCols ;// �-�� ��������
	double *mparrArg; // ������ �������� ��������
	double *mparrVal; // ������ ���������� ��������



	__fastcall ~TTable_1D() ;

	TTable_1D();
 // ����������� �����������
 TTable_1D (const TTable_1D &R) ;
 // �������� ������������
 TTable_1D &operator=(const TTable_1D  &R);
// ����� ������ 1
TTable_1D(const  double *parrArg,const  double *parrVal, const int NumCols);
// ����� ������ 2
TTable_1D(const  double *parrPoints,  const int NumCols);

double calcValue(const double VAlArg);

void ShowMe(wchar_t *FileName);

static TURPolyLine createPln(const TTable_1D table) ;

//void calcBoundBox();

};
#endif
