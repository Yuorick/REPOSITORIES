//---------------------------------------------------------------------------

#ifndef Table_2DH
#define Table_2DH

// ����� ��������� 2-� ������ �������
// ������� ������������ �� ���� ������ (������) ���������� ������ � ������
// �������� ��������� mparrArg
// �������� ������� �������� � �������  mparrVal
// �������� ���������� � ������� mparrArg
// ������ mparrArg ������������ �� ���� ���������� ������������������ �����������
// mNumCols - ����� ����� � �������� ������ �������
class TTable_1D;
//
class TTable_2D
{
public:

	int mNumCols ;//
	double *mparrArg; // ������ �������� ��������
	TTable_1D *mpArrTable_1D; // ������ ���������� ������-��������
	double Box[2];


	__fastcall ~TTable_2D() ;

	TTable_2D();
 // ����������� �����������
 TTable_2D (const TTable_2D &R) ;
 // �������� ������������
 TTable_2D &operator=(const TTable_2D  &R);
// ����� ������
TTable_2D( double *parrArg, TTable_1D *pArrTable_1D, const int NumCols);
TTable_2D( double *parrArgTab1, const int NumColsTab1, double *parrArgTab2, const int NumColsTab2, double *parrVal );

double calcValue(const double VAlRowArg, const double VAlColArg);

void calcBoundBox();

void ShowMe(wchar_t *FoldName);

};
#endif
