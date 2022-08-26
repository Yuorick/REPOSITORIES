//---------------------------------------------------------------------------

#ifndef Table_3DH
#define Table_3DH
// ����� ��������� 2-� ������ �������
// ������� ������������ �� ���� ������ (������) 2-� ������ ������ � ������
// �������� ��������� mparrArg
// �������� ������� �������� � �������  mparrVal
// �������� ���������� � ������� mparrArg
// ������ mparrArg ������������ �� ���� ���������� ������������������ �����������
// mNumCols - ����� ����� � �������� ������ �������
class TTable_2D;
//
class TTable_3D
{
public:

	int mNumCols ;//
	double *mparrArg; // ������ �������� ��������
	TTable_2D *mpArrTable_2D; // ������ ���������� ������-��������
	double Box[2];


    ~TTable_3D() ;

	TTable_3D();
 // ����������� �����������
 TTable_3D (const TTable_3D &R) ;
 // �������� ������������
 TTable_3D &operator=(const TTable_3D  &R);
// ����� ������
TTable_3D( double *parrArg, TTable_2D *pArrTable_2D, const int NumCols);
// ����� ������
TTable_3D( double *parrArg
	,TTable_2D Table_2_0
	,TTable_2D Table_2_1
	,TTable_2D Table_2_2
	) ;
	///
TTable_3D( double *parrArgTab1, const int NumColsTab1, double *parrArgTab2, const int NumColsTab2
 ,double *parrArgTab3, const int NumColsTab3, double *parrVal );

double calcValue(const double VAlTab3Arg,const double VAlRowArg, const double VAlColArg);

void calcBoundBox();

};
#endif