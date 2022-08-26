//---------------------------------------------------------------------------

#ifndef MeasureH
#define MeasureH
class TMeasure
{
public:
	// ����������� ������� ��������� �������
	int mDimY;


//  1.����� �������� ������
	double mTYZv;
//  2.������ ������
	double *mparrYZv;

 // 3.�������������� ������� ���
	double *mparrBMO_K;
	// 3.�������������� ������� ���
	double *mparrMMO_K;





   ~TMeasure() ;

	TMeasure () ;

	// ����������� �����������
	TMeasure  (const TMeasure  &R) ;
	// �������� ������������
	TMeasure  operator=(TMeasure   R2) ;

	TMeasure(const int DimY, const double ValDispBMO );

	TMeasure(const int DimY, const double ValDispBMO , const double ValDispMMO )  ;


};
#endif
