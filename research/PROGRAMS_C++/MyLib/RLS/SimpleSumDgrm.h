//---------------------------------------------------------------------------

#ifndef SimpleSumDgrmH
#define SimpleSumDgrmH
class TComp;
class TSimpleSumDgrm
{
public:
 // ���������� ����� �������� ���������
 double mKSum;
 //���� ������������ ��������� ���������
 double mScnSum;
 // ������� ���������� ��������� ���������
 double mTension;


 __fastcall  TSimpleSumDgrm::TSimpleSumDgrm() ;
// ����������� �����������
__fastcall  TSimpleSumDgrm (const TSimpleSumDgrm &R2) ;
 // ����� ������
 __fastcall TSimpleSumDgrm(const double KSum
	   , const double ScnSum, const double Tension ) ;
// �������� ������������
TSimpleSumDgrm   operator=(TSimpleSumDgrm  R2) ;

void  _fastcall calcPartial_vectG_and_mtrxH(TComp cmpSZv, double alfTrg,double  alfAnp,double valK11, double valK12
	 ,double valK21,double valK22,double* arr_gradK11, double* arr_gradK12, double* arr_gradK21
	 ,double* arr_gradK22, double*arr_HessK11, double* arr_HessK12, double* arr_HessK21
	 ,double* arr_HessK22,double*   arr_FGreek, double*  arr_dFGreek );
};
#endif
