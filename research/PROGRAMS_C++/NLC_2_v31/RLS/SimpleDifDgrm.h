//---------------------------------------------------------------------------

#ifndef SimpleDifDgrmH
#define SimpleDifDgrmH
class TComp;
class TSimpleDifDgrm
{
public:

// ������� ���� ������� ���������� ���������
 double mTang;
 // ���� ������������ ���������� ���������
 double	 mScnDif;


 __fastcall  TSimpleDifDgrm ::TSimpleDifDgrm() ;
// ����������� �����������
__fastcall  TSimpleDifDgrm (const TSimpleDifDgrm &R2) ;
 // ����� ������
 __fastcall TSimpleDifDgrm(const double Tang, const double ScnDif) ;
// �������� ������������
TSimpleDifDgrm   operator=(TSimpleDifDgrm  R2) ;




};
#endif