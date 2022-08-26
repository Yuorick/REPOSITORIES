//---------------------------------------------------------------------------

#ifndef PartDataH
#define PartDataH
//---------------------------------------------------------------------------
class TPartData
{
public:
// ��� ��������. 0 - ����������� �������������, 1 - �����,
// 2 - ������������, 3 - ��������������� ��������
// 4- ������������� ������� ���������������� ���������
	int miTypePart;
//   ������������ �������
	double mTimePart ;
// ��� ���������� ��������
	double mSigW ;
// ������ ����������, ����������� ��������
// ��� ������������ ��������  miTypePart = 0
  //  marrData[i]  =0 i = 0,...,9

// ��� �����   miTypePart =1
// marrData[0]  - ���� ��������� ������� ����� ��������� � ��������� � ����������
// marrData[1],marrData[2]- �����, ���������� 1-�� �������
// marrData[3],marrData[4]- �����, ���������� 2-��-�� �������
// marrData[5],marrData[6]- �����, ���������� 3-�� �������

// ������ ����������, ����������� �������� ��� ������� "������������"  miTypePart = 2
// marrData[0]  - ���� ��������� ������� ����� ��������� � ��������� � ����������
// marrData[1],marrData[2]- ����, ���������
// marrData[3],marrData[4]- ����
// marrData[5],marrData[6]- ����

// ������ ����������, ����������� ��������������� ��������  miTypePart = 3
// marrData[0], marrData[1],marrData[2]  - ��������� �� ������ ��� ���
// ��������� ����� ����

// ������ ����������, ����������� �������������  �� ������ ��������� ��������� miTypePart = 4
//  marrData[0] - ����������� ���������
//  marrData[1] - ���������� ���� ����������
//  marrData[2] - ����������� ��������� ���������� ����������, � �/�2

// ��������� ����� ����
	double marrData[10];



	// ����������� �� ���������
	TPartData () ;
	// ����������� �����������
	TPartData  (const TPartData  &R) ;
	// �������� ������������
	TPartData  operator=(TPartData   R2) ;

	// ����� �����������
	TPartData  (const int  iTypePart  ,const double TimePart,const double SigW, double *arrData 	) ;





}  ;
#endif