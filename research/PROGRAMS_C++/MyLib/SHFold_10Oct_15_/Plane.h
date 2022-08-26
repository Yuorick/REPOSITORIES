//---------------------------------------------------------------------------

#ifndef PlaneH
#define PlaneH

class TURPointXY;
class TPlane
{
public:

 double  marrS0[3]; // ����� �� ��������� - ������ ��������� ������������� ������� ��������� �� ���������  (���)
 double  marrF[9]; //�������������  ������� �������� �� ������� ��������� ��������� � �������� ������� ��������� (����������)
									// ������� �� �������� ��������� ����� ������� ��������� ���������. ������ ������� -
									// ��� ��� X ���, ���������� ������� - ������ ��������� ������� � ���������
									// �������������� ��������� ����� �� ��� � �������� �� �������� ���:
									// S = arrF * S����� +  marrS0
	 TPlane();
 // ����������� �����������
 TPlane (const TPlane &R) ;
 // �������� ������������
 TPlane operator=(TPlane  R);
 // �������� ������
 TPlane( double* arrS0, double* arrF);

bool findIntersectingPoint_with_Line(double *arrPosWorking, double *arrVeloWorking, TURPointXY *ppntIntersect) ;

void fillNormalVect(double *arrN);

void transform_xyzSSK_to_xyzSKP(double *arrPointINtersect_PrSK, double *arrPointINtersect_SKP);

void transform_xyzSKP_to_xyzSSK(double *arrPoint_SKP, double *arrPoint_PrSK);


};
#endif
