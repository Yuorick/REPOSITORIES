//---------------------------------------------------------------------------

#ifndef TransmitAntH
#define TransmitAntH
// ���������� �������
class TTransmitAnt
{
public:

 // �������� �� ��������
	double mPowerPrd;
	// �� �� ��������
	double mKYPrd;




 __fastcall  TTransmitAnt() ;
// ����������� �����������
__fastcall  TTransmitAnt (const TTransmitAnt &R2) ;

 // �������� ������������
 TTransmitAnt   &operator=(const TTransmitAnt  &R2) ;

  // ����� ������
 __fastcall TTransmitAnt(const double PowerPrd,const double KYPrd);





};
#endif
