//---------------------------------------------------------------------------

#ifndef ArtCannonH
#define ArtCannonH

enum enumCannonType {A_192M,A_190_01, AK_630,AK_630_2M,AK_176,CANNON_UNKNOWN};
class TArtCannon
{
public:
 // ���
  enumCannonType menumCannonType ;
// ����������������. �
	double mRateOfFire;
// �������� �� ����, ���, �����
	double mAngGroupedFire;
// ��� ����������������� �������
	double mSigRabT;
	// �������� �����
	double mRabT;
	// ����������� �� ���������
	TArtCannon () ;
	// ����������� �����������
	TArtCannon  (const TArtCannon  &R) ;
	// �������� ������������
	TArtCannon  &operator=(const TArtCannon   &R2) ;

	// ����� �����������2
   TArtCannon (const double RateOfFire, const double AngGroupedFire
  ,const double SigmaDelayT, const double RabT);

   TArtCannon (enumCannonType EnumCannonType, const double SigmaDelayT
   ,const double RabT);


}  ;
#endif