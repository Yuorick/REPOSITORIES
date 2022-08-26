//---------------------------------------------------------------------------

#ifndef EtalonSignH
#define EtalonSignH
class TEtalonSign
{
public:
 // ��������
 double mEtalonAmp;
 //���������
 double mEtalonDist;
 //���
 double mEtalonEPR;
 //��� ����� ���� ��������� ��������� 5�10
 double mNoiseSKZ_5P10;
 // ��� �������� ������� �������� ��������� ��������� 5�10
 double mEtalonSigAmplFact_5P10;
 // �������� �� ��������
 double mEtalonPowerPrd;
 // �� �� ��������
 double mEtalonKYPrd;
 // �������� �� �����
 double mEtalonKYPriem;


 __fastcall  TEtalonSign() ;
// ����������� �����������
__fastcall  TEtalonSign (const TEtalonSign &R2) ;
 // ����� ������
  __fastcall TEtalonSign(const double EtalonAmp,const double EtalonDist, const double EtalonAPR
   , const double NoiseSKZ_5P10, const double EtalonSigAmplFact_5P10) ;

  // ����� ������ 2
 __fastcall TEtalonSign(const double EtalonAmp,const double EtalonDist, const double EtalonAPR
   , const double NoiseSKZ_5P10, const double EtalonSigAmplFact_5P10, const double EtalonPowerPrd,  const double  mEtalonKYPrd
   , const double);


 // �������� ������������
 TEtalonSign   &operator=(const TEtalonSign  &R2) ;

 double calcNWaveEtalon();

};
#endif
