//---------------------------------------------------------------------------

#ifndef EbombH
#define EbombH
class TEbomb
{
public:
 //  ���� ������ ���������
 long double mFi;
 // ��������� ������ ���������
 long double mR ;
 // ����� ��������

 long double mW;



   TEbomb();
 // ����������� �����������
 TEbomb (const TEbomb &R) ;
 // �������� ������������
 TEbomb operator=(TEbomb  R);
 // �������� ������
  TEbomb(const long  double Fi,const long  double R,const long  double W );
 long double TEbomb::fncCalcPower ();

static void createGraphR_of_Fi (wchar_t *pwchOutFile, double valW, double valP);

 };
#endif
