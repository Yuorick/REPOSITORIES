//---------------------------------------------------------------------------

#ifndef SingleSignH
#define SingleSignH
// ��������� ������
class TSingleSign
{
public:

   // ���� �����
	double mBet;
	// ���� �����
	double mEps;
	// ���������
	double mAmpl;
	// ����
	double mPhase;




 __fastcall  TSingleSign() ;
// ����������� �����������
__fastcall  TSingleSign (const TSingleSign &R2) ;
 // ����� ������
 __fastcall TSingleSign(const double Bet,const double Eps, const double Amp, const double Phase);

 // �������� ������������
 TSingleSign   &operator=(const TSingleSign  &R2) ;




};
#endif