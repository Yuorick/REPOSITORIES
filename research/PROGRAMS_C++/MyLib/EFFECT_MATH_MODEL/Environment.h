//---------------------------------------------------------------------------

#ifndef EnvironmentH
#define EnvironmentH
//---------------------------------------------------------------------------
// 9 ����� �������� ������������� 20 �/� �������� �����
class TEnvironment
{
public:
// �������� �� ���� �� ����� Significance Wave Height (SWH). �� 0 �� 9 ������
	int  mBallWave;
//  ���� ����� �/�
	 double mWind_V  ;
	 // ����������� ������ ���� ������
	 double mWind_Alf ;
	 //
	 double mWind_VertV;


	// ����������� �� ���������
	 __fastcall TEnvironment () ;
	// ����������� �����������
	 __fastcall TEnvironment  (const TEnvironment  &R) ;
	// �������� ������������
      TEnvironment  &operator=(const TEnvironment   &R2) ;

	 // ����� �����������
	 __fastcall  TEnvironment (const int BallWave, const double Wind_V);

	  __fastcall  TEnvironment (const double Wind_V, const double Wind_Alf, const double  Wind_VertV)
;
}  ;
#endif
