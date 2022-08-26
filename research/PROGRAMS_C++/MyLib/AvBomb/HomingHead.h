//---------------------------------------------------------------------------

#ifndef HomingHeadH
#define HomingHeadH
// ����� ���������� ��������������� �������� �� �������� ������� ����������
class THomingHead
{
public:
 //  ���� ���������
 long double mFi;
 // ��������� ���������
 long double mR ;
 // ����� �� ���������, �
 long double mSigFl_R ;
 // ����� �� ����, ����
 long double mSigFl_U ;
 // �������� �������� ����� ��������
 long double m_h ;

 THomingHead();
 // ����������� �����������
 THomingHead (const THomingHead &R) ;
 // �������� ������������
 THomingHead operator=(THomingHead  R);
 // �������� ������
 THomingHead( const long  double Fi, const long  double R
   ,  const long  double SigFl_R, long  double SigFl_U , long  double h   );
 };
#endif
