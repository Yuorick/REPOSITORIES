//---------------------------------------------------------------------------

#ifndef RectRectWaveGuideH
#define RectRectWaveGuideH
// �������� ����� ��������� �������������� �������
class TRectWaveGuide
{
public:


 double ma ;  // ������� ������� ��������������
 double mb ; // ������� ������� ��������������
 double mFi1 ; //  �����������


   TRectWaveGuide();
 // ����������� �����������
 TRectWaveGuide (const TRectWaveGuide &R) ;
 // �������� ������������
 TRectWaveGuide operator=(TRectWaveGuide  R);
 // �������� ������
  TRectWaveGuide(  const double a,  const double b );

 // void ShowMe(wchar_t *FileName);
 double fncDiagr(const double VAlLambda, const double VAlTetta);

};
#endif
