//---------------------------------------------------------------------------

#ifndef GlonassH
#define GlonassH
class TGlonass
{
public:
 //  C�� �� ���������
 long double mSKOPos;
 // C�� �� C�������
 long double mSKOVeloc;




   TGlonass();
 // ����������� �����������
 TGlonass (const TGlonass &R) ;
 // �������� ������������
 TGlonass operator=(TGlonass  R);
 // �������� ������
 TGlonass(long const double SKOPos,long const double SKOVeloc );

 };
#endif
