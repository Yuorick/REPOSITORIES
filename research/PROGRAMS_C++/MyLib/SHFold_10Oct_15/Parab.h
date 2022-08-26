//---------------------------------------------------------------------------

#ifndef ParabH
#define ParabH
class TURPointXY ;
class TParab
{
public:

 double ma; // ����� �� ��� X
 double ml ;  //  ����� ����������



   TParab();
 // ����������� �����������
 TParab (const TParab &R) ;
 // �������� ������������
 TParab operator=(TParab  R);
 // �������� ������
 TParab ( const double a, const double b);

 void  clcCoeff ( double &vala, double &valb, double &valc) ;
 double  TParab::clcFncVal ( const double valx) ;

 int  fncLineIntersectParab (TURPointXY *pPontsInp,TURPointXY *pPontsOut) ;
 int fncSegmentCutParab( TURPointXY pntSegm0,  TURPointXY pntSegm1
							,TURPointXY *arrPntRez // ����� ����� � �����������
								) ;



};
#endif
