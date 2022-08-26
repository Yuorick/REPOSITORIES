//---------------------------------------------------------------------------

#ifndef ArcParabH
#define ArcParabH

// �������� ��� ������� ����� ����� ������������� �� ������ ������ (���������� ����������� ��������)
//  � ������ ����� (���������� ������� ��������).
// ������ � �������� � ����������, �������� �������� ���������� ��������.
// ��� ����� ���� ���������� ��� ���������� ������� � ��������� ����������������.
///////////////////////////////////////////////////////////////////////////////////
// �������� ����������� ���� �������
// ml - ��������� ��������
// ml/2 - �������� ����������
// ma - ���������� ������� �������� �� ��� OX
// �������� ���� X*X = -2PY
// ��������� ��������
// Y=-(X-ma)*(X-ma)/(2*ml)+ ml/2
// �� ���� ����� ���������� �� ��� OX
// ���������� ������� (ma, ml/2)
// ����� ���������� ����
class TURPointXY ;
class TArcParab
{
public:

 double ma; // ����� �� ��� X
 double ml ;  //  ����� ����������
 double m_x0 ;
 double m_x1 ;


   TArcParab();
 // ����������� �����������
 TArcParab (const TArcParab &R) ;
 // �������� ������������
 TArcParab &operator=(const TArcParab  &R);
 // �������� ������
 TArcParab ( const double a, const double b);
 	// ����� ������ 2
TArcParab( const double a, const double l, const int isign) ;


 void  clcCoeff ( double &vala, double &valb, double &valc) ;
 double  calcFncValue ( const double valx) ;

 int  fncLineIntersectParab (TURPointXY pPontsInp0,TURPointXY pPontsInp1,TURPointXY *pPontsOut);

 int fncSegmentCutArcParab( TURPointXY pntSegm0, TURPointXY pntSegm1 // �������
							,TURPointXY *arrPntRez // ����� ����� � �����������
								);
  void ShowMe(wchar_t *FileName);
  int fncSegmentCutParabLine(TURPointXY pntSegm0, TURPointXY pntSegm1 // �������
							,TURPointXY *arrPntRez // ����� ����� � �����������
								) ;
  void ShowFullGraph(wchar_t *FileName);



};

 bool Is_X_BelongeSegm( const double valx, const double vala,const double valb);
 double max_( const double x0, const double x1);
 double min_( const double x0, const double x1);


#endif
