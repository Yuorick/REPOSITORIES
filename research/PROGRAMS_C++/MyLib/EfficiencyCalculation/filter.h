//File filter.h
#ifndef filterH
#define filterH
extern int workFilter;
//���� �������� ������ ��� ��������� ����������
typedef struct Filter_Input
{
 double tk;//����� 
 double Xkg,Ykg,Hkg;//��������� ������ ��������� ���� � ����-��
 double sbet,seps,sD;//��������� ��������� ���������
 double Vcx,Vcy,Vch;//������ �������� ������ �������
}Filter_Input;//��������� ������� ������ � �������� ����������

typedef struct Filter_Output
{
 double X, Y, H;//��������� ������ ��������� ���� � ����-��
 double Vx,Vy,Vh;//������ ������ �������� ����
}Filter_Output;//��������� �������� ������ �� ��������� ����������

//�������� �������, ����������� �������� ����������
void filtr_ukf(Filter_Input *TFilter_Input, Filter_Output *TFilter_Output, bool bSliga, const double VAlSigW, const double VAlSigMMO);
void  fncExtrapolateCorMtrx (const double  p11h, const double p12h, const double p22h
		,const double  t,const double  VAlDispW, double *sp11h, double *sp12h, double *sp22h );

void  calcCoefUsilenia (const double  sp11h, const double sp12h
		,const double  VAlSigBMO, const double  VAlSigMMO, double *pUs0, double *pUs1, double *pvalP, double *pvalAlf);
void calcCorMtrxCompleted(const double VAlKExtr00, const double VAlKExtr01, const double VAlKExtr11
	,const double VAlP,const double VAlAlf,const double VAlP1
	, double* pvalK00, double* pvalK01, double* pvalK11);
#endif


 //End filter.h
