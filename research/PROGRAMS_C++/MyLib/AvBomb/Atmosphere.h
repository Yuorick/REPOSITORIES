//---------------------------------------------------------------------------

#ifndef AtmosphereH
#define AtmosphereH
#define ATM_PN0  1000.  //���������� �������� �������� �������
#define ATM_RoN0  1.2054  // ���������� ������� ��������� �������
#define ATM_AN0  340.8   //���������� ������� �������� ����� � �������
#define ATM_TAYN0  289   //���������� ����������� ����������� �� 0 ������
#define ATM_R_UNIVER  287.05287   //�������� ������� ���������� �������
#define R_ZEMLI 6371000.
#define G_ZEMLI 9.80665
void fncCalcNormTemperature(const long double valy,long double &valTay, long double &valGradTay)  ;
long double fncCalcOtnositP(const long double H) ;


#endif
