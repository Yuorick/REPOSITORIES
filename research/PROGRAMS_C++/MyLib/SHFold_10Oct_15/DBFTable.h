//---------------------------------------------------------------------------

#ifndef DBFTableH
#define DBFTableH
class TDBFTable
{
public:
// �-�� ������� � �������
	int m_NumRec;
// �-�� ����� � �������
	int m_NumFields;
// ��������� ������� ���� ����� � ������� (= �-�� �������� � ������� �������)
	int m_lenRow ;
// ������� ������� �����
	char *m_pcharrFieldName;
// ������ ���� ����
	char *m_pcharrFieldType;
// ������ ������ ������ i- �� ���� � ������ �������  m_pcharrFieldName
	int *m_piarrParts;
	// ����������
	   //	__fastcall ~TURPolygon() ;
  //	__fastcall ~TDBFTable() ;
	// ����������� �� ���������
	TDBFTable () ;
	// ����������� �����������
	TDBFTable  (const TDBFTable  &R) ;
	// �������� ������������
	TDBFTable  &operator=(const TDBFTable   &R2) ;
	// ����� ������
   //	TDBFTable  ( double *arrBox, const int iNumParts,const int iNumPoints,int *iarrParts
	//,TURPointXY *arrPoints);

}  ;
#endif
