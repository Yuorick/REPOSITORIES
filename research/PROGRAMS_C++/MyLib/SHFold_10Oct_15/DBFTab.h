//---------------------------------------------------------------------------

#ifndef DBFTabH
#define DBFTabH
class TDBFTab
{
public:
// �-�� ������� � �������
	int m_NumRec;
// �-�� ����� � �������
	int m_NumFields;
// ��������� ������� ���� ����� � ������� (= �-�� �������� � ������� �������)
	int m_lenRow ;
// ������������ ����� ������ � ��������� ����
	int m_lenFieldName;
// ������ � ���������� �����
	char *m_pcharrFieldName;
// ������� �� �������� �������� �������
	char *m_pcharrRecordsContent;
// ������ ���� ����
	char *m_pcharrFieldType;
// ������ ������ ������ i- �� ���� � ������ �������  m_pcharrFieldName
	int *m_piarrParts;
// ������ ���������� ������ ����� ������� � ����
	int *m_piarrDecCount ;
	// ����������

	__fastcall ~TDBFTab() ;
	// ����������� �� ���������
	TDBFTab () ;
	// ����������� �����������
	TDBFTab  (const TDBFTab  &R) ;
	// �������� ������������
	TDBFTab  &operator=(const TDBFTab   &R2) ;
	// ����� ������
	TDBFTab  ( wchar_t *wFilename);
	// ������ DBF �����
void WriteDBFFile( wchar_t *wFilename);


} ;
#endif
