//---------------------------------------------------------------------------

#ifndef DBFTabH
#define DBFTabH
class TDBFTab
{
public:
// к-во записей в таблице
	int m_NumRec;
// к-во полей в таблице
	int m_NumFields;
// суммадлин записей всех полей в таблице (= к-во столбцов в матрице записей)
	int m_lenRow ;
// максимальная длина строки с названием поля
	int m_lenFieldName;
// строка с названиями полей
	char *m_pcharrFieldName;
// матрица со строками контента записей
	char *m_pcharrRecordsContent;
// вектор типа поля
	char *m_pcharrFieldType;
// вектор номера начала i- го поля в строке матрицы  m_pcharrFieldName
	int *m_piarrParts;
// вектор количества знаков после запятой у поля
	int *m_piarrDecCount ;
	// деструктор

	__fastcall ~TDBFTab() ;
	// конструктор по умолчанию
	TDBFTab () ;
	// конструктор копирования
	TDBFTab  (const TDBFTab  &R) ;
	// оператор присваивания
    TDBFTab  & operator=(const TDBFTab   &R2) ;
	// парам констр
    TDBFTab  (  wchar_t *wFilename);
	// запись DBF файла
void WriteDBFFile( wchar_t *wFilename);


} ;
#endif
