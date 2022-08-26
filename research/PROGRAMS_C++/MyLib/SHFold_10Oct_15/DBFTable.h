//---------------------------------------------------------------------------

#ifndef DBFTableH
#define DBFTableH
class TDBFTable
{
public:
// к-во записей в таблице
	int m_NumRec;
// к-во полей в таблице
	int m_NumFields;
// суммадлин записей всех полей в таблице (= к-во столбцов в матрице записей)
	int m_lenRow ;
// матрица записей полей
	char *m_pcharrFieldName;
// вектор типа поля
	char *m_pcharrFieldType;
// вектор номера начала i- го поля в строке матрицы  m_pcharrFieldName
	int *m_piarrParts;
	// деструктор
	   //	__fastcall ~TURPolygon() ;
  //	__fastcall ~TDBFTable() ;
	// конструктор по умолчанию
	TDBFTable () ;
	// конструктор копирования
	TDBFTable  (const TDBFTable  &R) ;
	// оператор присваивания
	TDBFTable  &operator=(const TDBFTable   &R2) ;
	// парам констр
   //	TDBFTable  ( double *arrBox, const int iNumParts,const int iNumPoints,int *iarrParts
	//,TURPointXY *arrPoints);

}  ;
#endif
