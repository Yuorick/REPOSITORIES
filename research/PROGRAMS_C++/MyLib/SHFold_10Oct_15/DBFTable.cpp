//---------------------------------------------------------------------------


//#pragma hdrstop
 //#include <vcl.h>
#include "DBFTable.h"




#include <vcl.h>

#pragma hdrstop
 #include <stdio.h>
 #include <math.h>

 #include <stdlib.h>
//---------------------------------------------------------------------------

#pragma package(smart_init)

 
//---------------------------------------------------------------------------


TDBFTable::TDBFTable()
{
	m_NumRec = 0 ;
	m_NumFields = 0 ;
	m_lenRow = 0 ;
	m_pcharrFieldName = NULL ;
	m_pcharrFieldType = NULL ;
	m_piarrParts = NULL ;
}

//---------------------------------------------------------------------------


// конструктор копирования
 TDBFTable ::TDBFTable (const TDBFTable &R)
 {
  //
	m_NumRec  = R.m_NumRec ;
	m_lenRow = R.m_lenRow ;
	m_NumFields = R.m_NumFields ;
  // ьатрица записей
	m_pcharrFieldName = NULL;

	if(R.m_pcharrFieldName != NULL)
	{

		m_pcharrFieldName = new char[m_lenRow * m_NumRec];

		if(m_pcharrFieldName == NULL)
		{
		ShowMessage(L"Not memory for m_pcharrFieldName") ;
		Abort() ;
		}

		memcpy( m_pcharrFieldName,R.m_pcharrFieldName, m_lenRow * m_NumRec * sizeof(char));
	}
  //
  // // вектор типа поля

	m_pcharrFieldType = NULL;

	if(R.m_pcharrFieldType != NULL)
	{

		m_pcharrFieldType = new char[m_NumFields];

		if(m_pcharrFieldType == NULL)
		{
		ShowMessage(L"Not memory for m_pcharrFieldType") ;
		Abort() ;
		}

		memcpy( m_pcharrFieldType,R.m_pcharrFieldType, m_NumFields * sizeof(char));
	}
   // вектор номера начала i- го поля в строке матрицы  m_pcharrFieldName

	m_piarrParts = NULL;

	if(R.m_piarrParts != NULL)
	{

		m_piarrParts = new int[m_NumFields];

		if(R.m_piarrParts == NULL)
		{
		ShowMessage(L"Not memory for m_pcharrFieldName") ;
		Abort() ;
		}

		memcpy(m_piarrParts,R.m_piarrParts, m_NumFields * sizeof(int));
	}

 }
 // оператор присваивания
 TDBFTable &TDBFTable::operator=(const TDBFTable  &R)
 {
  //
	m_NumRec  = R.m_NumRec ;
	m_lenRow = R.m_lenRow ;
	m_NumFields = R.m_NumFields ;
  // ьатрица записей
	m_pcharrFieldName = NULL;

	if(R.m_pcharrFieldName != NULL)
	{

		m_pcharrFieldName = new char[m_lenRow * m_NumRec];

		if(m_pcharrFieldName == NULL)
		{
		ShowMessage(L"Not memory for m_pcharrFieldName") ;
		Abort() ;
		}

		memcpy( m_pcharrFieldName,R.m_pcharrFieldName, m_lenRow * m_NumRec * sizeof(char));
	}
  //
  // // вектор типа поля

	m_pcharrFieldType = NULL;

	if(R.m_pcharrFieldType != NULL)
	{

		m_pcharrFieldType = new char[m_NumFields];

		if(m_pcharrFieldType == NULL)
		{
		ShowMessage(L"Not memory for m_pcharrFieldType") ;
		Abort() ;
		}

		memcpy( m_pcharrFieldType,R.m_pcharrFieldType, m_NumFields * sizeof(char));
	}
 // вектор номера начала i- го поля в строке матрицы  m_pcharrFieldName

	m_piarrParts = NULL;

	if(R.m_piarrParts != NULL)
	{

		m_piarrParts = new int[m_NumFields];

		if(R.m_piarrParts == NULL)
		{
		ShowMessage(L"Not memory for m_pcharrFieldName") ;
		Abort() ;
		}

		memcpy( m_piarrParts,R.m_piarrParts, m_NumFields * sizeof(int));
	}

	return *this ;
 }

 // функции-члены
 /*
int  TDBFTable::LenString()
{
  int len = 0;
  for (int i =0; i < m_NumFields; i++)
  {
	int jj = m_pcharrFieldLeng[i] ;
	len +=jj;
  }
}  */


//#pragma package(smart_init)
