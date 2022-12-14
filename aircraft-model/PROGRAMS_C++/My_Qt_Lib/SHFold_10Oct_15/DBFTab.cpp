//---------------------------------------------------------------------------




 #include <stdlib.h>
 #include <string.h>
#include "DBFTab.h"
#include "stdio.h"

__fastcall TDBFTab::~TDBFTab()
{
	if(m_pcharrFieldName) delete m_pcharrFieldName ;
	m_pcharrFieldName = NULL ;
	if(m_pcharrFieldType) delete m_pcharrFieldType ;
	m_pcharrFieldType = NULL ;
	if(m_piarrParts) delete m_piarrParts ;
	m_piarrParts = NULL ;
	if(m_pcharrRecordsContent) delete m_piarrParts ;
	m_pcharrRecordsContent = NULL ;
	if(m_piarrDecCount) delete m_piarrDecCount ;
	m_piarrDecCount = NULL ;


}
//---------------------------------------------------------------------------


TDBFTab::TDBFTab()
{
	m_NumRec = 0 ;
	m_NumFields = 0 ;
	m_lenRow = 0 ;
	m_lenFieldName = 11;
	m_pcharrFieldName = NULL ;
	m_pcharrFieldType = NULL ;
	m_piarrParts = NULL ;
	m_pcharrRecordsContent = NULL ;
	m_piarrDecCount = NULL ;
}

//---------------------------------------------------------------------------


// конструктор копирования
 TDBFTab ::TDBFTab (const TDBFTab &R)
 {
  //
	m_NumRec  = R.m_NumRec ;
	m_lenRow = R.m_lenRow ;
	m_NumFields = R.m_NumFields ;
	m_lenFieldName = R.m_lenFieldName ;
  // ьатрица записей
	m_pcharrRecordsContent = NULL;

	if(R.m_pcharrRecordsContent != NULL)
	{

		m_pcharrRecordsContent = new char[m_lenRow * m_NumRec];

		if(m_pcharrRecordsContent == NULL)
		{
        //ShowMessage(L"Not memory for m_pcharrRecordsContent") ;
		//Abort() ;
		}

		memcpy( m_pcharrRecordsContent,R.m_pcharrRecordsContent, m_lenRow * m_NumRec * sizeof(char));
	}
  //
  // строка с названиями полей
	m_pcharrFieldName = NULL;

	if(R.m_pcharrFieldName != NULL)
	{

		m_pcharrFieldName = new char [m_lenFieldName * m_NumFields];

		if(m_pcharrFieldName == NULL)
		{
        //ShowMessage(L"Not memory for m_pcharrFieldName") ;
		//Abort() ;
		}

		memcpy( m_pcharrFieldName,R.m_pcharrFieldName, m_lenFieldName * m_NumFields* sizeof(char));
	}
  //
  // // вектор типа поля

	m_pcharrFieldType = NULL;

	if(R.m_pcharrFieldType != NULL)
	{

		m_pcharrFieldType = new char[m_NumFields];

		if(m_pcharrFieldType == NULL)
		{
        //ShowMessage(L"Not memory for m_pcharrFieldType") ;
		//Abort() ;
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
        //ShowMessage(L"Not memory for m_pcharrFieldName") ;
		//Abort() ;
		}

		memcpy(m_piarrParts,R.m_piarrParts, m_NumFields * sizeof(int));
	}

	// вектор знаков после запятой
	m_piarrDecCount = NULL;

	if(R.m_piarrDecCount != NULL)
	{

		m_piarrDecCount = new int[m_NumFields];

		if(R.m_piarrDecCount == NULL)
		{
        //ShowMessage(L"Not memory for m_piarrDecCount") ;
		//Abort() ;
		}

		memcpy(m_piarrDecCount,R.m_piarrDecCount, m_NumFields * sizeof(int));
	}

 }
 // оператор присваивания
 TDBFTab &TDBFTab::operator=(const TDBFTab  &R)
 {
  //
	m_NumRec  = R.m_NumRec ;
	m_lenRow = R.m_lenRow ;
	m_NumFields = R.m_NumFields ;
	m_lenFieldName = R.m_lenFieldName ;
  // ьатрица записей
	m_pcharrRecordsContent = NULL;

	if(R.m_pcharrRecordsContent != NULL)
	{

		m_pcharrRecordsContent = new char[m_lenRow * m_NumRec];

		if(m_pcharrRecordsContent == NULL)
		{
        //ShowMessage(L"Not memory for m_pcharrRecordsContent") ;
		//Abort() ;
		}

		memcpy( m_pcharrRecordsContent,R.m_pcharrRecordsContent, m_lenRow * m_NumRec * sizeof(char));
	}
  //
  // строка с названиями полей
	m_pcharrFieldName = NULL;

	if(R.m_pcharrFieldName != NULL)
	{

		m_pcharrFieldName = new char [m_lenFieldName * m_NumFields];

		if(m_pcharrFieldName == NULL)
		{
        //ShowMessage(L"Not memory for m_pcharrFieldName") ;
		//Abort() ;
		}

		memcpy( m_pcharrFieldName,R.m_pcharrFieldName, m_lenFieldName * m_NumFields* sizeof(char));
	}
  //
  // // вектор типа поля

	m_pcharrFieldType = NULL;

	if(R.m_pcharrFieldType != NULL)
	{

		m_pcharrFieldType = new char[m_NumFields];

		if(m_pcharrFieldType == NULL)
		{
        //ShowMessage(L"Not memory for m_pcharrFieldType") ;
		//Abort() ;
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
        //ShowMessage(L"Not memory for m_pcharrFieldName") ;
		//Abort() ;
		}

		memcpy(m_piarrParts,R.m_piarrParts, m_NumFields * sizeof(int));
	}
    // вектор знаков после запятой
	m_piarrDecCount = NULL;

	if(R.m_piarrDecCount != NULL)
	{

		m_piarrDecCount = new int[m_NumFields];

		if(R.m_piarrDecCount == NULL)
		{
        //ShowMessage(L"Not memory for m_piarrDecCount") ;
		//Abort() ;
		}

		memcpy(m_piarrDecCount,R.m_piarrDecCount, m_NumFields * sizeof(int));
	}
	return *this ;
 }

// параметр конструктор создания из файла
TDBFTab::TDBFTab( wchar_t *wFilename)
{
   m_lenFieldName =  11 ;
   FILE *fr ;
   if (	(fr=_wfopen(wFilename,L"rb")) == NULL)
   {
	   //String s_22 =  L"TDBFTab::TDBFTab( wchar_t *wFilename)\nNot possible to open file" ;
	  // s_22 = s_22 +  wFilename ;
      // //ShowMessage(s_22) ;
   }
  long offset =4;

 fseek(fr,offset,SEEK_SET);
   fread(&m_NumRec ,4, 1,fr) ;

  short lenHdr = 0 ;
   fread(&lenHdr ,2, 1,fr) ;
  const int ilenHdr =  lenHdr ;
  m_NumFields = ilenHdr/32 -1 ;

  m_piarrParts = new int [m_NumFields ] ;

  m_piarrDecCount = new int [m_NumFields ] ;

  m_pcharrFieldName = new char[m_lenFieldName * m_NumFields] ;//!!!!!!! m_NumFields

  m_pcharrFieldType = new char [m_NumFields ];
  short lenRecord = 0 ;
 // st = fread(&lenRecord ,2, 1,fr) ; 11.01.2016
  m_lenRow = 0 ;
  const int ilenRecord = lenRecord ;
  for (int i  =0 ; i < m_NumFields; i++)
  {
	m_piarrParts [i] = m_lenRow ;
	offset = 32 + i * 32 ;
	fseek(fr,offset,SEEK_SET);

	 fread(& m_pcharrFieldName[i * m_lenFieldName] ,1, m_lenFieldName,fr) ;
	 fread(& m_pcharrFieldType[i ] ,1, 1,fr) ;

	offset = 32 + i * 32 + 16;
	fseek(fr,offset,SEEK_SET);
	char lenField = 0 ;
	 fread(& lenField ,1, 1,fr) ;
	int ilenField =  lenField;
	m_lenRow +=  ilenField;
	char chDecCount = 0;
	 fread(& chDecCount ,1, 1,fr) ;
	m_piarrDecCount[ i ] =  chDecCount ;
  }

   // Создание таблицы(матрицы) записей
   m_pcharrRecordsContent = new char [ m_NumRec *  m_lenRow ] ;
  // Extract records
   for (int i = 0; i < m_NumRec; i++)
   {
	 offset = ilenHdr + i * ilenRecord ;
	fseek(fr,offset,SEEK_SET);
	 char delFlag = 0 ;
	  fread(& delFlag,1,1,fr);
	 if (delFlag == '\x20')
	 {
		fread(&m_pcharrRecordsContent[ i * m_lenRow],1,m_lenRow,fr);
	 }
   }
   fclose(fr) ;
}
// запись DBF файла
void TDBFTab::WriteDBFFile( wchar_t *wFilename)
{
	FILE  *fw ;
	const int Magic_Number1 = 0 ;
	fw=_wfopen(wFilename,L"wb");
  //	if(!fw) //////ShowMessage (L"TYrWrite::PutPointsToCsvFile\nFile is not opened !") ;
	//1. version number and Date of last update  - первые 0-3  байта
	// char ch = '\x03';
	int i0 = 201421059;
	fwrite(&i0,sizeof(int),1 ,fw) ;
	// 2. Number of records  in data file   4-7 байты
	fwrite(&m_NumRec,sizeof(int),1 ,fw) ;
	// 3. длина header
	int ilenHdr = 32 * ( 1 + m_NumFields ) +1;
	short lenHdr = ilenHdr;
	fwrite(&lenHdr,sizeof(short),1 ,fw) ;
	// 4. длина записи
	short lenRow = m_lenRow + 1 ;
	fwrite(&lenRow,sizeof(short),1 ,fw) ;
	// 5. волшебная строка 1
	int iarrMagic1[] = {0,0,0,0,22272};
	fwrite(iarrMagic1,sizeof(int),5 ,fw) ;

	for (int i = 0; i < m_NumFields ; i++)
	{
		fwrite(&m_pcharrFieldName [ i * m_lenFieldName],1,m_lenFieldName ,fw) ;
		fwrite(&m_pcharrFieldType [ i ],1, 1,fw) ;
		int numb1 = Magic_Number1;
		fwrite(&numb1 ,4, 1,fw) ;
		int iFieldLen = (i == (m_NumFields - 1))? (m_lenRow - m_piarrParts [ m_NumFields -1]): (m_piarrParts [i + 1]- m_piarrParts [ i]) ;
		char chFieldLen = iFieldLen;
		fwrite(&chFieldLen,1, 1,fw) ;
		char chDecCount = m_piarrDecCount[ i ] ;
		fwrite(&chDecCount,1, 1,fw) ;
		char charrMagic1[] = { '\x00','\x00'} ;
		fwrite(charrMagic1,1, 2,fw) ;
		int iarrMagic2[] = {0,0,1459617792};
		fwrite(iarrMagic2,4, 3,fw) ;

	}
	int chTerminator = '\x0D';
	fwrite(&chTerminator,sizeof(char),1 ,fw) ;

	// запмись контента
	char delFlag =  '\x20';
	for (int i = 0; i < m_NumRec; i++)
	{
	fwrite(&delFlag ,1, 1,fw) ;
	fwrite(&m_pcharrRecordsContent[ i * m_lenRow] ,1, m_lenRow,fw) ;
	}


	fclose(fw) ;
}

