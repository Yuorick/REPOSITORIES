//---------------------------------------------------------------------------

#ifndef DetonatorH
#define DetonatorH
enum enumDetonatorType {DVM,AR32A, D4MRM, CONTACT, MFIVU,AR51_LM,BARRIER, DETONATOR_UNKNOWN};

#define LEN_DOUBLE_ARR_DETONATOR_PARAMS 50
#define LEN_INT_ARR_DETONATOR_PARAMS 50
// ����� ��������� �������������� ����������

// ������ �������������� �� ������ ��� 76 ��
// marrDetonatorParams[0] - ������� ��������� ������� ������������
// marrDetonatorParams[1] - ���������� ���������� ������ ������������ ��� �������� �������


// ������ �������� ��51-��
// marrDetonatorParams[0] - ������ ���������������
class TDetonator
{
public:
 // ����������� ����� �������, �
	 enumDetonatorType  mEnumDetonatorType ;
	 double marrDetonatorParams [LEN_DOUBLE_ARR_DETONATOR_PARAMS];
	 int miarrDetonatorParams [LEN_INT_ARR_DETONATOR_PARAMS];

	// ����������� �� ���������
	TDetonator () ;
	// ����������� �����������
	TDetonator  (const TDetonator  &R) ;

	// �������� ������������
	TDetonator  &operator=(const TDetonator   &R2) ;
   // ���������������� �����������
	TDetonator(enumDetonatorType  EnumDetonatorType, double *arrDetonatorParams
 , int *iarrDetonatorParams) ;

   TDetonator(enumDetonatorType  EnumDetonatorType);


}  ;
#endif