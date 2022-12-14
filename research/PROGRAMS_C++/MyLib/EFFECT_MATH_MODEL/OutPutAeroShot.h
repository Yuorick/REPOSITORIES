//---------------------------------------------------------------------------

#ifndef OutPutAeroShotH
#define OutPutAeroShotH
// ????? ????????? ?????? ???????? ?????????? ?? ?????????? ???????? ?? ????????? ????
// ??? ?????? ??? ???????????? ???????????
class TOutPutAeroShot
{
public:
//   1 ????? ????????
double mTshot;
//   2 ???????? ?????
double mTfly;
//  3 ?????? ????????? ??????? ?? ?????? ???????? ? ???
double marrVS_GSK_Vessel[6];
// 4 ?????? ????????? ???? ??????? ?? ?????? ???????? ? ???
double marrVS_GSK_Targ[6];
//  5 ???? ? ????
double mEpsKGSK;
// 6 ???? ? ????
double mBetKGSK;
// 7 ???????? ??????? ???????? ??????? ????????? ???? ? ??? ?? ?????? ???????
double marrKTarg[36];
// 8 ???????? ??????? ???????? ??????? ????????? ??????? ? ???  ?? ?????? ???????
double marrKShell[36];

// 9 ?????? ???????
double marrMiss[6];


	// ??????????? ?? ?????????
	TOutPutAeroShot () ;
	// ??????????? ???????????
	TOutPutAeroShot  (const TOutPutAeroShot  &R) ;

	// ???????? ????????????
	TOutPutAeroShot  &operator=(const TOutPutAeroShot   &R2) ;
   // ???????????????? ???????????
   TOutPutAeroShot( const double Tshot
		,const double Tfly, double *arrVS_GSK_Vessel, double *arrVS_GSK_Targ
		,const double EpsKGSK,const double BetKGSK, double *arrKTarg
		,double *arrKShell, double *arrMiss);

   static void   swap( TOutPutAeroShot *a0,  TOutPutAeroShot *a1);

   static void   flipArray( TOutPutAeroShot *parr, const int len) ;


}  ;

#endif
