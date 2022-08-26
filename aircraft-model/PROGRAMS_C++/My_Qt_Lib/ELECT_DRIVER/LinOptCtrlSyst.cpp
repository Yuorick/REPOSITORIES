#include "LinOptCtrlSyst.h"
#include "MatrixProccess.h"
#include  <string.h>
#include <math.h>
#include <float.h>
#include "Equations.h"


#include "Comp.h"



QLinOptCtrlSyst::~QLinOptCtrlSyst()
{
    if(mpA) delete []mpA ;
    mpA = NULL ;
    if(mpB) delete []mpB ;
    mpB = NULL ;
    if(mpC) delete []mpC ;
    mpC = NULL ;
    if(marrPhVect0) delete []marrPhVect0 ;
    marrPhVect0 = NULL ;
    if(marrPhVect) delete []marrPhVect ;
    marrPhVect = NULL ;
}
//---------------------------------------------------------------------------
QLinOptCtrlSyst::QLinOptCtrlSyst()
{
  mQuantXVar = 0;
  mQuantUVar = 0;
  mpA = NULL;
  mpB= NULL;
  mpC= NULL;
  mT0 = 0.;
  marrPhVect0= NULL;
  mTCur = 0.;
  marrPhVect= NULL;
}


// конструктор копирования
 QLinOptCtrlSyst :: QLinOptCtrlSyst (const  QLinOptCtrlSyst &R)
 {
         mQuantXVar = R.mQuantXVar;
         mQuantUVar = R.mQuantUVar;
         mT0 = R.mT0;
         mTCur = R.mTCur;

         if (R.mpA != NULL)
         {
           if (!(mpA = new double  [mQuantXVar * mQuantXVar]))
               abort();
         }

         if (R.mpB != NULL)
         {
           if (!(mpB = new double  [mQuantXVar * mQuantUVar]))
               abort();
         }

         if (R.mpC != NULL)
         {
           if (!(mpC = new double  [mQuantXVar]))
               abort();
         }

         if (R.marrPhVect0 != NULL)
         {
           if (!(marrPhVect0 = new double  [mQuantXVar]))
               abort();
         }

         if (R.marrPhVect != NULL)
         {
           if (!(marrPhVect = new double  [mQuantXVar]))
               abort();
         }

         memcpy(mpA, R.mpA,  mQuantXVar * mQuantXVar * sizeof(double));
         memcpy(mpB, R.mpB,  mQuantXVar * mQuantUVar * sizeof(double));
         memcpy(mpC, R.mpC,              mQuantXVar  * sizeof(double));
         memcpy(marrPhVect0, R.marrPhVect0, mQuantXVar * sizeof(double));
         memcpy(marrPhVect,  R.marrPhVect,  mQuantXVar * sizeof(double));

 }

 // оператор присваивания
  QLinOptCtrlSyst  &QLinOptCtrlSyst::operator=( const QLinOptCtrlSyst  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      mQuantXVar = R.mQuantXVar;
      mQuantUVar = R.mQuantUVar;
      mT0 = R.mT0;
      mTCur = R.mTCur;

      if (R.mpA != NULL)
      {
        if (!(mpA = new double  [mQuantXVar * mQuantXVar]))
            abort();
      }

      if (R.mpB != NULL)
      {
        if (!(mpB = new double  [mQuantXVar * mQuantUVar]))
            abort();
      }

      if (R.mpC != NULL)
      {
        if (!(mpC = new double  [mQuantXVar]))
            abort();
      }

      if (R.marrPhVect0 != NULL)
      {
        if (!(marrPhVect0 = new double  [mQuantXVar]))
            abort();
      }

      if (R.marrPhVect != NULL)
      {
        if (!(marrPhVect = new double  [mQuantXVar]))
            abort();
      }

      memcpy(mpA, R.mpA,  mQuantXVar * mQuantXVar * sizeof(double));
      memcpy(mpB, R.mpB,  mQuantXVar * mQuantUVar * sizeof(double));
      memcpy(mpC, R.mpC,              mQuantXVar  * sizeof(double));
      memcpy(marrPhVect0, R.marrPhVect0, mQuantXVar * sizeof(double));
      memcpy(marrPhVect,  R.marrPhVect,  mQuantXVar * sizeof(double));


     return *this ;
 }


  // парам конструктор
     QLinOptCtrlSyst:: QLinOptCtrlSyst (const   int QuantXVar,
     const int QuantUVar,
     double *pA,
     double *pB,
     double *pC,
     const double T0,
     double *arrPhVect0,
     const double TCur,
     double *arrPhVect)
 {
         mQuantXVar = QuantXVar;
         mQuantUVar = QuantUVar;
         mT0 = T0;
         mTCur = TCur;
         mpA = NULL;
         mpB = NULL;
         mpA = NULL;
         marrPhVect0 = NULL;
         marrPhVect = NULL;

         if (pA != NULL)
         {
           if (!(mpA = new double  [mQuantXVar * mQuantXVar]))
               abort();
         }

         if (pB != NULL)
         {
           if (!(mpB = new double  [mQuantXVar * mQuantUVar]))
               abort();
         }

         if (pC != NULL)
         {
           if (!(mpC = new double  [mQuantXVar]))
               abort();
         }

         if (arrPhVect0 != NULL)
         {
           if (!(marrPhVect0 = new double  [mQuantXVar]))
               abort();
         }

         if (arrPhVect != NULL)
         {
           if (!(marrPhVect = new double  [mQuantXVar]))
               abort();
         }

         memcpy(mpA, pA,  mQuantXVar * mQuantXVar * sizeof(double));
         memcpy(mpB, pB,  mQuantXVar * mQuantUVar * sizeof(double));
         memcpy(mpC, pC,              mQuantXVar  * sizeof(double));
         memcpy(marrPhVect0, arrPhVect0, mQuantXVar * sizeof(double));
         memcpy(marrPhVect,  arrPhVect,  mQuantXVar * sizeof(double));

 }

//-----------------------------------------
