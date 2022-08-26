#include "CtrlPos.h"
#include "CtrlFollow.h"

QCtrlPos::QCtrlPos()
{

}
#include <math.h>
#include  <string.h>
#include "MatrixProccess.h"
#include "Equations.h"
#include "Comp.h"

// конструктор копирования
 QCtrlPos::QCtrlPos(const QCtrlPos &R):QCtrlFollow( R)
 {

 }

 // оператор присваивания
  QCtrlPos  &QCtrlPos::operator=( const QCtrlPos  &R)
 {
      if(this == &R)
      {
          return *this;
      }
     QCtrlFollow:: operator= (R);

     return *this ;
 }

//-----------------------------------
// парам конструктор
QCtrlPos::QCtrlPos( const QElectMotor ElectMotor ,  QLoad*  Load, const double *arrSpreadParams
                    , const double MomOut, const double VAlTettaBegin,const double VAlOmegaBegin
                    , const double T0,const double TCur
                    ,const double valh,TComp *pCmpArrLamb)
  :QCtrlFollow (ElectMotor ,  Load, arrSpreadParams, MomOut, VAlTettaBegin
                ,VAlOmegaBegin ,T0,  TCur,valh,pCmpArrLamb)
  {
   // marrObjective[1] = 0.;
   // marrObjective[0] = arrObjective[0];
  }



