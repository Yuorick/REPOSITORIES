#include "Lnsgm.h"
   const double EPS = 0.00000000001;
// класс описывает отрезок , лежащий на прямой
QLnSgm::QLnSgm()
{
    ma = 0.;
    mb = 0.;
}
//-------------------------------

// Конструктор копирования
 QLnSgm::QLnSgm (const QLnSgm &R2)
 {

     ma= R2.ma ;
      mb= R2.mb;

 }
 // парам констр
 QLnSgm::QLnSgm(const double a1,const double b1)
 {
     ma = a1;
     mb = b1;
 }



 // оператор присваивания
 QLnSgm &QLnSgm::operator=(const QLnSgm  &R2)
{
     ma= R2.ma ;
      mb= R2.mb;
     return *this ;
}

//длина
double QLnSgm::calcLeng()
{
    return mb - ma;
}

// пересечение 2-х сегментов
// если пересекаются, то возвращает true, результат в sgmRez
bool QLnSgm::productOfSgms( QLnSgm sgm0, QLnSgm sgm1, QLnSgm* psgmRez)
{
    // сегменты не пересекаются
    if ((sgm1.mb <= sgm0.ma + EPS)||(sgm1.ma >= sgm0.mb -EPS))
    {
     return false;
    }
    ///

    //

    if (sgm0.isValInsideSgm(sgm1.ma))
    {
        psgmRez->ma = sgm1.ma;
        if(sgm0.isValInsideSgm(sgm1.mb))
        {
          psgmRez->mb = sgm1.mb;
        }
        else
        {
            psgmRez->mb = sgm0.mb;
        }
    }
    else
    {
        psgmRez->ma = sgm0.ma;
        if(sgm0.isValInsideSgm(sgm1.mb))
        {
          psgmRez->mb = sgm1.mb;
        }
        else
        {
          psgmRez->mb = sgm0.mb;
        }
    }
    return true;
}

//-----------------------------------------


// вычитание из сегмента sgm0 сегмента sgm1
// Возвращает к-во сегментов разности, результат в parrSgmRez
int QLnSgm::SegmMinusSegm( QLnSgm sgm0, QLnSgm sgm1, QLnSgm* parrSgmRez)
{
    QLnSgm   sgmProduct;
   if (!productOfSgms(  sgm0, sgm1, &sgmProduct))
   {
      // сегменты не пересекаются
       parrSgmRez[0] = sgm0;
       return 1;
   }

   parrSgmRez[0] = QLnSgm(sgm0.ma, sgmProduct.ma);
   parrSgmRez[1] = QLnSgm(sgmProduct.mb, sgm0.mb);
   int ireturn =2;
   if (parrSgmRez[0].calcLeng() < EPS )
   {
     parrSgmRez[0] =parrSgmRez[1];
     if (parrSgmRez[0].calcLeng() < EPS )
     {
         return 0;
     }
     return 1;
   }

   if (parrSgmRez[1].calcLeng() < EPS )
   {
       return 1;
   }
  return ireturn;
}

//-------------------------------------------------
//если число x принадлежит отрезхку, то возвращает true
// в противном случае false
bool  QLnSgm::isValInsideSgm(const double x)
{

   if (( x-ma) * ( x-mb) < -0.00000001)
   {
       return true;
   }
   return false;
}

