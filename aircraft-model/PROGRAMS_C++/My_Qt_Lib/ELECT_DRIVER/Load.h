#ifndef LOAD_H
#define LOAD_H

class TEnvironment;
class QSegment;
class TPlane;

class QLoad
{
public:
    QLoad();

    // конструктор копирования
      QLoad (const  QLoad &R);

     // оператор присваивания
      QLoad  &operator=( const QLoad  &R);



      // парам конструктор 1
      QLoad ( const double  JPayLoad, const double  COmx);
      // парам конструктор 2
      QLoad ( const double  JRotor, const double  Cx, const double  Cv,const double  Max_dOm_po_dt);
      // парам конструктор 3




      //момент инерции нагрузки
      double mJPayLoad;
      // коэффиц сопротивления воздуха
      double mCx;
      // коэффиц трения воздуха
      double mCv;

      // ограничение по угловой скорости
      double mMax_dOm_po_dt;




      virtual void  setLength(const double  VAlLength);

      virtual double  getLength();

      virtual   double  calcAirResistMom(  TEnvironment &Environment,const double ValTet, const double VAlOm);

      inline double sign_(const double x)
      {
          return (x > 0.)?1.:-1;
      }

      virtual double calc_dAirResistMom_po_dOmega( TEnvironment &Environment,const double ValTet
                                                   , const double VAlOm);

      virtual double calc_dMa_po_Cx(const double VAlOm);

      double calc_dMa_po_Cv(const double VAlOm);

      double calc_d2Ma_po_dOm_po_Cv(const double VAlOm);

      virtual double calc_d2Ma_po_dOm_po_Cx(const double VAlOm);

      virtual int createInputDataReport(wchar_t*FileName, const bool bHeader);
};


#endif // LOAD_H
