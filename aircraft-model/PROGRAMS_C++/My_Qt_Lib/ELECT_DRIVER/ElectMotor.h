#ifndef ELECTMOTOR_H
#define ELECTMOTOR_H
// 1.вычисление индуктивности обмотки:
// L = Tэл.маг * R
// Tэл.маг - электромагнитная постоянная времени (=0.003), R - сопротвление


class QElectMotor
{
protected:
    // индуктивность катушки статора
     double mInductL;
     // момент инерции ротора
     double mJ0;
 private:
    // момент сопротивления при обесточенных обмотках, не более
    // относительно остаточного момента принята следующая гипотеза:
    // остаточный момент является суммой трех слагаемых
    // 1. момент сухого трения - 0.01 н/м
    // 2. знакопеременная функция (гармоническая)
    // - Fgarmon(Tetta, Omega) = Agarmon(Omega) * cos (Zp * Tetta + GrmnPh0)
    // Agarmon(Omega) = (0.9 * mMomResidual) * Omega /OmegaMax - амплитуда
    // GrmnPh0 - фаза
    // здесь Tetta - угл поворота ротора
    // Fgarmon(Tetta) достигает максимума, когда угол ротора совпадает с одной из катушек статора
    // амплитуда Fgarmon(Tetta) пропорциональна угловой скорости вращения ротора
    // 3. белый шум с интенсивностью:
    // Sig2W(Omega) = ( 1. - 0.9 * 0.9) * (mMomResidual * Omega /OmegaMax)* (mMomResidual * Omega /OmegaMax) / 9
    // Таким образом, при нулевой угловой скорости присутствует только лишь момент сухого трения в подшипнике
    //

    // максимальная амплитуда гармонической составляющей остаточного момента
    double mMaxAmpGrmn;
    // фаза гармонической составляющей
    double mGrmnPh0;
    // максимальная дисперсия шума остаточного момента теоретическая
    double mMaxDispRezid;

    // максимальная дисперсия шума остаточного момента практическая
    double mPracticalMaxDispRezid;

    // момент сухого трения
    double mMomDryFriction;


public:
     QElectMotor();

    // конструктор копирования
      QElectMotor (const  QElectMotor &R);

     // оператор присваивания
      QElectMotor  &operator=( const QElectMotor  &R);



      // парам конструктор
     // QElectMotor (const double  InductL , const double  Resist ,
                  // const double  Psi_f , const double   JRotor , const double MaxAmpGrmn
                 //  , const double  GrmnPh0,const double MaxDispMomRezid, const double PracticalMaxDispRezid);

      QElectMotor (const double  InductL , const double  Resist ,
               const double  Psi_f , const double   JRotor , const double MaxAmpGrmn
               , const double  GrmnPh0,const double MaxDispMomRezid, const double PracticalMaxDispRezid, const double MomDryFriction);



     // сопротивление катушки статора
      double mResist;

      // индуктивный поток
      double mPsi_f;

      // к-во пар
      double mZp;

      // величина, обратная индуктивности катушки статора
       double mInvL;

      // предельно допустимое напряжение
       double mUMax;
      // предельно допустимая амплитуда тока
       double mAMax;
      // предельно допустимая частота вращения
       double mOmegaMax;



       static double     calc_Psi_f (const double  VAlZp, const double   VAlCoeff );

       double     calcSumMomNoise (const double VAlTetta, const double VAlOmega);

       double SIGNUM_(const double a);

       double    calcDisp_MomNoise ( const double VAlOmega);

       double    getMaxHarmonAmp ();

       double    getHarmonPh ();

      // void   setMomResidual (const double MaxHarmonAmp);

       void    setHarmonPh (const double GrmnPh0);

        double  calc_dHarmMom_po_dOm (const double VAlTetta);

        double  calc_dHarmMom_po_dTetta (const double VAlTetta, const double VAlOm);

        double  calc_dHarmMom_po_dAmp (const double VAlTetta, const double VAlOm);

        double  calc_dHarmMom_po_dPh (const double VAlTetta, const double VAlOm);

        double  calc_HarmMom (const double VAlTetta, const double VAlOm);

        double   getLInv ();

        void   setLInv (const double Linv);

        double  getL ();

        double  getJ0 ();

        void   setMaxAmpGrmn (const double MaxAmpGrmn);

        static double calc_MaxDispMomRezid (const double MaxAmpGrmn);

        void   setMaxDispMomRezid (const double MaxDispMomRezid);

        double getMomDryFriction ();

        double calc_AmpHarmMom (const double VAlTetta, const double VAlOm);

        double calcDisp_FluctMomNoise ( const double VAlOmega);

        int createInputDataReport(wchar_t*FileName, const bool bHeader);


};
#endif // ELECTMOTOR_H
