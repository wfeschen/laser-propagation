#ifndef TABULATED_RATE_H_
#define TABULATED_RATE_H_

#include "../nonlinear/ionization.h"
#include "../util/interpolate.h"
#include "../linear/medium.h"

class TabulatedRate : public Ionization {
public:
  TabulatedRate(const std::string& filename, double density_of_neutrals, double z,
             double ionizing_fraction, Medium::Pressure pres);

  ~TabulatedRate();

  double rate(double I) {
    return gsl_spline_eval(spline, I, acc);
  }

  void update_density_of_neutrals(double z);

  void calculate_electron_density(const Radial& electric_field,
                                  Array2D<double>& ionization_rate,
                                  Array2D<double>& electron_density, double z) override;

  double density_of_neutrals, ionizing_fraction;
  Medium::Pressure p;
  double density_of_neutrals_0;
  std::vector<double> intensity_values, rate_values;
  gsl_spline* spline;
  gsl_interp_accel* acc;
};


#endif // TABULATED_RATE_H_
