#include "tabulated_rate.h"
#include "../util/io.h"
#include "../util/constants.h"
#include "../core/radial.h"

TabulatedRate::TabulatedRate(const std::string& filename, double density_of_neutrals,
                       double z, double ionizing_fraction, Medium::Pressure pres)
  :p(pres), density_of_neutrals(0), density_of_neutrals_0(density_of_neutrals), ionizing_fraction(ionizing_fraction){
  IO::read(filename, intensity_values, rate_values);
   
  density_of_neutrals = p.get_pressure(z) * density_of_neutrals_0; 

  // interpolation
  spline = gsl_spline_alloc(gsl_interp_linear, intensity_values.size());
  gsl_spline_init(spline, intensity_values.data(), rate_values.data(), intensity_values.size());
  acc = gsl_interp_accel_alloc();
}

TabulatedRate::~TabulatedRate() {
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline);
}

void TabulatedRate::update_density_of_neutrals(double z)
{
  density_of_neutrals = p.get_pressure(z) * density_of_neutrals_0; 
}

void TabulatedRate::calculate_electron_density(const Radial& electric_field,
                                            Array2D<double>& ionization_rate,
                                            Array2D<double>& electron_density, double z) {
  update_density_of_neutrals(z);
  double dt = electric_field.time[1] - electric_field.time[0];
  double eta = ionizing_fraction * density_of_neutrals * dt / 2;
  for (int i = 0; i < electric_field.Nradius; ++i) {
    double E = electric_field.rt(i, 0).real();
    double I = 0.5 * Constants::epsilon_0*Constants::c * std::pow(E, 2);
    ionization_rate(i, 0) = rate(I);
    electron_density(i, 0) = eta * ionization_rate(i, 0);
    for (int j = 1; j < electric_field.Ntime; ++j) {
      double E = electric_field.rt(i, j).real();
      double I = 0.5 * Constants::epsilon_0*Constants::c * std::pow(E, 2);
      ionization_rate(i, j) = rate(I);
      double exp = std::exp(-(ionization_rate(i, j) + ionization_rate(i, j-1)) / 2 * dt);
      electron_density(i, j) = exp * (electron_density(i, j-1) + eta * ionization_rate(i, j-1)) + eta * ionization_rate(i, j);
    }
  }
}

