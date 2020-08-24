#ifndef LINEAR_H_
#define LINEAR_H_

#include <functional>
#include <complex>
#include <iostream>

#include "medium.h"

namespace Linear {
  
class Base {
public:
  Base(Medium::IndexFunction linear_index, Medium::Pressure pres);
  virtual std::complex<double> kz(double kperp, double omega, double z) = 0;
  virtual double group_velocity(double kperp, double omega, double z);
  virtual double gvd(double kperp, double omega, double z);

  Medium::IndexFunction n;
  Medium::Pressure p;
};


class FreeSpace : public Base {
public:
  FreeSpace(Medium::IndexFunction linear_index, Medium::Pressure pres)
    :Base(linear_index, pres) {}
  
  std::complex<double> kz(double kperp, double omega, double z) override;
};


class DiffractionLess : public Base {
public:
  DiffractionLess(Medium::IndexFunction linear_index, Medium::Pressure pres)
    :Base(linear_index, pres) {}

  std::complex<double> kz(double, double omega, double z) override;
};


class Capillary : public Base {
public:
  Capillary(Medium::IndexFunction linear_index, double radius, double cladding_index,
	    Medium::Pressure pres)
    :Base(linear_index, pres), R(radius), nclad(cladding_index) {}

  std::complex<double> kz(double kperp, double omega, double z) override;

private:
  double R, nclad, pressure;
};

} // end namespace Linear
#endif // LINEAR_H_
