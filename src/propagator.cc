#include "propagator.h"
#include "constants.h"
#include "field.h"
#include "radial.h"
#include "util.h"
#include "linear.h"
#include "io.h"
#include <iostream>
#include <iomanip>
#include <limits>

Propagator::Propagator(int Nt, double time_min, double time_max,
                       double wave_min, double wave_max,
		       int Nr, double R, int Nk,
		       double abs_err, double rel_err, double first_step)
  :field(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   Pkerr(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   Jnon(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   Jplasma(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   Rho(Nr, Nt),
   RhoE(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   Ntime(Nt), Nradius(Nr), vg(0), n2(0),
   fraction(1),
   kz(field.Nkperp, field.Nomega),
   coef(field.Nkperp, field.Nomega),
   A(field.Nkperp, field.Nomega),
   z(0), abserr(abs_err), relerr(rel_err), first_step(first_step) {


  Nomega = field.Nomega;
  Nkperp = field.Nkperp;

  std::vector<double> wavelengths;
  for (auto o : field.omega) wavelengths.push_back(2*Constants::pi*Constants::c / o);

  system = {RHSfunction, nullptr, 2*A.vec().size(), this};

  driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, first_step,
  					 abserr, relerr);

  
}

Propagator::~Propagator() {
  gsl_odeiv2_driver_free(driver);
}

std::string Propagator::log_grid_info() {
  // log message about computation box size and parameters
  double wave_min = 2*Constants::pi*Constants::c / field.omega.back();
  double wave_max = 2*Constants::pi*Constants::c / field.omega.front();
  std::stringstream ss;
  ss << "*** Computational Grid ***\n";
  ss << "Supported Wavelengths: (" << wave_min << ", " << wave_max << ")\n";
  ss << "Ntime    =  " << Ntime << "\n";
  ss << "Nomega   = " << Nomega << "\n";
  ss << "Nradius  =  " << Nradius << "\n";
  ss << "Nkperp   = " << Nkperp << "\n";

  // ode solver
  ss << "ODE size = " << A.vec().size() << " complex<double>\n\n";
  return ss.str();
}

void Propagator::initialize_linear(const Linear& linear, double omega0) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    double kperp = field.kperp[i];
    for (int j = 0; j < Nomega; ++j) {
      double omega = field.omega[j];
      kz(i, j) = linear.kz(kperp, omega);
    }
  }
  vg = linear.group_velocity(field.kperp[0], omega0);

  // nonlinear coupling coefficient
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      coef(i, j) = imagi / (2*Constants::epsilon_0*kz(i, j).real()) * std::pow(field.omega[j]/Constants::c, 2);
    }
  }

  IO::write("kz.dat", kz.vec(), Nkperp, Nomega);
  IO::write("coef.dat", coef.vec(), Nkperp, Nomega);
}


void Propagator::initialize_field(const Field::Field& Efield) {
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Ntime; ++j) {
      field.rt(i, j) = Efield(field.radius[i], field.time[j]);
    }
  }
  
  // fill spectral array with initialized field
  field.transform_to_spectral();

  // copy spectral field to auxillary A which is passed to ODE solver
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      A(i, j) = field.kw(i, j);
    }
  }

  field.transform_to_temporal();
}

void Propagator::initialize_kerr(double n2value) {
  n2 = n2value;
}

void Propagator::initialize_pressure(double pressure_atm) {
  pressure = pressure_atm;
}

void Propagator::initialize_rate(const std::string& filename, double frac, 
				 double scale) {
  ionization_rate = std::make_unique<Ionization::Rate>(filename);
  fraction = frac;
  scaling = scale;
}


void Propagator::initialize_filters(double time_filter_min, double time_filter_max,
                                    double wave_filter_min, double wave_filter_max) {
  std::vector<std::reference_wrapper<Radial>> fields = {field, Pkerr, Jnon, Jplasma, RhoE};
  for (auto f : fields) {
    f.get().initialize_temporal_filter(time_filter_min, time_filter_max);
    f.get().initialize_spectral_filter(wave_filter_min, wave_filter_max);
  }
}

void Propagator::linear_step(Radial& radial, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < radial.Nkperp; ++i) {
    for (int j = 0; j < radial.Nomega; ++j) {
      auto arg = kz(i, j) - radial.omega[j] / vg;
      radial.kw(i, j) *= std::exp(-imagi * arg * dz);
    }
  }
}

void Propagator::linear_step(const std::complex<double>* A, Radial& radial, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < radial.Nkperp; ++i) {
    for (int j = 0; j < radial.Nomega; ++j) {
      auto arg = kz(i, j) - radial.omega[j] / vg;
      radial.kw(i, j) = A[i*radial.Nomega + j] * std::exp(-imagi * arg * dz);
    }
  }
}

void Propagator::nonlinear_step(double& z, double zi) {
  int status = gsl_odeiv2_driver_apply(driver, &z, zi,
				       reinterpret_cast<double*>(A.get_data_ptr()));
  if (status != GSL_SUCCESS) {
    throw std::runtime_error("gsl_ode error: " + std::to_string(status));
  }
}

void Propagator::calculate_electron_density() {
  const double rhoneut = 2.5e25 * pressure;
  const double dt = field.time[1] - field.time[0];
  for (int i = 0; i < Nradius; ++i) {
    Util::IntegratorTrapz integrator(dt);
    for (int j = 0; j < Ntime; ++j) {
      double E = field.rt(i, j).real();
  
      // optical field ionization
      double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(E, 2);
      double W = scaling * ionization_rate->operator()(2*I);
      double prob = integrator.add(W);
      double rho = fraction * prob * rhoneut;
      Rho(i, j) = rho;
    }
  }
}

void Propagator::calculate_rhs(double z, const std::complex<double>* A, std::complex<double>* dA) {
  // 1: shift to current z
  linear_step(A, field, z);

  // 2: transform A to E
  field.transform_to_temporal();

  // 3: calculate nonlinearities
  const double rhoneut = 2.5e25 * pressure;
  const double Ui = 2.5250303e-18;
  const double chi3 = 4.0/3.0 * Constants::epsilon_0*Constants::c * n2 * pressure;
  const double dt = field.time[1] - field.time[0];
  for (int i = 0; i < Nradius; ++i) {
    Util::IntegratorTrapz integrator(dt);
    for (int j = 0; j < Ntime; ++j) {
      double E = field.rt(i, j).real();

      // Kerr effect
      double kerr = Constants::epsilon_0 * chi3 * std::pow(E, 3);
      Pkerr.rt(i, j) = kerr;

      // optical field ionization
      double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(E, 2);
      double W = scaling * ionization_rate->operator()(2*I);
      double prob = integrator.add(W);
      double rho = fraction * prob * rhoneut;
      Rho(i, j) = rho;
      RhoE.rt(i, j) = rho * E;
      Jnon.rt(i, j) = W/(I+1) * Ui * (fraction*rhoneut - rho) * Constants::epsilon_0*Constants::c*E;
    }
  }

  // 4: Pkerr(r, t) -> Pkerr(kperp, omega)
  // add other responses to 4-6
  Pkerr.transform_to_spectral();
  RhoE.transform_to_spectral();
  Jnon.transform_to_spectral();

  // 5: linear_step(Pkerr.spectral, -z)
  linear_step(Pkerr, -z);
  linear_step(RhoE, -z);
  linear_step(Jnon, -z);
  
  const double tau = 350e-15 / pressure;
  // 6: dy[] = i*k/(2*epsilon_0) * Pkerr
  std::complex<double> imagi(0, 1);
  double pl1 = std::pow(Constants::e, 2) * tau / Constants::m_e;
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      const double omega = field.omega[j];
      const double ot = omega * tau;
      const std::complex<double> pl2 = (1.0 + imagi*ot) / (1.0 + std::pow(ot, 2));
      Jplasma.kw(i, j) = pl1 * pl2 * RhoE.kw(i, j);
      
      dA[i*Nomega + j] = coef(i, j) * (Pkerr.kw(i, j) + imagi/omega * (Jplasma.kw(i, j) + Jnon.kw(i, j)));
    }
  }
}


SimulationData Propagator::get_data() {
  SimulationData data = {field, Rho};
  return data;
}

int RHSfunction(double z, const double y[], double dy[], void* p) {
  Propagator& prop = *static_cast<Propagator*>(p);
  const std::complex<double>* A = reinterpret_cast<const std::complex<double>*>(y);
  std::complex<double>* dA = reinterpret_cast<std::complex<double>*>(dy);
  prop.calculate_rhs(z, A, dA);
  return GSL_SUCCESS;
}

