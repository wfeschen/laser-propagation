#include "nonlinear_response.h"

class Kerr : public NonlinearResponse {
public:
  Kerr(double n2, double n0, double z, Medium::Pressure pres): NonlinearResponse(pres), n_0(n0), n_2(n2) {
     double pressure = p.get_pressure(z);
     chi3 = 4.0/3.0 * Constants::epsilon_0 * Constants::c * n_2 * pressure * std::pow(n_0, 2);
  }

  void update_chi3(double z){
     double pressure = p.get_pressure(z);
     chi3 = 4.0/3.0 * Constants::epsilon_0 * Constants::c * n_2 * pressure * std::pow(n_0, 2);
  }

  void calculate_response(const std::vector<double>& radius,
                          const std::vector<double>& time,
                          const Array2D<std::complex<double>>& electric_field,
                          const Array2D<double>&,
                          Array2D<std::complex<double>>& response, double z) override {
    update_chi3(z);
    for (std::size_t i = 0; i < radius.size(); ++i) {
      for (std::size_t j = 0; j < time.size(); ++j) {
        double Ereal = electric_field(i, j).real();
        response(i, j) += Constants::epsilon_0 * chi3 * std::pow(Ereal, 3);
      }
    }
  }

private:
  double chi3;
  double n_0;
  double n_2;
};
