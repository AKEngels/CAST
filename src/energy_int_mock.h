

#ifndef CAST_ENERGY_INT_MOCK_H
#define CAST_ENERGY_INT_MOCK_H

#include "energy.h"
#include "coords.h"

namespace energy::interfaces {
  class mock : public interface_base {
  public:
    mock(coords::Coordinates* coordPointer) : interface_base(coordPointer){
    }

    void swap(interface_base& other) override {
      interface_base::swap(dynamic_cast<mock&>(other));
    }

    interface_base* clone(coords::Coordinates* coord_object) const override {
      return new mock(coord_object);
    }

    interface_base* move(coords::Coordinates* coord_object) override {
      return clone(coord_object);
    }

    void update(bool const) override {}

    coords::float_type e(void) override {
      if (coords->size() != 1)
        throw std::runtime_error("Mock interface can only be used for 1 atom.");
      auto x = coords->xyz(0).x();
      auto y = coords->xyz(0).y();
      auto z = coords->xyz(0).z();
      return  5*( - exp(-(pow((x-4)*0.25, 2) + pow((y-3)*0.3, 2))) - 0.9 *exp(-(pow((x+3)*0.25, 2) + pow((y-4)*0.3, 2))) - exp(-(pow((x+4)*0.25, 2) + pow((y+2)*0.3, 2))) - 0.9* exp(-(pow((x-3)*0.25, 2) + pow((y+3)*0.3, 2)))) + 0.1*z*z;
    }

    coords::float_type g(void) override {
      auto x = coords->xyz(0).x();
      auto y = coords->xyz(0).y();
      auto z = coords->xyz(0).z();
      auto grad_x = -4.5*(0.375 - 0.125*x)*exp(-0.5625*pow(0.333333333333333*x - 1, 2) - 0.81*pow(0.333333333333333*y + 1, 2)) - 5*(0.5 - 0.125*x)*exp(-pow(0.25*x - 1.0, 2) - 0.81*pow(0.333333333333333*y - 1, 2)) - 5*(-0.125*x - 0.5)*exp(-pow(0.25*x + 1.0, 2) - 0.36*pow(0.5*y + 1, 2)) - 4.5*(-0.125*x - 0.375)*exp(-0.5625*pow(0.333333333333333*x + 1, 2) - 1.44*pow(0.25*y - 1, 2));
      auto grad_y = -5*(0.54 - 0.18*y)*exp(-pow(0.25*x - 1.0, 2) - 0.81*pow(0.333333333333333*y - 1, 2)) - 4.5*(0.72 - 0.18*y)*exp(-0.5625*pow(0.333333333333333*x + 1, 2) - 1.44*pow(0.25*y - 1, 2)) - 4.5*(-0.18*y - 0.54)*exp(-0.5625*pow(0.333333333333333*x - 1, 2) - 0.81*pow(0.333333333333333*y + 1, 2)) - 5*(-0.18*y - 0.36)*exp(-pow(0.25*x + 1.0, 2) - 0.36*pow(0.5*y + 1, 2));
      auto grad_z = 0.2*z;
      coords::Representation_3D grad{coords::Cartesian_Point (grad_x, grad_y, grad_z)};
      coords->set_g_xyz(grad);
      return e();
    }

    coords::float_type h(void) override {
      // No hessian for now, sorry
      return e();
    }

    coords::float_type o(void) override {
      return e();
    }

    std::vector<coords::float_type> charges() const override {
      return {0.};
    }

    coords::Gradients_3D get_g_ext_chg() const override {
      return coords::Gradients_3D();
    }

    void print_E(std::ostream& S) const override {
      S << "Total Energy:      ";
      S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << energy;
    }

    void print_E_head(std::ostream& S, bool const endline) const override {
      S << "Energies\n";
      S << std::right << std::setw(24) << "E_bs";
      S << std::right << std::setw(24) << "E_coul";
      S << std::right << std::setw(24) << "E_rep";
      S << std::right << std::setw(24) << "SUM\n";
      if (endline) S << '\n';
    }

    void print_E_short(std::ostream& S, bool const endline) const override {
      S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
      S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
      S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
      S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << energy << '\n';
      if (endline) S << '\n';
    }

    void to_stream(std::ostream& S) const override {}
  };
}

#endif //CAST_ENERGY_INT_MOCK_H
