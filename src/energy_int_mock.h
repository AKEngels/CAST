

#ifndef CAST_ENERGY_INT_MOCK_H
#define CAST_ENERGY_INT_MOCK_H

#include "energy.h"
#include "coords.h"

#include "mock_parser.h"

namespace energy::interfaces {
  class mock : public interface_base {
  public:
    explicit mock(coords::Coordinates* coordPointer) : interface_base(coordPointer){
      using namespace ::mock;
      auto tokens = tokenize(Config::get().energy.mock_function);
      func_ = parseTokens(tokens);

      for (std::size_t i = 0; i<3; ++i)
        gradients_.emplace_back(func_->derivative(i)->simplify());
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

    coords::float_type e() override {
      if (coords->size() != 1)
        throw std::runtime_error("Mock interface can only be used for 1 atom.");
      auto x = coords->xyz(0).x();
      auto y = coords->xyz(0).y();
      auto z = coords->xyz(0).z();
      return (*func_)({x, y, z});
    }

    coords::float_type g() override {
      auto x = coords->xyz(0).x();
      auto y = coords->xyz(0).y();
      auto z = coords->xyz(0).z();
      auto grad_x = (*gradients_[0])({x, y, z});
      auto grad_y = (*gradients_[1])({x, y, z});
      auto grad_z = (*gradients_[2])({x, y, z});
      coords::Representation_3D grad{coords::Cartesian_Point (grad_x, grad_y, grad_z)};
      coords->set_g_xyz(grad);
      return e();
    }

    coords::float_type h() override {
      throw std::runtime_error("Mock interface does not support Hessian calculation");
    }

    coords::float_type o() override {
      throw std::runtime_error("Mock interface dows not support geometry optimization");
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

  private:
    std::unique_ptr<::mock::elements::Base> func_;
    std::vector<std::unique_ptr<::mock::elements::Base>> gradients_;
  };
}

#endif //CAST_ENERGY_INT_MOCK_H
