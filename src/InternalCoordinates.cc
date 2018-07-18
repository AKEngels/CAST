#include "InternalCoordinates.h"

#include "ic_util.h"

namespace InternalCoordinates {

  coords::float_type BondDistance::val(coords::Representation_3D const& xyz) const {

    auto const& a = xyz.at(index_a_ - 1u);
    auto const& b = xyz.at(index_b_ - 1u);

    return ic_util::euclid_dist<double>(a, b);
  }

  std::pair<coords::r3, coords::r3> BondDistance::der(coords::Representation_3D const& xyz) const {

    auto const& a = xyz.at(index_a_ - 1u);
    auto const& b = xyz.at(index_b_ - 1u);

    auto bond = ic_util::normalize(a - b);
    return std::make_pair(bond, -bond);
  }

  std::vector<coords::float_type>
    BondDistance::der_vec(coords::Representation_3D const& xyz) const {
    using scon::c3;

    auto firstder = der(xyz);
    std::vector<coords::r3> der_vec(xyz.size(), coords::r3(0.0, 0.0, 0.0));
    der_vec.at(index_a_ - 1) = firstder.first;
    der_vec.at(index_b_ - 1) = firstder.second;
    return ic_util::flatten_c3_vec(der_vec);
  }

  coords::float_type BondDistance::hessian_guess(coords::Representation_3D const& xyz) const {
    using ic_atom::element_period;
    using ic_atom::period;

    auto el_a = element_period(elem_a_);
    auto el_b = element_period(elem_b_);

    auto B_val{ 0.0 };
    if (el_a == period::one && el_b == period::one) {
      B_val = -0.244;
    }
    else if ((el_a == period::one && el_b == period::two) ||
      (el_a == period::two && el_b == period::one)) {
      B_val = 0.352;
    }
    else if ((el_a == period::one && el_b == period::three) ||
      (el_a == period::three && el_b == period::one)) {
      B_val = 0.660;
    }
    else if (el_a == period::two && el_b == period::two) {
      B_val = 1.085;
    }
    else if ((el_a == period::two && el_b == period::three) ||
      (el_a == period::three && el_b == period::two)) {
      B_val = 1.522;
    }
    else {
      B_val = 2.068;
    }
    return 1.734 / std::pow(val(xyz) - B_val, 3);
  }

  std::string BondDistance::info(coords::Representation_3D const & xyz) const
  {
    std::ostringstream oss;
    oss << "Bond: " << val(xyz) << " || " << index_a_ << " || " << index_b_ << " ||";
    return oss.str();
  }

  coords::float_type BondAngle::val(coords::Representation_3D const& xyz) const {
    using scon::cross;
    using scon::dot;
    using scon::len;

    auto const& a = xyz.at(index_a_ - 1u);
    auto const& b = xyz.at(index_b_ - 1u);
    auto const& c = xyz.at(index_c_ - 1u);

    auto u = a - b;
    auto v = c - b;
    auto uXv = cross(u, v);
    auto uDv = dot(u, v);
    auto l = len(uXv);
    return std::atan2(l, uDv);
  }

  std::tuple<coords::r3, coords::r3, coords::r3>
    BondAngle::der(coords::Representation_3D const& xyz) const {
    using coords::Cartesian_Point;
    using scon::cross;
    using scon::dot;
    using scon::len;

    auto const& a = xyz.at(index_a_ - 1u);
    auto const& b = xyz.at(index_b_ - 1u);
    auto const& c = xyz.at(index_c_ - 1u);

    auto u = a - b;
    auto v = c - b;
    auto lu = len(u);
    auto lv = len(v);
    u = ic_util::normalize(u);
    v = ic_util::normalize(v);
    auto cp1 = ic_util::normalize(Cartesian_Point(1.0, -1.0, 1.0));
    auto cp2 = ic_util::normalize(Cartesian_Point(-1.0, 1.0, 1.0));
    Cartesian_Point w_p(0.0, 0.0, 0.0);
    auto constexpr epsilon{ 0.01 };
    if (std::fabs(dot(u, v)) > (1. - epsilon)) {
      if (std::fabs(dot(u, cp1)) > (1. - epsilon)) {
        w_p = cross(u, cp2);
      }
      else {
        w_p = cross(u, cp1);
      }
    }
    else {
      w_p = cross(u, v);
    }
    auto w = ic_util::normalize(w_p);
    auto ad0 = cross(u, w) / lu;
    auto ad1 = cross(w, v) / lv;
    //                      a     b           c
    //                      m     o           n
    return std::make_tuple(ad0, -ad0 - ad1, ad1);
  }

  std::vector<coords::float_type>
    BondAngle::der_vec(coords::Representation_3D const& xyz) const {
    using scon::c3;

    auto firstder = BondAngle::der(xyz);
    std::vector<coords::r3> der_vec(xyz.size(), coords::r3(0.0, 0.0, 0.0));
    der_vec.at(index_a_ - 1) = std::get<0>(firstder);
    der_vec.at(index_b_ - 1) = std::get<1>(firstder);
    der_vec.at(index_c_ - 1) = std::get<2>(firstder);
    return ic_util::flatten_c3_vec(der_vec);
  }

  coords::float_type BondAngle::hessian_guess(coords::Representation_3D const& /*xyz*/) const {
    using ic_atom::element_period;
    using ic_atom::period;

    auto el_a = element_period(elem_a_);
    auto el_b = element_period(elem_b_);
    auto el_c = element_period(elem_c_);
    if (el_a == period::one ||
      el_b == period::one ||
      el_c == period::one) {
      return 0.160;
    }
    else {
      return 0.250;
    }
  }

  std::string BondAngle::info(coords::Representation_3D const & xyz) const {
    std::ostringstream oss;
    oss << "Angle: " << val(xyz) * SCON_180PI << " || " << index_a_ << " || " << index_b_ << " || " << index_c_ << " ||";
    return oss.str();
  }


  coords::float_type DihedralAngle::val(coords::Representation_3D const& xyz) const {
    using scon::cross;
    using scon::dot;
    using scon::len;


    //Eigentlich sollte das hier richtig sein:
    auto const & a = xyz.at(index_a_ - 1u);
    auto const & b = xyz.at(index_b_ - 1u);
    auto const & c = xyz.at(index_c_ - 1u);
    auto const & d = xyz.at(index_d_ - 1u);

    auto b1 = b - a;
    b1 /= len(b1);
    auto b2 = c - b;
    b2 /= len(b2);
    auto b3 = d - c;
    b3 /= len(b3);
    auto n1 = cross(b1, b2);
    auto n2 = cross(b2, b3);
    return std::atan2(dot(b1, n2), dot(n1, n2));

  }
 
  std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
    DihedralAngle::der(coords::Representation_3D const& xyz) const {
    using scon::cross;
    using scon::dot;
    using scon::len;

    auto const & a = xyz.at(index_a_ - 1u);
    auto const & b = xyz.at(index_b_ - 1u);
    auto const & c = xyz.at(index_c_ - 1u);
    auto const & d = xyz.at(index_d_ - 1u);

    auto u_p = a - b;
    auto w_p = c - b;
    auto v_p = d - c;
    auto u = ic_util::normalize(u_p);
    auto w = ic_util::normalize(w_p);
    auto v = ic_util::normalize(v_p);
    auto sin2_u = 1. - std::pow(dot(u, w), 2.);
    auto sin2_v = 1. - std::pow(dot(v, w), 2.);

    auto t1 = cross(u, w) / (len(u_p) * sin2_u);
    auto t2 = cross(v, w) / (len(v_p) * sin2_v);
    auto cuwd = cross(u, w) * dot(u, w);
    auto t3 = cuwd / (len(w_p) * sin2_u);
    //rechanged it: Lee-Ping pointed out that his numeric evaluated derivatives match with analytics without this sign ---> changed v's sign according to J. Chem. Ohys Vol 117 No. 20, 22 2002 p. 9160-9174
    auto cvwd = cross(v, w) * dot(v, w);//auto cvwd = cross(v, w) * dot(-v, w);
    auto t4 = cvwd / (len(w_p) * sin2_v);
    //exchanged +t2 with +t3 in o's derivative
    //                      a   b               c             d
    //                      m   o               p             n
    return std::make_tuple(t1, -t1 + t3 - t4, t2 - t3 + t4, -t2);
  }

  std::vector<coords::float_type>
    DihedralAngle::der_vec(coords::Representation_3D const& xyz) const {
    using scon::c3;

    auto firstder = der(xyz);
    coords::r3 temp(0.0, 0.0, 0.0);
    std::vector<coords::r3> der_vec(xyz.size(), temp);
    der_vec.at(index_a_ - 1) = std::get<0>(firstder);
    der_vec.at(index_b_ - 1) = std::get<1>(firstder);
    der_vec.at(index_c_ - 1) = std::get<2>(firstder);
    der_vec.at(index_d_ - 1) = std::get<3>(firstder);
    return ic_util::flatten_c3_vec(der_vec);
  }

  coords::float_type DihedralAngle::hessian_guess(coords::Representation_3D const & /*xyz*/) const{
    return 0.023;
  }

  std::string DihedralAngle::info(coords::Representation_3D const & xyz) const{
    std::ostringstream oss;
    oss << "Dihedral: " << val(xyz) * SCON_180PI << " || " << index_a_ << " || " << index_b_ << " || " << index_c_ << " || " << index_d_ << " ||";
    return oss.str();
  }

  coords::float_type OutOfPlane::hessian_guess(coords::Representation_3D const & xyz) const
  {
    auto const& a = xyz.at(index_a_ - 1u);
    auto const& b = xyz.at(index_b_ - 1u);
    auto const& c = xyz.at(index_c_ - 1u);
    auto const& d = xyz.at(index_d_ - 1u);
    auto r1 = b - a;
    auto r2 = b - c;
    auto r3 = b - d;
    auto r2Xr3 = cross(r2, r3);
    auto rd = dot(r1, r2Xr3);
    auto t2 = rd / (len(r1) * len(r2) * len(r3));
    auto dd = 1. - t2;
    auto d_pow = std::pow(dd, 4);
    return 0.045 * d_pow;
  }

  std::string OutOfPlane::info(coords::Representation_3D const & xyz) const
  {
    std::ostringstream oss;
    oss << "Out of plane: " << val(xyz) * SCON_180PI << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "||" << index_d_ << "\n";
    return oss.str();
  }

  std::vector<coords::float_type>
    TranslationX::der_vec(coords::Representation_3D const& xyz) const {
    using cp = coords::Cartesian_Point;

    return ic_util::flatten_c3_vec(der(xyz.size(), [](auto const & s) {
      return cp(1. / static_cast<coords::float_type>(s), 0., 0.);
    }));
  }

  std::string TranslationX::info(coords::Representation_3D const & xyz) const
  {
    std::ostringstream oss;
    oss << "Trans X: " << val(xyz);
    return oss.str();
  }

  std::vector<coords::float_type>
    TranslationY::der_vec(coords::Representation_3D const& xyz) const {

    using cp = coords::Cartesian_Point;

    return ic_util::flatten_c3_vec(der(xyz.size(), [](auto const & s) {
      return cp(0., 1. / static_cast<coords::float_type>(s), 0.);
    }));
  }

  std::string TranslationY::info(coords::Representation_3D const & xyz) const
  {
    std::ostringstream oss;
    oss << "Trans Y: " << val(xyz);
    return oss.str();
  }

  std::vector<coords::float_type>
    TranslationZ::der_vec(coords::Representation_3D const& xyz) const {
    using cp = coords::Cartesian_Point;

    return ic_util::flatten_c3_vec(der(xyz.size(), [](auto const & s) {
      return cp(0., 0., 1. / static_cast<coords::float_type>(s));
    }));
  }

  std::string TranslationZ::info(coords::Representation_3D const & xyz) const
  {
    std::ostringstream oss;
    oss << "Trans Z: " << val(xyz);
    return oss.str();
  }
}