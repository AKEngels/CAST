#ifndef optimization_dimer_header

#define optimization_dimer_header

#include <fstream>
#include <algorithm>
#include <cstddef>
#include <utility>

#include "Scon/function_trait.h"
#include "representation.h"
#include "optimization.h"

#include "Scon/scon_traits.h"
#include "Scon/scon_utility.h"
#include "Scon/scon_vect.h"

#include "cg.h"

/*

Lit.
[1] Kaestner, Sherwood, J. Chem. Phys., 128 (2003), doi 10.1063/1.2815812
[2] Henkelmann, Jonsson, J. Chem. Phys., 111 (1999), doi 10.1063/1.480097

*/


namespace optimization
{

  namespace dimer
  {
    template<class CallbackT>
    class dimer
    {

    public:

      using callback_type = CallbackT;
      using float_type = function_trait_detail::return_type<callback_type>;
      using rep_type = function_trait_detail::decayed_argument_type<callback_type, 0U>;
      using grad_type = function_trait_detail::decayed_argument_type<callback_type, 1U>;
      using point_type = Callback_Point<callback_type>;
      using minimum_type = Callback_Minimum < callback_type >;
      using ortho_vector = typename minimum_type::directions_type;

      callback_type force;
      ortho_vector passed_directions;

    private:

      point_type s0, s1;
      rep_type tau, init_x, last_tau;
      grad_type f_rot, f_trans;
      float_type C, D, rot_convergence;
      std::size_t max_rot_iterations;
      bool integrity, req_ortho;

    public:

      bool log;

      template<class CallArgs>
      dimer(CallArgs&& callback_arguments, minimum_type const & start,
        float_type const distance = 0.1, std::size_t max_rot_iter = 20U,
        float_type rot_conv = 0.01) 
        : force(std::forward<CallArgs>(callback_arguments)), 
        passed_directions(start.directions), s0(start),
        s1(start), tau(scon::randomized<rep_type>(start.x)), init_x(start.x), 
        last_tau(start.x), f_rot(), f_trans(), C(0), D(distance), 
        rot_convergence(rot_conv), max_rot_iterations(max_rot_iter), 
        integrity(true), req_ortho(false), log(false)
      {
        using std::sqrt;
        tau /= sqrt(dot(tau, tau));
        s1.x = s0.x + tau*D;
        update(true, true);
        check_ortho_directions();
        if (orthogonalization())
        {
          orthogonalize(tau);
          tau /= sqrt(dot(tau, tau));
          s1.x = s0.x + tau*D;
          update(false, true);
        }

        if (rotate_to_min_curve_cg() > 0U)
        {
          check_inversion();
          last_tau = tau;
          push_on();
        }

      }

      float_type rot_conv() const { return rot_convergence; }
      void rot_conv(float_type const val) { rot_convergence = val; }

      float_type operator() (rep_type const &x, grad_type &g, 
        std::size_t const, bool &go_on)
      {
        //std::cout << "Moving dimer to " << x << " with tau = " << tau << "\n";
        move_to(x);
        //std::cout << "New g0 " << s0.g << "\n";
        check_ortho_directions();
        //std::cout << "Orthogonalized tau = " << tau << "\n";
        std::size_t r = rotate_to_min_curve_cg();
        auto const ttd = dot(last_tau, tau);
        if ((1. - ttd) > 1.9)
        {
          tau = -tau;
          s1.x = s0.x + tau*D;
          update();
        }
        //std::cout << " C = " << C << " << " << tau << "\n";
        if (r == 0U) go_on = false;
        g = f_trans;
        //std::cout << "f_trans = " << f_trans << "\n";
        return -s0.f;
      }

      point_type const & middle() const { return s0; }
      point_type & middle() { return s0; }

      rep_type const & direction() const { return tau; }
      float_type curve() const { return C; }
      float_type distance() const { return D; }

      void push_over()
      {
        using std::sqrt;
        std::size_t i = 0;
        while (C < 0 && ++i < 10)
        {
          move_to(s1.x);
          //std::cout << " Pushing " << C << ", t: " << tau << "\n";
          D = D*1.2;
        }
        move_to(s1.x);
        //std::cout << " Pushing " << C << "\n";
      }

      void flip()
      {
        tau = -tau;
        s1.x = s0.x + tau*D;
        update();
      }

    private:

      void check_inversion()
      {
        auto s2x = s0.x - tau*D;
        grad_type s2g;
        bool b;
        auto s2f = force(s2x, s2g, 0u, b);
        if (s2f < s1.f)
        {
          using std::swap;
          swap(s2x, s1.x);
          swap(s2f, s1.f);
          swap(s2g, s1.g);
          tau = -tau;
          update(false, false);
        }
      }

      void check_ortho_directions()
      {
        std::vector<rep_type> ortho_required_directions;
        bool proceed;
        for (auto & v : passed_directions)
        {
          grad_type nf1;
          auto const nx1 = s0.x + v*D;
          force(nx1, nf1, 0, proceed);
          auto const f01(nf1 - s0.g);
          auto const f01_dot_v(dot(v, f01));
          auto const nC = f01_dot_v / D;
          if (nC < C)
          { // if curvature of passed direction lower, keep it
            ortho_required_directions.emplace_back(std::move(v));
          }
        }
        passed_directions.swap(ortho_required_directions);
      }

      void orthogonalize(rep_type & x) const
      {
        for (auto & v : passed_directions)
        {
          x -= v*dot(x, v);
        }
      }

      bool orthogonalization() const { return !passed_directions.empty(); }
      
      void update(bool const update_x0 = false, 
        bool const update_x1 = true)
      {
        using std::sqrt;
        bool proceed(true);
        if (update_x1)
        {
          s1.f = force(s1.x, s1.g, 1, proceed);
          orthogonalize(s1.g);
        }
        if (proceed && update_x0)
        {
          s0.f = force(s0.x, s0.g, 0, proceed);
          orthogonalize(s0.g);
        }
        if (!proceed) integrity = false;
        else
        {
          //if (log) s0.log();
          //if (log) s1.log();
          auto const f01(s1.g - s0.g);
          auto const minus_two_f01(-(f01 * 2));
          auto const f01_dot_tau(dot(tau, f01));
          f_rot = minus_two_f01 
            + tau * (f01_dot_tau * 2);  // ([1] Eq. (3))
          C = f01_dot_tau / D;  // ([1] Eq. (4))
          //f_trans = s0.g - tau * (dot(s0.g, tau) * 2); // ([1] Eq. (2))
          //f_trans = tau * (-dot(s0.g, tau));
          if (C > 0.0)
          {
            f_trans = tau * (-dot(s0.g, tau));
          }
          else
          {
            f_trans = s0.g - tau * (dot(s0.g, tau) * 2); // ([1] Eq. (2))
            f_trans = tau*dot(f_trans, tau);
          }

        }
      }

      void rotate(grad_type const &O, 
        float_type const phi)
      {
        using std::cos;
        using std::sin;
        using std::sqrt;
        auto const cosine = cos(phi)*D, sine = sin(phi)*D;
        s1.x = s0.x + (tau*cosine + O*sine); // ([2] Eq. (4))
        tau = s1.x - s0.x;
        if (orthogonalization())
        {
          orthogonalize(tau);
          tau = tau / sqrt(dot(tau, tau));
          s1.x = s0.x + tau*D;
        }
        else
        {
          tau = tau / sqrt(dot(tau, tau));
        }
        update();
      }

      void push_on()
      {
        using std::swap;
        std::swap(s0, s1);
        s1.x = s0.x + tau*D;
        update();
      }

      void move_to(rep_type const &new_x0)
      {
        using std::sqrt;
        s0.x = new_x0;
        s1.x = s0.x + tau*D;
        update(true, true);
      }

      std::size_t rotate_to_min_curve_cg()
      {
        using std::abs;
        using std::atan;
        using std::cos;
        using std::sin;
        using std::sqrt;
        //std::cout << "Init tau: " << tau << "\n";
        grad_type G_old, F_old, O_old;
        std::size_t i(0U);

				// this line was added on 25.7.19 in order to avoid errors of type "Unqually sized vectors given for dot product." (issue #3)
				// which is thrown in line 312 because G = f_rot (line 307) is empty if proceed = false in line 215 
				if (integrity == false) return 0U;    

        for ( ; i < max_rot_iterations; )
        {
          grad_type G;
          if (i > 0U)
          {
            // ([2] Eq. (17))
            auto const gamma = dot((f_rot - F_old), f_rot) 
              / dot(f_rot, f_rot);
            // ([2] Eq. (16))
            G = f_rot + O_old * gamma*sqrt(dot(G_old, G_old));
          }
          else
          {
            G = f_rot;
            //std::cout << "Init G: " << G << "\n";
          }

          G_old = G;
          G -= tau*dot(G, tau);
          G /= sqrt(dot(G, G));
          F_old = f_rot;
          grad_type f01(s1.g - s0.g), twof01(f01 * 2);
          // ([1] Eq. (6))
          float_type const dC_dphi(dot(twof01, G) / D);
          // ([1] Eq. (5))
          float_type const phi1(-atan(dC_dphi / (abs(C) * 2)) / 2);
          // out
          //std::cout << "phi1 " << phi1 << "\n";
          // tolerance angle
          if (abs(phi1) < rot_convergence)
          {
            //rotate(G, phi1);
            return i + 1U;
          }
          dimer dimer_tmp(*this);
          dimer_tmp.rotate(G, phi1);
          // new tau
          //std::cout << "New tau: " << dimer_tmp.tau << "\n";
          if (dimer_tmp.integrity == false) return 0U;
          // ([1] Eq. (8))
          float_type const b1(dC_dphi / 2);
          // ([1] Eq. (9))
          float_type const a1((C - dimer_tmp.C
            + b1 * sin(phi1 * 2)) / (1 - cos(phi1 * 2)));
          // ([1] Eq. (10))
          float_type const a0((C - a1) * 2);
          float_type phi(atan(b1 / a1) / 2);
          //out
          //std::cout << "phi " << phi << "\n";
          // ([1] Eq. (7))
          float_type const nC(a0 / 2 + a1 * cos(phi * 2)
            + b1 * sin(phi * 2));
          if (nC > C)
          {
            phi += float_type(3.14159265358979323846 / 2);
          }
          if (abs(phi) < rot_convergence)
          {
            //rotate(G, phi);
            return i + 1U;
          }
          else
          {
            rotate(G, phi);
            if (integrity == false) return 0U;
            O_old = G - tau * dot(G, tau);
            O_old /= sqrt(dot(O_old, O_old));
          }
          ++i;
        }
        return i;
      }

    };

    template<class CallbackT>
    struct linesearch
    {

      using callback_type = CallbackT;
      using float_type = function_trait_detail::return_type<callback_type>;
      using rep_type = function_trait_detail::decayed_argument_type<callback_type, 0U>;
      using grad_type = function_trait_detail::decayed_argument_type<callback_type, 1U>;

      using point_type = Point < rep_type, grad_type, float_type >;

      callback_type callback;

      local::status state;

      grad_type last_m;
      float_type max_step, min_step;
      

      linesearch(callback_type callback_object,
        float_type const max_stepsize = 1.0,
        float_type const min_stepsize = 0.02)
        : callback(callback_object), max_step(max_stepsize), 
        min_step(min_stepsize)
        { }

      local::status operator()(grad_type const & d,
        float_type & step, point_type & p,
        rep_type const & xp, std::size_t const iter)
      {

        using std::max;
        using std::min;
        using std::sqrt;
        grad_type m = d*step;
        m = callback.direction() * dot(callback.direction(), m);
        auto const mlen = sqrt(dot(m, m));
        auto steplen = max(min_step, min(mlen, max_step));
        m /= mlen;
        if (iter > 0)
        {
          auto const mmd = dot(m, last_m);
          if ((1. - mmd) > 1.9)
          {
            //callback.flip();
            return local::status::UNDEFINED;
          }
        }
        last_m = m;
        m *= steplen;
        p.x = xp + m;
        bool go_on(true);
        p.f = callback(p.x, p.g, iter, go_on);
        if (!go_on)
        {
          return local::status::ERR_CALLBACK_STOP;
        }
        return local::status::SUCCESS;

      }

    };

  }

  template<class CallbackT>
  struct CG_DimerMover
  {

    using callback_type = CallbackT;
    using dimer_type = dimer::dimer < callback_type >;
    using linesearch_type = dimer::linesearch < dimer_type >;
    using Optimizer = local::conjugate_gradient < linesearch_type >;
    using minimum_type = typename dimer_type::minimum_type;
    using point_type = typename dimer_type::point_type;
    using float_type = typename dimer_type::float_type;

    float_type distance;
    std::size_t iter;
    local::status state;

    CG_DimerMover(float_type const dimer_distance = 0.1, 
      std::size_t const max_iterations = 100U)
      : distance(dimer_distance), iter(max_iterations), state()
    { }

    callback_type operator() (minimum_type &p, callback_type c)
    {
      dimer_type d(std::move(c), p, distance);
      // push current dimer direction into minimum p
      p.directions.push_back(d.direction());
      //linesearch_type l(d);
      // Create optimizer
      Optimizer o(d);
      // Set config
      o.config.reset_offset = 5u;
      o.config.epsilon = 1.0e-3;
      o.config.max_iterations = iter;
      // Optimize dimer midpoint
      o(d.middle());
      // Dimer is the callback of the linesearch
      // translate dimer across maximum
      o.ls.callback.push_over();
      // Move dimer middle point into p
      p.point_type::operator=(std::move(o.ls.callback.middle()));
      // return callback which is "force" member 
      // of dimer which is callback member 
      // of linesearch which is ls member of optimizer
      return o.ls.callback.force;
    }

  };

}

#endif

