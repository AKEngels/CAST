#ifndef linesearch_header

#define linesearch_header

/*

 * The implementation has been adapted from:

 *      C library of Limited memory BFGS (L-BFGS).
 *
 * Copyright (c) 1990, Jorge Nocedal
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <cstddef> 
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <tuple>
#include <type_traits>

#include "function_trait.h"
#include "representation.h"

extern std::vector<std::vector<double>> locpx, locx, locd;

template<class T>
inline std::vector<double> change(T const &x)
{
  std::vector<double> r;
  r.reserve(x.size() * 3u);
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    r.push_back(x[i].x());
    r.push_back(x[i].y());
    r.push_back(x[i].z());
  }
  return r;
}

namespace optimization
{

  namespace local
  {

    enum status
    {
      ERR_INVALID_STEPSIZE = -100,
      ERR_GRAD_INCREASE,
      ERR_CALLBACK_STOP,
      ERR_ROUNDING,
      ERR_MAXSTEP,
      ERR_MINSTEP,
      ERR_WIDTHTOOSMALL,
      ERR_MAX_LS_ITERATIONS,
      ERR_MAX_ITERATIONS = -1,
      SUCCESS = 0,
      UNDEFINED
    };

    struct empty_void_functor
    {
      template<class ...T> void operator() (T&&...) { }
    };

    namespace linesearch
    {

      template<class CallbackT>
      struct is_valid_callback
      {
        typedef function_trait_detail::argument_type<CallbackT, 0U> argument_0;
        typedef function_trait_detail::argument_type<CallbackT, 1U> argument_1;
        static const bool value = function_trait_detail::function_traits<CallbackT>::arity == std::size_t(4u) // 4 arguments total
          //&& detail::has_dot<argument_0, argument_0>::value // dot of argument 0 and 0
          //&& detail::has_dot<argument_1, argument_1>::value // dot of argument 1 and 1
          //&& detail::has_dot<argument_0, argument_1>::value // dot of argument 0 and 1
          //&& detail::has_dot<argument_1, argument_0>::value // dot of argument 1 and 0
          && std::is_same<function_trait_detail::decayed_argument_type<CallbackT, 2U>, std::size_t>::value // third argument is size_t 
          && std::is_same<function_trait_detail::argument_type<CallbackT, 3u>, bool&>::value; // last argument bool reference
      };

      namespace _detail_
      {
        template < class T, class M > M get_member_type(M T:: *);
      }

      template<class T> using callback_type = decltype(_detail_::get_member_type(&T::callback));

      template<class L>
      struct is_valid_linesearch
      {
        static bool const value = is_valid_callback<callback_type<L>>::value;
      };


      template<class CallbackT>
      struct none
      {

        using callback_type = CallbackT;
        using float_type = function_trait_detail::return_type < callback_type >;
        using rep_type = function_trait_detail::decayed_argument_type < callback_type, 0U >;
        using grad_type = function_trait_detail::decayed_argument_type < callback_type, 1U >;

        using point_type = Point < rep_type, grad_type, float_type > ;

        callback_type callback;
        float_type max_step, min_step;
        status state;

        none(callback_type callback_object,
          float_type const max_stepsize = 1.e20,
          float_type const min_stepsize = 1.e-20)
          : callback(callback_object), max_step(max_stepsize), 
          min_step(min_stepsize), state(status::UNDEFINED)
        { }

        status operator() (grad_type const &d, float_type & step, 
          point_type & p, rep_type const &x, std::size_t const iter)
        {
          using std::max;
          using std::min;
          step = max(min_step, min(max_step, step));
          if (step > max_step) return status::ERR_MAXSTEP;
          if (step < min_step) return status::ERR_MINSTEP;
          p.x = x + d*step;
          bool go_on(true);
          p.f = callback(p.x, p.g, iter, go_on);
          if (!go_on) return status::ERR_CALLBACK_STOP;
          return status::SUCCESS;
        }

      };


      template<class F>
      bool signdiff(F const x, F const y)
      {
        using std::abs;
        return ((x*(y / abs(y))) < F(0));
      }

      namespace minimizers
      {
        /**
         * Find a minimizer of an interpolated cubic function.
         *  @return         The minimizer of the interpolated cubic.
         *  @param  u       The value of one point, u.
         *  @param  fu      The value of f(u).
         *  @param  du      The value of f'(u).
         *  @param  v       The value of another point, v.
         *  @param  fv      The value of f(v).
         *  @param  du      The value of f'(v).
         */
        template<class F>
        inline F cubic(F const u, F const fu, F const du, F const v, F const fv, F const dv)
        {
          using std::abs;
          using std::sqrt;
          using std::max;
          F const d(v - u), theta((fu - fv)*F(3) / d + du + dv);
          F const s(max(max(abs(theta), abs(du)), abs(dv)));
          F const a(theta / s);
          F gamma = s * sqrt(a * a - (du / s) * (dv / s));
          if (v < u) gamma = -gamma;
          F const r((gamma - du + theta) / (gamma - du + gamma + dv));
          return u + r*d;
        }
        /**
         * Find a minimizer of an interpolated cubic function.
         *  @return         The minimizer of the interpolated cubic.
         *  @param  u       The value of one point, u.
         *  @param  fu      The value of f(u).
         *  @param  du      The value of f'(u).
         *  @param  v       The value of another point, v.
         *  @param  fv      The value of f(v).
         *  @param  du      The value of f'(v).
         *  @param  xmin    The maximum value.
         *  @param  xmin    The minimum value.
         */
        template<class F>
        inline F cubic(F const u, F const fu, F const du,
          F const v, F const fv, F const dv,
          F const min_value, F const max_value)
        {
          using std::abs;
          using std::sqrt;
          using std::max;
          F const d(v - u);
          F const theta((fu - fv) * F(3) / d + du + dv);
          F const s(max(max(abs(theta), abs(du)), abs(dv)));
          F const a(theta / s);
          F gamma = s * sqrt(max(F(0), a*a - (du / s) * (dv / s)));
          if (u < v) gamma = -gamma;
          F const r((gamma - dv + theta) / (gamma - dv + gamma + du));
          if (r < F(0) && abs(gamma) > F(0))
          {
            return v - r*d;
          }
          else if (a < F(0))
          {
            return max_value;
          }
          else
          {
            return min_value;
          }
        }
        /**
         * Find a minimizer of an interpolated quadratic function.
         *  @return         The minimizer of the interpolated quadratic.
         *  @param  u       The value of one point, u.
         *  @param  fu      The value of f(u).
         *  @param  du      The value of f'(u).
         *  @param  v       The value of another point, v.
         *  @param  fv      The value of f(v).
         */
        template<class F>
        inline F quad(F const u, F const fu, F const du, F v, F const fv)
        {
          v -= u;
          return u + du / ((fu - fv) / v + du) / F(2) * v;
        }
        /**
         * Find a minimizer of an interpolated quadratic function.
         *  @return         The minimizer of the interpolated quadratic.
         *  @param  u       The value of one point, u.
         *  @param  du      The value of f'(u).
         *  @param  v       The value of another point, v.
         *  @param  dv      The value of f'(v).
         */
        template<class F>
        inline F quad(F const u, F const du, F const v, F const dv)
        {
          return v + dv / (dv - du) * (u - v);
        }

      }

      /*

        float_type = floating point type
        rep_type = representation class type
        grad_type = gradient class type
        callback_type callback for gradients from rep
        -> float_type operator() (rep_type const&, grad_type &, std::size_t const, bool&);

      */

      template<class float_type, class rep_type, class grad_type = rep_type>
      using callback_function_type = float_type(*)(rep_type const &, grad_type&, std::size_t const, bool&);

      template<class CallbackT>
      struct more_thuente
      {

        using callback_type = CallbackT;
        using float_type = function_trait_detail::return_type < callback_type >;
        using rep_type = function_trait_detail::decayed_argument_type < callback_type, 0U >;
        using grad_type = function_trait_detail::decayed_argument_type < callback_type, 1U >;
        using point_type = Point < rep_type, grad_type, float_type >;

        using F = float_type;

        //callback instances
        callback_type callback;

        struct configuration
        {
          float_type max_step, min_step, ftol, gtol, xtol;
          std::size_t max_iterations;
          bool ignore_callback_stop;
          configuration() : max_step(F(1.e20)), 
            min_step(F(1.e-20)), ftol(F(1.e-4)), 
            gtol(F(0.9)), xtol(F(1.e-16)), 
            max_iterations(50U), ignore_callback_stop(true){ }
        } config;

        status state;

        more_thuente(callback_type callback_object)
          : callback(std::move(callback_object)),
          config(), state(status::SUCCESS)
        { }

        status operator()(grad_type const & d, float_type & step, 
          point_type & p, rep_type const & xp, std::size_t const iter)
        { return line(d, step, p.x, p.g, p.f, xp, iter); }

      private:

        status line(grad_type const & d, float_type & step, 
          rep_type & p_x, grad_type & g, float_type & f,
          rep_type const & xp, std::size_t const iter)
        {

          using std::max;
          using std::min;
          using std::abs;

          std::size_t const N(config.max_iterations);
          bool brackt(false);

          float_type const dginit(dot(d, g)), finit(f), dgtest(config.ftol*dginit);
          float_type const dg_limit(min(config.ftol, config.gtol)*dginit);
          float_type dg = F();
          float_type stx = F(), fx(finit), dgx(dginit);
          float_type sty = F(), fy(finit), dgy(dginit);
          float_type stmin = F(), stmax = F(), ftest1 = F();
          float_type width(config.max_step - config.min_step);
          float_type prev_width((F)2 * width);
          
          int uinfo(0);
          bool stage1(true);

          state = status::SUCCESS;

          if (step < 0.)
          {
            state = status::ERR_INVALID_STEPSIZE;
            return state;
          }

          if (dginit > 0.)
          {
            state = status::ERR_GRAD_INCREASE;
            return state;
          }

          for (std::size_t i(0U); i < N; ++i)
          {

            if (brackt)
            {
              stmin = min(stx, sty);
              stmax = max(stx, sty);
            }
            else
            {
              stmin = stx;
              stmax = step + F(4) * (step - stx);
            }

            // Let step be in the range 
            // from config.min_step to config.max_step
            step = max(config.min_step, min(config.max_step, step));
            /*
              If an unusual termination is to occur then let
              stp be the lowest point obtained so far.
            */
            if ((brackt && ((step <= stmin || stmax <= step) || uinfo != 0)) ||
              (brackt && (stmax - stmin <= config.xtol * stmax)))
            {
              step = stx;
            }
            // Advance  p_x along d about amount step
            p_x = xp + d*step;
            // Grab new function value and gradients
            bool go_on(true);
            f = callback(p_x, g, iter, go_on);
            if (!go_on && !config.ignore_callback_stop)
            {
              return status::ERR_CALLBACK_STOP;
            }
            // d \cdot g
            dg = dot(d, g);
            //float_type const dg(dot(d, g));
            // f test value
            ftest1 = finit + step*dgtest;

            //printf("LSA:\n");
            //std::cout << "step\n" << step << "\n";
            //std::cout << "xp\n" << xp << "\n\n";
            //std::cout << "d\n" << d << "\n\n";
            //std::cout << "x\n" << p_x << "\n\n";
            //
            //printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", dg, stx, fx, dgx, sty, fy, dgy);
            //printf("%.6f %.6f %.6f %.6f\n", finit, ftest1, dginit, dgtest);
            //printf("%.6f %.6f %.6f %.6f\n", step, stmin, stmax, f);

            /* Test for errors and convergence. */
            if (brackt && ((step <= stmin || stmax <= step) || uinfo != 0))
            { /* Rounding errors prevent further progress. */
              state = status::ERR_ROUNDING;
              //std::cout << stmin << " !< " << step << " !< " << stmax << ", " << uinfo << "\n";
              break;
            }
            if (abs(step - config.max_step) < F(1.e-6)
              && f <= ftest1 && dg <= dgtest)
            { /* The step is the maximum value. */
              state = status::ERR_MAXSTEP;
              break;
            }
            if (abs(step - config.min_step) < F(1.e-6)
              && (ftest1 < f || dgtest <= dg))
            { /* The step is the minimum value. */
              //std::cout << "Stepl = " << step << "\n";
              state = status::ERR_MINSTEP;
              break;
            }
            if (brackt && (stmax - stmin) <= config.xtol*stmax)
            { /* Relative width of the interval 
                 of uncertainty is at most xtol. */
              state = status::ERR_WIDTHTOOSMALL;
              break;
            }
            if (f <= ftest1 && std::abs(dg) <= config.gtol * (-dginit))
            { /* The sufficient decrease condition 
                 and the directional derivative condition hold. */
              state = status::SUCCESS;
              break;
            }

            if (stage1 && f <= ftest1 && dg_limit <= dg)
            {
              stage1 = 0;
            }

            /*
            A modified function is used to predict the step only if
            we have not obtained a step for which the modified
            function has a nonpositive function value and nonnegative
            derivative, and if a lower function value has been
            obtained but the decrease is not sufficient.
            */
            if (stage1 && ftest1 < f && f <= fx)
            {
              /* Define the modified function and derivative values. */
              float_type fm = f - step*dgtest;
              float_type fxm = fx - stx * dgtest;
              float_type fym = fy - sty * dgtest;
              float_type dgm = dg - dgtest;
              float_type dgxm = dgx - dgtest;
              float_type dgym = dgy - dgtest;
              /*
              Call update_trial_interval() to update the interval of
              uncertainty and to compute the new step.
              */
              uinfo = more_thuente::update_trial_interval(stx, fxm, dgxm, sty, fym, dgym, step, fm, dgm, stmin, stmax, brackt);
              /* Reset the function and gradient values for f. */
              fx = fxm + stx * dgtest;
              fy = fym + sty * dgtest;
              dgx = dgxm + dgtest;
              dgy = dgym + dgtest;
            }
            else
            {
              /*
              Call update_trial_interval() to update the interval of
              uncertainty and to compute the new step.
              */
              uinfo = more_thuente::update_trial_interval(stx, fx, dgx, sty, fy, dgy, step, f, dg, stmin, stmax, brackt);
            }
            /*
            Force a sufficient decrease in the interval of uncertainty.
            */
            if (brackt)
            {
              if (F(0.66)*prev_width <= std::abs(sty - stx))
              {
                step = stx + F(0.5)*(sty - stx);
              }
              prev_width = width;
              width = abs(sty - stx);
            }
            if (i == N - 1) state = status::ERR_MAX_LS_ITERATIONS;
          }
          return state;
        }

        static int update_trial_interval(
          float_type &x, float_type &fx, float_type &dx,
          float_type &y, float_type &fy, float_type &dy,
          float_type &t, float_type const ft, float_type const dt,
          float_type const tmin, float_type const tmax, bool &brackt)
        {
          bool bound(false);
          bool const dsign(signdiff(dt, dx));
          float_type mc = F(); /* minimizer of an interpolated cubic. */
          float_type mq = F(); /* minimizer of an interpolated quadratic. */
          float_type newt = F();   /* new trial value. */
          using std::min;
          using std::max;
          if (brackt)
          {
            // t not in [x,y]
            if (t <= (min(x, y)) || (max(x, y)) <= t)
              return -3;
            // function increases
            if (F(0) <= dx*(t - x))
              return -2;
            // max < min?
            if (tmax < tmin)
              return -1;
          }
          if (fx < ft)
          {
            brackt = true;
            bound = true;
            mc = minimizers::cubic(x, fx, dx, t, ft, dt);
            mq = minimizers::quad(x, fx, dx, t, ft);
            if (abs(mc - x) < abs(mq - x)) newt = mc;
            else newt = mc + (mq - mc) / F(2);
          }
          else if (dsign)
          {
            brackt = true;
            bound = false;
            mc = minimizers::cubic(x, fx, dx, t, ft, dt);
            mq = minimizers::quad(x, dx, t, dt);
            newt = (abs(mc - t) > abs(mq - t)) ? mc : mq;
          }
          else if (abs(dt) < abs(dx))
          {
            bound = true;
            mc = minimizers::cubic(x, fx, dx, t, ft, dt, tmin, tmax);
            mq = minimizers::quad(x, dx, t, dt);
            if (brackt)
              newt = (abs(t - mc) < abs(t - mq)) ? mc : mq;
            else
              newt = (abs(t - mc) > abs(t - mq)) ? mc : mq;
          }
          else
          {
            bound = false;
            if (brackt) newt = minimizers::cubic(t, ft, dt, y, fy, dy);
            else if (x < t) newt = tmax;
            else newt = tmin;
          }

          if (fx < ft)
          {
            y = t;
            fy = ft;
            dy = dt;
          }
          else
          {
            if (dsign)
            {
              y = x;
              fy = fx;
              dy = dx;
            }
            x = t;
            fx = ft;
            dx = dt;
          }
          // clip newt in [tmin, tmax]
          newt = min(tmax, max(tmin, newt));
          // redefine newt if close to upper bound
          if (brackt && bound)
          {
            float_type const tmp_t = x + F(0.66)*(y - x);
            if (x < y)
            {
              if (tmp_t < newt)
              {
                newt = tmp_t;
              }
            }
            else if (newt < tmp_t)
            {
              newt = tmp_t;
            }
          }
          // Return the new trial value.
          t = newt;
          return 0;
        }
      };

    }

    template<class T, class...Args>
    linesearch::none<typename std::remove_reference<T>::type> 
      make_empty_ls(T && callback_object, Args && ... args)
    {
      return linesearch::none< typename std::remove_reference<T>::type >
        (std::forward<T>(callback_object), std::forward<Args>(args)...);
    }

    template<class T>
    linesearch::more_thuente<typename std::remove_reference<T>::type> 
      make_more_thuente(T && callback_object)
    {
      return linesearch::more_thuente<typename std::remove_reference<T>::type>
        (std::forward<T>(callback_object));
    }



  }

}

#endif
