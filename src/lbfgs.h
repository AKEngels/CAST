#ifndef lbfgs_header

#define lbfgs_header

/*


 * The implementation has been adapted from:

 *
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
 *


*/

#include <utility>
#include <vector>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include "ls.h"


namespace optimization
{
  namespace local
  {

    /*

    Requirements regarding representation type V and gradient type G:

    - function dot(V const &, V const &) returning the scalar product V^T times V
    - function dot(V const &, G const &) returning the scalar product V^T times G
    - Member function V::invert() Inverting every
    - Defined operator- to obtain V1-V2 andf G1-G2
    - Element type E of V needs to have operator* defined for multiplication with scalar values
    - Element type E of V needs to have operator+= (E const&) defined for adding one E to another

    */


    template<class LineSearchT, class LoggerT = empty_void_functor>
    class lbfgs
    {

    public:

      using linesearch_type = typename std::enable_if<
        linesearch::is_valid_linesearch<LineSearchT>::value, LineSearchT
      >::type;

      using float_type = typename linesearch_type::float_type;
      using rep_type = typename linesearch_type::rep_type;
      using grad_type = typename linesearch_type::grad_type;
      using callback_type = typename linesearch_type::callback_type;

      using point_type = Point < rep_type, grad_type, float_type >;
      using state_type = State < rep_type, grad_type >;

      using logger_type = LoggerT;

      // LineSearch instance
      linesearch_type ls;
      logger_type log;

      struct configuration
      {
        // Number of Hessian Corrections
        std::size_t m;
        // Number of History steps for delta convergence test
        std::size_t k;
        // maximum number of iterations
        std::size_t max_iterations;
        // Delta for convergence test.
        float_type delta;
        // Convergence epsilon
        float_type epsilon;
        // Constructor with initializer list for default values
        configuration() :
          m(6u), k(), max_iterations(500u),
          delta(), epsilon(F(1.e-4))
        { }
      } config;

      // Storage struct for point and norms
      struct xgf
      {

        point_type p;
        float_type xnorm, gnorm;

        xgf() : p(), xnorm(), gnorm() { }
        xgf(point_type const& v)
          : p(v), xnorm(), gnorm() { }

        void update(callback_type& _cbo,
          std::size_t _iteration, bool& _proceed)
        {
          p.f = _cbo(p.x, p.g, _iteration, _proceed);
          norms();
        }

        void norms()
        {
          using std::sqrt;
          xnorm = sqrt(dot(p.x, p.x));
          gnorm = sqrt(dot(p.g, p.g));
        }

      } xg, xg_p;

    private:

      using F = float_type;

      struct hesscorrection
      {
        state_type d;
        float_type dXdG, dGdG, alpha;
        hesscorrection() :
          d(), dXdG(), dGdG(), alpha()
        { }
        void assign(xgf const& current, xgf const& stored)
        {
          d.x = current.p.x - stored.p.x;
          d.g = current.p.g - stored.p.g;
          dXdG = dot(d.x, d.g);
          dGdG = dot(d.g, d.g);
        }
      };

      std::vector<hesscorrection> correctH;

      grad_type d;
      float_type dXdG, dGdG, step;
      std::size_t iteration, last_hess;
      status rstate;
      bool go_on;

      status optimize()
      {
        using std::sqrt;
        if (convergence())
        {
          return rstate = status::SUCCESS;
        }
        last_hess = 0U;
        step = float_type(1) / sqrt(dot(d, d));  // d = forces 
        log(xg.p);
        for (std::size_t i(1U); i <= config.max_iterations; ++i)
        {
          store_xg();
          auto const lsr = ls(d, step, xg.p, xg_p.p.x, i);  // goes to ls.h, line 251
          if (lsr != status::SUCCESS)
          {
            restore_xg();
            return lsr;
          }
          log(xg.p);
          xg.norms();
          ++iteration;
          if (convergence())
          {
            return rstate = status::SUCCESS;
          }
          assignHessianCorrection();
          updateHessian(i - 1);
          step = 1.0;
        }
        return rstate = status::ERR_MAX_ITERATIONS;
      }

      void init_()
      {
        rstate = status::UNDEFINED;
        correctH.assign(config.m,
          hesscorrection());
      }

      bool convergence()
      {
        using std::max;
        return ((xg.gnorm / max(xg.xnorm, F(1)))
          < config.epsilon);
      }

      void assignHessianCorrection()
      {
        if (last_hess < correctH.size())
        {
          correctH[last_hess].assign(xg, xg_p);
          dXdG = correctH[last_hess].dXdG;
          dGdG = correctH[last_hess].dGdG;
        }
      }

      void updateHessian(std::size_t const iter)
      {
        std::size_t const N = correctH.size(),
          M = std::min(N, iter);
        last_hess = (last_hess + 1U) % N;
        d = -xg.p.g;
        std::size_t j(last_hess);
        for (std::size_t i(0U); i < M; ++i)
        {
          j = (j + N - 1U) % N;
          correctH[j].alpha = dot(correctH[j].d.x, d) / correctH[j].dXdG;
          d += correctH[j].d.g * (-correctH[j].alpha);
        }
        d *= dXdG / dGdG;
        for (std::size_t i(0U); i < M; ++i)
        {
          float_type const beta(correctH[j].alpha
            - dot(correctH[j].d.g, d) / correctH[j].dXdG);
          d += correctH[j].d.x * beta;
          j = (j + 1U) % N;
        }
      }

      void store_xg() { xg_p = xg; }
      void restore_xg() { xg = xg_p; }

    public:

      lbfgs(linesearch_type linesearch_object,
        logger_type logger_object = logger_type())
        : ls(std::move(linesearch_object)),
        log(logger_object),
        xg(), xg_p(), correctH(),
        d(), dXdG(), dGdG(), step(),
        iteration(), last_hess(),
        rstate(), go_on(true)
      {
        init_();
      }

      void init(point_type const& p)
      {
        init_();
        xg.p = p;
        xg_p.p = p;
        iteration = 0u;
        xg.update(ls.callback, iteration, go_on);
        d = -xg.p.g;
      }

      point_type const& p() const { return xg.p; }
      point_type& p() { return xg.p; }

      status state() const { return rstate; }
      std::size_t iter() const { return iteration; }

      callback_type&& operator() (point_type& p)
      {
        init(p);
        rstate = optimize();
        return std::move(ls.callback);
      }

      callback_type operator() (point_type& p,
        callback_type callback)
      {
        ls.callback = std::move(callback);
        init(p);
        rstate = optimize();
        return ls.callback;
      }

    };

    template<class LinesearchT>
    inline lbfgs<typename std::remove_reference<LinesearchT>::type>
      make_lbfgs(LinesearchT&& ls)
    {
      return lbfgs<typename std::remove_reference<LinesearchT>::type>
        (std::forward<LinesearchT>(ls));
    }

    template<class LinesearchT, class LoggerT>
    inline lbfgs<typename std::remove_reference<LinesearchT>::type,
      typename std::remove_reference<LoggerT>::type>
      make_lbfgs(LinesearchT&& ls, LoggerT&& log)
    {
      return lbfgs<typename std::remove_reference<LinesearchT>::type,
        typename std::remove_reference<LoggerT>::type>
        (std::forward<LinesearchT>(ls), std::forward<LoggerT>(log));
    }

  }

}

#endif
