#ifndef conjugate_gradient_header

#define conjugate_gradient_header

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

#include <vector>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include "ls.h"

namespace optimization
{
  namespace local
  {

    template<class LineSearchT, class LoggerT = empty_void_functor>
    class conjugate_gradient
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

      linesearch_type ls;
      logger_type log;
     
      struct configuration
      {
        // maximum number of iterations
        std::size_t max_iterations, reset_offset;
        // Convergence epsilon
        float_type epsilon;
        configuration() :
          max_iterations(1000), reset_offset(10u), 
          epsilon(float_type(1.0e-5))
        { }
      } config;

      conjugate_gradient(linesearch_type linsearch_object, 
        logger_type logger_object = logger_type())
        : ls(std::move(linsearch_object)), log(std::move(logger_object)), 
        config(), stored(), iteration(0u), rstate(status::UNDEFINED)
      { }

      void init(point_type const & p = point_type())
      {
        rstate = status::UNDEFINED;
        stored = p;
        iteration = 0u;
      }

      point_type const & p() const { return stored; }
      point_type & p() { return stored; }

      status state() const { return rstate; }
      std::size_t iter() const { return iteration; }

      callback_type operator() (point_type & p)
      {
        init(p);
        rstate = optimize();
        return ls.callback;
      }

      callback_type operator() (point_type & p,
        callback_type callback)
      {
        ls.callback = std::move(callback);
        init(p);
        rstate = optimize();
        return ls.callback;
      }

    private:

      using F = float_type;

      point_type stored;
      std::size_t iteration;
      status rstate;

      bool convergence()
      {
        using std::max;
        using std::sqrt;
        auto const rdg = sqrt(dot(stored.g, stored.g));
        auto const rdx = sqrt(dot(stored.x, stored.x));
        return ((rdg / max(rdx, F(1))) < config.epsilon);
      }

      status optimize()
      {
        using std::sqrt;
        using std::min;
        iteration = 0U;
        point_type p(stored);
        grad_type G, G_old, d;
        bool go_on = true;
        p.f = ls.callback(p.x, p.g, iteration, go_on);
        stored = p;
        if (!go_on) return rstate = status::ERR_CALLBACK_STOP;
        float_type step = F(1) / sqrt(dot(p.g, p.g));
        for (; iteration < config.max_iterations;)
        {
          // Get step direction from conjugate gradient formulas
          if (iteration < 1U || (iteration%config.reset_offset) == 0)
          {
            G = p.g;
            d = -G;
            G_old = G;
          }
          else
          {
            auto const ddg = dot((p.g - stored.g), p.g);
            auto const dgg = dot(p.g, p.g);
            auto const gamma = ddg / dgg;
            auto const rdgo = sqrt(dot(G_old, G_old));
            G = p.g + G_old * gamma*rdgo;
            d = -G;
            G_old = G;
            step = 1.0;
          }
          // save current state before linesearch
          stored = p;
          // linesearch
          auto const lsr = ls(d, step, p, stored.x, iteration);
          // return error if linesearch fails
          if (lsr != status::SUCCESS) return lsr;
          // save current state and return if converged
          if (convergence())
          {
            stored = p;
            return rstate = status::SUCCESS;
          }
          // reset steplength
          step = 1.0;
          ++iteration;
        }
        return rstate = status::ERR_MAX_ITERATIONS;
      }

    };

    template<class LinesearchT>
    inline conjugate_gradient<LinesearchT> make_cg(LinesearchT && ls)
    {
      return conjugate_gradient<LinesearchT>(std::forward<LinesearchT>(ls));
    }

  }

}

#endif
