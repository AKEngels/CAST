#ifndef optimization_local_header
#define optimization_local_header

#include <vector>
#include <stdexcept>
#include <iterator>
#include <cstddef>
#include <cmath>
#include <utility>
#include <algorithm>

#include "lbfgs.h"
#include "ls.h"
#include "scon_vect.h"
#include "scon_traits.h"

namespace optimization
{

  namespace local
  {

    template<class underlying_callback>
    struct ortho_callback 
    {

      using callback_type = underlying_callback;
      using float_type = scon::return_type < callback_type > ;
      using rep_type = typename std::decay < scon::argument_type<callback_type, 0U> >::type;
      using grad_type = typename std::decay < scon::argument_type<callback_type, 1U> >::type;

      rep_type direction;
      underlying_callback callback;

      template<class... Args>
      ortho_callback(typename std::enable_if<std::is_constructible<underlying_callback,Args...>::value, 
                     rep_type>::type const &ortho_direction, Args && ... args)
        : direction(ortho_direction), callback(std::forward<Args>(args)...)
      {
        using std::sqrt;
        // normalize direction
        direction /= sqrt(dot(direction, direction));
      }

      float_type operator() (rep_type const & v, grad_type & g,
                             std::size_t const S, bool & go_on)
      {
        auto const E = callback(v, g, S, go_on);
        if (go_on) g -= direction * dot(g, direction);
        std::cout << "g = " << g << std::endl;
        return E;
      }

    };

    template<class callback_type, class rep_type>
    ortho_callback<callback_type> make_ortho_callback
      (callback_type && callback, rep_type && rep)
    {
      return ortho_callback<callback_type>(rep, callback);
    }

    template<class image_callback_type>
    struct neb_callback
    {

      using callback_type = image_callback_type;
      using float_type = scon::return_type < callback_type > ;
      using rep_type = typename std::decay < scon::argument_type<callback_type, 0U> >::type ;
      using grad_type = typename std::decay < scon::argument_type<callback_type, 1U> >::type;
      using rep_vec = scon::vector < rep_type > ;
      using grad_vec = scon::vector < grad_type > ;
      using float_vec = scon::vector < float_type > ;

      //static const float_type v = 1.5;

      callback_type image_callback;
      float_type distance;

      neb_callback(callback_type image_callback_obj)
        : image_callback(image_callback_obj), distance()
      { }

      rep_vec setup(std::size_t const num_images, rep_type const &start, rep_type const &end)
      {
        rep_vec images(num_images);
        using std::sqrt;
        auto const delta = (end - start) / (float_type)(num_images - 1);
        //std::cout << "delta = " << delta << std::endl;
        float_type const dist = sqrt(dot(delta, delta));
        //std::cout << "distance = " << dist << std::endl;
        distance = dist;
        for (std::size_t i(0u); i < num_images; ++i)
        {
          images[i] = start + delta * (float_type)i;
        }
        return images;
      }

      rep_vec setup(std::size_t const num_images, rep_vec const &vertices)
      {
        using std::sqrt;
        std::size_t const vs = vertices.size();
        auto dist_acc = float_type();
        rep_vec v_delta(vs - 1U);
        float_vec v_dist(vs - 1U);
        for (std::size_t i = 1u; i < vs; ++i)
        {
          std::size_t k = i - 1U;
          v_delta[k] = vertices[i] - vertices[k];
          v_dist[k] = sqrt(dot(v_delta[k], v_delta[k]));
          dist_acc += v_dist[k];
        }
        // mean inter image distance
        distance = dist_acc / (float_type)(num_images - 1U);
        rep_vec images;
        for (std::size_t i = 1u; i < vs; ++i)
        {
          std::size_t k = i - 1U;
          auto const pieces = std::round(v_dist[k] / distance);
          auto const pcs = static_cast<std::size_t>(pieces);
          v_delta[k] /= pieces;
          images.reserve(images.size() + pcs);
          for (std::size_t j = 0u; j < pcs; ++j)
          {
            images.push_back(vertices[k] + v_delta[k] * (float_type)j);
          }
        }
        if (vertices.size() > 1U) images.push_back(vertices.back());
        distance = dist_acc / (float_type)(images.size() - 1U);
        //std::cout << images.size() << " :: " << num_images << std::endl;
        return images;
      }

    };
    
    template<class rep_vec, class grad_vec, class float_type>
    struct neb_analyzer
    {

    };



    template<class image_callback_type>
    struct neb_quadratic_callback : neb_callback<image_callback_type>
    {

      // types
      using base_type = neb_callback < image_callback_type >;
      using callback_type = typename base_type::callback_type;
      using float_type = typename base_type::float_type;
      using rep_type = typename base_type::rep_type;
      using grad_type = typename base_type::grad_type;
      using rep_vec = typename base_type::rep_vec;
      using grad_vec = typename base_type::grad_vec;
      using float_vec = typename base_type::float_vec;
      
      // members
      float_type force;

      // ctor
      neb_quadratic_callback(callback_type images_callback, float_type const quad_force)
        : base_type(images_callback), force(quad_force)
      { }

      // function call operator
      float_type operator() (rep_vec const & v, grad_vec & g, 
                             std::size_t const S, bool & go_on)
      {
        using std::begin;
        using std::end;
        using std::sqrt;
        //using scon::dot;
        std::size_t const N = v.size();
        if (g.size() != N)
        {
          g.resize(N);
        }
        // Get energy of each image
        float_type E = float_type();
        for (std::size_t i(0U); i < N && go_on; ++i)
        {
          float_type const ie = this->image_callback(v[i], g[i], S, go_on);
          E += ie;
        }
        // Get inter-image energy
        if (go_on)
        {
          std::size_t const M = N - 1U;
          for (std::size_t i(1U); i < N; ++i)
          {
            std::size_t const k = i - 1u;
            rep_type const d_rep(v[i] - v[k]);
            float_type const d = sqrt(dot(d_rep, d_rep));
            float_type const delta = d - this->distance;
            if (delta > (float_type)0)
            {
              float_type dE = force * delta;
              E += dE * delta;
              dE *= (float_type)2 / d;
              grad_type const d_grad = d_rep*dE;
              if (i < M) g[i] += d_grad;
              if (i > 1U) g[k] -= d_grad;
            }
          }
          g[0U] *= (float_type)0.;
          g[N - 1U] *= (float_type)0.;
        }
        return E;
      }

    };

    template<class im_callb_ty>
    neb_quadratic_callback<im_callb_ty> 
      make_neb_quad_callback(im_callb_ty im_callback,
      scon::return_type<im_callb_ty> const force)
    {
      return neb_quadratic_callback<im_callb_ty>(im_callback, force);
    }

  }

}

#endif // optimization_local_header
