#ifndef histogram_774418c7_21fd_4101_a5f2_1c91e51ff62b_h_guard_
#define histogram_774418c7_21fd_4101_a5f2_1c91e51ff62b_h_guard_  
#pragma once

#include <vector>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <map>
#include <limits>



namespace statistics
{

  template<class T>
  using floif = typename std::enable_if<
    std::is_floating_point<T>::value, T>::type;

  namespace kernel
  {

    template<class T>
    inline floif<T> indicator(bool const x)
    {
      return x ? T{ 1 } : T{ 0 };
    }

    template<class T>
    inline floif<T> f_uniform(T u)
    {
      using std::abs;
      return 0.5*indicator<T>(abs(u) <= 1);
    }

    template<class T>
    floif<T> f_triangular(T u)
    {
      using std::abs;
      auto const au = abs(u);
      return (T{ 1 }-au)*indicator<T>(abs(u) <= 1);
    }

    template<class T>
    floif<T> f_epanechnikov(T u)
    {
      using std::abs;
      return (T{ 0.75 }-u*u)*indicator<T>(abs(u) <= 1);
    }

    template<class T>
    floif<T> f_quartic(T u)
    {
      using std::abs;
      return (T{ 0.9375 }-u*u)*indicator<T>(abs(u) <= 1);
    }

    template<class T>
    floif<T> f_gaussian(T u)
    {
      using std::exp;
      return T{ 0.3989422804014326779399460 }*exp(-T{ 0.5 }*u*u);
    }

  }

  namespace selector
  {

    namespace _detail
    {

      template<class T>
      struct accum
      {
        T sum;
        accum() : sum() {}
        void operator() (T const v)
        {
          sum += v;
        }
      };

      template<class T>
      struct ssd
      {
        T sum;
        T n;
        T av;
        T result;
        ssd(T average) : sum(), n(0), av(average), result() {}
        void operator() (T const v)
        {
          using std::sqrt;
          auto d = v - av;
          sum += d*d;
          n += T{ 1 };
          result = std::sqrt(sum / (n - T{ 1 }));
        }
      };

    }


    template<class T>
    floif<T> normal_distribution_approximation(std::vector<T> const &sample)
    {
      auto const av = std::for_each(sample.begin(), sample.end(),
        _detail::accum<T>{}).sum / static_cast<T>(sample.size());
      auto sd = std::for_each(sample.begin(), sample.end(), _detail::ssd<T>{av}).result;
      auto sd4 = sd*sd; // ^2
      sd4 = sd*sd; // ^4
      return std::pow(T{ 4 } *sd4 * sd / (T{ 3 } *static_cast<T>(sample.size())), T{ 1 } / T{ 5 });
    }
  }

  template<class T, class K, bool is_ari = std::is_floating_point<T>::value>
  class kernel_density_estimator;

  template<class T, class K>
  class kernel_density_estimator<T, K, true>
  {

    std::vector<T> sample;
    K kern;
    T bw;

  public:

    template<class U>
    kernel_density_estimator(U && k, std::vector<T> const & data = {}) 
      : sample(data), kern(std::forward<U>(k)), bw() {}

    template<class S>
    void select_bandwidth(S && selector_object)
    {
      using std::abs;
      bw = abs(selector_object(sample));
    }

    T bandwidth() const { return bw; }

    T operator() (T const x) const
    {
      T sum = T{};
      for (auto s : sample)
      {
        sum += kern( (x - s) / bw );
      }
      return sum / (static_cast<T>(sample.size())*bw);
    }

  };

  template<class T, class K, class S>
  inline kernel_density_estimator<T, K> make_kde(
    std::vector<T> const &values, K && kernel_object, S && selector_object)
  {
    auto kde = kernel_density_estimator<T, K>(std::forward<K>(kernel_object), values);
    kde.select_bandwidth(std::forward<S>(selector_object));
    return kde;
  }

  template<class T>
  inline auto make_default_kde(std::vector<T> const &values) -> 
    decltype(make_kde(values, &kernel::f_epanechnikov<T>, 
      &selector::normal_distribution_approximation<T>))
  {
    return make_kde(values, &kernel::f_epanechnikov<T>,
      &selector::normal_distribution_approximation<T>);
  }


  template<class T, bool isfp = std::is_floating_point<T>::value>
  class histogram;

  template<class T>
  class histogram<T, true>
  {

    T width, mini, maxi;
    std::map<std::size_t, std::size_t> histo_bins;
    
    std::map<std::size_t, std::vector<std::size_t>> bin_contains;

  public:

    using iterator = std::map<std::size_t, std::size_t>::const_iterator;

    histogram(std::vector<T> const &sample, T const bin_width) 
      : width(bin_width), mini(std::numeric_limits<T>::max()),
      maxi(std::numeric_limits<T>::min()), histo_bins(), bin_contains()
    {
      for (auto s : sample)
      {
        mini = std::min(mini, s);
        maxi = std::max(maxi, s);
      }
      std::size_t i = 0;
      for (auto s : sample)
      {
        auto bin_id = static_cast<std::size_t>(
          std::floor((s - mini) / width));
        ++histo_bins[bin_id];
        bin_contains[bin_id].push_back(i);
        ++i;
      }
    }

    struct bin_type
    {
      T const start, end;
      std::size_t const n;
      std::vector<std::size_t> const & elements;
    };

    bin_type bin(std::size_t const i)
    {
      auto strt = mini + i*width;
      return{ strt, strt + width, histo_bins[i], bin_contains[i] };
    }

    std::size_t bins() const 
    {
      return static_cast<std::size_t>(
        std::ceil((maxi - mini) / width));
    }

    std::map<std::size_t, std::size_t>::size_type size() const
    {
      return histo_bins.size();
    }

    iterator begin() const { return histo_bins.begin(); }
    iterator end() const { return histo_bins.end(); }

  };


}

namespace histo
{

  

  template<class T>
  class multidimensional
  {
    std::vector<T> m_bin_width, m_bin_shift,
      m_maximum, m_minimum;

    std::vector<std::size_t> m_bin_count;

    std::vector<std::vector<T>> m_values;
    std::vector<std::vector<std::size_t>> m_num_box_values;



  };
}

template<typename T>
class Multihistogram
{

private:

  // width and push value
  T w, p;
  // limits
  T m_max, m_min;
  // sum and mean
  std::vector<T> m_s, m_m;
  // values
  std::vector< std::vector<T> >  m_values;
  // boxes
  std::vector< std::vector<std::size_t> > m_boxes;
  // number of histograms
  std::size_t const m_histograms;
  // total number of values
  std::size_t m_valuecount, m_boxcount;

public:

  Multihistogram (std::size_t const number, T const width, T const push) :
    w(width), p(push), m_max(), m_min(), m_s(number), m_m(number), m_values(number), 
    m_boxes(number), m_histograms(number), m_valuecount(0U), m_boxcount(0U)
  { }

  void add_value (std::size_t const histogram, T const v)
  {
    using std::max;
    using std::min;
    using std::ceil;
    using std::floor;
    if (histogram < m_histograms)
    {
      if (m_valuecount < 1U) 
      {
        m_max = ceil(v);
        m_min = floor(v);
      }
      else
      {
        m_max = max(ceil(v), m_max);
        m_min = max(ceil(v), m_min);
      }
      m_s[histogram] += v;
      m_values[histogram].push_back(v);
      m_m[histogram] = m_s[histogram] / T(m_values[histogram].size());
      ++m_valuecount;
    }
  }

  T const& maximum () const { return m_max; }
  T const& minimum () const { return m_min; }

  T width () const { return (m_max-m_min); }
  T sum(std::size_t const i) const { return i < m_s.size() ? m_s[i] : T(); }
  T mean(std::size_t const i) const { return i < m_m.size() ? m_m[i] : T(); }
  T center (std::size_t const i) const { return i < boxes() ? minimum() + i*w+0.5*w+p : T(); }

  std::size_t histograms () const { return m_histograms; }
  std::size_t boxes () const { return m_boxcount; }

  std::size_t elements (std::size_t const histogram, std::size_t const box) const 
  { 
    return (histogram < m_boxes.size() && box < m_boxes[histogram].size()) ? m_boxes[histogram][box] : 0U; 
  }

  double std_deviation (std::size_t const i)
  {
    using std::sqrt;
    if (i < m_histograms && m_values[i].size() > 1)
    {
      T square_dev, square_dev_sum(0);
      for (std::size_t j(0u); j<m_values[i].size(); ++j)
      {
        square_dev = m_values[i][j]-m_m[i];
        square_dev *= square_dev;
        square_dev_sum += square_dev;
      }
      return sqrt(square_dev_sum / T(m_values[i].size()-1U));
    }
    return T();
  }

  void distribute (void)
  {
    using std::ceil;
    using std::floor;
    m_max += p;
    // Get number of distribution boxes
    m_boxcount = static_cast<std::size_t>(ceil(width() / w));
    for (std::size_t h(0U); h<m_histograms; ++h)
    {
      m_boxes[h].resize(m_boxcount);
      m_boxes[h].assign(m_boxcount, 0U);
      std::size_t const VALS(m_values[h].size());
      for (std::size_t i=0; i<VALS; ++i)
      {
        std::size_t const index = static_cast<std::size_t>(floor((m_values[h][i] - minimum()) / w + p));
        if (index < m_boxes[h].size()) ++m_boxes[h][index];
        else if (index > 0U && index == m_boxes[h].size()) ++m_boxes[h][index-1U];
      }
    }
  }

};

template<class T>
std::ostream& operator<< (std::ostream &stream, Multihistogram<T> const &histo)
{
  std::size_t const N(histo.boxes()), M(histo.histograms());
  stream << std::right << std::setw(10) << "Mean";
  stream << std::right << std::setw(13) << " ";
  for (std::size_t h(0U); h<M; ++h)
  {
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.mean(h);
  }
  stream << lineend;
  stream << std::right << std::setw(10) << "Sum";
  stream << std::right << std::setw(13) << " ";
  for (std::size_t h(0U); h<M; ++h)
  {
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.sum(h);
  }
  stream << lineend;
  for (std::size_t b(0U); b<N; ++b)
  {
    stream << std::right << std::setw(10) << b << " ";
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.center(b);
    for (std::size_t h(0U); h<M; ++h)
    {
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.elements(h, b);
    }
     stream << std::endl;
  }
  return stream;
}

#endif // histogram_774418c7_21fd_4101_a5f2_1c91e51ff62b_h_guard_
