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
  template<typename T>
  /**
   * This class enables histogramming of multidimensional data in multidimensional histograms.
   * The regular user should forget about all private member(-functions) and use the class as such:
   *
   * 1. Construct:
   * Either construct with a desired number of bins (adjusting width of bin automatically to suit data)
   * or with a defined width of a bin (adjusting number of bins automatically to suit data)
   *
   * 2. Add Values  via add_value
   * 3. call distribute() to perform histogramming
   * 4. write data and auxilary data.
   */
  class DimensionalHistogram //Single dimension multiple histograms
  {
    // Users should not care about the private members of this class.
    // If you need multidimensional histogramming, just use the public member functions,
    // and voila
  private:

    // width and push value
    std::vector<T> w, p;
    // upper and lower limits
    std::vector<T> m_max, m_min;
    // sum and mean
    std::vector<T> m_s, m_m;
    // values
    std::vector< std::vector<T> >  m_values;
    // boxes
    std::vector< std::size_t > m_boxes;
    // number of histograms
    std::size_t const m_dimensions;
    // total number of values
    std::size_t m_valuecount, m_boxcount;

    //FUNCTIONS

    // Maximum of all values in container for multiple histograms
    std::vector<T> const& maximum() const { return m_max; }
    T const& maximum(unsigned int i) const { return m_max[i]; }

    // Minimum of all values in container for multiple histograms
    std::vector<T> const& minimum() const { return m_min; }
    T const& minimum(unsigned int i) const { return m_min[i]; }

    // Range of all values in container for multiple histograms
    std::vector<T> width() const
    {
      std::vector<T> temp(m_min.size(), 0);
      for (unsigned int i = 0u; i < m_min.size(); i++)
      {
        temp[i] = (m_max[i] - m_min[i]);
      }
      return temp;
    }
    T width(unsigned int i) const { return (m_max[i] - m_min[i]); }

    // Sum of all values in single histogram
    std::vector<T> sum() const { return m_s; }
    T sum(unsigned int i) const { return m_s[i]; }

    // Mean of all values in single histogram
    std::vector<T> mean() const { return m_m; }
    T mean(unsigned int i) const { return m_m[i]; }

    // Center of histogram
    std::vector<T> center() const
    {
      std::vector<T> temp(m_dimensions, 0.);
      for (unsigned int i = 0u; i < m_dimensions; i++)
      {
        temp[i] = minimum(i) + (maximum() - minimum());
      }
      return temp;
    }
    T center(unsigned int dimension) const
    {
      return this->center()[dimension];
    }

    // Value associated with single bin (value at center of bin)
    std::vector<T> centerOfSingleBin(std::vector<size_t> const& bins)
    {
      std::vector<T> binsOut(bins.size());
      if (binsOut.size() == m_dimensions)
      {
        for (unsigned int i = 0u; i < binsOut.size(); i++)
        {
          binsOut[i] = this->minimum()[i] + bins[i] * w[i] + w[i] * 0.5;
        }
        return binsOut;
      }
      else
      {
        std::cerr << "fatal error in histogramming! -> centerOfSingleBin(std::vector<size_t> bins)\n";
        return std::vector<T>();
      }
    }

    // Sets Minimum and Maximum according to Values inside Histogramm
    // Called from within distribute()
    void setMinMax()
    {
      for (unsigned int i = 0u; i < m_values.size(); i++)
      {
        for (unsigned int j = 0u; j < m_values[i].size(); j++)
        {
          if (m_values[i][j] > m_max[j]) m_max[j] = m_values[i][j];
          if (m_values[i][j] < m_min[j]) m_min[j] = m_values[i][j];
        }
      }
    }

    // Access to the serially stored histogram bin values
    // Justification: Since this class should work for any number of dimensions,
    // bin data has to be stored serially. This is an N-Dimension to 1D-Mapping.
    inline std::size_t toIterator(std::vector<std::size_t> const& in)
    {
      if (in.size() == this->dimensions()) 
      {
        size_t accessIterator = 0;
        for (unsigned int i = 0u; i < in.size(); i++)
        {
          accessIterator += in[i] * (size_t) pow(m_boxcount, m_dimensions - i - 1);
        }
        return accessIterator;
      }
      else
      {
        std::cerr << "ERROR IN DIMENSIONS OF HISTOGRAM!\n" << std::flush;
        return 0u;
      }
    }

    //Covert serially stored histogram bin values to dimensional vector
    // Justification: Since this class should work for any number of dimensions,
    // bin data has to be stored serially. This is an N-Dimension to 1D-Mapping.
    inline std::vector<std::size_t> fromIterator(std::size_t in)
    {
      //in std::vec dim 0 is first bin, counting starts at 0.
      std::vector<std::size_t> out;
      out.resize(this->dimensions());

      if (in == 0) return std::vector<std::size_t>(this->dimensions(), 0u);

      size_t tempIn = in;
      for (int i = (int) out.size() - 1; i > -1; i--)
      {
        size_t tempCalc = tempIn % (int) this->numberOfBinsPerDimension();
        out[i] = tempCalc;
        tempIn = tempIn - tempCalc;
        if (tempIn < this->numberOfBinsPerDimension()) 
        { 
          tempIn = this->numberOfBinsPerDimension(); 
        }
        else
        { 
          tempIn = size_t(float(tempIn) / float(this->numberOfBinsPerDimension())); 
        }
      }
      return out;
    }

    inline std::size_t dimensions() const { return m_dimensions; }
    inline std::size_t numberOfBinsPerDimension() const { return m_boxcount; }

    // Returns number of elements in bin associated with (dimensionally stored) "in" 
    inline std::size_t element(std::vector<size_t> const& in) const
    {
      if (in.size() != m_dimensions || m_boxcount <= 0u)
      {
        std::cerr <<   "fatal error in histogramming, desired elements have wrong dimension. Talk to your admin.";
        return 0u;
      }
      else
      {
        return m_boxes[toIterator(in)];
      }
    }

    // Returns number of elements in bin associated with (1D-serial) "accessIterator" 
    // This should actually never be needed to be called exept fopr debug purposes.
    inline std::size_t element(std::size_t accessIterator)
    {
      return m_boxes[accessIterator];
    }


  public:

    // Construct a set of multiple histograms (number of histograms = "number")
    // where each bin has a width of "width". The number of bins is generated dynamicly
    // from this width of one bin and the total range of values later on.
    // "push" is appareantly a value by which the entire histogram is shifted.
    DimensionalHistogram(std::size_t const dimensions, T const width, T const push = T()) :
      w(dimensions, width), p(dimensions, push), m_max(dimensions, 0.), m_min(dimensions, 0.),
      m_s(dimensions, 0.), m_m(dimensions, 0u), m_values(0),
      m_boxes(dimensions, 0u), m_dimensions(dimensions), m_valuecount(0U), m_boxcount(0U) {};

    // Construct a set of multiple histograms (number of histograms = "number")
    // with a specified number of bins.
    // "push" is appareantly a value by which the entire histogram is shifted.
    DimensionalHistogram(std::size_t const dimensions, std::size_t const desiredNumberOfBoxes_, T const push = T()) :
      w(dimensions, 0.), p(dimensions, push), m_max( dimensions, 0.), m_min( dimensions, 0.),
      m_s(dimensions, 0.), m_m(dimensions, 0.), m_values(0),
      m_boxes(dimensions, 0u), m_dimensions(dimensions), m_valuecount(0U), m_boxcount(desiredNumberOfBoxes_) {};

    // Adds a value "newValue" to the histogram nr. "desiredHistogramInSet"
    inline void add_value(std::vector<T> const newValues)
    {
      // CHeck if value has valid dimensions
      if (newValues.size() <= m_dimensions)
      { //CHeck if value = firstValue
        if (m_valuecount < 1U)
        {
          for (unsigned int i = 0u; i < newValues.size(); i++)
          {//set min max
            m_max[i] = newValues[i];
            m_min[i] = newValues[i];
          }
        }
        else
        {
          //set min max
          for (unsigned int i = 0u; i < newValues.size(); i++)
          {
            m_max[i] = std::max(newValues[i], m_max[i]);
            m_min[i] = std::min(newValues[i], m_min[i]);
          }
        }
        //set sum and mean
        for (unsigned int i = 0u; i < newValues.size(); i++)
        {
          T tempsum = sqrt(newValues[i] * newValues[i]);
          m_s[i] += tempsum;
          m_m[i] = m_s[i] / T(m_values.size());
        }
        //save values
        m_values.push_back(newValues);
        ++m_valuecount;
      }
      else
      {
        std::cerr <<  "Fatal ERROR in multidimensional histogramming: add value\n" << std::flush;
      }
    }

    // this seems like it actually counts the values and thereby creates the histogram
    void distribute(void)
    {
      using std::ceil;
      using std::floor;
      for (unsigned int j = 0u; j < m_dimensions; j++)
      {
        m_max[j] += p[j];
      }

      // Get number of distribution boxes. 
      // Number of distribution boxes is equal to ceil()
      // of the total range of values devided by
      // the constructor specified width of one bin
      if (w[0] > T(0.))
      {
        T highestWidth = 0., absMax = 0., absMin=0.;
        for (unsigned int j = 0u; j < m_dimensions; j++)
        {
          highestWidth = std::max(highestWidth, this->width()[j]);
          if (absMax < this->maximum(j)) absMax = this->maximum(j);
          if (absMin > this->minimum(j)) absMin = this->minimum(j);
        }
        m_boxcount = static_cast<std::size_t>(std::ceil((absMax - absMin) / highestWidth));
      }
      else if (m_boxcount == 0)
      {
        std::cerr << "Severe Error in histogram.h. Stopping.";
        throw;
      }
      else if (m_boxcount > 0)
      {    
        for (unsigned int i = 0u; i < m_dimensions; i++)
        {
          w[i] = (m_max[i] - m_min[i]) / T(m_boxcount);
        }
      }

      // We need m_boxcount^dimension discrete histogram bins, so we resize to this
      //Adjust number of histogram bins
      m_boxes.resize( (size_t) pow(m_boxcount, m_dimensions));
      //Fills the std::vector with "m_boxcount^m_dimension" size_ts with value 0u.
      m_boxes.assign( (size_t) pow(m_boxcount, m_dimensions), 0u);

      // Now, the actual processing of values takes place.
      std::size_t const numberOfValues = m_values.size();
      for (std::size_t i = 0; i < (size_t) numberOfValues; ++i) // Iterate over the number of values. i is current value.
      {
        std::vector<std::size_t> index(m_dimensions, size_t());
        for (std::size_t h(0U); h < m_dimensions; ++h)
        {
          index[h] = static_cast<std::size_t>(floor((m_values[i][h] - minimum()[h]) / w[h] + p[h]));
          if (index[h] >= this->numberOfBinsPerDimension())
          {
            index[h] = this->numberOfBinsPerDimension() - 1u;
          }
          else if (index[h] < 0u)
          {
            index[h] = 0u;
          }
        }
        ++m_boxes[this->toIterator(index)];
      }
    }

    // Write histogrammed data in a format suitable for gnuplotting
    void write(std::string filename)
    {
      std::ofstream stream(filename, std::ios::out);
      stream << std::right << std::setw(13);
      size_t keeperForBlanklines = 0;
      for (unsigned int i = 0u; i < m_boxes.size(); i++)
      {

        std::vector<T> bins = centerOfSingleBin(fromIterator(i));

        if (fromIterator(i)[bins.size() - 2] != keeperForBlanklines)
        {
          keeperForBlanklines = fromIterator(i)[bins.size() - 2];
          stream << "\n";
        }

        for (unsigned int j = 0u; j < bins.size(); j++)
        {
          stream << std::right << std::setw(13) << bins[j] <<  " ";
        }
        stream << std::right << std::setw(13) << this->element((size_t)i) << "\n"; 
      }
    }

    // Write histogrammed data in a format suitable for gnuplotting
    // This data will be normated to be a probability density
    void writeProbabilityDensity(std::string filename)
    {
      std::ofstream stream(filename, std::ios::out);
      stream << std::right << std::setw(13);
      size_t keeperForBlanklines = 0;
      for (unsigned int i = 0u; i < m_boxes.size(); i++)
      {

        std::vector<T> bins = centerOfSingleBin(fromIterator(i));

        if (fromIterator(i)[bins.size() - 2] != keeperForBlanklines)
        {
          keeperForBlanklines = fromIterator(i)[bins.size() - 2];
          stream << "\n";
        }

        for (unsigned int j = 0u; j < bins.size(); j++)
        {
          stream << std::right << std::setw(13) << bins[j] << " ";
        }
        stream << std::right << std::setw(13) << (double)((double) this->element((size_t)i) / (double) this->m_valuecount) << "\n";
      }
    }

    //Write sums and means (->auxilary data) to file.
    void writeAuxilaryData(std::string filename)
    {
      std::ofstream stream(filename, std::ios::out);
      stream << "Number of Dimensions: " << this->dimensions() << "\n";
      for (unsigned int i = 0; i < this->dimensions(); i++)
      {
        stream << std::right << "Mean of dimension " << i+1 << ": " << this->mean(i) << "\n";
      }
      stream << "\n\n";

      for (unsigned int i = 0; i < this->dimensions(); i++)
      {
        stream << std::right << "Sum of dimension " << i + 1 << ": " << this->sum(i) << "\n";
      }
      stream << "\n\n";
    }
  };

  ///////////////////////////////////////////////////////////////
  //////////////////////////////
  ////////////////////////////////////////////////////////

  template<typename T>
  class Histograms //Single dimension multiple histograms
  {

  private:

    // width and push value
    T w, p;
    // upper and lower limits
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

    // Construct a set of multiple histograms (number of histograms = "number")
    // where each bin has a width of "width". The number of bins is generated dynamicly
    // from this width of one bin and the total range of values later on.
    // "push" is appareantly a value by which the entire histogram is shifted.
    Histograms(std::size_t const number, T const width, T const push = T()) :
      w(width), p(push), m_max(), m_min(), m_s(number), m_m(number), m_values(number),
      m_boxes(number), m_histograms(number), m_valuecount(0U), m_boxcount(0U) {};

    // Construct a set of multiple histograms (number of histograms = "number")
    // with a specified number of bins.
    // "push" is appareantly a value by which the entire histogram is shifted.
    Histograms(std::size_t const number, std::size_t const desiredNumberOfBoxes_, T const push = T()) :
      w(0.), p(push), m_max(), m_min(), m_s(number), m_m(number), m_values(number),
      m_boxes(number), m_histograms(number), m_valuecount(0U), m_boxcount(desiredNumberOfBoxes_) {};

    // Adds a value "newValue" to the histogram nr. "desiredHistogramInSet"
    void add_value(std::size_t const desiredHistogramInSet, T const newValue)
    {
      using std::max;
      using std::min;
      using std::ceil;
      using std::floor;
      if (desiredHistogramInSet < m_histograms)
      {
        if (m_valuecount < 1U)
        {
          m_max = ceil(newValue);
          m_min = floor(newValue);
        }
        else
        {
          m_max = max(ceil(newValue), m_max);
          m_min = max(floor(newValue), m_min);
        }
        m_s[desiredHistogramInSet] += newValue;
        m_values[desiredHistogramInSet].push_back(newValue);
        m_m[desiredHistogramInSet] = m_s[desiredHistogramInSet] / T(m_values[desiredHistogramInSet].size());
        ++m_valuecount;
      }
    }

    // Maximum of all values in container for multiple histograms
    T const& maximum() const { return m_max; }
    // Minimum of all values in container for multiple histograms
    T const& minimum() const { return m_min; }

    // Sets Minimum and Maximum according to Values inside Histogramm
    void setMinMax()
    {
      for (unsigned int i = 0u; i < m_histograms; i++)
      {
        for (unsigned int j = 0u; j < m_values[i].size(); j++)
        {
          if (m_values[i][j] > m_max) m_max = m_values[i][j];
          if (m_values[i][j] < m_min) m_min = m_values[i][j];
        }
      }
      if (m_boxcount > 0)
      {
        w = (m_max - m_min) / T(m_boxcount);
      }
    }

    // Range of all values in container for multiple histograms
    T width() const { return (m_max - m_min); }

    // Sum of all values in single histogram
    T sum(std::size_t const desiredHistogramInSet) const { return desiredHistogramInSet < m_s.size() ? m_s[desiredHistogramInSet] : T(); }

    // Mean of all values in single histogram
    T mean(std::size_t const desiredHistogramInSet) const { return desiredHistogramInSet < m_m.size() ? m_m[desiredHistogramInSet] : T(); }

    // Center of histogram bin nr. "numberOfHistogramBin"
    T center(std::size_t const numberOfHistogramBin) const { return numberOfHistogramBin < boxes() ? minimum() + numberOfHistogramBin*w + 0.5*w + p : T(); }


    std::size_t histograms() const { return m_histograms; }
    std::size_t boxes() const { return m_boxcount; }

    // Returns number of elements in bin nr. "box" of histogram nr. "desiredHistogramInSet"
    std::size_t elements(std::size_t const desiredHistogramInSet, std::size_t const box) const
    {
      return (desiredHistogramInSet < m_values.size() && box < m_boxes[desiredHistogramInSet].size()) ? m_boxes[desiredHistogramInSet][box] : 0U;
    }

    // Returns the standard deviation of histogram "desiredHistogramInSet" in set.
    double std_deviation(std::size_t const desiredHistogramInSet)
    {
      using std::sqrt;
      // Check if input is valid in range.
      if (m_values[desiredHistogramInSet].size() > 1)
      {
        T square_dev, square_dev_sum(0);
        for (std::size_t j(0u); j < m_values[desiredHistogramInSet].size(); ++j)
        {
          square_dev = m_values[desiredHistogramInSet][j] - m_m[desiredHistogramInSet];
          square_dev *= square_dev;
          square_dev_sum += square_dev;
        }
        return sqrt(square_dev_sum / T(m_values[desiredHistogramInSet].size() - 1U));
      }
      return T();
    }

    // this seems like it actually counts the values (again?)
    // and thereby creates the histogram.e Not sur though.
    void distribute(void)
    {
      this->setMinMax();
      using std::ceil;
      using std::floor;
      m_max += p;

      // Get number of distribution boxes. 
      // Number of distribution boxes is equal to ceil()
      // of the total range of values devided by
      // the constructor specified width of one bin
      if (w > T(0.))
      {
        m_boxcount = static_cast<std::size_t>(ceil(width() / w));
      }
      else if (m_boxcount == 0)
      {
        std::cerr << "Severe Error in histogram.h. Stopping.";
        throw;
      }
      for (std::size_t h(0U); h < m_histograms; ++h)
      {
        m_boxes[h].resize(m_boxcount); //Adjust number of histogram bins
        m_boxes[h].assign(m_boxcount, 0U); //Fills the std::vector with "m_boxcount" size_ts with value 0u.
        std::size_t const VALS = m_values[h].size();

        for (std::size_t i = 0; i < VALS; ++i) // Iterate over the number of values. i is current value.
        {
          // index = floor() of ( current Value - global minimum value in set of all histograms) over ( size of one bin) + push )
          std::size_t const index = static_cast<std::size_t>(floor((m_values[h][i] - minimum()) / w + p));

          if (index < m_boxes[h].size()) //If index is smaller than current number of histogram bins.
          {
            ++m_boxes[h][index]; // Increment count for the bin corresponding to "index"
          }
          else if (index > 0U && index == m_boxes[h].size()) //Else if index is euqal to current number of histogram bins.
          {
            ++m_boxes[h][index - 1U]; // Increment count for the bin corresponding to "index - 1u"
          }
        }
      }
    }

  };

  // Formatted output
  template<class T>
  std::ostream& operator<< (std::ostream &stream, Histograms<T> const &histo)
  {
    std::size_t const N(histo.boxes()), M(histo.histograms());
    stream << std::right << std::setw(10) << "Mean";
    stream << std::right << std::setw(13) << " ";
    for (std::size_t h(0U); h < M; ++h)
    {
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.mean(h);
    }
    stream << lineend;
    stream << std::right << std::setw(10) << "Sum";
    stream << std::right << std::setw(13) << " ";
    for (std::size_t h(0U); h < M; ++h)
    {
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.sum(h);
    }
    stream << lineend;
    for (std::size_t b(0U); b < N; ++b)
    {
      stream << std::right << std::setw(10) << b << " ";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.center(b);
      for (std::size_t h(0U); h < M; ++h)
      {
        stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << histo.elements(h, b);
      }
      stream << std::endl;
    }
    return stream;
  }
}
#endif // histogram_774418c7_21fd_4101_a5f2_1c91e51ff62b_h_guard_
