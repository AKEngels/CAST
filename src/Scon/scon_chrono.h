#pragma once

#include <cstdint>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "scon.h"
#include "scon_utility.h"

#if defined(SCON_CC11_CHRONO)
  #include <chrono>
#else
  #include <ctime>
#endif

#if defined (_MSC_VER) && _MSC_VER < 1900
#define WIN32_LEAN_AND_MEAN
#define NOATOM
#define NOGDI
#define NOGDICAPMASKS
#define NOMETAFILE
#define NOMINMAX
#define NOMSG
#define NOOPENFILE
#define NORASTEROPS
#define NOSCROLL
#define NOSOUND
#define NOSYSMETRICS
#define NOTEXTMETRIC
#define NOWH
#define NOCOMM
#define NOKANJI
#define NOCRYPT
#define NOMCX
//#define UNICODE
#define STRICT
#include <windows.h>
#include <ratio>
#undef WIN32_LEAN_AND_MEAN
#undef NOATOM
#undef NOGDI
#undef NOGDICAPMASKS
#undef NOMETAFILE
#undef NOMINMAX
#undef NOMSG
#undef NOOPENFILE
#undef NORASTEROPS
#undef NOSCROLL
#undef NOSOUND
#undef NOSYSMETRICS
#undef NOTEXTMETRIC
#undef NOWH
#undef NOCOMM
#undef NOKANJI
#undef NOCRYPT
#undef NOMCX
#undef UNICODE
#undef STRICT
#endif


namespace scon 
{

  namespace chrono
  {

#if defined(_MSC_VER) && _MSC_VER < 1900

    struct high_resolution_clock
    {	// wraps QueryPerformanceCounter
      typedef long long rep;
      typedef std::nano period;
      typedef std::chrono::nanoseconds duration;
      typedef std::chrono::time_point<high_resolution_clock> time_point;
      static const bool is_steady = true;

      static time_point now() _NOEXCEPT
      {	// get current time
        static const long long _Freq
          = _Query_perf_frequency();	// doesn't change after system boot
        const long long _Ctr = _Query_perf_counter();
        static_assert(period::num == 1, "This assumes period::num == 1.");
        const long long _Whole = (_Ctr / _Freq) * period::den;
        const long long _Part = (_Ctr % _Freq) * period::den / _Freq;
        return (time_point(duration(_Whole + _Part)));
      }
    };
#else
    using high_resolution_clock = std::chrono::high_resolution_clock;
#endif

    template<class Clock = high_resolution_clock>
    class high_resolution_timer_templ
    {

    public:

      using duration = typename Clock::duration;
      using time_point = typename Clock::time_point;

    protected:

      time_point start;

    public:

      high_resolution_timer_templ() : start(Clock::now()) { }
      duration operator() () const
      {
        return Clock::now() - start;
      }
      void restart() { start = Clock::now(); }

    }; // high_resolution_timer

    using high_resolution_timer = high_resolution_timer_templ<high_resolution_clock>;

    template<class Rep, class Duration>
    static inline double to_seconds (Duration && t)
    {
      using dura = std::chrono::duration<Rep>;
      auto sec_dura = std::chrono::duration_cast<dura>(std::forward<Duration>(t));
      return sec_dura.count();
    }

    template<class Clock>
    inline std::ostream& operator<< (std::ostream &strm, high_resolution_timer_templ<Clock> const &timer)
    {
      auto d = timer();
      strm << to_seconds<double>(d) << " s (" << d.count() << " ticks)";;
      return strm;
    }


  } // namespace chrono

}
