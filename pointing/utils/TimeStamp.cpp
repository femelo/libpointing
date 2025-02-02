 /* -*- mode: c++ -*-
 *
 * pointing/utils/TimeStamp.cpp --
 *
 * Initial software
 * Authors: Nicolas Roussel, Gery Casiez
 * Copyright © Inria
 *
 * http://libpointing.org/
 *
 * This software may be used and distributed according to the terms of
 * the GNU General Public License version 2 or any later version.
 *
 */

#define __STDC_LIMIT_MACROS 1

#include <pointing/utils/TimeStamp.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <time.h>

#ifndef _MSC_VER // Visual Studio C++
#include <sys/time.h>
#else
#define timegm _mkgmtime
#endif
#ifndef timegm
#define timegm(arg) (mktime(arg) - timezone)
#endif

#ifdef __APPLE__
#include <mach/mach_time.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

namespace pointing {

  const TimeStamp::inttime TimeStamp::one_nanosecond = 1LL ;
  const TimeStamp::inttime TimeStamp::one_microsecond = 1000LL*TimeStamp::one_nanosecond ;
  const TimeStamp::inttime TimeStamp::one_millisecond = 1000LL*TimeStamp::one_microsecond ;
  const TimeStamp::inttime TimeStamp::one_second = 1000LL*TimeStamp::one_millisecond ;
  const TimeStamp::inttime TimeStamp::one_minute = 60LL*TimeStamp::one_second ;
  const TimeStamp::inttime TimeStamp::one_hour = 60LL*TimeStamp::one_minute ;
  const TimeStamp::inttime TimeStamp::one_day = 24LL*TimeStamp::one_hour ;
  const TimeStamp::inttime TimeStamp::one_week = 7LL*TimeStamp::one_day ;

  const TimeStamp::inttime TimeStamp::undef = INT64_MIN ;
  const TimeStamp::inttime TimeStamp::min = -2147483648LL*TimeStamp::one_second ;
  const TimeStamp::inttime TimeStamp::max = 2147483647LL*TimeStamp::one_second ;

  // ----------------------------------------------------------------------

#ifdef _MSC_VER // Visual Studio C++
#include <Windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
 
  struct timezone 
  {
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
  };
 
  int gettimeofday(struct timeval *tv, struct timezone *tz)
  {
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;
 
    if (NULL != tv)
      {
	GetSystemTimeAsFileTime(&ft);
 
	tmpres |= ft.dwHighDateTime;
	tmpres <<= 32;
	tmpres |= ft.dwLowDateTime;
 
	/*converting file time to unix epoch*/
	tmpres /= 10;  /*convert into microseconds*/
	tmpres -= DELTA_EPOCH_IN_MICROSECS; 
	tv->tv_sec = (long)(tmpres / 1000000UL);
	tv->tv_usec = (long)(tmpres % 1000000UL);
      }
 
    if (NULL != tz)
      {
	if (!tzflag)
	  {
	    _tzset();
	    tzflag++;
	  }
	tz->tz_minuteswest = _timezone / 60;
	tz->tz_dsttime = _daylight;
      }
 
    return 0;
  }
#endif

  // ----------------------------------------------------------------------

  TimeStamp::TimeStamp(const TimeStamp::inttime t) {
    if ((t>=TimeStamp::min && t<=TimeStamp::max) || t==TimeStamp::undef)
      this->t = t ;
    else
      throw std::range_error("TimeStamp value out of range") ;
  }

  TimeStamp::inttime
  TimeStamp::operator=(const TimeStamp::inttime t) {
    if ((t>=TimeStamp::min && t<=TimeStamp::max) || t==TimeStamp::undef)
      this->t = t ;
    else
      throw std::range_error("TimeStamp value out of range") ;
    return t ;
  }

  TimeStamp::inttime 
  TimeStamp::now(void) {
#ifdef __APPLE__
    static uint64_t epoch, epoch_mach ;
    static mach_timebase_info_data_t sTimebaseInfo ;
    if (sTimebaseInfo.denom == 0) {
      struct timeval stamp ;
      gettimeofday(&stamp, 0) ;
      epoch_mach = mach_absolute_time() ;
      epoch = (inttime)stamp.tv_sec*one_second +  (inttime)stamp.tv_usec*one_microsecond ;
      (void) mach_timebase_info(&sTimebaseInfo) ;
      return epoch ;
    }
    uint64_t elapsedNano = (mach_absolute_time()-epoch_mach) * sTimebaseInfo.numer / sTimebaseInfo.denom ;
    return epoch + elapsedNano*TimeStamp::one_nanosecond ;
#elif _WIN32
    static uint64_t epoch ;
    static double freq;
    static LARGE_INTEGER li;
    static bool once = true;
    LARGE_INTEGER li2;
    if (once) {
        once = false;
         if(!QueryPerformanceFrequency(&li))
             std::runtime_error("QueryPerformanceFrequency failed!");
         freq = double(li.QuadPart);
         struct timeval stamp ;
         gettimeofday(&stamp, 0) ;
         epoch = (inttime)stamp.tv_sec*one_second +  (inttime)stamp.tv_usec*one_microsecond ;
         QueryPerformanceCounter(&li);
         return epoch ;
    }
    QueryPerformanceCounter(&li2);
    uint64_t elapsedMicros = double(li2.QuadPart-li.QuadPart)/freq*1000000;
    return epoch + elapsedMicros*one_microsecond ;
#else
    struct timeval stamp ;
    gettimeofday(&stamp, 0) ;
    return (inttime)stamp.tv_sec*one_second +  (inttime)stamp.tv_usec*one_microsecond ;
#endif
  }


  TimeStamp::inttime 
  TimeStamp::string2int(std::string s) {
    bool isPureInt = true ;
    inttime t = 0 ;
    for (unsigned int i=0; i<s.size(); ++i) {
      if (s[i]<'0' || s[i]>'9') {isPureInt = false ; break ;}
      t = 10*t + (s[i]-'0') ;
    }
    if (isPureInt) return t ;

    // 2006-02-25T11:59:12.113449
    // 2006_02_25_11_59_12_113449
    // 012345678901234567890
    //           1         2

    inttime frac = 0 ;
    struct tm aTm ;
    memset(&aTm, 0, sizeof(aTm)) ;
    /*int nbitems = */ sscanf(s.c_str(),"%4d-%2d-%2dT%2d:%2d:%2d.%9ldZ",
			      &aTm.tm_year,&aTm.tm_mon,&aTm.tm_mday,
			      &aTm.tm_hour,&aTm.tm_min,&aTm.tm_sec,
			      &frac) ;
#if 0
    std::cerr << nbitems << " items parsed: "
	      << std::setfill('0') << std::setw(4) << aTm.tm_year
	      << "-" << std::setfill('0') << std::setw(2) << aTm.tm_mon
	      << "-" << std::setfill('0') << std::setw(2) << aTm.tm_mday
	      << "T" << std::setfill('0') << std::setw(2) << aTm.tm_hour 
	      << ":" << std::setfill('0') << std::setw(2) << aTm.tm_min 
	      << ":" << std::setfill('0') << std::setw(2) << aTm.tm_sec 
	      << "." << std::setfill('0') << std::setw(9) << frac
	      << std::endl ;
#endif
    aTm.tm_year -= 1900 ;
    aTm.tm_mon-- ;

    return (inttime)timegm(&aTm)*one_second + frac ;
  }

  std::string
  TimeStamp::int2string(TimeStamp::inttime t) {
    if(t==TimeStamp::undef) return "undef" ;

    time_t sec = (time_t)(t/one_second) ;
    inttime frac = t - (inttime)(sec)*one_second ;
    if (t<0 && frac!=0) {
      sec -- ;
      frac += one_second ;
    }

    struct tm *tm = gmtime(&sec) ;
    if (tm==NULL) return "undef[NULL]" ;

    std::stringstream result ;
    result << std::setfill('0') << std::setw(4) << 1900+tm->tm_year
	   << "-" << std::setfill('0') << std::setw(2) << 1+tm->tm_mon
	   << "-" << std::setfill('0') << std::setw(2) << tm->tm_mday
	   << "T" << std::setfill('0') << std::setw(2) << tm->tm_hour
	   << ":" << std::setfill('0') << std::setw(2) << tm->tm_min
	   << ":" << std::setfill('0') << std::setw(2) << tm->tm_sec
	   << "." << std::setfill('0') << std::setw(9) << frac
	   << "Z" ;
    return result.str() ;
  }

  TimeStamp::inttime
  TimeStamp::ext2int(int year, int month, int day, 
		     int hour, int min, int sec, int frac) {
    struct tm tm ;
    tm.tm_sec = sec ;
    tm.tm_min = min ;
    tm.tm_hour = hour ;
    tm.tm_mday = day ;
    tm.tm_mon = month-1 ;
    tm.tm_year = year - 1900 ;
    tm.tm_wday = tm.tm_yday = 0 ; // ignored anyway
    tm.tm_isdst = 0 ; // ignored anyway
#ifndef _WIN32
    tm.tm_gmtoff = 0; // ignored anyway
    tm.tm_zone = 0 ;
#endif

    return (inttime)timegm(&tm)*one_second + (inttime)frac ;
  }

  void
  TimeStamp::int2ext(TimeStamp::inttime t, 
		     int *year, int *month, int *day,
		     int *hour, int *min, int *sec, int *frac,
		     int *wday, int *yday,
		     bool gmt) {
    time_t seconds = (time_t)(t/one_second) ;
    if (frac) *frac = (int)(t-seconds*one_second) ;

    struct tm *tm = gmt ? gmtime(&seconds) : localtime(&seconds) ;
    if (year) *year = 1900+tm->tm_year ;
    if (month) *month = 1+tm->tm_mon ;
    if (day) *day = tm->tm_mday ;
    if (hour) *hour = tm->tm_hour ;
    if (min) *min = tm->tm_min ;
    if (sec) *sec = tm->tm_sec ;
    if (wday) *wday = tm->tm_wday ;
    if (yday) *yday = tm->tm_yday ;
  }

  TimeStamp::inttime
  TimeStamp::getLocalUTCOffset(void) {
#ifdef _WIN32
    std::cerr << "Warning: TimeStamp::getLocalUTCOffset is not implemented" << std::endl ;
    return 0 ;
#else
    time_t seconds = (time_t)(now()/one_second) ;
    struct tm *tm = localtime(&seconds) ;
    return tm->tm_gmtoff*one_second ;
#endif
  }

}
