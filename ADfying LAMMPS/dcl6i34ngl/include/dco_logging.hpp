/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

// This file implements logging capability.
#pragma once

#include "dco_configuration.hpp"

namespace dco {

//* this is used as enumeration class for a specific logger; this
//* class defines the global logger
struct dco_global_logging {
  //** this function return identification string (used in logger
  //** output)
  static std::string id() { return "global"; }
};

//* returns current time as string
inline std::string get_current_time();

//* definition of the logger level. Max level can also be set with
//* DCO_LOG_MAX_LEVEL at compile time.
enum log_level_enum {logDISABLED = -1, logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4};

//* this class handles the formating for the logger
template <typename log_env>
  class logging {
  log_level_enum _level;
  std::string _file;
  int _line;
public:
  //** logging can be used with and without file / line numbering in output
  DCO_ALWAYS_INLINE logging(log_level_enum level=logERROR, std::string file="", int line=-1) : _level(level), _file(file), _line(line) {}

  //** the destructor actually writes to logger stream
  DCO_ALWAYS_INLINE ~logging() {
    if (_line == -1)
      os << std::endl;
    else
      os << " : " << _file << ":" << _line << std::endl;

    //*** this is responsible for enabling the compiler to completely
    //*** optimize logging away at compilation time, if
    //*** DCO_LOG_MAX_LEVEL < 0
    if (DCO_LOG_MAX_LEVEL >= 0 &&
        _level <= DCO_LOG_MAX_LEVEL &&
        _level <= level() &&
        stream()) {
      output(os.str());
    }
  }

  //** write default information to stream and return; in addition,
  //** care for correct indentation.
  DCO_ALWAYS_INLINE std::ostringstream &get_stream() {
    os << "dco/c++ [" << log_env::id() << "] - " << get_current_time();
    os << " " << std::setw(7) << to_string(_level) << ": ";
    os << std::string(_level > logINFO ? _level - logINFO : 0, '\t');
    return os;
  }

public:
  //** public routine for setting the level at run time.
  static inline log_level_enum &level() {
    static log_level_enum reporting_level = static_cast<log_level_enum>(DCO_LOG_MAX_LEVEL);
    return reporting_level;
  }

  //** public routine for converting a given level into a string.
  static inline std::string to_string(log_level_enum level) {
    static const char *const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
    return buffer[level];
  }

  
  //** public routine for converting a given string into a level.
  static inline log_level_enum from_string(const std::string &level) {
    if (level == "DEBUG4")  return logDEBUG4;
    if (level == "DEBUG3")  return logDEBUG3;
    if (level == "DEBUG2")  return logDEBUG2;
    if (level == "DEBUG1")  return logDEBUG1;
    if (level == "DEBUG")   return logDEBUG;
    if (level == "INFO")    return logINFO;
    if (level == "WARNING") return logWARNING;
    if (level == "ERROR")   return logERROR;
    logging(logWARNING).get_stream() << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return logINFO;
  }

  //** public interface to set stream of this logger.
  static inline FILE *&stream() {
    static FILE *pStream = stderr;
    return pStream;
  }

protected:
  std::ostringstream os;

private:
  logging(const logging &);
  logging &operator =(const logging &);

  //** flush stream
  static inline void output(const std::string &msg) {
    FILE *pStream = stream();
    if (!pStream)
      return;
    fprintf(pStream, "%s", msg.c_str());
    fflush(pStream);
  }
};

//* declare global logger: dco::logger
typedef logging<dco_global_logging> logger;

//* logging defines. Usage:
//**   DCO_LOGL(level) << ...; // this will print line and file as well
//**   DCO_LOG(level) << ...;  // this will NOT print line and file.
#define DCO_LOGL(level) if (level > DCO_LOG_MAX_LEVEL) {} dco::logger(level, __FILE__, __LINE__).get_stream()
#define DCO_LOG(level) if (level > DCO_LOG_MAX_LEVEL) {} else dco::logger(level).get_stream()

//* implementatino of get_current_time for Windows and Linux
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
inline std::string get_current_time() {
  const size_t SIZE=50;

  time_t ltime;
  time( &ltime );

  wchar_t wchar_result[SIZE];
  _wctime_s( wchar_result, SIZE, &ltime );

  std::wstring wstring_result(wchar_result);
  std::string result(wstring_result.begin(), wstring_result.end());
  result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());

  return result;
}
#else
inline std::string get_current_time() {
  char buffer[11];
  time_t t;
  time(&t);
  tm r;
  strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
  struct timeval tv;
  gettimeofday(&tv, 0);
  char result[100] = {0};
  std::sprintf(result, "%s.%03ld", buffer, static_cast<long>(tv.tv_usec) / 1000);
  return result;
}
#endif //WIN32

//** helper struct to log creation and remove of objects. Can be used
//** by inheriting from this one.
struct object_logger {
  const std::string _name;
  object_logger(std::string name) : _name(name) {
    DCO_LOG(logINFO) << "Creating " << _name;
  }
  ~object_logger() {
    DCO_LOG(logINFO) << "Removing " << _name;
  }
};

} // end namespace dco


