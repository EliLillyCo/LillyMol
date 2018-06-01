#ifndef IW_TIME_H
#define IW_TIME_H

#include <stdlib.h>
#include <time.h>

class IW_Time
{
  private:
    time_t _time;

  public:
    IW_Time() { _time = ::time(NULL);}

    time_t time() const { return _time;}

};

extern std::ostream &
operator << (std::ostream & os, const IW_Time & iwt)
{
  os << iwt.time();

  return os;
}

extern IWString &
operator << (IWString & os, const IW_Time & iwt)
{
  os << static_cast<unsigned int>(iwt.time());

  return os;
}

extern time_t
operator - (const IW_Time & t1, const IW_Time & t2)
{
  return t1.time() - t2.time();
}

#endif
