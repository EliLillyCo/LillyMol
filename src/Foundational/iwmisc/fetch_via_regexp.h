#if ! defined(FETCH_VIA_REGEXP_H)
#define FETCH_VIA_REGEXP_H

#include <memory>
#include <ostream>

#include "re2/re2.h"

#include "Foundational/cmdline/cmdline.h"
#include "misc.h"

class Fetch_via_Regexp
{
  private:
    int _n;
    std::unique_ptr<re2::RE2> * _rx;
    int * _matches_found;         // not thread safe

  public:
    Fetch_via_Regexp();
    ~Fetch_via_Regexp();

    int build (Command_Line & cl, const char flag);

    int active() const { return _n;}

    int report (std::ostream & output) const;

    int matches (const IWString & s);
};

#endif
