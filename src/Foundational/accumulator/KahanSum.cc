#include <stdlib.h>

#include "accumulator.h"


KahanSum::KahanSum ()
{
  _c = 0.0;
  _sum = 0.0;

  return;
}

KahanSum &
KahanSum::operator = (double d)
{
  _sum = d;
  _c = 0.0;

  return *this;
}

KahanSum &
KahanSum::operator += (double d)
{
  double y = d - _c;
  double t = _sum + y;
  _c = (t - _sum) - y;
  _sum = t;

  return *this;
}

KahanSum &
KahanSum::operator += (const KahanSum & rhs)
{
  _c += rhs._c;

  return this->operator+=(rhs._sum);
}
