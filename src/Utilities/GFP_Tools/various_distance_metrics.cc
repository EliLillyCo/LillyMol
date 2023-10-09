#include <stdlib.h>
#include <assert.h>

#include "Foundational/iwstring/iwstring.h"

#include "various_distance_metrics.h"

using std::cerr;

/*
  FVB define the similarity measure as a linear combination of the on-bit Tversky and the off-bit tversky.
  The ratio is determined by these parameters

  r1 + r2 = r3
*/

static similarity_type_t r1 = 2.0;
static similarity_type_t r2 = 1.0;
static similarity_type_t r1r2 = 3.0;

int
set_fvb_ratios (similarity_type_t b1, similarity_type_t b2)
{
  assert (b1 >= static_cast<similarity_type_t> (0.0));
  assert (b2 >= static_cast<similarity_type_t> (0.0));

  r1 = b1;
  r2 = b2;

  r1r2 = r1 + r2;

  cerr << "FVB ratios set to " << r1 << " " << r2 << " sum " << r1r2 << '\n';

  return 1;
}

int 
set_fvb_ratios (const const_IWSubstring & fvb)
{
  const_IWSubstring s1, s2;
  if (! fvb.split (s1, ',', s2))
  {
    cerr << "Invalid FVB ratio specifier '" << fvb << "'\n";
    return 0;
  }

  similarity_type_t f1;
  if (! s1.numeric_value (f1) || f1 < static_cast<similarity_type_t> (0.0))
  {
    cerr << "INvalid FVB ratio '" << s1 << "'\n";
    return 0;
  }

  similarity_type_t f2;
  if (! s2.numeric_value (f2) || f2 < static_cast<similarity_type_t> (0.0))
  {
    cerr << "INvalid FVB ratio '" << s1 << "'\n";
    return 0;
  }

  return set_fvb_ratios (f1, f2);
}

//#define DEBUG_FVB

similarity_type_t
fligner_verducci_blower (int nb,
                         int nset1,
                         int nset2,
                         int n00,
                         int n11)
{
  assert (n00 <= nb);
  assert (n11 <= nb);

  if (n11 == nb)    // all bits in common, they are identical
    return static_cast<similarity_type_t> (1.0);

  if (n00 == nb)    // both must be all zero
    return static_cast<similarity_type_t> (1.0);

  similarity_type_t phat = static_cast<similarity_type_t> (nset1 + nset2) / static_cast<similarity_type_t> (nb + nb);

  similarity_type_t t0 = static_cast<similarity_type_t> (n00) / static_cast<similarity_type_t> (nb - n11);
  similarity_type_t t1 = static_cast<similarity_type_t> (n11) / static_cast<similarity_type_t> (nb - n00);

  similarity_type_t rc = (r1 - phat) / r1r2 * t1 + (r2 + phat) / r1r2 * t0;

#ifdef DEBUG_FVB
  cerr << "t0 " << t0 << " and t1 " << t1 << " phat " << phat << '\n';
  cerr << "Returning " << rc << '\n';
#endif

  if (rc < static_cast<similarity_type_t> (0.0) || rc > static_cast<similarity_type_t> (1.0))
    cerr << "INvalid similarity " << rc << " nb = " << nb << " nset1 " << nset1 << " nset2 " << nset2 << " n00 " << n00 << " n11 " << n11 << '\n';
  assert (rc >= static_cast<similarity_type_t> (0.0) && rc <= static_cast<similarity_type_t> (1.0));

  return rc;
}

