// Experiment with different ways of storing sparse bitvectors.

#include <ctime>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "gfp.h"

using std::cerr;
using std::endl;

extern const unsigned char bic_table[256*256];

namespace {

int verbose = 0;

IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");

void usage(int rc) {
  exit(rc);
}

class AsVector {
  private:
    std::vector<uint32_t> _bits;
    int _nset;

    float _tanimoto(const AsVector& rhs) const;

  public:
    AsVector() {
    }

    int build(const Sparse_Fingerprint& sfp);

    float tanimoto(const AsVector& rhs) const;
};

int
AsVector::build(const Sparse_Fingerprint& sfp) {
  _bits.resize(sfp.nbits() * 2);

  int i = 0;
  unsigned int b;
  int c;
  int ndx = 0;
  _nset = 0;
  while (sfp.next_bit_set(i, b, c)) {
    _bits[ndx] = b;
    _bits[ndx + 1] = c;
    ndx += 2;
//  cerr << " i " << i << " bit " << b << " count " << c << endl;
    _nset += c;
  }
//cerr << " _nset " << _nset << endl;

  return _nset;
}

#ifdef SLOWER_QQ
float
AsVector::_tanimoto(const AsVector& rhs) const
{
  const int n1 = _bits.size();

  const uint32_t * d1 = _bits.data();
  const uint32_t * d2 = rhs._bits.data();

//cerr << "AsVector::_tanimoto: nbits " << n1 << " nset " << _nset << " RHS nb " << (rhs._bits.size() / 2) << " nset " << rhs._nset << '\n';

  int bits_in_common = 0;

  int i1 = 0;
  int i2 = 0;
  for (int i1 = 0; i1 < n1; i1 += 2)
  {
    while (d2[i2] < d1[i1])
      i2 += 2;

    if (d1[i1] == d2[i2]) {
//    cerr << " BIT in common " << d1[i1] << " counts " << d1[i1 + 1] << ' ' << d2[i2 + 1] << '\n';
//    bits_in_common += std::min(d1[i1 + 1], d2[i2 + 1]);
//    bits_in_common += ::bic_table[d1[i1 + 1] * 256 + d2[i2 + 1]];
      if (d1[i1 + 1] <= d2[i2 + 1])
        bits_in_common += d1[i1 + 1];
      else
        bits_in_common += d2[i2 + 1];
      i2 += 2;
    }
  }

//cerr << _nset << " and " << rhs._nset << " bic " << bits_in_common << endl;

  return static_cast<float> (bits_in_common) / static_cast<float>(_nset + rhs._nset - bits_in_common);
}
#endif

float
AsVector::_tanimoto(const AsVector& rhs) const
{
  const int n1 = _bits.size();

  const uint32_t * d1 = _bits.data();
  const uint32_t * d2 = rhs._bits.data();

//cerr << "AsVector::_tanimoto: nbits " << n1 << " nset " << _nset << " RHS nb " << (rhs._bits.size() / 2) << " nset " << rhs._nset << '\n';

  int bits_in_common = 0;

  int i1 = 0;
  int i2 = 0;
  while (i1 < n1) {
    if (d1[i1] < d2[i2]) {
      i1 += 2;
    } else if (d1[i1] > d2[i2]) {
      i2 += 2;
    } else {
//    cerr << " BIT in common " << d1[i1] << " counts " << d1[i1 + 1] << ' ' << d2[i2 + 1] << '\n';
//    bits_in_common += std::min(d1[i1 + 1], d2[i2 + 1]);
//    bits_in_common += ::bic_table[d1[i1 + 1] * 256 + d2[i2 + 1]];
      if (d1[i1 + 1] <= d2[i2 + 1])
        bits_in_common += d1[i1 + 1];
      else
        bits_in_common += d2[i2 + 1];
      i1 += 2;
      i2 += 2;
    }
  }

//cerr << _nset << " and " << rhs._nset << " bic " << bits_in_common << endl;

  return static_cast<float> (bits_in_common) / static_cast<float>(_nset + rhs._nset - bits_in_common);
}

float
AsVector::tanimoto(const AsVector& rhs) const
{
  if (_bits[_bits.size() - 2] < rhs._bits[rhs._bits.size() - 2])
    return _tanimoto(rhs);
  return rhs._tanimoto(*this);
}

float
separate_vectors(Sparse_Fingerprint* fp, const int number_fingerprints)
{
  const auto t0 = std::clock();

  float min_dist = 0.0;
  double tot_dist = 0.0;
  for (int i = 0; i < number_fingerprints; ++i) {
    for (int j = i + 1; j < number_fingerprints; ++j) {
      const float d = fp[i].tanimoto(fp[j]);
      if (d < min_dist)
        min_dist = d;
      tot_dist += d;
    }
  }

  cerr << "separate_vectors: time " << (std::clock() - t0) << " dist " << min_dist << " tpt " << tot_dist << '\n';

  return min_dist;
}

float
one_vector(Sparse_Fingerprint * fp, const int number_fingerprints) {
  AsVector * as_vector = new AsVector[number_fingerprints];
  std::unique_ptr<AsVector[]> free_as_vector(as_vector);

  for (int i = 0; i < number_fingerprints; ++i) {
    as_vector[i].build(fp[i]);
  }

  const auto t0 = std::clock();

  float min_dist = 1.0f;
  double tot_dist = 0.0;

  for (int i = 0; i < number_fingerprints; ++i) {
    for (int j = i + 1; j < number_fingerprints; ++j) {
      const float d = as_vector[i].tanimoto(as_vector[j]);
      if (d < min_dist)
        min_dist = d;
      tot_dist += d;
    }
  }

  cerr << "one_vector: time " << (std::clock() - t0) << " dist " << min_dist << " tpt dist " << tot_dist << '\n';

  return min_dist;
}

int
read_fingerprints(const char * fname,
                  IW_General_Fingerprint* & pool)
{
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  const int number_fingerprints = input.grep("^PCN<");
  if (number_fingerprints == 0) {
    cerr << "No fingerprints in '" << fname << "'\n";
    return 0;
  }

  pool = new IW_General_Fingerprint[number_fingerprints];

  IW_TDT tdt;
  for (int poolptr = 0; tdt.next(input); poolptr++)
  {
    int fatal;
    if (! pool[poolptr].construct_from_tdt(tdt, fatal))
    {
      if (fatal)
      {
        cerr << "Cannot parse tdt " << tdt;
        return 0;
      }

      continue;
    }
  }

  return number_fingerprints;
}

int
fpobj_spread (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vs:n:I:i:A:r:p:t:F:P:W:Q:O:N:V:b:S:M:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    return 1;
  }

  verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Must specify gfp file\n";
    return 1;
  }

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise GFP options\n";
      usage(23);
    }
  }
  else if (! initialise_fingerprints(cl[0], verbose))
  {
    cerr << "Cannot initialise fingerprints from '" << cl[0] << "'\n";
    return 11;
  }

  IW_General_Fingerprint  * pool = nullptr;

  const int number_fingerprints = read_fingerprints(cl[0], pool);

  if (verbose)
    cerr << "Read " << number_fingerprints << " from '" << cl[0] << "'\n";

  if (pool[0].number_sparse_fingerprints() != 1)
  {
    cerr << "Must have 1 sparse fingerprint in input\n";
    return 1;
  }

  Sparse_Fingerprint * fp = new Sparse_Fingerprint[number_fingerprints];
  std::unique_ptr<Sparse_Fingerprint[]> free_fp(fp);
  for (int i = 0; i < number_fingerprints; ++i) {
    fp[i] = pool[i].sparse_fingerprint(0);
  }

  delete [] pool;

  separate_vectors(fp, number_fingerprints);

  one_vector(fp, number_fingerprints);

  return 0;
}

}  // namespace

int
main (int argc, char ** argv)
{
  int rc = ::fpobj_spread(argc, argv);

  return rc;
}
