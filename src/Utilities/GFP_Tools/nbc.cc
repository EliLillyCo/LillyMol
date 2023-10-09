// Naieve Bayesian Classifier
#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"

#include "nbc.h"

using std::cerr;
using std::endl;

Naive_Bayesian_Classifier::Naive_Bayesian_Classifier()
{
  return;
}

int
Naive_Bayesian_Classifier::build(const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Naive_Bayesian_Classifier::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input);
}

int
Naive_Bayesian_Classifier::build(iwstring_data_source & input)
{
  assert(_fsbw.arrays_allocated());

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;

    if (! _build (buffer))
    {
      cerr << "Invalid weight specification '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

/*
  Unfortunately, these models can come from multiple sources.
  gfp_profile_activity_by_bits produces output like

SFP 0 bit 0 N = 80 in class 1 22 NB 0.746014 in class -1 58
FP 0 bit 2 N = 142 in class 1 54 NB 1.03873 in class -1 88

  Whereas naive_bayesian_classifier_group produces output like

FP 0 bit 0 class FOO NB 1.003
SFP 0 bit 0 class FOO NB 0.649
*/

static std::unique_ptr<RE2> build1_rx = std::make_unique<RE2>("(SFP|FP) ([0-9]+) bit ([0-9]+) N = [0-9]+ in class ");
static std::unique_ptr<RE2> build2_rx = std::make_unique<RE2>("(SFP|FP) ([0-9]+) bit ([0-9]+) NB ([^ ]+)");

int
Naive_Bayesian_Classifier::_build(const const_IWSubstring & buffer)
{
  re2::StringPiece tmp(buffer.data(), buffer.length());

  std::string fptype;
  uint32_t fpnum;
  uint32_t bitnum;
  if (RE2::PartialMatch(tmp, *build1_rx, &fptype, &fpnum, &bitnum)) {
    return _build_from_gfp_profile_activity_by_bits(buffer, fptype, fpnum, bitnum);
  }

  float nb;
  if (RE2::PartialMatch(tmp, *build2_rx, &fptype, &fpnum, &bitnum, &nb)) {
    return _build_from_naive_bayesian_classifier_group(buffer, fptype, fpnum, bitnum, nb);
  }

  cerr << "Naive_Bayesian_Classifier::_build:invalid input\n";
  return 0;
}

int
Naive_Bayesian_Classifier::_build_from_gfp_profile_activity_by_bits(const const_IWSubstring& buffer,
    const std::string& fptype,
    uint32_t fpnum,
    uint32_t bitnum)
{

  if (fpnum < 0 || fpnum > (NBC_FIXED_ARRAY_SIZE - 1))
  {
    cerr << "Naive_Bayesian_Classifier::_build:fixed array size overflow, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    return 0;
  }


  int i = 0;
  const_IWSubstring token;
  for (int j = 0; j < 7; j++)
  {
    buffer.nextword(token, i);
  }

  int prev_was_in = 0;
  int prev_was_class = 0;

  while (buffer.nextword(token, i))
  {
//  cerr << "Token is '" << token << "'\n";
    if ("in" == token)
    {
      prev_was_in = 1;
      continue;
    }

    if ("class" == token)
    {
      if (! prev_was_in)
      {
        prev_was_in = 0;
        continue;
      }

      prev_was_class = 1;
      continue;
    }

    if (prev_was_class)
    {
      _zclass = token;
      prev_was_class = 0;
      prev_was_in = 0;
      continue;
    }

    if ("NB" != token)
      continue;

    buffer.nextword(token, i);

    double w;
    if (! token.numeric_value(w) || w <= 0.0)
    {
      cerr << "Naive_Bayesian_Classifier::_build:invalid weight '" << token << "'\n";
      return 0;
    }

    if ("SFP" == fptype)
      _fsbw._sparsefp[fpnum].set_weight(bitnum, w);
    else if ("FP" == fptype)
      _fsbw._fp[fpnum].set_weight(bitnum, w);
    else
      abort();

//  cerr << "ndx " << ndx << " sizes " << _sparsefp[ndx].size() << " and " << _fp[ndx].size() << endl;
    return 1;
  }

  cerr << "Naive_Bayesian_Classifier::_build:never found Bayesian directive\n";
  return 0;
}

int
Naive_Bayesian_Classifier::_build_from_naive_bayesian_classifier_group(const const_IWSubstring& buffer,
        const std::string& fptype,
        uint32_t fpnum,
        uint32_t bitnum,
        float nb)
{
  if (fpnum < 0 || fpnum > (NBC_FIXED_ARRAY_SIZE - 1)) {
    cerr << "Naive_Bayesian_Classifier::_build:fixed array size overflow, contact LillyMol on github (https://github.com/EliLillyCo/LillyMol)\n";
    return 0;
  }

  if ("SFP" == fptype)
    _fsbw._sparsefp[fpnum].set_weight(bitnum, nb);
  else if ("FP" == fptype)
    _fsbw._fp[fpnum].set_weight(bitnum, nb);
  else {
    cerr << "Naive_Bayesian_Classifier::_build_from_naive_bayesian_classifier_group:huh?\n";
    abort();
  }

  return 1;
}

int
Naive_Bayesian_Classifier::_classify(int ndx,
                                     const IWDYFP & fp,
                                     double & rc) const
{
  const Bit_and_Weight & bw = _fsbw._fp[ndx];

  if (0 == bw.size())
    return 1;

  int nb = fp.nbits();

  for (int i = 0; i < nb; i++)
  {
    if (! fp.is_set(i))
      continue;

    bw.increment_weight(i, rc);
  }

  return 1;
}

int
Naive_Bayesian_Classifier::_classify(int ndx,
                                     const Sparse_Fingerprint & sfp,
                                     double & rc) const
{
  const Bit_and_Weight & bw = _fsbw._sparsefp[ndx];

//cerr << "ndx " << ndx << " bw contains " << bw.size() << " items\n";
  if (0 == bw.size())
    return 1;

  unsigned int b;
  int i = 0;
  int notused;

  while (sfp.next_bit_set (i, b, notused))
  {
    bw.increment_weight(b, rc);
  }

  return 1;
}

double
Naive_Bayesian_Classifier::classify(const IW_General_Fingerprint & fp) const
{
  int n = fp.nfingerprints();

  double rc = 0.0;

  for (int i = 0; i < n; i++)
  {
    _classify(i, fp[i], rc);
  }

  n = fp.number_sparse_fingerprints();

  for (int i = 0; i < n; i++)
  {
    _classify(i, fp.sparse_fingerprint(i), rc);
  }

  return rc;
}

int
Naive_Bayesian_Classifier::set_sfp_weight(int f,
                                          unsigned int b,
                                          double w)
{
  assert (f >= 0 && f < NBC_FIXED_ARRAY_SIZE);

  _fsbw._sparsefp[f].set_weight(b, w);

  return 1;
}

int
Naive_Bayesian_Classifier::set_fp_weight(int f,
                                          unsigned int b,
                                          double w)
{
  assert (f >= 0 && f < NBC_FIXED_ARRAY_SIZE);

  _fsbw._fp[f].set_weight(b, w);

  return 1;
}
