#include <cmath>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/cmdline/cmdline.h"
#include "fb_bits_and_weights.h"

FB_Bits_and_Weights_Fixed_Width::FB_Bits_and_Weights_Fixed_Width()
{
  _nbits = -1;
  _wmatch = nullptr;
  _wnmatch = nullptr;

  _include_non_matching_bits_in_similarity = 1;

  return;
}

FB_Bits_and_Weights_Fixed_Width::~FB_Bits_and_Weights_Fixed_Width()
{
  if (nullptr != _wmatch)
  {
    delete [] _wmatch;
    delete [] _wnmatch;
  }

  _nbits = -77;

  return;
}

int
FB_Bits_and_Weights_Fixed_Width::build_from_file (const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::build_from_file:cannot open '" << fname << "'\n";
    return 0;
  }

  return build_from_file (input);
}

int
FB_Bits_and_Weights_Fixed_Width::construct_from_command_line (Command_Line & cl,
                                                char flag,
                                                int verbose)
{
  int i = 0;
  const_IWSubstring f;

  IWString fname;
  while (cl.value(flag, f, i++))
  {
    if ("all" == f)
    {
      _include_non_matching_bits_in_similarity = 1;
      if (verbose)
        cerr << "Tanimoto measure will include weights for non-matching bits\n";
    }
    else if ("matches" == f)
    {
      _include_non_matching_bits_in_similarity = 0;
      if (verbose)
        cerr << "Tanimoto measure will only include weights for matching bits\n";
    }
    else if ("help" == f)
    {
      cerr << " - " << flag << " all        Tanimoto includes weights for matches and non matches\n";
      cerr << " - " << flag << " matches    Tanimoto only includes weights from matched bits \n";
      cerr << " - " << flag << " fname      file name with weights\n";
      cerr << " - " << flag << " help       this message\n";
      exit (0);
    }
    else if (fname.length())
    {
      cerr << "Can only specify one file for weighted fingerprints\n";
      return 0;
    }
    else
      fname = f;
  }

  if (0 == fname.length())
  {
    cerr << "Must specify bit weight file name via the -" << flag << " option\n";
    return 0;
  }

  return build_from_file(fname.null_terminated_chars());
}

int
FB_Bits_and_Weights_Fixed_Width::build_from_file (iwstring_data_source & input)
{
  input.set_translate_tabs(1);
  input.set_dos(1);

  _nbits = input.records_remaining() - 1;

  if (_nbits < 1)
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::build_from_file:empty file\n";
    return 0;
  }

  _wmatch = new double[_nbits];
  _wnmatch = new double[_nbits];

  const_IWSubstring buffer;

  input.next_record(buffer);   // skip header record

  for (int ndx = 0; input.next_record(buffer); ndx++)
  {
    if (! _parse_weight_record (buffer, ndx))
    {
      cerr << "FB_Bits_and_Weights_Fixed_Width::build_from_file:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

#ifdef ECHO_WEIGHTS
  for (int i = 0; i < _nbits; i++)
  {
    cerr << " bit " << i << " weights " << _wmatch[i] << " and " << _wnmatch[i] << endl;
  }
#endif

  return _nbits;
}

int
FB_Bits_and_Weights_Fixed_Width::_parse_weight_record (const const_IWSubstring & buffer,
                                                int ndx)
{
  if (buffer.nwords() < 3)
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::_parse_weight_record:must be at least three tokens\n";
    return 0;
  }

  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  int b;
  if (! token.numeric_value(b) || b < 0)
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::_parse_weight_record:invalid bit number '" << token << "'\n";
    return 0;
  }

  if (b != ndx)
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::_parse_weight_record:bit number mismatch, expected " << ndx << " got " << b << endl;
    return 0;
  }

  buffer.nextword(token, i);

  if (! token.numeric_value(_wmatch[ndx]) || _wmatch[ndx] < 0.0)
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::_parse_weight_record:invalid weight '" << token << "'\n";
    return 0;
  }

  buffer.nextword(token, i);

  if (! token.numeric_value(_wnmatch[ndx]) || _wnmatch[ndx] < 0.0)
  {
    cerr << "FB_Bits_and_Weights_Fixed_Width::_parse_weight_record:invalid weight '" << token << "'\n";
    return 0;
  }

  return 1;
}

double
FB_Bits_and_Weights_Fixed_Width::tanimoto (const int * fp1,
                                           const int * fp2) const
{
  assert (_nbits > 0);

  if (_include_non_matching_bits_in_similarity)
    return _tanimoto_include_non_matching_bits(fp1, fp2);
  else
    return _tanimoto_only_consider_matching_bits(fp1, fp2);
}

double
FB_Bits_and_Weights_Fixed_Width::_tanimoto_include_non_matching_bits (const int * fp1,
                                        const int * fp2) const
{
  double sumk = 0.0;
  double suml = 0.0;
  double summ = 0.0;

#ifdef DEBUG_FBEWT_TANIMOTO
  cerr << "Processing " << _nbits << " bits\n";
#endif

  for (int ndx = 0; ndx < _nbits; ndx++)
  {
#ifdef DEBUG_FBEWT_TANIMOTO
    cerr << "Bit " << ndx << " weights " << _wmatch[ndx] << " and " << _wnmatch[ndx] << endl;
#endif

    double h = fp1[ndx] * _wmatch[ndx] + (1 - fp1[ndx]) * _wnmatch[ndx];
    double i = fp2[ndx] * _wmatch[ndx] + (1 - fp2[ndx]) * _wnmatch[ndx];
    double j = (1.0 - fabs(static_cast<double>(fp1[ndx] - fp2[ndx]))) * h;

    assert (h >= 0.0);
    assert (i >= 0.0);
    assert (j >= 0.0);

    double maxfg;
    if (fp1[ndx] || fp2[ndx])
      maxfg = 1.0;
    else
      maxfg = 0.0;

    double k = maxfg * h;
    double l = maxfg * i;
    double m = maxfg * j;

#ifdef DEBUG_FBEWT_TANIMOTO
    cerr << " k = " << k << " l = " << j << " m = " << m << endl;
#endif

    sumk += k;
    suml += l;
    summ += m;
  }

#ifdef DEBUG_FBEWT_TANIMOTO
  cerr << "summ " << summ << " sumk " << sumk << " suml " << suml << " result " << (summ / (sumk + suml - summ)) << endl;
#endif

  return summ / (sumk + suml - summ);
}

static double
sum_non_zero_items (const double * w,
                    const int * consider,
                    int n)
{
  double rc = 0.0;

  for (int i = 0; i < n; i++)
  {
    if (1 == consider[i])
      rc += w[i];
  }

  return rc;
}

double
FB_Bits_and_Weights_Fixed_Width::_tanimoto_only_consider_matching_bits (const int * fp1,
                                        const int * fp2) const
{
  double o = sum_non_zero_items(_wmatch, fp1, _nbits);
  double p = sum_non_zero_items(_wmatch, fp2, _nbits);
  double summ = 0.0;

  for (int ndx = 0; ndx < _nbits; ndx++)
  {
    if (fp1[ndx] && fp2[ndx])
      summ += _wmatch[ndx];
  }

  return summ / (o + p - summ);
}
