#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "gfp.h"
#include "sparse_collection.h"

using std::cerr;
using std::endl;

static int singleton_threshold = 0;

Sparse_Fingerprint_Collection_Profile::Sparse_Fingerprint_Collection_Profile()
{
}

int
Sparse_Fingerprint_Collection_Profile::report(std::ostream & output) const
{
  output << "Sparse fingerprint set contains " << _xref.size() << " bits\n";
  if (_nset.n())
    output << "Profiled " << _nset.n() << " fingerprints, between " << _nset.minval() << " and " << _nset.maxval() << " bits set, ave " << _nset.average_if_available_minval_if_not() << endl;

  return output.good();
}

int
Set_of_Sparse_Fingerprint_Collection_Profile::report(std::ostream & output) const
{
  output << "Sparse fingerprint summary on " << _number_sparse_fingerprints << " sparse fingerprints\n";

  for (int i = 0; i < _number_sparse_fingerprints ; i++)
  {
    _sfcp[i].report(output);
  }

  return output.good();
}

Set_of_Sparse_Fingerprint_Collection_Profile::Set_of_Sparse_Fingerprint_Collection_Profile()
{
  _number_sparse_fingerprints = 0;

  _sfcp = nullptr;

  return;
}

Set_of_Sparse_Fingerprint_Collection_Profile::~Set_of_Sparse_Fingerprint_Collection_Profile()
{
  if (nullptr != _sfcp)
    delete [] _sfcp;

  _number_sparse_fingerprints = -2;

  return;
}

int
Set_of_Sparse_Fingerprint_Collection_Profile::resize(int s)
{
  assert (s > 0);

  if (nullptr != _sfcp)
    delete [] _sfcp;

  _sfcp = new Sparse_Fingerprint_Collection_Profile[s];

  if (nullptr == _sfcp)
  {
    cerr << "Set_of_Sparse_Fingerprint_Collection_Profile::resize:cannot allocate " << s << " profiles\n";
    return 0;
  }

  _number_sparse_fingerprints = s;

  return 1;
}

int
Set_of_Sparse_Fingerprint_Collection_Profile::build_profile(const IW_General_Fingerprint & fp)
{
  if (0 == _number_sparse_fingerprints)
  {
    int n = fp.number_sparse_fingerprints();
    if (0 == n)
    {
      cerr << "Sparse_Fingerprint_Collection_Profile::build_profile:no sparse fingerprints\n";
      return 0;
    }

    if (! resize(n))
      return 0;
  }

  for (int i = 0; i < _number_sparse_fingerprints; i++)
  {
    const Sparse_Fingerprint & fpi = fp.sparse_fingerprint(i);

    _sfcp[i].build_profile(fpi);
  }

  return 1;
}

int
Sparse_Fingerprint_Collection_Profile::build_profile(const Sparse_Fingerprint & fp)
{
  _nset.extra(fp.nset());

  int rc = 0;    // we return the number of new bits we detect

  int i = 0;
  unsigned int b;
  int c;

  while (fp.next_bit_set(i, b, c))
  {
    xref_t::const_iterator f = _xref.find(b);

    if (f != _xref.end())
    {
      unsigned int j = (*f).second;
      _count[j]++;
    }
    else
    {
      unsigned int s = _xref.size();
      _xref[b] = s;
      _count[s]++;
      rc++;
    }
  }

  return 1;
}

int
Set_of_Sparse_Fingerprint_Collection_Profile::finished_profiling(int verbose)
{
  if (0 != singleton_threshold)
  {
    for (int i = 0; i < _number_sparse_fingerprints; i++)
    {
      _sfcp[i].remove_singletons(singleton_threshold, verbose);
    }
  }
  else if (verbose)
  {
    for (int i = 0; i < _number_sparse_fingerprints; i++)
    {
      _sfcp[i].report(cerr);
    }
  }

  return 1;
}

int
Sparse_Fingerprint_Collection_Profile::remove_singletons(int threshold,
                                                         int verbose)
{
  resizable_array<unsigned int> bits_to_remove;

  bits_to_remove.resize(_xref.size());

  for (xref_t::const_iterator i = _xref.begin(); i != _xref.end(); ++i)
  {
    unsigned int j = (*i).second;    // the index in the _count array

//  cerr << _count[j] << " instances of bit " << (*i).first << endl;

    if (_count[j] <= threshold)
      bits_to_remove.add((*i).first);
  }

  int rc = bits_to_remove.number_elements();

  if (verbose)
    cerr << "Will remove " << rc << " of " << _xref.size() << " bits less than " << threshold << endl;

  if (0 == rc)
    return 0;

  for (int i = 0; i < rc; i++)
  {
    unsigned int b = bits_to_remove[i];

    _xref.erase(b);
  }

// Renumber the bits

  unsigned int * newcount = new unsigned int[_xref.size()]; std::unique_ptr<unsigned int[]> free_newcount(newcount);

  unsigned int b = 0;

  for (xref_t::iterator i = _xref.begin(); i != _xref.end(); ++i)
  {
    unsigned oldb = (*i).second;

    unsigned int oldc = _count[oldb];

    (*i).second = b;
    newcount[b] = oldc;
    b++;
  }

  _count.resize(0);        // discard existing data
  _count.resize(_xref.size());

  for (unsigned int i = 0; i < _xref.size(); i++)
  {
    _count.add(newcount[i]);
  }

  return rc;
}

int
Sparse_Fingerprint_Collection_Profile::convert_to_fixed_width(const Sparse_Fingerprint & fpfrom,
                                                              IWDYFP & fpto,
                                                              int & nextra) const
{
  nextra = 0;

  if (static_cast<unsigned int>(fpto.nbits()) != _xref.size())
    fpto.allocate_space_for_bits(_xref.size());

  int i = 0;
  unsigned int b;
  int c;

  while (fpfrom.next_bit_set(i, b, c))
  {
    xref_t::const_iterator f = _xref.find(b);

    if (f == _xref.end())
      nextra++;
    else
      fpto.set((*f).second);
  }

  return fpto.compute_nset();
}

//#define DEBUG_CONVERT_TO_FIXED_WIDTH

int
Sparse_Fingerprint_Collection_Profile::convert_to_fixed_width(const Sparse_Fingerprint & fpfrom,
                                        Fixed_Size_Counted_Fingerprint_uchar & fpto,
                                        int & extra_bits,
                                        int & extra_count) const
{
  extra_bits = 0;
  extra_count = 0;

  fpto.resize(_xref.size());

  int i = 0;
  unsigned int b;
  int c;

#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
  cerr << "Converting sparse fingerprint with " << fpfrom.nbits() << " bits\n";
#endif

  while (fpfrom.next_bit_set(i, b, c))
  {
    xref_t::const_iterator f = _xref.find(b);

    if (f == _xref.end())
    {
      extra_bits++;
      extra_count += c;
#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
      cerr << "Bit " << b << " not in xref\n";
#endif
    }
    else
    {
      int j = (*f).second;
#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
      cerr << "Bit " << b << " becomes bit " << j << " count " << c << endl;
#endif
      fpto.set(j, c);
    }
  }

#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
  cerr << extra_bits << " extra_bits, and " << extra_count << " extra count\n";
#endif

  return fpto.compute_nset();
}

int
Sparse_Fingerprint_Collection_Profile::convert_to_fixed_width(const Sparse_Fingerprint & fpfrom,
                                        Fixed_Size_Counted_Fingerprint_uint & fpto,
                                        int & extra_bits,
                                        int & extra_count) const
{
  extra_bits = 0;
  extra_count = 0;

  fpto.resize(_xref.size());

  int i = 0;
  unsigned int b;
  int c;

#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
  cerr << "Converting sparse fingerprint with " << fpfrom.nbits() << " bits\n";
#endif

  while (fpfrom.next_bit_set(i, b, c))
  {
    xref_t::const_iterator f = _xref.find(b);

    if (f == _xref.end())
    {
      extra_bits++;
      extra_count += c;
#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
      cerr << "Bit " << b << " not in xref\n";
#endif
    }
    else
    {
      int j = (*f).second;
#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
      cerr << "Bit " << b << " becomes bit " << j << " count " << c << endl;
#endif
      fpto.set(j, c);
    }
  }

#ifdef DEBUG_CONVERT_TO_FIXED_WIDTH
  cerr << extra_bits << " extra_bits, and " << extra_count << " extra count\n";
#endif

  return fpto.compute_nset();
}


static int
display_standard_sparse_to_dense_fingerprint_options(char flag,
                                                     std::ostream & os)
{
  os << " -" << flag << " bits     convert fixed width fingerprints to bits\n";
  os << " -" << flag << " count    convert to fixed width counted fingerprints\n";
  os << " -" << flag << " sgl=<nn> bits hit <nn> or fewer times are singletons and eliminated\n";
  os << " -" << flag << endl;

  return os.good();
}

int
parse_sparse_to_dense_fingerprint_specifications(Command_Line & cl,
                        char flag,
                        int verbose)
{
  int i = 0;
  const_IWSubstring k;

  while (cl.value(flag, k, i++))
  {
    if ("bits" == k)
    {
      set_convert_sparse_fingerprint_to_bits(1);
      if (verbose)
        cerr << "Sparse fingerprints converted to bits\n";
    }
    else if ("count" == k)
    {
      set_convert_sparse_fingerprint_to_fixed_width_counted(1);
      if (verbose)
        cerr << "Sparse fingerprints converted to counted bytes\n";
    }
    else if (k.starts_with("sgl="))
    {
      k.remove_leading_chars(4);
      if (! k.numeric_value(singleton_threshold) || singleton_threshold < 1)
      {
        cerr << "Invalid singleton threshold '" << k << "'\n";
        return 4;
      }
    }
    else if ("help" == k)
    {
      display_standard_sparse_to_dense_fingerprint_options(flag, cerr);
      exit(0);
    }
    else
    {
      cerr << "Unrecognised sparse fingerprint qualifier '" << k << "'\n";
      display_standard_sparse_to_dense_fingerprint_options(flag, cerr);
      return 0;
    }
  }

  return 1;
}


#ifdef NOW_IN_HFILE
int
Set_of_Sparse_Fingerprint_Collection_Profile::convert_to_fixed_width(int which_fingerprint,
                                        Sparse_Fingerprint const & fpfrom, 
                                        IWDYFP & fpto,
                                        int & extra_bits) const
{
  return _sfcp[which_fingerprint].convert_to_fixed_width(fpfrom, fpto, extra_bits);
}

int
Set_of_Sparse_Fingerprint_Collection_Profile::convert_to_fixed_width(int which_fingerprint,
                                Sparse_Fingerprint const & fpfrom,
                                Fixed_Size_Counted_Fingerprint_uchar & fpto, 
                                int & extra_bits,
                                int & extra_count) const
{
  return _sfcp[which_fingerprint].convert_to_fixed_width(fpfrom, fpto, extra_bits, extra_count);
}

int
Set_of_Sparse_Fingerprint_Collection_Profile::convert_to_fixed_width(int which_fingerprint,
                                Sparse_Fingerprint const & fpfrom,
                                Fixed_Size_Counted_Fingerprint_uint & fpto, 
                                int & extra_bits,
                                int & extra_count) const
{
  return _sfcp[which_fingerprint].convert_to_fixed_width(fpfrom, fpto, extra_bits, extra_count);
}
#endif
