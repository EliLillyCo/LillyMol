/*
  Echo's fingerpints from a fingerprint file
*/

#include <stdlib.h>
#include <memory>
using namespace std;

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "iwdigits.h"
#include "iw_stl_hash_map.h"
#include "accumulator.h"
#include "misc.h"

#include "dyfp.h"
#include "sparsefp.h"

const char * prog_name = NULL;

static int verbose = 0;

static IWString fptag;

static IWString identifier_tag ("PCN<");

static unsigned int fingerprints_read = 0;

static int non_colliding_form = 1;

static int input_is_hex_encoded = 0;

static unsigned int nbits = 0;    // for fixed width fingerprints

static const float zero = static_cast<float> (0.0);

static float lower_population_cutoff = zero;

static unsigned int lower_support_cutoff_number = 0;

static float upper_population_cutoff = zero;

static unsigned int upper_support_cutoff_number = 0;

static int min_non_zero_values_needed_on_each_row = 0;

static IWString_and_File_Descriptor stream_for_not_enough_non_zero_values;

static int write_identifiers_suppressed_by_dash_y = 0;

static int records_discarded_for_not_enough_non_zero_values = 0;

static int gsub_multi_token_names = 0;

static int produce_sparse_form_for_nikil = 0;

static IWString descriptor_prefix (" F");

static IWDigits iwdigits;

static int interpret_fixed_width_as_numeric = 0;

static int fold_to_constant_width = 0;

/*
  By default, we create arbitrary numbers
*/

static int write_actual_bit_numbers = 0;

/*
  It can be convenient to sort the bits so the most common occur as the first descriptors
*/

static int sort_by_frequency = 0;

static int write_zero_count_as_negative_one = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Converts a fingerprint set from a gfp file to a descriptor file\n";
  cerr << " -F <tag>       specify tag of the fingerprint to be processed\n";
  cerr << " -f             fixed width fingerprints (use the underlying programme instead)\n";
  cerr << " -h             fixed width hex encoded fingerprints\n";
//cerr << " -m             molecular properties\n";
  cerr << " -x <fraction>  discard bits that hit less than <fraction> of the time\n";
  cerr << " -n <number>    discard bits hit in less than <number> of the fingerprints\n";
  cerr << " -X <fraction>  discard bits that hit more than <fraction> of the time\n";
  cerr << " -N <number>    discard bits hit in more than <number> of the fingerprints\n";
  cerr << "                    use either -x or -n, and either -X or -N\n";
  cerr << " -y <number>    discard records unless they have at least <number> non zero values\n";
  cerr << " -s             sort by frequency\n";
  cerr << " -P <string>    prefix for descriptors produced\n";
  cerr << " -b             use actual bit numbers instead of sequential numbers for descriptors\n";
  cerr << " -w <n>         interpret fixed width fp's as numbers of width <n> bytes\n";
  cerr << " -u             gsub space to _ in multi-token names\n";
  cerr << " -d <nbits>     fold sparse fingerprints to constant width <nbits>\n";
  cerr << " -q             if a bit is absent write as -1 rather than 0\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

class Bit_and_Count
{
  private:
    unsigned int _bit;
    unsigned int _count;
  public:
    unsigned int bit() const { return _bit;}
    void set_bit(unsigned int s) { _bit = s;}

    unsigned int count() const { return _count;}
    void set_count(unsigned int s) { _count = s;}
};

/*static int
not_enough_non_zero_values(const int * tmp,
                           int number_descriptors)
{
  int number_non_zero_values = 0;
  for (int i = 0; i < number_non_zero_values; i++)
  {
    if (0 == tmp[i])
      continue;

    number_non_zero_values++;

    if (number_non_zero_values >= min_non_zero_values_needed_on_each_row)
      return 0;    // this record is OK
  }

  if (number_non_zero_values < min_non_zero_values_needed_on_each_row)
    return 1;   // too few

  return 0;    // we have enough
}*/

static int
write_descriptors(const IWString & pcn,
                  const Bit_and_Count * bac, 
                  int number_descriptors, 
                  IWString & output_buffer)
{
  output_buffer << '"' << pcn << '"';

  for (int i = 0; i < number_descriptors; i++)
  {
    if (0 == bac[i].count())
    {
      number_descriptors = i + 1;
      break;
    }
  }

  output_buffer << ',' << number_descriptors;

  for (int i = 0; i < number_descriptors; i++)
  {
    const Bit_and_Count & baci = bac[i];

    if (0 == baci.count())
      break;

    output_buffer << ',' << baci.bit() << ',' << baci.count();
  }

  output_buffer << '\n';

  return 1;
}

static int
discard_because_of_dash_y (const IWString & pcn,
                           int nbits)
{
  if (0 == min_non_zero_values_needed_on_each_row)
    return 0;

  if (0 == nbits)
    ;
  else if (nbits >= min_non_zero_values_needed_on_each_row)
    return 0;

  if (write_identifiers_suppressed_by_dash_y)
    stream_for_not_enough_non_zero_values << pcn << '\n';

  records_discarded_for_not_enough_non_zero_values++;

  return 1;
}

/*
  One could argue that it would be more efficient to detect the
  too-sparse fingerprints earlier, but putting it here makes for
  very simple code
*/

static int
write_descriptors (const IWString & pcn,
                   const int * tmp, 
                   int number_descriptors, 
                   IWString & output_buffer)
{
  output_buffer << pcn;

  for (int i = 0; i < number_descriptors; i++)
  {
    if (tmp[i] >= 0)
      iwdigits.append_number (output_buffer, tmp[i]);
    else if (write_zero_count_as_negative_one)
      output_buffer << " -1";
    else
    {
      cerr << "Unexpected negative value, i = " << i << " value " << tmp[i] << endl;
      output_buffer << ' ' << tmp[i];
    }
  }

  output_buffer.add ('\n');

  return 1;
}

static int
fill_non_colliding_fingerprint_output (const const_IWSubstring & buffer,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         int * tmp)
{
  Sparse_Fingerprint fp;
  (void) fp.construct_from_daylight_ascii_representation (buffer);

  int istart = 0;
  unsigned int zbit;
  int zcount;

  int rc = 0;

  while (fp.next_bit_set (istart, zbit, zcount))
  {
    IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = col.find (zbit);
    if (f == col.end())
      continue;

//  cerr << "Bit " << zbit << " found in column " << f->second << endl;
    tmp[(*f).second] = zcount;
    rc++;
  }

  return rc;
}

static int
fill_fixed_width_fingerprint_output_numeric2(const unsigned short * s,
                         int nbits,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         int * tmp)
{
  int nshort = nbits / 16;

  int rc = 0;

  for (int i = 0; i < nshort; i++)
  {
    if (0 == s[i])
      continue;

    IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = col.find(i);
    if (f == col.end())
      continue;

    tmp[(*f).second] = s[i];
    rc++;
  }

  return rc;
}

static int
fill_fixed_width_fingerprint_output_numeric1(const unsigned char * s,
                         int nbits,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         int * tmp)
{
  int nbytes = nbits / 8;

  int rc = 0;

  for (int i = 0; i < nbytes; i++)
  {
    if (0 == s[i])
      continue;

    IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = col.find(i);
    if (f == col.end())
      continue;

    tmp[(*f).second] = s[i];
    rc++;
  }

  return rc;
}

static int
fill_fixed_width_fingerprint_output_numeric(const IWDYFP & fp,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         int * tmp)
{
  if (0 != fp.nbits() % 8)
  {
    cerr << "Sorry, cannot handle fingerprint with " << fp.nbits() << " bits\n";
    return 0;
  }

  if (1 == interpret_fixed_width_as_numeric)
    return fill_fixed_width_fingerprint_output_numeric1(fp.bits(), fp.nbits(), col, tmp);
  else if (2 == interpret_fixed_width_as_numeric)
    return fill_fixed_width_fingerprint_output_numeric2(reinterpret_cast<const unsigned short *>(fp.bits()), fp.nbits(), col, tmp);
  else
    return 0;
}

static int
fill_fixed_width_fingerprint_output (const const_IWSubstring & buffer,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         int * tmp)
{
  IWDYFP fp;
  (void) fp.construct_from_daylight_ascii_representation (buffer);

  int rc = 0;

  if (interpret_fixed_width_as_numeric)
    return fill_fixed_width_fingerprint_output_numeric(fp, col, tmp);

  for (unsigned int i = 0; i < nbits; i++)
  {
    if (! fp.is_set (i))
      continue;

    IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = col.find (i);
    if (f == col.end())
      continue;

    tmp[(*f).second] = 1;
    rc++;
  }

  return rc;
}

static int
fill_fingerprint_output(const const_IWSubstring & buffer,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         Bit_and_Count * bac)
{
  const_IWSubstring tmp (buffer);
  tmp.remove_leading_chars (fptag.length());
  tmp.chop();

  Sparse_Fingerprint fp;
  (void) fp.construct_from_daylight_ascii_representation (tmp);

  int istart = 0;
  unsigned int zbit;
  int zcount;

  int ndx = 0;

  while (fp.next_bit_set (istart, zbit, zcount))
  {
    IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = col.find (zbit);
    if (f == col.end())
      continue;

    if (write_actual_bit_numbers)
      bac[ndx].set_bit(zbit);
    else
      bac[ndx].set_bit((*f).second);

    bac[ndx].set_count(zcount);
    ndx++;
  }

  bac[ndx].set_count(0);

  return 1;
}

static int
fill_fingerprint_output (const const_IWSubstring & buffer,
                         const IW_Hash_Map<unsigned int, unsigned int> & col,
                         int * tmp)
{
  const_IWSubstring fp (buffer);
  fp.remove_leading_chars (fptag.length());
  fp.chop();

  int n = col.size();

  if (fold_to_constant_width > 0 && fold_to_constant_width < n)
    n = fold_to_constant_width;

  if (write_zero_count_as_negative_one)
    std::fill_n(tmp, n, -1);
  else
    std::fill_n(tmp, n, 0);

  if (non_colliding_form)
    return fill_non_colliding_fingerprint_output (fp, col, tmp);
  else
    return fill_fixed_width_fingerprint_output (fp, col, tmp);
}

static int
rejected_by_upper_or_lower_population (unsigned int zbit,
                                       unsigned int zcount,
                                       int & bits_below_threshold,
                                       int & bits_above_threshold)
{
  float ratio = static_cast<float> (zcount) / static_cast<float> (fingerprints_read);

  if (verbose > 2)
    cerr << "Fingerprint " << zbit << " in " << zcount << " molecules, ratio " << ratio << endl;

  if (lower_support_cutoff_number > 0 && zcount < lower_support_cutoff_number)
  {
    bits_below_threshold++;
    return 1;
  }

  if (upper_support_cutoff_number > 0 && zcount > upper_support_cutoff_number)
  {
    bits_above_threshold++;
    return 1;
  }

  return 0;
}

static int
report_bits_above_and_below_threshold (ostream & os,
                                       int bits_below_threshold,
                                       int bits_above_threshold)
{
  if (lower_support_cutoff_number > 0)
    os << bits_below_threshold << " bits below threshold " << lower_support_cutoff_number << '\n';
    
  if (upper_support_cutoff_number > 0)
    os << bits_above_threshold << " bits above threshold " << upper_support_cutoff_number << '\n';

  return os.good();
}

static int
bit_and_count_comparitor (const void * v1, const void * v2)
{
  const pair<unsigned int, unsigned int> * p1 = (const pair<unsigned int, unsigned int> *) v1;
  const pair<unsigned int, unsigned int> * p2 = (const pair<unsigned int, unsigned int> *) v2;

  if (p1->second < p2->second)
    return 1;

  if (p1->second > p2->second)
    return -1;

  return 0;
}

static int
do_sort_by_frequency (const IW_Hash_Map<unsigned int, unsigned int> & fpcount,
                      IW_Hash_Map<unsigned int, unsigned int> & col,
                      IWString & header,
                      pair<unsigned int, unsigned int> * bit_and_count)
{
  int number_descriptors = 0;
  int bits_below_threshold = 0;
  int bits_above_threshold = 0;

//for (IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = fpcount.begin(); f != fpcount.end(); ++f)
  for (auto f = fpcount.cbegin(); f != fpcount.cend(); ++f)
  {
    if (0 == lower_support_cutoff_number && 0 == upper_support_cutoff_number)
      ;
    else if (rejected_by_upper_or_lower_population ((*f).first, (*f).second, bits_below_threshold, bits_above_threshold))
      continue;

    pair<unsigned int, unsigned int> & pi = bit_and_count[number_descriptors];
    number_descriptors++;

    pi.first  = (*f).first;
    pi.second = (*f).second;
  }

  if (verbose)
    report_bits_above_and_below_threshold (cerr, bits_below_threshold, bits_above_threshold);

  if (0 == number_descriptors)
    return 0;

  qsort (bit_and_count, number_descriptors, sizeof (pair<unsigned int, unsigned int>), bit_and_count_comparitor);

  for (int i = 0; i < number_descriptors; i++)
  {
    const pair<unsigned int, unsigned int> & pi = bit_and_count[i];

    unsigned int b = pi.first;

    if (verbose > 2)
      cerr << "Bit " << b << " found " << pi.second << " times\n";

    col[b] = i;

    header << descriptor_prefix << b;
  }

  return number_descriptors;
}

static int
do_sort_by_frequency (const IW_Hash_Map<unsigned int, unsigned int> & fpcount,
                      IW_Hash_Map<unsigned int, unsigned int> & col,
                      IWString & header)
{
  pair<unsigned int, unsigned int> * bit_and_count = new pair<unsigned int, unsigned int>[fpcount.size()];

  int rc = do_sort_by_frequency (fpcount, col, header, bit_and_count);

  delete [] bit_and_count;

  return rc;
}

template <typename T>
int
fingerprints_to_descriptors (iwstring_data_source & input,
                             const IW_Hash_Map<unsigned int, unsigned int> & col,
                             T * tmp,
                             IWString_and_File_Descriptor & output)
{
  const_IWSubstring input_buffer;

  int ncols;
  if (fold_to_constant_width > 0)
    ncols = fold_to_constant_width;
  else
    ncols = col.size();

  IWString pcn;
  int nbits;

  while (input.next_record (input_buffer))
  {
    if ('|' == input_buffer)
    {
      if (0 == pcn.length())
      {
        cerr << "Yipes, TDT with no '" << identifier_tag << "' tag\n";
        return 0;
      }

      if (discard_because_of_dash_y(pcn, nbits))
        ;
      else
        write_descriptors (pcn, tmp, ncols, output);

      pcn.resize_keep_storage (0);

      output.write_if_buffer_holds_more_than(32768);
    }
    else if (input_buffer.starts_with (fptag))
    {
      nbits = fill_fingerprint_output (input_buffer, col, tmp);
    }
    else if (input_buffer.starts_with (identifier_tag))
    {
      pcn = input_buffer;
      pcn.remove_leading_chars (identifier_tag.length());
      pcn.chop();
      if (gsub_multi_token_names)
        pcn.gsub (' ', '_');
      else
        pcn.truncate_at_first(' ');
    }
  }

  output.flush();

  return 1;
}

static int
scan_fixed_width_fingerprint_numeric2(const unsigned short * s,
                                      int nbits,
                                      IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  int nshort = nbits / 16;

  for (int i = 0; i < nshort; i++)
  {
    if (s[i])
      fpcount[i]++;
  }

  return 1;
}

static int
scan_fixed_width_fingerprint_numeric1(const unsigned char * s,
                                      int nbits,
                                      IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  int nbytes = nbits / 8;

  for (int i = 0; i < nbytes; i++)
  {
    if (s[i])
      fpcount[i]++;
  }

  return 1;
}

static int
scan_fixed_width_fingerprint_numeric(const IWDYFP & fp,
                          IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  if (1 == interpret_fixed_width_as_numeric)
    return scan_fixed_width_fingerprint_numeric1(fp.bits(), fp.nbits(), fpcount);
  else if (2 == interpret_fixed_width_as_numeric)
    return scan_fixed_width_fingerprint_numeric2(reinterpret_cast<const unsigned short *>(fp.bits()), fp.nbits(), fpcount);
  else
    return 0;
}

static int
scan_fixed_width_fingerprint (const IWDYFP & fp,
                  IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  unsigned int nb = fp.nbits();

  if (0 == nbits)
  {
    nbits = nb;
    if (verbose)
      cerr << "Fixed width fingerprint " << fptag << " has " << nbits << " bits\n";
  }
  else if (nb != nbits)
  {
    cerr << "Bit count mismatch, expected " << nbits << " got " << nb << endl;
    return 0;
  }

  if (interpret_fixed_width_as_numeric)
    return scan_fixed_width_fingerprint_numeric(fp, fpcount);

  for (unsigned int i = 0; i < nb; i++)
  {
    if (fp.is_set(i))
      fpcount[i]++;
  }

  return 1;
}

static int
scan_fixed_width_fingerprint (const const_IWSubstring & buffer,
                  IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  IWDYFP fp;

  if (! fp.construct_from_daylight_ascii_representation (buffer))
  {
    cerr << "Invalid fixed width fingerprint '" << buffer << "'\n";
    return 0;
  }

  return scan_fixed_width_fingerprint (fp, fpcount);
}

static int
scan_hex_encoded_fixed_width_fingerprint (const const_IWSubstring & buffer,
                  IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  IWDYFP fp;

  if (! fp.construct_from_hex (buffer))
  {
    cerr << "Invalid hex encoding\n";
    return 0;
  }

  return scan_fixed_width_fingerprint (fp, fpcount);
}

static int
scan_non_colliding_fingerprint (const const_IWSubstring & buffer,
                  IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  if (0 == buffer.length())    // empty fingerprint
    return 1;

  Sparse_Fingerprint fp;
  if (! fp.construct_from_daylight_ascii_representation (buffer))
  {
    cerr << "Invalid sparse fingerprint specification '" << buffer << "'\n";
    return 0;
  }

  int istart = 0;
  unsigned int zbit;
  int zcount;

  while (fp.next_bit_set (istart, zbit, zcount))
  {
    fpcount[zbit]++;
  }

  return 1;
}

static int
scan_fingerprint (const const_IWSubstring & buffer,
                  IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  assert (buffer.ends_with ('>'));

  const_IWSubstring fp (buffer);
  fp.remove_leading_chars (fptag.length());
  fp.chop();

  if (non_colliding_form)
    return scan_non_colliding_fingerprint (fp, fpcount);
  else if (input_is_hex_encoded)
    return scan_hex_encoded_fixed_width_fingerprint (fp, fpcount);
  else
    return scan_fixed_width_fingerprint (fp, fpcount);
}

static int
fingerprints_to_descriptors (iwstring_data_source & input,
                             IWString_and_File_Descriptor & output)
{
  IW_Hash_Map<unsigned int, unsigned int> fpcount;

  int got_fingerprint_this_tdt = 0;

  const_IWSubstring buffer;
  while (input.next_record (buffer))
  {
    if ('|' == buffer)
    {
      if (! got_fingerprint_this_tdt)
      {
        cerr << "Fingerprint without '" << fptag << "', near line " << input.lines_read() << endl;

        if (0 == fingerprints_read)    // For example the Daylight fingerprint with a header TDT
          continue;
        else
          return 0;
      }

      fingerprints_read++;

      got_fingerprint_this_tdt = 0;
    }
    else if (buffer.starts_with (fptag))
    {
      if (got_fingerprint_this_tdt)
      {
        cerr << "Duplicate '" << fptag << "' within TDT, impossible\n";
        return 0;
      }

      if (! scan_fingerprint (buffer, fpcount))
      {
        cerr << "Fatal error scanning fingerprint on line " << input.lines_read() << endl;
        cerr << buffer << endl;
        return 0;
      }

      got_fingerprint_this_tdt = 1;
    }
  }

  if (0 == fpcount.size())
  {
    cerr << "Very strange, no fingerprints counted\n";
    return 0;
  }

  if (verbose)
    cerr << "Found " << fpcount.size() << " bits set in " << fingerprints_read << " fingerprints\n";

  if (verbose > 1)
  {
    for (const auto i : fpcount)
    {
      cerr << " bit " << i.first << " hit " << i.second << " time\n";
    }
  }

// Convert any ratios specified for support levels to numbers

  if (zero != lower_population_cutoff)
  {
    lower_support_cutoff_number = static_cast<unsigned int> (lower_population_cutoff * static_cast<double>(fingerprints_read) + 0.4999);
    if (0 == lower_population_cutoff)
      lower_population_cutoff = 1;

    if (verbose)
      cerr << "Will discard bits that occur in fewer than " << lower_population_cutoff << " molecules\n";
  }

  if (zero != upper_population_cutoff)
  {
    upper_support_cutoff_number = static_cast<unsigned int> (upper_population_cutoff * static_cast<double>(fingerprints_read) + 0.4999);
    assert (upper_support_cutoff_number > 1);
    if (fingerprints_read == upper_support_cutoff_number)
      upper_support_cutoff_number = fingerprints_read - 1;

    if (verbose)
      cerr << "Will discard bits that occur in more than " << upper_support_cutoff_number << " molecules\n";
  }

  int number_descriptors = 0;

  IWString header ("Name");
  header.resize (10 * fpcount.size());

  IW_Hash_Map<unsigned int, unsigned int> col;

//cerr << "Building header, lower_support_cutoff_number " << lower_support_cutoff_number << " upper_support_cutoff_number " << upper_support_cutoff_number << ", fold_to_constant_width " << fold_to_constant_width << endl;

  if (sort_by_frequency)
    number_descriptors = do_sort_by_frequency (fpcount, col, header);
  else if (lower_support_cutoff_number > 0 || upper_support_cutoff_number > 0)
  {
    int bits_below_threshold = 0;
    int bits_above_threshold = 0;

    for (auto f = fpcount.cbegin(); f != fpcount.cend(); ++f)
    {
      if (verbose > 2)
        cerr << "Fingerprint " << (*f).first << " in " << (*f).second << " molecules\n";

      if (lower_support_cutoff_number > 0 && (*f).second < lower_support_cutoff_number)
        bits_below_threshold++;
      else if (upper_support_cutoff_number > 0 && (*f).second > upper_support_cutoff_number)
        bits_above_threshold++;
      else
      {
        header << descriptor_prefix << (*f).first;
        col[(*f).first] = number_descriptors;
        number_descriptors++;
      }
    }

    if (verbose)
      report_bits_above_and_below_threshold (cerr, bits_below_threshold, bits_above_threshold);
  }
  else if (fold_to_constant_width > 0)
  {
    int * to_col = new_int(fold_to_constant_width); unique_ptr<int[]> free_to_col(to_col);

    for (auto f = fpcount.cbegin(); f != fpcount.cend(); ++f)
    {
      const auto b = f->first % fold_to_constant_width;
      to_col[b]++;
      col[f->first] = b;
    }

    for (auto i = 0; i < fold_to_constant_width; ++i)
    {
      header << descriptor_prefix << i;
    }

    if (verbose)
    {
      cerr << "After folding to " << fold_to_constant_width << " bits, how many of the initial bits hit each fixed width bit\n";
      int zero_bits_map_here = 0;
      Accumulator_Int<int> acc;

      for (auto i = 0; i < fold_to_constant_width; ++i)
      {
        if (0 == to_col[i])
          zero_bits_map_here++;
        else
          acc.extra(to_col[i]);

        if (verbose > 1)
          cerr << " " << to_col[i] << " non colliding bits hashed to bit " << i << endl;
      }

      cerr << zero_bits_map_here << " of " << fold_to_constant_width << " fixed width bits never set\n";
      cerr << "Of bits set, set between " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;
    }

    number_descriptors = fold_to_constant_width;
  }
  else
  {
    number_descriptors = 0;

    for (auto i = fpcount.cbegin(); i != fpcount.cend(); ++i)
    {
      unsigned int b = (*i).first;
      if (write_actual_bit_numbers)
        header << descriptor_prefix << b;
      else
        header << descriptor_prefix << number_descriptors;

      col[b] = number_descriptors;
      number_descriptors++;
    }
  }

  if (0 == number_descriptors)
  {
    cerr << "Bad news, support levels " << lower_support_cutoff_number << " and " << upper_support_cutoff_number << " leave no fingerprints\n";
    return 0;
  }

  header << '\n';

  if (! produce_sparse_form_for_nikil)
    output << header;

  if (verbose)
    cerr << "Will produce " << number_descriptors << " descriptors\n";

  if (! input.seekg (static_cast<off_t> (0)))
  {
    cerr << "Yipes, cannot seek back to start of file\n";
    return 0;
  }

  int rc;

  if (produce_sparse_form_for_nikil)
  {
    Bit_and_Count * bac = new Bit_and_Count[number_descriptors];
    rc = fingerprints_to_descriptors(input, col, bac, output);
    delete [] bac;
  }
  else
  {
    int * tmp = new_int (number_descriptors); std::unique_ptr<int[]> free_tmp (tmp);
    if (NULL == tmp)
    {
      cerr << "Bad news, cannot allocate " << number_descriptors << " objects\n";
      return 0;
    }

    rc = fingerprints_to_descriptors(input, col, tmp, output);
  }

  return rc;
}


static int
fingerprints_to_descriptors (const char * fname,
                             IWString_and_File_Descriptor & output)
{
  iwstring_data_source input (fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return fingerprints_to_descriptors (input, output);
}

static int
fingerprints_to_descriptors (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vF:fhx:X:P:sbmy:Y:un:N:kw:d:q");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (! cl.option_present ('F'))
  {
    cerr << "Must specify fingerprint to echo via the -F option\n";
    usage (7);
  }

  if (cl.option_present ('F'))
  {
    fptag = cl.string_value ('F');

    if (! fptag.ends_with ('<'))
      fptag += '<';

    if (verbose)
      cerr << "Will process fingerprint '" << fptag << "'\n";
  }

  if (cl.option_present ('f'))
  {
    non_colliding_form = 0;

    if (verbose)
      cerr << "Will process as fixed width fingerprint\n";

    if (cl.option_present('w'))
    {
      if (! cl.value('w', interpret_fixed_width_as_numeric) || interpret_fixed_width_as_numeric < 1)
      {
        cerr << "The interpret fixed width as numeric (-w) must be a whole +ve number\n";
        usage(3);
      }

      if (verbose)
        cerr << "Will interpret fixed width fingerprints as numeric of width " << interpret_fixed_width_as_numeric << " bytes\n";
    }
  }

  if (cl.option_present ('h'))
  {
    non_colliding_form = 0;

    input_is_hex_encoded = 1;

    if (verbose)
      cerr << "Hex encoded fixed width form\n";
  }

  if (cl.option_present('d'))
  {
    if (cl.option_present('b'))
    {
      cerr << "Sorry, the -b option cannot be used with the -d option\n";
      usage(1);
    }

    if (! cl.value('d', fold_to_constant_width) || fold_to_constant_width < 1)
    {
      cerr << "The fold to constant width option (-d) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will fold fingerprints to a constant width of " << fold_to_constant_width << " bits\n";
  }

  if (cl.option_present('q'))
  {
    write_zero_count_as_negative_one = 1;

    if (verbose)
      cerr << "Will write zero count bits as negative one\n";
  }


  int support_level_specified = 0;
  if (cl.option_present('x'))
    support_level_specified = 1;
  else if (cl.option_present('X'))
    support_level_specified++;
  else if (cl.option_present('n'))
    support_level_specified++;
  else if (cl.option_present('N'))
    support_level_specified++;

  if (cl.option_present('y'))
  {
    if (0 == support_level_specified)
    {
      cerr << "The -y option only makes sense when used with the -x,-n and/or -X,-N options\n";
      usage(4);
    }

    if (! cl.value('y', min_non_zero_values_needed_on_each_row) || min_non_zero_values_needed_on_each_row < 1)
    {
      cerr << "The minimum number of non zero columns (-y) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard records that don't have " << min_non_zero_values_needed_on_each_row << " or more non-zero values\n";

    if (cl.option_present('Y'))
    {
      const char * y = cl.option_value('Y');

      stream_for_not_enough_non_zero_values.open(y);

      if (! stream_for_not_enough_non_zero_values.good())
      {
        cerr << "Cannot open -Y file '" << y << "'\n";
        return 3;
      }

      if (verbose)
        cerr << "Identifiers discarded by -Y value written to '" << y << "'\n";

      write_identifiers_suppressed_by_dash_y = 1;
    }
  }

  if (cl.option_present('x') && cl.option_present('n'))
  {
    cerr << "Sorry, the -x and -n options are mutually exclusive\n";
    usage(4);
  }

  if (cl.option_present('X') && cl.option_present('N'))
  {
    cerr << "Sorry, the -X and -N options are mutually exclusive\n";
    usage(4);
  }

  if (cl.option_present ('x'))
  {
    if (! cl.value ('x', lower_population_cutoff) || lower_population_cutoff < 0.0 || lower_population_cutoff >= 1.0)
    {
      cerr << "The lower population cutoff must be a valid fraction\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will discard bits set less than " << lower_population_cutoff << " of the time\n";
  }
  else if (cl.option_present('n'))
  {
    if (! cl.value('n', lower_support_cutoff_number) || lower_support_cutoff_number < 1)
    {
      cerr << "The lower support level number (-n) option must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard bits with fewer than " << lower_support_cutoff_number << " molecules hit\n";
  }

  if (cl.option_present ('X'))
  {
    if (! cl.value ('X', upper_population_cutoff) || upper_population_cutoff < 0.0 || upper_population_cutoff >= 1.0)
    {
      cerr << "The upper population cutoff must be a valid fraction\n";
      usage (7);
    }

    if (verbose)
      cerr << "Will discard bits set more than " << upper_population_cutoff << " of the time\n";
  }
  else if (cl.option_present('N'))
  {
    if (! cl.value('N', upper_support_cutoff_number) || upper_support_cutoff_number < 1)
    {
      cerr << "The upper support level number (-N) option must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard bits with more than " << upper_support_cutoff_number << " molecules hit\n";

    if (lower_support_cutoff_number > 0 && lower_support_cutoff_number > upper_support_cutoff_number)
    {
      cerr << "Inconsistent -n (" << lower_support_cutoff_number << ") and -N (" << upper_support_cutoff_number << ") option values\n";
      return 4;
    }
  }

  if (0 == upper_population_cutoff)
    ;
  else if (lower_population_cutoff > upper_population_cutoff)
  {
    cerr << "Impossible lower (" << lower_population_cutoff << ") and upper (" << upper_population_cutoff << ") specification\n";
    return 9;
  }

  if (cl.option_present ('P'))
  {
    descriptor_prefix = cl.string_value ('P');

    if (verbose)
      cerr << "Descriptors will have the prefix '" << descriptor_prefix << "'\n";

    if (! descriptor_prefix.starts_with (' '))
      descriptor_prefix.insert_before (0, ' ');
  }

  if (cl.option_present('u'))
  {
    gsub_multi_token_names = 1;
    if (verbose)
      cerr << "Will gsub to _ any multi-token names\n";
  }

  if (cl.option_present ('s'))
  {
    sort_by_frequency = 1;

    if (verbose)
      cerr << "Bits will be output in order of frequency of occurrence in the set\n";
  }

  if (cl.option_present ('b'))
  {
    write_actual_bit_numbers = 1;

    if (verbose)
      cerr << "Will write bit numbers rather than sequential names\n";
  }

  if (cl.option_present('k'))
  {
    produce_sparse_form_for_nikil = 1;
    if (verbose)
      cerr << "Will produce sparse form for Nikil\n";
  }

  iwdigits.set_include_leading_space (1);
  iwdigits.initialise (100);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! fingerprints_to_descriptors(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    if (records_discarded_for_not_enough_non_zero_values)
      cerr << records_discarded_for_not_enough_non_zero_values << " discarded for not enough non-zero values\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fingerprints_to_descriptors (argc, argv);

  return rc;
}
