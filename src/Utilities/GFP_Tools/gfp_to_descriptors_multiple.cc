// Converts multiple fingerprints to descriptor form.

#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static IWString identifier_tag("PCN<");

static const float zero = static_cast<float>(0.0);

static float lower_population_cutoff = zero;

static unsigned int lower_support_cutoff_number = 0;

static float upper_population_cutoff = zero;

static unsigned int upper_support_cutoff_number = std::numeric_limits<int>::max();

/*
  We can quickly check whether or not lower or upper limits have been specified
*/

static int bit_count_limits_set = 0;

static int min_non_zero_values_needed_on_each_row = 0;

static IWString_and_File_Descriptor stream_for_not_enough_non_zero_values;

static int write_identifiers_suppressed_by_dash_y = 0;

static int records_discarded_for_not_enough_non_zero_values = 0;

static int gsub_multi_token_names = 0;

static IWString descriptor_prefix;

static IWDigits iwdigits;

static int fold_to_constant_width = 0;

static int * collision_count = nullptr;

#ifdef DEPRECATE_DASH_B_OPTION
/*
  By default, we create arbitrary numbers
*/

static int write_actual_bit_numbers = 1;
#endif

/*
  It can be convenient to sort the bits so the most common occur as the first descriptors
*/

static int sort_by_frequency = 0;

static int write_zero_count_as_negative_one = 0;

static char output_separator = ' ';

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "Converts a fingerprint set from a gfp file to a descriptor file\n";
  cerr << " -F <tag>       specify tag of the fingerprint to be processed, standard .gfp syntax\n";
//cerr << " -m             molecular properties\n";
  cerr << " -x <fraction>  discard bits that hit less than <fraction> of the time\n";
  cerr << " -n <number>    discard bits hit in less than <number> of the fingerprints\n";
  cerr << " -X <fraction>  discard bits that hit more than <fraction> of the time\n";
  cerr << " -N <number>    discard bits hit in more than <number> of the fingerprints\n";
  cerr << "                    use either -x or -n, and either -X or -N\n";
  cerr << " -y <number>    discard records unless they have at least <number> non zero values\n";
  cerr << " -s             sort columns by frequency (within fingeprint)\n";
  cerr << " -I <string>    prefix for descriptors produced\n";
  cerr << " -u             gsub space to _ in multi-token names\n";
  cerr << " -d <nbits>     fold sparse fingerprints to constant width <nbits>\n";
  cerr << " -q             if a bit is absent write as -1 rather than 0\n";
  cerr << " -o <sep>       output separator (def ' ')\n";
  cerr << " -v             verbose output\n";

  exit(rc);
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

static int
discard_because_of_dash_y(const IWString & pcn,
                          const int nbits)
{
//cerr << "cmp " << nbits << " with " << min_non_zero_values_needed_on_each_row << endl;

  if (nbits >= min_non_zero_values_needed_on_each_row)   // do not discard
    return 0;

  if (write_identifiers_suppressed_by_dash_y)
    stream_for_not_enough_non_zero_values << pcn << '\n';

  records_discarded_for_not_enough_non_zero_values++;

  return 1;
}

static int
write_descriptors(const IWString & pcn,
                  const int * tmp, 
                  const int number_descriptors, 
                  IWString_and_File_Descriptor & output)
{
  output << pcn;

  for (int i = 0; i < number_descriptors; i++)
  {
    if (tmp[i] >= 0)
      iwdigits.append_number(output, tmp[i]);
    else if (write_zero_count_as_negative_one)
      output << output_separator << "-1";
    else
    {
      cerr << "Unexpected negative value, i = " << i << " value " << tmp[i] << endl;
      output << ' ' << tmp[i];
    }
  }

  output.add('\n');

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
fill_output_array(const Sparse_Fingerprint & fp,
                  const IW_Hash_Map<unsigned int, unsigned int> & col,
                  int * tmp)
{
  int istart = 0;
  unsigned int zbit;
  int zcount;

  int rc = 0;

  while (fp.next_bit_set(istart, zbit, zcount))
  {
    IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = col.find(zbit);
    if (f == col.end())
      continue;

//  cerr << "Bit " << zbit << " found in column " << f->second << endl;
    tmp[(*f).second] = zcount;
    rc++;
  }

  return rc;
}

static int
fill_output_array(const IWDYFP & fp,
                  const std::unordered_map<unsigned int, unsigned int> & col,
                  int * tmp)
{
  int i = 0;
  int b;

  int rc = 0;

  while((b = fp.next_on_bit(i)) >= 0)
  {
    const auto f = col.find(b);

    if (f == col.end())
      continue;

    tmp[f->second]++;
    rc++;
  }

  return rc;
}

static int
rejected_by_upper_or_lower_population(unsigned int zbit,
                                       unsigned int zcount,
                                       const int pool_size,
                                       int & bits_below_threshold,
                                       int & bits_above_threshold)
{
  if (! bit_count_limits_set)
    return 0;

  float ratio = static_cast<float>(zcount) / static_cast<float>(pool_size);

  if (verbose > 2)
    cerr << "Fingerprint " << zbit << " in " << zcount << " molecules, ratio " << ratio << endl;

  if (zcount < lower_support_cutoff_number)
  {
    bits_below_threshold++;
    return 1;
  }

  if (zcount > upper_support_cutoff_number)
  {
    bits_above_threshold++;
    return 1;
  }

  return 0;
}

static int
report_bits_above_and_below_threshold(std::ostream & os,
                                       int bits_below_threshold,
                                       int bits_above_threshold)
{
  if (! bit_count_limits_set)
    return 1;

  if (0 != lower_support_cutoff_number)
    os << bits_below_threshold << " bits below threshold " << lower_support_cutoff_number << '\n';
    
  if (std::numeric_limits<int>::max() != upper_support_cutoff_number)
    os << bits_above_threshold << " bits above threshold " << upper_support_cutoff_number << '\n';

  return os.good();
}

static int
establish_header_and_columns_folded(const int old_to_constant_width,
                      const IW_Hash_Map<unsigned int, unsigned int> * fpcount,
                      const int nfp,
                      IW_Hash_Map<unsigned int, unsigned int> * col,
                      const IWString * dprefix,
                      int & number_descriptors,
                      IWString & header)
{
  int * to_col = new_int(fold_to_constant_width); std::unique_ptr<int[]> free_to_col(to_col);

  int bits_below_threshold = 0;
  int bits_above_threshold = 0;

  for (int i = 0; i < nfp; ++i)
  {
    for (auto f : fpcount[i])
    {
      if (f.second < lower_support_cutoff_number)
        bits_below_threshold++;
      else if (f.second > upper_support_cutoff_number)
        bits_above_threshold++;
      else
      {
        const auto b = f.first % fold_to_constant_width;
        if (0 != to_col[b])
          collision_count[b]++;
        to_col[b]++;
        col[i][f.first] = b;
      }
    }
  }

  for (auto i = 0; i < fold_to_constant_width; ++i)
  {
    header << dprefix[0] << i;
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

  if (verbose)
    report_bits_above_and_below_threshold(cerr, bits_below_threshold, bits_above_threshold);

  return 1;
}

static int
establish_header_and_columns_with_thresholds(const IW_Hash_Map<unsigned int, unsigned int> * fpcount,
                      const int nfp,
                      IW_Hash_Map<unsigned int, unsigned int> * col,
                      const IWString * dprefix,
                      int & number_descriptors,
                      IWString & header)
{
  int bits_below_threshold = 0;
  int bits_above_threshold = 0;

  for (int i = 0; i < nfp; ++i)
  {
    for (auto f : fpcount[i])
    {
      if (verbose > 2)
        cerr << "bit " << f.first << " in " << f.second << " molecules\n";

      if (f.second < lower_support_cutoff_number)
        bits_below_threshold++;
      else if (f.second > upper_support_cutoff_number)
        bits_above_threshold++;
      else
      {
        header << dprefix[i] << f.first;
        col[i][f.first] = number_descriptors;
        number_descriptors++;
      }
    }
  }

  if (verbose)
    report_bits_above_and_below_threshold(cerr, bits_below_threshold, bits_above_threshold);

  return 1;
}

static int
bit_and_count_comparitor(const void * v1, const void * v2)
{
  const std::pair<unsigned int, unsigned int> * p1 = (const std::pair<unsigned int, unsigned int> *) v1;
  const std::pair<unsigned int, unsigned int> * p2 = (const std::pair<unsigned int, unsigned int> *) v2;

  if (p1->second < p2->second)
    return 1;

  if (p1->second > p2->second)
    return -1;

  return 0;
}

static int
do_sort_by_frequency(const IW_Hash_Map<unsigned int, unsigned int> & fpcount,
                      IW_Hash_Map<unsigned int, unsigned int> & col,
                      const IWString & descriptor_prefix,
                      int & number_descriptors,
                      IWString & header,
                      const int pool_size)
{
  std::pair<unsigned int, unsigned int> * bit_and_count = new std::pair<unsigned int, unsigned int>[fpcount.size()]; std::unique_ptr<std::pair<unsigned int, unsigned int>[]> free_bit_and_count(bit_and_count);

  int bits_below_threshold = 0;
  int bits_above_threshold = 0;

  int ndx = 0;

//for (IW_Hash_Map<unsigned int, unsigned int>::const_iterator f = fpcount.begin(); f != fpcount.end(); ++f)
  for (auto f = fpcount.cbegin(); f != fpcount.cend(); ++f)
  {
    if (rejected_by_upper_or_lower_population((*f).first, (*f).second, pool_size, bits_below_threshold, bits_above_threshold))
      continue;

    std::pair<unsigned int, unsigned int> & pi = bit_and_count[ndx];
    ndx++;

    pi.first  = (*f).first;
    pi.second = (*f).second;
  }

  if (verbose)
    report_bits_above_and_below_threshold(cerr, bits_below_threshold, bits_above_threshold);

  if (0 == ndx)
    return 0;

  qsort(bit_and_count, ndx, sizeof(std::pair<unsigned int, unsigned int>), bit_and_count_comparitor);

  for (int i = 0; i < ndx; i++)
  {
    const std::pair<unsigned int, unsigned int> & pi = bit_and_count[i];

    unsigned int b = pi.first;

    if (verbose > 2)
      cerr << "Bit " << b << " found " << pi.second << " times\n";

    col[b] = i;

    header << descriptor_prefix << b;
  }

  number_descriptors += ndx;

  return number_descriptors;
}

int
fingerprints_to_descriptors(const IW_General_Fingerprint * pool,
                            const int pool_size,
                            const int nfp,
                            const int * fpnumber,
                            const std::unordered_map<unsigned int, unsigned int> * col,
                            const int number_descriptors,
                            IWString_and_File_Descriptor & output)
{
  int * tmp = new int[number_descriptors]; std::unique_ptr<int[]> free_tmp(tmp);
  if (nullptr == tmp)
  {
    cerr << "Cannot allocate " << number_descriptors << " columns\n";
    return 0;
  }

  int fvalue = 0;
  if (write_zero_count_as_negative_one)
    fvalue = -1;

  for (int i = 0; i < pool_size; ++i)
  {
    const int initial_size = output.number_elements();

    std::fill_n(tmp, number_descriptors, fvalue);

    int nbits = 0;

    for (int j = 0; j < nfp; ++j)
    {
      if (fpnumber[j] < 1000)
        nbits += fill_output_array(pool[i][fpnumber[j]], col[j], tmp);
      else
        nbits += fill_output_array(pool[i].sparse_fingerprint(fpnumber[j]-1000), col[j], tmp);
    }

    if (discard_because_of_dash_y(pool[i].id(), nbits))
      output.resize_keep_storage(initial_size);
    else 
      write_descriptors(pool[i].id(), tmp, number_descriptors, output);
  }

  output.flush();

  return 1;
}

static int
scan_fixed_width_fingerprint(const IWDYFP & fp,
                             std::unordered_map<unsigned int, unsigned int> & fpcount,
                             int & nbits,
                             const IWString & fptag)
{
  unsigned int nb = fp.nbits();

  if (0 == nbits)
  {
    nbits = nb;
    if (verbose)
      cerr << "Fixed width fingerprint " << fptag << " has " << nbits << " bits\n";
  }
  else if (nb != static_cast<unsigned int>(nbits))
  {
    cerr << "Bit count mismatch, expected " << nbits << " got " << nb << endl;
    return 0;
  }

  for (unsigned int i = 0; i < nb; i++)
  {
    if (fp.is_set(i))
      fpcount[i]++;
  }

  return 1;
}

static int
scan_non_colliding_fingerprint(const Sparse_Fingerprint & fp,
                  IW_Hash_Map<unsigned int, unsigned int> & fpcount)
{
  int istart = 0;
  unsigned int zbit;
  int zcount;

  while (fp.next_bit_set(istart, zbit, zcount))
  {
    fpcount[zbit]++;
  }

  return 1;
}

static int
fingerprints_to_descriptors(const IW_General_Fingerprint * pool,
                            const int pool_size,
                            const resizable_array_p<IWString> & fptag,
                            const int * fpnumber,
                            IWString_and_File_Descriptor & output)
{
  const int nfp = fptag.number_elements();

  std::unordered_map<unsigned int, unsigned int> * fpcount = new std::unordered_map<unsigned int, unsigned int>[nfp]; std::unique_ptr<std::unordered_map<unsigned int, unsigned int>[]> free_fpcount(fpcount);

  int * nbits = new_int(nfp); std::unique_ptr<int[]> free_nbits(nbits);

  for (int i = 0; i < pool_size; ++i)
  {
    for (int j = 0; j < nfp; ++j)
    {
      if (fpnumber[j] < 1000)
        scan_fixed_width_fingerprint(pool[i][fpnumber[j]], fpcount[j], nbits[j], *fptag[j]);
      else
        scan_non_colliding_fingerprint(pool[i].sparse_fingerprint(fpnumber[j] - 1000), fpcount[j]);
    }
  }

  for (int i = 0; i < nfp; ++i)
  {
    if (0 == fpcount[i].size())
    {
      cerr << "Very strange, no fingerprints counted for fingerprint " << i << " across " << pool_size << " fingerprints\n";
      return 0;
    }

    if (verbose)
      cerr << "Found " << fpcount[i].size() << " bits set in " << *fptag[i] << " fingerprint\n";
  }

  IWString header("Name");
  header.resize(10 * fpcount[0].size());

  std::unordered_map<unsigned int, unsigned int> * col = new std::unordered_map<unsigned int, unsigned int>[nfp]; std::unique_ptr<std::unordered_map<unsigned int, unsigned int>[]> free_col(col);

  IWString * dprefix = new IWString[nfp]; std::unique_ptr<IWString[]> free_dprefix(dprefix);
  if (1 == nfp)
  {
    dprefix[0] = descriptor_prefix;
  }
  else
  {
    for (int i = 0; i < nfp; ++i)
    {
      dprefix[i] = descriptor_prefix;
      dprefix[i] << i << '.';
    }
  }

  int number_descriptors = 0;

//cerr << "Building header, lower_support_cutoff_number " << lower_support_cutoff_number << " upper_support_cutoff_number " << upper_support_cutoff_number << ", fold_to_constant_width " << fold_to_constant_width << endl;

  if (sort_by_frequency)
  {
    for (int i = 0; i < nfp; ++i)
    {
      do_sort_by_frequency(fpcount[i], col[i], dprefix[i], number_descriptors, header, pool_size);
    }
  }
  else if (fold_to_constant_width > 0)
    establish_header_and_columns_folded(fold_to_constant_width, fpcount, nfp, col, dprefix, number_descriptors, header);
  else if (bit_count_limits_set)
    establish_header_and_columns_with_thresholds(fpcount, nfp, col, dprefix, number_descriptors, header);
  else
  {
    for (int i = 0; i < nfp; ++i)
    {
      for (auto j : fpcount[i])
      {
        const unsigned int b = j.first;

        header << dprefix[i] << b;
        col[i][b] = number_descriptors;
        number_descriptors++;
      }
    }
  }

  if (0 == number_descriptors)
  {
    cerr << "Bad news, support levels " << lower_support_cutoff_number << " and " << upper_support_cutoff_number << " leave no fingerprints\n";
    return 0;
  }

  if (verbose)
    cerr << "Will produce " << number_descriptors << " descriptors\n";

  header << '\n';

  output << header;

  return fingerprints_to_descriptors(pool, pool_size, nfp, fpnumber, col, number_descriptors, output);
}

static IW_General_Fingerprint *
read_fingerprints(iwstring_data_source & input,
                  int & pool_size)
{
  pool_size = input.grep(identifier_tag);

  if (0 == pool_size)
  {
    cerr << "No fingerprints in input\n";
    return nullptr;
  }

  IW_General_Fingerprint * rc = new IW_General_Fingerprint[pool_size];

  IW_TDT tdt;
  for (int ndx = 0; tdt.next(input); ndx++)
  {
    int fatal;
    if (! rc[ndx].construct_from_tdt(tdt, fatal))
    {
      cerr << "Error building fingerprint\n";
      return 0;
    }
  }

  return rc;
}

static IW_General_Fingerprint *
read_fingerprints(const char * fname,
                  int & pool_size)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "read_fingerprints:cannot open '" << fname << "'\n";
    return nullptr;
  }

  return read_fingerprints(input, pool_size);
}

static int
fingerprints_to_descriptors(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vF:x:X:P:sbmy:Y:un:N:qo:d:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise fingerprints\n";
      return 1;
    }
  }

  if (cl.option_present('o'))
  {
    IWString o = cl.string_value('o');

    if (! char_name_to_char(o))
    {
      cerr << "unrecoginsed output file separator '" << o << "'\n";
      return 1;
    }

    output_separator= o[0];

    if (verbose)
      cerr << "Output separator '" << output_separator << "'\n";
  }

  resizable_array_p<IWString> fptag;

#ifdef NOT_DONE_ANY_MORE
  if (cl.option_present('F'))
  {
    IWString f;
    for (int i = 0; cl.value('F', f, i); ++i)
    {
      int j = 0;
      IWString token;

      while (f.nextword(token, j, ','))
      {
        if (! token.ends_with('<'))
          token += '<';
        fptag.add(new IWString(token));
        if (verbose)
          cerr << "Will process fingerprint '" << token << "'\n";
      }
    }
  }
#endif

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

  if (cl.option_present('x'))
  {
    if (! cl.value ('x', lower_population_cutoff) || lower_population_cutoff < 0.0 || lower_population_cutoff >= 1.0)
    {
      cerr << "The lower population cutoff must be a valid fraction\n";
      usage(7);
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

    bit_count_limits_set = 1;
  }

  if (cl.option_present('X'))
  {
    if (! cl.value('X', upper_population_cutoff) || upper_population_cutoff < 0.0 || upper_population_cutoff >= 1.0)
    {
      cerr << "The upper population cutoff must be a valid fraction\n";
      usage(7);
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

    if (0 != lower_support_cutoff_number && lower_support_cutoff_number > upper_support_cutoff_number)
    {
      cerr << "Inconsistent -n (" << lower_support_cutoff_number << ") and -N (" << upper_support_cutoff_number << ") option values\n";
      return 4;
    }

    bit_count_limits_set = 1;
  }

  if (0 == upper_population_cutoff)
    ;
  else if (lower_population_cutoff > upper_population_cutoff)
  {
    cerr << "Impossible lower (" << lower_population_cutoff << ") and upper (" << upper_population_cutoff << ") specification\n";
    return 9;
  }

  if (cl.option_present('I'))
  {
    descriptor_prefix = cl.string_value('I');

    if (verbose)
      cerr << "Descriptors will have the prefix '" << descriptor_prefix << "'\n";

    if (! descriptor_prefix.starts_with(output_separator))
      descriptor_prefix.insert_before(0, output_separator);
  }
  else
   descriptor_prefix << output_separator << 'F';

  if (cl.option_present('u'))
  {
    gsub_multi_token_names = 1;
    if (verbose)
      cerr << "Will gsub to _ any multi-token names\n";
  }

  if (cl.option_present('s'))
  {
    sort_by_frequency = 1;

    if (verbose)
      cerr << "Bits will be output in order of frequency of occurrence in the set\n";
  }

#ifdef DEPRECATE_DASH_B_OPTION
  if (cl.option_present('b'))
  {
    write_actual_bit_numbers = 1;

    if (verbose)
      cerr << "Will write bit numbers rather than sequential names\n";
  }
#endif

  if (1)
  {
    IWString o(output_separator);
    iwdigits.set_leading_string(o);
    iwdigits.initialise(100);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, cannot do multiple input files, try cat'ing them together\n";
    return 1;
  }

  int pool_size = 0;
  IW_General_Fingerprint * pool = read_fingerprints(cl[0], pool_size);

  if (nullptr == pool || 0 == pool_size)
  {
    cerr << "Cannot read fingerprints '" << cl[0] << "'\n";
    return 1;
  }

  if (verbose)
    cerr << "Read " << pool_size << " fingerprints from '" << cl[0] << "'\n";

  if (0 == fptag.number_elements())    // consume all of them
  {
    for (int i = 0; i < number_fingerprints(); ++i)
    {
      fptag.add(new IWString(fixed_fingerprint_tag(i)));
    }
    for (int i = 0; i < number_sparse_fingerprints(); ++i)
    {
      fptag.add(new IWString(sparse_fingerprint_tag(i)));
    }
  }

  if (fptag.number_elements() > 1 && sort_by_frequency)
  {
    cerr << "WARNING, sorting by frequency with multiple fingerprints does not sort across fingeprints!\n";
  }

  int * fpnumber = new_int(fptag.number_elements(), -1); std::unique_ptr<int[]> free_fpnumber(fpnumber);

  for (int i = 0; i < fptag.number_elements(); ++i)
  {
    const IWString & t = *fptag[i];

    for (int j = 0; j < number_fingerprints(); ++j)
    {
      if (t == fixed_fingerprint_tag(j))
      {
        fpnumber[i] = j;
        break;
      }
    }

    if (fpnumber[i] >= 0)
      continue;

    for (int j = 0; j < number_sparse_fingerprints(); ++j)
    {
      if (t == sparse_fingerprint_tag(j))
      {
        fpnumber[i] = 1000 + j;
        break;
      }
    }

    if (fpnumber[i] < 0)
    {
      cerr << "No fingerprint match for '" << t << "'\n";
      return 1;
    }
  }

// Convert any ratios specified for support levels to numbers

  if (zero != lower_population_cutoff)
  {
    lower_support_cutoff_number = static_cast<unsigned int>(lower_population_cutoff * static_cast<double>(pool_size) + 0.4999);
    if (0 == lower_support_cutoff_number)
      lower_support_cutoff_number = 1;

    if (verbose)
      cerr << "Will discard bits that occur in fewer than " << lower_population_cutoff << " molecules\n";

    bit_count_limits_set = 1;
  }

  if (zero != upper_population_cutoff)
  {
    upper_support_cutoff_number = static_cast<unsigned int>(upper_population_cutoff * static_cast<double>(pool_size) + 0.4999);
    assert (upper_support_cutoff_number > 1);
    if (static_cast<unsigned int>(pool_size) == upper_support_cutoff_number)
      upper_support_cutoff_number = pool_size - 1;

    if (verbose)
      cerr << "Will discard bits that occur in more than " << upper_support_cutoff_number << " molecules\n";

    bit_count_limits_set = 1;
  }

  if (cl.option_present('d'))
  {
    if (cl.option_present('b'))
    {
      cerr << "Sorry, the -b option cannot be used with the -d option\n";
      usage(1);
    }

    if (sort_by_frequency)
    {
      cerr << "Folded fingerprints and sort by frequency cannot both be used, see Ian\n";
      return 1;
    }

    if (! cl.value('d', fold_to_constant_width) || fold_to_constant_width < 1)
    {
      cerr << "The fold to constant width option (-d) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will fold fingerprints to a constant width of " << fold_to_constant_width << " bits\n";

    collision_count = new_int(fold_to_constant_width);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! fingerprints_to_descriptors(pool, pool_size, fptag, fpnumber, output))
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

  delete [] pool;

  if (nullptr != collision_count)
  {
    if (verbose)
    {
      int bits_with_collisions = 0;
      int max_collisions = 0;

      for (int i = 0; i < fold_to_constant_width; ++i)
      {
        if (0 == collision_count[i])
          continue;

        bits_with_collisions++;
        if (verbose > 1)
          cerr << " bit " << i << ' ' << collision_count[i] << " collisions\n";

        if (collision_count[i] > max_collisions)
          max_collisions = collision_count[i];
      }
      cerr << bits_with_collisions << " of " << fold_to_constant_width << " bits with collisions (fraction " << (static_cast<float>(bits_with_collisions)/static_cast<float>(fold_to_constant_width)) << ") max collision count " << max_collisions << endl;
    }
    delete [] collision_count;
  }

  return rc;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = fingerprints_to_descriptors(argc, argv);

  return rc;
}
