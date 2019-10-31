/*
  We often need to extract random records from a file
*/

#include <memory>
#include <random>

#ifdef FOOBAR
#if (__GNUC__ >= 3)
#include <ext/hash_map>
using namespace __gnu_cxx;
#else
#include <hash_map>
#endif
#endif

using namespace std;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "iw_stl_hash_map.h"
#include "iwcrex.h"

#include "cmdline.h"
#include "iwmmap.h"
#include "accumulator.h"

static int verbose = 0;

static unsigned int records_to_select = 1;

static unsigned int max_iterations = 0;

static int header_records_to_keep = 0;

static int columns_in_input = 0;

static IWString missing_value('.');

static int with_replacement = 0;

static int unbiased_selection = 0;

static int sorted_record_order = 1;

/*
  When doing stratified sampling, we need to count the number of non-zero
  occurrences in each column
*/

static int * non_zero_samples = NULL;

/*
  We may want to also choose a number of records that follow the selected records
*/

static int extra_records_to_select = 0;

static int sequential_scan = 0;

static IW_Regular_Expression rx;

static IW_Regular_Expression end_of_record_rx;

static resizable_array<int> specific_records_to_select;

static int like_grep = 0;

static int write_records_as_they_are_selected = 0;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Extracts random records from a file\n";

  cerr << " -R <pattern>    selected records must match <pattern>\n";
  cerr << " -n <number>     the number of records to select (default " << records_to_select << ")\n";
  cerr << " -e <number>     extra records to also select after each randomly chosen record\n";
  cerr << " -j              is a descriptor file - write the first record (same as -J 1)\n";
  cerr << " -J <number>     number of header records to preserve\n";
  cerr << " -a <number>     number of attempts per record to find a random record (default 10)\n";
  cerr << " -r              sample with replacement - may get the same record multiple times\n";
  cerr << " -u              unbiased selection - every record has equal probability of being\n";
  cerr << "                 chosen. Must scan the input file once however\n";
  cerr << " -f              fast sorted selection. Works best for large files\n";
  cerr << " -z              by default records are sorted in natural order, do random order instead\n";
  cerr << " -E <pattern>    end of record pattern - for dealing with multi-line records\n"; 
  cerr << " -s <number>     specific record(s) to select (first is 1)\n";
  cerr << " -g              work like grep - sequentially scan the file (only appropriate with -E option)\n";
  cerr << " -q              sequential scan - good for large files\n";
  cerr << " -p <prob>       read the file line at a time. Write lines with probability <prob>\n";
  cerr << " -v              verbose output\n";

  exit(rc);
}

static int
random_records_sequentual_scan(IWString_Data_Source_MMAP & input,
                                const double p,
                                IWString_and_File_Descriptor & output)
{
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_real_distribution<double> u(0.0, 1.0);

  const_IWSubstring buffer;

  unsigned int records_selected = 0;

  while (input.next_record(buffer) && records_selected < records_to_select)
  {
    if (u(rng) >= p)
      continue;

    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
    records_selected++;

    if (extra_records_to_select)
    {
      for (int i = 0; i < extra_records_to_select; ++i)
      {
        if (! input.next_record(buffer))
        {
          cerr << "EOF getting " << extra_records_to_select << " extra records\n";
          return 0;
        }
        output << buffer << '\n';
        output.write_if_buffer_holds_more_than(4096);
      }
    }
  }

  if (records_selected < records_to_select)
    cerr << "random_records_sequentual_scan:requested " << records_to_select << " but only wrote " << records_selected << endl;

  return 1;
}

class streampos_hash_fn : public hash<off_t>
{
  private:
  public:
    size_t operator()(const off_t &) const;
};

size_t
streampos_hash_fn::operator()(const off_t & s) const
{
  if (8 == sizeof(off_t))
    return static_cast<size_t>(s);
  else
    return static_cast<size_t>(s);
}

static int
choose_records_sequtentially(IWString_Data_Source_MMAP & input,
                              off_t start_position,
                              IWString_and_File_Descriptor & output)
{
  off_t fsize = input.file_size() - input.tellg();

  off_t chunk_size = fsize / records_to_select - 2;

  if (chunk_size < 2)
  {
    cerr << "Cannot sequentially select " << records_to_select << " from file of size " << fsize << endl;
    return 0;
  }

  if (verbose)
    cerr << "For selecting " << records_to_select << " from file of size " << fsize << " initial chunk size " << chunk_size << endl;

  unsigned int nsel = 0;
  const_IWSubstring buffer;   // scope here for efficiency

  Accumulator_Int<int> record_length;

  while (1)
  {
    if (! input.seekg(chunk_size, SEEK_CUR))
    {
      cerr << "Cannot seek " << chunk_size << " from " << input.tellg() << endl;
      break;
    }

    if (! input.next_record(buffer))
    {
      cerr << "Could not get record after " << input.tellg() << endl;
      break;
    }

    if (! input.next_record(buffer))
    {
      cerr << "Cannot fetch record after " << input.tellg() << endl;
      break;
    }

    record_length.extra(buffer.length());

    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(32768);

    nsel++;

    if (nsel >= records_to_select)
      break;

    int ave_record_length = static_cast<int>(record_length.average_if_available_minval_if_not()) + 1;

    chunk_size = (fsize - input.tellg()) / (records_to_select - nsel) - ave_record_length;

//  cerr << "Offset now at " << input.tellg() << ", chunk size " << chunk_size << endl;
  }

  if (nsel < records_to_select)
    cerr << "Warning, only selected " << nsel << " records\n";
    
  return 1;
}

/*
  Does SELECTED already contain P. If not, add it.
*/

static int
can_be_output(IW_Hash_Map<off_t, int, streampos_hash_fn> & selected,
               const off_t & p,
               const off_t & start_position)
{
  if (p < start_position)
    return 0;

  IW_Hash_Map<off_t, int, streampos_hash_fn>::iterator f = selected.find(p);

  if (f == selected.end())    // new record
  {
    selected[p] = 1;
    return 1;
  }

  if (with_replacement)
  {
    selected[p]++;
    return 1;
  }

  return 0;
}

/*
  Since a "record" is determined by end_of_record_rx, we need a special function
  for examining the next record.
*/

static int
fetch_next_record(IWString_Data_Source_MMAP & input,
                   int & rx_matched)
{
  if (rx.active())
    rx_matched = 0;
  else
    rx_matched = 1;

  const_IWSubstring buffer;

  int lines_this_record = 0;   

  while (input.next_record(buffer))
  {
    lines_this_record++;

    if (rx.active() && 0 == rx_matched && rx.matches(buffer))
      rx_matched = 1;

    if (! end_of_record_rx.active())    // one line per record
      return lines_this_record;

    if (end_of_record_rx.matches(buffer))
      return lines_this_record;
  }

  return 0;
}

static int
choose_a_record(IWString_Data_Source_MMAP & input,
                 off_t start_position,
                 off_t & p,
                 int & rx_matched,
                 int & fatal)
{
  const off_t file_size = input.file_size() - static_cast<off_t>(3);

  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<off_t> u(start_position, file_size);

  while(1)
  {
    auto o = u(rng);

    if (! input.seekg(o))
    {
      cerr << "Yipes, cannot seek to " << o << " in file\n";
      fatal = 1;
      return 0;
    }
  
    if (o > 0)     // most probably in the middle of a record. Skip to end of record
    {
      int notused;
      if (! fetch_next_record(input, notused))
      {
        cerr << "Very strange, cannot read from " << o << endl;
        fatal = 1;
        return 0;
      }

      o = input.tellg();
    }
  
    if (fetch_next_record(input, rx_matched))   // will fail if we seek'd into the last record
    {
      p = o;
      return 1;
    }

//  Must have seek'd into last record. Since we have very little change of getting the first record otherwise, do it here

    p = start_position;

    input.seekg(p);

    return fetch_next_record(input, rx_matched);
  }
}


static int
get_extra_records(IWString_Data_Source_MMAP & input,
                   int & ndx,
                   off_t * ov)
{
  off_t o = input.tellg();

  int rc = 0;
  for (int i = 0; i < extra_records_to_select; i++)
  {
    int notused;
    if (! fetch_next_record(input, notused))
    {
//    cerr << "Only got " << i << " of " << extra_records_to_select << " extra records\n";
      break;
    }

    ov[ndx] = o;
    ndx++;

    o = input.tellg();
    rc++;
  }

  if (rc == extra_records_to_select)
    return 1;

  ndx -= rc;

  return 0;
}

static int
streampos_comparitor(const void * pp1, const void * pp2)
{
//streampos p1 = * (reinterpret_cast<streampos *> (pp1));
//streampos p2 = * (reinterpret_cast<streampos *> (pp2));
  off_t p1 = * ((off_t *) pp1);
  off_t p2 = * ((off_t *) pp2);

  if (p1 < p2)
    return -1;

  if (p1 > p2)
    return 1;

  return 0;
}

static int
do_write_records_as_they_are_selected(IWString_Data_Source_MMAP & input,
                                       const int average_record_length,
                                       IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;
  for (auto i = 0; i < header_records_to_keep; ++i)
  {
    if (! input.next_record(buffer))
    {
      cerr << "Cannot echo header record " << i << endl;
      return 0;
    }

    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);
  }

  off_t start_position = input.tellg();

  off_t file_size = input.file_size();

  int buffer_size = records_to_select * 4;

  off_t * p = new off_t[buffer_size];

  std::random_device rd;

  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<off_t> u(start_position, file_size - average_record_length);

  for (auto i = 0; i < buffer_size; ++i)
  {
    p[i] = u(rng);
  }

  qsort(p, buffer_size, sizeof(off_t), streampos_comparitor);    // we always sort

// get rid of dups

  int ndx = 0;
  for (auto i = 1; i < buffer_size; ++i)
  {
    if (p[i] != p[ndx])
    {
      ndx++;
      p[ndx] = p[i];
    }
  }
 
  buffer_size = ndx;

  unsigned int rc = 0;
  for (auto i = 0; i < buffer_size; ++i)
  {
    auto o = p[i];

    if (o > average_record_length)
    {
      o -= average_record_length;

      if (! input.seekg(p[i]))
      {
        cerr << "Huh, within file if size " << file_size << " cannot seek to " << p[i] << endl;
        return 0;
      }

      if (! input.next_record(buffer))
      {
        cerr << "Huh, could not get next record, position " << p[i] << " file size " << file_size << endl;
        return 0;
      }
    }
    else
      o = start_position;

    if (! input.next_record(buffer))
    {
      cerr << "Huh, cannot fetch record near " << p[i] << " file size " << file_size << endl;
      return 0;
    }

    if (rx.active() && ! rx.matches(buffer))
      continue;

    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);

    rc++;

    for (auto j = 0; j < extra_records_to_select; ++j)
    {
      if (! input.next_record(buffer))
      {
        cerr << "Cannot fetch extra record '" << j << " from near " << p[i] << ", file size " << file_size << endl;
        break;
      }

      output << buffer << '\n';
      output.write_if_buffer_holds_more_than(8192);
    }

    if (rc >= records_to_select)
      return 1;

    o = input.tellg();

    while (p[i] <= o)
      ++i;
  }

  if (rc < records_to_select)
    cerr << "Warning, only retrieved " << rc << " of " << records_to_select << " requested\n";

  return rc;
}

static int
choose_like_grep(IWString_Data_Source_MMAP & input,
                  off_t & start_position,
                  off_t * ov)
{
  int records_found = 0;

  int ndx = 0;

  off_t o = input.tellg();

  int rx_matched;
  while (fetch_next_record(input, rx_matched))
  {
    if (rx_matched)
    {
      ov[ndx] = o;
      ndx++;
      records_found++;
      if (extra_records_to_select)
      {
        if (! get_extra_records(input, ndx, ov))
          records_found--;
      }
    }

    o = input.tellg();
  }

  return records_found;
}

static int
choose_specific_records(IWString_Data_Source_MMAP & input,
                         off_t & start_position,
                         off_t * ov)
{
  int records_read = 0;
  int records_found = 0;   // index into OV array

  off_t o = input.tellg();

  int rx_matched;

  while (fetch_next_record(input, rx_matched))
  {
    records_read++;

    if (! rx_matched)
      cerr << "Warning, record " << records_read << " does not match regular expression\n";

    if (specific_records_to_select.contains(records_read))
    {
      ov[records_found] = o;
      records_found++;
      if (records_found == specific_records_to_select.number_elements())
        return records_found;
    }

    o = input.tellg();
  }

  return 0;
}

#ifdef IWSGI
template class resizable_array<off_t>;
template class resizable_array_base<off_t>;
#endif

static int
choose_the_records_uniformly(IWString_Data_Source_MMAP & input,
                              off_t & start_position,
                              off_t * ov)
{
  resizable_array<off_t> o;    // the offset for every record in the file
  o.resize(16000);

  o.add(input.tellg());

  int rx_matched;
  while (fetch_next_record(input, rx_matched))
  {
    o.add(input.tellg());
  }

  o.pop();     // the last entry in the array should be the file size

  if (verbose)
    cerr << "Input contains " << o.number_elements() << " records\n";

  if (static_cast<unsigned int>(o.number_elements()) <= records_to_select && 0 == with_replacement)    // should be more careful about is_descriptor_file, where one less record is to be selected
  {
    cerr << "Only " << o.number_elements() << " records in the file, requested " << records_to_select << ", fetching all records\n";
    for (int i = 0; i < o.number_elements() - 1; i++)
    {
      ov[i] = o[i];
    }

    return o.number_elements();
  }

  IW_Hash_Map<off_t, int, streampos_hash_fn> selected;

  unsigned int number_selected = 0;
  unsigned int iterations = 0;
  unsigned int ndx = 0;

  int istart;

  if (header_records_to_keep > 0)
    istart = header_records_to_keep;
  else
    istart = 0;

  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<int> u(istart, o.number_elements() - 1);

  while (number_selected < records_to_select && iterations < max_iterations)
  {
    int j = u(rng);

    iterations++;

    if (! can_be_output(selected, o[j], start_position))
      continue;

//  Should do something about the regular expres

    ov[ndx] = o[j];
    ndx++;

    number_selected++;

    if (0 == extra_records_to_select)
      continue;

    input.seekg(o[j]);    // go back to start position

    int extra_records_selected = 0;

    for (int i = 0; i < extra_records_to_select; i++)
    {
      if (j + i + 1 >= o.number_elements() - 1)
        break;

      off_t p = o[j + i + 1];
      selected[p]++;
      ov[ndx] = p;
      ndx++;
      extra_records_selected++;
    }

    if (extra_records_selected < extra_records_to_select)
    {
      number_selected--;
      ndx -= (extra_records_selected + 1);
    }
  }

  if (verbose)
    cerr << "Selected " << selected.size() << " of " << o.number_elements() << " records for output\n";

  return ndx;
}


static int
choose_the_records(IWString_Data_Source_MMAP & input,
                    off_t start_position,
                    off_t * ov)
{
//cerr << "Choosing between " << start_position << " and " << input.file_size() << endl;

  IW_Hash_Map<off_t, int, streampos_hash_fn> selected;

  unsigned int iterations = 0;
  unsigned int number_selected = 0;
  unsigned int ndx = 0;
  while (number_selected < records_to_select && iterations < max_iterations)
  {
    iterations++;

//  cerr << "Iteration " << iterations << endl;

    off_t p;
    int rx_matched;
    int fatal = 0;
    if (! choose_a_record(input, start_position, p, rx_matched, fatal))
    {
      if (fatal)
        return 0;

      continue;
    }

//  cerr << "Checking record at " << p << endl;

    if (rx.active() && ! rx_matched)
      continue;

    if (! can_be_output(selected, p, start_position))
      continue;

    ov[ndx] = p;
    ndx++;
    number_selected++;
    int extra_records_selected = 0;
    for (int i = 0; i < extra_records_to_select; i++)
    {
      p = input.tellg();

      int notused;
      if (! fetch_next_record(input, notused))
        break;

      selected[p]++;     // may already be in SELECTED
      ov[ndx] = p;
      ndx++;
      extra_records_selected++;
    }

    if (extra_records_selected < extra_records_to_select)
    {
      number_selected -= 1;
      ndx -= (extra_records_selected + 1);
    }
  }

  return ndx;
}

/*static int
echo_header_record(IWString_Data_Source_MMAP & input, 
                    IWString_and_File_Descriptor & output)
{
  return input.echo_records(output, 1);
}*/

static int
echo_the_records(IWString_Data_Source_MMAP & input,
                  const off_t * p,
                  int records_selected,
                  IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < records_selected && output.good(); i++)
  {
    if (! input.seekg(p[i]))
    {
      cerr << "echo_the_records: huh, cannot seek to " << p[i] << endl;
      return 0;
    }

    const_IWSubstring buffer;

    int records_written = 0;

    while (input.next_record(buffer))
    {
      records_written++;

      output << buffer << '\n';

      output.write_if_buffer_holds_more_than(32768);

      if (! end_of_record_rx.active())
        break;

      if (end_of_record_rx.matches(buffer))
        break;
    }

    if (0 == records_written)
    {
      cerr << "echo_the_records: huh, cannot read record from " << p[i] << endl;
      return 0;
    }

  }

  return output.good();
}

static int
echo_file(IWString_Data_Source_MMAP & input,
           IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer) && output.good())
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(32768);
  }

  return output.good();
}

static int
_at_least_X_records_remaining(IWString_Data_Source_MMAP & input,
                              int records_to_select,
                              int & average_record_length)
{

  const_IWSubstring buffer;

  int bytes_read = 0;

  for (auto i = 0; i < records_to_select; ++i)
  {
    if (! input.next_record(buffer))
      return 0;

    bytes_read += buffer.length() + 1;
  }

  average_record_length = bytes_read / records_to_select + 1;

  return 1;
}

static int
at_least_X_records_remaining(IWString_Data_Source_MMAP & input,
                              int records_to_select,
                              int & average_record_length)
{
  const auto sbegin = input.tellg();

  const auto rc = _at_least_X_records_remaining(input, records_to_select, average_record_length);

  input.seekg(sbegin);

  return rc;
}

static int
random_records(const char * fname, 
                const double prob,
                IWString_and_File_Descriptor & output)
{
  assert (records_to_select > 0);

  IWString_Data_Source_MMAP input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open input file '" << fname << "'\n";
    return 0;
  }

// If the file is a descriptor file, we need to make sure we never select the first record again

  off_t start_position;

  if (header_records_to_keep > 0)
  {
    if (! input.echo_records(output, header_records_to_keep))
    {
      cerr << "Cannot echo " << header_records_to_keep << " header records\n";
      return 0;
    }

    start_position = input.tellg();
  }
  else
  {
    start_position = static_cast<off_t>(0);
  }

  if (prob > 0.0)
    return random_records_sequentual_scan(input, prob, output);

  int average_record_length;

  if (! at_least_X_records_remaining(input, records_to_select, average_record_length))
  {
    cerr << "Too few records in the input file\n";
    return echo_file(input, output);
  }

/*if (! input.at_least_X_records_remaining (records_to_select))
  {
    cerr << "Too few records in the input file\n";
    return echo_file (input, output);
  }*/

  if (verbose)
    cerr << "Input file has " << input.file_size() << " bytes\n";

  if (sequential_scan)
    return choose_records_sequtentially(input, start_position, output);

  if (write_records_as_they_are_selected)
    return do_write_records_as_they_are_selected(input, average_record_length, output);

  off_t * p = new off_t[records_to_select + records_to_select * extra_records_to_select]; std::unique_ptr<off_t[]> free_p(p);

  int records_selected;

  if (specific_records_to_select.number_elements())
    records_selected = choose_specific_records(input, start_position, p);
  else if (like_grep)
    records_selected = choose_like_grep(input, start_position, p);
  else if (unbiased_selection)
    records_selected = choose_the_records_uniformly(input, start_position, p);
  else
    records_selected = choose_the_records(input, start_position, p);

  if (0 == records_selected)
  {
    cerr << "Cannot choose " << records_to_select;
    if (0 == with_replacement)
      cerr << " different";
    cerr << " records from '" << fname << "'\n";

    return 0;
  }

  if (verbose)
    cerr << "Selected " << records_selected << " records for output\n";

  if (records_selected > 1 && sorted_record_order)
    qsort(p, records_selected, sizeof(off_t), streampos_comparitor);

  if (verbose > 2)
  {
    for (int i = 0; i < records_selected; i++)
    {
      cerr << "Will echo " << p[i] << endl;
    }
  }

  return echo_the_records(input, p, records_selected, output);
}

static int
random_records (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vn:jJ:R:e:a:ruzs:E:gqfp:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    if (! cl.value('n', records_to_select) || records_to_select < 1)
    {
      cerr << "The records to select option (-n) must be followed by a whole positive number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will select " << records_to_select << " records\n";
  }
  else
  {
    records_to_select = 1;

    if (verbose)
      cerr << "By default will select 1 record\n";
  }

  if (cl.option_present('j'))
  {
    header_records_to_keep = 1;

    if (verbose)
      cerr << "Will treat as a descriptor file\n";
  }

  if (cl.option_present('J'))
  {
    if (! cl.value('J', header_records_to_keep) || header_records_to_keep < 0)
    {
      cerr << "The header records to keep option (-J) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will retain the first " << header_records_to_keep << " records\n";
  }

  if (cl.option_present('r'))
  {
    with_replacement = 1;

    if (verbose)
      cerr << "Will select with replacement\n";
  }

  if (cl.option_present('R'))
  {
    IWString r = cl.string_value('R');

    if (! rx.set_pattern(r))
    {
      cerr << "Invalid regular expression '" << r << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Will only retrieve records matching '" << rx.source() << "'\n";
  }

  if (cl.option_present('E'))
  {
    IWString r = cl.string_value('E');

    if (! end_of_record_rx.set_pattern(r))
    {
      cerr << "Invalid end of record regular expression '" << r << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "End of record regular expression'" << end_of_record_rx.source() << "'\n";
  }

  if (cl.option_present('e'))
  {
    if (! cl.value('e', extra_records_to_select) || extra_records_to_select < 0)
    {
      cerr << "The extra records to choose option (-e) must be followed by a whole number\n";
      usage(6);
    }

    if (verbose)
      cerr << "Will select " << extra_records_to_select << " extra records along with each randomly chosen record\n";
  }

  if (cl.option_present('u'))
  {
    unbiased_selection = 1;

    if (verbose)
      cerr << "Strictly uniform record probability\n";
  }

  if (cl.option_present('z'))
  {
    sorted_record_order = 0;

    if (verbose)
      cerr << "Records will be in random order on output\n";
  }

  if (cl.option_present('f'))
  {
    write_records_as_they_are_selected = 1;

    if (verbose)
      cerr << "Fast sequential record selection\n";
  }

  double p = 0.0;

  if (cl.option_present('p'))
  {
    if (! cl.value('p', p) || p <= 0.0 || p >= 1.0)
    {
      cerr << "The probability of record selection (-p) must be a valid probability\n";
      usage(1);
    }

    if (verbose)
      cerr << "Sequential reading of file, write records with probability " << p << endl;
  }

  if (cl.option_present('a'))
  {
    int a;
    if (! cl.value('a', a) || a < 1)
    {
      cerr << "The maximum iterations per record flag (-a) must be followed by a whole positive number\n";
      usage(6);
    }

    max_iterations = a * records_to_select;

    if (verbose)
      cerr << "A maximum of " << a << " attempts to identify each record will be performed\n";
  }
  else if (p > 0.0)
    ;
  else
  {
    max_iterations = records_to_select * 10;

    if (verbose)
      cerr << "By default, max attempts per record set to 10\n";
  }

  if (cl.option_present('s'))
  {
    if (cl.option_present('n') || cl.option_present('g'))
    {
      cerr << "The -n option is mutually inconsistent with the -s option\n";
      usage(5);
    }

    int i = 0;
    const_IWSubstring s;
    while (cl.value('s', s, i++))
    {
      int j;
      if (! s.numeric_value(j) || j < 1)
      {
        cerr << "INvalid record to select '" << s << "'\n";
        return 8;
      }

      specific_records_to_select.add(j);
    }

    records_to_select = specific_records_to_select.number_elements();
  }

  if (cl.option_present('g'))
  {
    if (! rx.active())
    {
      cerr << "In order to work like grep, there must be a regular expression\n";
      usage(5);
    }

    like_grep = 1;

    if (verbose)
      cerr << "Will work like grep\n";
  }

  if (cl.option_present('q'))
  {
    sequential_scan = 1;

    if (verbose)
      cerr << "Will work as a sequential scan\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Must specify input file on command line\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, only works on a single file\n";
    usage(3);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  for (int i = 0; i < cl.number_elements(); i++)    // but only works for 1 file
  {
    if (! random_records(cl[i], p, output))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = random_records(argc, argv);

  return rc;
}
