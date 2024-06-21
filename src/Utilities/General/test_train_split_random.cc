/*
  Splits a set of molecules into test and train splits, while
  ensuring uniform distribution
*/

#include <stdlib.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <random>
#include <valarray>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

using std::cerr;
using std::endl;

const char * prog_name = NULL;

static int verbose = 0;

static float fraction_in_training_set = static_cast<float>(0.5);

static int number_splits = 1;

/*
  After how many splits will we start refining them
*/

static int refine_splits = 0;

/*
  And how many items in each split will we check for adjustment
*/

static int refinement_steps = 0;

static int total_refinements_done = 0;

static IWString output_stem("ttRsplit");

static IWString suffix;

static IWString identifier_tag("PCN<");

static int write_training_set = 1;

static int write_test_set = 1;

#define INPUT_TYPE_ID 1
#define INPUT_TYPE_SMI 2
#define INPUT_TYPE_TDT 3
#define INPUT_TYPE_DM 4
#define INPUT_TYPE_DSC 5

static int input_type = 0;

static int write_header_to_split_files = 0;

static void
usage (int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Splits a set of molecules into test and train\n";
  cerr << " -p <pct>       percent of molecules in the training set (default 50)\n";
  cerr << " -n <number>    number of splits to create\n";
  cerr << " -S <stem>      write results to files starting with <stem> (default '" << output_stem << ")\n";
  cerr << " -k <...>       what kind of input file\n";
  cerr << " -k id          file of ID's, one per line (default)\n";
  cerr << " -k smi         smiles file\n";
  cerr << " -k tdt         TDT file\n";
  cerr << " -k dm          distance matrix\n";
  cerr << " -k dsc         descriptor file\n";
  cerr << " -w test        write the test     set identifiers\n";
  cerr << " -w train       write the training set identifiers\n";
  cerr << " -h             write a header record to each of the split files\n";
  cerr << " -F <suffix>    append <suffix> to all files produced\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit (rc);
}

class TTItem
{
  private:
    IWString _id;

    int _times_in_training_set;

    int _in_training_set;

    static std::random_device _rd;
    static std::default_random_engine _rng;

    std::uniform_real_distribution<double> _u;

//  private functions

    int _volunteer_lower_probability (int splits_remaining, int times_should_have_been_chosen);
    int _volunteer_increased_probability (int splits_remaining, int times_should_have_been_chosen);

  public:
    TTItem();

    const IWString & id() const { return _id;}
    void set_id (const const_IWSubstring & s) { _id = s;}

    int times_in_training_set() const { return _times_in_training_set;}

    int in_training_set() const { return _in_training_set;}
    void set_in_training_set (int s);

    int volunteer (int, int);

    void switch_training_set_membership();
};

std::default_random_engine TTItem::_rng;

TTItem::TTItem() : _u(0.0, 1.0)
{
  _times_in_training_set = 0;

  _in_training_set = 0;

  return;
}

void
TTItem::set_in_training_set (int s)
{
  _in_training_set = s;

  if (_in_training_set)
    _times_in_training_set++;

  return;
}

/*
  When choosing whether or not to volunteer, we need to always try
  to move towards the average. The degree of movement will be a
  function of the "urgency" of getting there.
*/

int
TTItem::volunteer (int splits_remaining,
                   int times_should_have_been_chosen)
{
  int rc;

  if (0 == _times_in_training_set || 0 == times_should_have_been_chosen)
    rc = 1;
  else if (_times_in_training_set > times_should_have_been_chosen)
    rc = _volunteer_lower_probability(splits_remaining, times_should_have_been_chosen);
  else
    rc = _volunteer_increased_probability(splits_remaining, times_should_have_been_chosen);

  if (rc)
    set_in_training_set(1);

  return rc;
}

int
TTItem::_volunteer_lower_probability (int splits_remaining,
                                      int times_should_have_been_chosen)
{
  assert (_times_in_training_set > times_should_have_been_chosen);

  double p = static_cast<double>(times_should_have_been_chosen) / static_cast<double>(_times_in_training_set);

  if (_u(_rng) < p)
    return 1;
  else
    return 0;
}

int
TTItem::_volunteer_increased_probability (int splits_remaining,
                                          int times_should_have_been_chosen)
{
  assert (_times_in_training_set <= times_should_have_been_chosen);

  double p = static_cast<double>(_times_in_training_set) / static_cast<double>(times_should_have_been_chosen);

  if (_u(_rng) < p)
    return 1;
  else
    return 0;
}


void
TTItem::switch_training_set_membership()
{
  if (_in_training_set)
  {
    _in_training_set = 0;
    _times_in_training_set--;
  }
  else
  {
    _in_training_set = 1;
    _times_in_training_set++;
  }

  return;
}

static int number_items = 0;

TTItem * ttitem = NULL;

static void
set_training_set_membership (int s)
{
  for (int i = 0; i < number_items; i++)
  {
    ttitem[i].set_in_training_set(s);
  }
  
  return;
}

class Times_in_Training_Set_Comparitor
{
  private:
  public:
    int operator() (const TTItem *, const TTItem *) const;
};

int
Times_in_Training_Set_Comparitor::operator() (const TTItem * pt1,
                                               const TTItem * pt2) const
{
  int t1 = pt1->times_in_training_set();
  int t2 = pt2->times_in_training_set();

  if (t1 < t2)
    return -1;
  if (t1 > t2)
    return 1;
  return 0;
}

static int
adjust_split_toward_average()
{
  TTItem ** tmp = new TTItem *[number_items]; std::unique_ptr<TTItem *[]> free_tmp(tmp);

  for (int i = 0; i < number_items; i++)
  {
    tmp[i] = &ttitem[i];
  }

  Times_in_Training_Set_Comparitor idtsc;

  iwqsort(tmp, number_items, idtsc);

  assert (tmp[0]->times_in_training_set() <= tmp[number_items - 1]->times_in_training_set());

  int ilow = 0;
  int ihigh = number_items - 1;
  int refinements_done = 0;

  while (ilow < ihigh && refinements_done < refinement_steps)
  {
    TTItem * low  = tmp[ilow];
    TTItem * high = tmp[ihigh];

#ifdef DEBUG_SWAPPING
    cerr << "Looking for swap low " << ilow << " now " << low->in_training_set() << " times " << low->times_in_training_set() << " and high " << ihigh << " now " << high->in_training_set() << " times " << high->times_in_training_set() << endl;
#endif

    if (high->times_in_training_set() == low->times_in_training_set())
      break;

    if (! high->in_training_set())    // we need to move these out of the training set
    {
      ihigh--;
      continue;
    }

    if (low->in_training_set())    // we need to increase traing set membership, low already in this set
    {
      ilow++;
      continue;
    }

    low->switch_training_set_membership();
    high->switch_training_set_membership();

#ifdef DEBUG_SWAPPING
    cerr << "Swapped " << low->id() << " now " << low->in_training_set() << " times " << low->times_in_training_set() << " and " << high->id() << " now " << high->in_training_set() << " times " << high->times_in_training_set() << endl;
#endif

    refinements_done++;

    ilow++;
    ihigh--;
  }

  total_refinements_done += refinements_done;

  return 1;
}

int
generate_random_split(const int* ndx,
                      int number_items,
                      int items_needed,
                      int split_number,
                      float average_times_chosen) {
  int items_selected = 0;

  for (int i = 0; i < number_items; ++i) {
    TTItem & t = ttitem[ndx[i]];

    if (t.in_training_set())
      continue;

    if (t.times_in_training_set() > 0)
      continue;

    t.set_in_training_set(1);

    items_selected++;
    if (items_selected >= items_needed)
      return items_needed;
  }

  for (int i = 0; i < number_items; ++i) {
#ifdef DEBUG_FUNCTION_OBJECTS
    cerr << "IN for loop, i = " << i << endl;
#endif

    TTItem & t = ttitem[ndx[i]];

    if (t.in_training_set()) {
      continue;
    }

    if (t.volunteer(number_splits - split_number, average_times_chosen)) {
      items_selected++;
      if (items_selected >= items_needed)
        return items_selected;
    }
  }

  return items_selected;
}

/*
  Since DELTA can be either positive of negative, we use function objects
  to control the stopping criterion
*/

template <typename O>
int
generate_random_split (int istart,
                       int delta,
                       int istop,
                       const O & o,
                       int items_needed,
                       int split_number,
                       int average_times_chosen)
{
  int items_selected = 0;

#ifdef DEBUG_FUNCTION_OBJECTS
  cerr << "Looking from " << istart << " to " << istop << " by " << delta << endl;
#endif

// First do a scan in which we just pick up things that have not been
// selected before

  for (int i = istart; o(i, istop); i += delta)
  {
    TTItem & t = ttitem[i];

    if (t.in_training_set())
      continue;

    if (t.times_in_training_set() > 0)
      continue;

    t.set_in_training_set(1);

    items_selected++;
    if (items_selected >= items_needed)
      return items_needed;
  }

  for (int i = istart; o(i, istop); i += delta)
  {
#ifdef DEBUG_FUNCTION_OBJECTS
    cerr << "IN for loop, i = " << i << endl;
#endif

    TTItem & t = ttitem[i];

    if (t.in_training_set())
      continue;

    if (t.volunteer(number_splits - split_number, average_times_chosen))
    {
      items_selected++;
      if (items_selected >= items_needed)
        return items_selected;
    }
  }

  return items_selected;
}

static int
generate_random_split (int items_in_training_set,
                       int split_number)
{
  set_training_set_membership(0);

  // how often should an individual have been chosen
  int average_times_chosen = static_cast<int>(fraction_in_training_set * split_number - 1);

  int items_selected = 0;

  static std::random_device rd;
  static std::default_random_engine rng;
  std::uniform_int_distribution<int> u0(0, number_items / 2);
  std::uniform_int_distribution<int> u1(1, number_items / 2);
  static std::bernoulli_distribution b(0.50);

  std::unique_ptr<int[]> ndx = std::make_unique<int[]>(number_items);
  std::iota(ndx.get(), ndx.get() + number_items, 0);
  std::shuffle(ndx.get(), ndx.get() + number_items, rng);

  while (1)
  {
    int items_needed = items_in_training_set - items_selected;
    items_selected += generate_random_split(ndx.get(), number_items, items_needed, split_number, average_times_chosen);

    if (items_selected >= items_in_training_set)
      return 1;
  }
}

template <typename C>
int
write_split (std::ostream & output,
             const C & c)
{
  int items_written = 0;

  if (write_header_to_split_files)
    output << "ID\n";

  for (int i = 0; i < number_items; i++)
  {
    const TTItem & t = ttitem[i];

    if (c(t.in_training_set()))
    {
      output << t.id() << '\n';
//    cerr << "Just wrote i = " << i << " '" << t.id() << "'\n";
      items_written++;
    }
  }

  if (verbose > 1)
    cerr << "Wrote " << items_written << " items\n";

  return output.good();
}

template <typename C>
int
write_split (IWString & fname,
             const C & c)
{
  std::ofstream output(fname.null_terminated_chars(), std::ios::out);

  if (! output.good())
  {
    cerr << "Cannot open output file '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    cerr << "Writing '" << fname << "'\n";

  return write_split(output, c);
}

/*
  I am sure there must be a way of doing this with standard library fns...
*/

class IWEqual
{
  private:
    int _n;

  public:
    IWEqual (int);

    void set (int s) { _n = s;}

    int operator() (int i) const { return _n == i;}
};

IWEqual::IWEqual (int i) : _n (i)
{
}

#ifdef __GNUG__
template int write_split (IWString &, const IWEqual &);
//template int write_split (std::ostream &, const IWEqual &);
#endif

static int
write_split (int s)
{
  IWString fname;

  IWEqual eq(1);

  if (write_training_set)
  {
    fname << output_stem << s << ".train";
    if (suffix.length())
      fname << suffix;

    if (! write_split(fname, eq))
      return 0;
  }

  if (write_test_set)
  {
    fname = output_stem;
    fname << s << ".test";
    if (suffix.length())
      fname << suffix;

    eq.set(0);

    return write_split(fname, eq);
  }

  return 1;
}

static void
initial_random_assignment (int items_in_training_set)
{
  int items_selected = 0;

  set_training_set_membership(0);

  std::random_device rd;
  std::default_random_engine rng;
  std::uniform_int_distribution<int> u(0, number_items - 1);
  while (items_selected < items_in_training_set)
  {
    const int i = u(rng);

    if (ttitem[i].in_training_set())
      continue;

    ttitem[i].set_in_training_set(1);
    items_selected++;
  }

  return;
}

static int
test_train_split_random()
{
  int items_in_training_set = static_cast<int>(static_cast<float>(number_items) * fraction_in_training_set + 0.01);

  if (0 == items_in_training_set)
  {
    cerr << "At " << fraction_in_training_set << " in the training set, only one molecule in training set\n";
    items_in_training_set = 1;
  }
  else if (items_in_training_set == number_items)
  {
    cerr << "At " << fraction_in_training_set << " in the training set, all " << number_items << " items in training set\n";
    items_in_training_set = number_items - 1;
  }

  if (verbose)
    cerr << items_in_training_set << " of " << number_items << " items in training set\n";

  initial_random_assignment(items_in_training_set);

  int splits_created = 0;

  if (! write_split(splits_created))
    return 0;

  splits_created++;

  while (splits_created < number_splits)
  {
    generate_random_split(items_in_training_set, splits_created);

    if (splits_created > refine_splits)
      adjust_split_toward_average();

    if (! write_split(splits_created))
      return 0;

    splits_created++;
  }

  return 1;
}

static int
allocate_arrays (int n)
{
  if (0 == n)
  {
    cerr << "No data, n = 0\n";
    return 0;
  }

  ttitem = new TTItem[n];

  if (NULL == ttitem)
  {
    cerr << "Cannot allocate " << n << " items\n";
    return 0;
  }

  number_items = n;

  return 1;
}

static int
read_identifier_tdt (const_IWSubstring buffer,    // not a reference
                     TTItem & t)
{
  assert (buffer.ends_with('>'));

  buffer.remove_leading_chars(identifier_tag.length());
  buffer.chop();

  if (buffer.nwords() > 1)
    buffer.truncate_at_first(' ');

  t.set_id(buffer);
  
  return 1;
}

static int
read_identifiers_tdt (iwstring_data_source & input)
{
  IWString tmp;
  tmp << "^" << identifier_tag;

  int n = input.grep(tmp);

  if (! allocate_arrays(n))
    return 0;

  const_IWSubstring buffer;

  n = 0;
  while (input.next_record(buffer))
  {
    if ("|" == buffer)
      n++;
    else if (! buffer.starts_with(identifier_tag))
      continue;
    else if (! read_identifier_tdt(buffer, ttitem[n]))
    {
      cerr << "Cannot process TDT record '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
read_identifier_dm (const const_IWSubstring & buffer,
                    TTItem & t)
{
  const_IWSubstring id;

  if (! buffer.word(1, id))
  {
    cerr << "Cannot extract identifier from DM record\n";
    return 0;
  }

  t.set_id(id);

  return 1;
}

static int
read_identifiers_dm (iwstring_data_source & input)
{
  int n = 0;
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with("ID"))
      n++;
    else if ("|" == buffer)
      break;
  }

  if (! allocate_arrays(n))
    return 0;

  if (! input.seekg(static_cast<off_t>(0)))
  {
    cerr << "Cannot seek back to beginning of dm file\n";
    return 0;
  }

  n = 0;
  while (input.next_record(buffer))
  {
    if ("|" == buffer)
      break;

    if (! buffer.starts_with("ID"))
      continue;

    if (! read_identifier_dm(buffer, ttitem[n]))
    {
      cerr << "Fatal error processing distance matrix record '" << buffer << "'\n";
      return 0;
    }

    n++;
  }

  return 1;
}

static int
read_identifier_smi (const const_IWSubstring & buffer,
                     TTItem & t)
{
  const_IWSubstring id;
  if (! buffer.word(1, id))
  {
    cerr << "Cannot extract smiles identifier\n";
    return 0;
  }

  t.set_id(id);

  return 1;
}

static int
read_identifiers_smi (iwstring_data_source & input)
{
  int n = input.records_remaining();

  if (! allocate_arrays(n))
    return 0;

  const_IWSubstring buffer;

  int i = 0;
  while (input.next_record(buffer))
  {
    if (! read_identifier_smi(buffer, ttitem[i]))
    {
      cerr << "Fatal error reading smiles on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    i++;
  }

  return 1;
}

static int
read_identifier_id (const_IWSubstring buffer,    // note, not a reference
                    TTItem & t)
{
  buffer.truncate_at_first(' ');

  t.set_id(buffer);

  return 1;
}

static int
read_identifiers_id (iwstring_data_source & input)
{
  int n = input.records_remaining();

  if (! allocate_arrays(n))
    return 0;

  const_IWSubstring buffer;

  int i = 0;

  while (input.next_record(buffer))
  {
    if (! read_identifier_id(buffer, ttitem[i]))
    {
      cerr << "Fatal error on line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }

    i++;
  }

  return 1;
}

static int
read_identifiers_dsc (iwstring_data_source & input)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "Cannot read descriptor file header record\n";
    return 0;
  }

  return read_identifiers_id(input);
}

static int
read_identifiers (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (INPUT_TYPE_ID == input_type)
    return read_identifiers_id(input);
  if (INPUT_TYPE_SMI == input_type)
    return read_identifiers_smi(input);
  if (INPUT_TYPE_TDT == input_type)
    return read_identifiers_tdt(input);
  if (INPUT_TYPE_DM == input_type)
    return read_identifiers_dm(input);
  if (INPUT_TYPE_DSC == input_type)
    return read_identifiers_dsc(input);

  const_IWSubstring myfname(fname);

  if (myfname.ends_with(".smi"))
    return read_identifiers_smi(input);

  if (myfname.ends_with(".gfp") || myfname.ends_with(".tdt"))
    return read_identifiers_tdt(input);

  if (myfname.ends_with(".dm"))
    return read_identifiers_dm(input);

  return read_identifiers_id(input);
}

static int
test_train_split_random (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vp:n:S:k:R:r:w:hF:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n'))
  {
    if (! cl.value('n', number_splits) || number_splits < 1)
    {
      cerr << "The number of splits to create option (-n) must be a whole positive number\n";
      usage (4);
    }

    if (verbose)
      cerr << "Will create " << number_splits << " random splits\n";
  }

  if (cl.option_present('p'))
  {
    int p;
    if (! cl.value ('p', p) || p < 1 || p > 99)
    {
      cerr << "The percent in training set option (-p) must be a valid percentage\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will put " << p << " percent of the items in the training set\n";

    fraction_in_training_set = static_cast<float> (p) / static_cast<float> (100.0);
  }

  if (cl.option_present('h'))
  {
    write_header_to_split_files = 1;

    if (verbose)
      cerr << "Will write a header record to each split\n";
  }

  if (cl.option_present('S'))
  {
    output_stem = cl.string_value ('S');

    if (verbose)
      cerr << "Will create file(s) will stem '" << output_stem << "'\n";
  }

  if (cl.option_present('F'))
  {
    suffix = cl.string_value ('F');

    if (verbose)
      cerr << "Will create file(s) will suffix '" << suffix << "'\n";

    if (! suffix.starts_with('.'))
      suffix.insert_at_beginning('.');
  }

  input_type = 0;

  if (cl.option_present('k'))
  {
    const_IWSubstring k = cl.string_value ('k');

    if ("id" == k)
      input_type = INPUT_TYPE_ID;
    else if ("smi" == k)
      input_type = INPUT_TYPE_SMI;
    else if ("tdt" == k)
      input_type = INPUT_TYPE_TDT;
    else if ("dm" == k)
      input_type = INPUT_TYPE_DM;
    else if ("dsc" == k)
      input_type = INPUT_TYPE_DSC;
    else
    {
      cerr << "Unrecognised input type specifier '" << k << "'\n";
      usage(8);
    }
  }

  if (cl.option_present('w'))
  {
    write_test_set = 0;
    write_training_set = 0;

    int i = 0;
    const_IWSubstring w;
    while (cl.value ('w', w, i++))
    {
      if ("test" == w)
      {
        write_test_set = 1;
        if (verbose)
          cerr << "Will write the test set identifiers\n";
      }
      else if ("train" == w)
      {
        write_training_set = 1;
        if (verbose)
          cerr << "Will write the training set identifiers\n";
      }
      else
      {
        cerr << "Unrecognised -w qualifier '" << w << "'\n";
        usage(5);
      }
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (1 != cl.number_elements())
  {
    cerr << "Sorry, only know how to process one file at a time\n";
    return 6;
  }

  if (! read_identifiers(cl[0]))
  {
    cerr << "Cannot read identifiers from '" << cl[0] << "'\n";
    return 5;
  }

  if (verbose)
    cerr << "Read " << number_items << " identifiers\n";

  if (cl.option_present('R'))
  {
    if (! cl.value('R', refine_splits) || refine_splits < 1)
    {
      cerr << "The refinement start (-R) option must be a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will refine all splits after the " << refine_splits << " split\n";
  }

  if (cl.option_present('r'))
  {
    if (0 == refine_splits)
      refine_splits = 2;

    if (! cl.value('r', refinement_steps) || refinement_steps < 1)
    {
      cerr << "The number of items to check for refinement option (-r) must be a whole positive number\n";
      usage(5);
    }

    if (refinement_steps > number_items / 2)
    {
      cerr << "Too many refinement steps, max " << (number_items / 2) << endl;
      refinement_steps = number_items / 2;
    }

    if (verbose)
      cerr << "Will refine up to " << refinement_steps << " items at each refinement\n";
  }
  else if (refine_splits)
  {
    refinement_steps = number_items / 10;
    if (0 == refinement_steps)
      refinement_steps = 2;
  }

  int rc = 1;
  for (int i = 0; i < number_items; i++)
  {
    if (0 == ttitem[i].id().length())
    {
      cerr << "Zero length identifier, impossible. i = " << i << endl;
      rc = 0;
    }
  }

  if (0 == rc)
    return 5;

  rc = test_train_split_random();

  if (0 == rc)
    return 3;

  if (verbose)
  {
    Accumulator_Int<int> acc;
    extending_resizable_array<int> times_in_training_set;

    for (int i = 0; i < number_items; i++)
    {
      const TTItem & ti = ttitem[i];

      int t = ti.times_in_training_set();

      acc.extra (t);
      times_in_training_set[t]++;

      if (verbose > 2)
        cerr << " i = " << i << " '" << ti.id() << " in training set " << ti.times_in_training_set() << endl;
    }

    cerr << "Items in training set between " << acc.minval() << " and " << acc.maxval() << " ave " << acc.average() << endl;
    for (int i = 0; i < times_in_training_set.number_elements(); i++)
    {
      if (times_in_training_set[i])
        cerr << times_in_training_set[i] << " items in training set " << i << " times\n";
    }

    if (total_refinements_done)
      cerr << total_refinements_done << " refinements done\n";
  }

  if (rc)
    return 0;
  else
    return 3;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = test_train_split_random (argc, argv);

  return rc;
}
