/*
  Group identical molecules together and study their activities
*/

#include <iostream>
#include <memory>
#include <functional>
#include <random>
using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "cmdline.h"
#include "iw_stl_hash_map.h"
#include "iw_stl_hash_set.h"
#include "accumulator.h"
#include "iwqsort.h"
#include "numeric_data_from_file.h"
#include "misc.h"

#include "istream_and_type.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "etrans.h"
#include "rmele.h"
#include "allowed_elements.h"

const char * prog_name = NULL;

static int verbose = 0;

static unsigned int molecules_read = 0;

static int remove_leading_zeros_from_identifiers = 0;

static Chemical_Standardisation chemical_standardisation;

static Allowed_Elements allowed_elements;

static int reduce_to_largest_fragment = 0;

static int ignore_chirality = 0;

static int reduce_to_graph_form = 0;

static int ignore_cis_trans_bonds = 0;

static Element_Transformations element_transformations;

static IWString_and_File_Descriptor stream_for_non_duplicates;

static float max_difference_for_merging = static_cast<float>(-1.0);

static int check_tolerances_by_ratio_method = 0;

static IWString_and_File_Descriptor stream_for_inconsistent_groups;

static int activity_column = 1;

static int take_first_of_multi_valued_data = 0;
static int take_average_of_multi_valued_data = 0;
static int take_max_of_multi_valued_data = 0;
static int take_min_of_multi_valued_data = 0;
static int vote_for_multi_valued_data = 0;   // not implemented
static int remove_multi_valued_data = 0;
static int make_replicate_molecules_of_multi_valued_data = 0;
static int take_random_value_from_range_of_multi_valued_data = 0;
static int discard_inconsistent_class_assignments = 0;

static int discarded_for_inconsistent_class_assignments = 0;

static int molecules_discarded_for_non_allowed_elements = 0;

static int molecules_discarded_for_covalent_non_organics = 0;

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = std::numeric_limits<int>::max();

static int molecules_discarded_for_too_few_atoms = 0;
static int molecules_discarded_for_too_many_atoms = 0;

static IWString activity_name;

/*
  As an optional behaviour we can only display the first line of a 
  common structure grouping
*/

static int only_display_common_group = 0;

/*
  When creating the -M file, we can write either the max value or the item
  closest to the average
*/

static int write_average_for_each_group = 0;
static int write_randomly_selected_item_for_each_group = 0;
static int write_median_item_for_each_group = 0;
static int write_random_value_from_range_of_each_group = 0;

static Elements_to_Remove elements_to_remove;

/*
  We may read in multi-valued data - same identifier, different values
*/

typedef IW_STL_Hash_Map<IWString, resizable_array<float> *> ID_to_Activity;

/*
  If we read in class data, we need to be able to write the same
  class labels when we are done
*/

static IW_STL_Hash_Map_float class_name_to_number;
static IWString negative_class, positive_class;

static std::random_device rd;

static IWString_and_File_Descriptor stream_for_multi_valued_data;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Groups identical molecules and compares activities\n";
  cerr << "  -X <fname>    experimental data\n";
  cerr << "  -e <col>      activity is column <col> in the -X file\n";
  cerr << "  -e <n>        activity is token <n> in the molecule name (no -X)\n";
  cerr << "  -a            reduce to graph form\n";
  cerr << "  -c            ignore chirality\n";
  cerr << "  -x            ignore cis-trans bonds\n";
  cerr << "  -m            only display common structure record of groups\n";
  cerr << "  -t <tol>      activity difference for items to be merged\n";
  cerr << "  -r            interpret the -t option as a ratio between activities\n";
  cerr << "  -O <fname>    write groups outside tolerance to <fname>\n";
  cerr << "  -N <fname>    write non duplicates to <fname> - cannot use with -t,-O\n";
  cerr << "  -M <stem>     write new smiles and activity files with max activity \n";
  cerr << "  -M max        write max activity to the -M file\n";
  cerr << "  -M ave        write average activity to the -M file\n";
  cerr << "  -M median     write median activity to the -M file\n";
  cerr << "  -M rand       write a random example to the -M file\n";
  cerr << "  -M range      write a random value from within the range to the -M file\n";
  cerr << "  -M rmdup      discard all duplicate structures\n";
  cerr << "  -M rminc      discard structures with inconsistent class assignments\n";
  cerr << "  -T ...        standard element transformation options\n";
  cerr << "  -V ...        what to do with multi-valued activities (-V help for info)\n";
  cerr << "  -K <ele>      remove elements of type <ele>\n";
  cerr << "  -W ...        specification of allowed elements\n";
  cerr << "  -U <fname>    tabular output format\n";
  cerr << "  -b <natoms>   lower atom count cutoff\n";
  cerr << "  -B <natoms>   lower atom count cutoff\n";
  cerr << "  -z            remove leading 0's from identifiers\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

class Smiles_ID_Activity
{
  private:
    IWString _smiles;
    IWString _id;
    float _activity;

  public:
    Smiles_ID_Activity(const IWString & smi, const IWString & id, float act) :
                        _smiles(smi), _id(id), _activity(act) {}


    const IWString & id() const { return _id;}

    float activity() const { return _activity;}
    void  set_activity (float s) { _activity = s;}

    int do_write (IWString_and_File_Descriptor & output) const;
};

int
Smiles_ID_Activity::do_write(IWString_and_File_Descriptor & output) const
{
  output << _smiles << ' ' << _id << ' ';

  if (0 == negative_class.length())   // not a classification problem
    output << _activity;
  else if (_activity < 0.0f)
    output << negative_class;
  else
    output << positive_class;

  output << '\n';

  return 1;
}

class Smiles_ID_Activity_Comparator
{
  private:
  public:
    int operator() (const Smiles_ID_Activity *, const Smiles_ID_Activity *) const;
};

int
Smiles_ID_Activity_Comparator::operator() (const Smiles_ID_Activity * sida1,
                                           const Smiles_ID_Activity * sida2) const
{
  if (sida1->activity() < sida2->activity())
    return -1;

  if (sida1->activity() > sida2->activity())
    return 1;

  return 0;
}

class Group_of_Molecules : public IWString
{
  public:
    typedef float activity_type_t;

  private:
    Accumulator<activity_type_t> _activity;

    resizable_array_p<Smiles_ID_Activity> _sida;

    std::mt19937_64 _rng;

//  private functions

    int _identify_item_with_max_activity () const;
    int _identify_item_with_activity_closest_to_average () const;
    int _identify_most_common_item (int & tied) const;

    int _write_max_activity (const Smiles_ID_Activity * sida, 
                             IWString_and_File_Descriptor & stream_for_smiles,
                             IWString_and_File_Descriptor & stream_for_activity) const;

  public:
    Group_of_Molecules (const IWString & s) : IWString(s), _rng(rd()) {};

    int n() const { return _activity.n();}

    int extra (const IWString & smi, const IWString & id, activity_type_t act);

    float max_difference() const;

    int is_within_tolerance (float t) const;
    int min_and_max_within_ratio (float t) const;

    int write_structure_group (IWString_and_File_Descriptor & output) const;

    int write_first_member (IWString_and_File_Descriptor & output) const;

    int write_merged_data (IWString_and_File_Descriptor & output) const;

    int write_max_activity (IWString_and_File_Descriptor &, IWString_and_File_Descriptor &) const;
    int write_random_item (IWString_and_File_Descriptor &, IWString_and_File_Descriptor &);
    int write_median_item (IWString_and_File_Descriptor & stream_for_smiles, IWString_and_File_Descriptor & stream_for_activity);
    int write_item_closest_to_average (IWString_and_File_Descriptor & stream_for_smiles,
                                       IWString_and_File_Descriptor & stream_for_activity) const;
    int write_most_common_item (IWString_and_File_Descriptor & stream_for_smiles,
                                IWString_and_File_Descriptor & stream_for_activity) const;
    int write_if_classes_consistent (IWString_and_File_Descriptor & stream_for_smiles,
                                     IWString_and_File_Descriptor & stream_for_activity) const;
    int write_first_member (IWString_and_File_Descriptor & stream_for_smiles,
                            IWString_and_File_Descriptor & stream_for_activity) const;

    int write_random_value_from_range_of_each_group (IWString_and_File_Descriptor & stream_for_smiles, IWString_and_File_Descriptor & stream_for_activity);

    int write_tabular_output (IWString_and_File_Descriptor &) const;
};

template class resizable_array_p<Smiles_ID_Activity>;
template class resizable_array_base<Smiles_ID_Activity *>;

int
Group_of_Molecules::extra (const IWString & smi,
                           const IWString & id,
                           activity_type_t act)
{
  _activity.extra(act);

  Smiles_ID_Activity * t = new Smiles_ID_Activity(smi, id, act);

  _sida.add(t);

  return 1;
}

float
Group_of_Molecules::max_difference() const
{
  if (_activity.n() < 2)
    return static_cast<float>(0.0);

  return _activity.maxval() - _activity.minval();
}

int
Group_of_Molecules::is_within_tolerance (float t) const
{
  return max_difference() <= t;
}

int
Group_of_Molecules::min_and_max_within_ratio (float t) const
{
  float mi = _activity.minval();
  float ma = _activity.maxval();

  if (fabs(mi - ma) < 1.0e-05)
    return 1;

  if (0.0F == ma)
    std::swap(mi, ma);

  float ratio = (ma - mi) / ma;

  return fabs(ratio) <= static_cast<double>(t);
}

int
Group_of_Molecules::write_structure_group (IWString_and_File_Descriptor & output) const
{
  output << (*this) << ' ';

  if (only_display_common_group)
  {
    for (int i = 0; i < _sida.number_elements(); i++)
    {
      if (i > 0)
        output << ':';

      output << _sida[i]->id();
    }

    output << " N = " << _activity.n() << " min " << _activity.minval() << " max " << _activity.maxval() << " ave " << static_cast<float>(_activity.average()) << " diff " << (static_cast<float>(_activity.maxval() - _activity.minval())) << '\n';

    return 1;
  }

  output << " N = " << _activity.n() << " min " << _activity.minval() << " max " << _activity.maxval() << " ave " << static_cast<float>(_activity.average()) << " diff " << (static_cast<float>(_activity.maxval() - _activity.minval())) << '\n';

  for (int i = 0; i < _sida.number_elements(); i++)
  {
    _sida[i]->do_write(output);
  }

  return 1;
}

int
Group_of_Molecules::write_first_member(IWString_and_File_Descriptor & output) const
{
  assert (_sida.number_elements() > 0);

  return _sida[0]->do_write(output);
}

class Group_of_Molecules_Comparator : public std::binary_function<const Group_of_Molecules *, 
                                                                  const Group_of_Molecules *,
                                                                  int>
{
  private:
  public:
    int operator() (const Group_of_Molecules *, const Group_of_Molecules *) const;
};

/*
  We have decided that the multiple values we have are within an
  acceptable tolerance. Write a consensus score
*/

int
Group_of_Molecules::write_merged_data (IWString_and_File_Descriptor & output) const
{
  assert (_sida.number_elements() > 1);

  output << (*this) << ' ';

  Accumulator<float> acc;

  for (int i = 0; i < _sida.number_elements(); i++)
  {
    const Smiles_ID_Activity * s = _sida[i];

    if (i > 0)
      output << ':';

    output << s->id();

    acc.extra(s->activity());
  }

  output << ' ' << static_cast<float>(acc.average()) << " N = " << _sida.number_elements() << '\n';

  return 1;
}

int
Group_of_Molecules::_write_max_activity (const Smiles_ID_Activity * sida, 
                                         IWString_and_File_Descriptor & stream_for_smiles,
                                         IWString_and_File_Descriptor & stream_for_activity) const
{
  sida->do_write(stream_for_smiles);

//cerr << "class_name_to_number " << class_name_to_number.size() << endl;
  if (0 == class_name_to_number.size())
    stream_for_activity << sida->id() << ' ' << sida->activity();
  else if (sida->activity() < 0.0F)
    stream_for_activity << sida->id() << ' ' << negative_class;
  else
    stream_for_activity << sida->id() << ' ' << positive_class;

  stream_for_activity << '\n';

  return 1;
}

int
Group_of_Molecules::_identify_item_with_max_activity() const
{
  int n = _sida.number_elements();

  if (1 == n)
    return 0;

  float max_activity = _sida[0]->activity();
  int rc = 0;

  for (int i = 1; i < n; i++)
  {
    if (_sida[i]->activity() > max_activity)
    {
      max_activity = _sida[i]->activity();
      rc = i;
    }
  }

  return rc;
}

int
Group_of_Molecules::_identify_most_common_item (int & tied) const   // vote
{
  tied = 0;

  int n = _sida.number_elements();

  if (1 == n)
    return 0;

  if (2 == n)
    return 0;

  resizable_array<float> v;
  resizable_array<int> c;
  resizable_array<int> f;     //first item with that value

  v.add(_sida[0]->activity());
  c.add(1);
  f.add(0);

  for (int i = 1; i < n; i++)
  {
    float a = _sida[i]->activity();

    int added_to_existing_value = 0;

    for (int j = 0; j < v.number_elements(); j++)
    {
      if (fabs(a - v[j]) < 1.0e-05)
      {
        c[j]++;
        added_to_existing_value = 1;
        break;
      }
    }

    if (! added_to_existing_value)
    {
      v.add(a);
      c.add(1);
      f.add(i);
    }
  }

  float most_common = v[0];
  int max_count = c[0];
  int item_with_most_common_value = 0;
  int items_with_highest_count = 1;

  for (int i = 1; i < v.number_elements(); i++)
  {
    if (c[i] > max_count)
    {
      most_common = v[i];
      max_count = c[i];
      item_with_most_common_value = f[i];
      items_with_highest_count = 1;
    }
    else if (c[i] == max_count)
      items_with_highest_count++;
  }

  if (1 == items_with_highest_count)
    return item_with_most_common_value;

  tied = 1;

  return item_with_most_common_value;
}

int
Group_of_Molecules::_identify_item_with_activity_closest_to_average() const
{
  int n = _sida.number_elements();

  if (1 == n)   // won't happen
  {
    return 0;
  }

  if (2 == n)    // if just two items, we take the one with highest activity
  {
    if (_sida[0]->activity() > _sida[1]->activity())
      return 0;
    else
      return 1;
  }

  Accumulator<float> acc;

  for (int i = 0; i < n; i++)
  {
    acc.extra(_sida[i]->activity());
  }

  double ave = acc.average();
  double min_d = fabs(_sida[0]->activity() - ave);

  int rc = 0;

  for (int i = 1; i < n; i++)
  {
    double d = fabs(_sida[i]->activity() - ave);

    if (d < min_d)
    {
      min_d = d;
      rc = i;
    }
  }

  return rc;
}


int
Group_of_Molecules::write_max_activity (IWString_and_File_Descriptor & stream_for_smiles,
                                        IWString_and_File_Descriptor & stream_for_activity) const
{
  int n = _sida.number_elements();

  if (1 == n)
    return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);

  int item_with_max_activity = _identify_item_with_max_activity();

  return _write_max_activity(_sida[item_with_max_activity], stream_for_smiles, stream_for_activity);
}

int
Group_of_Molecules::write_item_closest_to_average (IWString_and_File_Descriptor & stream_for_smiles,
                                                   IWString_and_File_Descriptor & stream_for_activity) const
{
  if (1 == _sida.number_elements())
    return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);

  int item_closest_to_average = _identify_item_with_activity_closest_to_average();

  float asave = _sida[item_closest_to_average]->activity();

  _sida[item_closest_to_average]->set_activity(_activity.average());

  int rc = _write_max_activity(_sida[item_closest_to_average], stream_for_smiles, stream_for_activity);

  _sida[item_closest_to_average]->set_activity(asave);

  return rc;
}

int
Group_of_Molecules::write_median_item (IWString_and_File_Descriptor & stream_for_smiles,
                                       IWString_and_File_Descriptor & stream_for_activity)
{
  int n = _sida.number_elements();

  if (n <= 2)
    return write_item_closest_to_average(stream_for_smiles, stream_for_activity);

  Smiles_ID_Activity_Comparator sidac;
  _sida.iwqsort(sidac);

  if (1 == n % 2)    // odd number of values
    return _write_max_activity(_sida[n / 2], stream_for_smiles, stream_for_activity);

// When an even number, take the average of the middle two

  activity_type_t a1 = _sida[n/2-1]->activity();
  activity_type_t a2 = _sida[n/2]->activity();

  activity_type_t ave = (a1 + a2) * 0.5f;

  activity_type_t asave = _sida[n/2]->activity();

  _sida[n/2]->set_activity(ave);

  int rc = _write_max_activity(_sida[n/2], stream_for_smiles, stream_for_activity);

  _sida[n/2]->set_activity(asave);

  return rc;
}

int
Group_of_Molecules::write_random_item (IWString_and_File_Descriptor & stream_for_smiles,
                                       IWString_and_File_Descriptor & stream_for_activity)
{
  int n = _sida.number_elements();

  if (1 == n)
    return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);

  std::uniform_int_distribution<int> u(0, n-1);

  int i = u(_rng);

  return _write_max_activity(_sida[i], stream_for_smiles, stream_for_activity);
}

int
Group_of_Molecules::write_most_common_item (IWString_and_File_Descriptor & stream_for_smiles,
                                            IWString_and_File_Descriptor & stream_for_activity) const
{
  int tied;

  int m = _identify_most_common_item(tied);

  return _write_max_activity(_sida[m], stream_for_smiles, stream_for_activity);
}

int
Group_of_Molecules::write_first_member (IWString_and_File_Descriptor & stream_for_smiles,
                                        IWString_and_File_Descriptor & stream_for_activity) const
{
  return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);
}

int
Group_of_Molecules::write_if_classes_consistent (IWString_and_File_Descriptor & stream_for_smiles,
                                                 IWString_and_File_Descriptor & stream_for_activity) const
{
  int n = _sida.number_elements();

  if (1 == n)
    return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);

  float class_assigned = _sida[0]->activity();

  for (int i = 1; i < n; i++)
  {
    float c = _sida[i]->activity();

    if (fabs(c - class_assigned) > 0.001)   // inconsistent, do NOT write
    {
      cerr << "Classes " << c << " and " << class_assigned << endl;
      discarded_for_inconsistent_class_assignments++;
      return 1;
    }
  }

  return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);
}

int
Group_of_Molecules::write_random_value_from_range_of_each_group (IWString_and_File_Descriptor & stream_for_smiles,
                                                IWString_and_File_Descriptor & stream_for_activity)
{
  if (1 == _sida.number_elements())
    return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);

  activity_type_t min_activity = _sida[0]->activity();
  activity_type_t max_activity = min_activity;

  for (int i = 1; i < _sida.number_elements(); ++i)
  {
    const auto a = _sida[i]->activity();

    if (a < min_activity)
      min_activity = a;
    else if (a > max_activity)
      max_activity = a;
  }

  if (min_activity == max_activity)
    return _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);

  std::uniform_real_distribution<double> u(min_activity, max_activity);

  min_activity = _sida[0]->activity();      //  save the value
  _sida[0]->set_activity(u(_rng));
  auto rc = _write_max_activity(_sida[0], stream_for_smiles, stream_for_activity);
  _sida[0]->set_activity(min_activity);
  return rc;
}

int
Group_of_Molecules::write_tabular_output (IWString_and_File_Descriptor & output) const
{
  int n = _sida.number_elements();

  if (1 == n)
    return 1;

  output << n << '\t' <<_activity.minval() << '\t' << _activity.maxval() << '\t' << (_activity.maxval() - _activity.minval()) << '\t' << static_cast<float>(_activity.average());
  for (int i = 0; i < n; i++)
  {
    output << '\t' << _sida[i]->id() << '\t' << _sida[i]->activity();
  }
  output << '\n';

  return 1;
}

int
Group_of_Molecules_Comparator::operator() (const Group_of_Molecules * g1, 
                                           const Group_of_Molecules * g2) const
{
  float d1 = g1->max_difference();
  float d2 = g2->max_difference();

  if (d1 < d2)
    return 1;

  if (d1 > d2)
    return -1;

  return 0;
}

static void
preprocess (Molecule & m)
{
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  return;
}

static int
fetch_activity_from_token_of_name(IWString & id,
                                  const ID_to_Activity & id_to_activity,
                                  resizable_array<float> & activity)
{
  const_IWSubstring token;
  if (! id.word(activity_column, token))
  {
    cerr << "Not enough tokens in name '" << id << "'\n";
    return 0;
  }

  float a;
  if (! token.numeric_value(a))
  {
    cerr << "Invalid activity token '" << token << "'\n";
    return 0;
  }

  activity.add(a);

  id.truncate_at_first(' ');

  return 1;
}

static int
fetch_activity (const ID_to_Activity & id_to_activity,
                IWString & id,
                resizable_array<float> & activity)
{
  if (0 == id_to_activity.size())
    return fetch_activity_from_token_of_name(id, id_to_activity, activity);

  id.truncate_at_first(' ');

  ID_to_Activity::const_iterator f = id_to_activity.find(id);

  if (f != id_to_activity.end())
  {
    activity = *(f->second);
    return 1;
  }

  if (! remove_leading_zeros_from_identifiers)
    return 0;
  else if (! id.starts_with('0'))
    return 0;

  IWString tmp(id);
  tmp.remove_leading_chars('0');

  f = id_to_activity.find(tmp);

  if (f == id_to_activity.end())
    return 0;

  activity = *(f->second);

  return 1;
}

static int
perform_transformations(Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (ignore_chirality)
    m.remove_all_chiral_centres();

  if (ignore_cis_trans_bonds)
    m.revert_all_directional_bonds_to_non_directional();

  if (reduce_to_graph_form)
    m.change_to_graph_form();

  if (element_transformations.active())
    element_transformations.process(m);

  if (elements_to_remove.active())
    elements_to_remove.process(m);

  return 1;
}

static int
discarded_for_atom_count(const Molecule & m)
{
  const auto matoms = m.natoms();
  if (matoms < lower_atom_count_cutoff)
  {
    if (verbose)
      cerr << m.name() << " rejected for too few atoms " << m.natoms() << endl;
    molecules_discarded_for_too_few_atoms++;
    return 1;
  }

  if (matoms > upper_atom_count_cutoff)
  {
    if (verbose)
      cerr << m.name() << " rejected for too many atoms " << m.natoms() << endl;
    molecules_discarded_for_too_many_atoms++;
    return 1;
  }

  return 0;     // not discarded
}

static int
contains_non_allowed_elements(const Molecule & m,
                              const Allowed_Elements & allowed_elements)
{
  if (allowed_elements.contains_non_allowed_atoms(m))
  {
    cerr << m.name() << " rejected for non allowed elements\n";
    molecules_discarded_for_non_allowed_elements++;
    return 1;
  }

  return 0;
}

static int
discarded_for_covalently_bonded_non_organics(const Molecule & m)
{
  if (m.organic_only())
    return 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    if (a->element()->organic())
      continue;

    if (a->ncon() > 0)
      return 1;
  }

  return 0;
}

/*
  We always return 1
*/

static int
group_molecules(Molecule & m,
                const ID_to_Activity & id_to_activity,
                IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                IW_STL_Hash_Set & chirality_removed)
{
  if (contains_non_allowed_elements(m, allowed_elements))   // check before any frgment reduction
    return 1;

  if (discarded_for_covalently_bonded_non_organics(m))
    return 1;

  IWString mname = m.name();

  resizable_array<float> activity;

  if (! fetch_activity(id_to_activity, mname, activity))
  {
    cerr << "No activity data for '" << mname << "'\n";
    return 0;
  }

  mname.truncate_at_first(' ');

  if (ignore_chirality && m.chiral_centres())
    chirality_removed.insert(mname);

  int need_to_make_copy = 0;

  if (reduce_to_graph_form)
    need_to_make_copy = 1;
  else if (ignore_chirality && m.chiral_centres())
    need_to_make_copy = 1;
  else if (ignore_cis_trans_bonds && m.cis_trans_bonds_present())
    need_to_make_copy = 1;
  else if (element_transformations.active())
    need_to_make_copy = 1;
  else if (reduce_to_largest_fragment && m.number_fragments() > 1)
    need_to_make_copy = 1;

  IWString usmi;

  if (need_to_make_copy)
  {
    int hcount;
    if (reduce_to_graph_form)
      hcount = m.implicit_hydrogens();
    else
      hcount = 0;

    Molecule mcopy(m);
    perform_transformations(mcopy);

    if (discarded_for_atom_count(mcopy))
      return 1;

    usmi = mcopy.unique_smiles();
    if (hcount > 0)
      usmi << ":H" << hcount;
  }
  else
  {
    if (discarded_for_atom_count(m))
      return 1;

    usmi = m.unique_smiles();
  }

  const auto f = structure_group.find(usmi);

  Group_of_Molecules * g;

  if (f == structure_group.end())
  {
    g = new Group_of_Molecules(usmi);
    structure_group[usmi] = g;
  }
  else
    g = (*f).second;

  if (! need_to_make_copy)
    m.invalidate_smiles();

  for (int i = 0; i < activity.number_elements(); ++i)
  {
    g->extra(m.smiles(), mname, activity[i]);
  }

  return 1;
}

static int
group_molecules(data_source_and_type<Molecule> & input,
                const ID_to_Activity & id_to_activity,
                IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                IW_STL_Hash_Set & chirality_removed)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

//  cerr << "Read " << m->name() << endl;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! group_molecules(*m, id_to_activity, structure_group, chirality_removed))
      return 0;
  }

  return 1;
}

static int
group_molecules(const char * fname, int input_type, 
                const ID_to_Activity & id_to_activity,
                IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                IW_STL_Hash_Set & chirality_removed)
{
  assert (NULL != fname);

  if (0 == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return group_molecules(input, id_to_activity, structure_group, chirality_removed);
}

static int
check_for_within_tolerance(const Group_of_Molecules & g,
                           const float max_difference_for_merging)
{
  if (check_tolerances_by_ratio_method)
    return g.min_and_max_within_ratio(max_difference_for_merging);

  return g.is_within_tolerance(max_difference_for_merging);
}

static int
do_merge_output (const IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                 IWString_and_File_Descriptor & output)
{
  int molecules_written = 0;
  int groups_within_tolerance = 0;
  int groups_outside_tolerance = 0;

  Accumulator<float> acc;

  for (auto i : structure_group)
  {
    const Group_of_Molecules * g = i.second;

    if (1 == g->n())
    {
      g->write_first_member(output);
      molecules_written++;
      continue;
    }

    acc.extra(g->max_difference());

    if (g->is_within_tolerance(max_difference_for_merging))
    {
      g->write_merged_data(output);
      groups_within_tolerance++;
      molecules_written++;
    }
    else
    {
      groups_outside_tolerance++;
      if (stream_for_inconsistent_groups.is_open())
        g->write_structure_group(stream_for_inconsistent_groups);
    }
  }

  if (verbose)
  {
    cerr << "Wrote " <<  molecules_written << " from " << stream_for_inconsistent_groups.size() << " structure groups\n";
    cerr << groups_within_tolerance << " groups within tolerance, " << groups_outside_tolerance << " outside tolerance\n";
    if (acc.n())
      cerr << "Within group activity differences between " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average_if_available_minval_if_not()) << '\n';
  }

  return 1;
}

/*
  The -M file stuff
*/

static int
write_reconciled_data (const IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                       IWString_and_File_Descriptor & stream_for_smiles,
                       IWString_and_File_Descriptor & stream_for_activity)
{
  int groups_outside_tolerance = 0;

  for (auto i : structure_group)
  {
    Group_of_Molecules * g = i.second;   // not const because the median method does a sort

    if (1 == g->n())
      g->write_first_member(stream_for_smiles, stream_for_activity);
    else if (remove_multi_valued_data)
      continue;
    else if (max_difference_for_merging >= 0.0 && ! check_for_within_tolerance(*g, max_difference_for_merging))
    {
      if (stream_for_inconsistent_groups.is_open())
        g->write_structure_group(stream_for_inconsistent_groups);

      groups_outside_tolerance++;

      continue;
    }
    else if (write_average_for_each_group)
      g->write_item_closest_to_average(stream_for_smiles, stream_for_activity);
    else if (write_median_item_for_each_group)
      g->write_median_item(stream_for_smiles, stream_for_activity);
    else if (discard_inconsistent_class_assignments)
      g->write_if_classes_consistent(stream_for_smiles, stream_for_activity);
    else if (take_first_of_multi_valued_data)
      g->write_first_member(stream_for_smiles, stream_for_activity);
    else if (vote_for_multi_valued_data)
      g->write_most_common_item(stream_for_smiles, stream_for_activity);
    else if (write_randomly_selected_item_for_each_group)
      g->write_random_item(stream_for_smiles, stream_for_activity);
    else if (write_random_value_from_range_of_each_group)
      g->write_random_value_from_range_of_each_group(stream_for_smiles, stream_for_activity);
    else
      g->write_max_activity (stream_for_smiles, stream_for_activity);

    stream_for_smiles.write_if_buffer_holds_more_than(32768);
    stream_for_activity.write_if_buffer_holds_more_than(32768);
  } 

  if (verbose && max_difference_for_merging > 0.0)
    cerr << groups_outside_tolerance << " groups outside tolerance\n";

  return 1;
}

static int
count_items_per_structure(const IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                          extending_resizable_array<int> & items_per_structure)
{
  for (auto i : structure_group)
  {
    const Group_of_Molecules * g = i.second;

    int n = g->n();

    items_per_structure[n]++;
  }

  if (verbose)
  {
    for (int i = 0; i < items_per_structure.number_elements(); i++)
    {
      if (items_per_structure[i] > 0)
        cerr << items_per_structure[i] << " structures had " << i << " examples\n";
    }
  }

  return 1;
}

static int
write_reconciled_data (const IW_STL_Hash_Map<IWString, Group_of_Molecules *> & structure_group,
                       const IWString & stem)
{
  IWString_and_File_Descriptor stream_for_smiles;
  IWString_and_File_Descriptor stream_for_activity;

  IWString fname;
  fname << stem << ".smi";

  if (! stream_for_smiles.open(fname.c_str()))
  {
    cerr << "Cannot open maximum activity smiles file '" << fname << "'\n";
    return 0;
  }

  fname = stem;
  fname << ".activity";

  if (! stream_for_activity.open(fname.c_str()))
  {
    cerr << "Cannot open maximum activity file '" << fname << "'\n";
    return 0;
  }

  stream_for_activity << "ID " << activity_name << '\n';

  return write_reconciled_data(structure_group, stream_for_smiles, stream_for_activity);
}

static float
average(const resizable_array<float> & v)
{
  const auto n = v.number_elements();

  float sum = 0.0f;
  for (auto x : v)
  {
    sum += x;
  }

  return sum / static_cast<float>(n);
}

static void
update_global_accumulators(const resizable_array<float> & v,
                           Accumulator<float> & acc_diff,
                           int & number_zero_differences,
                           const float arange)
{
//const auto ave = average(v);

  const auto n = v.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    const auto vi = v[i];

    for (unsigned int j = i + i; j < n; ++j)
    {
      const auto d = fabsf(v[j] - vi);

      if (d < (0.01f * arange))
        number_zero_differences++;
      else
        acc_diff.extra(d);
    }
  }

  return;
}

static int
write_multi_valued_information(const IWString & id, 
                               resizable_array<float> & acc, 
                               IWString_and_File_Descriptor & output)
{
  const char sep = ' ';

  output << id << sep << acc.number_elements() << sep << acc.min_val() << sep << acc.max_val() << sep << static_cast<float>(average(acc));

  acc.iwqsort_lambda([](const float f1, const float f2)
  {
    if (f1 < f2)
      return -1;

    if (f1 > f2)
      return 1;

    return 0;
  });

  for (auto x : acc)
  {
    output << sep << x;
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
handle_multi_valued_activities(ID_to_Activity & id_to_activity)
{
  float min_activity = std::numeric_limits<float>::max();
  float max_activity = - min_activity;

  for (auto i : id_to_activity)
  {
    const auto * acc = i.second;

    const auto n = acc->size();

    if (1 == n)
    {
      const auto v0 = acc->item(0);
      if (v0 < min_activity)
        min_activity = v0;
      if (v0 > max_activity)
        max_activity = v0;
    }
    else
    {
      for (unsigned int j = 0; j < n; ++j)
      {
        const auto vj = acc->item(j);
        if (vj < min_activity)
          min_activity = vj;
        if (vj > max_activity)
          max_activity = vj;
      }
    }
  }

  std::mt19937_64 rng(rd());

  Accumulator<float> acc_diff;
  int zero_diff = 0;

  for (ID_to_Activity::iterator i = id_to_activity.begin(); i != id_to_activity.end(); ++i)
  {
    auto * acc = (*i).second;

    if (1 == acc->number_elements())
      continue;

    update_global_accumulators(*acc, acc_diff, zero_diff, max_activity - min_activity);

    if (stream_for_multi_valued_data.is_open())
      write_multi_valued_information(i->first, *acc, stream_for_multi_valued_data);

    if (take_average_of_multi_valued_data)
    {
      double a = average(*acc);
      acc->resize_keep_storage(0);
      acc->add(a);
    }
    else if (take_max_of_multi_valued_data)
    { 
      double a = acc->max_val();
      acc->resize_keep_storage(0);
      acc->add(a);
    }
    else if (take_min_of_multi_valued_data)
    {
      double a = acc->min_val();
      acc->resize_keep_storage(0);
      acc->add(a);
    }
    else if (take_random_value_from_range_of_multi_valued_data)
    {
      std::uniform_real_distribution<double> u(acc->min_val(), acc->max_val());
      acc->resize_keep_storage(0);
      acc->add(u(rng));
    }
    else if (vote_for_multi_valued_data)
    {
      if (acc->min_val() != acc->max_val())
        acc->resize_keep_storage(0);
    }
  }

  cerr << acc_diff.n() << " measurements (same ID, different value), diffs btw " << acc_diff.minval() << " and " << acc_diff.maxval() << " ave " << static_cast<float>(acc_diff.average());
  if (acc_diff.n() > 2)
    cerr << " std " << static_cast<float>(sqrt(acc_diff.variance()));
  cerr << endl;

  if (zero_diff > 0)
  {
    acc_diff.extra(0.0f, zero_diff);
    cerr << zero_diff << " measurements (same ID, same value), ave diff including these " << static_cast<float>(acc_diff.average()) << "\n";
  }

  return 1;
}

static int
read_activity_record(const const_IWSubstring & buffer,
                     ID_to_Activity & id_to_activity)
{
  static bool first_call = true;

  IWString id;

  int i = 0;
  buffer.nextword(id, i);

  if (remove_leading_zeros_from_identifiers)
    id.remove_leading_chars('0');

  const_IWSubstring token;

  if (1 == activity_column)
  {
    if (! buffer.nextword(token, i))
    {
      cerr << "NO experimental data available\n";
      return 0;
    }
//  cerr << "Buffer '" << buffer << "', extracted '" << token << "'\n";
  }
  else
  {
    if (! buffer.word(activity_column, token))
    {
      cerr << "Cannot fetch token for experimental data, activity_column " << activity_column << endl;
      return 0;
    }
  }

  if (0 == token.length())   // how could that happen?
  {
    cerr << "Cannot extract activity token, activity_column " << activity_column << endl;
    return 0;
  }

  float a;

  if (token.numeric_value(a))    // great
    first_call = false;
  else if (first_call)   // header record in activity file
  {
    activity_name = token;
    first_call = false;
    return 1;
  }
  else if ('.' == token)
  {
    cerr << "Ignoring missing value for '" << id << "'\n";
    return 1;
  }
  else
  {
    unsigned int s = class_name_to_number.size();

    const auto f = class_name_to_number.find(token);

    if (f != class_name_to_number.end())
      a = (*f).second;
    else if (0 == s)
    {
      class_name_to_number[token] = -1.0f;
      negative_class = token;
      a = -1.0f;
    }
    else if (1 == s)
    {
      class_name_to_number[token] =  1.0f;
      positive_class = token;
      a = 1.0f;
    }
    else
    {
      cerr << "Non numeric experimental value '" << token << "', maybe too many classes\n";
      return 0;
    }
  }

  ID_to_Activity::iterator f = id_to_activity.find(id);

  if (f == id_to_activity.end())
  {
    resizable_array<float> * tmp = new resizable_array<float>;
    tmp->add(a);

    id_to_activity[id] = tmp;

    return 1;
  }

  if (verbose > 1)
    cerr << "Duplicate data for '" << id << "'\n";

  if (take_first_of_multi_valued_data)
    return 1;

  (*f).second->add(a);

  return 1;
}

static int 
read_activity_data(iwstring_data_source & input,
                   ID_to_Activity & id_to_activity)
{
  input.set_translate_tabs(1);
  input.set_strip_trailing_blanks(1);

  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;

    if (! read_activity_record(buffer, id_to_activity))
    {
      cerr << "Bad activity data record '" << buffer << "'\n";
      return 0;
    }
  }

  return static_cast<int>(id_to_activity.size());
}

static int 
read_activity_data(const char * fname,
                   ID_to_Activity & id_to_activity)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open activity file '" << fname << "'\n";
    return 0;
  }

  return read_activity_data(input, id_to_activity);
}

static int
do_write_tabular_output(const Group_of_Molecules * const * gm,
                        const int n,
                        IWString_and_File_Descriptor & output)
{
  output << "N\tMinval\tMaxval\tRange\tAverage\tID\tActivity\n";

  for (int i = 0; i < n; i++)
  {
    gm[i]->write_tabular_output(output);

    output.write_if_buffer_holds_more_than(32768);
  }

  output.flush();

  return 1;
}

static void
display_dash_V_options(std::ostream & os)
{
  os << " -V first         take first of multi-valued activity data\n";
  os << " -V ave           take average of multi-valued activity data\n";
  os << " -V max           take max value of multi-valued activity data\n";
  os << " -V WRITE=<fname> write duplicate info to <fname>\n";
  os << " -V append        append all values associated with each id\n";
//os << " -V rm            remove all multi-valued activity data\n";

  exit(0);
}

static int
activity_consistency (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:lX:caT:szmN:O:t:rxM:V:e:K:W:U:b:B:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('T'))
  {
    if (! element_transformations.construct_from_command_line(cl, verbose, 'T'))
    {
      cerr << "Cannot initialise element transformations (-T)\n";
      usage(6);
    }

    if (verbose)
      cerr << "Initialised " << element_transformations.number_elements() << " element transformations\n";
  }

  if (cl.option_present('K'))
  {
    if (! elements_to_remove.construct_from_command_line(cl, verbose, 'K'))
    {
      cerr << "Cannot initialise elements to remove (-K) option\n";
      usage(3);
    }
  }

  if (cl.option_present('W'))
  {
    if (! allowed_elements.build_from_command_line(cl, 'W', verbose))
    {
      cerr << "Cannot initialise allowed elements (-W) option\n";
      usage(3);
    }
  }

  if (cl.option_present('c'))
  {
    ignore_chirality = 1;

    if (verbose)
      cerr << "Will ignore chirality\n";
  }

  if (cl.option_present('x'))
  {
    ignore_cis_trans_bonds = 1;

    if (verbose)
      cerr << "Will ignore cis trans bonds\n";
  }

  if (cl.option_present('a'))
  {
    reduce_to_graph_form = 1;

    if (verbose)
      cerr << "Will compare molecules in their graph form\n";
  }

  if (cl.option_present('V'))
  {
    int specifications_given = 0;

    IWString v;
    for (int i = 0; cl.value('V', v, i); ++i)
    {
      IWString lcv = v;     // lowercase version of V
      lcv.to_lowercase();

      if ("max" == lcv)
      {
        take_max_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will take the max of multi-valued activity data\n";
        specifications_given++;
      }
      else if ("min" == lcv)
      {
        take_min_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will take the min of multi-valued activity data\n";
        specifications_given++;
      }
      else if ("ave" == lcv)
      {
        take_average_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will take the ave of multi-valued activity data\n";
        specifications_given++;
      }
      else if ("first" == lcv)
      {
        take_first_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will take the first of multi-valued activity data\n";
        specifications_given++;
      }
      else if ("vote" == lcv)
      {
        vote_for_multi_valued_data = 1;
        if (verbose)
          cerr << "Assuming classification model, will vote for multi-valued data\n";
        specifications_given++;
      }
      else if ("extra" == lcv || "append" == lcv)
      {
        make_replicate_molecules_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Duplicate data associated with same ID gets counted as extra sample\n";
        specifications_given++;
      }
      else if ("range" == lcv)
      {
        take_random_value_from_range_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will make a random selection from within the range if any ID has multiple values\n";
        specifications_given++;
      }
      else if (v.starts_with("WRITE="))
      {
        v.remove_leading_chars(6);
        if (! stream_for_multi_valued_data.open(v.null_terminated_chars()))
        {
          cerr << "Cannot open stream for multi valued data '" << v << "'\n";
          return 1;
        }

        if (verbose)
          cerr << "Multi valued data written to '" << v << "'\n";

        stream_for_multi_valued_data << "ID N Min Max Ave ...\n";
      }
      else if ("help" == lcv)
      {
        display_dash_V_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -V qualifier '" << v << "'\n";
        display_dash_V_options(cerr);
      }
    }

    if (specifications_given > 1)
    {
      cerr << "Can only specify one course of action when multiple identifiers present\n";
    }
  }

  if (cl.option_present('e'))
  {
    if (! cl.value('e', activity_column) || activity_column < 2)
    {
      cerr << "INvalid activity column specification(-x)\n";
      usage(5);
    }

    if (verbose)
      cerr << "Activity data in column " << activity_column << "\n";

    activity_column--;
  }

  if (cl.option_present('X'))
    ;
  else if (activity_column > 0)
    ;
  else
  {
    cerr << "Must specify activity data via the -X option or\n";
    usage(4);
  }

  if (cl.option_present('z'))
  {
    remove_leading_zeros_from_identifiers = 1;

    if (verbose)
      cerr << "Will strip leading 0's from identifiers\n";
  }

  ID_to_Activity id_to_activity;

  if (cl.option_present('X'))
  {
    const char * x = cl.option_value('X');

    if (! read_activity_data (x, id_to_activity))
    {
      cerr << "Cannot read activity data (-X), file '" << x << "'\n";
      return 8;
    }

    if (verbose)
      cerr << "Read " << id_to_activity.size() << " activity dataitems\n";

    if (1 == class_name_to_number.size())
    {
      cerr << "Got just one class label, something is wrong with the activity file, look for non numeric values..\n";
      return 2;
    }

    extending_resizable_array<int> values_per_id;

    int multi_valued_data_present = 0;   // same ID with multiple data values

    for (auto i : id_to_activity)
    {
      const auto * acc = i.second;

      int n = acc->number_elements();

      values_per_id[n]++;

      if (n > 1)
        multi_valued_data_present++;
    }

    if (multi_valued_data_present)
    {
      cerr << "Multi valued data present for " << multi_valued_data_present << " molecules\n";
      for (int i = 1; i < values_per_id.number_elements(); i++)
      {
        if (values_per_id[i] > 0)
          cerr << values_per_id[i] << " items had " << i << " values\n";
      }

      if (take_max_of_multi_valued_data)
        ;
      else if (take_min_of_multi_valued_data)
        ;
      else if (take_average_of_multi_valued_data)
        ;
      else if (take_first_of_multi_valued_data)   // will not happen here
        ;
      else if (vote_for_multi_valued_data)
        ;
      else if (make_replicate_molecules_of_multi_valued_data)
        ;
      else
      {
        cerr << "Must specify some handling of multi-valued data via the -V option\n";
        display_dash_V_options(cerr);
      }

      handle_multi_valued_activities(id_to_activity);
    }
  }

//if (class_name_to_number.size() > 0)
//  setup_positive_and_negative_class_labels();

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage (6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('N'))
  {
    const char * n = cl.option_value('N');

    if (! stream_for_non_duplicates.open(n))
    {
      cerr << "Cannot open stream for non duplicates '" << n << "'\n";
      return 6;
    }

    if (verbose)
      cerr << "Non-duplicate molecules written to '" << n << "'\n";
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', max_difference_for_merging) || max_difference_for_merging < 0.0)
    {
      cerr << "The max difference for merging (-t) must be a valid non negative value\n";
      usage(6);
    }

    if (verbose)
      cerr << "Differences less than " << max_difference_for_merging << " will be merged\n";

    if (cl.option_present('r'))
    {
      check_tolerances_by_ratio_method = 1;
      if (verbose)
        cerr << "Tolerances for grouping done by ratio\n";
    }
  }

  if (cl.option_present('O'))
  {
    if (! cl.option_present('t'))
    {
      cerr << "When using a inconsistent stream (-O) must also specify max tolerance (-T)\n";
      usage(4);
    }

    IWString o = cl.string_value('O');

    if (! o.ends_with(".smi"))
      o << ".smi";

    if (! stream_for_inconsistent_groups.open(o.null_terminated_chars()))
    {
      cerr << "Cannot open stream for inconsistent data '" << o << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Differences beyond " << max_difference_for_merging << " written to '" << o << "'\n";
  }

  IWString_and_File_Descriptor stream_for_tabular_output;

  if (cl.option_present('U'))
  {
    const char * u = cl.option_value('U');

    if (! stream_for_tabular_output.open(u))
    {
      cerr << "Cannot open stream for tabular output '" << u << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "Tabular output written to '" << u << "'\n";
  }

  if (cl.option_present('b'))
  {
    if (! cl.value('b', lower_atom_count_cutoff) || lower_atom_count_cutoff < 1)
    {
      cerr << "THe lower atom count cutoff (-b) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will discard molecules with fewer than " << lower_atom_count_cutoff << " atoms\n";
  }

  if (cl.option_present('B'))
  {
    if (! cl.value('B', upper_atom_count_cutoff) || upper_atom_count_cutoff < lower_atom_count_cutoff)
    {
      cerr << "THe upper atom count cutoff (-B) must be a whole +ve number, greater than " << lower_atom_count_cutoff << endl;
      usage(2);
    }

    if (verbose)
      cerr << "Will discard molecules with more than " << upper_atom_count_cutoff << " atoms\n";
  }

  set_copy_name_in_molecule_copy_constructor(1);

  IW_STL_Hash_Map<IWString, Group_of_Molecules *> structure_group;

  IW_STL_Hash_Set chirality_removed;

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! group_molecules(cl[i], input_type, id_to_activity, structure_group, chirality_removed))
    {
      cerr << "Cannot read molecules from '" << cl[i] << "'\n";
      return 4;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, " << structure_group.size() << " unique_structures\n";
  }

  if (molecules_read == structure_group.size())
    cerr << "No duplicate structures encountered\n";

  if (cl.option_present('m') && cl.option_present('t'))
  {
    cerr << "Sorry, the -m and -t options cannot be used together\n";
    usage(5);
  }

  if (cl.option_present('s') && cl.option_present('t'))
  {
    cerr << "Sorry, the -s and -t options cannot be used together\n";
    usage(5);
  }

  extending_resizable_array<int> items_per_structure;

  count_items_per_structure(structure_group, items_per_structure);

  if (max_difference_for_merging > static_cast<float>(0.0) && ! cl.option_present('M'))
  {
    IWString_and_File_Descriptor output(1);
    return do_merge_output(structure_group, output);
  }

  if (cl.option_present('m'))
  {
    only_display_common_group = 1;
  }

  if (cl.option_present('M'))
  {
    IWString fname;

    int i = 0;
    const_IWSubstring m;

    int write_specifications = 0;

    while (cl.value('M', m, i++))
    {
      if ("ave" == m || "AVE" == m)
      {
        write_average_for_each_group = 1;
        if (verbose)
          cerr << "Will write the average activity for each group\n";
        write_specifications++;
      }
      else if ("median" == m)
      {
        write_median_item_for_each_group = 1;
        if (verbose)
          cerr << "Will write the median for odd numbered groups of molecules\n";
        write_specifications++;
      }
      else if ("rand" == m)
      {
        write_randomly_selected_item_for_each_group = 1;
        if (verbose)
          cerr << "Will write a randomly selected molecule from each group\n";
        write_specifications++;
      }
      else if ("range" == m)
      {
        write_random_value_from_range_of_each_group = 1;
        if (verbose)
          cerr << "Will choose a random value from within the range of each group\n";
        write_specifications++;
      }
      else if ("max" == m || "MAX" == m)
      {
        take_max_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will take the largest activity of a multi-structure group\n";
        write_specifications++;
      }
      else if ("min" == m || "MIN" == m)
      {
        take_min_of_multi_valued_data = 1;
        if (verbose)
          cerr << "Will take the smallest activity value of a multi-structure group\n";
        write_specifications++;
      }
      else if ("rmdup" == m || "rm" == m)
      {
        remove_multi_valued_data = 1;

        if (verbose)
          cerr << "Will remove multi valued data\n";
      }
      else if ("rminc" == m || "RMINC" == m)
      {
        discard_inconsistent_class_assignments = 1;
        if (verbose)
          cerr << "In a classification scenario, will discard structures with conflicting class assignments\n";
      }
      else if (0 == fname.length())
        fname = m;
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        usage(3);
      }
    }

    if (0 == fname.length())
    {
      cerr << "NO file name specified for -M option\n";
      usage(1);
    }

    if (write_specifications > 1)
    {
      cerr << "Sorry, must specify just one writing specification for multi-valued data\n";
      usage(5);
    }

    if (verbose)
      cerr << "Writing to '" << fname << "'\n";

    if (! write_reconciled_data(structure_group, fname))
      return 3;

    if (class_name_to_number.size() > 0)
    {
      take_max_of_multi_valued_data = 0;
      take_min_of_multi_valued_data = 0;
      vote_for_multi_valued_data = 1;
    }
  }

  int single_page_output = cl.option_present('s');

  int max_examples_per_structure = 0;

  int items_to_display = 0;

  for (int i = 0; i < items_per_structure.number_elements(); i++)
  {
    if (0 == items_per_structure[i])  // no instances of I examples of a structure
      continue;

    if (items_per_structure[i] > max_examples_per_structure)
      max_examples_per_structure = items_per_structure[i];

    items_to_display += i * items_per_structure[i];
  }

  for (auto i: structure_group)
  {
    const Group_of_Molecules * g = i.second;

    int n = g->n();

    items_per_structure[n]++;

    if (n > max_examples_per_structure)
      max_examples_per_structure = n;

    if (n > 1)
      items_to_display++;
  }

  if (single_page_output)
  {
    max_examples_per_structure++;    // don't forget parent structure

    if (0 == max_examples_per_structure % 2)
      ;
    else
      max_examples_per_structure += 1;
  }

  if (0 == items_to_display)
    cerr << "No multi-valued structures encountered???\n"; 

  Group_of_Molecules ** gm = new Group_of_Molecules*[items_to_display]; std::unique_ptr<Group_of_Molecules*[]> free_gm(gm);

  items_to_display = 0;

  Accumulator<float> range_acc;

  for(auto i: structure_group)
  {
    Group_of_Molecules * g = i.second;

    int n = g->n();

    if (1 == n)
    {
      if (stream_for_non_duplicates.active())
        g->write_first_member(stream_for_non_duplicates);
      continue;
    }

    if (verbose)
      range_acc.extra(g->max_difference());

    gm[items_to_display] = g;
    items_to_display++;
  }

  if (0 == verbose)
    ;
  else if (range_acc.n() > 0)
    cerr << items_to_display << " groups, ranges between " << range_acc.minval() << " to " << range_acc.maxval() << " ave " << static_cast<float>(range_acc.average_if_available_minval_if_not()) << endl;

  Group_of_Molecules_Comparator gmc;
  if (items_to_display > 1)
    iwqsort (gm, items_to_display, gmc);

  if (stream_for_tabular_output.is_open())
    do_write_tabular_output(gm, items_to_display, stream_for_tabular_output);

  IWString_and_File_Descriptor output(1);

  for (int i = 0; i <items_to_display; i++)
  {
    const Group_of_Molecules * g = gm[i];

    g->write_structure_group(output);

    if (single_page_output)
    {
      for (int j = g->n() + 1; j < max_examples_per_structure; j++)
      {
        output << "*\n";
      }
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  if (verbose)
  {
    if (discard_inconsistent_class_assignments)
      cerr << discarded_for_inconsistent_class_assignments << " molecules discarded for inconsistent class assignments\n";
    if (molecules_discarded_for_non_allowed_elements)
      cerr << molecules_discarded_for_non_allowed_elements << " molecules_discarded_for_non_allowed_elements\n";
    if (molecules_discarded_for_covalent_non_organics)
      cerr << molecules_discarded_for_covalent_non_organics << " discarded for covalently bonded non organics\n";
  }

#ifdef CARE_ABOUT_FREEING_MEMORY
  for (auto i : id_to_activity)
  {
    delete i.second;
  }

  for (auto i : structure_group)
  {
    delete i.second;
  }
#endif

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = activity_consistency (argc, argv);

  return rc;
}
