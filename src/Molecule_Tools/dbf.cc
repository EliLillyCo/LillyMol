/*
  We often want descriptors that are the distances between features.

  We can compute either topological or through space distances

  Lots of hard-coded stuff here for donors/acceptors and charges
*/

#include <stdlib.h>
#include <memory>
#include <iostream>

using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/histogram/iwhistogram.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/minmaxspc.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#ifdef DBF_MOLVOL
#include "surface_area_molvol.h"
#endif

static Elements_to_Remove elements_to_remove;

static Chemical_Standardisation chemical_standardisation;

static Element_Transformations element_transformations;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static Charge_Assigner charge_assigner;

static int verbose = 0;

static int reduce_to_largest_fragment = 0;

static int molecules_read = 0;

static int do_topological_distances = 0;
static int do_spatial_distances = 0;

static Accumulator<coord_t> distance_stats;

#ifdef DBF_MOLVOL
static Surface_Area_Molvol surface_area_molvol;
#endif

/*
  If someone specifies a histogram on the command line, we store the specification here, and 
  then use this to initialise all the histograms
*/

static IWString histogram_specification;

/*
  The user may want to specify their own features
*/

static resizable_array_p<Substructure_Query> queries;

/*
  Once we introduced the positive_negative_acceptor_donor concept, user specified
  queries will be placed after those, so each query now needs a feature number
*/

static resizable_array<int> query_to_feature_number;

/*
  But the feature name array is across all features, including positive_negative_acceptor_donor
*/

static resizable_array_p<IWString> feature_name;

static int nq = 0;

/*
  There are some features that are too hard to implement with queries. If
  they are active, their value will be set to their putative query number
*/

static int hard_coded_query_centre_of_aromatic_5_membered_ring = -1;
static int hard_coded_query_centre_of_aromatic_6_membered_ring = -1;

static int number_hard_coded_queries = 0;

static extending_resizable_array<int> matches_per_molecule;
static extending_resizable_array<int> molecules_matching_query;
static extending_resizable_array<int> matches_per_query;

/*
  What to do with multiple hits to a query
*/

static int take_first_of_multiple_hits = 0;

static int ignore_queries_hitting_multiple_times = 0;

/*
  Sometimes we want to ignore interactions that may be closer
  than a given distance. I've implemented a filter for 2d, do 
  I need one for 3d?
*/

static int ignore_atoms_closer_than_2d = 0;

static int consider_bonded_features = 0;

static coord_t ignore_atoms_closer_than_3d = 0.0;

static Molecule_Output_Object stream_for_labelled_atoms;

/*
  Sometimes we may want the NKEEP shortest distances as descriptors
*/

static int nkeep = 0;

static IWString missing_value('.');

/*
  If a query doesn't hit any atoms, we need a string for everything missing
*/

static IWString all_values_missing;

static IWDigits iwdigits;

static Fraction_as_String fraction_as_string;

static IWString_and_File_Descriptor stream_for_all_distances;

static IWString tag;
static IWString smiles_tag("$SMI<");

/*
  We can speed things up and avoid problems with isotopes by hard coding the common pharmacaphore things
*/

static int positive_negative_acceptor_donor = 0;

static float bit_count_scaling_factor = 1.0f;

#ifdef MAX_DBF_NOT_USED
static float max_dbf [] = {2,
2,
7,
5,
10,
10,
10,
9,
10,
9,
10,
16,
12.21,
11,
15,
12,
7,
8,
7,
7,
15,
10.25,
10,
14,
11,
12,
19,
13,
10,
18,
11.33,
12,
16,
12,
9.507,
10.12,
9.765,
8.58,
9.399,
8.88,
9.108,
15.05,
11.13,
10.26,
13.98,
11.23,
7.068,
7.284,
7.191,
6.963,
14.58,
9.603,
9.837,
12.81,
10.85,
10.58,
17.39,
11.75,
8.937,
17.51,
10.68,
10.96,
15.16,
11.26,
1.2,
1.224,
1.2,
1.065,
1.155,
1.104,
1.14,
1.245,
1.161,
1.23,
1.239,
1.23,
1.002,
1.041,
1.017,
1.11,
1.218,
1.137,
1.17,
1.227,
1.182,
1.236,
1.263,
1.236,
1.2,
1.329,
1.218,
1.23,
1.245,
1.23};
#endif


template <typename T>
class Distances : public Accumulator_Base<T, T>, private resizable_array<T>
{
#ifdef IW_TWO_PHASE_TEMPLATES
  protected:
    using resizable_array<T>::_number_elements;
    using resizable_array<T>::_things;
    using resizable_array<T>::_elements_allocated;
#endif
  private:

//  There are two histograms. _global_histogram is accumulated across the whole dataset.
//  _current_molecule_histogram is for the current molecule

    IWHistogram * _global_histogram;

    IWHistogram * _current_molecule_histogram;

  public:
    Distances ();
    ~Distances ();

    int activate_current_molecule_histogram ();

    void extra (T);

    void reset ();

    void append_results (IWString &) const;
    void append_current_molecule_histogram_results (IWString &) const;
    void append_results_integer (IWString &) const;
    void append_results_float   (IWString &) const;
    void write_all_results (IWString_and_File_Descriptor & buffer);

    void set_bits (int b, float bit_count_scaling_factor, Sparse_Fingerprint_Creator &) const;

    const IWHistogram * global_histogram () const { return _global_histogram;}
};

template <typename T>
Distances<T>::Distances ()
{
  if (histogram_specification.length() > 0)
  {
    _global_histogram = new IWHistogram;
    (void) _global_histogram->initialise(histogram_specification);
  }
  else
  {
    _global_histogram = nullptr;
    _current_molecule_histogram = nullptr;
  }

  return;
}

template <typename T>
int
Distances<T>::activate_current_molecule_histogram ()
{
  assert (nullptr == _current_molecule_histogram);
  assert (0 != histogram_specification.length());

  _current_molecule_histogram = new IWHistogram;

  return _current_molecule_histogram->initialise(histogram_specification);
}

template <typename T>
Distances<T>::~Distances ()
{
  if (nullptr != _global_histogram)
    delete [] _global_histogram;

  if (nullptr != _current_molecule_histogram)
    delete [] _current_molecule_histogram;

  return;
}

/*
  Note that _global_histogram does not get reset during a normal reset. That histogram is designed to 
  accumulate data from all molecules
*/

template <typename T>
void
Distances<T>::reset ()
{
  Accumulator_Base<T, T>::reset();

  if (nullptr != _current_molecule_histogram)
  {
    if (nullptr != _global_histogram)
      _global_histogram->add(*_current_molecule_histogram);

    _current_molecule_histogram->reset();
  }

  this->resize_keep_storage(0);

  return;
}

template <typename T>
void
Distances<T>::extra (T v)
{
  Accumulator_Base<T, T>::extra(v);

  if (nullptr != _current_molecule_histogram)
    _current_molecule_histogram->extra(static_cast<double>(v));

  if (0 == nkeep)
    return;

  if (0 == _elements_allocated)
  {
    this->resize(nkeep);
    this->add(v);

    return;
  }

  if (0 == _number_elements || (_number_elements < nkeep && v >= this->last_item()))
  {
    this->add(v);

    return;
  }

  for (int i = 0; i < _number_elements; i++)
  {
    if (v < _things[i])
    {
      this->insert_before(i, v);
      break;
    }
  }

  if (_number_elements > nkeep)
    this->resize(nkeep);

  return;
}

template <typename T>
void
Distances<T>::append_results (IWString & buffer) const
{
  int number_samples = Accumulator_Base<T, T>::n();

//cerr << "distribution has " << _number_elements << " items, n = " << number_samples << " between " << minval() << " and " << maxval() << endl;

  if (0 == number_samples)
  {
    buffer << all_values_missing;
    return;
  }

  buffer << ' ' << this->minval() << ' ' << this->maxval();

  buffer << ' ';

  if (1 == number_samples)
    buffer << this->minval();
  else
    buffer << static_cast<float>(this->average());

  if (0 == nkeep)
    return;

  assert (_number_elements <= nkeep);

  for (int i = 1; i < _number_elements; i++)    // the first sample stored is the minimum value - we already wrote that
  {
    buffer << ' ' << _things[i];
  }

  int nextra = nkeep - _number_elements + 1;

  if (nextra > 0)
  {
    for (int i = 0; i < nextra; i++)
    {
      buffer << ' ' << missing_value << i;
    }
  }

  return;
}

template <typename T>
class Float_Comparator
{
  private:
  public:
    int operator() (T, T) const;
};

template <typename T>
int
Float_Comparator<T>::operator() (T f1, T f2) const
{
  if (f1 < f2)
    return -1;

  if (f1 > f2)
    return 1;

  return 0;
}

template <typename T>
void
Distances<T>::write_all_results (IWString_and_File_Descriptor & buffer)
{
  Float_Comparator<T> c;

  this->iwqsort(c);

  for (int i = 0; i < _number_elements; i++)
  {
    buffer << ' ';

    buffer << _things[i];
  }

  return;
}

template <typename T>
void
Distances<T>::append_results_integer (IWString & buffer) const
{
  int number_samples = Accumulator_Base<T, T>::n();

#ifdef DEBUG_APPEND_RESULTS_INTEGER
  cerr << "distribution has " << _number_elements << " items, n = " << number_samples << " between " << this->minval() << " and " << this->maxval() << endl;
#endif

  if (0 == number_samples)
  {
    buffer << all_values_missing;
    return;
  }

  iwdigits.append_number(buffer, static_cast<int>(this->minval()));
  iwdigits.append_number(buffer, static_cast<int>(this->maxval()));

  if (1 == number_samples)
    iwdigits.append_number(buffer, static_cast<int>(this->minval()));
  else
    buffer << ' ' << static_cast<float>(this->average());

  if (0 == nkeep)
    return;

  assert (_number_elements <= nkeep);

  for (int i = 1; i < _number_elements; i++)    // the first sample stored is the minimum value - we already wrote that
  {
    iwdigits.append_number(buffer, static_cast<int>(_things[i]));
  }

  int nextra = nkeep - _number_elements + 1;

  if (nextra > 0)
  {
    for (int i = 0; i < nextra; i++)
    {
      buffer << ' ' << missing_value;
    }
  }

  return;
}

template <typename T>
void
Distances<T>::append_results_float (IWString & buffer) const
{
  int number_samples = Accumulator_Base<T, T>::n();

#ifdef DEBUG_APPEND_RESULTS_FLOAT
  cerr << "distribution has " << _number_elements << " items, n = " << number_samples << " between " << this->minval() << " and " << this->maxval();
  if (this->n() > 0)
    cerr << " ave " << this->average();
  cerr << endl;
#endif

  if (0 == number_samples)
  {
    buffer << all_values_missing;
    return;
  }

  fraction_as_string.append_number(buffer, this->minval());
  fraction_as_string.append_number(buffer, this->maxval());

  if (1 == number_samples)
    fraction_as_string.append_number(buffer, this->minval());
  else
   fraction_as_string.append_number(buffer, this->average());

  if (0 == nkeep)
    return;

  assert (_number_elements <= nkeep);

  for (int i = 1; i < _number_elements; i++)    // the first sample stored is the minimum value - we already wrote that
  {
    fraction_as_string.append_number(buffer, _things[i]);
  }

  int nextra = nkeep - _number_elements + 1;

  if (nextra > 0)
  {
    for (int i = 0; i < nextra; i++)
    {
      buffer << ' ' << missing_value;
    }
  }

  return;
}

template <typename T>
void 
Distances<T>::set_bits (int b, float bit_count_scaling_factor, Sparse_Fingerprint_Creator & sfc) const
{
  const int number_samples = Accumulator_Base<T, T>::n();

  if (0 == number_samples)
    return;

  sfc.hit_bit(b,   static_cast<int>(this->minval()  * bit_count_scaling_factor) + 1);
  sfc.hit_bit(b+1, static_cast<int>(this->maxval()  * bit_count_scaling_factor) + 1);
  sfc.hit_bit(b+2, static_cast<int>(this->average() * bit_count_scaling_factor) + 1);

  return;
}

template <typename T>
void
Distances<T>::append_current_molecule_histogram_results (IWString & buffer) const
{
  assert (nullptr != _current_molecule_histogram);

  const unsigned int * c = _current_molecule_histogram->raw_counts();
  int nbuckets = _current_molecule_histogram->nbuckets();

  for (int i = 0; i < nbuckets; i++)
  {
    int ci = c[i];

    iwdigits.append_number(buffer, ci);
  }

  return;
}

template class Distances<int>;
template class Distances<coord_t>;

static Distances<int> * bond_distances = nullptr;
static Distances<coord_t> * spatial_distances = nullptr;

/*
  We can optionally do the ratio between the spatial and bond distances
*/

static int compute_spatial_topological_ratio = 0;

static Distances<coord_t> * spatial_bond_ratio = nullptr;

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

  cerr << "  -p 2d          compute topological distances\n";
  cerr << "  -p 3d          compute spatial     distances\n";
  cerr << "  -q <query>     identify features as queries\n";
  cerr << "  -s <smarts>    identify features as smarts\n";
  cerr << "  -e             implement the default positive, negative, acceptor, donor queries (no smarts needed)\n";
  cerr << "  -k <nkeep>     report the <nkeep> shortest distances for each type\n";
  cerr << "  -z first       take the first match when multiple matches to a query\n";
  cerr << "  -z ignore      ignore any query that matches multiple times\n";
  cerr << "  -z oknomatch   ignore molecules where the query doesn't match\n";
  cerr << "  -b <number>    ignore atom pairs <number> bonds apart or closer\n";
  cerr << "  -d             allow bonded features\n";
  cerr << "  -T <hist>      specify histogram for profiling -T min,max,nbuckets\n";
  cerr << "  -C <file>      detect all spatial distances and write to <file>\n";
  cerr << "  -r             compute spatial/topological distance ratios\n";
  cerr << "  -W ...         various hard-coded queries, enter '-W help' for info\n";
  cerr << "  -L <fname>     write labelled molecules to <fname> - last query match wins!\n";
  cerr << "  -i <type>      input type\n";
  (void) display_standard_aromaticity_options(cerr);
  cerr << "  -t ...         element transformation options, enter '-t help'\n";
  (void) display_standard_chemical_standardisation_options(cerr, 'g');
  (void) display_standard_charge_assigner_options(cerr, 'N');
  cerr << "  -H <...>       donor/acceptor options, enter '-H help' for details\n";
  cerr << "  -u <string>    missing value specification\n";
  cerr << "  -n             suppress normal output\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -E ...         standard element options, enter '-E help' for info\n";
#ifdef DBF_MOLVOL
  (void) display_standard_vdw_radius_types(cerr, 'V');
#endif
  cerr << "  -v             verbose output\n";

  exit(rc);
}

/*
  When recording matched atoms, we need to differentiate between matches
  that match just a single atom and those that match multiple atoms. For
  those that match multiple atoms, we create a dummy atom at the geometric
  centre of the atoms that matched.
*/

class Coords_and_Stuff : public Coordinates
{
  private:

//  If this is the result of a single atom match, then _zatom will be set

    atom_number_t _zatom;

//  the query number from which we are derived

    const int _query_number;

  public:
    Coords_and_Stuff (int);

    int query_number () const { return _query_number;}

    int single_atom_match () const { return INVALID_ATOM_NUMBER != _zatom; }

    atom_number_t matched_atom () const { return _zatom;}

    void set_matched_atom (atom_number_t a) { _zatom = a;}
};

Coords_and_Stuff::Coords_and_Stuff (int qnum) : _query_number(qnum)
{
  _zatom = INVALID_ATOM_NUMBER;

  return;
}

template class resizable_array_p<Coords_and_Stuff>;
template class resizable_array_base<Coords_and_Stuff *>;

static int
write_atom_hit_info (Molecule & m,
                     const int * atom_label,
                     const resizable_array_p<Coords_and_Stuff> & c,
                     Molecule_Output_Object & output)
{
  IWString save_name(m.name());

  IWString tmp(m.name());

  for (int i = 0; i < c.number_elements(); ++i)
  {
    tmp << ' ' << c[i]->matched_atom() << ':' << c[i]->query_number();
  }

  m.set_name(tmp);

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    m.set_isotope(i, i);
  }

  int rc = output.write(m);

  m.set_name(save_name);

  return rc;
}

static int
add_matched_atom (const Molecule & m,
                  const atom_number_t zatom,
                  const int query_number,
                  resizable_array_p<Coords_and_Stuff> & c)
{
  assert (query_number >= 0);

  matches_per_query[query_number]++;

  Coords_and_Stuff * tmp = new Coords_and_Stuff(query_number);

  tmp->set_matched_atom(zatom);

  const Atom * a = m.atomi(zatom);

  tmp->setxyz(*a);

  c.add(tmp);

  return 1;
}

/*
  Identify the atoms that match each query
*/

static int
identify_matched_atoms (Molecule & m,
                        int query_number,
                        const Set_of_Atoms & embedding,
                        resizable_array_p<Coords_and_Stuff> & c)
{
  assert (query_number >= 0);

  int atoms_in_embedding = embedding.number_elements();

// If the smarts contained atoms that are to be excluded from the embedding, there will be
// INVALID_ATOM_NUMBER items in the embedding. Note we just check the 2nd atom, which is
// not really correct

  if (atoms_in_embedding > 1 && INVALID_ATOM_NUMBER == embedding[1])
    atoms_in_embedding = 1;

  if (1 == atoms_in_embedding)
    return add_matched_atom(m, embedding[0], query_number, c);

  Coords_and_Stuff * tmp = new Coords_and_Stuff(query_number);

// Multiple atoms matched, we must form the geometric centre

  matches_per_query[query_number]++;

  tmp->setxyz(static_cast<coord_t>(0.0), static_cast<coord_t>(0.0), static_cast<coord_t>(0.0));

  for (int i = 0; i < atoms_in_embedding; i++)
  {
    atom_number_t j = embedding[i];

    if (INVALID_ATOM_NUMBER == j)
      continue;

    const Atom * ai = m.atomi(embedding[i]);

    tmp->operator += (*ai);
  }

  tmp->operator /= (static_cast<coord_t>(atoms_in_embedding));

  c.add(tmp);

  return 1;
}

static int
do_hard_coded_query_centre_of_aromatic_ring (Molecule & m,
                                             int ring_size,
                                             int query_number,
                                             resizable_array_p<Coords_and_Stuff> & c)
{
  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->number_elements() != ring_size)
      continue;

    if (! ri->is_aromatic())
      continue;

    identify_matched_atoms(m, query_number, *ri, c);
    rc++;
  }

  if (rc)
    molecules_matching_query[query_number]++;

  return rc;
}

static int
do_hard_coded_features (Molecule & m,
                        resizable_array_p<Coords_and_Stuff> & c,
                        int & matches_this_molecule)
{
  m.compute_aromaticity_if_needed();

  if (hard_coded_query_centre_of_aromatic_6_membered_ring >= 0)
    matches_this_molecule += do_hard_coded_query_centre_of_aromatic_ring(m, 6, hard_coded_query_centre_of_aromatic_6_membered_ring, c);

  if (hard_coded_query_centre_of_aromatic_5_membered_ring >= 0)
    matches_this_molecule += do_hard_coded_query_centre_of_aromatic_ring(m, 5, hard_coded_query_centre_of_aromatic_5_membered_ring, c);

  return 1;
}

static int
do_positive_negative_acceptor_donor (const Molecule & m,
                                     int * atom_label,
                                     resizable_array_p<Coords_and_Stuff> & c)
{
  const int matoms = m.natoms();

  int pnad_matches[4];
  std::fill_n(pnad_matches, 4, 0);

  for (int i = 0; i < matoms; ++i)
  {
    const formal_charge_t fc = m.formal_charge(i);
    if (0 == fc)
      ;
    else if (1 == fc)
    {
      if (nullptr != atom_label)
        atom_label[i] = 1;
      pnad_matches[0] += add_matched_atom(m, i, 0, c);
    }
    else if (-1 == fc)
    {
      if (nullptr != atom_label)
        atom_label[i] = 2;
      pnad_matches[1] += add_matched_atom(m, i, 1, c);
    }

    const isotope_t iso = m.isotope(i);
    if (0 == iso)
      continue;

    if (1 == iso || 2 == iso)
    {
      if (nullptr != atom_label)
        atom_label[i] = 3;
      pnad_matches[2] += add_matched_atom(m, i, 2, c);
    }
    if (2 == iso || 3 == iso)
    {
      if (nullptr != atom_label)
        atom_label[i] = 4;
      pnad_matches[3] += add_matched_atom(m, i, 3, c);
    }
  }

  int rc = 0;

  for (int i = 0; i < 4; ++i)
  {
    if (0 == pnad_matches[i])
      continue;

    rc += pnad_matches[i];

    molecules_matching_query[i]++;
  }

  return rc;
}

static int
identify_matched_atoms (Molecule & m,
                        const isotope_t * isotope,
                        resizable_array_p<Coords_and_Stuff> & c)
{
  int matches_this_molecule = 0;     // the number of queries finding a match

  int * atom_label = nullptr;

  if (stream_for_labelled_atoms.is_open())
    atom_label = new_int(m.natoms());

  std::unique_ptr<int[]> free_atom_label(atom_label);

  if (positive_negative_acceptor_donor)
  {
    matches_this_molecule += do_positive_negative_acceptor_donor(m, atom_label, c);
    if (nullptr != isotope)
    {
      m.set_isotopes(isotope);
    }
  }

  Molecule_to_Match target(&m);

  for (int i = 0; i < nq; i++)
  {
    Substructure_Results sresults;

//  cerr << " Query " << i << " (" << queries[i]->comment() << ") is query_to_feature_number " << query_to_feature_number[i] << " name " << feature_name[query_to_feature_number[i]] << endl;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (verbose > 2)
      cerr << nhits << " hits to query " << i << endl;

    if (0 == nhits)
      continue;

    molecules_matching_query[query_to_feature_number[i]]++;

    matches_this_molecule++;

    if (1 == nhits)
      ;
    else if (take_first_of_multiple_hits)
      nhits = 1;
    else if (ignore_queries_hitting_multiple_times)
      continue;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      if (nullptr != atom_label)
        e->set_vector(atom_label, query_to_feature_number[i] + 1);

      identify_matched_atoms(m, query_to_feature_number[i], *e, c);
    }
  }

  do_hard_coded_features(m, c, matches_this_molecule);

  matches_per_molecule[matches_this_molecule]++;

  if (stream_for_labelled_atoms.is_open())
    return write_atom_hit_info(m, atom_label, c, stream_for_labelled_atoms);

  return 1;
}

template <typename T>
void
reset_distances (Distances<T> * d,
                 const int nfeatures)
{
  for (int i = 0; i < nfeatures; i++)
  {
    for (int j = i; j < nfeatures; j++)
    {
      d[i * nfeatures + j].reset();
    }
  }

  return;
}

template void reset_distances(Distances<int> *, int);
template void reset_distances(Distances<coord_t> *, int);

/*
  We are processing a pair of matches, and need to know whether or not
  they violate any topological rules. Note that we can only apply topology
  rules to matches with one atom
*/

static int
ok_topology (Molecule & m,
             const Coords_and_Stuff & ci,
             const Coords_and_Stuff & cj,
             int & b)
{
  b = 0;

  if (! ci.single_atom_match())
    return 1;

  if (! cj.single_atom_match())
    return 1;

  atom_number_t anumi = ci.matched_atom();
  atom_number_t anumj = cj.matched_atom();

  if (anumi == anumj)      // not needed, will be caught by the zero spatial distance check. Just paranoid...
    return 0;

  b = m.bonds_between(anumi, anumj);
  assert (b > 0);

  if (b > 1)
    ;
  else if (consider_bonded_features)
    return 1;
  else
    return 0;

  if (ignore_atoms_closer_than_2d && b <= ignore_atoms_closer_than_2d)
    return 0;

  return 1;
}

static int
do_spatial_and_topological_distances_computation (Molecule & m,
                                                  const int nfeatures,
                                                  const resizable_array_p<Coords_and_Stuff> & coords)
{
  reset_distances(bond_distances, nfeatures);
  reset_distances(spatial_distances, nfeatures);
  reset_distances(spatial_bond_ratio, nfeatures);

  int rc = 0;

  int nc = coords.number_elements();

  for (int i = 0; i < nc; i++)
  {
    Coords_and_Stuff * ci = coords[i];

    int iqn = ci->query_number();

    for (int j = i + 1; j < nc; j++)
    {
      Coords_and_Stuff * cj = coords[j];

      int jqn = cj->query_number();

      coord_t d = ci->distance(*cj);

      if (d <= ignore_atoms_closer_than_3d)
        continue;
      
//    If these are single atom matches, we can do checks on connectivity, and we can
//    compute topological parameters

      int b;
      if (! ok_topology(m, *ci, *cj, b))
        continue;

      int qni = iqn;
      int qnj = jqn;
      if (qni > qnj)
        std::swap(qni, qnj);

      Distances<int>   & topolog_dist = bond_distances[qni * nfeatures + qnj];
      Distances<float> & spatial_dist = spatial_distances[qni * nfeatures + qnj];
      Distances<float> & sptopratio   = spatial_bond_ratio[qni * nfeatures + qnj];

      if (b > 0)     // bond distance determined, these must be single atom matches
      {
        topolog_dist.extra(b);

        if (compute_spatial_topological_ratio)
          sptopratio.extra(static_cast<float>(d) / static_cast<float>(b));
      }

      spatial_dist.extra(d);

      rc++;
    }
  }

  return rc;
}

/*
  The business with lstart is to avoid duplicating things in the I set of atoms
*/

static int
do_topological_distances_computation (Molecule & m,
                                      const int nfeatures,
                                      const resizable_array_p<Coords_and_Stuff> & coords)
{
  reset_distances(bond_distances, nfeatures);

  const int nc = coords.number_elements();
//cerr << "in " << m.name() << " " << nc << " features\n";

  if (0 == nc)
    return 0;

  for (int i = 0; i < nc; i++)
  {
    const Coords_and_Stuff * ci = coords[i];

    if (! ci->single_atom_match())     // cannot do topological matches on multi-atom-matches
      continue;

    const atom_number_t ai = ci->matched_atom();

    const int iqn = ci->query_number();

    for (int j = i + 1; j < nc; j++)
    {
      const Coords_and_Stuff * cj = coords[j];

      if (! cj->single_atom_match())
        continue;

      atom_number_t aj = cj->matched_atom();

      if (ai == aj)
        continue;

      int d = m.bonds_between(ai, aj);

      if (1 == d)      // are bonded
        continue;

      if (d <= ignore_atoms_closer_than_2d)
        continue;

      int jqn = cj->query_number();

//    cerr << "Btw " << iqn << " and " << jqn << " dist " << d << endl;

      if (iqn <= jqn)
        bond_distances[iqn * nfeatures + jqn].extra(d);
      else
        bond_distances[jqn * nfeatures + iqn].extra(d);
    }
  }

  return 1;
}

static int
do_spatial_distances_computation (Molecule & m,
                                  const int nfeatures,
                                  const resizable_array_p<Coords_and_Stuff> & coords)
{
  reset_distances(spatial_distances, nfeatures);

  const int nc = coords.number_elements();

  int rc = 0;

  for (int i = 0; i < nc; i++)
  {
    Coords_and_Stuff * ci = coords[i];

    int iqn = ci->query_number();

    for (int j = i + 1; j < nc; j++)
    {
      Coords_and_Stuff * cj = coords[j];

      int jqn = cj->query_number();

      coord_t d = ci->distance(*cj);

      if (d < ignore_atoms_closer_than_3d)
        continue;

      int b;     // not used here
      if (! ok_topology(m, *ci, *cj, b))
        continue;

      if (iqn <= jqn)
        spatial_distances[iqn * nfeatures + jqn].extra(d);
      else
        spatial_distances[jqn * nfeatures + iqn].extra(d);

      rc++;
    }
  }

  return rc;
}

static int
gather_results (const resizable_array_p<Coords_and_Stuff> & coords,
                Sparse_Fingerprint_Creator & sfc,
                const int nfeatures,
                IWString & buffer)
{
  int * tmp = new_int(nfeatures); std::unique_ptr<int[]> free_tmp(tmp);

  int nc = coords.number_elements();

  for (int i = 0; i < nc; i++)
  {
    const Coords_and_Stuff * ci = coords[i];

    tmp[ci->query_number()]++;

//  cerr << "Item " << i << " is a match to query " << ci->query_number() << endl;
  }

  if (tag.length())
  {
    for (int i = 0; i < nfeatures; ++i)
    {
      if (0 == tmp[i])
        continue;

      int c = static_cast<int>(tmp[i] * bit_count_scaling_factor) + 1;
        sfc.hit_bit(i, c);
    }
  }
  else
  {
    for (int i = 0; i < nfeatures; i++)
    {
      buffer << ' ' << tmp[i];
    }
  }

  int bit_offset = nfeatures;

  if (do_topological_distances)
  {
    if (tag.length())
    {
      for (int i = 0; i < nfeatures; i++)
      {
        for (int j = i; j < nfeatures; j++)
        {
          const int ndx = i * nfeatures + j;
          const Distances<int> & d = bond_distances[ndx];
          d.set_bits(bit_offset + 3 * ndx, bit_count_scaling_factor, sfc);
        }
      }
    }
    else
    {
      for (int i = 0; i < nfeatures; i++)
      {
        for (int j = i; j < nfeatures; j++)
        {
          const Distances<int> & d = bond_distances[i * nfeatures + j];
          d.append_results_integer(buffer);
        }
      }
    }
  }

  bit_offset += (3 * nfeatures * nfeatures);

  if (do_spatial_distances)
  {
    if (tag.length())
    {
      for (int i = 0; i < nfeatures; i++)
      {
        for (int j = i; j < nfeatures; j++)
        {
          const int ndx = i * nfeatures + j;
          Distances<coord_t> & d = spatial_distances[ndx];
          d.set_bits(bit_offset + 3 * ndx, bit_count_scaling_factor, sfc);
        }
      }
    }
    else
    {
      for (int i = 0; i < nfeatures; i++)
      {
        for (int j = i; j < nfeatures; j++)
        {
          Distances<coord_t> & d = spatial_distances[i * nfeatures + j];
          d.append_results_float(buffer);
        }
      }
    }
  }

  bit_offset += (3 * nfeatures * nfeatures);

  if (compute_spatial_topological_ratio)
  {
    if (tag.length())
    {
      for (int i = 0; i < nfeatures; i++)
      {
        for (int j = i; j < nfeatures; j++)
        {
          const int ndx = i * nfeatures + j;
          Distances<float> & d = spatial_bond_ratio[ndx];
          d.set_bits(bit_offset + 3 * ndx, bit_count_scaling_factor * 15.0f, sfc);   // these are ratios, so we need to expand them
        }
      }
    }
    else
    {
      for (int i = 0; i < nfeatures; i++)
      {
        for (int j = i; j < nfeatures; j++)
        {
          Distances<float> & d = spatial_bond_ratio[i * nfeatures + j];
          d.append_results_float(buffer);
        }
      }
    }
  }

  return 1;
}

static int
do_compute_the_descriptors (Molecule & m,
                            const int nfeatures,
                            const resizable_array_p<Coords_and_Stuff> & coords)
{
  if (do_topological_distances && do_spatial_distances)
    return do_spatial_and_topological_distances_computation(m, nfeatures, coords);

  if (do_topological_distances)
    return do_topological_distances_computation(m, nfeatures, coords);

  if (do_spatial_distances)
    return do_spatial_distances_computation(m, nfeatures, coords);

  cerr << "Huh, what am I supposed to do?\n";

  return 0;
}

template <typename T>
int
write_all_distances (Molecule & m,
                     const int nfeatures,
                     Distances<T> * distances,
                     IWString_and_File_Descriptor & output)
{
  if (nullptr == distances)
    return 0;

  int name_written = 0;

  for (int i = 0; i < nfeatures; i++)
  {
    for (int j = i; j < nfeatures; j++)
    {
      Distances<T> & d = distances[i * nfeatures + j];

      if (d.n())
      {
        if (! name_written)
        {
          append_first_token_of_name(m.name(), output);
          output << '\n';
          name_written = 1;
        }

        output << i << ' ' << j;

        d.write_all_results(output);
        output << "\n";
      }
    }
  }

  if (name_written)
    output << "|\n";

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

static int
distance_between_features (Molecule & m,
                           const isotope_t * isotopes,
                           const int nfeatures,
                           IWString & output)
{
  (void) m.ring_membership();    // force ring perception

  if (do_topological_distances || ignore_atoms_closer_than_2d)
    m.recompute_distance_matrix();

  resizable_array_p<Coords_and_Stuff> coords;
  coords.resize(nfeatures * 10);

  if (! identify_matched_atoms(m, isotopes, coords))
    return 0;

// Compute the descriptors

  int nd = do_compute_the_descriptors(m, nfeatures, coords);

  if (stream_for_all_distances.is_open())
  {
    if (0 == nd)   // no distances found
      return 1;

    if (nullptr != spatial_distances)
      return write_all_distances(m, nfeatures, spatial_distances, stream_for_all_distances);
    if (nullptr != bond_distances)
      return write_all_distances(m, nfeatures, bond_distances, stream_for_all_distances);

    return 1;
  }

  Sparse_Fingerprint_Creator sfc;    // might not be used

  if (0 == tag.length())
    append_first_token_of_name(m.name(), output);

  gather_results(coords, sfc, nfeatures, output);

  if (tag.length())
  {
    IWString tmp;
    sfc.daylight_ascii_form_with_counts_encoded(tag, output);

    output << tmp << '\n';
  }
  else
    output << '\n';

  return 1;
}

// If isotopes are encountered, they will be stored in `isotopes`.
static int
preprocess(Molecule & m,
           std::unique_ptr<isotope_t[]>& isotopes) {
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (element_transformations.active())
    (void) element_transformations.process(m);

  if (charge_assigner.active())
    charge_assigner.process(m);

  if (donor_acceptor_assigner.active())
  {
    if (queries.number_elements())    // may contain an isotopic query, so shield the isotopes
    {
      isotopes = m.GetIsotopes();
      m.transform_to_non_isotopic_form();
    }

    donor_acceptor_assigner.process(m);
  }

  return 1;
}

static void
check_for_coordinates (const Molecule & m)
{
  int d = m.highest_coordinate_dimensionality();

  if (d < 3)
  {
    cerr << "Spatial distances requested, but first molecule has <3 dimensions.\n";
    cerr << "Spatial distances not computed\n";

    do_spatial_distances = 0;
  }

  return;
}

static int
distance_between_features (data_source_and_type<Molecule> & input,
                           const int nfeatures,
                           IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (verbose > 1)
      cerr << molecules_read << " '" << m->name() << "'\n";

    if (do_spatial_distances && 1 == molecules_read)
      check_for_coordinates(*m);

    std::unique_ptr<isotope_t[]> isotopes;
    preprocess(*m, isotopes);

    if (! distance_between_features(*m, isotopes.get(), nfeatures, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
distance_between_features (const char * fname,
                           FileType input_type,
                           const int nfeatures,
                           IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return distance_between_features(input, nfeatures, output);
}

static int
distance_between_features_filter_process (Molecule & m,
                                          const int nfeatures,
                                          IWString_and_File_Descriptor & output)
{
  std::unique_ptr<isotope_t[]> isotopes;

  preprocess(m, isotopes);

  return distance_between_features(m, isotopes.get(), nfeatures, output);
}

static int
distance_between_features_filter_process (const const_IWSubstring & buffer,
                                          const int nfeatures,
                                          IWString_and_File_Descriptor & output)
{
  Molecule m;

  if (! m.build_from_smiles(buffer))
  {
    cerr << "Cannot parse smiles '" << buffer << "'\n";
    return 0;
  }

  return distance_between_features_filter_process(m, nfeatures, output);
}

static int
distance_between_features_filter (iwstring_data_source & input,
                                  const int nfeatures,
                                  IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(8192);

    if (! buffer.starts_with(smiles_tag))
      continue;

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();

    if (! distance_between_features_filter_process(buffer, nfeatures, output))
      return 0;
  }

  output.flush();

  return 1;
}

static int
distance_between_features_filter (const char * fname,
                                  const int nfeatures,
                                  IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return distance_between_features_filter(input, nfeatures, output);
}

static int
create_from_smarts (const_IWSubstring & smarts,
                    resizable_array_p<Substructure_Query> & queries)
{
  Substructure_Query * q = new Substructure_Query;

  if (! q->create_from_smarts(smarts))
  {
    delete q;

    cerr << "Cannot parse smarts '" << smarts << "'\n";
    return 0;
  }

  queries.add(q);

  return 1;
}

/*
  Appends all the info for a descriptor combination
*/

static void
append_keeper_header (IWString & header,
                      const char * stem,
                      int n,
                      char c[2])
{
  header << stem << "min" << c[0] << c[1];
  header << stem << "max" << c[0] << c[1];
  header << stem << "ave" << c[0] << c[1];

//cerr << "Processing '" << c[0] << c[1] << "' n = " << nkeep << endl;

  for (int i = 0; i < n; i++)
  {
    header << stem << c[0] << c[1] << i;
  }

  return;
}

static char alphabet [] = {'a',
                          'b',
                          'c',
                          'd',
                          'e',
                          'f',
                          'g',
                          'h',
                          'i',
                          'j',
                          'k',
                          'l',
                          'm',
                          'n',
                          'o',
                          'p',
                          'q',
                          'r',
                          's',
                          't',
                          'u',
                          'v',
                          'w',
                          'x',
                          'y',
                          'z',
                          '_',
                          '0',
                          '1',
                          '2',
                          '3',
                          '4',
                          '5',
                          '6',
                          '7',
                          '8',
                          '9'};

static int
build_header (IWString & header,
              const int nfeatures,
              const char * stem)
{
  assert (' ' == stem[0]);

  for (int i = 0; i < nfeatures; i++)
  {
    char s[2];
    s[0] = alphabet[i];

    for (int j = i; j < nfeatures; j++)
    {
      s[1] = alphabet[j];
      append_keeper_header(header, stem, nkeep, s);
    }
  }

  return 1;
}

static int
build_header (IWString & header,
              const int nfeatures)
{
  for (int i = 0; i < nfeatures; i++)
  {
    if (i > 0)
      header << ' ';

    header << "dbf_n" << i;
  }

  if (do_topological_distances)
  {
    build_header(header, nfeatures, " dbf_");
  }

  if (do_spatial_distances)
  {
    build_header(header, nfeatures, " dbf_3");
  }

  if (compute_spatial_topological_ratio)
    build_header(header, nfeatures, " dbf_r");

  return 1;
}

static void
display_hard_coded_features(char flag,
                            std::ostream & os)
{
  os << "  -" << flag << " arom5     perceive 5 membered aromatic rings\n";
  os << "  -" << flag << " arom6     perceive 6 membered aromatic rings\n";
  os << "  -" << flag << " arom      perceive 5 and 6 membered aromatic rings\n";

  return;
}

static int
distance_between_features (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vi:E:A:t:g:lp:q:s:z:N:H:b:du:nk:T:V:C:rW:J:y:eL:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl))
  {
    usage(2);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
  {
    cerr << "Cannot process aromaticity options (-A)\n";
    usage(5);
  }

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('t'))
  {
    if (! element_transformations.construct_from_command_line(cl, verbose, 't'))
    {
      cerr << "Cannot parse -t option\n";
      usage(11);
    }
  }

  if (cl.option_present('N'))
  {
    if (! charge_assigner.construct_from_command_line(cl, verbose > 1, 'N'))
    {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      usage(4);
    }
  }

  if (cl.option_present('H'))
  {
    if (! donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose))
    {
      cerr << "Cannot initialise donor/acceptor assigner (-H option)\n";
      usage(6);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (! cl.option_present('p'))
  {
    do_topological_distances = 1;
    if (verbose)
      cerr << "Doing topological distances by default\n";
  }
  else
  {
    const_IWSubstring p;
    int i = 0;
    while (cl.value('p', p, i++))
    {
      if ("2d" == p)
      {
        do_topological_distances = 1;
        if (verbose)
          cerr << "Will compute topological distances\n";
      }
      else if ("3d" == p)
      {
        do_spatial_distances = 1;
        if (verbose)
          cerr << "Will compute spatial distances\n";
      }
      else
      {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(9);
      }
    }
  }

// Since this programme might be wrapped in a shel script, take the last value
// for the -u option, as that's something people may want to change

  if (cl.option_present('u'))
  {
    int i = 0;
    while (cl.value('u', missing_value, i++))
    {
    }

    if (verbose)
      cerr << "Missing values written as '" << missing_value << "'\n";
  }

  if (cl.option_present('z'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('z', m, i++))
    {
      if ("first" == m)
      {
        take_first_of_multiple_hits = 1;
        if (verbose)
          cerr << "Will take the first of multiple hits for user query\n";
      }
      else if ("ignore" == m)
      {
        ignore_queries_hitting_multiple_times = 1;
        if (verbose)
          cerr << "Will ignore multiple query matches in the user query\n";
      }
      else if ("oknomatch" == m)
      {
      }
      else
      {
        cerr << "Unrecognised -z qualifier '" << m << "'\n";
        usage(17);
      }
    }
  }

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring smarts;
    while (cl.value('s', smarts, i++))
    {
      if (! create_from_smarts(smarts, queries))
      {
        cerr << "Cannot parse -s option '" << smarts << "'\n";
        return 1;
      }
    }
  }

  if (cl.option_present('q'))
  {
    int i = 0;
    const_IWSubstring q;

    while (cl.value('q', q, i++))
    {
      if (! process_cmdline_token('q', q, queries, verbose))
      {
        cerr << "Cannot parse -q option '" << q << "'\n";
        return 4;
      }
    }
  }

  nq = queries.number_elements();

  if (verbose)
    cerr << "Defined " << nq << " user defined queries\n";

  if (cl.option_present('e'))
  {
    if (! charge_assigner.active() && ! donor_acceptor_assigner.active())
    {
      cerr << "Warning, request for default features, but charge assigner and donor acceptor assignor both inactive! Repeat -e to allow\n";
      if (1 == cl.option_count('e'))
        return 1;
    }

    positive_negative_acceptor_donor = 1;

    if (verbose)
      cerr << "Will do default set of queries, positive, negative, acceptor, donor\n";
  }

  if (positive_negative_acceptor_donor)
  {
    feature_name.add(new IWString("POS"));
    feature_name.add(new IWString("NEG"));
    feature_name.add(new IWString("ACC"));
    feature_name.add(new IWString("DON"));
  }

  if (positive_negative_acceptor_donor && nq)
  {
    query_to_feature_number.resize(nq);
    if (positive_negative_acceptor_donor)
    {
      for (int i = 0; i < nq; ++i)
      {
        query_to_feature_number.add(4 + i);
        feature_name.add(new IWString(queries[i]->comment()));
      }
    }
  }
  else if (nq)
  {
    for (int i = 0; i < nq; ++i)
    {
      query_to_feature_number.add(i);
      feature_name.add(new IWString(queries[i]->comment()));
    }
  }

  if (cl.option_present('W'))
  {
    int i = 0;
    const_IWSubstring w;
    while (cl.value('W', w, i++))
    {
      if ("arom5" == w)
      {
        hard_coded_query_centre_of_aromatic_5_membered_ring = (nq + 4 * positive_negative_acceptor_donor) + number_hard_coded_queries;
        if (verbose)
          cerr << "Will perceive aromatic 5 membered rings, feature " << number_hard_coded_queries << endl;
        number_hard_coded_queries++;
        feature_name.add(new IWString("AROM5"));
      }
      else if ("arom6" == w)
      {
        hard_coded_query_centre_of_aromatic_6_membered_ring = nq + 4 * positive_negative_acceptor_donor + number_hard_coded_queries;
        if (verbose)
          cerr << "Will perceive aromatic 6 membered rings, feature " << number_hard_coded_queries << endl;
        number_hard_coded_queries++;
        feature_name.add(new IWString("AROM6"));
      }
      else if ("arom" == w)
      {
        hard_coded_query_centre_of_aromatic_5_membered_ring = nq + 4 * positive_negative_acceptor_donor + number_hard_coded_queries;
        number_hard_coded_queries++;
        feature_name.add(new IWString("AROM5"));
        hard_coded_query_centre_of_aromatic_6_membered_ring = nq + 4 * positive_negative_acceptor_donor + number_hard_coded_queries;
        number_hard_coded_queries++;
        feature_name.add(new IWString("AROM6"));
        if (verbose)
          cerr << "Will perceive 5 and 6 sized aromatic rings, features " << hard_coded_query_centre_of_aromatic_5_membered_ring << " and " << hard_coded_query_centre_of_aromatic_6_membered_ring << endl;
      }
      else if ("help" == w)
      {
        display_hard_coded_features('W', cerr);
        return 0;
      }
      else
      {
        cerr << "Unrecognised -W qualifier '" << w << "'\n";
        display_hard_coded_features('W', cerr);
        return 5;
      }
    }
  }

  int nfeatures = nq + number_hard_coded_queries + (4 * positive_negative_acceptor_donor);

  if (0 == nfeatures)
  {
    cerr << "No queries defined, cannot continue\n";
    usage(8);
  }

  if (verbose)
    cerr << "Defined " << nfeatures << " features\n";

  for (int i = 0; i < nq; i++)
  {
    queries[i]->set_find_unique_embeddings_only(1);
  }

// We initialise a temporary histogram to make sure the histogram specification is OK

  if (cl.option_present('T'))
  {
    histogram_specification = cl.string_value('T');

    IWHistogram tmp;
    if (! tmp.initialise(histogram_specification))
    {
      cerr << "INvalid histogram specification '" << histogram_specification << "'\n";
      return 0;
    }

    if (verbose)
    {
      cerr << "Histograms initialised \n";
      tmp.debug_print(cerr);
    }
  }

  iwdigits.set_include_leading_space(1);

  if (! iwdigits.initialise(100))
  {
    cerr << "Cannot initialise digits\n";
    return 4;
  }

  fraction_as_string.set_include_leading_space(1);

  if (! fraction_as_string.initialise(0.0, 30.0, 4))
  {
    cerr << "Cannot initialise string floating point representations\n";
    return 3;
  }

  if (cl.option_present('r'))
  {
    compute_spatial_topological_ratio = 1;
    do_topological_distances = 1;
    do_spatial_distances = 1;

    if (verbose)
      cerr << "Will compute the spatial/topological dstance ratio\n";
  }

#ifdef DBF_MOLVOL
  int vdw_type = 0;
  if (cl.option_present('V'))
  {
    if (! set_default_van_der_waals_radius_type(cl, 'V', vdw_type, verbose > 1))
    {
      cerr << "Cannot set vdw radius type (-V option)\n";
      usage(8);
    }

    surface_area_molvol.set_vdw_type(vdw_type);
  }
#endif

// Since we don't want to repeat the shortest distance, we increase nkeep silently
// and then don't bother writing out the shortest value

  if (cl.option_present('k'))
  {
    if (! cl.value('k', nkeep) || nkeep < 0)
    {
      cerr << "The number of descriptors option (-k) must be >= 0\n";
      usage(13);
    }

    if (verbose)
      cerr << "Descriptors will be the " << nkeep << " shortest distances\n";

    nkeep++;
  }

// Now that nkeep is known, construct the various missing value strings

  all_values_missing << ' ' << missing_value << ' ' << missing_value << ' ' << missing_value;
  for (int i = 0; i < nkeep; i++)
  {
    all_values_missing << ' ' << missing_value;
  }

  if (do_topological_distances)
    bond_distances = new Distances<int> [nfeatures * nfeatures];

  std::unique_ptr<Distances<int>[]> free_bond_distances(bond_distances);

  if (do_spatial_distances)
    spatial_distances = new Distances<coord_t> [nfeatures * nfeatures];

  std::unique_ptr<Distances<coord_t>[]> free_spatial_distances(spatial_distances);

  if (do_spatial_distances && do_topological_distances)
    spatial_bond_ratio = new Distances<float> [nfeatures * nfeatures];     // actually only used if the -r option is set.

  std::unique_ptr<Distances<float>[]> free_spatial_bond_ratio(spatial_bond_ratio);

  if (cl.option_present('b'))
  {
    if (! cl.value('b', ignore_atoms_closer_than_2d) || ignore_atoms_closer_than_2d < 1)
    {
      cerr << "The ignore atom pairs cutoff value (-b option) must be a whole positive number\n";
      usage(5);
    }

    if (verbose)
      cerr << "Will ignore atom pairs " << ignore_atoms_closer_than_2d << " or fewer bonds apart\n";
  }

  if (cl.option_present('d'))
  {
    consider_bonded_features = 1;

    if (verbose)
      cerr << "Will consider bonded features\n";
  }

  if (cl.option_present('J'))
  {
    cl.value('J', tag);
    if (! tag.ends_with('<'))
      tag << '<';

    if (verbose)
      cerr << "Will act as a TDT filter\n";

    if (cl.option_present('y'))
    {
      if (! cl.value('y', bit_count_scaling_factor) || bit_count_scaling_factor <= 0.0f)
      {
        cerr << "The bit count scaling factor (-y) must be +ve floating point number\n";
        usage(2);
      }

      if (verbose)
        cerr << "Counts will be scaled by " << bit_count_scaling_factor << " for bit counts in the fingerprint\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (tag.length())
    ;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  set_default_iwstring_float_concatenation_precision(4);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('C'))
  {
    if (0 == nkeep)
    {
      nkeep = 20;
      cerr << "Request for all distances -C, but no -k specified, default value " << nkeep << endl;
    }

    const char * c = cl.option_value('C');

    if (! stream_for_all_distances.open(c))
    {
      cerr << "Cannot open stream for all distances '" << c << "'\n";
      return 3;
    }

    if (verbose)
      cerr << "All distances written to '" << c << "'\n";
  }

  if (cl.option_present('L'))
  {
    const_IWSubstring l = cl.option_value('L');
    stream_for_labelled_atoms.add_output_type(FILE_TYPE_SMI);
    stream_for_labelled_atoms.add_output_type(FILE_TYPE_SDF);
    if (stream_for_labelled_atoms.would_overwrite_input_files(cl, l))
    {
      cerr << "Labelled atoms file -L " << l << " cannot overwrite input file(s)\n";
      return 1;
    }

    if (! stream_for_labelled_atoms.new_stem(l))
    {
      cerr << "Cannot open stream for labelled molecules '" << l << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Labelled molecules written to " << l << "'\n";
  }

  IWString_and_File_Descriptor output(1);

// Write the header for descriptor files

  if (cl.option_present('C'))
    ;
  else if (tag.length())
    ;
  else
  {
    IWString header;

    build_header(header, nfeatures);

    output << "Name " << header << '\n';
  }

  int rc = 0;
  if (tag.length())
  {
    if (! distance_between_features_filter(cl[0], nfeatures, output))
      rc = 2;
  }
  else
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! distance_between_features(cl[i], input_type, nfeatures, output))
      {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";

    for (int i = 0; i < matches_per_molecule.number_elements(); i++)
    {
      if (matches_per_molecule[i])
        cerr << matches_per_molecule[i] << " molecules matched " << i << " queries\n";
    }

    for (int i = 0; i < matches_per_query.number_elements(); ++i)
    {
      if (i >= feature_name.number_elements())
        break;

      cerr << " query " << i << " " << *(feature_name[i]) << ": " << molecules_matching_query[i] << " molecules matched, " << matches_per_query[i] << " sites\n";
    }
  }

  if (histogram_specification.length() && nullptr != spatial_distances)
  {
    for (int i = 0; i < nfeatures; i++)
    {
      for (int j = i; j < nfeatures; j++)
      {
        const Distances<float> & dist = spatial_distances[i * nfeatures + j];

        const IWHistogram * h = dist.global_histogram();

        assert (nullptr != h);

        cerr << "i = " << i << " j = " << j << endl;
        h->debug_print(cerr);
      }
    }
  }

  if (stream_for_all_distances.is_open())
    stream_for_all_distances.flush();

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = distance_between_features(argc, argv);

  return rc;
}
