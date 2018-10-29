/*
  Random mutations of smiles
*/

#include <iostream>
#include <time.h>
#include <memory>
#include <limits>
#include <random>
using std::cerr;
using std::endl;

#define REPORT_PROGRESS_IMPLEMENTATION

#include "cmdline.h"
#include "iw_stl_hash_set.h"
#include "accumulator.h"
#include "misc.h"
#include "report_progress.h"

#include "istream_and_type.h"
#include "molecule.h"
#include "target.h"
#include "smiles.h"
#include "path.h"
#include "substructure.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "numass.h"
#include "random_reactions.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static Number_Assigner number_assigner;

static int reduce_to_largest_fragment = 0;

static float probability_inter_molecule_change = 0.5;

static float probability_intra_molecule_change = 0.5;

static int max_chars_to_remove = 1;

static int random_replicates = 1;

static int number_times_single_character_swap = 1;

static int number_iterations = 0;

static int run_for = 0;

static time_t tzero;    // not initialised

static float probability_refresh_from_initial = 0.0;
static float probability_refresh_from_most_recent_valid = 0.0;

static int complete_refresh = 0;

typedef long long large_integer;

static large_integer smiles_produced = 0;

static large_integer molecular_interpretation_attempted = 0;

static large_integer stop_after_producing = std::numeric_limits<int>::max();

static large_integer valid_molecules_produced = 0;

static large_integer invalid_valence_produced = 0;

static int reject_for_non_periodic_table_elements = 1;

static int non_periodic_table_molecules_rejected = 0;

static int discard_non_organic_molecules = 0;

static int non_organic_molecules_discarded = 0;

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = std::numeric_limits<int>::max();

static large_integer rejected_for_too_small = 0;
static large_integer rejected_for_too_large = 0;

static int min_nrings = 0;
static int max_rings = std::numeric_limits<int>::max();

static large_integer rejected_for_too_few_rings = 0;
static large_integer rejected_for_too_many_rings = 0;

static int largest_ring_size_allowed = 0;

static int rejected_for_recreating_initial_structure = 0;

static large_integer smiles_written = 0;

static Report_Progress_Template<large_integer> report_progress;

static Accumulator_Int<int> acc_natoms;
static Accumulator_Int<int> acc_nrings;

static resizable_array_p<Substructure_Query> queries_to_match;
static resizable_array_p<Substructure_Query> queries_to_avoid;

static large_integer rejected_for_not_matching_needed_query = 0;
static large_integer rejected_for_matching_avoid_query = 0;

static int translate_non_organics = 0;

static extending_resizable_array<int> valid_molecules_produced_per_generation;

static IWString_and_File_Descriptor stream_for_all_valid_smiles;

/*
  Don't write out any of the starting structures
*/

static IW_STL_Hash_Set initial_structures;

static IWString name_stem;

static Random_Reactions rxn;

static double probability_do_reaction = 0.0;

static std::random_device(rd);

static std::mt19937_64 rng(rd());

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Random mutations of smiles strings\n";
  cerr << "  -N <n>        number of iterations to run\n";
  cerr << "  -n <n>        complete refresh from initial smiles every <n> iterations\n";
  cerr << "  -x <nchar>    max number of characters to remove during excisions\n";
  cerr << "  -p <n>        generate <n> random replicates of each starting molecule\n";
  cerr << "  -w <n>        number of single character swaps to do each time\n";
  cerr << "  -F <prob>     probability of refreshing each smiles from initial smiles\n";
  cerr << "  -f <prob>     probability of refreshing each smiles from most recent valid smiles\n";
  cerr << "  -c <n>        minimum number of atoms in molecules produced\n";
  cerr << "  -C <n>        maximum number of atoms in molecules produced\n";
  cerr << "  -r <n>        minimum number of rings in molecules produced\n";
  cerr << "  -R <n>        maximum number of rings in molecules produced\n";
  cerr << "  -z <n>        maximum ring size permitted in molecules produced\n";
  cerr << "  -S <fname>    sidechain library - file of smiles\n";
  cerr << "  -Y <stem>     name stem for molecules created\n";
  cerr << "  -q            queries that must be present in products\n";
  cerr << "  -Q            queries that must NOT be present in products\n";
  cerr << "  -b <n>        report progress every <n> molecules formed\n";
  cerr << "  -O ...        miscellaneous options, enter '-O help' for info\n";
  cerr << "  -a ...        number assigner options, enter '-a help' for info\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static double
random_number_between_01 ()
{
  std::uniform_real_distribution<double> u(0.0, 1.0);

  return u(rng);
}

/*
  Sidechains are inserted into smiles
*/

static resizable_array_p<IWString> sidechain_library;

static int
read_sidechain_library(iwstring_data_source & input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#') || 0 == buffer.length())
      continue;

    buffer.truncate_at_first(' ');

    IWString * t = new IWString(buffer);
    sidechain_library.add(t);
  }

  return sidechain_library.number_elements();
}

static int
read_sidechain_library(const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open sidechain library '" << fname << "'\n";
    return 0;
  }

  return read_sidechain_library(input);
}

static int
insert_member_of_sidechain_library(IWString & s,
                                   const IWString & sidechain)
{
  const int n = s.length();

  std::uniform_int_distribution<int> u(0, n - 1);

  int istart = u(rng);

  s.insert(sidechain, istart);

  return 1;
}

static int
insert_member_of_sidechain_library(IWString & s)
{
  std::uniform_int_distribution<int> u(0, sidechain_library.number_elements() - 1);

  const int n = u(rng);

  return insert_member_of_sidechain_library(s, *(sidechain_library[n]));
}

template <typename T>
void
swap_items (T * v, int n)
{
  if (n < 2)
    return;

  if (2 == n)
  {
    std::swap(v[0], v[1]);
    return;
  }

  std::uniform_int_distribution<int> u(0, n - 1);

  const int i1 = u(rng);

  int i2 = u(rng);
  while (i2 == i1)
  {
    i2 = u(rng);
  }

  std::swap(v[i1], v[i2]);

  return;
}

static void
shuffle_chars_within_smiles (IWString & s)
{
  const int n = s.length();

  if (1 == n)
    return;

  char * c = const_cast<char *>(s.rawchars());

  double r = random_number_between_01();

  if (r < 0.40)
  {
    swap_items(c, n);
  }
  else if (r < 0.70)
  {
    unsigned short * s = reinterpret_cast<unsigned short *>(c);
    swap_items(s, n / 2);
  }
  else if (r < 0.90)
  {
    unsigned int * s = reinterpret_cast<unsigned int *>(c);
    swap_items(s, n / 4);
  }
  else
  {
    double * s = reinterpret_cast<double *>(c);
    swap_items(s, n / 8);
  }

  return;
}

static int
excise_characters(IWString & s)
{
  int chars_to_remove = max_chars_to_remove;

  const int n = s.length();

  if (1 == n)
    return 0;

  if (chars_to_remove >= n)
    chars_to_remove = n - 1;

  std::uniform_int_distribution<int> u(0, n - 1);

  int cstart = u(rng);

  if (cstart + chars_to_remove >= n)
    s.iwtruncate(cstart + 1);
  else
    s.erase(cstart, cstart + chars_to_remove - 1);

  return 1;
}

static void
breed (IWString & s1,
       IWString & s2)
{
  std::uniform_int_distribution<int> u(0, s1.length() - 1);

  int start1 = u(rng);

  std::uniform_int_distribution<int>::param_type p(0, s2.length() - 1);

  u.param(p);
  const int start2 = u(rng);

  IWString tmp1;
  tmp1.strncpy(s1.rawchars(), start1 + 1);
  tmp1.strncat(s2.rawchars() + start2, s2.length() - start2 - 1);

  IWString tmp2;
  tmp2.strncpy(s2.rawchars(), start2 + 1);
  tmp2.strncat(s1.rawchars() + start1, s1.length() - start1 - 1);

  s1 = tmp1;
  s2 = tmp2;

  return;
}

static void
swap_some_characters (IWString & s1,
                      IWString & s2)
{
  int n = number_times_single_character_swap;
  if (number_times_single_character_swap > 1)
  {
    std::uniform_int_distribution<int> u(1, number_times_single_character_swap);
    n = u(rng);
  }

  std::uniform_int_distribution<int> u1(0, s1.length() - 1);
  std::uniform_int_distribution<int> u2(0, s2.length() - 1);

  for (int i = 0; i < n; i++)
  {
    int j1 = u1(rng);
    int j2 = u2(rng);

    char t = s1[j1];
    if ('[' == t || ']' == t)  // should also check s2
      continue;

    s1[j1] = s2[j2];
    s2[j2] = t;
  }

  return;
}

static int
matches_any_of_these (Molecule & m,
                      resizable_array_p<Substructure_Query> & query)
{
  int n = query.number_elements();

  Molecule_to_Match target(&m);

  for (int i = 0; i < n; i++)
  {
    if (query[i]->substructure_search(target))
      return 1;
  }

  return 0;
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();

// Truncate name to first token only

  IWString mname(m.name());

  mname.truncate_at_first(' ');

  if (mname != m.name())
    m.set_name(mname);

  return;
}

//static resizable_array_p<IWString> smiles;
static resizable_array_p<IWString> initial_name;

static int
fill_smiles_array (data_source_and_type<Molecule> & input,
                   resizable_array_p<IWString> & smiles)
{
  IW_STL_Hash_Set already_seen;

  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (m->natoms() < 2)
    {
      cerr << "Cannot process molecules with 1 atom, ignoring '" << m->name() << "'\n";
      continue;
    }

    initial_structures.insert(m->unique_smiles());

    for (int i = 0; i < random_replicates; i++)
    {
      IWString * s = new IWString(m->random_smiles());
      if (already_seen.contains(*s))
        delete s;
      else
      {
        already_seen.insert(*s);
        smiles.add(s);
        IWString * n = new IWString(m->name());
        initial_name.add(n);
      }
    }
  }

  assert (initial_name.number_elements() == smiles.number_elements());

  return smiles.number_elements();
}

static int
fill_smiles_array (const char * fname, int input_type,
                   resizable_array_p<IWString> & smiles)
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

  return fill_smiles_array(input, smiles);
}

static int
do_translate_non_organics (Molecule & m)
{
  if (m.organic_only())
    return 0;

  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    const Element * e = m.elementi(i);

    if (e->organic())
      continue;

    if (0.5 < random_number_between_01())
      m.set_atomic_number(i, 6);
    else
      m.set_atomic_number(i, 7);

    rc++;
  }

  return rc;
}

static int
contains_multi_connected_halogen (const Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    const atomic_number_t z = a->atomic_number();

    if (6 == z || 7 == z || 8 == z)
      continue;

    if (9 == z || 17 == z || 35 == z || 53 == z)
    {
      if (a->ncon() > 1)
        return 1;
    }
  }

  return 0;    // no multi-connected halogens found
}

static int
contains_out_of_range_formal_charge (const Molecule & m)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    formal_charge_t fc = m.formal_charge(i);

    if (0 == fc)
      ;
    else if (1 == fc)
      ;
    else if (-1 == fc)
      ;
    else
      return 1;
  }

  return 0;   // no strange charges found
}

static int
try_various_subsets_of_the_smiles (Molecule & m,
                                   IWString & s)
{
  if (s.length() < 30)
    return 0;

  const_IWSubstring h;
  s.from_to(0, s.length() / 2, h);
  if (m.build_from_smiles(h))
  {
    s.iwtruncate(s.length() / 2 + 1);
    return 1;
  }

  s.from_to(s.length() / 2, s.length() - 1, h);
  if (m.build_from_smiles(h))
  {
    s.remove_leading_chars(s.length() / 2);
    return 1;
  }

  s.remove_leading_chars(s.length() / 4);

  return m.build_from_smiles(s);
}

/*
  For ring opening and closing matching, there must be an equal
  number of digits in the smiles, so the number of each digit
  must be equal
*/

static void
digits_must_match (IWString & s)
{
  int n = s.number_elements();

  int dcount[10];

  set_vector(dcount, 10, 0);

  const char * r = s.rawchars();

  int inside_square_bracket = 0;

  for (int i = 0; i < n; i++)
  {
    char c = r[i];

    if ('[' == c)
      inside_square_bracket = 1;
    else if (']' == c)
      inside_square_bracket = 0;
    else if (inside_square_bracket)
      ;
    else if (isdigit(c))
    {
      int d = c - '0';
      dcount[d]++;
    }
  }

  for (int i = 1; i < 10; i++)
  {
    if (0 == dcount[i])
      continue;

    if (0 == dcount[i] % 2)
      continue;

    s << i;
  }

  return;
}

/*
  Look for two character strings. When encountered, change the
  first character to sto
*/

static void
do_translate2 (IWString & s,
               const char * sfrom,
               const char sto)
{
  char * r = const_cast<char *>(s.rawchars());

  int n = s.length();

  char f1 = sfrom[0];
  char f2 = sfrom[1];

  for (int i = 1; i < n; i++)
  {
    if (f2 != r[i])
      continue;

    if (f1 == r[i - 1])
      r[i - 1] = sto;
  }
}

static void
parenthesis_munging (IWString & s,
                     char oparen, char cparen)
{
  const char * r = s.rawchars();

  int n = s.length();

  int paren_level = 0;

  for (int i = 0; i < n; i++)
  {
    if (oparen == r[i])
    {
      if (i > 0 && oparen == r[i-1])
        s[i] = 'C';
      else
        paren_level++;
    }
    else if (cparen == r[i])
    {
      if (paren_level <= 0)
      {
        s[i] = oparen;
        paren_level++;
      }
      else if (i > 0 && oparen == r[i - 1])
        s[i] = 'N';
      else
        paren_level--;
    }
  }

  for (int i = 0; i < paren_level; i++)
  {
    s += cparen;
  }

  return;
}

/*
  Checks for common errors, and makes attempts at fixing some obvious errors
*/

static int
valid_smiles(IWString & s,
             IWString & smiles)
{
  smiles_produced++;

  int n = s.length();

  if (lower_atom_count_cutoff > 0 && n < lower_atom_count_cutoff)
    return 0;

  const char * sraw = s.rawchars();

  for (int i = 1; i < n; i++)
  {
    if ('l' == sraw[i])    // Only chlorine
    {
      if ('C' != sraw[i - 1])
        s[i - 1] = 'C';
    }
    else if ('r' == sraw[i])
    {
      if ('B' != sraw[i - 1])   // Only Br allowed
        s[i - 1] = 'B';
    }
  }

  char s0 = sraw[0];

  if (! isalpha(s0))
    s[0] = 'C';
  else if ('(' == s0)
    s[0] = 'N';
  else if ('#' == s0)
    s[0] = 'C';
  else if ('=' == s0)
    s[0] = 'C';
  else if ('r' == s0)
    s[0] = 'O';
  else if ('l' == s0)
    s[0] = 'C';

// Fix any ==

  for (int i = 1; i < n; i++)
  {
    if ('=' != sraw[i])
      continue;

    if ('=' == sraw[i - 1])
      s[i] = 'C';
  }

  do_translate2(s, "[(", 'C');

  do_translate2(s, "[)", 'C');

  do_translate2(s, "=)", 'C');

  do_translate2(s, "=(", 'O');

  do_translate2(s, "#=", 'N');

  do_translate2(s, "=#", 'C');

  do_translate2(s, "(1", 'C');

  do_translate2(s, ")1", 'C');

  do_translate2(s, "[=", 'c');

  do_translate2(s, "[#", 'N');

  do_translate2(s, "[+", 'C');

  do_translate2(s, "[-", 'O');

  parenthesis_munging(s, '(', ')');
  parenthesis_munging(s, '[', ']');

  digits_must_match(s);

  molecular_interpretation_attempted++;

// the logic is somewhat wrong here, we might reject something
// here that is a valid smiles, but it will not be written to
// stream_for_all_valid_smiles. But we want the speed of this...

  int matoms = count_atoms_in_smiles(s);

  if (lower_atom_count_cutoff > 0 && matoms < lower_atom_count_cutoff)
  {
    rejected_for_too_small++;
    return 0;
  }

  if (matoms > upper_atom_count_cutoff)
  {
    rejected_for_too_large++;
    return 0;
  }

  Molecule m;
  if (m.build_from_smiles(s))   // great
    ;
  else if (try_various_subsets_of_the_smiles(m, s))
    ;
  else
    return 0;

  valid_molecules_produced++;

  if (stream_for_all_valid_smiles.is_open())
  {
    stream_for_all_valid_smiles << m.smiles() << '\n';
    stream_for_all_valid_smiles.write_if_buffer_holds_more_than(4096);
  }

  if (translate_non_organics)
    do_translate_non_organics(m);

  int nr = m.nrings();

  if (min_nrings > 0 && nr < min_nrings)
  {
    rejected_for_too_few_rings++;
    return 0;
  }

  if (nr > max_rings)
  {
    rejected_for_too_many_rings++;
    return 0;
  }

  if (! m.valence_ok())
  {
    invalid_valence_produced++;
    return 0;
  }

  if (discard_non_organic_molecules && ! m.organic_only())
  {
    non_organic_molecules_discarded++;
    return 0;
  }

  if (reject_for_non_periodic_table_elements && m.contains_non_periodic_table_elements())
  {
    non_periodic_table_molecules_rejected++;
    return 0;
  }

  if (contains_out_of_range_formal_charge(m))
  {
    invalid_valence_produced++;
    return 0;
  }

  if (contains_multi_connected_halogen(m))
  {
    invalid_valence_produced++;
    return 0;
  }

  if (largest_ring_size_allowed > 0 && nr > 0)
  {
    for (int i = nr - 1; i >= 0; i--)
    {
      const Ring * ri = m.ringi(i);

      if (ri->number_elements() > largest_ring_size_allowed)
        return 0;
    }
  }

  if (queries_to_match.number_elements() && ! matches_any_of_these(m, queries_to_match))
  {
    rejected_for_not_matching_needed_query++;
    return 0;
  }

  if (queries_to_avoid.number_elements() && matches_any_of_these(m, queries_to_avoid))
  {
    rejected_for_matching_avoid_query++;
    return 0;
  }

  const IWString & usmi = m.unique_smiles();

  if (initial_structures.contains(usmi))
  {
    rejected_for_recreating_initial_structure++;
    return 0;
  }

  m.invalidate_smiles();

  smiles = m.smiles();

  if (verbose)
  {
    acc_natoms.extra(m.natoms());
    acc_nrings.extra(nr);
  }

  return 1;
}

static int
do_intra_molecular_change(IWString & s)
{
  random_number_t r = random_number_between_01();

  if (r < 0.33)
    excise_characters(s);
  else if (r < 0.66)
    shuffle_chars_within_smiles(s);
  else if (sidechain_library.number_elements())
    insert_member_of_sidechain_library(s);

  return 1;
}

static int
do_inter_molecular_change(IWString & s1,
                          IWString & s2)
{
  if (s1.length() < 4 || s2.length() < 4)   // 4 is just arbitrary
    swap_some_characters(s1, s2);
  else if (random_number_between_01() < 0.5)
    breed(s1, s2);
  else
    swap_some_characters(s1, s2);

  return 1;
}

static int
do_refresh (resizable_array_p<IWString> & smiles,
            const IWString * r,
            int * times_since_last_valid_smiles,
            random_number_t p)
{
  int n = smiles.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    if (0 == times_since_last_valid_smiles[i])
      continue;

    if (random_number_between_01() > p)
    {
      *(smiles[i]) = r[i];
      times_since_last_valid_smiles[i] = 0;
      rc++;
    }
  }

  if (verbose > 2)
    cerr << "Refresh processed " << rc << " molecules\n";

  return rc;
}

static int
something_other_than(int n,
                     int i)
{
  std::uniform_int_distribution<int> u(0, n - 1);

  int j = u(rng);

  while (j == i)
  {
    j = u(rng);
  }

  return j;
}

static int
need_to_stop_for_any_reason()
{
  if (run_for > 0)
  {
    time_t tnow = time(NULL);

    if (tnow - tzero > run_for)
    {
      if (verbose)
        cerr << "Computation complete after " << run_for << " seconds\n";
      return 1;
    }
  }

  if (valid_molecules_produced > stop_after_producing)
  {
    if (verbose)
      cerr << "Computation complete having produced " << valid_molecules_produced << " molecules, " <<stop_after_producing << endl;
    return 1;
  }

  return 0;
}

static int
smiles_mutation (resizable_array_p<IWString> & smiles,
                 IWString_and_File_Descriptor & output)
{
  int n = smiles.number_elements();

  assert (n > 0);

  if (1 == n && probability_inter_molecule_change > 0.0)
  {
    cerr << "Only one string present, cannot do inter molecular changes\n";
    probability_inter_molecule_change = 0.0;
  }

  IWString * initial_smiles = new IWString[n]; std::unique_ptr<IWString[]> free_initial_smiles(initial_smiles);
  IWString * last_valid_smiles = new IWString[n]; std::unique_ptr<IWString[]> free_last_valid_smiles(last_valid_smiles);
  int * times_since_last_valid_smiles = new_int(n); std::unique_ptr<int[]> free_times_since_last_valid_smiles(times_since_last_valid_smiles);
  int * changed_this_iteration = new int[n]; std::unique_ptr<int[]> free_changed_this_iteration(changed_this_iteration);

  for (int i = 0; i < n; i++)
  {
    initial_smiles[i] = *(smiles[i]);
    last_valid_smiles[i] = *(smiles[i]);
  }

// For useful info we need to keep track of iterations from restarting

  int generations_since_last_complete_refresh = 0;

  int iteration;

  for (iteration = 0; iteration < number_iterations; iteration++)
  {
    if (verbose > 1 && iteration > 0)
      cerr << "Begin iteration " << iteration << ", produced " << smiles_produced << " smiles, " << valid_molecules_produced << " valid molecules, " << smiles_written << " written\n";

    set_vector(changed_this_iteration, n, 0);

    if (0 == iteration)
      ;
    else if (probability_refresh_from_initial > 0.0 && probability_refresh_from_initial > random_number_between_01())
    {
      do_refresh(smiles, initial_smiles, times_since_last_valid_smiles, probability_refresh_from_initial);
      generations_since_last_complete_refresh = 0;
    }
    else if (probability_refresh_from_most_recent_valid > 0.0 && probability_refresh_from_most_recent_valid > random_number_between_01())
      do_refresh(smiles, last_valid_smiles, times_since_last_valid_smiles, probability_refresh_from_most_recent_valid);
    else if (complete_refresh > 0 && 0 == iteration % complete_refresh)
    {
      do_refresh(smiles, initial_smiles, times_since_last_valid_smiles, 0.0);
      generations_since_last_complete_refresh = 0;
    }

    for (int j = 0; j < n; j++)
    {
      if (0 == smiles[j]->length())
        continue;

      if (random_number_between_01() > probability_intra_molecule_change)
      {
        if (smiles[j]->length() < 3)
          continue;

        do_intra_molecular_change(*(smiles[j]));
        times_since_last_valid_smiles[j]++;
        changed_this_iteration[j] = 1;
      }

      if (random_number_between_01() > probability_inter_molecule_change)
      {
        int k = something_other_than(n, j);

        if (0 == smiles[k]->length())
          continue;

        do_inter_molecular_change(*(smiles[j]), *(smiles[k]));
        times_since_last_valid_smiles[j]++;
        times_since_last_valid_smiles[k]++;
        changed_this_iteration[j] = 1;
        changed_this_iteration[k] = 1;
      }

      if (rxn.active() && random_number_between_01() < probability_do_reaction)
      {
        IWString & t = *(smiles[j]);
        rxn.perform_random_reaction(t);
      }
    }

    generations_since_last_complete_refresh++;

    int produced_this_generation = 0;
    int valid_molecules_this_generation = 0;

    IWString smiles_from_molecule;   // scope here for efficiency

    for (int j = 0; j < n; j++)
    {
      if (! changed_this_iteration[j])
        continue;

      produced_this_generation++;

      if (! valid_smiles(*(smiles[j]), smiles_from_molecule))
        continue;

      output << smiles_from_molecule << ' ';

      if (name_stem.length())
      {
        output << name_stem;
        output << smiles_written << '\n';
      }
      else if (number_assigner.active())
      {
        IWString tmp(*(initial_name[j]));
        number_assigner.process(tmp);
        output << tmp << '\n';
      }
      else
        output << *(initial_name[j]) << '\n';

      smiles_written++;
      times_since_last_valid_smiles[j] = 0;
      valid_molecules_this_generation++;

      if (report_progress())
        cerr << "Generated " << smiles_produced << " wrote " << smiles_written << " valid smiles\n";
    }

    if (verbose > 1)
      cerr << "Iteration " << iteration << " generated " << produced_this_generation << " molecules, " << valid_molecules_this_generation << " valid\n";

    if (valid_molecules_this_generation > 0)
      valid_molecules_produced_per_generation[generations_since_last_complete_refresh] += valid_molecules_this_generation;

    output.write_if_buffer_holds_more_than(4096);

    if (verbose > 1)
      cerr << " iteration = " << iteration << " number_iterations = " << number_iterations << endl;

    if (need_to_stop_for_any_reason())
      break;
  }

  if (verbose)
    cerr << "Returning from " << iteration << " outer loop iterations\n";

  return 1;
}

static void
display_dash_o_options (std::ostream & os)
{
  os << " -O pia=<prob>         probability of an intra molecular mutation\n";
  os << " -O pie=<prob>         probability of an inter molecular mutation\n";
  os << " -O trans              translate non organics to something organic\n";
  os << " -O run=<secs>         stop after <secs> seconds\n";
  os << " -O prd=<n>            stop after producing <n> valid molecules\n";
  os << " -O oknp               ok to produce non periodic table elements\n";
  os << " -O organic            reject any non-organic molecule produced\n";

  exit (1);
}

static int
smiles_mutation (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lS:x:p:w:N:n:F:f:c:C:r:R:q:Q:M:z:Y:b:O:a:X:");

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

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (! read_sidechain_library(s))
    {
      cerr << "Cannot read sidechains from '" << s << "'\n";
      return 0;
    }

    if (verbose)
      cerr << "Read " << sidechain_library.number_elements() << "' from '" << s << "'\n";
  }

  if (cl.option_present('N'))
  {
    if (! cl.value('N', number_iterations) || number_iterations < 0)
    {
      cerr << "The number of iterations option (-N) must be a non negative number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will run " << number_iterations << " iterations\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', complete_refresh) || complete_refresh < 0)
    {
      cerr << "The number of iterations for complete refresh (-n) must be a whole +ve number\n";
      usage(4);
    }

    if (0 == number_iterations)
      ;
    else if (complete_refresh >= number_iterations)
    {
      cerr << "When doing " << number_iterations << " iterations (-N), the iterations to refresh (-n) must be less than -N value, " << complete_refresh << " invalid\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will revert to initial smiles every " << complete_refresh << " iterations\n";
  }

  if (cl.option_present('F'))
  {
    if (! cl.value('F', probability_refresh_from_initial) || probability_refresh_from_initial < 0.0 || probability_refresh_from_initial > 1.0)
    {
      cerr << "The probability of refresh from initial option (-F) must be a valid probability\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will refresh from initial smiles with probability " << probability_refresh_from_initial << '\n';
  }

  if (cl.option_present('f'))
  {
    if (! cl.value('f', probability_refresh_from_most_recent_valid) || probability_refresh_from_most_recent_valid < 0.0 || probability_refresh_from_most_recent_valid > 1.0)
    {
      cerr << "The probability of refresh from most recent valid smiles option (-f) must be a valid probability\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will refresh from most recent valid smiles with probability " << probability_refresh_from_most_recent_valid << '\n';
  }

  if (cl.option_present('Y'))
  {
    name_stem = cl.string_value('Y');

    if (verbose)
      cerr << "Molecules created will have names starting with '" << name_stem << "'\n";
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', max_chars_to_remove) || max_chars_to_remove < 0)
    {
      cerr << "The max chars to remove option (-x) must be a non negative number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will remove as many as " << max_chars_to_remove << " contiguous characters during excision operations\n";
  }

  if (cl.option_present('p'))
  {
    if (! cl.value('p', random_replicates) || random_replicates < 1)
    {
      cerr << "The number of random replicates option (-p) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will generate as many as " << random_replicates << " random replicates of each molecule\n";
  }

  if (cl.option_present('w'))
  {
    if (! cl.value('w', number_times_single_character_swap) || number_times_single_character_swap < 1)
    {
      cerr << "The number of single character swaps option (-w) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will swap " << number_times_single_character_swap << " single characters\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', lower_atom_count_cutoff) || lower_atom_count_cutoff < 1)
    {
      cerr << "The lower atom count option (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules with fewer than " << lower_atom_count_cutoff << " atoms\n";
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', upper_atom_count_cutoff) || upper_atom_count_cutoff < 1)
    {
      cerr << "The upper atom count option (-w) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules with more than " << upper_atom_count_cutoff << " atoms\n";
  
    if (upper_atom_count_cutoff < lower_atom_count_cutoff)
    {
      cerr << "upper_atom_count_cutoff " << upper_atom_count_cutoff << " less than lower_atom_count_cutoff " << lower_atom_count_cutoff << endl;
      return 8;
    }
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', min_nrings) || min_nrings < 0)
    {
      cerr << "The lower ring count option (-r) must be a whole non negative number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules with fewer than " << min_nrings << " rings\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', max_rings) || max_rings < 0)
    {
      cerr << "The max ring option (-R) must be a whole non negative number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will discard molecules with more than " << max_rings << " rings\n";
  
    if (max_rings < min_nrings)
    {
      cerr << "max_rings " << max_rings << " less than min_nrings " << min_nrings << endl;
      return 8;
    }
  }

  if (cl.option_present('z'))
  {
    if (! cl.value('z', largest_ring_size_allowed) || largest_ring_size_allowed < 3)
    {
      cerr << "The largest ring size option (-z) must be a valid ring size\n";
      usage(3);
    }

    if (verbose)
      cerr << "Molecules with rings of size above " << largest_ring_size_allowed << " will be discarded\n";
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries_to_match, verbose > 0, 'q'))
    {
      cerr << "Cannot read queries to match (-q) option\n";
      return 4;
    }

    if (verbose)
      cerr << "Read " << queries_to_match.number_elements() << " queries to match\n";

    for (int i = 0; i < queries_to_match.number_elements(); i++)
    {
      queries_to_match[i]->set_max_matches_to_find(1);
    }
  }

  if (cl.option_present('Q'))
  {
    if (! process_queries(cl, queries_to_avoid, verbose > 0, 'Q'))
    {
      cerr << "Cannot read queries to avoid (-Q) option\n";
      return 4;
    }

    if (verbose)
      cerr << "Read " << queries_to_avoid.number_elements() << " queries to avoid\n";

    for (int i = 0; i < queries_to_avoid.number_elements(); i++)
    {
      queries_to_avoid[i]->set_max_matches_to_find(1);
    }
  }

  if (cl.option_present('b'))
  {
    if (! report_progress.initialise(cl, 'b', verbose))
    {
      cerr << "The report every option (-b) must be a whole +ve number\n";
      usage(3);
    }
  }

  if (cl.option_present('O'))
  {
    int i = 0;
    const_IWSubstring o;
    while (cl.value('O', o, i++))
    {
      if (o.starts_with("pie="))
      {
        o.remove_leading_chars(4);

        if (! o.numeric_value(probability_inter_molecule_change) || probability_inter_molecule_change < 0.0 || probability_inter_molecule_change > 1.0)
        {
          cerr << "Invalid probability of inter molecular change 'pie=" << o << "'\n";
          display_dash_o_options(cerr);
        }

        if (verbose)
          cerr << "Probability of inter molecule change " << probability_inter_molecule_change << endl;
      }
      else if (o.starts_with("pia="))
      {
        o.remove_leading_chars(4);

        if (! o.numeric_value(probability_intra_molecule_change) || probability_intra_molecule_change < 0.0 || probability_intra_molecule_change > 1.0)
        {
          cerr << "Invalid probability of intra molecular change 'pia=" << o << "'\n";
          display_dash_o_options(cerr);
        }

        if (verbose)
          cerr << "Probability of intra molecular change " << probability_intra_molecule_change << endl;
      }
      else if ("trans" == o)
      {
        translate_non_organics = 1;

        if (verbose)
          cerr << "Will translate non-organic elements to organic types\n";
      }
      else if (o.starts_with("run="))
      {
        o.remove_leading_chars(4);

        if (! o.numeric_value(run_for) || run_for < 1)
        {
          cerr << "The '-O run=nnn' option must specify a valid +ve number of seconds\n";
          return 3;
        }

        if (verbose)
          cerr << "Will run for " << run_for << " seconds\n";

        time(&tzero);
        number_iterations = std::numeric_limits<int>::max();
      }
      else if (o.starts_with("prd="))
      {
        o.remove_leading_chars(4);

        if (! o.numeric_value(stop_after_producing) || stop_after_producing < 1)
        {
          cerr << "Will stop after producing " << stop_after_producing << " valid molecules\n";
          return 3;
        }

        if (verbose)
          cerr << "Will produce only " << stop_after_producing << " valid molecules\n";

        number_iterations = std::numeric_limits<int>::max();
      }
      else if ("oknp" == o)
      {
        reject_for_non_periodic_table_elements = 0;

        if (verbose)
          cerr << "Will allow production of non periodic table elements\n";
      }
      else if ("organic" == o)
      {
        discard_non_organic_molecules = 1;

        if (verbose)
          cerr << "Non organic molecules will be discarded\n";
      }
      else if (o.starts_with("prxn="))
      {
        o.remove_leading_chars(5);
        if (! o.numeric_value(probability_do_reaction) || probability_do_reaction < 0.0 || probability_do_reaction > 1.0)
        {
          cerr << "INvalid reaction probability 'prxn=" << o << "'\n";
          return 1;
        }
      }
      else if ("help" == o)
      {
        display_dash_o_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -O qualifier '" << o << "'\n";
        display_dash_o_options(cerr);
      }
    }
  }

  if (number_iterations)
    ;
  else if (run_for > 0)
    ;
  else
    number_iterations = 1;

  if (! cl.option_present('a'))
    ;
  else if (! number_assigner.initialise(cl, 'a', verbose))
  {
    cerr << "Cannot initialise number assigner (-a)\n";
    return 3;
  }

  if (cl.option_present('X'))
  {
    if (! rxn.build(cl, 'X', verbose))
    {
      cerr << "Cannot initialise reactions (-X)\n";
      return 2;
    }

    if (0.0 == probability_do_reaction)
      probability_do_reaction = 0.05;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  iw_random_seed();

  set_include_aromaticity_in_smiles(1);   // turned off below. Use to populate the pool
  set_input_aromatic_structures(1);
  set_add_same_bond_twice_fatal(0);
  set_display_abnormal_valence_messages(0);
  set_warn_aromatic_chain_atoms(0);
  set_display_no_kekule_form_message(0);
  set_display_already_bonded_error_message(0);
  set_display_messages_about_unable_to_compute_implicit_hydgogens(0);
  set_display_smiles_interpretation_error_messages(0);
  set_display_strange_chemistry_messages(0);
  set_display_unusual_hcount_warning_messages(0);
  set_include_implicit_hydrogens_on_aromatic_n_and_p(0);

  resizable_array_p<IWString> smiles;

  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! fill_smiles_array(cl[i], input_type, smiles))
    {
      cerr << "Cannot read molecules from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  if (verbose)
  {
    cerr << "Read " << smiles.number_elements() << " molecules\n";
  }

  if (cl.option_present('M'))
  {
    const char * m = cl.option_value('M');
    if (! stream_for_all_valid_smiles.open(m))
    {
      cerr << "Cannot open stream for all valid smiles '" << m << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "All valid smiles written to '" << m << "'\n";
  }

  set_include_aromaticity_in_smiles(0);

  IWString_and_File_Descriptor output(1);

  smiles_mutation(smiles, output);

  output.flush();

  if (verbose)
  {
    cerr << "produced " << smiles_produced << " smiles\n";
    cerr << molecular_interpretation_attempted << " attempts at molecular interpretation\n";
    cerr << valid_molecules_produced << " valid molecules produced\n";
    if (lower_atom_count_cutoff > 0)
      cerr << rejected_for_too_small << " rejected for fewer than " << lower_atom_count_cutoff << " atoms\n";
    if (std::numeric_limits<int>::max() != upper_atom_count_cutoff)
      cerr << rejected_for_too_large << " rejected for more than " << upper_atom_count_cutoff << " atoms\n";
    if (min_nrings > 0)
      cerr << rejected_for_too_few_rings << " rejected for fewer than " << min_nrings << " rings\n";
    if (std::numeric_limits<int>::max() != max_rings)
      cerr << rejected_for_too_many_rings << " rejected for more than " << max_rings << " rings\n";
    cerr << invalid_valence_produced << " molecules with invalid valences produced\n";
    cerr << rejected_for_recreating_initial_structure << " rejected for re-creating initial structures\n";
    cerr << non_organic_molecules_discarded << " non organic molecules discarded\n";
    cerr << non_periodic_table_molecules_rejected << " non periodic table molecules discarded\n";
    cerr << smiles_written << " smiles written\n";

    if (acc_natoms.n() > 0)
    {
      cerr << "Molecules had between " << acc_natoms.minval() << " and " << acc_natoms.maxval() << " atoms, ave " << static_cast<float>(acc_natoms.average()) << '\n';
      cerr << "Molecules had between " << acc_nrings.minval() << " and " << acc_nrings.maxval() << " rings, ave " << static_cast<float>(acc_nrings.average()) << '\n';
    }

    for (int i = 0; i < valid_molecules_produced_per_generation.number_elements(); i++)
    {
      if (valid_molecules_produced_per_generation[i])
        cerr << valid_molecules_produced_per_generation[i] << " molecules produced at generation " << i << '\n';
    }

    if (queries_to_avoid.number_elements())
      cerr << rejected_for_matching_avoid_query << " molecules rejected for matching an avoid query\n";

    if (queries_to_match.number_elements())
      cerr << rejected_for_not_matching_needed_query << " molecules rejected for not matching a needed query\n";
  }

  return 0;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = smiles_mutation (argc, argv);

  return rc;
}
