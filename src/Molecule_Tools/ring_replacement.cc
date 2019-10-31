/*
  companion programme to ring_extraction.
*/

#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "cmdline.h"
#include "misc.h"
#include "accumulator.h"

#define ISTREAM_AND_TYPE_IMPLEMENTATION
#include "istream_and_type.h"
#include "molecule.h"
#include "aromatic.h"
#include "output.h"
#include "iwstandard.h"
#include "path.h"
#include "substructure.h"
#include "target.h"
#include "reaction_duplicate.h"
#include "toggle_kekule_form.h"
#include "numass.h"
#include "smiles.h"

#include "ring_ext_rep.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int molecules_written = 0;

static int molecules_getting_replacements = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static Number_Assigner number_assigner;

static Ring_Extraction_Replacement_Conditions rerc;

static resizable_array_p<Substructure_Query> queries;
static resizable_array_p<Substructure_Query> products_must_have, products_must_not_have;

static int molecules_discarded_by_must_have_queries = 0;
static int molecules_discarded_by_must_not_have_queries = 0;

static int molecules_with_invalid_valences_suppressed = 0;

static int molecules_saved_by_finding_different_kekule_forms = 0;

static Accumulator_Int<int> acc_variants;

static Reaction_Duplicate reaction_duplicate;

static int unique_molecules_only = 0;
static int unique_molecules_names = 0;

static int preserve_aromaticity = 1;

static int molecules_rejected_for_aromaticity_changes = 0;

static Molecule_Output_Object stream_for_molecules_not_transformed;
static Molecule_Output_Object stream_for_molecules_with_bad_valences;

static int write_parent_molecule = 0;

static int min_examples_needed = 0;

static int replacement_rings_discarded_for_count = 0;

static IWString name_token_separator(' ');

static resizable_array_p<Substructure_Query> discard_library_rings;

static int library_rings_discarded_for_discard_library_rings = 0;

static int try_extra_substituents = 0;

static extending_resizable_array<int> molecules_changed_with_extra_substituents;

#define MAX_RING_SIZE 10

class Set_of_Ring_Sizes
{
  private:
    int _nr;
    int _aliphatic_ring_sizes[MAX_RING_SIZE];
    int _aromatic_ring_sizes[MAX_RING_SIZE];
    resizable_array<int> _fused_sizes;     // size1 * MAX_RING_SIZE + size2. Just ring systems with 2 rings

//  private functions

    void _fill_fused_system_array (Molecule & m, resizable_array<int> & fused_sizes);

  public:
    Set_of_Ring_Sizes ();

    int initialise(Molecule & m);

    int rings_match(Molecule & m) const;
};

Set_of_Ring_Sizes::Set_of_Ring_Sizes ()
{
  _nr = 0;

  std::fill_n(_aliphatic_ring_sizes, MAX_RING_SIZE, 0);
  std::fill_n(_aromatic_ring_sizes, MAX_RING_SIZE, 0);

  return;
}

/*
  Look at all fused systems that are of size 2 rings.
  Generate a unique number that describes that combination and place in FUSED_SIZES
*/

static void
fill_fused_system_array (Molecule & m,
                         resizable_array<int> & fused_sizes)
{
  const int nr = m.nrings();

  int * already_done = new_int(nr); std::unique_ptr<int[]> free_already_done(already_done);

  for (auto i = 0; i < nr; ++i)
  {
    if (already_done[i])
      continue;

    const auto ri = m.ringi(i);
    const Ring * second_ring = nullptr;

    if (0 == ri->fused_ring_neighbours())
      continue;

    int fss = 1;

    for (auto j = i + 1; j < nr; ++j)
    {
      if (already_done[j])
        continue;

      const auto rj = m.ringi(j);

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        continue;

      already_done[j] = 1;
      second_ring = rj;
      fss++;
    }

    if (2 == fss)     // we only deal with these
    {
      const Ring * rx = ri;
      const Ring * ry = second_ring;

      if (rx->number_elements() > ry->number_elements())
        std::swap(rx, ry);

      int q1 = rx->number_elements();
      int q2 = ry->number_elements();

      if (rx->is_aromatic())
        q1 = q1 + 20;
      if (ry->is_aromatic())
        q2 = q2 + 17;

      fused_sizes.add(40 * q1 + q2);
    }
  }

  return;
}


int
Set_of_Ring_Sizes::initialise(Molecule & m)
{
  m.compute_aromaticity_if_needed();

  _nr = m.nrings();

  for (int i = 0; i < _nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    const int s = ri->number_elements();

    if (s >= MAX_RING_SIZE)
      continue;

    if (ri->is_aromatic())
      _aromatic_ring_sizes[s]++;
    else
      _aliphatic_ring_sizes[s]++;
  }

  if (_nr < 2)
    return _nr;

  fill_fused_system_array(m, _fused_sizes);

  return _nr;
}

int
Set_of_Ring_Sizes::rings_match(Molecule & m) const
{
  const int nr = m.nrings();

#ifdef DEBUG_RINGS_MATCH
  cerr << "Checking molecule with " << nr << " rings vs " << _nr << " rings available\n";
  for (int i = 0; i < MAX_RING_SIZE; ++i)
  {
    if (_aliphatic_ring_sizes[i] || _aromatic_ring_sizes[i])
      cerr << " ring size " << i << " target has " << _aliphatic_ring_sizes[i] << " aliphatic and " << _aromatic_ring_sizes[i] << " aromatic\n";
  }
#endif

  if (nr < _nr)
    return 0;

  int candidate_rings = nr;    // we will reduce this as we eliminate possibilities

  m.compute_aromaticity_if_needed();

  int aliphatic_rings[MAX_RING_SIZE], aromatic_rings[MAX_RING_SIZE];
  std::fill_n(aliphatic_rings, MAX_RING_SIZE, 0);
  std::fill_n(aromatic_rings,  MAX_RING_SIZE, 0);

  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);

    const int s = ri->number_elements();

    if (s >= MAX_RING_SIZE)
    {
      candidate_rings--;
      if (candidate_rings < _nr)
        return 0;
      continue;
    }

    if ((ri->is_aromatic() && 0 == _aromatic_ring_sizes[s]) ||
        (! ri->is_aromatic() && 0 == _aliphatic_ring_sizes[s]))
    {
      candidate_rings--;
      if (candidate_rings < _nr)
        return 0;
    }

    if (ri->is_aromatic())
      aromatic_rings[s]++;
    else
      aliphatic_rings[s]++;
  }

  for (auto i = 0; i < MAX_RING_SIZE; ++i)
  {
    if (aromatic_rings[i] < _aromatic_ring_sizes[i])
      return 0;

    if (aliphatic_rings[i] < _aliphatic_ring_sizes[i])
      return 0;
  }

  if (0 == _fused_sizes.number_elements())
    return 1;

#ifdef DEBUG_RINGS_MATCH
  cerr << "At line " << __LINE__ << " still checking\n";
#endif

  resizable_array<int> fused_sizes;
  fill_fused_system_array(m, fused_sizes);

  for (unsigned int i = 0; i < _fused_sizes.size(); ++i)
  {
    if (! fused_sizes.contains(_fused_sizes[i]))
      return 0;
  }

  return 1;
}

static bool
matches_any_of_these (Molecule_to_Match & target,
                      resizable_array_p<Substructure_Query> & q)
{
  const auto n = q.size();

  for (unsigned int i = 0; i < n; ++i)
  {
    if (q[i]->substructure_search(target))
      return true;
  }

  return false;
}

/*
  We have a collection of information about the possible replacement rings

  we try to impose a clockwise ordering of the atoms. First the (optional)
  ring fusion points. Then any number of substitution points. The important
  thing is to preserve a clock-wise ordering
*/

class Replacement_Ring : public Molecule
{
  private:
    IWString _id;

    Set_of_Atoms _ring_fusion_points;

    Set_of_Atoms _substitution_points;

    Substructure_Query _query;

    Set_of_Ring_Sizes _sors;

    int _count;

//  private functions

    int _process (const Molecule & m, const IWString & initial_usmi, const Set_of_Atoms & e, std::ostream & output, const int current_cnt) const;
  public:
    Replacement_Ring();

    int build (const const_IWSubstring &, const Ring_Extraction_Replacement_Conditions &);

    int number_ring_fusion_points() const { return _ring_fusion_points.number_elements();}
    int number_substitution_points() const { return _substitution_points.number_elements();}

    int count () const { return _count;}

    int process (Molecule & m, const int * process_these_atoms, std::ostream & output);
};

static resizable_array_p<Replacement_Ring> replacement_rings;

template class resizable_array_p<Replacement_Ring>;
template class resizable_array_base<Replacement_Ring *>;
template class data_source_and_type<Replacement_Ring>;

Replacement_Ring::Replacement_Ring()
{
  _count = 0;

  return;
}

int
Replacement_Ring::build (const const_IWSubstring & buffer,
                         const Ring_Extraction_Replacement_Conditions & rerc)
{
  if (buffer.nwords() < 3)
  {
    cerr << "Replacement_Ring::build:ring specification must have at least <smarts> <smiles> <id>\n";
    return 0;
  }

  const_IWSubstring smiles;
  int i = 0;
  buffer.nextword(smiles, i);

  const_IWSubstring smarts;
  buffer.nextword(smarts, i);

  buffer.nextword(_id, i);

  if (! _query.create_from_smarts(smarts))
  {
    cerr << "Replacement_Ring::build:invalid smarts '" << smarts << "'\n";
    return 0;
  }

  _query.set_do_not_perceive_symmetry_equivalent_matches(1);
  
  if (! this->build_from_smiles(smiles))
  {
    cerr << "Replacement_Ring::build:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  const int matoms = natoms();

  assert (matoms > 0);

  for (int i = 0; i < matoms; ++i)
  {
    const int iso = isotope(i);
    if (0 == iso)
      continue;

    if (hcount(i))
      continue;


    unset_all_implicit_hydrogen_information(i);
    recompute_implicit_hydrogens(i);
  }

  if (0 == nrings())
  {
    cerr << "Replacement_Ring::build:no rings '" << _id << "'\n";
    return 0;
  }

#ifdef CARE_ABOUT_ISOTOPES

  assert (matoms > 0);

  for (int i = 0; i < matoms; i++)
  {
    int iso = isotope(i);

    if (0 == iso)
      continue;

    if (rerc.isotope_for_substitution_points() == iso)
      _substitution_points.add(i);
    else if (rerc.isotope_for_ring_fusion() == iso)
      _ring_fusion_points.add(i);
    else
    {
      cerr << "Replacement_Ring::initialise:unrecognised isotope " << iso << ", in '" << _id << "'\n";
    }
  }
#endif

  const_IWSubstring token;
  if (buffer.nextword(token, i))
  {
    if (! token.numeric_value(_count) || _count < 1)
    {
      cerr << "Replacement_Ring::build:invalid count '" << buffer << "'\n";
      return 0;
    }
  }

  _sors.initialise(*this);

  return 1;
}

int
Replacement_Ring::process (Molecule & m,
                           const int * process_these_atoms,
                           std::ostream & output)
{
  if (! _sors.rings_match(m))
    return 0;

  Molecule_to_Match target(&m);
  Substructure_Results sresults;

  int nhits = _query.substructure_search(m, sresults);

  if (0 == nhits)
    return 0;

  const IWString initial_usmi = m.unique_smiles();

  for (int i = 0; i < nhits; i++)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    if (! e->any_members_set_in_array(process_these_atoms))
      continue;

    _process(m, initial_usmi, *e, output, i);
  }

  return 1;
}


static int
can_find_different_kekule_form (Molecule & m,
                                int * aromatic_atom)
{
  set_display_no_kekule_form_message(0);

  int rc = m.find_kekule_form(aromatic_atom);

  set_display_no_kekule_form_message(1);

  if (0 == rc)         // failed to find Kekule form
  {
    molecules_with_invalid_valences_suppressed++;

    if (stream_for_molecules_with_bad_valences.active())
      stream_for_molecules_with_bad_valences.write(m);

    return 0;
  }

  if (m.valence_ok())    // found Kekule form and valence is OK
  {
    molecules_saved_by_finding_different_kekule_forms++;
    return 1;
  }

  molecules_with_invalid_valences_suppressed++;
  return 0;
}

int
Replacement_Ring::_process (const Molecule & initial_molecule,
                            const IWString & initial_usmi,
                            const Set_of_Atoms & e,
                            std::ostream & output,
                            const int current_cnt) const
{
  Molecule m(initial_molecule);
  m.set_name(initial_molecule.name());

  set_display_abnormal_valence_messages(0);

  int initial_aromaticity = -1;

  if (preserve_aromaticity)
    initial_aromaticity = m.aromatic_atom_count();

  int initial_matoms = m.natoms();

  m.add_molecule(this);

  int * aromatic_atom = new_int(m.natoms()); std::unique_ptr<int[]> free_aromatic_atom(aromatic_atom);

  m.aromaticity(aromatic_atom);

// We rely on the order of the matched atoms being the same as the order in the replacement piece

  int extra_atoms = natoms();

//  Because we don't know how the Molecule will change the Atom, we need to do
//  the changes in two passes. First to identify what changes, then to do the changes

  resizable_array_p<Bond> bonds_to_be_removed, bonds_to_be_added;

  for (int i = 0; i < extra_atoms; i++)
  {
    atom_number_t old_atom = e[i];

    const Atom * a = m.atomi(old_atom);

    int acon = a->ncon();

    m.set_isotope(initial_matoms + i, 0);
    m.set_implicit_hydrogens_known(initial_matoms + i, 0);

    for (int j = 0; j < acon; j++)
    {
      const Bond * b = a->item(j);

      atom_number_t k = b->other(old_atom);

      if (e.contains(k))
        continue;

      if (m.are_bonded(old_atom, k))
      {
        Bond * b = new Bond(old_atom, k, SINGLE_BOND);
        bonds_to_be_removed.add(b);
      }

      Bond * newbond = new Bond(k, initial_matoms + i, b->btype());
      bonds_to_be_added.add(newbond);
    }
  }

// Now do any bond breakages and changes

  int n = bonds_to_be_removed.number_elements();
  for (int i = 0; i < n; i++)
  {
    const Bond * b = bonds_to_be_removed[i];
    m.remove_bond_between_atoms(b->a1(), b->a2());
  }

  n = bonds_to_be_added.number_elements();
  for (int i = 0; i < n; i++)
  {
    const Bond * b = bonds_to_be_added[i];
    m.add_bond(b->a1(), b->a2(), b->btype());
  }

// Restore the original atom ordering

  for (int i = 0; i < extra_atoms; i++)
  {
    m.swap_atoms(initial_matoms + i, e[i]);
  }

  m.resize(initial_matoms);

  if (m.valence_ok())
    ;
  else if (! can_find_different_kekule_form(m, aromatic_atom))
  {
    set_display_abnormal_valence_messages(1);
    return 0;
  }

  set_display_abnormal_valence_messages(1);

  if (preserve_aromaticity && initial_aromaticity != m.aromatic_atom_count())
  {
    molecules_rejected_for_aromaticity_changes++;
    return 0;
  }

  if (initial_usmi == m.unique_smiles())
    return 0;

  if (unique_molecules_only && reaction_duplicate.is_duplicate(m))
    return 0;

  if (products_must_not_have.size() > 0 || products_must_have.size() > 0)
  {
    Molecule_to_Match target(&m);
    if (products_must_not_have.size() > 0 && matches_any_of_these(target, products_must_not_have))
    {
      molecules_discarded_by_must_not_have_queries++;
      return 0;
    }

    if (products_must_have.size() > 0 && ! matches_any_of_these(target, products_must_have))
    {
      molecules_discarded_by_must_have_queries++;
      return 0;
    }
  }

  if (number_assigner.active())
    number_assigner.process(m);

  if (unique_molecules_names) {
    output << m.smiles() << ' ' << m.name() << "_" << current_cnt << name_token_separator << "%%" << name_token_separator << _id << name_token_separator << _count << '\n';
  } 
  else {
    output << m.smiles() << ' ' << m.name() << name_token_separator << "%%" << name_token_separator << _id << name_token_separator << _count << '\n';
  }
  molecules_written++;

  return 1;
}

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Companion programme to ring_extraction\n";
//display_standard_ring_ext_rep_options (cerr);
  cerr << "  -P <fname>    file of labelled rings created by ring_extraction\n";
  cerr << "  -s <smarts>   only replace rings matched by <smarts>\n";
  cerr << "  -q <query>    only replace rings matched by <query>\n";
  cerr << "  -u            unique molecules only\n";
  cerr << "  -a            allow change of aromaticity on ring replacement\n";
  cerr << "  -B <fname>    write molecules not transformed to <fname>\n";
  cerr << "  -p            write parent molecule\n";
  cerr << "  -n <ex>       only process replacement rings with <ex> or more examples\n";
  cerr << "  -M <query>    product molecules MUST     match a query in <query>\n";
  cerr << "  -X <query>    product molecules must NOT match any queries in <query>\n";
  cerr << "  -D <query>    discard any replacement ring that matches <query>\n";
  cerr << "  -e ...        standard number assigner options\n";
  cerr << "  -x 'char'     name token separator (default ' ')\n";
  cerr << "  -m <n>        add <n> methyl groups to ring atoms to find more substituted replacements\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static int
read_ring_to_be_substituted (const const_IWSubstring & buffer,
                             resizable_array_p<Replacement_Ring> & rings)
{
  Replacement_Ring * r = new Replacement_Ring;

  if (! r->build(buffer, rerc))
  {
    delete r;
    return 0;
  }

  if (r->count() < min_examples_needed)
  {
    replacement_rings_discarded_for_count++;
    delete r;
    return 1;
  }

  if (discard_library_rings.size())
  {
    Molecule_to_Match target(r);
    if (matches_any_of_these(target, discard_library_rings))
    {
      library_rings_discarded_for_discard_library_rings++;
      delete r;
      return 1;
    }
  }

  rings.add(r);

  return 1;
}

static int
read_rings_to_be_substituted(iwstring_data_source & input,
                             resizable_array_p<Replacement_Ring> & rings)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (! read_ring_to_be_substituted(buffer, rings))
    {
      cerr << "Invalid replacement ring specification, line " << input.lines_read() << endl;
      cerr << buffer << endl;
      return 0;
    }
  }

  return 1;
}

static int
read_rings_to_be_substituted(const const_IWSubstring & fname,
                             resizable_array_p<Replacement_Ring> & rings)
{
  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_rings_to_be_substituted(input, rings);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (rerc.remove_chirality())
    m.remove_all_chiral_centres();

  return;
}

static int
identify_rings_to_be_processed (Molecule & m,
                                int * process_these_atoms)
{
  Molecule_to_Match target(&m);

  int nq = queries.number_elements();

  int rc = 0;

  for (int i = 0; i< nq; i++)
  {
    Substructure_Results sresults;

    int nhits = queries[i]->substructure_search(target, sresults);

    if (0 == nhits)
      continue;

    rc++;

    sresults.each_embedding_set_vector(process_these_atoms, 1);
  }

  return rc;
}

static int
do_try_extra_substituents(const int level,
                          Molecule & m,
                          const int * process_these_atoms,
                          std::ostream & output)
{
  const int n = replacement_rings.number_elements();

  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == process_these_atoms[i])
      continue;

    if (! m.is_ring_atom(i))
      continue;

    if (0 == m.hcount(i))
      continue;

    Atom * c = new Atom(6);
    m.add(c);
    m.add_bond(i, matoms, SINGLE_BOND);

//  cerr << "Checking " << m.smiles() << " level " << level << endl;

    for (int j = 0; j < n; ++j)
    {
      if (replacement_rings[j]->process(m, process_these_atoms, output))
        rc++;
    }

    if (level > 0)
      rc += do_try_extra_substituents(level - 1, m, process_these_atoms, output);

    m.resize(matoms);
  }

  if (rc)
    molecules_changed_with_extra_substituents[level]++;

  return rc;
}

/*
  Here's where you do whatever you want to do with the molecule
  In this case, we count the number of nitrogen atoms
*/

static int
ring_replacement(Molecule & m,
                  std::ostream & output)
{
  const int nr = m.nrings();

  if (0 == nr)
    return 1;

  if (! m.valence_ok())
  {
    cerr << "ring_replacement:cannot process molecules with invalid valences '" << m.name() << "'\n";
    return 1;
  }

  (void) reaction_duplicate.is_duplicate(m);

  if (write_parent_molecule)
    output << m.smiles() << ' ' << m.name() << '\n';

  const int matoms = m.natoms();

  int * process_these_atoms = new int[matoms + matoms]; std::unique_ptr<int[]> free_process_these_atoms(process_these_atoms);

  if (queries.number_elements())
  {
    std::fill_n(process_these_atoms, matoms + matoms, 0);
    if (! identify_rings_to_be_processed(m, process_these_atoms))
      return 1;
  }
  else
    std::fill_n(process_these_atoms, matoms + matoms, 1);

  const int n = replacement_rings.number_elements();

  int this_molecule_changed = 0;

  for (int i = 0; i < n; i++)
  {
    if (replacement_rings[i]->process(m, process_these_atoms, output))
      this_molecule_changed++;
  }

  if (try_extra_substituents)
    this_molecule_changed += do_try_extra_substituents(try_extra_substituents, m, process_these_atoms, output);

  if (this_molecule_changed > 0)
    acc_variants.extra(this_molecule_changed);

  if (verbose > 1)
    cerr << "Input molecule " << m.name() << " generated " << this_molecule_changed << " variants\n";

  if (this_molecule_changed)
    molecules_getting_replacements++;
  else if (stream_for_molecules_not_transformed.active())
    stream_for_molecules_not_transformed.write(m);

  return 1;
}

static int
ring_replacement (data_source_and_type<Molecule> & input,
                     std::ostream & output)
{
  Molecule * m;
  while (NULL != (m = input.next_molecule()))
  {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (! ring_replacement(*m, output))
      return 0;
  }

  return output.good();
}

static int
ring_replacement (const char * fname, int input_type, std::ostream & output)
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

  return ring_replacement(input, output);
}

static int
ring_replacement (int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:i:g:laP:s:q:uQB:pV:n:e:x:M:X:D:m:");

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
  else
    set_global_aromaticity_type(Daylight);

  set_unset_implicit_hydrogens_known_if_possible(1);

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

  if (! rerc.initialise(cl, verbose))
  {
    cerr << "Cannot initialise remove and extract option(s)\n";
    return 4;
  }

  int input_type = 0;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.option_present('n'))
  {
    if (! cl.value('n', min_examples_needed) || min_examples_needed < 1)
    {
      cerr << "The minimum number of examples needed (-n) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will not consider replacement rings that have fewer than " << min_examples_needed << " examples\n";
  }

  if (cl.option_present('e'))
  {
    if (! number_assigner.initialise(cl, 'e', verbose))
    {
      cerr << "Cannot initialise number assigner (-e)\n";
      return 2;
    }
  }

  if (cl.option_present('x'))
  {
    cl.value('x', name_token_separator);

    if (verbose)
      cerr << "Name token separator set to '" << name_token_separator << "'\n";
  }

  if (cl.option_present('D'))
  {
    if (! process_queries(cl, discard_library_rings, verbose, 'D'))
    {
      cerr << "Cannot read discard library rings queries (-D)\n";
      return 2;
    }

    if (verbose)
      cerr << "Read " << discard_library_rings.size() << " discard library rings (-D) queries\n";
  }

  if (! cl.option_present('P'))
  {
    cerr << "Must specify (-P) file of mared rings generated by ring_extraction\n";
    usage(4);
  }

  if (cl.option_present('P'))
  {
    const_IWSubstring p;
    for (auto i = 0; cl.value('P', p, i); ++i)
    {
      if (! read_rings_to_be_substituted(p, replacement_rings))
      {
        cerr << "Cannot read set of replacement rings from -P file '" << p << "'\n";
        return 8;
      }
    }

    if (0 == replacement_rings.number_elements())
    {
      cerr << "No replacement ring specifications read\n";
      return 1;
    }

    if (verbose)
      cerr << "Read " << replacement_rings.number_elements() << " replacement rings\n";

    if (verbose && min_examples_needed > 0)
      cerr << "Discarded " <<  replacement_rings_discarded_for_count << " rings for fewer than " << min_examples_needed << " examples\n";

    if (verbose && discard_library_rings.size())
      cerr << "Discarded " << library_rings_discarded_for_discard_library_rings << " rings that matched a -D query\n";
  }

  if (cl.option_present('s'))
  {
    int i = 0;
    const_IWSubstring smarts;
    while (cl.value('s', smarts, i++))
    {
      Substructure_Query * q = new Substructure_Query;
      if (! q->create_from_smarts(smarts))
      {
        delete q;
        return 60 + i;
      }

      if (verbose)
        cerr << "Created query from smarts '" << smarts << "'\n";

      queries.add(q);
    }
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries, verbose, 'q'))
    {
      cerr << "Cannot read queries (-q)\n";
      return 4;
    }
  }

  if (cl.option_present('M'))
  {
    if (! process_queries(cl, products_must_have, verbose, 'M'))
    {
      cerr << "Cannot read products must have queries (-M)\n";
      return 2;
    }
  }

  if (cl.option_present('X'))
  {
    if (! process_queries(cl, products_must_not_have, verbose, 'X'))
    {
      cerr << "Cannot read products must not have queries (-X)\n";
      return 2;
    }
  }

  if (verbose && queries.number_elements() > 0)
    cerr << "Read " << queries.number_elements() << " queries for rings to process\n";

  if (cl.option_present('m'))
  {
    if (! cl.value('m', try_extra_substituents) || try_extra_substituents < 1)
    {
      cerr << "The number of extra substituents to try (-m) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will try to add as many as " << try_extra_substituents << " extra substituents\n";
  }

  if (cl.option_present('u'))
  {
    unique_molecules_only = 1;

    if (verbose)
      cerr << "Will suppress duplicate products\n";
  }

  if (cl.option_present('Q'))
  {
    unique_molecules_names = 1;

    if (verbose)
      cerr << "Will generate unique product names\n";
  }

  if (cl.option_present('a'))
  {
    preserve_aromaticity = 0;

    if (verbose)
      cerr << "Aromaticity may be destroyed during a ring replacement\n";
  }

  if (cl.option_present('p'))
  {
    write_parent_molecule = 1;

    if (verbose)
      cerr << "The starting molecule will be written\n";
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('B'))
  {
    const_IWSubstring b = cl.string_value('B');

    stream_for_molecules_not_transformed.add_output_type(SMI);
    if (stream_for_molecules_not_transformed.would_overwrite_input_files(cl, b))
    {
      cerr << "Cannot overwrite input file(s), -B = " << b << "'\n";
      return 3;
    }

    if (! stream_for_molecules_not_transformed.new_stem(b))
    {
      cerr << "Cannot initialise stream for non-chaning molecules '" << b << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Molecules not changed written to '" << b << "'\n";
  }

  if (cl.option_present('V'))
  {
    const_IWSubstring v = cl.string_value('V');

    stream_for_molecules_with_bad_valences.add_output_type(SMI);
    if (stream_for_molecules_with_bad_valences.would_overwrite_input_files(cl, v))
    {
      cerr << "Cannot overwrite input file(s), -V = " << v << "'\n";
      return 3;
    }

    if (! stream_for_molecules_with_bad_valences.new_stem(v))
    {
      cerr << "Cannot initialise stream for bad valence molecules '" << v << "'\n";
      return 4;
    }

    if (verbose)
      cerr << "Molecules with bad valences written to '" << v << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! ring_replacement(cl[i], input_type, std::cout))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << molecules_getting_replacements << " molecules got replacement rings\n";
    if (acc_variants.n() > 0)
      cerr << "Created as few as " << acc_variants.minval() << " and as many as " << acc_variants.maxval() << " variants. Ave " << acc_variants.average() << endl;

    if (try_extra_substituents)
    {
      int last_non_zero = -1;

      for (int i = 0; i < molecules_changed_with_extra_substituents.number_elements(); ++i)
      {
        if (molecules_changed_with_extra_substituents[i] > 0)
          last_non_zero = i;
      }

      for (int i = 0; i <= last_non_zero; ++i)
      {
        cerr << molecules_changed_with_extra_substituents[i] << " molecules found extra replacements with " << i << " extra substituents\n";
      }
    }

    cerr << molecules_rejected_for_aromaticity_changes << " molecules discarded for aromaticity changes\n";
    cerr << molecules_with_invalid_valences_suppressed << " molecules discarded for invalid valences\n";
    if (products_must_have.size())
      cerr << molecules_discarded_by_must_have_queries << " molecules discarded for failing products must have queries\n";
    if (products_must_not_have.size())
      cerr << molecules_discarded_by_must_not_have_queries << " molecules discarded for failing products must NOT have queries\n";
    cerr << molecules_saved_by_finding_different_kekule_forms << " valence errors saved by finding different Kekule forms\n";
    cerr << molecules_written << " molecules written\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = ring_replacement(argc, argv);

  return rc;
}
