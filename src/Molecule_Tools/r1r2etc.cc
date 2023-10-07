/*
  Examines connection points.
*/

#include <stdint.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/ematch.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

/*
  Need to compute the following descriptors

  dtodo = ['natoms', 'nrings', 'amw', 'ro5_ohnh', 'ro5_on', 'nvrtspsa', 'halogen', 'fccsp3', 'ringsys', 'arring', 'alring', 'ringatom', 'mxdst']
*/

/*
  Lifted from iwdescr
*/

static int
compute_rule_of_five_stuff (Molecule & m,
                            int & on, int & ohnh)
{
  int matoms = m.natoms();

  on = 0;
  ohnh = 0;

  for (int i = 0; i < matoms; i++)
  {
    const atomic_number_t z = m.atomic_number(i);

    if (7 == z || 8 == z)
    {
      on++;
      if (m.hcount(i) > 0)
        ohnh++;
    }
  }

  return 1;
}

static int
compute_halogen (const Molecule & m)
{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (m.is_halogen(i))
      rc++;
  }

  return rc;
}

static float
compute_fccsp3 (const Molecule & m)
{
  const int matoms = m.natoms();

  int carbon = 0;
  int csp3 = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    if (6 != a->atomic_number())
      continue;

    carbon++;

    if (a->ncon() == a->nbonds())
      csp3++;
  }

  if (0 == carbon)
    return 0.0f;

  return static_cast<float>(csp3) / static_cast<float>(carbon);
}

static float
compute_fcsp3 (const Molecule & m)
{
  const int matoms = m.natoms();

  int csp3 = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    if (6 != a->atomic_number())
      continue;

    if (a->ncon() == a->nbonds())
      csp3++;
  }

  return static_cast<float>(csp3) / static_cast<float>(matoms);
}


static int
compute_mxdst (Molecule & m)
{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    for (int j = i + 1; j < matoms; ++j)
    {
      const auto b = m.bonds_between(i, j);

      if (b > rc)
        rc = b;
    }
  }

  return rc;
}

static int
compute_ring_related (Molecule & m,
                      int & ring_systems,
                      int & alring,
                      int & arring)
{
  ring_systems = 0;
  alring = 0;
  arring = 0;

  const int nrings = m.nrings();

  if (0 == nrings)
    return 1;

  if (1 == nrings)
  {
    ring_systems = 1;
    const Ring * r0 = m.ringi(0);
    if (r0->is_aromatic())
      arring = 1;
    else
      alring = 1;

    return 1;
  }

  if (2 == nrings)
  {
    const Ring * r0 = m.ringi(0);
    const Ring * r1 = m.ringi(1);

    if (r0->is_aromatic())
      arring++;
    else
      alring++;

    if (r1->is_aromatic())
      arring++;
    else
      alring++;

    if (r0->is_fused())
      ring_systems = 1;
    else
      ring_systems = 2;

    return 1;
  }

  int * ring_already_done = new_int(nrings);std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < nrings; ++i)
  {
    const Ring * ri = m.ringi(i);
    if (ri->is_aromatic())
      arring++;
    else
      alring++;

    if (ring_already_done[i])
      continue;

    if (! ri->is_fused())
      continue;

    for (int j = i + 1; j < nrings; ++j)
    {
      if (ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        continue;

      ring_already_done[j] = 1;
    }

    ring_systems++;
  }

  return 1;
}

static int
compute_ring_atoms (Molecule & m)
{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (m.is_ring_atom(i))
      rc++;
  }

  return rc;
}

static int
ensure_parent_written (const IWString & buffer,
                       int & parent_written,
                       IWString & output)
{
  if (parent_written)
    return 0;

  output << buffer;

  parent_written = 1;

  return 1;
}

/*
  For each query, we need to keep track of which matched atoms
  correspond to the various R atoms. 
  For each R number, we keep the corresponding matched query atoms
*/

class Query_and_R_Atoms : public Substructure_Query
{
  private:
    int * _r;              // arbitrary numbers, 3,1,2
    int * _m;              // corresponding query atom, 18,0,5
    int _nr;

// For efficiency, we pre-allocate an array to hold the atoms and bond types lost as
// the bonds are removed. This makes this class non-thread safe

    atom_number_t * _atom1;
    atom_number_t * _atom2;          // may be empty if a hydrogen substituent
    bond_type_t * _bond_lost;

    int _maximum_substituent_size;

    int _fcsp3;

    int _remove_isotopes_from_scaffold;

    int _translate_sidechain_isotopes_to;

    int _only_substituents_at_matched_atoms;
    
    int _allow_symmetry_related_matches;

    char _output_separator;
    
    int _brief_output;

// private functions

    int _identify_attachment_points (Molecule & m, const Set_of_Atoms & e, const int ndx);
    int _append_smiles_and_molecular_descriptors (Molecule & m, IWString & buffer) const;
    int _append_smiles_and_molecular_descriptors_Hygrogen (IWString & buffer) const;
    int _remove_explicit_hydrogens (Molecule & m, int preserve_isotopic_labelled) const;

  public:
    Query_and_R_Atoms();
    ~Query_and_R_Atoms();

    int number_R_atoms () const { return _nr;}

    void set_output_separator (char s) {_output_separator = s;}
    void set_brief_output (int b) {_brief_output = b;}
    void set_maximum_substituent_size (int s) {_maximum_substituent_size = s;}
    void set_remove_isotopes_from_scaffold (int s) { _remove_isotopes_from_scaffold = s;}
    void set_translate_sidechain_isotopes_to (int s) { _translate_sidechain_isotopes_to = s;}
    void set_only_substituents_at_matched_atoms (int s) { _only_substituents_at_matched_atoms = s;}
    void set_allow_symmetry_related_matches (int s) { _allow_symmetry_related_matches = s;}

    void set_fcsp3 (int s) { _fcsp3 = s;}

    int build (MDL_Molecule & m, Molecule_to_Query_Specifications & mqs, int & r_atoms_in_queries);

    int identify_attachment_points (Molecule & m, const Set_of_Atoms & e) const;

    int process (Molecule & m, const Set_of_Atoms & e, IWString & output);

    int process_for_jibo (Molecule & m, const Set_of_Atoms & e, const int preserve_isotopic_labelled,
                                     IWString_and_File_Descriptor & buffer,
                                     IWString mainName,
                                     IWString mainSmiles);

    int is_duplicate_embedding (const int matoms, const Set_of_Atoms & e, resizable_array<uint64_t> & a) const;
};

Query_and_R_Atoms::Query_and_R_Atoms()
{
  _r = nullptr;
  _m = nullptr;
  _nr = 0;
  _atom1 = nullptr;
  _atom2 = nullptr;
  _bond_lost = nullptr;

  _maximum_substituent_size = std::numeric_limits<int>::max();

  _fcsp3 = 1;
  
  _remove_isotopes_from_scaffold = 0;

  _translate_sidechain_isotopes_to = 0;

  _only_substituents_at_matched_atoms = 1;

  _allow_symmetry_related_matches = 1;

  _output_separator = ' ';
  
  _brief_output = 0;

  return;
}

Query_and_R_Atoms::~Query_and_R_Atoms()
{
  if (nullptr != _r)
    delete [] _r;

  if (nullptr != _m)
    delete [] _m;

  if (nullptr != _bond_lost)
    delete [] _bond_lost;

  if (nullptr != _atom1)
    delete [] _atom1;

  if (nullptr != _atom2)
    delete [] _atom2;

  return;
}

int
Query_and_R_Atoms::build (MDL_Molecule & m,
                          Molecule_to_Query_Specifications & mqs,
                          int & r_atoms_in_queries)
{
  const int matoms = m.natoms();

  auto rntmp = new std::pair<int, int>[matoms]; std::unique_ptr<std::pair<int, int>[]> free_rntmp(rntmp);
  int * atom_is_r = new_int(matoms, 0); std::unique_ptr<int[]> free_atom_is_r(atom_is_r);

  int ndx = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const IWString & s = m.atomic_symbol(i);
    if (2 != s.length())
      continue;

    if ('R' != s[0])
      continue;

    if (1 != m.ncon(i))
    {
      cerr << "R1R2Etc::_read_query:atom " << s << " has " << m.ncon(i) << " connections, impossible\n";
      return 0;
    }

    char s1 = s[1];

    int n = s1 - '0';

    if (n < 1 || n > 9)
      continue;

    atom_is_r[i] = 1;

    rntmp[ndx].first = n;
    rntmp[ndx].second = i;
    ndx++;
  }

  if (0 == ndx)
  {
    cerr << "Query_and_R_Atoms::build:no R atoms in molecule, cannot continue\n";
    return 0;
  }

  r_atoms_in_queries += ndx;    // update global counter

// Set max connectivity of all non R atoms to what they are now
// Remember that we will be searching against molecules with explicit hydrogen atoms

  if (_only_substituents_at_matched_atoms)
  {
    for (int i = 0; i < matoms; ++i)
    {
      if (atom_is_r[i])
        continue;

      if (1 == m.atomic_number(i))
        continue;

      auto a = m.mdl_atom_data(i);

//    In theory this is a good thing. But in the specific case of a query with a 3 connected Nitrogen atom,
//    that would prevent it matching a 4 connected, charged Nitrogen atom, so we suppress this condition.

//    const auto ih = m.implicit_hydrogens(i);
//    if (0 == ih) {
//      a->set_max_ncon(m.ncon(i));
//      continue;
//    }

      a->set_hcount(m.hcount(i) + 1);    // remember, ISIS does strange things with the hcount designator
//    cerr << "Atom " << i << " on " << m.smarts_equivalent_for_atom(i) << " set to " << (m.hcount(i)+1) << endl;
    }
  }
  else
  {
    for (int i = 0; i < matoms; ++i)   // remove this, don't need to do anything
    {
      if (atom_is_r[i])
        continue;

      if (1 == m.atomic_number(i))
        continue;
    }
  }

  Element_Matcher ematch;

  if (! ematch.construct_from_string("RX=^R[0-9]$"))
  {
    cerr << "Yipes, cannot initialise R element regular expression\n";
    return 0;
  }

  if (! m.change_R_groups_to_match_any_atom(ematch, _only_substituents_at_matched_atoms))
  {
    cerr << "Query_and_R_Atoms::build:did not convert any R groups\n";
    return 0;
  }

  if (ndx > 1) {
    std::sort(rntmp, rntmp + ndx, [] (const std::pair<int, int> & r1, const std::pair<int, int> & r2) { return r1.first < r2.first;});
  }

  _nr = ndx;
  _r = new int[_nr];
  _m = new int[_nr];

// At this stage we have a sorted list of R numbers and the corresponding atom number.
// Make sure the R numbers are unique

  for (int i = 0; i < ndx; ++i)
  {
    if (i > 0 && rntmp[i-1].first == rntmp[i].first)
    {
      cerr << "R1R2Etc::_read_query:duplicate R numbers " << rntmp[i].first << endl;
      return 0;
    }

    _r[i] = rntmp[i].first;            // R number
    _m[i] = rntmp[i].second;           // corresponding matched atom
//  cerr << "Query_and_R_Atoms:loaded _r array with " << _r[i] << endl;
  }

  if (! Substructure_Query::create_from_molecule(m, mqs))
  {
    cerr << "R1R2Etc::_read_query:cannot create query\n";
    return 0;
  }

  if (! _allow_symmetry_related_matches)
    this->set_do_not_perceive_symmetry_equivalent_matches(1);

  _bond_lost = new bond_type_t[ndx];
  _atom1 = new atom_number_t[ndx];
  _atom2 = new atom_number_t[ndx];

  return 1;
}

int
Query_and_R_Atoms::_append_smiles_and_molecular_descriptors (Molecule & m,
                                           IWString & buffer) const
{
  const auto matoms = m.natoms();

  if (0 == matoms)
    return _append_smiles_and_molecular_descriptors_Hygrogen(buffer);

  if (_translate_sidechain_isotopes_to)
  {
    for (int i = 0; i < matoms; ++i)
    {
      if (m.isotope(i) > 0)
        m.set_isotope(i, _translate_sidechain_isotopes_to);
    }
  }

  buffer << _output_separator << m.unique_smiles();

  if (!_brief_output)
  {
    buffer << _output_separator << matoms;
    buffer << _output_separator << m.nrings();
    buffer << _output_separator << static_cast<float>(m.molecular_weight_ignore_isotopes());

    int on, ohnh;
    compute_rule_of_five_stuff(m, on, ohnh);

    buffer <<_output_separator << on;
    buffer << _output_separator << ohnh;

  // for PSA attach a Carbon atom to stabilise the calculations

    Molecule mcopy(m);
    for (int i = 0; i < matoms; ++i)
    {
      if (0 == mcopy.isotope(i))
        continue;

      const auto e = get_element_from_atomic_number(6);

      mcopy.add(e);
      mcopy.add_bond(i, matoms, SINGLE_BOND);
      break;
    }

    buffer << _output_separator << static_cast<float>(novartis_polar_surface_area(mcopy));

    buffer << _output_separator << compute_halogen(m);

    if (_fcsp3)
      buffer << _output_separator << compute_fcsp3(m);
    else
      buffer << _output_separator << compute_fccsp3(m);

    int ring_systems, alring, arring;
    compute_ring_related(m, ring_systems, alring, arring);

    buffer << _output_separator << ring_systems;
    buffer << _output_separator << arring;
    buffer << _output_separator << alring;

    buffer << _output_separator << compute_ring_atoms(m);

    buffer << _output_separator << compute_mxdst(m);
  }
  return 1;
}

int
Query_and_R_Atoms::_append_smiles_and_molecular_descriptors_Hygrogen (IWString & buffer) const
{
  buffer << _output_separator << 'H';
  
  if (!_brief_output)
  {
    buffer << _output_separator << '1';           // natoms
    buffer << _output_separator << '0';           // nrings

    const auto h = get_element_from_atomic_number(1);

    buffer << _output_separator << h->atomic_mass();

    IWString sep0;
    sep0 << _output_separator << '0';

    for (int i = 0; i < 10; ++i)
    {
      buffer << sep0;
    }   
  }

  return 1;
}

/*
  Add the appropriate R group to every even numbered isotopic atom
*/

static int
add_r_atoms_to_scaffold (Molecule & m)
{
  m.remove_all(1);    // get rid of explicit hydrogens

  const auto matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    const auto iso = m.isotope(i);
    if (0 == iso)
      continue;

    if (0 != iso % 2)
      continue;

    IWString s;
    s << 'R' << iso / 2;

    const Element * e = get_element_from_symbol_no_case_conversion(s);

    m.add(e);

    m.add_bond(i, matoms + rc, SINGLE_BOND);
    rc++;
  }

  return rc;
}


/*
  We need to find a pair of isotopic atoms, iso1 and iso1+1
*/

static int
identify_atoms_with_isotopes (const Molecule & m,
                              const isotope_t iso1,
                              const isotope_t iso2,
                              atom_number_t & a1,
                              atom_number_t & a2)
{
  a1 = INVALID_ATOM_NUMBER;
  a2 = INVALID_ATOM_NUMBER;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    const auto iso = m.isotope(i);

    if (0 == iso)
      continue;

    if (iso == iso1)
    {
      a1 = i;
      if (INVALID_ATOM_NUMBER != a2)
        return 1;
    }
    else if (iso == iso2)
    {
      a2 = i;
      if (INVALID_ATOM_NUMBER != a1)
        return 1;
    }
  }

  return 0;    // should not happen
}

#ifdef NO_LONGER_USED
static int
reduce_to_scaffold (Molecule & m)
{
  const int matoms = m.natoms();

  Set_of_Atoms atoms_to_keep;

  for (int i = 0; i < matoms; ++i)
  {
    const isotope_t isoi = m.isotope(i);
    if (0 == isoi)
      continue;

    const Atom * a = m.atomi(i);
    const int acon = a->ncon();

    for (int j = 0; j < acon; ++j)
    {
      const atom_number_t k = a->other(i, j);

      const isotope_t isoj = m.isotope(j);

      if (0 == isoj)
        continue;

      m.remove_bond_between_atoms(i, k);

      if (isoi > isoj)
        atoms_to_keep.add(k);
      else
        atoms_to_keep.add(i);
    }
  }

  if (atoms_to_keep.empty()) {
    cerr << "reduce_to_scaffold:no labelled atoms found " << m.smiles() << endl;
    return 0;
  }

  const atom_number_t f = m.fragment_membership(atoms_to_keep[0]);

  m.delete_all_fragments_except(f);

  return 1;
}
#endif

atom_number_t
matched_atom_adjacent(const Molecule & m,
                      const Set_of_Atoms & e,
                      const atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);

    if (! e.contains(j))
      continue;

    return j;
  }

  return INVALID_ATOM_NUMBER;
}

int
Query_and_R_Atoms::_identify_attachment_points (Molecule & m,
                                                const Set_of_Atoms & e,
                                                const int ndx)
{
  _bond_lost[ndx] = INVALID_BOND_TYPE;

  const auto r = _r[ndx];

  const atom_number_t zatom = e[_m[ndx]];

  m.set_isotope(zatom, 2*r+1);
//cerr << "Set atom " << zatom << " type " << m.atomic_symbol(zatom) << " to " << (2*r+1) << endl;

  _atom1[ndx] = zatom;
  _atom2[ndx] = INVALID_ATOM_NUMBER;

  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const auto b = a->item(i);

    const auto j = b->other(zatom);

    if (! e.contains(j))
      continue;

    _bond_lost[ndx] = BOND_TYPE_ONLY(b->btype());
    _atom2[ndx] = j;
//  cerr << "Setting atom " << j << " type " << m.atomic_symbol(j) << " to " << (2*r) << endl;
    m.set_isotope(j, 2*r);

//  cerr << "Now " << m.smiles() << ' ' << m.name() << endl;
    return 1;
  }

  cerr << "Query_and_R_Atoms::_identify_attachment_points:no matched atom adjacent " << m.name() << endl;
  return 0;
}

/*
  We have matched. For each of our R atoms, we need to
  identify the bond being broken and the fragment
*/

int
Query_and_R_Atoms::process (Molecule & m,
                            const Set_of_Atoms & e,
                            IWString & output)
{
  cerr << "QQQQQQQQQ" << endl << "Calling process()" <<endl << "QQQQQQQQQ";
  Molecule mcopy(m);

  for (int i = 0; i < _nr; ++i)
  {
    if (! _identify_attachment_points(mcopy, e, i))    // ignore errors
      continue;
  }

  IWString buffer;

  buffer << m.smiles() << ' ' << m.name();    // always a space here

  buffer << _output_separator << mcopy.smiles() << ' ' << m.name();    // always a space

  int parent_written = 0;

  for (int i = 0; i < _nr; ++i)
  {
    mcopy.remove_bond_between_atoms(_atom1[i], _atom2[i]);
  }

//cerr << "After removing bonds " << mcopy.smiles() << endl;

  for (int i = 0; i < _nr; ++i)
  {
    const int r = _r[i];

    atom_number_t a1, a2;
    if (! identify_atoms_with_isotopes(mcopy, 2*r+1, 2*r, a1, a2))
    {
      cerr << "Query_and_R_Atoms::process:huh, no isotopes for break " << r << " in " << mcopy.smiles() << " " << m.name() << endl;
      return 0;
    }

    Molecule f;

    if (mcopy.number_fragments() > 1)
      mcopy.excise_fragment(a1, f);           // fragment f contains the fragment with atom A1

//  cerr << "i = " << i << " residual " << mcopy.smiles() << endl;
//  cerr << "frags " << f.smiles() << endl;

    f.remove_all(1);

    if (f.natoms() > _maximum_substituent_size)
      continue;

    ensure_parent_written(buffer, parent_written, output);

    if (0 == f.natoms())
    {
      output << _output_separator << 'H' << " R" << r;    // always a space
    }
    else
      output << _output_separator << f.unique_smiles() << " R" << r;    // always a space
  }

  if (! parent_written)
    return 1;

//Molecule mcopy(m);
//reduce_to_scaffold(mcopy);

  add_r_atoms_to_scaffold (mcopy);

  output << _output_separator << mcopy.unique_smiles() << ' ' << "SCAFFOLD " << m.name();

  output << '\n';

  return 1;
}

int
Query_and_R_Atoms::identify_attachment_points (Molecule & m,
                                               const Set_of_Atoms & e) const
{
  for (int i = 0; i < _nr; ++i)
  {
    const int ri = _r[i];
    const int mi = _m[i];

    const atom_number_t zatom = e[mi];

    m.set_isotope(zatom, 2*ri+1);

    atom_number_t j = matched_atom_adjacent(m, e, zatom);
    if (INVALID_ATOM_NUMBER == j)   // should not happen
      return 0;

    m.set_isotope(j, 2*ri);
  }

  return 1;
}

int
Query_and_R_Atoms::_remove_explicit_hydrogens (Molecule & m,
                                               int preserve_isotopic_labelled) const
{
  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; --i)
  {
    if (1 != m.atomic_number(i))
      continue;

    if (preserve_isotopic_labelled && m.isotope(i) > 0)
      continue;

    m.remove_atom(i);
    rc++;
  }

  return rc;
}

/*
  We just build up a large number that represents the matched atoms - in order
*/

int
Query_and_R_Atoms::is_duplicate_embedding (const int matoms, 
                                           const Set_of_Atoms & e,
                                           resizable_array<uint64_t> & a) const
{
  uint64_t s = 0;

  for (int i = 0; i < _nr; ++i)
  {
    const int r = _m[i];

    s += matoms * s + e[r];
  }

  if (a.contains(s))    // seen this before
    return 1;

  a.add(s);

  return 0;   // not a dup
}

static int
non_hydrogen_atoms (const Molecule & m)
{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; ++i)
  {
    if (1 != m.atomic_number(i))
      rc++;
  }

  return rc;
}

int
Query_and_R_Atoms::process_for_jibo (Molecule & m,
                                     const Set_of_Atoms & e,
                                     const int preserve_isotopic_labelled_hydrogen,
                                     IWString_and_File_Descriptor & output,
                                     IWString mainName,
                                     IWString mainSmiles)
{
  Molecule mcopy(m);

  for (int i = 0; i < _nr; ++i)
  {
    if (! _identify_attachment_points(mcopy, e, i))    // ignore errors
      continue;
  }

//#define DEBUG_PROCESS_FOR_JIBO
#ifdef DEBUG_PROCESS_FOR_JIBO
  cerr << "Starting smiles " << mcopy.smiles() << endl;
#endif

  IWString buffer;

  for (int i = 0; i < _nr; ++i)
  {
    const int r = _r[i];

    atom_number_t a1, a2;
    identify_atoms_with_isotopes(mcopy, 2*r, 2*r+1, a1, a2);
    int natoms = mcopy.natoms();
    
    if (INVALID_ATOM_NUMBER == a2)
    {
      cerr << "Query_and_R_Atoms::process_for_jibo:no isotopes for r = " << r << " in " <<  mcopy.smiles() << endl;
      return 0;
    }

    if (a1 < 0 || a1 >= natoms)
    {      
      cerr << "Query_and_R_Atoms::process_for_jibo:atom 1 (" << a1 << ") is invalid for r = " << r << " in " <<  mcopy.smiles() << endl;
      return 0;
    }

    if (a2 < 0 || a2 >= natoms)
    {      
      cerr << "Query_and_R_Atoms::process_for_jibo:atom 2 (" << a2 << ") is invalid for r = " << r << " in " <<  mcopy.smiles() << endl;
      return 0;
    }
    if (a1 == a2)
    {      
      cerr << "Query_and_R_Atoms::process_for_jibo:atoms 1 and 2 (" << a1 << ", " << a2 << ") are the same for r = " << r << " in " <<  mcopy.smiles() << " - the bond cannot be removed" << endl;
      return 0;
    }

    if (mcopy.are_bonded(a1, a2)) {
      cerr << "Atoms are bonded, should not happen " << mcopy.smiles() << '\n';
    } else {
      mcopy.remove_bond_between_atoms(a1, a2);
    }

    Molecule f;

    if (mcopy.number_fragments() > 1)
      mcopy.split_off_fragments(a1, f);             // fragment containing A1 is retained.

    _remove_explicit_hydrogens(f, preserve_isotopic_labelled_hydrogen);

    const int fatoms = non_hydrogen_atoms(f);

    if (fatoms > _maximum_substituent_size)
      continue;

//  if ( if we have a molecule to add to the sidechain, do it here. Need to pass as an argument...

    _append_smiles_and_molecular_descriptors(f, buffer);
  }

  add_r_atoms_to_scaffold(mcopy);

  mcopy.remove_all(1);
  if (_remove_isotopes_from_scaffold)
    mcopy.transform_to_non_isotopic_form();
 
 
  output << mainName << _output_separator << mainSmiles;


  output << _output_separator << mcopy.unique_smiles() << _output_separator << mcopy.natoms() << buffer;

  output << '\n';
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

class R1R2Etc
{
  private:
    int _verbose;
    int _molecules_read;
    Chemical_Standardisation _chemical_standardisation;
    int _reduce_to_largest_fragment;

    resizable_array_p<Query_and_R_Atoms> _queries;

    IWString_and_File_Descriptor _stream_for_molecules_hitting_no_queries;

    int _take_first_of_multiple_hits;
    int _ignore_molecules_hit_hitting_any_queries;
    int _unique_embeddings_only;
    int _molecules_not_hitting_any_queries;

    int  _maximum_substituent_size;
    float _maximum_substituent_fraction;   // not implemented

    int _jibo_output;

    int _fcsp3;             // or fccsp3
    
    int _brief_output;

    int _remove_isotopes_from_scaffold;

    int _translate_sidechain_isotopes_to;

    int _only_substituents_at_matched_atoms;

    int _allow_symmetry_related_matches;

    int _remove_all_chiral_centres;

    int _preserve_isotopic_labelled_hydrogen;

    Molecule _add_to_sidechains;

    char _output_separator;

// private functions

    void _default_values();
    void _usage (int rc);
    void _preprocess (Molecule & m);

    int _read_add_to_sidechains (iwstring_data_source & input);
    int _read_add_to_sidechains (const char * fname);

    int _remove_duplicate_embeddings (Molecule & m, Query_and_R_Atoms & q,
                                       Substructure_Results & sresults) const;

    int _unmatched_heavy_atom_available_instead_of_hydrogen (Molecule & m, const Set_of_Atoms & e) const;
    void _remove_extraeous_hydrogen_matches (Molecule & m, Query_and_R_Atoms * q, Substructure_Results & sresults) const;

    int _handle_zero_substructure_hits (Molecule & m, IWString_and_File_Descriptor & output);

    int _read_query (const char * fname, Molecule_to_Query_Specifications & mqs, int & r_atoms_in_queries);
    int _read_query (data_source_and_type<MDL_Molecule> &, Molecule_to_Query_Specifications & mqs, int & r_atoms_in_queries);
    int _read_query (MDL_Molecule &, Molecule_to_Query_Specifications & mqs, int & r_atoms_in_queries);

    int _write_jibo_header (const int nr, IWString_and_File_Descriptor &) const;
    int _output_for_jibo (Molecule & m, IWString_and_File_Descriptor & output);
    int _output_for_jibo (Molecule & m, Substructure_Results * sresults, const int * ndx, IWString_and_File_Descriptor & output);

    int _R1R2Etc (Molecule & m, const Set_of_Atoms & e, IWString_and_File_Descriptor & output);
    int _R1R2Etc (Molecule & m, IWString_and_File_Descriptor & output);
    int _R1R2Etc (data_source_and_type<Molecule> & input, IWString_and_File_Descriptor & output);
    int _R1R2Etc (const char * fname, FileType input_type, IWString_and_File_Descriptor & output);

  public:
    R1R2Etc();
    ~R1R2Etc();

    int operator() (int argc, char ** argv);
};

R1R2Etc::R1R2Etc ()
{
  _default_values();

  return;
}

void
R1R2Etc::_default_values()
{
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _take_first_of_multiple_hits = 0;
  _ignore_molecules_hit_hitting_any_queries = 0;
  _unique_embeddings_only = 0;
  _molecules_not_hitting_any_queries = 0;

  _maximum_substituent_size = std::numeric_limits<int>::max();
  _maximum_substituent_fraction = 1.0f;

  _jibo_output = 0;

  _fcsp3 = 1;
  
  _brief_output = 0;

  _remove_isotopes_from_scaffold = 0;

  _preserve_isotopic_labelled_hydrogen = 0;

  _translate_sidechain_isotopes_to = 0;

  _only_substituents_at_matched_atoms = 1;

  _allow_symmetry_related_matches = 1;

  _remove_all_chiral_centres = 0;

  _output_separator = ' ';

  return;
}

R1R2Etc::~R1R2Etc()
{
  return;
}

void
R1R2Etc::_usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Identifies substituents by query molecules labelled with 'R1', 'R2' atoms\n";
  cerr << " -q <fname>    query molecule - must have atoms labelled R1...\n";
  cerr << " -z i          ignore molecules not hitting any query\n";
  cerr << " -z f          in the case of multiple query matches, take the first\n";
  cerr << " -M sep        Maximum substituent size innumber of atoms.  Larger ones are ignored\n";
  cerr << " -o sep        defines the separator to be used. Values can be \"tab\" \"space\" \"comma\" \"newline\" or escaped single char\n";
  cerr << " -u            unique embeddings only\n";
  cerr << " -J <fname>    write output to <fname>\n";
  cerr << " -F            calculate fccsp3 rather than fcsp3\n";
  cerr << " -Y <fname>    echo first query to <fname>\n";
  cerr << " -I            remove isotopes from the scaffold structure\n";
  cerr << " -U 1          translate isotopes in fragments to 1\n";
  // -S is not implmented
  //cerr << " -S <fname>    molecule to add to sidechains (probably methane or benzene)\n";
  cerr << " -b            allow substituents at atoms other than the R labelled atoms\n";
  cerr << " -B            brief output - only the R1_smiles, R2_smiles, etc\n";
  cerr << " -Z <fname>    stream for smiles of molecules not reacting\n";
  cerr << " -k            turn OFF matching of symmetry equivalent matches\n";
  cerr << " -h            produce isotopically labelled Hydrogen attachment atoms\n";
  cerr << " -l            reduce to largest fragment\n";
  cerr << " -c            remove all chiral centers\n";
  cerr << " -i <type>     input specification\n";
  cerr << " -g ...        chemical standardisation options\n";
  cerr << " -E ...        standard element specifications\n";
  cerr << " -A ...        standard aromaticity specifications\n";
  cerr << " -v            verbose output\n";

  exit(rc);
}

void
R1R2Etc::_preprocess (Molecule & m)
{
  if (_reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (_chemical_standardisation.active())
    _chemical_standardisation.process(m);


  const int matoms = m.natoms();

  const Element * h = get_element_from_atomic_number(1);

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == m.implicit_hydrogens(i))
      continue;

    int x = m.natoms();
    m.add(h);
    m.add_bond(i, x, SINGLE_BOND);
  }

//m.make_implicit_hydrogens_explicit();

  return;
}

int
R1R2Etc::_handle_zero_substructure_hits (Molecule & m,
                                         IWString_and_File_Descriptor & output)
{
  _molecules_not_hitting_any_queries++;

  if (_verbose > 1 || ! _ignore_molecules_hit_hitting_any_queries)
    cerr << "None of " << _queries.number_elements() << " queries hit '" << m.name() << "'\n";

  if (_stream_for_molecules_hitting_no_queries.active())
  {
    _stream_for_molecules_hitting_no_queries << m.smiles() << ' ' << m.name() << '\n';
    _stream_for_molecules_hitting_no_queries.write_if_buffer_holds_more_than(4096);
  }

  if (_ignore_molecules_hit_hitting_any_queries)
    return 1;

  return 0;
}

/*
  in this particular invocation, we just process the first embedding in each SRESULTS object
*/

int
R1R2Etc::_output_for_jibo (Molecule & m,
                           Substructure_Results * sresults,
                           const int * ndx,
                           IWString_and_File_Descriptor & output)
{
  Molecule mcopy(m);
  mcopy.remove_all(1);

//  output << m.name() << _output_separator << mcopy.smiles();

  const int nq = _queries.size();

  for (int i = 0; i < nq; ++i)
  {
//  cerr << "Query " << i << " has " << sresults[i].number_embeddings() << " embeddings\n";
    if (0 == sresults[i].number_embeddings())
      continue;

    const auto e = sresults[i].embedding(ndx[i]);

    _queries[i]->process_for_jibo(m, *e, _preserve_isotopic_labelled_hydrogen, output,  m.name(), mcopy.smiles());
  }

  //output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

/*
  When we allow symmetric matches, we can get lots of matches, but the actual atoms appearing
  at the Rx positions may be duplicated
*/

int
R1R2Etc::_remove_duplicate_embeddings (Molecule & m,
                                       Query_and_R_Atoms & q,
                                       Substructure_Results & sresults) const
{
  const int matoms = m.natoms();

  resizable_array<uint64_t> a;

  for (int i = sresults.number_embeddings() - 1; i >= 0; --i)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    if (q.is_duplicate_embedding(matoms, *e, a))
    {
      sresults.remove_embedding(i);
    }
  }

  return 1;
}

//#define DEBUG_UNMATCHED_HEAVY_ATOM_AVAILABLE_INSTEAD_OF_HYDROGEN
/*
  We have multiple embeddings. If there are embeddings that include a Hydrogen atom, but
  if there is another unmatched atom at that site, discard that embedding.
  This arose from things like

  ClC1=CC2=C(NC(=O)C(=C2C2=CC=CC=C2)N2CCC(C(=O)O)CC2)C=C1 PBCHM16460913

  Where a query picked up two matches in the piperidine, one of which
  includes the explicit Hydrogen we created

  Our current implementation does not check that it is indeed an R atom
  that caused the Hydrogen match, may need to do that....
*/

int
R1R2Etc::_unmatched_heavy_atom_available_instead_of_hydrogen (Molecule & m,
                                      const Set_of_Atoms & e) const
{
  const int n = e.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const atom_number_t h = e[i];

#ifdef DEBUG_UNMATCHED_HEAVY_ATOM_AVAILABLE_INSTEAD_OF_HYDROGEN
    cerr << "Checking matched atom " << i << " atom " << h << " atomic number " << m.atomic_number(h) << endl;
#endif

    if (1 != m.atomic_number(h))
      continue;

    const atom_number_t x = m.other(h, 0);    // the atom to which the Hydrogen is attached

    const Atom * a = m.atomi(x);

    const int acon = a->ncon();

#ifdef DEBUG_UNMATCHED_HEAVY_ATOM_AVAILABLE_INSTEAD_OF_HYDROGEN
    cerr << "H bonded to " << x << " ncon " << acon << endl;
#endif

    for (int j = 0; j < acon; ++j)
    {
      const atom_number_t k = a->other(x, j);

      if (k == h)
        continue;

      if (1 == m.atomic_number(k))    // not really a heavy atom, just another H
        continue;

      if (e.contains(k))
        continue;

#ifdef DEBUG_UNMATCHED_HEAVY_ATOM_AVAILABLE_INSTEAD_OF_HYDROGEN
      cerr << "Attached heavy atom " << k << " not im embedding\n";
#endif

      return 1;    // there is an unmatched heavy atom here, ignore the Hydrogen match
    }
  }

  return 0;
}

void
R1R2Etc::_remove_extraeous_hydrogen_matches (Molecule & m,
                                             Query_and_R_Atoms * q,
                                             Substructure_Results & sresults) const
{
  for (int i = sresults.number_embeddings() - 1; i >= 0; --i)
  {
    const auto e = sresults.embedding(i);

    if (_unmatched_heavy_atom_available_instead_of_hydrogen(m, *e))
      sresults.remove_embedding(i);
  }

  return;
}
 
static int
next_index(const int n,
           const Substructure_Results * sresults,
           int * ndx)
{
  assert (n > 0);

  const int s = sresults[n-1].number_embeddings();

  if (0 == s)    // no matches to that query
  {
    if (n > 1)
      return next_index(n-1, sresults, ndx);
    return 0;
  }

  if (ndx[n-1] < s-1)
  {
    ndx[n-1]++;
    return 1;
  }

  if (1 == n)
    return 0;

  ndx[n-1] = 0;
  return next_index(n-1, sresults, ndx);
}

int
R1R2Etc::_output_for_jibo (Molecule & m,
                           IWString_and_File_Descriptor & output)
{
  Molecule_to_Match target(&m);

  int queries_hitting_this_molecule = 0;
  int multiple_matches_found = 0;

  const int nq = _queries.size();

  Substructure_Results * sresults = new Substructure_Results[nq]; std::unique_ptr<Substructure_Results[]> free_sresults(sresults);

  for (int i = 0; i < nq; ++i)
  {
    const int nhits = _queries[i]->substructure_search(target, sresults[i]);

    if (0 == nhits && _verbose > 2)
      cerr << "Processing " << m.name() << ", only matched " << _queries[i]->max_query_atoms_matched_in_search() << " query atoms\n";

    if (0 == nhits)
      continue;

    queries_hitting_this_molecule++;

//#define ECHO_QUERY_MATCHES_QWEQWEQWE
#ifdef ECHO_QUERY_MATCHES_QWEQWEQWE
    cerr << nhits << " hits to query " << i << endl;

    for (int k = 0; k < nhits; ++k)
    {
      const auto e = sresults[i].embedding(k);
      cerr << *e << endl;
    }
#endif

    if (nhits > 1)
    {
      _remove_extraeous_hydrogen_matches(m, _queries[i], sresults[i]);
//    cerr << "After _remove_extraeous_hydrogen_matches " << sresults[i].number_embeddings() << endl;

      _remove_duplicate_embeddings(m, *(_queries[i]), sresults[i]);

//    cerr << "After _remove_duplicate_embeddings " << sresults[i].number_embeddings() << endl;

      if (sresults[i].number_embeddings() > 1)
        multiple_matches_found++;
    }
  }

//cerr << "multiple_matches_found " << multiple_matches_found << endl;

  if (0 == queries_hitting_this_molecule)
    return _handle_zero_substructure_hits(m, output);

  int * ndx = new_int(nq); std::unique_ptr<int[]> free_ndx(ndx);

  if (0 == multiple_matches_found)
    return _output_for_jibo(m, sresults, ndx, output);

  do
  {
/*  for (int i = 0; i < nq; ++i)
    {
      cerr << " i = " << i << " ndx " << ndx[i];
    }
    cerr << endl;*/
    _output_for_jibo(m, sresults, ndx, output);
  } while (next_index(nq, sresults, ndx));

  return 1;
}

int
R1R2Etc::_R1R2Etc (Molecule & m,
                   IWString_and_File_Descriptor & output)
{
  Molecule_to_Match target(&m);

  int queries_hitting_this_molecule = 0;

  const int nq = _queries.size();

  for (int i = 0; i < nq; ++i)
  {
    Substructure_Results sresults;

    auto nhits = _queries[i]->substructure_search(target, sresults);

    if (_verbose > 2)
      cerr << "query " << i << ' ' << sresults.max_query_atoms_matched_in_search() << " query atoms matched\n";

    if (0 == nhits)
      continue;

    queries_hitting_this_molecule++;

    for (int j = 0; j < nhits; ++j)
    {
      const auto e = sresults.embedding(j);

      if (! _queries[i]->process(m, *e, output))
        return 0;
    }
  }

  if (queries_hitting_this_molecule)
    return 1;

  return _handle_zero_substructure_hits(m, output);
}

int
R1R2Etc::_R1R2Etc (data_source_and_type<Molecule> & input,
                   IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    _molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    _preprocess(*m);

    if (_jibo_output)
    {
      if (! _output_for_jibo(*m, output))
        return 0;
    }
    else if (! _R1R2Etc(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}


int
R1R2Etc::_R1R2Etc (const char * fname, FileType input_type,
                 IWString_and_File_Descriptor & output)
{
  assert (nullptr != fname);

  if (input_type == FILE_TYPE_INVALID)
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

  if (_verbose > 1)
    input.set_verbose(1);

  return _R1R2Etc(input, output);
}

int
R1R2Etc::_read_query (MDL_Molecule & m,
                      Molecule_to_Query_Specifications & mqs,
                      int & r_atoms_in_queries)
{
  Query_and_R_Atoms * q = new Query_and_R_Atoms();

  if (_remove_all_chiral_centres)
    m.remove_all_chiral_centres();

  q->set_only_substituents_at_matched_atoms(_only_substituents_at_matched_atoms);

  int tmp = 0;
  if (! q->build(m, mqs, tmp))
  {
    cerr << "R1R2Etc::_read_query:cannot build query from '" << m.name() << "'\n";
    delete q;
    return 0;
  }

  if (0 == r_atoms_in_queries)
    r_atoms_in_queries = tmp;

  if (_unique_embeddings_only)
    q->set_find_unique_embeddings_only(1);

  if (_take_first_of_multiple_hits)
    q->set_max_matches_to_find(1);

// Now copy over some variables we have

  q->set_output_separator(_output_separator);
  q->set_brief_output(_brief_output);
  q->set_maximum_substituent_size(_maximum_substituent_size);
  q->set_fcsp3(_fcsp3);
  q->set_remove_isotopes_from_scaffold(_remove_isotopes_from_scaffold);
  q->set_translate_sidechain_isotopes_to(_translate_sidechain_isotopes_to);
  q->set_allow_symmetry_related_matches(_allow_symmetry_related_matches);

  _queries.add(q);

  return 1;
}

int
R1R2Etc::_read_query (data_source_and_type<MDL_Molecule> & input,
                      Molecule_to_Query_Specifications & mqs,
                      int & r_atoms_in_queries)
{
  MDL_Molecule * m;

  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<MDL_Molecule> free_m(m);

    if (! _read_query(*m, mqs, r_atoms_in_queries))
    {
      cerr << "Fatal error processing query molecule " << input.molecules_read() << endl;
      return 0;
    }
  }

  return _queries.number_elements();
}

int
R1R2Etc::_read_query (const char * fname,
                      Molecule_to_Query_Specifications & mqs,
                      int & r_atoms_in_queries)
{
  const auto input_type = discern_file_type_from_name(fname);

  if (input_type == FILE_TYPE_INVALID)
  {
    cerr << "Cannot discern input type from '" << fname << "'\n";
    return 0;
  }

  data_source_and_type<MDL_Molecule> input (input_type, fname);

  if (! input.good())
  {
    cerr << "R1R2Etc::_read_query:cannot open '" << fname << "'\n";
    return 0;
  }

  return _read_query(input, mqs, r_atoms_in_queries);
}

int
R1R2Etc::_write_jibo_header (const int nr, IWString_and_File_Descriptor & output) const
{
  output << "CompoundName" << _output_separator << "SMILES" << _output_separator;
  output << "Scaffold" << _output_separator << "Scaffold_natoms";

  for (int i = 0; i < nr; ++i)
  {
    IWString rn;
    rn << _output_separator << 'R' << (i+1) << '_';
    output << rn << "SMILES";
    
    if (!_brief_output)
    {
      output << rn << "natoms";
      output << rn << "nrings";
      output << rn << "amw";
      output << rn << "ro5_ohnh";
      output << rn << "ro5_on";
      output << rn << "nvrtspsa";
      output << rn << "halogen";
      if (_fcsp3)
        output << rn << "fcsp3";
      else
        output << rn << "fccsp3";
      output << rn << "ringsys";
      output << rn << "arring";
      output << rn << "alring";
      output << rn << "ringatom";
      output << rn << "mxdst";
    }
  }

  output << '\n';

  return 1;
}

int
R1R2Etc::_read_add_to_sidechains (iwstring_data_source & input)
{
  const_IWSubstring buffer;
  if (! input.next_record(buffer))
  {
    cerr << "R1R2Etc::_read_add_to_sidechains:cannot read data\n";
    return 0;
  }

  if (! _add_to_sidechains.build_from_smiles(buffer))
  {
    cerr << "R1R2Etc::_read_add_to_sidechains:invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return _add_to_sidechains.natoms();
}

int
R1R2Etc::_read_add_to_sidechains (const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "R1R2Etc::_read_add_to_sidechains:cannot open '" << fname << "'\n";
    return 0;
  }

  return _read_add_to_sidechains(input);
}

int
R1R2Etc::operator() (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:i:g:lq:z:uo:M:Z:J:Y:FIU:bckhS:B");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    _usage(1);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, _verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      _usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, _verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g'))
  {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      _usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    _reduce_to_largest_fragment = 1;

    if (_verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('M'))
  {
    if (! cl.value('M', _maximum_substituent_size) || _maximum_substituent_size < 1)
    {
      cerr << "R1R2Etc::the maximum substituent size (-M) must be a whole +ve number\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "Will not write fragments with more than " << _maximum_substituent_size << " atoms\n";
  }

  if (cl.option_present('F'))
  {
    _fcsp3 = 0;

    if (_verbose)
      cerr << "Will compute fccsp3 rather than fcsp3\n";
  }
 
  if (cl.option_present('B'))
  {
    _brief_output = 1;

    if (_verbose)
      cerr << "Will output only the parent and fragment smiles\n";
  }



  if (cl.option_present('I'))
  {
    _remove_isotopes_from_scaffold = 1;

    if (_verbose)
      cerr << "Will remove isotopes from the scaffold\n";
  }

  if (cl.option_present('U'))
  {
    if (! cl.value('U', _translate_sidechain_isotopes_to) || _translate_sidechain_isotopes_to < 0)
    {
      cerr << "R1R2Etc::the common isotope for sidechains (-U) must be a valid isotope\n";
      _usage(1);
    }

    if (_verbose)
      cerr << "Will label all sidechain attachment points with isotope " << _translate_sidechain_isotopes_to << endl;
  }

  if (cl.option_present('b'))
  {
    _only_substituents_at_matched_atoms = 0;

    if (_verbose)
      cerr << "Will allow substituents at other locations in the query\n";
  }

  if (cl.option_present('c'))
  {
    _remove_all_chiral_centres = 1;

    if (_verbose)
      cerr << "Will remove chiral centres from query molecules\n";
  }

  if (cl.option_present('k'))
  {
    _allow_symmetry_related_matches = 0;

    if (_verbose)
      cerr << "Will NOT perceive symmetry equivalent matches\n";
  }

  if (cl.option_present('h'))
  {
    _preserve_isotopic_labelled_hydrogen = 1;

    if (_verbose)
      cerr << "Will produce isotopically labelled Hydrogen substituents\n";
  }

  if (cl.option_present('S'))
  {
    const char * fname = cl.option_value('S');
    if (! _read_add_to_sidechains(fname))
    {
      cerr << "Cannot read molecule to add to sidechains '" << fname << "'\n";
      return 1;
    }

    if (_verbose)
      cerr << "Will add " << _add_to_sidechains.smiles() << " to sidechains\n";
  }

  if (cl.option_present('z'))
  {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i)
    {
      if ('f' == z)
      {
        _take_first_of_multiple_hits = 1;

        if (_verbose)
          cerr << "If multiple query matches, will take the first match\n";
      }
      else if ('i' == z)
      {
        _ignore_molecules_hit_hitting_any_queries = 1;

        if (_verbose)
          cerr << "Will ignore molecules not hitting any queries\n";
      }
      else
      {
        cerr << "R1R2Etc::unrecognised -z qualifier '" << z << "'\n";
        _usage(1);
      }
    }
  }

  if (cl.option_present('u'))
  {
    _unique_embeddings_only = 1;

    if (_verbose)
      cerr << "Will only match unique embeddings\n";
  }

  if (! cl.option_present('q'))
  {
    cerr << "Must specify one or more query molecules via the -q option\n";
    _usage(1);
  }

  if (cl.option_present('o'))     // must be processed before -q options
  {
    const_IWSubstring o = cl.string_value('o');
    if ("tab" == o)
      _output_separator = '\t';
    else if ("space" == o)
      _output_separator = ' ';
    else if ("comma" == o)
      _output_separator = ',';
    else if ("newline" == o)
      _output_separator = '\n';
    else if (1 == o.length())
      _output_separator = o[0];
    else
    {
      cerr << "Unrecognised -o qualifier '" << o << "'\n";
      return 2;
    }
  }

  for (int i = 0; i < 10; ++i)
  {
    IWString r;
    r << 'R' << i;
    create_element_with_symbol(r.c_str(), r.length());
  }

  Molecule_to_Query_Specifications mqs;

  mqs.set_ring_atoms_conserve_ring_membership(1);

  int r_atoms_in_queries = 0;

  if (1)
  {
    IWString q;
    for (int i = 0; cl.value('q', q, i); ++i)
    {
      if (! _read_query (q.c_str(), mqs, r_atoms_in_queries))
      {
        cerr << "Cannot process query molecule '" << q << "'\n";
        return (i+1);
      }
    }

    for (int i = 1; i < _queries.number_elements(); ++i)
    {
      if (_queries[i-1]->number_R_atoms() != _queries[i]->number_R_atoms())
      {
        cerr << "R1R2Etc::all queries must have the same number of R atoms, mismatch " << _queries[i-1]->number_R_atoms() << " and " << _queries[i]->number_R_atoms() << endl;
        return i;
      }
    }

    if (cl.option_present('Y'))
    {
      IWString stem;
      cl.value('Y', stem);

      for (int i = 0; i < _queries.number_elements(); ++i)
      {
        IWString fname;
        fname << stem << i << ".qry";

//      _queries[i]->debug_print(cerr);
        _queries[i]->write_msi(fname);
      }

      if (_verbose)
        cerr << "Echo'd " << _queries.number_elements() << " queries to stem '" << stem << "'\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      _usage (6);
    }
  }
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    _usage(2);
  }

  if (cl.option_present('Z'))
  {
    IWString z = cl.option_value('Z');

    if (! z.ends_with(".smi"))
      z << ".smi";

    if (! _stream_for_molecules_hitting_no_queries.open(z.c_str()))
    {
      cerr << "R1R2Etc::Cannot open stream for molecules not hitting queries '" << z << "'\n";
      return 2;
    }

    if (_verbose)
      cerr << "Molecules not hitting any queries written to '" << z << "'\n";
  }

  int rc = 0;

  if (cl.option_present('J'))
  {
    const char * j = cl.option_value('J');
    IWString_and_File_Descriptor output;
    if (! output.open(j))
    {
      cerr << "Cannot open stream for Jibo output '" << j << "'\n";
      return 1;
    }

    _jibo_output =1;

    _write_jibo_header(r_atoms_in_queries, output);

    for (int i = 0; i < cl.number_elements(); ++i)
    {
      if (! _R1R2Etc(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }
    output.flush();
  }
  else
  {
    IWString_and_File_Descriptor output(1);

    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! _R1R2Etc(cl[i], input_type, output))
      {
        rc = i + 1;
        break;
      }
    }

    output.flush();
  }

  if (_verbose)
  {
    cerr << "Read " << _molecules_read << " molecules\n";
    if (0 == _molecules_not_hitting_any_queries)
      cerr << "All molecules hit a query\n";
    else
      cerr << _molecules_not_hitting_any_queries << " molecules did not hit any queries\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  R1R2Etc R1R2Etc;
  
  return R1R2Etc (argc, argv);
}
