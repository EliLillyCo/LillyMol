#include <iostream>
#include <memory>

#include "re2/re2.h"

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/path.h"

#include "tp_first_pass_lib.h"

namespace lilly_medchem_rules {

using std::cerr;
using std::endl;

template <typename F>
F
Fraction(int numerator, int denominator) {
  return static_cast<F>(numerator) / static_cast<F>(denominator);
}

// Returns true if `value` is in [0, 1].
template <typename T>
bool
OkFraction(T value) {
  if (value < static_cast<T>(0)) {
    return false;
  }
  if (value > static_cast<T>(1.0)) {
    return false;
  }
  return true;
}

int
CountOrFraction::Build(const const_IWSubstring& token) {
  re2::StringPiece tmp(token.data(), token.length());
  if (RE2::FullMatch(tmp, "^(\\d+)$", &count)) {
    return 1;
  }

  if (! token.numeric_value(fraction) || ! OkFraction(fraction)) {
    cerr << "CountOrFraction::Build:invalid fraction '" << token << "'\n";
    return 0;
  }

  return 1;
}

bool
CountOrFraction::GreaterThan(int c, int tot) const {
  if (count >= 0) {
    return c > count;
  }
  if (fraction >= 0.0f) {
    return Fraction<float>(c, tot) >= fraction;
  }

  return false;
}

// Return true if `m` contains both a Carbon atom as well
// as at least a N or O.
int
InterestingAtoms(const Molecule & m) {
  int carbon = 0;
  int nitrogen = 0;
  int oxygen = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);
    if (6 == z) {
      carbon = 1;
      if (nitrogen || oxygen)
        return 1;
    } else if (7 == z) {
      nitrogen = 1;
      if (carbon)
        return 1;
    } else if (8 == z) {
      oxygen = 1;
      if (carbon)
        return 1;
    }
  }
  
  return 0;
}


MCFirstPass::MCFirstPass() {
// initialise the array of allowable elements

  std::fill_n(_ok_elements, HIGHEST_ATOMIC_NUMBER + 1, 0);

  _ok_elements[6] = 1;
  _ok_elements[7] = 1;
  _ok_elements[8] = 1;
  _ok_elements[9] = 1;
  _ok_elements[15] = 1;
  _ok_elements[16] = 1;
  _ok_elements[17] = 1;
  _ok_elements[35] = 1;
  _ok_elements[53] = 1;
  _ok_elements[3]  = 1;    // Li
  _ok_elements[11] = 1;    // Na
  _ok_elements[12] = 1;    // Mg
  _ok_elements[19] = 1;    // K
  _ok_elements[20] = 1;    // Ca
}

// Chirality is currently not handled in this tool.
#ifdef DELETE_SOON
void
DisplayChiralityOptions(std::ostream& output) {
  output << "  -s discard     discard all chiral information in the input\n";
  output << "  -s good        ignore erroneous chiral input\n";
  output << "  -s 1           include chiral info in output (default)\n";
  output << "  -s 0           exclude chiral info from output\n";
}
#endif

//#define DEBUG_MC_FIRST_PASS

int
MCFirstPass::Build(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('e')) {  // OK elements.
    IWString e;
    for (int i = 0; cl.value('e', e, i); ++i) {
      const Element * o = get_element_from_symbol_no_case_conversion(e);
      if (nullptr == o) {
        cerr << "Sorry, non periodic table element '" << e << "', cannot be OK\n";
        return 0;
      }

      atomic_number_t z = o->atomic_number();

      assert (z >= 0 && z <= HIGHEST_ATOMIC_NUMBER);

      _ok_elements[z] = 1;

      if (_verbose)
        cerr << "Element " << e << " atomic number " << z << " allowed\n";

      if (5 == z || 14 == z)
      {
        Element * x = const_cast<Element *>(o);
        x->set_organic(1);
      }
    }
  }

  if (cl.option_present('n')) {
    IWString ele;
    for (int i = 0; cl.value('n', ele, i); ++i) {
      const Element * e = get_element_from_symbol_no_case_conversion(ele);
      if (e == NULL) {
        cerr << "Unrecognised element '" << ele << "'\n";
        return 1;
      }
      const_cast<Element*>(e)->set_organic(0);
      if (_verbose) {
        cerr << "Element '" << ele << "' marked non organic\n";
      }
    }
  }

  for (int i = 1; i < HIGHEST_ATOMIC_NUMBER; i++)
  {
    if (_ok_elements[i])
      continue;

    const Element * e = get_element_from_atomic_number(i);
    if (! e->organic())
      continue;

    _ok_elements[i] = 1;
    if (_verbose)
      cerr << "Element " << e->symbol() << " atomic number " << i << " allowed\n";
  }

  if (cl.option_present('t')) {
    if (! _element_transformations.construct_from_command_line(cl, _verbose, 't')) {
      return 0;
    }
  }

  if (cl.option_present('X')) {
    if (! _elements_to_remove.construct_from_command_line(cl, _verbose, 'X')) {
      cerr << "Cannot discern elements to remove from -X switch\n";
      return 0;
    }
  }

  if (cl.option_present('k')) {
    _exclude_molecules_with_no_interesting_atoms = 0;

    if (_verbose)
      cerr << "Will allow molecules with no interesting atoms to pass\n";
  }

  if (cl.option_present('f')) {
    if (! _exclude_molecules_with_no_interesting_atoms) {
      cerr << "The -f and -k options cannot be used together\n";
      return 0;
    }

    if (! cl.value('f', _min_fraction_interesting_atoms) ||
        ! OkFraction(_min_fraction_interesting_atoms)) {
      cerr << "The minimum fraction of interesting atoms needed (-f) must be a valid fraction\n";
      return 0;
    }

    if (_verbose)
      cerr << "Molecules must have a minimum fraction " << _min_fraction_interesting_atoms << " of interesting atoms\n";
  }

  if (cl.option_present('y')) {
    _allow_non_periodic_table_elements_if_not_connected = 1;

    if (_verbose)
      cerr << "Non periodic table elements allowed if not connected\n";
  }

  if (! cl.option_present('c'))
    ;
  else if (cl.value ('c', _lower_atom_count_cutoff) && _lower_atom_count_cutoff > 0) {
    if (_verbose)
      cerr << "Will exclude molecules with fewer than " << _lower_atom_count_cutoff << " atoms\n";
  } else {
    cerr << "Cannot discern lower atom count cutoff from " << cl.option_value ('c') << "'\n";
    return 0;
  }

  if (! cl.option_present('C'))
    ;
  else if (cl.value('C', _upper_atom_count_cutoff))
  {
    if (_upper_atom_count_cutoff < _lower_atom_count_cutoff)
    {
      cerr << "Upper atom count cutoff " << _upper_atom_count_cutoff << 
              " must be greater than lower atom count cutoff " << _lower_atom_count_cutoff << '\n';
      return 0;
    }
    if (_verbose)
      cerr << "Will exclude molecules with more than " << _upper_atom_count_cutoff << " atoms\n";
  } else {
    cerr << "Cannot discern upper atom count cutoff from '" << cl.option_value ('C') << "'\n";
    return 0;
  }

  if (cl.option_present('r')) {
    if (! cl.value ('r', _lower_ring_count_cutoff) ||
          _lower_ring_count_cutoff < 1)
    {
      cerr << "-r option needs a whole number > 0\n";
      return 0;
    }

    if (_verbose)
      cerr << "Molecules containing fewer than " << _lower_ring_count_cutoff <<
              " will be ignored\n";
  }

  if (cl.option_present('R')) {
    if (! cl.value('R', _upper_ring_count_cutoff) ||
          _upper_ring_count_cutoff < _lower_ring_count_cutoff) {
      cerr << "-R option needs a whole number > " << _lower_ring_count_cutoff << '\n';
      return 0;
    }

    if (_verbose)
      cerr << "Molecules containing more than " << _upper_ring_count_cutoff <<
              " rings will be ignored\n";
  }

  if (cl.option_present('Z')) {
    if (! cl.value('Z', _upper_ring_size_cutoff) || _upper_ring_size_cutoff < 3)
    {
      cerr << "The upper ring size cutoff value (-Z) must be a valid ring size\n";
      return 0;
    }

    if (_verbose)
      cerr << "Will discard molecules with ring sizes > " << _upper_ring_size_cutoff << '\n';
  }

  if (cl.option_present('V'))
  {
    _skip_molecules_with_abnormal_valences = 1;
    if (_verbose)
      cerr << "Molecules containing abnormal valences will be skipped\n";
  }

  if (cl.option_present('I')) {
    IWString tmp;
    for (int i = 0; cl.value('I', tmp, i); ++i) {
      if ('0' == tmp) {
        _exclude_isotopes = 1;
        if (_verbose)
          cerr << "Molecules containing isotopes will be excluded\n";
      } else if ('1' == tmp) {
        _exclude_isotopes = 0;
        if (_verbose)
          cerr << "No action taken on molecules containing isotopes\n";
      } else if ("convert" == tmp) {
        _convert_isotopes = 1;
        if (_verbose)
          cerr << "All isotopic atoms will be converted to their non isotopic form\n";
      } else {
        cerr << "Unrecognised -I qualifier '" << tmp << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('p')) {
    const_IWSubstring p;
    cl.value ('p', p);

    if ("remove" == p) {
      _remove_non_printing_chars = 1;
      if (_verbose)
        cerr << "Non printing characters removed from molecule names\n";
    } else if (1 == p.length ()) {
      _translate_non_printing_chars = p[0];
      if (_verbose)
        cerr << "Non printing characters translated to '" << _translate_non_printing_chars << "'\n";
    } else {
      cerr << "Invalid/unrecognised -p qualifier '" << p << "'\n";
      return 0;
    }
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', _max_ring_bond_ratio) || OkFraction(_max_ring_bond_ratio)) {
      cerr << "The lower ring bond ratio (-b) option must be followed by a valid fraction\n";
      return 0;
    }

    if (_verbose)
      cerr << "Molecules with ring bond ratio's of " << _max_ring_bond_ratio << " or less are rejected\n";
  }

  if (cl.option_present('H')) {
    if (! cl.value('H', _upper_ring_system_size) || _upper_ring_system_size < 2) {
      cerr << "The upper ring system size (-H) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will reject molecules containing ring systems with more than " << _upper_ring_system_size << " rings\n";
    }
  }

  // Count or fraction specifications look like nn:cf
  if (cl.option_present('F')) {
    const_IWSubstring token;
    for (int i = 0; cl.value('F', token, i); ++i) {
      const_IWSubstring nn, cf;
      if (! token.split(nn, ':', cf) || nn.empty() || cf.empty()) {
        cerr << "Invalid elemental fraction specification '" << token << "'\n";
        return 0;
      }
      atomic_number_t z;
      if (! nn.numeric_value(z) || z < 1 || z> HIGHEST_ATOMIC_NUMBER) {
        cerr << "Invalid atomic number '" << token << "'\n";
        return 0;
      }
      CountOrFraction count_or_fraction;
      if (! count_or_fraction.Build(cf)) {
        cerr << "Invalid count or fraction specification '" << token << "'\n";
        return 0;
      }
      _element_fraction[z] = count_or_fraction;
    }
  }

  return 1;
}

int
MCFirstPass::set_ok_element(const Element* e) {
  if (! e->is_in_periodic_table()) {
    cerr << "MCFirstPass::set_ok_element:non periodic table element '" << e->symbol() << "'\n";
    return 0;
  }
  _ok_elements[e->atomic_number()] = 1;
  return 1;
}

int
MCFirstPass::Rejected(Molecule& m,
                      MCFirstPassCounter& counter,
                      IWString& rejection_reason) {
  counter.molecules_processed++;

  rejection_reason.resize_keep_storage(0);

  _elements_to_remove.process(m);

  if (RejectedInner(m, counter, rejection_reason)) {
    counter.molecules_rejected++;
    return 1;
  }

  if (_convert_isotopes)
    ConvertIsotopes(m, counter);

  return 0;  // Not rejected.
}

int
MCFirstPass::RejectedInner(Molecule& m,
                      MCFirstPassCounter& counter,
                      IWString& rejection_reason) {
  const int matoms = m.natoms();
  if (matoms == 0) {
    counter.empty_molecule++;
    rejection_reason = "Empty molecule";
    return 1;
  }
  
  _element_transformations.process(m);

  if (ExcludeForAtomType(m, counter, rejection_reason)) {
    return 1;
  } 
  
  if (_skip_molecules_with_abnormal_valences && ! m.valence_ok()) {
    counter.molecules_with_abnormal_valences++;
    if (_verbose > 1)
      cerr << "Molecule contains abnormal valence(s)\n";
    rejection_reason << "abnormal_valence";
    return 1;
  }

  if (! FragmentsAreOK(m, counter, rejection_reason)) {
    return 1;
  }

  if (_exclude_molecules_with_no_interesting_atoms &&
            ExcludeForNoInterestingAtoms(m, counter, rejection_reason)) {
    return 1;
  } 
  
  if (_min_fraction_interesting_atoms > 0.0 &&
                ExcludeForTooFewInterestingAtoms(m, counter, rejection_reason)) {
    return 1;
  } 
  
  if (_lower_ring_count_cutoff &&
           m.nrings() < _lower_ring_count_cutoff) {
    counter.molecules_with_too_few_rings++;
    if (_verbose > 1)
      cerr << "Molecule contains " << m.nrings() << " rings, which is below cutoff\n";
    rejection_reason << "not_enough_rings";
    return 1;
  } 
  
  if (_upper_ring_count_cutoff &&
           m.nrings() > _upper_ring_count_cutoff) {
    counter.molecules_with_too_many_rings++;
    if (_verbose > 1)
      cerr << "Molecule contains " << m.nrings() << " rings, which is above cutoff\n";
    rejection_reason << "too_many_rings";
    return 1;
  } 
  
  if (_exclude_isotopes && m.number_isotopic_atoms()) {
    counter.molecules_containing_isotopes++;
    if (_verbose > 1)
      cerr << "Molecule contains isotopes\n";
    rejection_reason << "isotopes";
    return 1;
  } 
  
  if (RejectForRingBondRatio(m, counter, rejection_reason)) {
    return 1;
  } 
  
  if (_upper_ring_size_cutoff > 0 &&
             RejectForRingSizeCondition(m, counter, rejection_reason)) {
    return 1;
  }

  if (RejectedForElementFraction(m, counter, rejection_reason)) {
    return 1;
  }

  if (RejectedForRingSystemSize(m, counter, rejection_reason)) {
    return 1;
  }

  return 0;
}

int
MCFirstPass::RejectedForElementFraction(Molecule& m,
                MCFirstPassCounter& count,
                IWString& rejection_reason) const {
  if (_element_fraction.empty()) {
    return 0;
  }

  const int matoms = m.natoms();

  std::unordered_map<atomic_number_t, int> elements = m.ElementCount();
  for (const auto& [z, count_or_fraction] : _element_fraction) {
    const auto iter = elements.find(z);
    if (iter == elements.end()) {
      continue;
    }
    if (count_or_fraction.GreaterThan(iter->second, matoms)) {
      rejection_reason << "element_fraction";
      return 1;
    }
  }

  return 0;  // Not rejected.
}

int
MCFirstPass::RejectedForRingSystemSize(Molecule& m,
                MCFirstPassCounter& counter,
                IWString& rejection_reason) const {
  const int nr = m.nrings();
  if (nr < 2) {
    return 0;   // NOt rejected.
  }

  std::unique_ptr<int[]> ring_alread_done(new_int(nr));

  for (int i = 0; i < nr; ++i) {
    if (ring_alread_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);
    if (! ri->is_fused()) {
      continue;
    }

    int system_size = 1;
    for (int j = i + 1; j < nr; ++j) {
      if (ring_alread_done[j]) {
        continue;
      }
      const Ring* rj = m.ringi(j);
      if (ri->fused_system_identifier() != rj->fused_system_identifier()) {
        continue;
      }
      ring_alread_done[j] = 1;
      ++system_size;
    }
    if (system_size > _upper_ring_system_size) {
      counter.molecules_with_bad_ring_system_size++;
      rejection_reason << "ring_system_size";
      return 1;
    }
  }

  return 0;  // Not rejected.
}

bool
MCFirstPass::ExcludeForAtomType(const Molecule& m,
                             MCFirstPassCounter& counter,
                             IWString& rejection_reason) const {
#ifdef DEBUG_MC_FIRST_PASS
  cerr << "ExcludeForAtomType\n";
#endif
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    const Element * e = a->element();

    if (1 == e->atomic_number() && 1 != a->ncon())
    {
      if (_verbose > 1)
        cerr << "Contains two valent hydrogen\n";
      counter.molecules_with_two_connected_hydrogens++;
      rejection_reason = "two_valent_hydrogen";
      return 1;    // yes, exclude this molecule
    }

    if (e->organic())
      continue;

    if (! e->is_in_periodic_table()) {
      // If not excluded, they are OK.
      if (! _exclude_molecules_containing_non_periodic_table_elements) {
        continue;
      }

      if (_allow_non_periodic_table_elements_if_not_connected && 0 == a->ncon()) {
        if (_verbose > 1)
          cerr << "Allowing singly connected non-periodic table element '" << e->symbol() << "'\n";
        continue;
      }

      if (_verbose > 1)
        cerr << "Contains non periodic table atom '" << e->symbol() << "'\n";
      counter.molecules_containing_non_periodic_table_elements++;
      rejection_reason = "non_periodic_table_atom";
      return 1;    // yes, exclude this molecule
    }

    if (! _ok_elements[e->atomic_number()])
    {
      if (_verbose > 1)
        cerr << "Contains non-allowed atom '" << e->symbol() << "'\n";
      counter.molecules_containing_non_allowed_atom_types++;
      rejection_reason = "non_allowed_atom";
      return 1;    // yes, exclude this molecule
    }

    if (a->ncon() > 0)
    {
      if (_verbose > 1)
        cerr << "Contains covalently bound non-organic '" << e->symbol() << "'\n";
      counter.molecules_containing_colvalent_non_organics++;
      rejection_reason = "covalent_non-organic";
      return 1;    // yes, exclude this molecule
    }
  }

  return 0;    // not rejected
}

int
MCFirstPass::FragmentsAreOK(Molecule & m,
                        MCFirstPassCounter& counter,
                        IWString& rejection_reason) const {
#ifdef DEBUG_MC_FIRST_PASS
  cerr << "FragmentsAreOK\n";
#endif
  const int nf = m.number_fragments();
  if (nf == 1) {
    return OkLowerUpperAtomCountCutoff(m.natoms(), counter, rejection_reason);
  }

  int atoms_in_largest_fragment = 0;
  int number_large_fragments = 0;

  for (int i = 0; i < nf; i++)
  {
    const int aif = m.atoms_in_fragment(i);

    if (aif < atoms_in_largest_fragment) {
      continue;
    }

    if (aif > _reject_if_any_fragment_larger_than) {
      if (_verbose > 1) {
        cerr << "Too many atoms " << aif << " in fragment\n";
      }
      counter.molecules_with_fragment_too_large++;
      rejection_reason << "fragment too large " << aif;
      return 0;
    }

    if (aif > atoms_in_largest_fragment) {
      number_large_fragments = 1;
      atoms_in_largest_fragment = aif;
    } else if (aif == atoms_in_largest_fragment) {
      number_large_fragments++;
    }
  }

  if (number_large_fragments > 1)
  {
    if (_verbose > 1)
      cerr << "Mixture " << number_large_fragments << " large fragments\n";
    rejection_reason << "mixture";
    counter.mixtures_rejected++;
    return 0;
  }

  return OkLowerUpperAtomCountCutoff(atoms_in_largest_fragment, counter, rejection_reason);
}

int
MCFirstPass::OkLowerUpperAtomCountCutoff(const int natoms,
                                         MCFirstPassCounter& counter,
                                         IWString& rejection_reason) const {

  if (_upper_atom_count_cutoff > 0 && natoms > _upper_atom_count_cutoff)
  {
    if (_verbose > 1)
      cerr << "Too many atoms " << natoms << '\n';
    rejection_reason = "too_many_atoms";
    counter.molecules_above_atom_count_cutoff++;
    return 0;
  }

  if (natoms < _lower_atom_count_cutoff)
  {
    counter.molecules_below_atom_count_cutoff++;
    rejection_reason = "too_few_atoms";
    if (_verbose > 1)
      cerr << "Too few atoms " << natoms << '\n';
    return 0;
  }

  return 1;
}

double
Fraction(int numerator, int denominator) {
  return static_cast<double>(numerator) /
         static_cast<double>(denominator);
}

int
MCFirstPass::NonPrintingCharactersInName(Molecule & m) const {
  IWString mname = m.name();

  int rc = 0;
  for (int i = mname.length() - 1; i >= 0; i--) {
    char c = mname[i];

    if (isprint(c))
      continue;

    if ('\0' != _translate_non_printing_chars)
      mname[i] = _translate_non_printing_chars;
    else if (_remove_non_printing_chars) {
      mname.remove_item(i);
      i++;
    }

    rc++;
  }

  if (rc == 0) {
    return rc;
  }

  if ('\0' != _translate_non_printing_chars || _remove_non_printing_chars) {
    m.set_name(mname);
  }

  if (_verbose > 1) {
    cerr << "Removed/changed " << rc << " non printing chars in name\n";
  }

  return rc;
}

int
interesting_atoms(const Molecule & m)
{
  int carbon = 0;
  int nitrogen = 0;
  int oxygen = 0;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t z = m.atomic_number(i);
    if (6 == z)
    {
      carbon = 1;
      if (nitrogen || oxygen)
        return 1;
    }
    else if (7 == z)
    {
      nitrogen = 1;
      if (carbon)
        return 1;
    }
    else if (8 == z)
    {
      oxygen = 1;
      if (carbon)
        return 1;
    }
  }
  
  return 0;
}

static int
count_interesting_atoms(const Molecule & m,
                        const int * include_atom,
                        const int flag)
{
  int carbon = 0;
  int nitrogen = 0;
  int oxygen = 0;

  const auto matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (nullptr != include_atom && flag != include_atom[i])
      continue;

    atomic_number_t z = m.atomic_number(i);
    if (6 == z)
      carbon++;
    else if (7 == z)
      nitrogen++;
    
    else if (8 == z)
      oxygen++;
  }

  if (0 == oxygen && 0 == nitrogen)
    return 0;

  if (0 == carbon)
    return 0;
  
  return nitrogen + oxygen;
}

int
MCFirstPass::ExcludeForNoInterestingAtoms(Molecule & m,
                             MCFirstPassCounter& counter,
                             IWString& rejection_reason) const {
#ifdef DEBUG_MC_FIRST_PASS
  cerr << "ExcludeForNoInterestingAtoms\n";
#endif
  int nf = m.number_fragments();
  if (1 == nf)
  {
    if (InterestingAtoms(m))
      return 0;    // do not reject this molecule
    else
    {
      rejection_reason = "no_interesting_atoms";
      counter.molecules_with_no_interesting_atoms++;
      return 1;    // yes, reject this molecule
    }
  }

  resizable_array_p<Molecule> components;
  m.create_components(components);

  int largest_frag = 0;
  int interesting_atoms_in_largest_frag = 0;

  for (int i = 0; i < nf; i++)
  {
    Molecule * c = components[i];

    int catoms = c->natoms();

    assert (catoms == m.atoms_in_fragment(i));

    if (catoms < largest_frag)
      continue;

    interesting_atoms_in_largest_frag = interesting_atoms(*c);
    largest_frag = catoms;
  }

  if (interesting_atoms_in_largest_frag)
    return 0;     // do not exclude it
  else
  {
    rejection_reason = "no_interesting_atoms";
    counter.molecules_with_no_interesting_atoms++;
    return 1;     // yes, exclude this molecule
  }
}

int
MCFirstPass::ExcludeForTooFewInterestingAtoms(Molecule & m,
                                  MCFirstPassCounter& counter,
                                  IWString& rejection_reason) const {
#ifdef DEBUG_MC_FIRST_PASS
  cerr << "ExcludeForTooFewInterestingAtoms\n";
#endif
  const int nf = m.number_fragments();

  if (nf == 1) {
    const int interesting_atoms = count_interesting_atoms(m, nullptr, 0);
    if (Fraction(interesting_atoms, m.natoms()) >=
                _min_fraction_interesting_atoms)
      return 0;                // not rejected
    counter.molecules_with_too_few_interesting_atoms++;
    rejection_reason = "too few interesting atoms";
    return 1;          // yes, exclude this molecule
  }

  const int matoms = m.natoms();

  int * f = new int[matoms]; std::unique_ptr<int[]> free_f(f);

  m.fragment_membership(f);

  int atoms_in_largest_fragment = 0;
  int interesting_atoms_in_largest_fragment = 0;

  for (int i = 0; i < nf; ++i)
  {
    const auto interesting_atoms = count_interesting_atoms(m, f, i+1);
    const auto aif = m.atoms_in_fragment(i);

    if (aif > atoms_in_largest_fragment)
    {
      atoms_in_largest_fragment = aif;
      interesting_atoms_in_largest_fragment = interesting_atoms;
    }
    else if (aif == atoms_in_largest_fragment && interesting_atoms > interesting_atoms_in_largest_fragment)
      interesting_atoms_in_largest_fragment = interesting_atoms;
  }

  if (Fraction(interesting_atoms_in_largest_fragment, atoms_in_largest_fragment) >=
          _min_fraction_interesting_atoms)
    return 0;            // do not reject this molecule

  rejection_reason = "Too_few_interesting_atoms";
  counter.molecules_with_too_few_interesting_atoms++;

  return 1;        // yes, reject this molecule
}

int
MCFirstPass::RejectForRingBondRatio(Molecule & m,
                        MCFirstPassCounter& counter,
                        IWString& rejection_reason) const {
#ifdef DEBUG_MC_FIRST_PASS
  cerr << "RejectForRingBondRatio\n";
#endif
  if (_max_ring_bond_ratio > 1.0)
    return 0;

  (void) m.ring_membership();

  int nb = m.nedges();

  int ring_bonds = 0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);
    if (b->nrings())
      ring_bonds++;
  }

  if (Fraction(ring_bonds, nb) >= _max_ring_bond_ratio) {
    counter.molecules_with_bad_max_ring_bond_ratios++;
    rejection_reason << "bad_ring_bond_ratio";
    return 1;     // return 1 means reject this molecule
  }

  return 0;       // return 0 means don't reject this molecule
}

int
MCFirstPass::RejectForRingSizeCondition(Molecule & m,
                               MCFirstPassCounter& counter,
                               IWString& rejection_reason) const {
#ifdef DEBUG_MC_FIRST_PASS
  cerr << "RejectForRingSizeCondition\n";
#endif
  for (int i = m.nrings() - 1; i >= 0; i--) {
    const Ring * ri = m.ringi(i);

    if (ri->number_elements() <= _upper_ring_size_cutoff) {
      continue;
    }

    counter.molecules_with_ring_sizes_out_of_range++;
    if (_verbose > 1)
      cerr << "Ring size out of range\n";
    rejection_reason << "ring size";
    return 1;
  }

  return 0;
}

int
MCFirstPass::ConvertIsotopes(Molecule & m, MCFirstPassCounter& counter) const {
  int rc = m.transform_to_non_isotopic_form();

  if (rc)
    counter.molecules_with_isotopes_converted++;

  return rc;
}

int
MCFirstPass::Report(const MCFirstPassCounter& counter,
                    std::ostream& output) const {
  output << counter.molecules_processed << " molecules processed, ";

  if (counter.molecules_processed == 0)
    return 1;

  output << counter.molecules_rejected << " molecules rejected\n";

  if (counter.molecules_containing_non_allowed_atom_types) {
    output << "Skipped " << counter.molecules_containing_non_allowed_atom_types <<
      " molecules containing non allowed atoms\n";
  }

  if (counter.molecules_containing_non_periodic_table_elements) {
    output << "Skipped " << counter.molecules_containing_non_periodic_table_elements <<
      " molecules containing non periodic table atoms\n";
  }

  if (counter.molecules_containing_colvalent_non_organics) {
    output << "Skipped " << counter.molecules_containing_colvalent_non_organics <<
      " molecules containing covalently bonded non organics\n";
  }

  if (counter.molecules_with_no_interesting_atoms) {
    output << "Skipped " << counter.molecules_with_no_interesting_atoms <<
            " molecules with no interesting atoms\n";
  }

  if (counter.molecules_below_atom_count_cutoff) {
    output << "Skipped " << counter.molecules_below_atom_count_cutoff << 
            " molecules with atom count less than " << _lower_atom_count_cutoff << '\n';
  }

  if (counter.molecules_above_atom_count_cutoff) {
    output << "Skipped " << counter.molecules_above_atom_count_cutoff << 
            " molecules with atom count greater than " << _upper_atom_count_cutoff << '\n';
  }

  if (counter.molecules_with_too_few_rings) {
    output << "Skipped " << counter.molecules_with_too_few_rings << 
            " molecules having fewer than " << _lower_ring_count_cutoff << " rings\n";
  }

  if (counter.molecules_with_too_many_rings) {
    output << "Skipped " << counter.molecules_with_too_many_rings << 
            " molecules having more than " << _upper_ring_count_cutoff << " rings\n";
  }

  if (counter.molecules_with_ring_sizes_out_of_range) {
    output << "Skipped " << counter.molecules_with_ring_sizes_out_of_range << " molecules with rings containing more than " << _upper_ring_size_cutoff << " atoms\n";
  }
  
  if (counter.molecules_containing_isotopes) {
    output << counter.molecules_containing_isotopes << " molecules contained isotopic atoms\n";
  }

  if (counter.molecules_containing_isotopes) {
    output << counter.molecules_containing_isotopes << " molecules containing isotopes\n";
  }

  if (counter.molecules_with_isotopes_converted) {
    output << counter.molecules_with_isotopes_converted << " isotopes converted to natural form\n";
  }

  _elements_to_remove.report(output);

  if (_element_transformations.number_elements()) {
    _element_transformations.debug_print(output);
  }

  if (counter.molecules_with_abnormal_valences) {
    output << counter.molecules_with_abnormal_valences << " molecules containing abnormal valences\n";
  }

  if (counter.molecules_with_two_connected_hydrogens) {
    output << counter.molecules_with_two_connected_hydrogens << " molecules containing two connected hydrogens\n";
  }

  if (counter.molecules_with_bad_max_ring_bond_ratios) {
    output << counter.molecules_with_bad_max_ring_bond_ratios << " molecules with ring bond ratios out of range\n";
  }

  return output.good();
}

}  // namespace lilly_medchem_rules
