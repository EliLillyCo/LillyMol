/*
  Computes chirality fingerprints
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/extended_connectivity_fp.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString tag("NCCHR<");

static int function_as_tdt_filter = 0;

static int discard_existing_chirality = 0;

static int chirality_from_3d = 0;

static int remove_invalid_chiral_centres = 0;

static int skip_non_chiral_molecules = 0;

static int non_chiral_molecules_skipped = 0;

/*
  We differentiate between the number of chiral centres in the molecule, and the
  number we process.
*/

static extending_resizable_array<int> chiral_centre_count;

static extending_resizable_array<int> chiral_centres_processed;

static extending_resizable_array<int> expansions_done;

static resizable_array_p<Substructure_Query> only_process_atoms;

static int max_iterations = std::numeric_limits<int>::max();

static Atom_Typing_Specification atom_typing_specification;

static int include_EC_bits_central_atom = 0;

static int fingerprint_individual_layers = 0;

static int replicates = 1;

/*
  We can work as a tester by bringing in random variants of the same
  molecule and computing the fingerprint
*/

static int ntest = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Computes chirality fingerprints\n";
  cerr << "  -J <tag>      tag to procude\n";
  cerr << "  -P <atype>    atom typing to use\n";
  cerr << "  -s <smarts>   only process atoms matching <smarts>\n";
  cerr << "  -q <query>    only process atoms matching <query>\n";
  cerr << "  -f            function as a TDT filter\n";
  cerr << "  -u            remove invalid chiral centres\n";
  cerr << "  -3            discern chirality from 3D structure\n";
  cerr << "  -x            discard all existing chirality information\n";
  cerr << "  -n            discard molecules containing no chirality on input\n";
  cerr << "  -e <rep>      produce EC type bits for central atom\n";
  cerr << "  -r <rep>      how many bit replicates to generate\n";
  cerr << "  -y <rep>      fingerprint each layer as the chiral centre is discerned\n";
  cerr << "  -T <ntest>    work as a tester\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

/*
  For each atom around the chiral centre, we need to keep track of the expansion

  We keep track of the bonds that define our next layer.

*/

class Atom_Layer : public resizable_array<const Bond*> {
 private:
 public:
  void
  add_bond(const Bond* b) {
    this->add(b);
  }
};

class Atom_Expansion {
 private:
  atom_number_t _a;  // bonded to the chiral atom
  const int _matoms;
  int _layer;
  int _stopped;

  resizable_array_p<Atom_Layer> _layers;

  int* _atom_in_layer;

  uint64_t _hash;

  //  When we make the fingerprint, we need to be sorted by position in the chiral centre

  int _position;

 public:
  //  Atom_Expansion(const atom_number_t, const atom_number_t, const int matoms, const int
  //  * atype);
  Atom_Expansion(const atom_number_t, const Bond* b, const int matoms,
                 const unsigned int* atype);
  ~Atom_Expansion();

  int debug_print(std::ostream&) const;

  int
  stopped() const {
    return _stopped;
  }

  uint64_t
  hash() const {
    return _hash;
  }

  void
  set_hash(const uint64_t s) {
    _hash = s;
  }

  atom_number_t
  start_atom() const {
    return _a;
  }

  void
  set_position(const int s) {
    _position = s;
  }

  int
  position() const {
    return _position;
  }

  int expand(const Molecule& m, const unsigned int* atype);

  void
  set_stopped(const int s) {
    _stopped = s;
  }

  int
  atom_in_layer(const atom_number_t a) const {
    return _atom_in_layer[a];
  }

  int
  nlayers() const {
    return _layers.number_elements();
  }

  const Atom_Layer*
  atom_layer(const int s) const {
    return _layers[s];
  }
};

#ifdef OLLDDDDD
Atom_Expansion::Atom_Expansion(const atom_number_t chiral_atom, const atom_number_t a1,
                               const int matoms, const int* atype)
    : _a(a1), _matoms(matoms) {
  _atom_in_layer = new_int(matoms, -1);
  _atom_in_layer[chiral_atom] = 0;

  if (a1 >= 0) {
    _atom_in_layer[a1] = 1;
  }

  _layer = 1;
  _stopped = 0;
  if (a1 >= 0) {
    _hash = atype[a1];
  } else {
    _hash = 0;
  }

  _position = 0;
}
#endif

Atom_Expansion::Atom_Expansion(const atom_number_t chiral_atom, const Bond* b,
                               const int matoms, const unsigned int* atype)
    : _matoms(matoms) {
  _layer = 1;

  _atom_in_layer = new_int(matoms, -1);
  _atom_in_layer[chiral_atom] = 0;

  if (NULL == b) {
    _hash = 0;
    _stopped = 1;
    _a = INVALID_ATOM_NUMBER;
  } else {
    _stopped = 0;
    _a = b->other(chiral_atom);
    _atom_in_layer[_a] = 1;
    _hash = atype[_a];
    _layers.add(new Atom_Layer());
    _layers.last_item()->add_bond(b);
  }

  _position = 0;
}

Atom_Expansion::~Atom_Expansion() {
  delete[] _atom_in_layer;

  return;
}

int
Atom_Expansion::debug_print(std::ostream& output) const {
  output << "Atom_Expansion::debug_print:atom " << _a << " _layer " << _layer
         << " stopped " << _stopped << " nlayers " << _layers.number_elements()
         << " hash " << _hash << '\n';
  for (int i = 0; i < _matoms; ++i) {
    if (_atom_in_layer[i] >= 0) {
      output << " atom " << i << " in layer " << _atom_in_layer[i] << '\n';
    }
  }

  return 1;
}

static uint64_t
bond_hash(const Bond& b) {
  if (b.is_aromatic()) {
    return 4;
  }

  if (b.is_single_bond()) {
    return 1;
  }
  if (b.is_double_bond()) {
    return 2;
  }
  if (b.is_triple_bond()) {
    return 3;
  }

  return 5;
}

class Atom_Expansion_Sorter {
 public:
  int
  operator()(const Atom_Expansion* a1, const Atom_Expansion* a2) {
    const int nl1 = a1->nlayers();
    const int nl2 = a2->nlayers();
    if (nl1 < nl2) {
      return 1;
    }

    if (nl1 > nl2) {
      return -1;
    }

    if (a1->hash() < a2->hash()) {
      return 1;
    }

    if (a1->hash() > a2->hash()) {
      return -1;
    }

    return 0;
  }
};

int
Atom_Expansion::expand(const Molecule& m, const unsigned int* atype) {
  _layers.add(new Atom_Layer);

  int rc = 0;

  uint64_t h = 0;

  for (int i = 0; i < _matoms; ++i) {
    if (_atom_in_layer[i] != _layer) {
      continue;
    }

    const Atom* a = m.atomi(i);

    const int acon = a->ncon();

    for (int j = 0; j < acon; ++j) {
      const Bond* b = a->item(j);

      const atom_number_t k = b->other(i);

      if (_layer + 1 == _atom_in_layer[k]) {  // forming a ring
        ;
      } else if (_atom_in_layer[k] >= 0) {
        continue;
      }

      _atom_in_layer[k] = _layer + 1;

      _layers.last_item()->add_bond(b);

      h += atype[i] * bond_hash(*b) * atype[k];

      //    cerr << "from " << i << ' ' << m.smarts_equivalent_for_atom(i) << " to " << k
      //    << ' ' << m.smarts_equivalent_for_atom(k) << '\n';
      rc++;
    }
  }

  _layer++;

  _hash = _hash * _layer + h;

  if (0 == rc) {
    _stopped = 1;
  }

  return rc;
}

static int
not_chiral(Molecule& m, const Chiral_Centre&, const unsigned int* atype,
           Sparse_Fingerprint_Creator& sfc) {
  return 1;
}

static int
resolved(const resizable_array_p<Atom_Expansion>& expansion) {
  const int n = expansion.number_elements();

  const Atom_Expansion* prev = expansion[0];

  for (int i = 1; i < n; ++i) {
    const auto* curr = expansion[i];
    if (prev->hash() == curr->hash()) {
      return 0;
    }

    prev = curr;
  }

  return 1;
}

template void resizable_array_base<Atom_Expansion*>::iwqsort<Atom_Expansion_Sorter>(
    Atom_Expansion_Sorter&);

// #define DEBUG_CHIRALITY_FINGERPRINT

static int
chirality_fingerprint(Molecule& m, const unsigned int* atype, const int* include_atom,
                      const Chiral_Centre& c, Sparse_Fingerprint_Creator& sfc) {
  const int matoms = m.natoms();

  const auto zatom = c.a();

  const Atom* a = m.atomi(zatom);

  const int acon = a->ncon();

  Atom_Expansion_Sorter aes;

  resizable_array_p<Atom_Expansion> expansion;

  for (int i = 0; i < acon; ++i) {
    const Bond* b = a->item(i);

    expansion.add(new Atom_Expansion(zatom, b, matoms, atype));
  }

  if (3 == acon) {
    int ok_to_process = 0;
    if (1 == m.implicit_hydrogens(zatom)) {
      ok_to_process = 1;
    } else if (16 == a->atomic_number()) {  // we ignore 3 connected S atoms
      return 1;
    } else {
      int lp;
      m.lone_pair_count(zatom, lp);
      if (1 == lp) {
        ok_to_process = 1;
      }
    }

    if (!ok_to_process) {
      cerr << "chirality_fingerprint:three connected atom, but no imp H or LP "
           << m.name() << " atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom)
           << '\n';
      return 0;
    }

    Atom_Expansion* x = new Atom_Expansion(zatom, NULL, matoms, atype);
    x->set_stopped(1);
    x->set_hash(7);  // just some arbitrary, but small number

    expansion.add(x);
  }

  assert(4 == expansion.number_elements());

  expansion.iwqsort(aes);

  if (!resolved(expansion)) {
    for (int i = 0; i < max_iterations; ++i) {
      int expanded = 0;

      for (int j = 0; j < acon; ++j) {
        Atom_Expansion* a = expansion[j];

        if (a->stopped()) {
          continue;
        }

        if (a->expand(m, atype)) {
          expanded++;
        }
      }

#ifdef DEBUG_CHIRALITY_FINGERPRINT
      cerr << "expanded " << expanded << '\n';
#endif

      if (0 == expanded) {
        return not_chiral(m, c, atype, sfc);
      }

      expansion.iwqsort(aes);

      if (resolved(expansion)) {
        break;
      }
    }
  }

#ifdef DEBUG_CHIRALITY_FINGERPRINT
  cerr << "Resolved around atom " << zatom << '\n';
  for (int i = 0; i < 4; ++i) {
    const auto* ps = expansion[i];
    if (INVALID_ATOM_NUMBER != ps->start_atom()) {
      cerr << " atom " << ps->start_atom() << ' '
           << m.smarts_equivalent_for_atom(ps->start_atom()) << " hash " << ps->hash()
           << '\n';
    } else {
      cerr << " H\n";
    }
  }
  c.debug_print(cerr);
#endif

  int max_steps = 0;

  for (auto* psi : expansion) {
    if (psi->nlayers() > max_steps) {
      max_steps = psi->nlayers();
    }
  }

  expansions_done[max_steps]++;

  if (include_EC_bits_central_atom) {
    EC_Fingerprint_Generator ecfpg;

    ecfpg.set_max_radius(max_steps);

    resizable_array<unsigned int> bits;

    ecfpg.generate_fingerprint(m, zatom, atype, include_atom, 1, bits);

    for (auto b : bits) {
      //    cerr << "Shell expansion setting bit " << b << '\n';
      sfc.hit_bit(b, include_EC_bits_central_atom);
    }
  }

  // Now something that is sensitive to the order

  sfc.hit_bit(81923 * max_steps, 5);

  const auto x0 = expansion[0]->start_atom();
  const auto x1 = expansion[1]->start_atom();
  const auto x2 = expansion[2]->start_atom();
  // never check the 4th item because it may be Hydrogen
  if (x0 == c.top_front()) {
    if (x1 == c.top_back()) {
      if (x2 == c.right_down()) {
        expansion.swap_elements(2, 3);
      }
    } else if (x1 == c.left_down()) {
      if (x2 == c.top_back()) {
        expansion.swap_elements(2, 3);
      }
    } else  // right down
    {
      if (x2 == c.left_down()) {
        expansion.swap_elements(2, 3);
      }
    }
  } else if (x0 == c.top_back()) {
    if (x1 == c.top_front()) {
      if (x2 == c.left_down()) {
        expansion.swap_elements(2, 3);
      }
    } else if (x1 == c.left_down()) {
      if (x2 == c.right_down()) {
        expansion.swap_elements(2, 3);
      }
    } else  // right down
    {
      if (x2 == c.top_front()) {
        expansion.swap_elements(2, 3);
      }
    }
  } else if (x0 == c.left_down()) {
    if (x1 == c.top_front()) {
      if (x2 == c.right_down()) {
        expansion.swap_elements(2, 3);
      }
    } else if (x1 == c.top_back()) {
      if (x2 == c.top_front()) {
        expansion.swap_elements(2, 3);
      }
    } else  // right down
    {
      if (x2 == c.top_back()) {
        expansion.swap_elements(2, 3);
      }
    }
  } else if (x0 == c.right_down()) {
    if (x1 == c.top_front()) {
      if (x2 == c.top_back()) {
        expansion.swap_elements(2, 3);
      }
    } else if (x1 == c.top_back()) {
      if (x2 == c.left_down()) {
        expansion.swap_elements(2, 3);
      }
    } else  // left down
    {
      if (x2 == c.top_front()) {
        expansion.swap_elements(2, 3);
      }
    }
  }

#ifdef DEBUG_CHIRALITY_FINGERPRINT
  cerr << "After uniqueness\n";
  for (int i = 0; i < 4; ++i) {
    const auto* ps = expansion[i];
    if (INVALID_ATOM_NUMBER != ps->start_atom()) {
      cerr << " atom " << ps->start_atom() << ' '
           << m.smarts_equivalent_for_atom(ps->start_atom()) << " hash " << ps->hash()
           << '\n';
    } else {
      cerr << " H\n";
    }

    //  ps->debug_print(cerr);
  }

  cerr << "Computing fingerprint, max_steps " << max_steps << '\n';
#endif

  uint64_t h = atype[zatom];

  for (int layer = 1; layer <= max_steps; ++layer)  // layer 0 is the central atom
  {
    uint64_t hlayer = 0;

    for (int j = 0; j < 4; ++j) {
      const auto* psj = expansion[j];

#ifdef DEBUG_CHIRALITY_FINGERPRINT
      cerr << "layer " << layer << " expansion[" << j << "] has " << psj->nlayers()
           << " layers\n";
#endif

      if (layer > psj->nlayers()) {
        continue;
      }

      const Atom_Layer* al =
          psj->atom_layer(layer - 1);  // there is no zero layer in the sets of bonds

      const int nb = al->number_elements();

#ifdef DEBUG_CHIRALITY_FINGERPRINT
      cerr << "Atom layer contains " << nb << " bonds\n";
#endif

      for (int k = 0; k < nb; ++k) {
        const Bond* b = al->item(k);

#ifdef DEBUG_CHIRALITY_FINGERPRINT
        cerr << " k = " << k << " bond " << b->a1() << ' '
             << m.smarts_equivalent_for_atom(b->a1()) << " (layer "
             << psj->atom_in_layer(b->a1()) << ") " << b->a2() << ' '
             << m.smarts_equivalent_for_atom(b->a2()) << " (layer "
             << psj->atom_in_layer(b->a2()) << ")\n";
#endif

        if (layer == psj->atom_in_layer(b->a1())) {
          hlayer += (bond_hash(*b) * atype[b->a1()] * (j + 1));
        } else if (layer == psj->atom_in_layer(b->a2())) {
          hlayer += (bond_hash(*b) * atype[b->a2()] * (j + 1));
        } else {
          cerr << "HUH\n";
        }
      }
    }

    h = h * layer * 3149 + (layer * hlayer);

    if (fingerprint_individual_layers) {
      sfc.hit_bit(h, fingerprint_individual_layers);
    }
  }

#ifdef DEBUG_CHIRALITY_FINGERPRINT
  cerr << "Setting bit " << h << '\n';
#endif
  sfc.hit_bit(h, replicates);

  return 1;
}

static int
do_tests(Molecule& m, const int* include_atom, unsigned int* atype) {
  const int nchiral = m.chiral_centres();

  IWString fp;

  int rc = 1;

  for (int t = 0; t < ntest; ++t) {
    Sparse_Fingerprint_Creator sfc;

    const IWString& rsmi = m.random_smiles();
    Molecule x;
    if (!x.build_from_smiles(rsmi)) {
      return 0;
    }

    x.set_name(m.name());

    if (verbose > 2) {
      cerr << "Random variant " << t << " " << x.smiles() << '\n';
    }

    atom_typing_specification.assign_atom_types(x, atype);

    for (int i = 0; i < nchiral; i++) {
      const Chiral_Centre* c = x.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

      const atom_number_t a = c->a();

      if (0 == include_atom[a]) {
        continue;
      }

      chirality_fingerprint(x, atype, include_atom, *c, sfc);
    }

    IWString s;
    sfc.daylight_ascii_form_with_counts_encoded(tag, s);

    if (0 == t) {
      fp = s;
    } else if (fp != s) {
      cerr << "Fingerprint mismatch " << m.smiles() << ' ' << m.name() << '\n';
      cerr << "Random variant       " << x.smiles() << '\n';
      cerr << "Expected " << fp << '\n';
      cerr << "Got      " << s << '\n';
      rc = 0;
    }
  }

  return rc;
}

static int
chirality_fingerprint(Molecule& m, const int* process_atom,
                      IWString_and_File_Descriptor& output) {
  preprocess(m);

  molecules_read++;

  const int nchiral = m.chiral_centres();

  chiral_centre_count[nchiral]++;  // the number in the molecule, which might be different
                                   // from how many we process

  if (0 == nchiral && skip_non_chiral_molecules) {
    non_chiral_molecules_skipped++;

    return 1;
  }

  if (ntest || function_as_tdt_filter) {  // no need to write the tags
    ;
  } else  // write the tags
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  if (0 == nchiral)  // no chirality here
  {
    if (ntest) {
      ;
    } else if (!function_as_tdt_filter) {
      output << tag << ">\n";
      output << "|\n";
    }

    return 1;
  }

  const int matoms = m.natoms();

  unsigned int* atype = new unsigned int[matoms];
  std::unique_ptr<unsigned int[]> free_atype(atype);
  int* include_atom = new_int(matoms, 1);
  std::unique_ptr<int[]> free_include_atom(include_atom);

  atom_typing_specification.assign_atom_types(m, atype);

  if (ntest > 0) {
    return do_tests(m, include_atom, atype);
  }

  Sparse_Fingerprint_Creator sfc;

  int rc = 0;

  for (int i = 0; i < nchiral; i++) {
    const Chiral_Centre* c = m.chiral_centre_in_molecule_not_indexed_by_atom_number(i);

    const atom_number_t a = c->a();

    if (0 == process_atom[a]) {
      continue;
    }

    rc++;

    chirality_fingerprint(m, atype, include_atom, *c, sfc);
  }

  chiral_centres_processed[rc]++;  // the number we process, rather than the number in the
                                   // molecule

  IWString s;
  sfc.daylight_ascii_form_with_counts_encoded(tag, s);

  output << s << '\n';

  if (!function_as_tdt_filter) {
    output << "|\n";
  }

  return output.good();
}

static int
chirality_fingerprint(Molecule& m, IWString_and_File_Descriptor& output) {
  if (discard_existing_chirality) {
    m.remove_all_chiral_centres();
  }

  if (chirality_from_3d) {
    m.discern_chirality_from_3d_structure();
  }

  if (remove_invalid_chiral_centres) {
    lillymol::do_remove_invalid_chiral_centres(m);
  }

  const auto matoms = m.natoms();

  int* process_atom = new int[matoms];
  std::unique_ptr<int[]> free_process_atom(process_atom);

  if (0 == only_process_atoms.number_elements()) {
    std::fill_n(process_atom, matoms, 1);

    return chirality_fingerprint(m, process_atom, output);
  }

  std::fill_n(process_atom, matoms, 0);

  int queries_matching = 0;

  Molecule_to_Match target(&m);

  for (int i = 0; i < only_process_atoms.number_elements(); ++i) {
    Substructure_Results sresults;

    const auto nhits = only_process_atoms[i]->substructure_search(target, sresults);

    if (0 == nhits) {
      continue;
    }

    queries_matching++;

    sresults.each_embedding_set_vector(process_atom, 1);
  }

  if (0 == queries_matching) {
    cerr << "None of " << only_process_atoms.number_elements() << " queries matched "
         << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  return chirality_fingerprint(m, process_atom, output);
}

static int
chirality_fingerprint(data_source_and_type<Molecule>& input,
                      IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (NULL != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (!chirality_fingerprint(*m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
chirality_fingerprint_record(const_IWSubstring& buffer,
                             IWString_and_File_Descriptor& output) {
  assert(buffer.starts_with(smiles_tag));
  assert(buffer.ends_with('>'));

  buffer.remove_leading_chars(smiles_tag.length());
  buffer.chop();

  Molecule m;

  if (!m.build_from_smiles(buffer)) {
    cerr << "chirality_fingerprint_record:invalid smiles '" << buffer << "'\n";
    return 0;
  }

  preprocess(m);

  return chirality_fingerprint(m, output);
}

static int
chirality_fingerprint_filter(iwstring_data_source& input,
                             IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!chirality_fingerprint_record(buffer, output)) {
      cerr << "Fatal error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
chirality_fingerprint_filter(const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "chirality_fingerprint:cannot open '" << fname << "'\n";
    return 0;
  }

  return chirality_fingerprint_filter(input, output);
}

static int
chirality_fingerprint(const char* fname, FileType input_type,
                      IWString_and_File_Descriptor& output) {
  assert(NULL != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1) {
    input.set_verbose(1);
  }

  return chirality_fingerprint(input, output);
}

static int
chirality_fingerprint(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lJ:P:fx3uq:s:nT:e:r:y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('J')) {
    cl.value('J', tag);

    if (verbose) {
      cerr << "Fingerprints generated with tag '" << tag << "'\n";
    }

    if (!tag.ends_with('<')) {
      tag += '<';
    }
  }

  if (cl.option_present('x')) {
    discard_existing_chirality = 1;

    if (verbose) {
      cerr << "Will discard all existing chirality information\n";
    }
  }

  if (cl.option_present('u')) {
    remove_invalid_chiral_centres = 1;

    if (verbose) {
      cerr << "Will remove any invalid chirality information\n";
    }
  }

  if (cl.option_present('T')) {
    if (!cl.value('T', ntest) || ntest < 1) {
      cerr << "The number of test cycles to perform (-T) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will perform " << ntest << " random variant tests\n";
    }
  }

  if (cl.option_present('3')) {
    chirality_from_3d = 1;

    if (verbose) {
      cerr << "Will discern chirality from atom coordinates\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    if (!atom_typing_specification.build(p)) {
      cerr << "Invalid atom typing specification '" << p << "'\n";
      return 1;
    }
  } else {
    atom_typing_specification.build("UST:Z");
  }

  if (cl.option_present('n')) {
    skip_non_chiral_molecules = 1;

    if (verbose) {
      cerr << "Will skip non chiral molecules\n";
    }
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', include_EC_bits_central_atom) ||
        include_EC_bits_central_atom < 1) {
      cerr << "The number of central atom EC bit replicates (-e) must be a whole +ve "
              "number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will include EC type bits for central atom, "
           << include_EC_bits_central_atom << " bit replicates\n";
    }
  }

  if (cl.option_present('r')) {
    if (!cl.value('r', replicates) || replicates < 1) {
      cerr << "The number of fingerprint replicates (-r) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will produce " << replicates
           << " bit replicates of the chirality fingerprint\n";
    }
  }

  if (cl.option_present('y')) {
    if (!cl.value('y', fingerprint_individual_layers) ||
        fingerprint_individual_layers < 1) {
      cerr << "The fingerprint individual layers replicates (-y) option must be a whole "
              "+ve number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will produce " << fingerprint_individual_layers
           << " bit replicates of each layer\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (cl.option_present('f')) {
    if (cl.number_elements() > 1) {
      cerr << "Extra arguments ignored when working as filter\n";
      return 1;
    }

    function_as_tdt_filter = 1;

    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }

    if (!chirality_fingerprint_filter(cl[0], output)) {
      rc = 1;
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!chirality_fingerprint(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    for (int i = 0; i < chiral_centre_count.number_elements(); ++i) {
      if (chiral_centre_count[i] > 0) {
        cerr << chiral_centre_count[i] << " molecules had " << i << " chiral centres\n";
      }
    }

    for (int i = 0; i < expansions_done.number_elements(); ++i) {
      if (expansions_done[i]) {
        cerr << i << " expansions done " << expansions_done[i] << " times\n";
      }
    }

    if (skip_non_chiral_molecules) {
      cerr << "skipped " << non_chiral_molecules_skipped << " non chiral molecules\n";
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = chirality_fingerprint(argc, argv);

  return rc;
}
