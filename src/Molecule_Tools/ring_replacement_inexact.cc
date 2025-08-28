// Ring replacement where the replacement ring is not necessarunset_isotopes();
// an exact replacement for the ring being lost.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_set.h"

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#ifdef BUILD_BAZEL
#include "Molecule_Lib/reaction.pb.h"
#include "Molecule_Tools/replacement_ring.pb.h"
#else
#include "reaction.pb.h"
#include "replacement_ring.pb.h"
#endif

namespace ring_replacement_inexact {

using std::cerr;
using RplRing::ReplacementRing;
using combinations::Combinations;

constexpr int kYtterbium = 39;

// Can be used during debugging, change name...
#define REMOVE_BOND_BETWEEN_ATOMS_DEBUG(m, a1, a2) \
  { \
    if (! m.are_bonded(a1, a2)) { \
      cerr << m.smiles() << ' ' << m.name() << '\n'; \
      write_isotopically_labelled_smiles(m, false, cerr); \
      cerr << " atoms " << a2 << " and " << a2 << '\n'; \
    } else { \
      m.remove_bond_between_atoms(a1, a2); \
    } \
  }

#define REMOVE_BOND_BETWEEN_ATOMS(m, a1, a2) m.remove_bond_between_atoms(a1, a2)

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(
Ring replacement where ring sizes can change and substitution patterns are not necessarily
preserved.
Typical usage

ring_replacement_inexact -R rings_6a6a.smi -s '[/IWfss2r5r5]' -c file.smi

which processes two fused ring systems with two aromatic, five membered rings and
replaces those ring systems with fused 6 membered aromatics instead.

With some common, useful options
ring_replacement_inexact -z i -p -c -v -R rings_6a6a.smi -s '[/IWfss2r5r5]' file.smi

 -R <fname>     one or more RplRing::ReplacementRing textproto files, such as what is generated
                by ring_extraction.
                Use  '-R F:file.txt'  for a file containing a list of replacement rings.
 -n <n>         minimum support level for replacement rings, the 'n:' value in the proto.
 -j <n>         discard replacements with more than <n> points of attachment.
 -s <smarts>    specify one or more ring atoms that define the ring(s) to be removed
                and replaced.
 -q <query>     query file specification of the atoms to be removed - same as the -s option.
 -z i           ignore molecules that do not match any of the -s/-q queries.
 -x <n>         max number of products per starting molecule. Arbitrary variants will be generated.
 -o             disallow ortho substitions of rings.
 -I .           remove the isotopic labels from product molecules. Second arg is mandatory, but not parsed.
 -p             write the parent molecule before writing the variants.
 -V             discard any product with an invalid valence.
 -e             preserve the same-ring attachment patterns of the starting molecule.
 -Y <query>     queries that product molecules must contain.
                tsubstructure -q syntax, so smarts is '-Y SMARTS:n'
                textproto query file is '-Y PROTO:/path/to/file/qry`
 -N <query>     queries for features that must NOT be in product molecules.
 -X ...         miscellaneous options, enter '-X help' for info.
 -c             remove chirality from input molecules.
 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

// Examine the atoms attached to `zatom` and if any of them are zero in `spinach`
// return that atom number.
atom_number_t
BondedToScaffoldAtom(const Molecule& m, atom_number_t zatom, const int* spinach) {
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);

    // spinach 0 means in the scaffold.
    if (spinach[o] == 0) {
      return o;
    }
  }

  return kInvalidAtomNumber;
}

// We need to describe substituents
class Substituent {
  private:
    // This will be an isotopically labelled molecule.
    Molecule _m;

    // The atom number of the attachment point in the starting molecule.
    atom_number_t _atom_number;

    // The atom in the starting molecule to which we are attached.
    // This will be an atom in the scaffold.
    atom_number_t _attached_to;

    // The first ring number to which this substituent is joined.
    int _ring_number;

    // The fused system identifier to which this atom is attached
    int _fused_system_identifier;

    int _is_fused;

    isotope_t _isotope;

  public:
    // Not initialised, rely on caller to initialise.
    // Maybe should have 3 argument constructor.
    Substituent() {};

    Molecule& molecule() {
      return _m;
    }

    atom_number_t atom_number() const {
      return _atom_number;
    }
    void set_atom_number(atom_number_t s) {
      _atom_number = s;
    }

    atom_number_t attached_to() const {
      return _attached_to;
    }

    void set_attachment_point(atom_number_t s) {
      _attached_to = s;
    }

    int ring_number() const {
      return _ring_number;
    }

    void set_ring_number(int s) {
      _ring_number = s;
    }

    void set_is_fused(int s) {
      _is_fused = s;
    }
    int is_fused() const {
      return _is_fused;
    }

    const IWString& smiles() {
      return _m.smiles();
    }

    int fused_system_identifier() const {
      return _fused_system_identifier;
    }
    void set_fused_system_identifier(int s) {
      _fused_system_identifier = s;
    }

    isotope_t isotope() const {
      return _isotope;
    }
    void set_isotope(isotope_t s) {
      _isotope = s;
    }

    // `m` is assumed to be the molecule from which we were formed.
    // Break the bond joining this substituent to the ring system.
    int BreakBond(Molecule& m) const;
};

int
Substituent::BreakBond(Molecule& m) const {
#ifdef CHECK_BONDS_BEFORE_BREAKING
  if (! m.are_bonded(_atom_number, _attached_to)) {
    cerr << "Substituent::BreakBond:not bonded\n";
    return 0;
  }
#endif

  return m.remove_bond_between_atoms(_atom_number, _attached_to);
}

// There can be multiple Substituents associated with a given ring.

class SubstituentsForRing : public resizable_array_p<Substituent> {
  private:
    // The ring number.
    const int _ring_number;

    // Each substituent has applied an isotope to an atom in the ring
    // system and an atom in the substituent. For efficiency, we store
    // a list of the isotopes that are to be broken in order to separate
    // the substituents from the ring system.
    extending_resizable_array<isotope_t> _bonds_to_break;

  public:
    SubstituentsForRing(int r) : _ring_number(r) {
    }

    int AddSubstituent(Substituent* s);

    // For each Substituent, break the bonds joining the
    // substituent to the ring system.
    // Return the atom number of an atom in the ring
    // system - does not matter which one.
    atom_number_t BreakBonds(Molecule& m) const;

    int number_attachment_points() const {
      return _number_elements;
    }

    // We know the isotopic atoms used for attachment points,
    // look for those isotopes in `m` and return the corresponding
    // atom numbers.
    Set_of_Atoms AttachmentPoints(const Molecule& m) const;
    // Does the same thing, just different signature.
    int AttachmentPoints(const Molecule& m, Set_of_Atoms& result) const;
};

int
SubstituentsForRing::AddSubstituent(Substituent* s) {
  if (s->ring_number() != _ring_number) {
    return 0;
  }

  _bonds_to_break[s->isotope()] = 1;

  this->add(s);

  return 1;
}

// Each substituent breaks the bond it specifies.
// Return the atom number of an atom in the ring
// that will be removed.
atom_number_t
SubstituentsForRing::BreakBonds(Molecule& m) const {
  for (const Substituent* s : *this) {
    s->BreakBond(m);
  }

  return _things[0]->attached_to();
}

// And a ring system consists of multiple rings.
class SubstituentsForRingSystem : public 
                resizable_array_p<SubstituentsForRing> {
  private:

  public:
    int AddSubstituent(Substituent* s);

    int nrings() const {
      return _number_elements;
    }

    int number_attachment_points() const;

    int RemoveRingSystem(Molecule& m) const;

    // We know the isotopic atoms used for attachment points,
    // look for those isotopes in `m` and return the corresponding
    // atom numbers.
    Set_of_Atoms AttachmentPoints(const Molecule& m) const;
};

int
SubstituentsForRingSystem::AddSubstituent(Substituent* s) {
  int rnumber = s->ring_number();

  for (SubstituentsForRing* ring : *this) {
    if (ring->AddSubstituent(s)) {
      return 1;
    }
  }

  std::unique_ptr<SubstituentsForRing> new_ring = std::make_unique<SubstituentsForRing>(rnumber);
  new_ring->AddSubstituent(s);

  this->add(new_ring.release());

  return 1;
}

int
SubstituentsForRingSystem::number_attachment_points() const {
  int rc = 0;
  for (const SubstituentsForRing* s : *this) {
    rc += s->number_attachment_points();
  }

  return rc;
}

int
SubstituentsForRingSystem::RemoveRingSystem(Molecule& m) const {
  atom_number_t in_scaffold = kInvalidAtomNumber;

  for (const SubstituentsForRing* substituents : *this) {
    in_scaffold = substituents->BreakBonds(m);
  }

  if (in_scaffold == kInvalidAtomNumber) [[ unlikely ]] {
    return 0;
  }

  return m.remove_fragment_containing_atom(in_scaffold);
}

Set_of_Atoms
SubstituentsForRing::AttachmentPoints(const Molecule& m) const {
  Set_of_Atoms result;

  AttachmentPoints(m, result);

  return result;
}

int
SubstituentsForRing::AttachmentPoints(const Molecule& m,
                                      Set_of_Atoms& result) const {
  int nfound = 0;

  // The molecule has been formed by adding the replacement ring
  // to the existing substituents, so we need to find the last
  // occurrence of each of our isotopes.

  for (int i = m.natoms() - 1; i >= 0; --i) {
    const isotope_t iso = m.isotope(i);
    if (iso == 0) {
      continue;
    }

    if (_bonds_to_break[iso] == 0) {
      continue;
    }

    result << i;
    ++nfound;
    if (nfound == _number_elements) {
      return result.size();
    }
  }

  return result.size();
}

Set_of_Atoms
SubstituentsForRingSystem::AttachmentPoints(const Molecule& m) const {
  Set_of_Atoms result;

  for (const SubstituentsForRing* r : *this) {
    r->AttachmentPoints(m, result);
  }

  return result;
}

class Replacement {
  private:
    Molecule _m;

    uint32_t _count;

    IWString _id;

    uint32_t _nrings;

    Set_of_Atoms _attachment_points;

  public:
    // If `use_existing_attachment_points` is set and the molecule in the proto
    // has no isotopic atoms, that is a failure. But a harmless failure.
    // In that case, `failed_for_no_attachment_points` will be set, so the caller
    // knows not to complain.
    int Build(const ReplacementRing& proto, int use_existing_attachment_points,
              int& failed_for_no_attachment_points);

    const Molecule& molecule() const {
      return _m;
    }
    Molecule& molecule() {
      return _m;
    }

    uint32_t count() const {
      return _count;
    }

    uint32_t nrings() const {
      return _nrings;
    }

    const IWString& id() const {
      return _id;
    }

    const Set_of_Atoms& attachment_points() const {
      return _attachment_points;
    }

    Set_of_Atoms attachment_points() {
      return _attachment_points;
    }

    int number_attachment_points() const {
      return _attachment_points.number_elements();
    }
};

// If use_existing_attachment_points is set, the attachment
// points are restricted to the labelled atoms in the proto.
int
Replacement::Build(const ReplacementRing& proto,
                  int use_existing_attachment_points,
                  int& failed_for_no_attachment_points) {
  failed_for_no_attachment_points = 0;

  if (proto.smi().empty()) {
    cerr << "Replacement::Build:empty smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _m.build_from_smiles(proto.smi())) {
    cerr << "Replacement::Build:invalid smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  const int matoms = _m.natoms();

  if (use_existing_attachment_points) {
    for (int i = 0; i < matoms; ++i) {
      if (_m.isotope(i) > 0) {
        if (_m.ring_bond_count(i) == 0) {
          // Cases like [1CH2]=C1NC2=C(C=CC=C2)N1 Z13708695.5A6a
          // Silently ignore.
          // cerr << "Attached to non ring atom " << _m.smiles() << ' ' << proto.id() << '\n';
        } else {
          _attachment_points << i;
        }
      }
    }
  } else {  // anything with a H atom.
    for (int i = 0; i < matoms; ++i) {
      if (_m.hcount(i) == 0) {
        continue;
      }
      // Do not substitute things like =N exocyclic sites.
      if (_m.ring_bond_count(i) == 0) {
        continue;
      }

      _attachment_points << i;
    }
  }

  if (_attachment_points.empty()) {
    if (use_existing_attachment_points) {
      failed_for_no_attachment_points = 1;
    }
    return 0;
  }

  _m.unset_isotopes();

  // Attachment points must be sorted - so next_permutation works.
  if (_attachment_points.size() > 1) {
    std::sort(_attachment_points.begin(), _attachment_points.end(),
      []( atom_number_t a1, atom_number_t a2) {
        return a1 < a2;
      });
  }

  _count = proto.n();

  _id = proto.id();

  _nrings = _m.nrings();

  return 1;
}

class SameRing {
  private:
    resizable_array<Substituent*> _substituents;

  public:
};

// Build a reaction that has `ring_substitions` subsitions points
// in the ring, and `sidechains` to be attached sidechains.
int
BuildReaction(uint32_t ring_substitions, uint32_t sidechains,
              IWReaction& result) {
  ReactionProto::Reaction proto;

  IWString smarts;
  smarts << "[<" << (ring_substitions + 1) << ']';

  proto.mutable_scaffold()->set_id(0);
  proto.mutable_scaffold()->add_smarts(smarts.data(), smarts.length());

  IWString not_used;   // File name
  Sidechain_Match_Conditions smc;
  if (! result.ConstructFromProto(proto, not_used, smc)) {
    cerr << "BuildReaction:cannot build reaction from\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}
              


// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    resizable_array_p<Substructure_Query> _queries;

    // We only use replacement rings that satisfy a support requirement.
    uint32_t _min_support_requirement;

    resizable_array_p<Replacement> _replacements;

    // By default, we attach substituents to any open atom in the
    // replacement ring. If this is set, we only use the existing
    // attachment points.
    int _use_existing_attachment_points;

    isotope_t _isotope;

    int _write_parent_molecule;

    uint32_t _molecules_not_matching_queries;
    int _ignore_molecules_not_matching_queries;
    uint32_t _no_substituents;

    // If two substituents are attaced to the same ring in the
    // starting molecule, they must be attached to the same ring
    // in the product.
    // Note that this works even if we are compressing something like
    // a naphthalene to a benzene - all same ring relationships are
    // preserved, we never examine any different ring relationships.
    int _preserve_ring_substitutions;

    int _remove_isotopes_from_products;

    // We can limit the number of products generated per starting molecule.
    // This will be approximate since checking is not done at every step.
    uint32_t _max_products_per_starting_molecule;

    // We can disallow ortho replacements.
    int _allow_ortho_replacements;
    uint32_t _ortho_substitutions_suppressed;

    // Structural constraints on what gets generated.
    resizable_array_p<Substructure_Query> _products_must_contain;
    resizable_array_p<Substructure_Query> _products_cannot_contain;

    // Incremented when a molecule fails the _products_must_contain
    // or _products_cannot_contain constaints.
    uint32_t _discarded_for_query_mismatch;

    extending_resizable_array<uint32_t> _products_per_molecule;

    Chemical_Standardisation _chemical_standardisation;

    // The initial implementation did not handle the case of an
    // aliphatic ring with two connections on the same atom.
    // In order to keep things simple, we disconnect the two
    // substituents, insert a Y atom which is joined to both
    // the ring and the two substituents. In products molecules,
    // we reverse that transformation.
    int _two_substituents_ytterbium;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    // Do not write duplicates.
    // Duplicates can arise if there is symmetry in the system.
    // TODO:ianwatson Investigate if there is something we can do with the
    // symmetry of the replacement rings to try and cut down on the number of
    // redundant products generated. This might be hard...
    absl::flat_hash_set<IWString> _seen;

    int _discard_invalid_valence;
    uint32_t _invalid_valence_discarded;

    // The number of duplicates we detect with the _seen has.
    uint32_t _duplicates_discarded;

  // Private functions
    int ReadReplacementRings(IWString& fname, absl::flat_hash_set<IWString>& already_found);
    int ReadReplacementRings(iwstring_data_source& input, absl::flat_hash_set<IWString>& already_found);
    int ReadReplacementRings(const const_IWSubstring& buffer, absl::flat_hash_set<IWString>& already_found);
    int BuildReplacement(const const_IWSubstring& buffer,
                          absl::flat_hash_set<IWString>& already_found);
    int ReadFileOfReplacementRings(IWString& fname,
                iwstring_data_source& input,
                absl::flat_hash_set<IWString>& already_found);
    int ReadFileOfReplacementRings(IWString& fname,
                                    absl::flat_hash_set<IWString>& already_found);

    uint32_t ProcessInner(Molecule& m, IWString_and_File_Descriptor& output);

    int IdentifyMatchedAtoms(Molecule& m, int* matched_atoms);
    int OkSubstituent(Molecule& fragment);
    int MoleculeToSubstituents(Molecule& m,
           const int* matched_atoms,
           SubstituentsForRingSystem* substituents);
    int ReplaceSingleRing(Molecule& m,
                const SubstituentsForRing& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output);
    int ReplaceSingleRingInner(Molecule& m,
                int initial_matoms,
                const SubstituentsForRing& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output);
    int ReplaceSingleRingInner(Molecule& m,
                                int initial_matoms,
                                const Set_of_Atoms& substituent_atoms,
                                const Replacement& replacement,
                                const Set_of_Atoms& ra,
                                IWString_and_File_Descriptor& output);
    uint32_t ReplaceRingSystem(Molecule& m,
                const SubstituentsForRingSystem& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output);
    int ReplaceRingSystemSubstituentsAnywhere(Molecule& m,
                int initial_matoms,
                const SubstituentsForRingSystem& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output);
    int ReplaceRingSystemPreserveRingSubstitution(Molecule& m,
                int initial_matoms,
                const SubstituentsForRingSystem& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output);
    uint32_t Process(Molecule& m,
                 const SubstituentsForRingSystem& substituents,
                 IWString_and_File_Descriptor& output);
    uint32_t Process(Molecule& m,
                 const SubstituentsForRingSystem& substituents,
                 const Replacement& replacement,
                 IWString_and_File_Descriptor& output);
    int OkOrthoSubstituents(const Molecule& m, const Set_of_Atoms& ring_atoms, int n);
    int MaybeWrite(Molecule& m, const Replacement& replacement, IWString_and_File_Descriptor& output);

    int OkSubstructures(Molecule& m);
    int MatchesMustNotHaveQueries(Molecule& m);
    int MatchesMustHaveQueries(Molecule& m);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    uint32_t Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;

  _preserve_ring_substitutions = 0;

  _min_support_requirement = 0;

  _two_substituents_ytterbium = 1;

  _allow_ortho_replacements = 1;
  _ortho_substitutions_suppressed = 0;

  _max_products_per_starting_molecule = std::numeric_limits<uint32_t>::max();

  _remove_isotopes_from_products = 0;
  _invalid_valence_discarded = 1;

  _discarded_for_query_mismatch = 0;

  _write_parent_molecule = 0;

  _discard_invalid_valence = 0;

  _use_existing_attachment_points = 1;

  _ignore_molecules_not_matching_queries = 0;
  _molecules_not_matching_queries = 0;
  _no_substituents = 0;

  _duplicates_discarded = 0;

  _isotope = 1;  // We need this to be non zero
}

void
DisplayDashXOptions(std::ostream& output) {
  output << "-X any             any available ring atom can receive substituents - not just isotopically labelled\n";

  ::exit(0);
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('X')) {
    IWString x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "any") {
        _use_existing_attachment_points = 0;
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }


  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (!q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }

      _queries << q.release();
    }
  }

  if (cl.option_present('q')) {
    static constexpr int kFlag = ' ';  // not used.
    if (!process_queries(cl, _queries, _verbose, kFlag)) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  if (_queries.empty()) {
    cerr << "No queries specified, specify ring atom matches via the -s and/or -q options\n";
    return 0;
  }

  if (cl.option_present('Y')) {
    if (! process_queries(cl, _products_must_contain, _verbose, 'Y')) {
      cerr << "Cannot read product must contain queries (-Y)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Defined " << _products_must_contain.size() << 
              " products must contain queries\n";
    }
  }

  if (cl.option_present('N')) {
    if (! process_queries(cl, _products_cannot_contain, _verbose, 'N')) {
      cerr << "Cannot read product must not contain queries (-N)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Defined " << _products_cannot_contain.size() << 
              " products must NOT contain queries\n";
    }
  }

  if (_verbose) {
    cerr << "Defined " << _queries.size() << " queries to define replacement rings\n";
  }

  if (cl.option_present('e')) {
    _preserve_ring_substitutions = 1;
    if (_verbose) {
      cerr << "Same ring attachment patterns will NOT be preserved\n";
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_products_per_starting_molecule)) {
      cerr << "Invalid max products per starting molecule (-x)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will stop enumeration of a molecule once " << _max_products_per_starting_molecule <<
              " new molecules have been generated\n";
    }
  }

  if (cl.option_present('o')) {
    _allow_ortho_replacements = 0;
    if (_verbose) {
      cerr << "Will not generate ortho substitued variants\n";
    }
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matching_queries = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching any query\n";
        }
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('p')) {
    _write_parent_molecule = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('I')) {
    _remove_isotopes_from_products = 1;
    if (_verbose) {
      cerr << "Will transform products to non isotopic variants\n";
    }
  }

  if (cl.option_present('V')) {
    _discard_invalid_valence = 1;
    if (_verbose) {
      cerr << "Will discard products containing invalid valences\n";
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _min_support_requirement)) {
      cerr << "The minimum support requirement for a replacement ring must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only consider replacement rings with at least " << 
              _min_support_requirement << " examples\n";
    }
  }

  if (! cl.option_present('R')) {
    cerr << "Must specify one or more replacement rings via the -R option\n";
    return 0;
  }

  if (cl.option_present('R')) {
    absl::flat_hash_set<IWString> already_found;
    IWString r;
    for (int i = 0; cl.value('R', r, i); ++i) {
      if (! ReadReplacementRings(r, already_found)) {
        cerr << "Cannot read replacement rings from '" << r << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Read " << _replacements.size() << " replacement rings\n";
    }
  }

  return 1;
}

int
Options::ReadReplacementRings(IWString& fname,
                              absl::flat_hash_set<IWString>& already_found) {
  if (fname.starts_with("F:")) {
    fname.remove_leading_chars(2);
    return ReadFileOfReplacementRings(fname, already_found);
  }

  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadReplacementRings:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReplacementRings(input, already_found);
}

int
Options::ReadFileOfReplacementRings(IWString& fname,
                                    absl::flat_hash_set<IWString>& already_found) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFileOfReplacementRings:cannot open '" << fname << "'\n";
    return 0;
  }

  // The directory...
  if (fname.contains('/')) {
    fname.truncate_at_last('/');
  }

  return ReadFileOfReplacementRings(fname, input, already_found);
}

int
Options::ReadFileOfReplacementRings(IWString& dirname,
                iwstring_data_source& input,
                absl::flat_hash_set<IWString>& already_found) {
  IWString buffer;
  while (input.next_record(buffer)) {
    IWString fname;
    if (buffer.starts_with('/')) {
      fname = buffer;
    } else {
      fname << dirname << '/' << buffer;
    }
    if (! ReadReplacementRings(fname, already_found)) {
      cerr << "Options::ReadFileOfReplacementRings:cannot read '" << fname << "'\n";
      return 0;
    }
  }

  return _replacements.size();
}

int
Options::ReadReplacementRings(iwstring_data_source& input,
                              absl::flat_hash_set<IWString>& already_found) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (! BuildReplacement(buffer, already_found)) {
      cerr << "Options::ReadReplacementRings:invalid input\n";
      cerr << buffer << '\n';
      return 0;
    }
  }
  return 1;
}

int
Options::BuildReplacement(const const_IWSubstring& buffer,
                          absl::flat_hash_set<IWString>& already_found) {

  google::protobuf::io::ArrayInputStream zero_copy_array(buffer.data(), buffer.nchars());
  ReplacementRing proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "Options:ReadReplacementRings:cannot parse proto " << buffer << '\n';
    return 0;
  }

  // Not processing these right now - why not?
  if (proto.has_exo()) {
    return 1;
  }

  if (proto.n() < _min_support_requirement) {
    return 1;
  }

  std::unique_ptr<Replacement> r = std::make_unique<Replacement>();

  int failed_for_no_attachment_points = 0;
  if (r->Build(proto, _use_existing_attachment_points, failed_for_no_attachment_points)) {
    // great
  } else if (failed_for_no_attachment_points) {
    // safe to ignore.
    return 1;
  } else {
    cerr << "Options::ReadReplacementRings:invalid proto " << proto.ShortDebugString() << '\n';
    cerr << "ignored\n";
    return 1;
  }

  if (_use_existing_attachment_points) {
    const IWString& usmi = r->molecule().unique_smiles();
    if (auto iter = already_found.find(usmi); iter != already_found.end()) {
      return 1;
    }

    already_found.insert(usmi);
    r->molecule().invalidate_smiles();
  }

  _replacements << r.release();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _molecules_not_matching_queries << " molecules did not match any queries\n";
  output << _no_substituents << " molecuels had no substients at matched atoms\n";
  if (! _allow_ortho_replacements) {
    output << _ortho_substitutions_suppressed << " ortho substited products suppressed\n";
  }
  if (_discard_invalid_valence) {
    output << _invalid_valence_discarded << " products with invalid valences discarded\n";
  }
  if (_products_must_contain.size() || _products_cannot_contain.size() > 0) {
    output << _discarded_for_query_mismatch << " products discarded via -Y or -N queries\n";
  }

  for (int i = 0; i < _products_per_molecule.number_elements(); ++i) {
    if (_products_per_molecule[i] == 0) {
      continue;
    }
    output << _products_per_molecule[i] << " molecules generated " << i << " products\n";
  }

  output << "Generated " << _seen.size() << " molecules\n";
  output << _duplicates_discarded << " duplicates discarded\n";

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    m.revert_all_directional_bonds_to_non_directional();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

// For now, there are no constraints on substituents.
int
Options::OkSubstituent(Molecule& fragment) {
  return 1;
  // const int matoms = fragment.natoms();

  // const int nrings = fragment.nrings();

  return 1;
}

// If `zatom` is bonded to an unmatched atom, return that atom.
atom_number_t
BondedToUnmatchedAtom(const Molecule& m,
                      atom_number_t zatom,
                      const int* matched) {
  for (const Bond* b : m[zatom]) {
    if (b->nrings()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (! matched[o]) {
      return o;
    }
  }

  return kInvalidAtomNumber;
}

// `matched_atoms` has been labelled by ring system number.
int
Options::MoleculeToSubstituents(Molecule& m,
           const int* matched_atoms,
           SubstituentsForRingSystem* substituents) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> down_the_bond = std::make_unique<int[]>(matoms);

  // We put unique isotopes on the atoms defining a substituent.
  // One on the ring atom, one on the first atom in the substituent.
  isotope_t iso = 1;

  // Look for matched atoms with substituents.
  for (int i = 0; i < matoms; ++i) {
    if (! matched_atoms[i]) {
      continue;
    }

    if (m.ncon(i) <= 2) {  // No substituents here.
      continue;
    }

    atom_number_t unmatched = BondedToUnmatchedAtom(m, i, matched_atoms);
    if (unmatched == kInvalidAtomNumber) {
      continue;
    }

    std::fill_n(down_the_bond.get(), matoms, 0);
    auto dtb = m.DownTheBond(i, unmatched, down_the_bond.get());
    if (! dtb) {
      continue;
    }

    std::unique_ptr<Substituent> substituent = std::make_unique<Substituent>();
    m.create_subset(substituent->molecule(), down_the_bond.get());

    if (! OkSubstituent(substituent->molecule())) {
      continue;
    }

    substituent->set_atom_number(unmatched);
    substituent->set_attachment_point(i);

    m.set_isotope(i, iso);
    m.set_isotope(unmatched, iso);
    substituent->set_isotope(iso);

    ++iso;

    const Ring* r = m.ring_containing_atom(i);
    assert(r != nullptr);

    substituent->set_ring_number(r->ring_number());
    // substituent->set_fused_system_identifier(matched[i]);

    if (r->is_fused()) {
      substituent->set_is_fused(1);
    } else {
      substituent->set_is_fused(0);
    }

    substituents[matched_atoms[i]].AddSubstituent(substituent.release());
  }

  return 1;
}


// Look for the first match to `_queries` in `m` and mark all matched
// atoms in `matched_atoms`.
// Return the number of queries matching.
// Unclear if we should match the first query or maybe all of them...
int
Options::IdentifyMatchedAtoms(Molecule& m,
                              int* matched_atoms) {
  Molecule_to_Match target(&m);

  for (Substructure_Query* q : _queries) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    sresults.each_embedding_set_vector(matched_atoms, 1);
    return 1;
  }

  return 0;
}

// Unfortunately we don't have a ready means of detecting
// fatal errors, we always return 1.
uint32_t
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  if (m.nrings() == 0) {
    return 1;
  }

  const uint32_t rc = ProcessInner(m, output);

  ++_products_per_molecule[rc];

  return 1;
}

int
IncludesIsotopicAtoms(const Molecule& m, const Set_of_Atoms& r) {
  for (atom_number_t a : r) {
    if (m.isotope(a)) {
      return 1;
    }
  }

  return 0;
}


int
TwoSubstituentsInsertYtterbium(Molecule& m, const int* matched_atoms) {
  const int matoms = m.natoms();

  Set_of_Atoms in_ring;
  Set_of_Atoms in_substituent;

  Set_of_Atoms connections;  // scope here for efficiency.
  connections.reserve(2);

  for (int i = 0; i < matoms; ++i) {
    if (! matched_atoms[i]) {
      continue;
    }

    const Atom& a = m[i];

    if (a.atomic_number() != 6) {
      continue;
    }
    if (a.ncon() != 4) {
      continue;
    }
    if (m.ring_bond_count(i) == 0) [[ unlikely ]] {  // can this happen?
      continue;
    }

    connections.resize_keep_storage(0);
    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      if (matched_atoms[o]) {
        continue;
      }
      if (m.ring_bond_count(o)) {
        continue;
      }

      connections << o;
    }

    if (connections.size() != 2) {
      continue;
    }

    in_ring << i;
    in_substituent += connections;
  }

  if (in_substituent.empty()) {
    return 0;
  }

  // cerr << "in ring " << in_ring << " in_substituent " << in_substituent << '\n';
  assert(in_ring.size() * 2 == in_substituent.size());

  static const Element* ytterbium = get_element_from_atomic_number(kYtterbium);

  for (uint32_t i = 0; i < in_ring.size(); ++i) {
    atom_number_t r = in_ring[i];
    atom_number_t s1 = in_substituent[i + i];
    atom_number_t s2 = in_substituent[i + i + 1];
    //m.remove_bond_between_atoms(r, s1);
    REMOVE_BOND_BETWEEN_ATOMS(m, r, s1);
    //m.remove_bond_between_atoms(r, s2);
    REMOVE_BOND_BETWEEN_ATOMS(m, r, s2);

    int natoms = m.natoms();
    m.add(ytterbium);
    m.add_bond(natoms, r, SINGLE_BOND);
    m.add_bond(natoms, s1, SINGLE_BOND);
    m.add_bond(natoms, s2, SINGLE_BOND);
  }

  // cerr << m.smiles() << " ytterbium\n";
  return in_ring.number_elements();
}

// Note that this will also follow spiro fusions.
void
ExtendMatchedAtomsAcrossFusions(Molecule& m,
                atom_number_t zatom,
                int* matched_atoms,
                int flag,
                int* visited) {
  visited[zatom] = 1;
  matched_atoms[zatom] = flag;

  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (visited[o]) {
      continue;
    }

    // =O type spinach atom.
    if (b->is_double_bond() && m.ncon(o) == 1) {
      matched_atoms[o] = flag;
      continue;
    }

    if (! b->nrings()) {
      continue;
    }

    ExtendMatchedAtomsAcrossFusions(m, o, matched_atoms, flag, visited);
  }
}

int
ExtendMatchedAtomsAcrossFusions(Molecule& m,
                int* matched_atoms) {
  const int matoms = m.natoms();

  m.ring_membership();  // Force sssr

  std::unique_ptr<int[]> visited = std::make_unique<int[]>(matoms);
  std::fill_n(visited.get(), matoms, 0);

  int system_number = 0;
  for (int i = 0; i < matoms; ++i) {
    if (visited[i]) {
      continue;
    }
    if (matched_atoms[i] == 0) {
      continue;
    }
    if (m.ring_bond_count(i) == 0) {
      continue;
    }

    ++system_number;

    ExtendMatchedAtomsAcrossFusions(m, i, matched_atoms, system_number, visited.get());
  }

  return system_number;
}

// If we have added Ytterbium atoms to the molecule, we need to adjust
// the size of the matched atoms array.
void
ResizeMatchedAtomsArray(std::unique_ptr<int[]>& matched_atoms,
                        int old_size, int new_size) {
  assert(old_size < new_size);

  std::unique_ptr<int[]> new_array = std::make_unique<int[]>(new_size);
  std::copy_n(matched_atoms.get(), old_size, new_array.get());

  // Fill the new array members with zero.
  for (int i = old_size; i < new_size; ++i) {
    new_array[i] = 0;
  }

  matched_atoms.reset(new_array.release());
}

uint32_t
Options::ProcessInner(Molecule& m, IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> matched_atoms = std::make_unique<int[]>(matoms);
  std::fill_n(matched_atoms.get(), matoms, 0);

  if (! IdentifyMatchedAtoms(m, matched_atoms.get())) {
    ++_molecules_not_matching_queries;
    if (_ignore_molecules_not_matching_queries) {
      return 1;
    }
    cerr << "Cannot identify matched atoms " << m.name() << " none of "
         << _queries.size() << " queries matched\n";
    return 0;
  }

  if (_write_parent_molecule) {
    output << m.smiles() << ' ' << m.name() << " PARENT\n";
  }

  const int nsys = ExtendMatchedAtomsAcrossFusions(m, matched_atoms.get());

//#define DEBUG_RING_REPLACEMENT
#ifdef DEBUG_RING_REPLACEMENT
  cerr << m.smiles() << ' ' << m.name() << " begin processing\n";
  cerr << "After assignment of ring systems\n";
  cerr << m.smiles() << ' ' << m.name() << '\n';
#endif

  int contains_ytterbium;
  if (_two_substituents_ytterbium) {
    const int initial_matoms = m.natoms();
    contains_ytterbium = TwoSubstituentsInsertYtterbium(m, matched_atoms.get());
    if (contains_ytterbium) {
      ResizeMatchedAtomsArray(matched_atoms, initial_matoms, m.natoms());
    }
  } else {
    contains_ytterbium = 0;
  }

  std::unique_ptr<SubstituentsForRingSystem[]> substituents =
                std::make_unique<SubstituentsForRingSystem[]>(nsys + 1);
  if (! MoleculeToSubstituents(m, matched_atoms.get(), substituents.get())) {
    ++_no_substituents;
    return 0;
  }

  // Process each ring system.
  uint32_t rc = 0;
  for (int i = 1; i <= nsys; ++i) {
    rc += Process(m, substituents[i], output);
    if (rc > _max_products_per_starting_molecule) {
      return rc;
    }
  }

  return rc;
}

// Process a given ring system. A ring system consists of one of
// more rings.
uint32_t
Options::Process(Molecule& m,
                 const SubstituentsForRingSystem& substituents,
                 IWString_and_File_Descriptor& output) {
  Molecule mcopy(m);
  substituents.RemoveRingSystem(mcopy);

  uint32_t rc = 0;
  for (const Replacement* r : _replacements) {
    rc += Process(mcopy, substituents, *r, output);
    if (rc > _max_products_per_starting_molecule) {
      return rc;
    }
  }

  return rc;
}

uint32_t
Options::Process(Molecule& m,
                 const SubstituentsForRingSystem& substituents,
                 const Replacement& replacement,
                 IWString_and_File_Descriptor& output) {

  // New ring system cannot accommodate the substituents.
  if (substituents.number_attachment_points() > 
      replacement.number_attachment_points()) {
    return 0;
  }

  // cerr << "Rings in substituents " << substituents.nrings() << '\n';

  if (substituents.nrings() == 1) {
    return ReplaceSingleRing(m, *substituents[0], replacement, output);
  }

  return ReplaceRingSystem(m, substituents, replacement, output);
}

uint32_t
Options::ReplaceRingSystem(Molecule& m,
                const SubstituentsForRingSystem& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output) {
  const int initial_matoms = m.natoms();

  m.add_molecule(&replacement.molecule());

  int rc;
  if (_preserve_ring_substitutions) {
    rc = ReplaceRingSystemPreserveRingSubstitution(m, initial_matoms, substituents, replacement, output);
  } else {
    rc = ReplaceRingSystemSubstituentsAnywhere(m, initial_matoms, substituents, replacement, output);
  }

  m.resize(initial_matoms);

  return rc;
}

// For each ring, we need to add the substituents to that ring, and
// then iterate over rings.
int
Options::ReplaceRingSystemPreserveRingSubstitution(Molecule& m,
                int initial_matoms,
                const SubstituentsForRingSystem& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output) {
  return 1;
}

int
Options::ReplaceRingSystemSubstituentsAnywhere(Molecule& m,
                int initial_matoms,
                const SubstituentsForRingSystem& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output) {
  const Set_of_Atoms sa = substituents.AttachmentPoints(m);
  // cerr << m.smiles() <<  " ReplaceRingSystemSubstituentsAnywhere\n";
  const uint32_t n = sa.size();

  // THis must have been sorted.
  Set_of_Atoms ra = replacement.attachment_points();
  ra.EachAtomIncrement(initial_matoms);

  if (ra.size() > sa.size()) {
    return ReplaceSingleRingInner(m, initial_matoms, sa, replacement, ra, output);
  }

  // cerr << "Replacement atoms " << ra << '\n';
  do {
    bool ok_adjacent_atoms = 1;
    if (! _allow_ortho_replacements && ! OkOrthoSubstituents(m, ra, n)) {
      continue;
    }
    for (uint32_t i = 0; i < n; ++i) {
      const atom_number_t a1 = sa[i];
      const atom_number_t a2 = ra[i];
      if (m.atomic_number(a1) != 6 && m.atomic_number(a2) != 6) {
        ok_adjacent_atoms = false;
      }
      if (m.formal_charge(a1) || m.formal_charge(a2)) {
        ok_adjacent_atoms = false;
      }
      m.add_bond(sa[i], ra[i], SINGLE_BOND);
    }

    if (ok_adjacent_atoms) {
      MaybeWrite(m, replacement, output);
    }

    for (uint32_t i = 0; i < n; ++i) {
      //m.remove_bond_between_atoms(sa[i], ra[i]);
      REMOVE_BOND_BETWEEN_ATOMS(m, sa[i], ra[i]);
    }
  } while (std::next_permutation(ra.begin(), ra.end()));

  return 1;
}

int
Options::ReplaceSingleRing(Molecule& m,
                const SubstituentsForRing& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output) {
  const int initial_matoms = m.natoms();

  m.add_molecule(&replacement.molecule());

  int rc = ReplaceSingleRingInner(m, initial_matoms, substituents, replacement, output);

  m.resize(initial_matoms);

  return rc;
}

// 'ring_atoms' are atom numbers in `m` to which we are going to attach
// `n` substituents.
// Examine all 
int
Options::OkOrthoSubstituents(const Molecule& m, const Set_of_Atoms& ring_atoms, int n) {
  if (n == 1) {
    return 1;
  }

  for (int i = 0; i < n; ++i) {
    atom_number_t ai = ring_atoms[i];
    for (int j = i + 1; j < n; ++j) {
      if (m.are_bonded(ai, ring_atoms[j])) {
        ++_ortho_substitutions_suppressed;
        return 0;
      }
    }
  }

  return 1;
}

// this works even if the number of attachment points in
// the replacement is larger. Note that we permute all
// attachment points from the replacement, and just use
// the first n of them.
int
Options::ReplaceSingleRingInner(Molecule& m,
                int initial_matoms,
                const SubstituentsForRing& substituents,
                const Replacement& replacement,
                IWString_and_File_Descriptor& output) {
  const Set_of_Atoms sa = substituents.AttachmentPoints(m);
  const uint32_t n = sa.size();

  // THis must have been sorted.
  Set_of_Atoms ra = replacement.attachment_points();
  ra.EachAtomIncrement(initial_matoms);

  // If the replacement has more attachment points than what is in
  // substituents, we need a different n choose k algorithm.

  // If more sites available in the replacement, we need a different algorithm.
  if (ra.size() > n) {
    return ReplaceSingleRingInner(m, initial_matoms, sa, replacement, ra, output);
  }

  // cerr << "Replacement atoms " << ra << '\n';
  do {
    bool ok_adjacent_atoms = 1;
    if (! _allow_ortho_replacements && ! OkOrthoSubstituents(m, ra, n)) {
      continue;
    }
    for (uint32_t i = 0; i < n; ++i) {
      const atom_number_t a1 = sa[i];
      const atom_number_t a2 = ra[i];
      if (m.atomic_number(a1) != 6 && m.atomic_number(a2) != 6) {
        ok_adjacent_atoms = false;
      }
      m.add_bond(sa[i], ra[i], SINGLE_BOND);
    }

    if (ok_adjacent_atoms) {
      MaybeWrite(m, replacement, output);
    }

    for (uint32_t i = 0; i < n; ++i) {
      //m.remove_bond_between_atoms(sa[i], ra[i]);
      REMOVE_BOND_BETWEEN_ATOMS(m, sa[i], ra[i]);
    }
  } while (std::next_permutation(sa.begin(), sa.end()));

  return 1;
}

// Return a vector of length `s1`. The last entries must be the
// numbers 1,2,3...s2.
// The first entries will be leading 0's.
std::vector<uint32_t>
LeadingZerosAndIota(uint32_t s1, uint32_t s2) {
  assert(s1 > s2);

  std::vector<uint32_t> result;
  result.reserve(s1);

  for (uint32_t i = 0; i < (s1 - s2); ++i) {
    result.push_back(0);
  }

  for (uint32_t i = 1; i <= s2; ++i) {
    result.push_back(i);
  }

  return result;
}

int
IndicesToAtoms(const std::vector<uint32_t>& indices, const Set_of_Atoms& ra, 
               size_t size_needed,
               Set_of_Atoms& destination) {
  destination.resize_keep_storage(0);

  const size_t n = indices.size();
  for (size_t i = 0; i < n; ++i) {
    int ndx = indices[i];
    if (ndx == 0) {
      continue;
    }
    destination << ra[ndx - 1];
  }

#ifdef DEBUG_INDICES_TO_ATOMS
  cerr << "IndicesToAtoms";
  for (uint32_t i : indices) {
    cerr << ' ' << i;
  }
  cerr << '\n';
  cerr << "Result " << destination << " size match " << size_needed << '\n';
#endif

  return destination.size() == size_needed;
}

// This variant is the case where the number of sites in the number
// of sites in `replacement` is greater than the number of sites in
// the substituents (sa). We use the algorithm suggested at
// https://stackoverflow.com/questions/64037553/iterative-algorithm-for-n-choose-k-without-repetitions-order-matters
int
Options::ReplaceSingleRingInner(Molecule& m,
                                int initial_matoms,
                                const Set_of_Atoms& sa,
                                const Replacement& replacement,
                                const Set_of_Atoms& replacement_atoms,
                                IWString_and_File_Descriptor& output) {
  assert(replacement_atoms.size() > sa.size());

  const uint32_t n = sa.size();

  // Load an array with leading 0's and then the indices 1,2,3...
  std::vector<uint32_t> v = LeadingZerosAndIota(replacement_atoms.size(), n);

  // Scope here for efficiency.
  Set_of_Atoms ra;
  ra.reserve(n);

#ifdef DEBUG_SINGLE_RING_INNER_X
  cerr << "replacement_atoms " << replacement_atoms << '\n';
  cerr << "initial_matoms " << initial_matoms << " natoms " << m.natoms() << '\n';

  cerr << "indices ";
  for (auto i : v) {
    cerr << ' ' << i;
  }
  cerr << '\n';
#endif

  do {
    if (! IndicesToAtoms(v, replacement_atoms, n, ra)) {
      continue;
    }
    bool ok_adjacent_atoms = 1;
    if (! _allow_ortho_replacements && ! OkOrthoSubstituents(m, ra, n)) {
      continue;
    }
    for (uint32_t i = 0; i < n; ++i) {
      const atom_number_t a1 = sa[i];
      const atom_number_t a2 = ra[i];
      if (m.atomic_number(a1) != 6 && m.atomic_number(a2) != 6) {
        ok_adjacent_atoms = false;
      }
      m.add_bond(sa[i], ra[i], SINGLE_BOND);
    }

    if (ok_adjacent_atoms) {
      MaybeWrite(m, replacement, output);
    }

    for (uint32_t i = 0; i < n; ++i) {
      //m.remove_bond_between_atoms(sa[i], ra[i]);
      REMOVE_BOND_BETWEEN_ATOMS(m, sa[i], ra[i]);
    }
  } while (std::next_permutation(v.begin(), v.end()));

  return 1;
}

atom_number_t
SameIsotope(const Molecule& m, atom_number_t zatom) {
  const isotope_t iso = m.isotope(zatom);

  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (m.isotope(o) == iso) {
      return o;
    }
  }

  return kInvalidAtomNumber;
}

isotope_t
HighestIsotope(const Molecule& m) {
  const int matoms = m.natoms();

  isotope_t rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso = m.isotope(i);
    if (iso > rc) {
      rc = iso;
    }
  }

  return rc;
}


// We have added a Ytterbium atom to handle the case of two
// substiuents on the same ring atom. Remove the Ytterbium atom
// and reconnect the two substituents.
// Return 0 if we fail.
int
RemoveYtterbiumTemporaryAtoms(Molecule& m, int nytterbium) {
  const int matoms = m.natoms();

  Set_of_Atoms y;
  Set_of_Atoms in_ring;
  Set_of_Atoms in_substituent;

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m[i];
    if (a.atomic_number() != kYtterbium) {
      continue;
    }

    y << i;
    for (const Bond* b : a) {
      atom_number_t o = b->other(i);
      if (m.ring_bond_count(o)) {
        in_ring << o;
      } else {
        in_substituent << o;
      }
    }

    if (y.number_elements() >= nytterbium) {
      break;
    }
  }

  if (y.size() != in_ring.size() ||
      in_ring.size() * 2 != in_substituent.size()) {
    cerr << "RemoveYtterbiumTemporaryAtoms:size mismatch\n";
    cerr << "Y atoms " << y << '\n';
    cerr << "Ring atoms " << in_ring << '\n';
    cerr << "Substituent atoms " << in_substituent << '\n';
    cerr << m.smiles() << ' ' << m.name() << '\n';
    // Need to do this otherwise we get infinite recursion in MaybeWrite.
    m.remove_all(kYtterbium);
    return 0;
  }

  // cerr << m.smiles() << " Before fixing ytt\n";

  const int n = y.number_elements();
  for (int i = 0; i < n; ++i) {
    const atom_number_t r = in_ring[i];
    const atom_number_t s1 = in_substituent[i + i];
    const atom_number_t s2 = in_substituent[i + i + 1];

    //m.remove_bond_between_atoms(y[i], r);
    REMOVE_BOND_BETWEEN_ATOMS(m, y[i], r);

    if (m.hcount(r) < 2) {  // would generate a valence error.
      return 0;
    }

    //m.remove_bond_between_atoms(y[i], s1);
    REMOVE_BOND_BETWEEN_ATOMS(m, y[i], s1);
    //m.remove_bond_between_atoms(y[i], s2);
    REMOVE_BOND_BETWEEN_ATOMS(m, y[i], s2);

    m.add_bond(r, s1, SINGLE_BOND);
    m.add_bond(r, s2, SINGLE_BOND);
  }

  // cerr << m.smiles() << " before removing Y atoms\n";
  if (y.size() == 1) {
    m.remove_atom(y.front());
    return 1;
  }

  for (int i = n - 1; i >= 0; --i) {
    m.remove_atom(y[i]);
  }

  return 1;
}

// Return 1 if any of the _products_must_contain queries
// match `m`.
int
Options::MatchesMustHaveQueries(Molecule& m) {
  Molecule_to_Match target(&m);

  Substructure_Results sresults;
  for (Substructure_Query* q : _products_must_contain) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

// Return 1 if any of the _products_cannot_contain queries
// match `m`.
int
Options::MatchesMustNotHaveQueries(Molecule& m) {
  Molecule_to_Match target(&m);

  Substructure_Results sresults;
  for (Substructure_Query* q : _products_cannot_contain) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

// Return 1 if the molecule satisfies any substructure constraints
int
Options::OkSubstructures(Molecule& m) {
  if (_products_cannot_contain.size() > 0) {
    if (MatchesMustNotHaveQueries(m)) {
      ++_discarded_for_query_mismatch;
      return 0;
    }
  }

  if (_products_must_contain.size() > 0) {
    if (! MatchesMustHaveQueries(m)) {
      ++_discarded_for_query_mismatch;
      return 0;
    }
  }

  return 1;
}

int
Options::MaybeWrite(Molecule& m, const Replacement& replacement,
                    IWString_and_File_Descriptor& output) {
  if (_two_substituents_ytterbium) {
    if (const int nytterbium = m.natoms(kYtterbium); nytterbium) {
      Molecule mcopy(m);
      if (RemoveYtterbiumTemporaryAtoms(mcopy, nytterbium)) {
        return MaybeWrite(mcopy, replacement, output);
      } else {
        return 0;
      }
    }
  }

  if (_discard_invalid_valence && ! m.valence_ok()) {
    ++_invalid_valence_discarded;
    return 0;
  }

  if (auto iter = _seen.find(m.unique_smiles()); iter != _seen.end()) {
    ++_duplicates_discarded;
    return 0;
  }

  _seen.insert(m.unique_smiles());

  if (! OkSubstructures(m)) {
    return 0;
  }

  m.invalidate_smiles();

  std::unique_ptr<isotope_t[]> isotopes;
  if (_remove_isotopes_from_products) {
    isotopes = m.GetIsotopes();
    m.unset_isotopes();
  }

  output << m.smiles() << ' ' << m.name() << " %% " << replacement.id() <<
            ' ' << replacement.count() << '\n';

  if (_remove_isotopes_from_products) {
    m.set_isotopes(isotopes.get());
  }

  output.write_if_buffer_holds_more_than(4192);

  return 1;
}

int
RingReplacementInexact(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
RingReplacementInexact(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! RingReplacementInexact(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
RingReplacementInexact(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "RingReplacementInexact:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return RingReplacementInexact(options, input, output);
}

int
RingReplacementInexact(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:s:q:I:R:z:pn:ex:oVY:N:X:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! RingReplacementInexact(options, fname, input_type, output)) {
      cerr << "RingReplacementInexact::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace ring_replacement_inexact

int
main(int argc, char ** argv) {

  int rc = ring_replacement_inexact::RingReplacementInexact(argc, argv);

  return rc;
}
