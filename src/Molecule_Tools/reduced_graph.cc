// Tool for implementing reduce graph concepts.
// The reduced features are specified as protos.
#include <stdlib.h>

#include <filesystem>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/reduced_graph.pb.h"

namespace reduced_graph {

namespace fs = std::filesystem;

using std::cerr;

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
  cerr << "Generates reduced rings similar to Harper et al https://doi.org/10.1021/ci049860f\n";
  cerr << " -R <fname>       file containing GraphReduction protos describing features\n";
  cerr << "                  the file contains names of proto files\n";
  cerr << " -x               add doubly bonded exocyclic heteroatoms to rings\n";
  cerr << " -I <fname>       write isotopically labelled molecules to <fname> - for debugging\n";
  cerr << " -m <natoms>      discard molecules having fewer than <matoms> after preprocessing\n";
  cerr << " -l               reduce to largest fragment\n";
  cerr << " -c               remove chirality\n";
  //cerr << " -T e1=e2         element transformation options\n";  not implemented
  cerr << " -N ...           charge assigner specification\n";
  cerr << " -H ...           donor acceptor specification\n";
  cerr << " -g ...           chemical standardisation\n";
  cerr << " -A ...           aromaticity options\n";
  cerr << " -E ...           element options\n";
  cerr << " -v               verbose output\n";
  // clang-format on

  ::exit(rc);
}

// When discerning rings, we perceive a ring, which is then assigned an
// abstract atom type. Keep track of that info here.
class TypeAndRingAtoms {
  private:
    // The type of the ring - isotope to be applied.
    int _type;

    // The element for the abstract atom.
    const Element * _element;

    // The abstract atom number that was placed in the molecule.
    // This will not be known until the ring is emplaced.
    atom_number_t _atom;

    // The rings in this ring.
    Set_of_Atoms _atoms;

  public:
    TypeAndRingAtoms(int atype, const Element* e, const Set_of_Atoms& s);

    const Element* element() const {
      return _element;
    }

    const Set_of_Atoms& atoms() const {
      return _atoms;
    }

    void set_atom(atom_number_t s) {
      _atom = s;
    }
    atom_number_t atom() const {
      return _atom;
    }

    void add_atom(atom_number_t zatom) {
      _atoms << zatom;
    }

    int AtomsInCommon(const TypeAndRingAtoms& rhs) const;
};

TypeAndRingAtoms::TypeAndRingAtoms(int atype, const Element* e,
                const Set_of_Atoms& s) :
                _type(atype),
                _element(e),
                _atoms(s) {
  _atom = INVALID_ATOM_NUMBER;
  auto cmp = Int_Comparator_Larger();
  _atoms.iwqsort(cmp);
}

// Return the number of atoms in common between _atoms and rhs._atoms.
// Depends on both lists being sorted.
int
TypeAndRingAtoms::AtomsInCommon(const TypeAndRingAtoms& rhs) const {
  int rc = 0;
  const int n1 = _atoms.number_elements();
  const int n2 = rhs._atoms.number_elements();
  int i1 = 0;
  int i2 = 0;
  while (i1 < n1 && i2 < n2) {
    const atom_number_t a1 = _atoms[i1];
    const atom_number_t a2 = rhs._atoms[i2];
    if (a1 == a2) {
      ++rc;
      ++i1;
      ++i2;
    } else if (a1 < a2) {
      ++i1;
    } else {
      ++i2;
    }
  }

  return rc;
}

// Description of a specific graph reduction.
class GraphReductionSpec {
  private:
    // Queries that define the feature.
    resizable_array_p<Substructure_Query> _queries;

    // The dummy element assigned to this feature type.
    const Element* _element;

    // Does this describe a ring?
    bool _ring;

    // Each reduction is assigned a number that indicates its priority.
    int _my_number;

    // Whether or not doubly bonded exocyclics are included with rings.
    int _include_doubly_bonded_heteroatoms_with_rings;

  // private functions

    int SetElement(const GraphReduction& proto);
    int AssignChain(Molecule_to_Match& target, int * rgroup,
                               const resizable_array_p<Set_of_Atoms>& embeddings);
    int AssignRing(Molecule_to_Match& target, int * rgroup,
                               int * storage,
                               Substructure_Results& sresults,
                               resizable_array_p<TypeAndRingAtoms> & rings);
  public:
    GraphReductionSpec();

    void SetNumber(int s) {
      _my_number = s;
    }
    int Number() const {
      return _my_number;
    }

    const Element* TransformElement() const {
      return _element;
    }

    bool is_ring() const {
      return _ring;
    }

    void set_include_doubly_bonded_heteroatoms_with_rings(int s) {
      _include_doubly_bonded_heteroatoms_with_rings = s;
    }

    // Since a GraphReduction proto may have file names inside it, building
    // needs to know the name of the proto file so names can be found in
    // the same directory.
    int BuildFromProto(IWString& outer_file_name, const GraphReduction& proto);

    // Perform substructure searches in `target`.
    // If matches are found, update `rgroup`.
    // If we are a ring, then also update `rings`.
    // `storage` is just scratch space we use.
    int Assign(Molecule_to_Match& target, int * rgroup,
               int * storage, resizable_array_p<TypeAndRingAtoms>& rings);
};

GraphReductionSpec::GraphReductionSpec() {
  _my_number = 1;
}

int
GraphReductionSpec::Assign(Molecule_to_Match& target,
                           int * rgroup,
                           int * storage,
                           resizable_array_p<TypeAndRingAtoms>& rings) {
  int rc = 0;
  Substructure_Results sresults;
  for (Substructure_Query * q : _queries) {
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
#ifdef DEBUG_ASSIGN
    cerr << "Got " << nhits << " to query " << _my_number << " " <<  _ring << '\n';
#endif

    ++rc;
    if (_ring) {
      AssignRing(target, rgroup, storage, sresults, rings);
    } else {
      AssignChain(target, rgroup, sresults.embeddings());
    }
  }

  return rc;
}

int
GraphReductionSpec::AssignChain(Molecule_to_Match& target,
                                int * rgroup,
                                const resizable_array_p<Set_of_Atoms>& embeddings) {
  for (const Set_of_Atoms * e : embeddings) {
    for (atom_number_t i : *e) {
      if (rgroup[i] == 0) {
        rgroup[i] = _my_number;
      }
    }
  }

  return 1;
}

//#define DEBUG_ASSIGN_RING
// Ring queries need to be handled carefully.
// Each embedding will contain 1 atom, so we assemble a list of all matched atoms.
// Then for each ring in the molecule, we identify those rings that have every member
// set in that list. Those are the rings.
int
GraphReductionSpec::AssignRing(Molecule_to_Match& target,
                               int * rgroup,
                               int * storage,
                               Substructure_Results& sresults,
                               resizable_array_p<TypeAndRingAtoms> & rings) {
#ifdef DEBUG_ASSIGN_RING
  cerr << "Ring abstraction " << _my_number << " ring " << _ring << " has " << sresults.number_embeddings() << " hits\n";
#endif
  Molecule & m = *target.molecule();
  const int matoms = m.natoms();
  std::fill_n(storage, matoms, 0);
  sresults.each_embedding_set_vector(storage, 1);

  const int nr = m.nrings();
  for (int i = 0; i < nr; ++i) {
    const Set_of_Atoms* ri = m.ringi(i);
#ifdef DEBUG_ASSIGN_RING
    cerr << "Ring " << *ri << " all it " << ri->all_members_set_in_array(storage, 1) << '\n';
#endif
    if (! ri->all_members_set_in_array(storage, 1)) {
      continue;
    }
    int any_atoms_hit = 0;  // Maybe another query has already covered this ring.
    for (auto j : *ri) {
      if (rgroup[j] == 0) {
        rgroup[j] = _my_number;
        ++any_atoms_hit;
      }
    }
    if (! any_atoms_hit) {
      continue;
    }

    TypeAndRingAtoms * ratn = new TypeAndRingAtoms(_my_number, _element, *ri);
    rings << ratn;
    if (_include_doubly_bonded_heteroatoms_with_rings) {
      for (auto ring_atom : *ri) {
        const Atom * a = m.atomi(ring_atom);
        if (a->ncon() == 2) {
          continue;
        }
        for (const Bond * b : *a) {
          if (b->is_single_bond()) {
            continue;
          }
          atom_number_t o = b->other(ring_atom);
          if (m.atomic_number(o) == 6) {
            continue;
          }
          if (rgroup[o]) {  // also checks in `ri`.
            continue;
          }
          rgroup[o] = _my_number;
          ratn->add_atom(o);
        }
      }
    }
  }

  return 1;
}

// A collection of GraphReductionSpec.
class SetOfGraphReductions {
  private:
    resizable_array_p<GraphReductionSpec> _graph_reductions;

    // We need to identify the isolating carbon query.
    int _isolating_carbon_type;

  // private functions.
    int Initialise(IWString& fname, iwstring_data_source& input);

  public:
    SetOfGraphReductions();

    int Initialise(IWString& fname);

    // Call the Assign method of each of the _graph_reductions.
    int Assign(Molecule& m, int * rgraph,
               int * storage, resizable_array_p<TypeAndRingAtoms>& rings);

    resizable_array_p<GraphReductionSpec>& Reductions() {
      return _graph_reductions;
    }

    size_t size() const { return _graph_reductions.size();}

    int isolating_carbon_type() const {
      return _isolating_carbon_type;
    }

    // Passed to each GraphReductionSpec.
    void set_include_doubly_bonded_heteroatoms_with_rings(int s);
};

SetOfGraphReductions::SetOfGraphReductions() {
  _isolating_carbon_type = -1;
}

int
SetOfGraphReductions::Initialise(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "SetOfGraphReductions::Initialise:cannot open '" << fname << "'\n";
    return 0;
  }

  return Initialise(fname, input);
}

void
SetOfGraphReductions::set_include_doubly_bonded_heteroatoms_with_rings(int s) {
  for (GraphReductionSpec * grs : _graph_reductions) {
    grs->set_include_doubly_bonded_heteroatoms_with_rings(s);
  }
}

int
SetOfGraphReductions::Initialise(IWString& outer_file_name, iwstring_data_source& input) {
  IWString fname;
  input.set_skip_blank_lines(1);
  while (input.next_record(fname)) {
    std::optional<IWString> proto_fname = iwmisc::FileOrPath(outer_file_name, fname);
    if (! proto_fname) {
      cerr << "SetOfGraphReductions::Initialise:cannot find '" << fname << "'\n";
      return 0;
    }
    std::optional<reduced_graph::GraphReduction> proto = iwmisc::ReadTextProtoCommentsOK<reduced_graph::GraphReduction>(*proto_fname);
    if (! proto) {
      cerr << "SetOfGraphReductions::Initialise:cannot read proto file " << *proto_fname << "\n";
      return 0;
    }
    std::unique_ptr<GraphReductionSpec> grs = std::make_unique<GraphReductionSpec>();
    if (! grs->BuildFromProto(outer_file_name, *proto)) {
      cerr << "SetOfGraphReductions::Initialise:cannot parse " << proto->ShortDebugString() << '\n';
      return 0;
    }

    grs->SetNumber(_graph_reductions.number_elements() + 1);
    _graph_reductions << grs.release();
  }

  // Determine which one is the isolating carbon query, Zinc.
  const Element* zinc = get_element_from_atomic_number(30);

  for (const GraphReductionSpec* spec : _graph_reductions) {
    if (spec->TransformElement() == zinc) {
      _isolating_carbon_type = spec->Number();
    }
  }

  return _graph_reductions.number_elements();
}

int
GraphReductionSpec::BuildFromProto(IWString& outer_file_name, const GraphReduction& proto) {
  if (proto.query_size() == 0) {
    cerr << "GraphReductionSpec:Initialise:no queries " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! SetElement(proto)) {
    cerr << "GraphReductionSpec::BuildFromProto:missing or invalid element\n";
    return 0;
  }

  for (const std::string& fname : proto.query()) {
    std::optional<IWString> proto_fname = iwmisc::FileOrPath(outer_file_name, fname);
    if (! proto_fname) {
      cerr << "GraphReductionSpec::Initialise:cannot find '" << fname << "'\n";
      cerr << "outer_file_name " << outer_file_name << '\n';
      return 0;
    }

    std::optional<SubstructureSearch::SubstructureQuery> proto =
          iwmisc::ReadTextProtoCommentsOK<SubstructureSearch::SubstructureQuery>(*proto_fname);
    if (! proto) {
      cerr << "GraphReductionSpec::BuildFromProto:cannot read " << *proto_fname << '\n';
      return 0;
    }

    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    if (! qry->ConstructFromProto(*proto)) {
      cerr << "GraphReductionSpec::Initialise:cannot parse " << proto->ShortDebugString() << '\n';
      return 0;
    }
    _queries << qry.release();
  }

  return 1;
}

int
GraphReductionSpec::SetElement(const GraphReduction& proto) {
  if (proto.element().empty()) {
    cerr << "GraphReduction::SetElement:no element\n";
    return 0;
  }

  IWString ele(proto.element());
  _element = get_element_from_symbol_no_case_conversion(ele);
  if (_element == nullptr) {
    cerr << "GraphReduction::SetElement:invalid element '" << ele << "'\n";
    return 0;
  }
  _ring = proto.ring();
  return 1;
}

int
SetOfGraphReductions::Assign(Molecule& m,
                             int * rgroup,
                             int * storage,
                             resizable_array_p<TypeAndRingAtoms>& rings) {
  Molecule_to_Match target(&m);
  for (GraphReductionSpec * grf : _graph_reductions) {
    grf->Assign(target, rgroup, storage, rings);
  }

  return 1;
}

// Holds most of the configurable things associated with this.
// Primarily instantiated from the command line.
struct Options {
    int verbose = 0;

    FileType input_type = FILE_TYPE_INVALID;

    int reduce_to_largest_fragment = 0;

    int remove_chirality = 0;

    Chemical_Standardisation chemical_standardisation;

    Element_Transformations element_transformations;

    Charge_Assigner charge_assigner;

    Donor_Acceptor_Assigner donor_acceptor_assigner;

    SetOfGraphReductions graph_reductions;

    int molecules_read = 0;

    IWString_and_File_Descriptor stream_for_isotopically_labelled;

    // As an alternative to isotopic labels. But cactvs seems to
    // not display atom map numbers.
    int label_atoms_by_atom_map_number = 0;

    // It is useful to add doubly bonded heteroatoms to rings.
    int include_doubly_bonded_heteroatoms_with_rings = 0;

    // The Gillet paper recursively removes terminal groups.
    int remove_terminal_groups = 1;
    // Some molecules will be destroyed by terminal group removal.
    int empty_after_terminal_removal = 0;
    // For those atomic numbers that are removed if they are
    // found singly connected. Carbon, halogens...
    extending_resizable_array<int> remove_if_terminal;

    // Removing terminal groups may leave a molecule that is too small.
    int min_atoms_needed = 1;
    int molecules_with_too_few_atoms = 0;

    // The Gillet paper transforms unclassified heteroatoms to Carbon.
    int transform_unclassified_heteroatoms_to_carbon = 0;

    // Statistics on coverage by the queries.
    extending_resizable_array<int> atoms_matched;
    extending_resizable_array<int> atoms_not_matched;
    Accumulator<double> coverage;

    // A mapping from isotopes applied to elements to assign.
    //std::unordered_map<int, const Element*> iso_to_element;
    const Element** iso_to_element;

    // Private functions

  private:
    int InitialiseIsoToElement();

    int MaybeRemoveAtom(Molecule& m, atom_number_t zatom,
                int * is_terminal, Set_of_Atoms& examine_next_step);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int RemoveTerminalGroups(Molecule& m);

    // For all atoms in `m` that are singly connected, and set in remove_if_terminal
    // mark the corresponding entry in `is_terminal`.
    int IdentifyTerminalGroups(Molecule& m, int * is_terminal);

    int MaybeWriteIsotopicallyLabelled(Molecule& m, const int * rgroup);

    int TransformUnclassifiedHeteroatomsToCarbon(Molecule& m, int * reduction);

    int UpdateMatchedAtomsCounters(const Molecule& m, const int * rgroup);

    int Report(std::ostream& output) const;
};

Options::Options() {
  iso_to_element = nullptr;

  remove_if_terminal[53] = 1;
  remove_if_terminal[35] = 1;
  remove_if_terminal[17] = 1;
  remove_if_terminal[9] = 1;
  remove_if_terminal[6] = 1;
}

Options::~Options() {
  if (iso_to_element != nullptr) {
    delete [] iso_to_element;
  }
}

int
Options::Initialise(Command_Line& cl) {

  verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('N')) {
    if (!charge_assigner.construct_from_command_line(cl, verbose, 'N')) {
      cerr << "Cannot determine charge assigner from command line\n";
      Usage(1);
    }

    if (Daylight != global_aromaticity_type()) {
      set_global_aromaticity_type(Daylight);
      if (verbose)
        cerr << "Global aromaticity type set to Dayligth for Charge Assigner";
    }
  }

  if (cl.option_present('H')) {
    if (!donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose)) {
      cerr << "Cannot initialise donor/acceptor assigner (-H option)\n";
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!element_transformations.construct_from_command_line(cl, verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    remove_chirality = 1;
    if (verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (cl.option_present('x')) {
    include_doubly_bonded_heteroatoms_with_rings = 1;
    if (verbose) {
      cerr << "Will include_doubly_bonded_heteroatoms_with_rings\n";
    }
  }

  if (! cl.option_present('R')) {
    cerr << "Must specify graph reductions via the -R option\n";
    return 0;
  }

  IWString fname;
  for (int i = 0; cl.value('R', fname, i); ++i) {
    if (! graph_reductions.Initialise(fname)) {
      cerr << "Cannot initialise graph reduction '" << fname << "'\n";
      return 0;
    }
  }

  if (verbose) {
    cerr << "Read " << graph_reductions.size() << " graph reductions\n";
  }

  if (include_doubly_bonded_heteroatoms_with_rings) {
    graph_reductions.set_include_doubly_bonded_heteroatoms_with_rings(1);
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', min_atoms_needed) || min_atoms_needed < 1) {
      cerr << "The min atoms needed option (-m) must be a whole +ve number\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will skip molecules having fewer than " << min_atoms_needed << " atoms\n";
    }
  }

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

  if (cl.option_present('I')) {
    IWString fname;
    IWString token;
    for (int i = 0; cl.value('I', token, i); ++i) {
      if (token == "amap") {
        label_atoms_by_atom_map_number = 1;
      } else {
        fname = token;
      }
    }
    if (fname.empty()) {
      cerr << "Must specify file name for isotopically labelled smiles\n";
      return 0;
    }

    if (! fname.ends_with(".smi")) {
      fname << ".smi";
    }
    if (! stream_for_isotopically_labelled.open(fname.null_terminated_chars())) {
      cerr << "Options::initialise:cannot open stream for isotopically labelled molecules '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "Isotopically labelled molecules written to " << fname << '\n';
    }
  }

  return InitialiseIsoToElement();
}

int
Options::InitialiseIsoToElement() {
  int highest_isotope = 0;
  for (const GraphReductionSpec* grdspec : graph_reductions.Reductions()) {
    if (grdspec->Number() > highest_isotope) {
      highest_isotope = grdspec->Number();
    }
  }

  if (highest_isotope == 0) {
    return 0;
  }

  iso_to_element = new const Element*[highest_isotope + 1];
  std::fill_n(iso_to_element, highest_isotope + 1, nullptr);

  for (const GraphReductionSpec* grdspec : graph_reductions.Reductions()) {
    int number = grdspec->Number();
    const Element* ele = grdspec->TransformElement();
    iso_to_element[number] = ele;
  }

//#define DEBUG_ISO_TO_ELEMENT
#ifdef DEBUG_ISO_TO_ELEMENT
  for (int i = 1; i <= highest_isotope; ++i) {
    cerr << i << " iso_to_element " << iso_to_element[i]->symbol() << '\n';
  }
#endif

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (remove_chirality) {
    m.remove_all_chiral_centres();
  }

#ifdef DO_NOT_DO_THIS_HERE
  if (remove_terminal_groups) {
    RemoveTerminalGroups(m);
  }
#endif

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (m.natoms() < min_atoms_needed) {
    ++molecules_with_too_few_atoms;
    return 0;
  }

  if (charge_assigner.active()) {
    charge_assigner.process(m);
  }

  if (donor_acceptor_assigner.active()) {
    donor_acceptor_assigner.process(m);
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << molecules_read << " molecules\n";
  output << molecules_with_too_few_atoms << " molecules with too few atoms\n";

  for (int i = 0; i < atoms_matched.number_elements(); ++i) {
    if (atoms_matched[i]) {
      output << atoms_matched[i] << " molecules had " << i << " matched atoms\n";
    }
  }

  output << "Coverage btw " << coverage.minval() << " and " << coverage.maxval() << " mean " << coverage.average() << '\n';

  return 1;
}

int
Options::UpdateMatchedAtomsCounters(const Molecule& m,
                                    const int * rgroup) {
  const int matoms = m.natoms();
  int matched = count_non_zero_occurrences_in_array(rgroup, matoms);

  ++atoms_matched[matched];
  ++atoms_not_matched[matoms - matched];
  coverage.extra(iwmisc::Fraction<double>(matched, matoms));

  return 1;
}

int
Options::MaybeWriteIsotopicallyLabelled(Molecule& m, const int * rgroup) {
  if (! stream_for_isotopically_labelled.is_open()) {
    return 0;
  }

  const int matoms = m.natoms();
  int matched = 0;
  for (int i = 0; i < matoms; ++i) {
    if (rgroup[i] == 0) {
      continue;
    }

    ++matched;
    if (label_atoms_by_atom_map_number) {
      m.set_atom_map_number(i, rgroup[i]);
    } else {
      m.set_isotope(i, rgroup[i]);
    }
  }

  float cov = iwmisc::Fraction<double>(matched, matoms);
  coverage.extra(cov);

  stream_for_isotopically_labelled << m.smiles() << ' ' << m.name() << ' ' << cov << '\n';
  stream_for_isotopically_labelled.write_if_buffer_holds_more_than(2048);

  m.transform_to_non_isotopic_form();

  return 1;
}

// Identify atoms in `m` that are in the same `rgroup` as `zatom`.
// Atoms in the group are set in `in_group`, atoms that are attached
// but in different groups are set in `attached_to_group`.
int
IdentifyContiguousSection(Molecule& m,
                          const int * rgroup,
                          atom_number_t zatom,
                          int * in_group,
                          int * attached_to_group) {
  const Atom * atom = m.atomi(zatom);
  int rc = 1;
  in_group[zatom] = 1;
  for (const Bond * b : *atom) {
    atom_number_t j = b->other(zatom);
    if (in_group[j]) {
      continue;
    }
    if (attached_to_group[j]) {
      continue;
    }

#ifdef DEBUG_IDENTIFYCONTIGUOUSSECTION
    cerr << "From " << zatom << " to " << j << " same group " << (rgroup[j] == rgroup[zatom]) << '\n';
#endif
    if (rgroup[j] == rgroup[zatom]) {
      rc += IdentifyContiguousSection(m, rgroup, j, in_group, attached_to_group);
    } else {
      attached_to_group[j] = 1;
    }
  }

  return rc;
}

// Rings matches are stored in `rings`. For each ring, add a new
// abstract atom to `reduced_form` and update both `new_group` and `xref`.
int
AddAbstractRings(const resizable_array_p<TypeAndRingAtoms>& rings,
                 int * new_group,
                 int * xref,
                 Molecule& reduced_form) {
  int uid = 1;  // unique id for each group.
#ifdef DEBUG_ADD_ABSTRACT_RINGS
  cerr << "Adding " << rings.size() << " rings\n";
#endif
  for (TypeAndRingAtoms* ring : rings) {
    atom_number_t new_atom = reduced_form.natoms();
    ring->set_atom(new_atom);
    reduced_form.add(ring->element());
    for (atom_number_t a : ring->atoms()) {
      reduced_form.add_bond(a, new_atom, SINGLE_BOND);
      xref[a] = new_atom;
      new_group[a] = uid;
    }
    ++uid;
  }

  return uid;
}

// If there have been multiple rings discerned, see if any of them
// share atoms and therefore need to be bonded, but adding bonds to `reduced_form`.
int
MakeBondsBetweenFusedRings(const resizable_array_p<TypeAndRingAtoms>&rings,
                Molecule& reduced_form) {
  const int nr = rings.number_elements();
  if (nr < 2) {
    return 1;
  }

  for (int i = 0; i < nr; ++i) {
    const TypeAndRingAtoms* r1 = rings[i];
    for (int j = i + 1; j < nr; ++j) {
      const TypeAndRingAtoms* r2 = rings[j];
      const int aic = r1->AtomsInCommon(*r2);
      //cerr << "i " << i << " j " << j << " aic " << aic << '\n';
      if (aic == 0) {
        continue;
      }

      if (aic == 1) {  // Spiro fusion
        reduced_form.add_bond(r1->atom(), r2->atom(), SINGLE_BOND);
      } else if (aic == 2) {  // Normal fusion
        reduced_form.add_bond(r1->atom(), r2->atom(), DOUBLE_BOND);
      } else {
        reduced_form.add_bond(r1->atom(), r2->atom(), TRIPLE_BOND);
      }
    }
  }

  return 1;
}

// Recursively remove terminal groups, as determined by `remove_if_terminal`.
int
Options::RemoveTerminalGroups(Molecule& m) {
  Set_of_Atoms removed_this_step;
  do {
    removed_this_step.resize_keep_storage(0);
    const int matoms = m.natoms();
    for (int i = 0; i < matoms; ++i) {
      const Atom * a = m.atomi(i);
      if (a->ncon() != 1) {
        continue;
      }
      if (remove_if_terminal[a->atomic_number()]) {
        removed_this_step << i;
      }
    }

    if (removed_this_step.size() > 0) {
      m.remove_atoms(removed_this_step);
    }
  } while (removed_this_step.size() > 0);

  if (m.natoms() == 0) {
    ++empty_after_terminal_removal;
  }

  return m.natoms();
}

// If `zatom` looks like a terminal atom, mark it in `is_terminal`
// and add the atom to which it connected to `examine_next_step`.
int
Options::MaybeRemoveAtom(Molecule& m,
                atom_number_t zatom,
                int * is_terminal,
                Set_of_Atoms& examine_next_step) {
  const Atom * a = m.atomi(zatom);
  if (a->ncon() != 1) {
    return 0;
  }
  if (! remove_if_terminal[a->atomic_number()]) {
    return 0;
  }
  is_terminal[zatom] = 1;
  atom_number_t j = a->other(zatom, 0);
  examine_next_step << j;
  m.remove_bond_between_atoms(zatom, j);
  m.unset_all_implicit_hydrogen_information(j);
  return 1;
}

int
Options::IdentifyTerminalGroups(Molecule& m,
                                int * is_terminal) {
  // Only atoms from which a terminal group has been removed need to be
  // examined on each iteration.
  Molecule mcopy(m);  // Expensive, but too hard to keep track of otherwise.
  Set_of_Atoms examine_next_step;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    MaybeRemoveAtom(mcopy, i, is_terminal, examine_next_step);
  }

  while (examine_next_step.size() > 0) {
    auto to_shift = examine_next_step.number_elements();
    for (int i = 0; i < to_shift; ++i) {
      atom_number_t zatom = examine_next_step[i];
      MaybeRemoveAtom(mcopy, zatom, is_terminal, examine_next_step);
    }
    // If the array is same size as when we started, we did not find anything.
    if (examine_next_step.number_elements() == to_shift) {
      return m.natoms();
    }
    examine_next_step.erase(0, to_shift);
  }

  return 1;
}

int
Options::TransformUnclassifiedHeteroatomsToCarbon(Molecule& m,
                                        int * reduction) {
  int seemingly_carbon = graph_reductions.isolating_carbon_type();
  if (seemingly_carbon <= 0) {
    return 0;
  }

  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (reduction[i] > 0) {
      continue;
    }
    reduction[i] = seemingly_carbon;
    ++rc;
  }

  return rc;
}

//#define DEBUG_REDUCED_FORM

// Generate a reduced form from `m`.
// Individual reductions are marked in `reduction`.
// `rings` holds information about pseudo atoms defined by rings.
int
ReducedForm(Options& options,
            Molecule& m,
            int * reduction,
            resizable_array_p<TypeAndRingAtoms>& rings,
            IWString_and_File_Descriptor& output) {
  Molecule reduced_form(m);
  const int matoms = m.natoms();

  // One temporary rather than many.
  int * in_group = new int[matoms + matoms + matoms + matoms + matoms]; std::unique_ptr<int[]> free_in_group(in_group);
  int * attached_to_group = in_group + matoms;
  int * new_group = in_group + matoms + matoms;
  std::fill_n(new_group, matoms, 0);
  int * xref = in_group + matoms + matoms + matoms;
  std::iota(xref, xref + matoms, 0);
  int * is_terminal = in_group + matoms + matoms + matoms + matoms;
  std::fill_n(is_terminal, matoms, 0);

  AddAbstractRings(rings, new_group, xref, reduced_form);
  MakeBondsBetweenFusedRings(rings, reduced_form);
  options.IdentifyTerminalGroups(m, is_terminal);
  // terminal atoms will be removed, so let's not do any abstractions for them.
  for (int i = 0; i < matoms; ++i) {
    if (is_terminal[i]) {
      reduction[i] = 0;
    }
  }

  // Process all non-ring contiguous groups.
  // Note that we do not check for single atom groups.

  int uid = rings.number_elements() + 1;  // Ensure unique id for new groups.
  for (int i = 0; i < matoms; ++i) {
    if (new_group[i] > 0) {  // Already assigned.
      continue;
    }
    if (reduction[i] == 0) {  // Not changing, or don't care.
      continue;
    }
    if (m.is_ring_atom(i)) {  // Abstract ring atoms already added.
      continue;
    }
    std::fill_n(in_group, matoms + matoms, 0);
    IdentifyContiguousSection(m, reduction, i, in_group, attached_to_group);
    atom_number_t new_atom = reduced_form.natoms();
    reduced_form.add(options.iso_to_element[reduction[i]]);  // Even if just one atom in the group.
#ifdef DEBUG_REDUCED_FORM
    cerr << "Begin group starting with atom " << i << '\n';
#endif
    for (int j = 0; j < matoms; ++j) {
      if (in_group[j]) {
        new_group[j] = uid;
        xref[j] = new_atom;
      } else if (attached_to_group[j]) {
        reduced_form.add_bond(new_atom, j, SINGLE_BOND);
      }
    }
    ++uid;
  }

#ifdef DEBUG_REDUCED_FORM
  cerr << "xref\n";
  for (int i = 0; i < matoms; ++i) {
    cerr << i << ' ' << " xref " << xref[i] << " reduction " << reduction[i] << ' ' << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

  // Capture bonds between rings and unchanged atoms.
  // Note that this over-enthusiastically forms bonds within certain
  // fused systems, C(=O)C1=NNC2=C1CCC1=C2N=C(N)N=C1 CHEMBL538195.begin
  for (const TypeAndRingAtoms* ring : rings) {
    for (atom_number_t i : ring->atoms()) {
      const Atom * a = reduced_form.atomi(i);
      for (const Bond * b : *a) {
        atom_number_t j = b->other(i);
        if (j >= matoms) {  // Bond to pseudo atom...
          continue;
        }
        if (xref[j] != j) {  // mapped to a pseudo atom
          continue;
        }
        if (! reduced_form.are_bonded(ring->atom(), j)) {
          reduced_form.add_bond(ring->atom(), j, b->btype());
        }
      }
    }
  }

  // Any atoms that are not transforming to a reduced form need to have their bonding adjusted

  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = reduced_form.bondi(i);
    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
#ifdef DEBUG_REDUCED_FORM
    cerr << "Bond is " << *b << " xref " << xref[a1] << " and " << xref[a2] << '\n';
#endif
    int x1 = xref[a1];
    int x2 = xref[a2];

    if (x1 == x2) {  // Both atoms mapped to the same pseudo atom.
      continue;
    }
    if (x1 == a1 && x2 == a2) {  // Atoms are unchanged.
      continue;
    }

    // Both atoms are being replaced by pseudo atoms, make a bond between them.
    if (x1 >= matoms && x2 >= matoms) {  // Bond between abstract atoms
      if (! reduced_form.are_bonded(x1, x2)) {
        reduced_form.add_bond(x1, x2, SINGLE_BOND);
      }
      continue;
    }
  }

#ifdef DEBUG_REDUCED_FORM
  cerr << "Before removing atoms\n";
  reduced_form.debug_print(cerr);
#endif
  for (int i = matoms - 1; i >= 0; --i) {
    if (new_group[i] > 0 || is_terminal[i]) {
      reduced_form.remove_atom(i);
    }
  }
#ifdef DEBUG_REDUCED_FORM
  reduced_form.debug_print(cerr);
  reduced_form.invalidate_fragment_membership();
#endif

  output << m.smiles() << ' ' << m.name() << ".reduced\n";
  output << reduced_form.smiles() << ' ' << m.name() << ".abstract\n";
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
ReducedGraph(Options& options,
             Molecule& m,
             IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  int * reduction = new int[matoms + matoms];
  std::unique_ptr<int[]> free_reduction(reduction);
  std::fill_n(reduction, matoms, 0);
  int * storage = reduction + matoms;

  Molecule reduced_form(m);
  resizable_array_p<TypeAndRingAtoms> rings;
  if (! options.graph_reductions.Assign(m, reduction, storage, rings)) {
    cerr << "Graph reduction failed\n";
    return 0;
  }
#ifdef DEBUG_REDUCED_GRAPH
  cerr << "Returned with " << rings.size() << " rings\n";
#endif

  if (options.transform_unclassified_heteroatoms_to_carbon) {
    options.TransformUnclassifiedHeteroatomsToCarbon(m, reduction);
  }

  options.MaybeWriteIsotopicallyLabelled(m, reduction);

  options.UpdateMatchedAtomsCounters(m, reduction);

  return ReducedForm(options, m, reduction, rings, output);
}

int
ReducedGraph(Options& options,
             data_source_and_type<Molecule>& input,
             IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);
    ++options.molecules_read;

    if (! options.Preprocess(*m)) {
      continue;
    }
#ifdef DEBUG_REDUCED_GRAPH
    output << m->smiles() << " begin\n";
#endif

    if (! ReducedGraph(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ReducedGraph(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  if (options.input_type == FILE_TYPE_INVALID) {
    options.input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(options.input_type, fname);
  if (! input.good()) {
    cerr << "ReducedGraph: cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose > 1) {
    input.set_verbose(1);
  }

  return ReducedGraph(options, input, output);
}

int
ReducedGraph(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:H:N:T:A:lcg:i:I:R:m:x");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! ReducedGraph(options, fname, output)) {
      cerr << "ReducedGraph::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace reduced_graph

int
main(int argc, char ** argv) {

  int rc = reduced_graph::ReducedGraph(argc, argv);

  return rc;
}
