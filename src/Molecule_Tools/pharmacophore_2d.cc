// Generate pharmacophore-like queries for 2D molecules

#include <iostream>
#include <memory>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"
#include "Molecule_Lib/target.h"

namespace pharacaphore_2d {
using std::cerr;
using std::cout;

int verbose = 0;

int molecules_read = 0;
int no_functional_groups = 0;
int only_one_functional_group = 0;

resizable_array_p<Substructure_Query> atoms_to_ignore;

Accumulator_Int<int> functional_group_count;
Accumulator<float> fraction_atoms_in_functional_groups;

// When generating distance constraints, ignore atoms outside the
// range specified here.
int min_separation = 1;
int max_separation = std::numeric_limits<int>::max();

// When two atoms are separated by `d` bonds, the min_bonds_between
// will be `d - delta_shorter` (if positive) and the max_bonds_between
// will be `d + delta_longer`.
int delta_shorter = 0;
int delta_longer = 0;

int ncon_becomes_min_ncon = 0;

IWString_and_File_Descriptor stream_for_labelled_smiles;

void
Usage(int rc) {
  cerr << " -n               ncon values become min_ncon values in the query\n";
  cerr << " -d <dist>        min separation between atoms\n";
  cerr << " -D <dist>        max separation between atoms\n";
  cerr << " -L <fname>       write labelled molecule to <fname>\n";
  cerr << " -v               verbose output\n";
  exit(rc);
}

void
Preprocess(Molecule& m) {
}

void
AddDistanceConstraints(Molecule& m,
                       const int * functional_group,
                       SubstructureSearch::SingleSubstructureQuery& query) {

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (functional_group[i] < 0) {
      continue;
    }
    for (int j = 0; j < i - 1; ++j) {
      if (functional_group[j] < 0) {
        continue;
      }
      if (functional_group[i] == functional_group[j]) {
        continue;
      }
      // Atoms in two different functional groups.
      const int d = m.bonds_between(i, j);
      if (d < min_separation || d > max_separation) {
        continue;
      }
      SubstructureSearch::SeparatedAtoms* separated_atoms = query.add_separated_atoms();
      separated_atoms->set_a1(i);
      separated_atoms->set_a2(j);
      if (delta_longer == 0 && delta_shorter == 0) {
        separated_atoms->add_bonds_between(d);
      } else {
        if (d - delta_shorter > 0) {
          separated_atoms->set_min_bonds_between(d - delta_shorter);
        }
        separated_atoms->set_max_bonds_between(d + delta_longer);
      }
    }
  }
}

void
AddDistanceConstraints(Molecule& m,
                       const int * functional_group,
                       int g1,
                       int g2,
                       SubstructureSearch::SingleSubstructureQuery& query) {

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }
      if (functional_group[i] != g1 || functional_group[j] != g2) {
        continue;
      }
      const int d = m.bonds_between(i, j);
      if (d < min_separation || d > max_separation) {
        continue;
      }
      SubstructureSearch::SeparatedAtoms* separated_atoms = query.add_separated_atoms();
      separated_atoms->set_a1(i);
      separated_atoms->set_a2(j);
      if (delta_longer == 0 && delta_shorter == 0) {
        separated_atoms->add_bonds_between(d);
      } else {
        if (d - delta_shorter > 0) {
          separated_atoms->set_min_bonds_between(d - delta_shorter);
        }
        separated_atoms->set_max_bonds_between(d + delta_longer);
      }
    }
  }
}

void
CopyAttributes(Molecule& m,
               atom_number_t zatom,
               SubstructureSearch::SubstructureAtom& query_atom) {
  SubstructureSearch::SubstructureAtomSpecifier* qas = query_atom.add_atom_properties();
  const Atom& atom = m.atom(zatom);
  qas->add_atomic_number(atom.atomic_number());

  if (ncon_becomes_min_ncon) {
    qas->set_min_ncon(atom.ncon());
  } else {
    qas->add_ncon(atom.ncon());
    qas->add_nbonds(m.nbonds(zatom));
  }

  // An open question as to whether ring membership should be allowed to change.
  if (m.is_aromatic(zatom)) {
    qas->set_aromatic(true);
  } else {
    qas->add_ring_bond_count(m.ring_bond_count(zatom));
  }
}

int
AddBonds(Molecule& m,
         const resizable_array<atom_number_t>& atom_order_in_smiles,
         const atom_number_t zatom,
         const int * functional_group,
         int group_number,
         SubstructureSearch::SubstructureAtom& query_atom) {
  const Atom& atom = m.atom(zatom);
  int rc = 0;
  for (const Bond * b : atom) {
    const atom_number_t other = b->other(zatom);
    if (functional_group[other] != group_number) {
      continue;
    }
    if (atom_order_in_smiles[other] > atom_order_in_smiles[zatom]) {
      continue;
    }
//  cerr << "AddBonds:adding bond between " << zatom << " fg " << functional_group[zatom] << " and " << other << " fg " << functional_group[other] << '\n';
    SubstructureSearch::SubstructureBond* query_bond = query_atom.add_query_bond();
    query_bond->set_other_end(other);
    if (b->is_aromatic()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_AROMATIC_BOND);
    } else if (b->is_single_bond()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_SINGLE_BOND);
    } else if (b->is_double_bond()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_DOUBLE_BOND);
    } else if (b->is_triple_bond()) {
      query_bond->add_btype(SubstructureSearch::BondType::SS_TRIPLE_BOND);
    }
    rc++;
  }

  return rc;
}

// Return the atoms for which `functional_group[i] == group_number` ordered
// by `atom_order_in_smiles`.
// natoms is the size of the `functional_group` array.
Set_of_Atoms
AtomsInFunctionalGroup(int natoms,
                       const int * functional_group,
                       int group_number,
                       const resizable_array<atom_number_t>& atom_order_in_smiles) {
  Set_of_Atoms atoms_in_group;
  for (int i = 0; i < natoms; ++i) {
    if (functional_group[i] == group_number) {
      atoms_in_group << i;
    }
  }
  atoms_in_group.iwqsort_lambda([&atom_order_in_smiles](int a1, int a2) {
    if (atom_order_in_smiles[a1] < atom_order_in_smiles[a2])
      return -1;
    if (atom_order_in_smiles[a1] > atom_order_in_smiles[a2])
      return 1;
    return 0;  // Should never happen here.
  });
  
  return atoms_in_group;
}

int
BuildQuery(Molecule& m,
           const int * functional_group,
           int group_number,
           const resizable_array<atom_number_t>& atom_order_in_smiles,
           SubstructureSearch::SingleSubstructureQuery& query) {
  const int matoms = m.natoms();
  const resizable_array<atom_number_t> atoms_in_functional_group = AtomsInFunctionalGroup(matoms, functional_group, group_number, atom_order_in_smiles);
  for (atom_number_t a : atoms_in_functional_group) {
    SubstructureSearch::SubstructureAtom * atom = query.add_query_atom();
    atom->set_id(a);
    CopyAttributes(m, a, *atom);
    AddBonds(m, atom_order_in_smiles, a, functional_group, group_number, *atom);
  }

  return 1;
}

#ifdef NOT_NEEDED_ASDASD
Set_of_Atoms
GetFunctionalGroup(const int * functional_group,
                   int group_number,
                   int n) {
  Set_of_Atoms result;
  for (int i = 0; i < n; ++i) {
    if (functional_group[i] == group_number) {
      result << i;
    }
  }
  return result;
}
#endif

int
Pharacaphore2d(Molecule& m,
               const int * functional_group,
               int number_functional_groups,
               IWString_and_File_Descriptor& output) {
  if (verbose > 1)
    cerr << m.name() << " has " << number_functional_groups << " functional groups\n";

  // Force a smiles computation, and copy the atom order.
  m.smiles();
  resizable_array<atom_number_t> atom_order_in_smiles = m.atom_order_in_smiles();

  SubstructureSearch::SubstructureQuery composite_query;
  composite_query.set_name(m.name().data(), m.name().length());
  SubstructureSearch::SingleSubstructureQuery * query = composite_query.add_query();
  query->set_respect_initial_atom_numbering(true);
  int group_number = 0;
  for (int i = 1; i <= number_functional_groups; ++i, ++group_number) {
    BuildQuery(m, functional_group, i, atom_order_in_smiles, *query);
  }
  
  AddDistanceConstraints(m, functional_group, *query);

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(composite_query, zero_copy_output.get())) {
    cerr << "Pharacaphore2d:cannot write\n";
    return 0;
  }

  return 1;
}

int
IdentifyAtomsToIgnore(Molecule& m,
                      resizable_array_p<Substructure_Query>& atoms_to_ignore,
                      int * ignore_atom) {
  Molecule_to_Match target(&m);
  for (Substructure_Query* q : atoms_to_ignore) {
    Substructure_Results sresults;
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    sresults.each_embedding_set_vector(ignore_atom, 1);
  }

  return 1;  // Not sure what should be returned, ignored for now.
}

int
IdentifyFunctionalGroup(Molecule& m,
                        const atom_number_t zatom,
                        int * fg,
                        int group_number,
                        const int * ignore) {
  fg[zatom] = group_number;
  int rc = 1;
  const Atom& a = m.atom(zatom);
  for (const Bond * b : a) {
    const atom_number_t other = b->other(zatom);
    if (ignore[other]) {
      continue;
    }
    if (fg[other] > 0) {  // Might be our group or another group.
      continue;
    }

    const atomic_number_t zother = m.atomic_number(other);
    if (b->is_single_bond() && zother == 6) {
      continue;
    }
    if (b->is_aromatic() && zother == 6) {
      continue;
    }
    // Always break at an aromatic ring.
    if (b->is_single_bond() && b->nrings() == 0 && (m.is_aromatic(zatom) || m.is_aromatic(other))) {
      continue;
    }
    rc += IdentifyFunctionalGroup(m, other, fg, group_number, ignore);
  }

  return rc;
}

std::tuple<int, int>
IdentifyFunctionalGroups(Molecule& m,
                         int * fg,
                         const int * ignore) {
  const int matoms = m.natoms();

  m.compute_aromaticity_if_needed();

  int atoms_in_functional_groups = 0;
  int group_number = 0;
  for (int i = 0; i < matoms; ++i) {
    if (fg[i] > 0) {
      continue;
    }
    if (ignore[i]) {
      continue;
    }
    if (m.atomic_number(i) == 6) {
      continue;
    }
    group_number++;
    atoms_in_functional_groups += IdentifyFunctionalGroup(m, i, fg, group_number, ignore);
  }

  return {group_number, atoms_in_functional_groups};
}

int
WriteLabelledSmiles(const Molecule & m,
                    const int * functional_group,
                    IWString_and_File_Descriptor& stream_for_labelled_smiles) {
  Molecule mcopy(m);
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (functional_group[i] < 0) {
      continue;
    }
    mcopy.set_atom_map_number(i, functional_group[i]);
  }
  stream_for_labelled_smiles << mcopy.smiles() << ' ' << m.name() << '\n';
  stream_for_labelled_smiles.write_if_buffer_holds_more_than(8192);
  return 1;
}

int
Pharacaphore2d(Molecule& m,
               IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  if (matoms == 0) {
    cerr << "Ignoring empty molecule " << m.name() << '\n';
    return 1;
  }

  std::unique_ptr<int[]> ignore_atoms(new_int(matoms));
  if (atoms_to_ignore.number_elements() > 0) {
    IdentifyAtomsToIgnore(m, atoms_to_ignore, ignore_atoms.get());
  }

  std::unique_ptr<int[]> functional_group(new_int(matoms, -1));
  const auto [number_functional_groups, atoms_in_functional_groups] =
              IdentifyFunctionalGroups(m, functional_group.get(), ignore_atoms.get());
  if (verbose) {
    functional_group_count.extra(number_functional_groups);
    fraction_atoms_in_functional_groups.extra(static_cast<float>(atoms_in_functional_groups) /
                static_cast<float>(matoms));
  }
  if (number_functional_groups == 0) {
    no_functional_groups++;
    return 1;
  }
  
  if (stream_for_labelled_smiles.active()) {
    WriteLabelledSmiles(m, functional_group.get(), stream_for_labelled_smiles);
  }

  if (number_functional_groups == 1) {
    only_one_functional_group++;
    return 1;
  }

  return Pharacaphore2d(m, functional_group.get(), number_functional_groups, output);
}

int
Pharacaphore2d(Molecule& m,
               IWString& output_fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(output_fname)) {
    cerr << "Cannot open " << output_fname << '\n';
    return 0;
  }

  return Pharacaphore2d(m, output);
}

int
Pharacaphore2d (data_source_and_type<Molecule> & input,
                const IWString& output_stem) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    molecules_read++;
    std::unique_ptr<Molecule> free_m(m);
    Preprocess(*m);
    IWString output_fname;
    output_fname << output_stem << molecules_read << ".proto";
    if (! Pharacaphore2d(*m, output_fname)) {
      cerr << "Fatal error processing " << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
Pharacaphore2d(const char * fname,
               FileType input_type,
               const IWString& output_stem) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  return Pharacaphore2d(input, output_stem);
}

int
Pharacaphore2d(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:S:L:d:D:n");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('d')) {
    if (! cl.value('d', min_separation) || min_separation < 1) {
      cerr << "Invalid minimum separation (-d)\n";
      Usage(1);
    }
    if (verbose) {
      cerr << "Will ignore atoms separated by less than " << min_separation << " bonds\n";
    }
  }

  if (cl.option_present('D')) {
    if (! cl.value('D', max_separation) || max_separation < min_separation) {
      cerr << "Invalid maximum separation (-D)\n";
      Usage(1);
    }
    if (verbose) {
      cerr << "Will ignore atoms separated by more than " << max_separation << " bonds\n";
    }
  }

  if (cl.option_present('n')) {
    ncon_becomes_min_ncon = 1;
    if (verbose) {
      cerr << "Ncon values become min_ncon query attributes\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot process -i option\n";
      return 1;
    }
  } else if (! all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 4;
  }

  IWString output_stem;
  if (! cl.option_present('S')) {
    cerr << "Must specify the output stem (-S)\n";
    Usage(1);
  }

  if (cl.option_present('L')) {
    IWString fname = cl.string_value('L');
    if (! fname.ends_with(".smi")) {
      fname << ".smi";
    }
    if (! stream_for_labelled_smiles.open(fname)) {
      cerr << "Cannot open stream for labelled molecules '" << fname << "'\n";
      return 1;
    }
    if (verbose)
      cerr << "Labelled molecules written to '" << fname << "'\n";
  }

  cl.value('S', output_stem);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  for (const auto* fname : cl) {
    if (! Pharacaphore2d(fname, input_type, output_stem)) {
      cerr << "Error processing " << fname << "\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << no_functional_groups << " molecules had no functional groups\n";
    cerr << only_one_functional_group << " molecules had only one functional group\n";
  }

  return 0;
}

}  // namespace pharacaphore_2d

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return pharacaphore_2d::Pharacaphore2d(argc, argv);
}
