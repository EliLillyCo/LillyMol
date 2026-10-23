#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#define SPARSE_FP_CREATOR_IMPLEMENTATION

#include "Molecule_Tools/iwecfp_lib.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/target.h"

namespace iwecfp {

using std::cerr;

static constexpr int kProcessingFinished = 1;
static constexpr int kReadyToProcess = 2;
static constexpr int kNextTime = 3;
static constexpr int kMatchedByStartingAtomQuery = 4;

Iwecfp::Iwecfp() = default;

Iwecfp::~Iwecfp() {
  Finalise();
}

void
Usage(int rc) {
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif

  cerr << R"(Generates EC type fingerprints. Usually invoked by gfp_make.sh.
This is usually invoked from gfp_make.sh.
iwecfp -l -g all -R 3 -P C -J NCEC3 file.smi > file.gfp
 -r <rad>         min shell radius.
 -R <rad>         max shell radius.
 -m               multiplicative formation of bits - loses precision.
 -s               each radius gets its own fingerprint.
 -c               generate R/S chiriality sensitive fingerprint based on isotopic label on input.
 -J <tag>         output tag name, usually NCEC3 for a radius 3 fingerprint.
 -P ...           atom typing specification, enter '-P help' for info.
 -b               all bond type info lost - useful with -P UST:N no atom type.
 -q               query file specifying atoms from which shells are grown
 -M ...           read or write bits found. Enter '-M help' for info.
 -B ...           all bits generated can be written. Enter '-B help' for info.
 -Q ...           options controlling equalisation of times atoms are fingerprinted.
 -X ...           obscure options. Enter '-X help' for info.
 -f               filter existing TDT/fingerprint file.
 -i               input type
 -l               reduce to largest fragment
 -v               verbose output.
)";

  ::exit(rc);
}

void
DisplayDashGOptions(std::ostream& os) {
  os << " -G presence record presence of bits only, not count\n";
  os << " -G WRITE=fname write the global fingerprint to <fname>\n";
  ::exit(0);
}

int
Iwecfp::ReadBitsToInvestigate(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#') || buffer.empty()) {
      continue;
    }

    buffer.truncate_at_first(' ');

    unsigned int b;
    if (!buffer.numeric_value(b)) {
      cerr << "Invalid bit number '" << buffer << "', line " << input.lines_read() << '\n';
      return 0;
    }

    _bits_to_investigate[b] = 0;
  }

  return static_cast<int>(_bits_to_investigate.size());
}

int
Iwecfp::ReadBitsToInvestigate(const char* fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open bits to be examined file name '" << fname << "'\n";
    return 0;
  }

  return ReadBitsToInvestigate(input);
}

int
Iwecfp::CheckAgainstList(Molecule& m, const IWString& smarts,
                         const atom_number_t centre_of_shell,
                         unsigned int sum_so_far, const int radius) {
  auto f = _bits_to_investigate.find(sum_so_far);
  if (f == _bits_to_investigate.end()) {
    return 0;
  }

  ++_bits_found;
  ++(*f).second;

  if (radius > 0) {
    m.set_atom_map_number(centre_of_shell, radius);
  }

  _stream_for_bit_meanings << m.smiles() << ' ' << m.name() << " bit "
                           << sum_so_far << " atom " << centre_of_shell << ' '
                           << smarts << " radius " << radius << '\n';

  if (radius > 0) {
    m.set_atom_map_number(centre_of_shell, 0);
  }

  _stream_for_bit_meanings.write_if_buffer_holds_more_than(32768);
  return 1;
}

// Note that this is quite inefficient since setting the isotope will result
// in destruction of the distance matrix. This is not used frequently, so
// ignore for now.
void
Iwecfp::WriteLabelledSmiles(const Molecule& m, int centre_of_shell, int radius,
                            IWString_and_File_Descriptor& output) {
  Molecule mcopy(m);
  mcopy.recompute_distance_matrix();

  const int matoms = mcopy.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == centre_of_shell) {
      continue;
    }

    const int d = mcopy.bonds_between(centre_of_shell, i);
    if (d > radius) {
      continue;
    }

    if (_write_smiles_with_bit_meanings == 1) {
      mcopy.set_isotope(i, d);
    } else {
      mcopy.set_atom_map_number(i, d);
    }
  }

  output << mcopy.smiles();
}

void
Iwecfp::WriteBit(const int centre_of_shell,
                 const IWString& smarts_for_centre_of_shell,
                 const int radius, unsigned int b, Molecule& m,
                 IWString_and_File_Descriptor& output) {
  _iwdigits_center.append_number(output, centre_of_shell);
  _iwdigits.append_number(output, radius);
  _iwdigits.append_number(output, b);

  constexpr char sep = ' ';
  output << sep << smarts_for_centre_of_shell;

  if (_write_smiles_with_bit_meanings) {
    output << sep;
    WriteLabelledSmiles(m, centre_of_shell, radius, output);
  }

  output << '\n';
  output.write_if_buffer_holds_more_than(32768);
}

int
Iwecfp::BondConstant(const Bond* bond) const {
  if (_all_bonds_same_type) {
    return 1;
  }

  if (bond->is_aromatic()) {
    return 11;
  }

  if (bond->is_triple_bond()) {
    return 7;
  }

  if (bond->is_double_bond()) {
    return 5;
  }

  return 3;
}

void
Iwecfp::Increment(unsigned int& sum_so_far, const int bc,
                  const atype_t atom_constant) const {
  if (_additive) {
    sum_so_far += bc * atom_constant;
  } else {
    sum_so_far *= bc * atom_constant;
  }
}

int
Iwecfp::GenerateShells(const int matoms, int radius, const int max_radius,
                       const Atom* const* atoms, const atype_t* atom_constant,
                       int* processing_status, unsigned int sum_so_far,
                       Molecule& m, Sparse_Fingerprint_Creator* sfc) {
  ++radius;

  if (_additive) {
    sum_so_far *= 7879;
  }

  const int add_tails_here = (_add_tails > 0 && radius <= _add_tails);

  for (int i = 0; i < matoms; ++i) {
    if (kReadyToProcess != processing_status[i]) {
      continue;
    }

    const Atom* ai = atoms[i];
    const int acon = ai->ncon();
    const Bond* const* bonds = ai->rawdata();

    for (int j = 0; j < acon; ++j) {
      const Bond* b = bonds[j];
      atom_number_t k = b->other(i);

      if (kProcessingFinished == processing_status[k]) {
        const int bc = BondConstant(b);
        Increment(sum_so_far, bc, atom_constant[i]);
        if (add_tails_here) {
          sfc->hit_bit(sum_so_far);
        }
      } else if (kReadyToProcess == processing_status[k]) {
        ;
      } else {
        processing_status[k] = kNextTime;
      }
    }
  }

  if (radius >= _min_shell_radius) {
    sfc->hit_bit(sum_so_far);
    if (_centre_atom_isotope) {
      sfc->hit_bit(sum_so_far + _centre_atom_isotope);
    }

    if (_looking_for_bit_meanings) {
      CheckAgainstList(*_current_molecule, _smarts_for_centre_of_shell,
                       _centre_of_shell, sum_so_far, radius);
    } else if (_stream_for_all_bits.is_open()) {
      WriteBit(_centre_of_shell, _smarts_for_centre_of_shell, radius,
               sum_so_far, m, _stream_for_all_bits);
    }
  }

  if (max_radius > 0 && radius >= max_radius) {
    return 1;
  }

  int continue_processing = 0;
  for (int i = 0; i < matoms; ++i) {
    if (kReadyToProcess == processing_status[i]) {
      processing_status[i] = kProcessingFinished;
    } else if (kNextTime == processing_status[i]) {
      processing_status[i] = kReadyToProcess;
      continue_processing = 1;
    }
  }

  if (!continue_processing) {
    return 1;
  }

  if (_each_shell_gets_own_fingerprint) {
    return GenerateShells(matoms, radius, max_radius, atoms, atom_constant,
                          processing_status, sum_so_far, m, sfc + 1);
  } else {
    return GenerateShells(matoms, radius, max_radius, atoms, atom_constant,
                          processing_status, sum_so_far, m, sfc);
  }
}

int
Iwecfp::GenerateShellsSubset(const int matoms, int radius, const int max_radius,
                       const Atom* const* atoms, const atype_t* atom_constant,
                       const int* include_atom, int* processing_status,
                       unsigned int sum_so_far, Molecule& m,
                       Sparse_Fingerprint_Creator* sfc) {
  ++radius;

  if (_additive) {
    sum_so_far *= 7879;
  }

  const int add_tails_here = (_add_tails > 0 && radius <= _add_tails);

  for (int i = 0; i < matoms; ++i) {
    if (kReadyToProcess != processing_status[i]) {
      continue;
    }

    const Atom* ai = atoms[i];
    const int acon = ai->ncon();
    const Bond* const* bonds = ai->rawdata();

    for (int j = 0; j < acon; ++j) {
      const Bond* b = bonds[j];
      atom_number_t k = b->other(i);
      if (!include_atom[k]) {
        continue;
      }

      if (kProcessingFinished == processing_status[k]) {
        const int bc = BondConstant(b);
        Increment(sum_so_far, bc, atom_constant[i]);
        if (add_tails_here) {
          sfc->hit_bit(sum_so_far);
        }
      } else if (kReadyToProcess == processing_status[k]) {
        ;
      } else {
        processing_status[k] = kNextTime;
      }
    }
  }

  if (radius >= _min_shell_radius) {
    sfc->hit_bit(sum_so_far);
    if (_centre_atom_isotope) {
      sfc->hit_bit(sum_so_far + _centre_atom_isotope);
    }

    if (_looking_for_bit_meanings) {
      CheckAgainstList(*_current_molecule, _smarts_for_centre_of_shell,
                       _centre_of_shell, sum_so_far, radius);
    } else if (_stream_for_all_bits.is_open()) {
      WriteBit(_centre_of_shell, _smarts_for_centre_of_shell, radius,
               sum_so_far, m, _stream_for_all_bits);
    }
  }

  if (max_radius > 0 && radius >= max_radius) {
    return 1;
  }

  int continue_processing = 0;
  for (int i = 0; i < matoms; ++i) {
    if (kReadyToProcess == processing_status[i]) {
      processing_status[i] = kProcessingFinished;
    } else if (kNextTime == processing_status[i]) {
      processing_status[i] = kReadyToProcess;
      continue_processing = 1;
    }
  }

  if (!continue_processing) {
    return 1;
  }

  if (_each_shell_gets_own_fingerprint) {
    return GenerateShellsSubset(matoms, radius, max_radius, atoms, atom_constant,
                          include_atom, processing_status, sum_so_far, m, sfc + 1);
  } else {
    return GenerateShellsSubset(matoms, radius, max_radius, atoms, atom_constant,
                          include_atom, processing_status, sum_so_far, m, sfc);
  }
}

void
Iwecfp::FormBit(Molecule& m, const atype_t* atom_constant,
                const Atom* const* atoms, const atom_number_t zatom, int max_r,
                int* processing_status, Sparse_Fingerprint_Creator* sfc) {
  const int matoms = m.natoms();

  std::fill_n(processing_status, matoms, 0);
  processing_status[zatom] = kProcessingFinished;

  const auto e = atom_constant[zatom];
  sfc[0].hit_bit(e);

  const Atom* a = atoms[zatom];
  const int acon = a->ncon();

  if (max_r == 0) {
    return;
  }

  for (int i = 0; i < acon; ++i) {
    const Bond* b = a->item(i);
    const auto j = b->other(zatom);
    processing_status[j] = kReadyToProcess;
  }

  if (_each_shell_gets_own_fingerprint) {
    GenerateShells(matoms, 0, _max_shell_radius, atoms, atom_constant,
                   processing_status, e, m, sfc + 1);
  } else {
    GenerateShells(matoms, 0, _max_shell_radius, atoms, atom_constant,
                   processing_status, e, m, sfc);
  }
}

void
Iwecfp::FormBitSubset(Molecule& m, const atype_t* atom_constant,
                const Atom* const* atoms, const int* include_atom,
                const atom_number_t zatom, int max_r, int* processing_status,
                Sparse_Fingerprint_Creator* sfc) {
  const int matoms = m.natoms();

  std::fill_n(processing_status, matoms, 0);
  processing_status[zatom] = kProcessingFinished;

  const auto e = atom_constant[zatom];
  sfc[0].hit_bit(e);

  const Atom* a = atoms[zatom];
  const int acon = a->ncon();

  if (max_r == 0) {
    return;
  }

  for (int i = 0; i < acon; ++i) {
    const Bond* b = a->item(i);
    const auto j = b->other(zatom);
    if (include_atom[j]) {
      processing_status[j] = kReadyToProcess;
    }
  }

  if (_each_shell_gets_own_fingerprint) {
    GenerateShellsSubset(matoms, 0, _max_shell_radius, atoms, atom_constant,
                   include_atom, processing_status, e, m, sfc + 1);
  } else {
    GenerateShellsSubset(matoms, 0, _max_shell_radius, atoms, atom_constant,
                   include_atom, processing_status, e, m, sfc);
  }
}

void
Iwecfp::IdentifyAtomsWithinRange(Molecule& m, Set_of_Atoms* atoms_within_range) {
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    Set_of_Atoms* ai = atoms_within_range + i * (_max_shell_radius + 1);
    for (int r = 0; r <= _max_shell_radius; ++r) {
      ai[r].add(i);
    }
  }

  for (int i = 0; i < matoms; ++i) {
    Set_of_Atoms* ai = atoms_within_range + i * (_max_shell_radius + 1);
    for (int j = i + 1; j < matoms; ++j) {
      const int d = m.bonds_between(i, j);
      if (d > _max_shell_radius) {
        continue;
      }

      Set_of_Atoms* aj = atoms_within_range + j * (_max_shell_radius + 1);
      for (int r = 1; r <= _max_shell_radius; ++r) {
        if (d <= r) {
          ai[r].add(j);
          aj[r].add(i);
        }
      }
    }
  }
}

static int
OurOwnCustomVariance(const int matoms, const int sumv, const int sumv2) {
  const int rc = matoms * sumv2 - (sumv * sumv);
  if (rc >= 0) {
    return rc;
  }

  cerr << "Negative variance: matoms " << matoms << " sumv " << sumv
       << " sumv2 " << sumv2 << '\n';
  return std::numeric_limits<int>::max();
}

static int
ComputeVariance(const Set_of_Atoms& s, const int* visited, const int matoms,
                int sumv, int sumv2) {
  const int n = s.number_elements();

  for (int i = 0; i < n; ++i) {
    const atom_number_t j = s[i];
    const int v = visited[j];

    sumv += 1;
    sumv2 = sumv2 - (v * v) + (v + 1) * (v + 1);
  }

  return OurOwnCustomVariance(matoms, sumv, sumv2);
}

static void
RecomputeSums(const int* visited, const int matoms, int& sumv, int& sumv2) {
  sumv = 0;
  sumv2 = 0;

  for (int i = 0; i < matoms; ++i) {
    sumv += visited[i];
    sumv2 += visited[i] * visited[i];
  }
}

int
Iwecfp::DoEqualiseAtomCoverage(Molecule& m, const Atom* const* atoms,
                               const atype_t* atom_constant,
                               int* processing_status,
                               Sparse_Fingerprint_Creator* sfc) {
  const int matoms = m.natoms();

  Set_of_Atoms* atoms_within_range =
      new Set_of_Atoms[matoms * (_max_shell_radius + 1)];
  std::unique_ptr<Set_of_Atoms[]> free_atoms_within_range(atoms_within_range);

  IdentifyAtomsWithinRange(m, atoms_within_range);

  int* visited = new_int(matoms, 1);
  std::unique_ptr<int[]> free_visited(visited);

  for (int i = 0; i < matoms; ++i) {
    for (int j = 0; j < matoms; ++j) {
      if (j == i) {
        continue;
      }

      const int d = m.bonds_between(i, j);
      if (d > _max_shell_radius) {
        continue;
      }

      visited[j] += _max_shell_radius - d + 1;
    }
  }

  int max_visited = visited[0];
  int min_visited = visited[0];

  for (int i = 0; i < matoms; ++i) {
    if (visited[i] > max_visited) {
      max_visited = visited[i];
    } else if (visited[i] < min_visited) {
      min_visited = visited[i];
    }
  }

  if ((static_cast<float>(max_visited - min_visited) /
       static_cast<float>(max_visited)) > 0.8f) {
    return 1;
  }

  int sumv = 0;
  int sumv2 = 0;

  const int narv = matoms * (_max_shell_radius + 1);
  int* var = new int[narv];
  std::unique_ptr<int[]> free_var(var);

  RecomputeSums(visited, matoms, sumv, sumv2);

  if (_label_by_visited && _verbose > 1) {
    for (int i = 0; i < matoms; ++i) {
      m.set_atom_map_number(i, visited[i]);
    }
    cerr << m.smiles() << ' ' << m.name() << " before equalisation, var "
         << OurOwnCustomVariance(matoms, sumv, sumv2) << '\n';
  }

  for (int i = 0; i < _equalise_atom_coverage; ++i) {
    for (int j = 0; j < matoms; ++j) {
      for (int r = 0; r <= _max_shell_radius; ++r) {
        const Set_of_Atoms& s =
            atoms_within_range[j * (_max_shell_radius + 1) + r];
        var[j * (_max_shell_radius + 1) + r] =
            ComputeVariance(s, visited, matoms, sumv, sumv2);
      }
    }

    int min_variance = var[0];
    resizable_array<int> mins;

    for (int j = 1; j < narv; ++j) {
      if (j % (_max_shell_radius + 1) == 0) {
        continue;
      }

      if (var[j] > min_variance) {
        continue;
      }

      if (var[j] < min_variance) {
        mins.resize_keep_storage(0);
        mins.add(j);
        min_variance = var[j];
      } else {
        mins.add(j);
      }
    }

    for (int j = 0; j < mins.number_elements(); ++j) {
      const Set_of_Atoms& s = atoms_within_range[mins[j]];
      s.increment_vector(visited, 1);
      FormBit(m, atom_constant, atoms, s[0], s.number_elements() - 1,
              processing_status, sfc);
    }

    RecomputeSums(visited, matoms, sumv, sumv2);
  }

  if (_label_by_visited) {
    for (int i = 0; i < matoms; ++i) {
      m.set_atom_map_number(i, visited[i]);
    }

    if (_verbose > 1) {
      cerr << m.smiles() << ' ' << m.name() << " after equalisation, var "
           << OurOwnCustomVariance(matoms, sumv, sumv2) << '\n';
    }
  }

  return 1;
}

int
Iwecfp::IdentifyStartAtoms(Molecule& m, int* processing_status, int matched_flag) {
  Molecule_to_Match target(&m);

  int rc = 0;
  for (Substructure_Query* q : _start_atom_query) {
    Substructure_Results sresults;
    if (!q->substructure_search(target, sresults)) {
      continue;
    }

    sresults.each_embedding_set_vector(processing_status, matched_flag);
    ++rc;
  }

  return rc;
}


FingerprintResult
Iwecfp::Fingerprint(Molecule& m, const atype_t* atom_constant,
                const int* include_atom, Sparse_Fingerprint_Creator* sfc) {
  if (include_atom == nullptr) {
    return Fingerprint(m, atom_constant, sfc);
  }

  if (_equalise_atom_coverage) {
    cerr << "Iwecfp::Fingerprint:subset fingerprints do not support atom coverage equalisation\n";
    return FingerprintResult::kFatal;
  }

  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  Atom** atoms = new Atom*[matoms];
  std::unique_ptr<Atom*[]> free_atoms(atoms);
  m.atoms(const_cast<const Atom**>(atoms));

  int set_centre_atom_member_variables = 0;
  if (_looking_for_bit_meanings || _stream_for_all_bits.is_open()) {
    set_centre_atom_member_variables = 1;
    _current_molecule = &m;
    if (_stream_for_all_bits.is_open()) {
      _stream_for_all_bits << m.name() << '\n';
    }
  }

  if (_central_atom_possible_chiral) {
    set_centre_atom_member_variables = 1;
  }

  int* processing_status = new int[matoms];
  std::unique_ptr<int[]> free_processing_status(processing_status);

  std::unique_ptr<int[]> start_atom_query_match;
  if (_start_atom_query.size() > 0) {
    start_atom_query_match = std::make_unique<int[]>(matoms);

    if (!IdentifyStartAtoms(m, start_atom_query_match.get(), 1)) {
      if (_verbose) {
        cerr << "Iwecfp:no start atoms defined for " << _start_atom_query.size() << " queries\n";
      }
      return FingerprintResult::kNoStartAtoms;
    }
  }

  int shell_centres = 0;
  for (int i = 0; i < matoms; ++i) {
    if (!include_atom[i]) {
      continue;
    }

    if (start_atom_query_match && !start_atom_query_match[i]) {
      continue;
    }

    ++shell_centres;

    const atype_t e = atom_constant[i];

    if (set_centre_atom_member_variables) {
      _smarts_for_centre_of_shell = m.smarts_equivalent_for_atom(i);
      _centre_of_shell = i;
      _centre_atom_isotope = m.isotope(i);
    } else {
      _centre_atom_isotope = 0;
    }

    if (_min_shell_radius == 0) {
      sfc[0].hit_bit(e);
    }

    if (_looking_for_bit_meanings) {
      CheckAgainstList(m, _smarts_for_centre_of_shell, i, e, 0);
    } else if (_stream_for_all_bits.is_open()) {
      WriteBit(i, _smarts_for_centre_of_shell, 0, e, m, _stream_for_all_bits);
    }

    std::fill_n(processing_status, matoms, 0);
    processing_status[i] = kProcessingFinished;

    const Atom* ai = atoms[i];
    const int acon = ai->ncon();
    for (int j = 0; j < acon; ++j) {
      atom_number_t k = ai->other(i, j);
      if (include_atom[k]) {
        processing_status[k] = kReadyToProcess;
      }
    }

    if (_each_shell_gets_own_fingerprint) {
      GenerateShellsSubset(matoms, 0, _max_shell_radius, const_cast<const Atom* const*>(atoms),
                     atom_constant, include_atom, processing_status, e, m, sfc + 1);
    } else {
      GenerateShellsSubset(matoms, 0, _max_shell_radius, const_cast<const Atom* const*>(atoms),
                     atom_constant, include_atom, processing_status, e, m, sfc);
    }
  }

  if (_stream_for_all_bits.is_open()) {
    _stream_for_all_bits << "|\n";
  }

  if (shell_centres == 0) {
    return FingerprintResult::kNoStartAtoms;
  }

  return FingerprintResult::kOk;
}

FingerprintResult
Iwecfp::Fingerprint(Molecule& m, const atype_t* atom_constant,
                Sparse_Fingerprint_Creator* sfc) {
  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  Atom** atoms = new Atom*[matoms];
  std::unique_ptr<Atom*[]> free_atoms(atoms);
  m.atoms(const_cast<const Atom**>(atoms));

  int set_centre_atom_member_variables = 0;
  if (_looking_for_bit_meanings || _stream_for_all_bits.is_open()) {
    set_centre_atom_member_variables = 1;
    _current_molecule = &m;
    if (_stream_for_all_bits.is_open()) {
      _stream_for_all_bits << m.name() << '\n';
    }
  }

  if (_central_atom_possible_chiral) {
    set_centre_atom_member_variables = 1;
  }

  int* processing_status = new int[matoms];
  std::unique_ptr<int[]> free_processing_status(processing_status);

  std::unique_ptr<int[]> start_atom_query_match;
  if (_start_atom_query.size() > 0) {
    start_atom_query_match = std::make_unique<int[]>(matoms);

    if (!IdentifyStartAtoms(m, start_atom_query_match.get(), 1)) {
      if (_verbose) {
        cerr << "Iwecfp:no start atoms defined for " << _start_atom_query.size() << " queries\n";
      }
      return FingerprintResult::kNoStartAtoms;
    }
  }

  for (int i = 0; i < matoms; ++i) {
    if (start_atom_query_match && !start_atom_query_match[i]) {
      continue;
    }

    const atype_t e = atom_constant[i];

    if (set_centre_atom_member_variables) {
      _smarts_for_centre_of_shell = m.smarts_equivalent_for_atom(i);
      _centre_of_shell = i;
      _centre_atom_isotope = m.isotope(i);
    } else {
      _centre_atom_isotope = 0;
    }

    if (_min_shell_radius == 0) {
      sfc[0].hit_bit(e);
    }

    if (_looking_for_bit_meanings) {
      CheckAgainstList(m, _smarts_for_centre_of_shell, i, e, 0);
    } else if (_stream_for_all_bits.is_open()) {
      WriteBit(i, _smarts_for_centre_of_shell, 0, e, m, _stream_for_all_bits);
    }

    std::fill_n(processing_status, matoms, 0);
    processing_status[i] = kProcessingFinished;

    const Atom* ai = atoms[i];
    const int acon = ai->ncon();
    for (int j = 0; j < acon; ++j) {
      atom_number_t k = ai->other(i, j);
      processing_status[k] = kReadyToProcess;
    }

    if (_each_shell_gets_own_fingerprint) {
      GenerateShells(matoms, 0, _max_shell_radius, const_cast<const Atom* const*>(atoms),
                     atom_constant, processing_status, e, m, sfc + 1);
    } else {
      GenerateShells(matoms, 0, _max_shell_radius, const_cast<const Atom* const*>(atoms),
                     atom_constant, processing_status, e, m, sfc);
    }
  }

  if (_equalise_atom_coverage) {
    DoEqualiseAtomCoverage(m, const_cast<const Atom* const*>(atoms), atom_constant,
                           processing_status, sfc);
  }

  if (_stream_for_all_bits.is_open()) {
    _stream_for_all_bits << "|\n";
  }

  return FingerprintResult::kOk;
}

void
DisplayDashMOptions(int rc) {
  cerr << R"(The -M option specifies files of bit numbers to either write or read.
 -M WRITE=<fname> write bits found to <fname>.
 -M READ=<fname>  identify bits to be examined - one per line.
)";
  ::exit(rc);
}

void
DisplayDashBOptions(int rc) {
  cerr << R"(Writes info about all bits generated.
 -B write         info on all bits produced to <fname> (large!)
 -B smiles        include isotopically labelled smiles in the -B file
 -B smilesm       include atom map number labelled smiles in the -B file
)";

  ::exit(rc);
}

void
DisplayDashXOptions(int rc) {
  cerr << R"(The following -X qualifiers are recognised.
 -X addtails            At each shell, single atom "tails" are fingerprinted.
)";

  ::exit(rc);
}

int
Iwecfp::Initialise(Command_Line& cl) {
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  _verbose = cl.option_count('v');

  if (cl.option_present('r')) {
    if (!cl.value('r', _min_shell_radius) || _min_shell_radius < 0) {
      cerr << "The min shell radius (-r) must be a whole +ve number\n";
      Usage(2);
    }

    if (_verbose) {
      cerr << "Will only fingerprint paths larger than " << _min_shell_radius << " bonds\n";
    }
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', _max_shell_radius) || _max_shell_radius < 0) {
      cerr << "The max shell radius (-R) must be a whole +ve number\n";
      Usage(2);
    }

    if (_max_shell_radius < _min_shell_radius) {
      cerr << "Inconsistent min " << _min_shell_radius << " and max "
           << _max_shell_radius << " shell radius\n";
      Usage(6);
    }

    if (_verbose) {
      cerr << "Max radius " << _max_shell_radius << '\n';
    }
  }

  if (cl.option_present('s')) {
    if (!cl.option_present('R')) {
      cerr << "Sorry, the -s option only works when the -R option is specified\n";
      Usage(4);
    }

    _each_shell_gets_own_fingerprint = 1;
    if (_verbose) {
      cerr << "Each radius will get its own fingerprint\n";
    }
  }

  if (cl.option_present('b')) {
    _all_bonds_same_type = 1;
    if (_verbose) {
      cerr << "All bonds considered identical\n";
    }
  }

  if (cl.option_present('m')) {
    _additive = 0;
    if (_verbose) {
      cerr << "Fingerprints formed with multiplication operations\n";
    }
  }

  if (cl.option_present('c')) {
    _central_atom_possible_chiral = 1;
    if (_verbose) {
      cerr << "Central atoms will only be CD4 or CD3H types\n";
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, _start_atom_query, _verbose, 'q')) {
      cerr << "iwecfp:cannot process start atom queries (-q)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _start_atom_query.size() << " start atom queries\n";
    }
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x.starts_with("addtails=")) {
        x.remove_leading_chars(9);
        if (! x.numeric_value(_add_tails) || _add_tails < 1) {
          cerr << "The addtails directive must be a whole +ve number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will fingerprint tails out to radius " << _add_tails << '\n';
        }
      } else if (x == "help") {
        DisplayDashXOptions(0);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(1);
      }
    }
  }

  if (cl.option_present('M')) {
    IWString read_fname, write_fname;
    int i = 0;
    const_IWSubstring m;

    while (cl.value('M', m, i++)) {
      if (m.starts_with("READ=")) {
        m.remove_leading_chars(5);
        read_fname = m;
      } else if (m.starts_with("WRITE=")) {
        m.remove_leading_chars(6);
        write_fname = m;
      } else if (m == "help") {
        DisplayDashMOptions(0);
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        DisplayDashMOptions(1);
      }
    }

    if (read_fname.empty() || write_fname.empty()) {
      cerr << "For identifying features, must have both read and write file names (-M)\n";
      Usage(3);
    }

    if (read_fname == write_fname) {
      cerr << "Cannot use same file for reading and writing (-M)\n";
      Usage(2);
    }

    if (!ReadBitsToInvestigate(read_fname.null_terminated_chars())) {
      cerr << "Cannot read bits to investigate from '" << read_fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _bits_to_investigate.size()
           << " bits to investigate from '" << read_fname << "'\n";
    }

    if (!_stream_for_bit_meanings.open(write_fname.null_terminated_chars())) {
      cerr << "Cannot open stream for bit interpretations '" << write_fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Bit interpretations written to '" << write_fname << "'\n";
    }

    _looking_for_bit_meanings = 1;
  }

  if (cl.option_present('B')) {
    const_IWSubstring b;
    IWString fname;

    for (int i = 0; cl.value('B', b, i); ++i) {
      if (b == "smiles") {
        _write_smiles_with_bit_meanings = 1;
      } else if (b == "smilesm") {
        _write_smiles_with_bit_meanings = 2;
      } else if (b == "help") {
        DisplayDashBOptions(0);
      } else {
        fname = b;
      }
    }

    if (fname.empty()) {
      cerr << "No file for bit meanings (-B)\n";
      return 0;
    }

    if (!_stream_for_all_bits.open(fname)) {
      cerr << "Cannot open stream for all bits '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Info on all bits written to '" << fname << "'\n";
    }

    _iwdigits_center.initialise(100);
    _iwdigits.set_include_leading_space(1);
    _iwdigits.initialise(100);
  }

  if (cl.option_present('Q')) {
    const_IWSubstring q;
    for (int i = 0; cl.value('Q', q, i); ++i) {
      if (q == "label") {
        _label_by_visited = 1;
      } else if (q == "help") {
        // Original source had no detailed help here.
      } else if (!q.numeric_value(_equalise_atom_coverage) ||
                 _equalise_atom_coverage < 1) {
        cerr << "The number of visited equalisation steps (-Q) must be a whole +ve number\n";
        Usage(1);
      }
    }
  }

  return 1;
}

int
Iwecfp::Report(std::ostream& output) const {
  if (_looking_for_bit_meanings) {
    output << "Found " << _bits_found << " bits in lookup file\n";
    for (const auto& [bit, count] : _bits_to_investigate) {
      output << "Found " << count << " instances of " << bit << '\n';
    }
  }

  return output.good();
}

int
Iwecfp::Finalise() {
  return 1;
}

}  // namespace iwecfp
