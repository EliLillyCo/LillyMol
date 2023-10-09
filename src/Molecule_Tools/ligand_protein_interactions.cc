// Implemention of Niraj's ideas around specific ligand/protein interactions

#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/ligand_protein_interactions.pb.h"

namespace ligand_protein {

using std::cerr;
using std::setw;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << " -P <fname>       pdb file containing protein\n";
  cerr << " -S <stem>        file name step for output file(s)\n";
  cerr << " -v               verbose output\n";
  ::exit(rc);
}

// Within the protein we store each of the rings and the center atom.
struct RingAndCenter {
    Set_of_Atoms atoms;
    Space_Vector<float> center;
    // This will be known for the ligand, and inferred in the protein.
    // Initialised to -1 to indicate not initialised.
    int aromatic = -1;

  public:
    int Initialise(const Molecule& m, const Ring* r);
};

// Copy the atoms in `r` to this->atoms and compute the center.
int
RingAndCenter::Initialise(const Molecule& m, const Ring* r) {
  atoms = *r;

  uint32_t rsize = r->size();
  for (uint32_t i = 0; i < rsize; ++i) {
    center += *m.atomi(r->item(i));
  }

  if (r->is_aromatic()) {
    aromatic = 1;
  }

  center /= static_cast<float>(rsize);
  cerr << "Center initialised to " << center << ' ' << atoms << '\n';

  return 1;
}

// We cannot instantiate a distance matrix for the protein, so in order
// to identify atoms in range, we do a search each time.
// This class does these radius-defined collecting of nearby atoms.
class IdentifyConnectedAtoms {
  private:
    Molecule* _m;

    int _max_radius;

    int * _visited;

  // Private functions
    void ConnectedAtoms(atom_number_t zatom, int radius,
                                       uint32_t max_to_find,
                                       Set_of_Atoms& result);

  public:
    IdentifyConnectedAtoms();
    ~IdentifyConnectedAtoms();

    int Initialise(Molecule* m);

    void set_max_radius(int s) {
      _max_radius = s;
    }

    Set_of_Atoms ConnectedAtoms(atom_number_t zatom,
                                uint32_t max_to_find = std::numeric_limits<uint32_t>::max());
};

// Performs the interaction computation.
class Options {
  private:
    int _verbose = 0;

    FileType _input_type;

    // Data on both the ligand and the protein.
    // Perhaps it would make sense to create a struct to hold
    // this common information.
    Molecule * _ligand;
    uint32_t _atoms_in_ligand;
    int* _ligand_ring_system;
    // The number of ring systems.
    int _rings_in_ligand;
    RingAndCenter* _ligand_rings;

    Molecule * _protein;
    uint32_t _atoms_in_protein;
    int* _protein_ring_system;
    // The number of ring systems.
    int _rings_in_protein;
    RingAndCenter* _protein_rings;

    // A scratch array sizes to _protein->natoms();
    int* _protein_scratch;

    // The distance between each atom in _ligand and each atom in _protein.
    float* _distance;

    // For each atom in the ligand, the protein atoms that are in range.
    resizable_array<uint32_t>* _neighbours;

    //Parsed from the -C option.
    LigandProtein::InteractionConfig _config;

    IdentifyConnectedAtoms _id_connected_atoms;

    // Variables that keep track of what we find.
    int _acceptor_donor_found;
    int _donor_acceptor_found;
    int _pipi_found;
    int _donor_to_pi_found;

    // File name stem for outputs, from the -S option.
    IWString _stem;
    int _file_name_index = 0;

  // private functions

    IWString ResidueName(atom_number_t zatom) const;
    int ReadProtein(const char * fname);
    int ReadLigand(const char * fname);
    // Operates on either the ligand or the protein.
    int IdentifyRingSystems(Molecule& m, int* ring_system);
    int InitialiseDistanceMatrix();
    float DistanceBetweenLPAtoms(int ligand_atom, int protein_atom) const;
    int IdentifyHBonds();
    int IdentifyHBonds(atom_number_t zatom);
    int IdentifyHBondsAcceptorInLigand(atom_number_t zatom);
    int IdentifyHBondsDonorInLigand(atom_number_t zatom, const Set_of_Atoms& hydrogens_attached);
    int IdentifyPiPi();

    int OkHBondDistance(distance_t distance) const;
    int OkHBondAngle(distance_t distance) const;
    int OkHPiBond(const Molecule& m, atom_number_t hydrogen, const RingAndCenter& rc);
    int OkHPiAngle(const Molecule& m1,
                   atom_number_t zatom,
                   const Molecule& m2,
                   const RingAndCenter& rc) const;
    int OkPiPiDistance(const Space_Vector<float>& p1, const Space_Vector<float>& p2) const;
    int OkPiPiAngle(const Molecule& m1,
                     const RingAndCenter& rc1,
                     const Molecule& m2,
                     const RingAndCenter& rc2) const;

    int InitialiseRingsInLigand();
    int InitialiseRingsInProtein();
    int IdentifyAcceptorsInProtein(atom_number_t h);
    int IdentifyDonorsInProtein(atom_number_t zatom);
    int ReconnectHydrogens(Molecule& m) const;
    int ReconnectHydrogen(Molecule& m, atom_number_t zatom) const;
    int MakeBondsFromDistances(Molecule& m) const;
    int AddProximalAtoms(const Molecule& m,
                          int* xref,
                          atom_number_t zatom,
                          Molecule& destination);

    IWString MakePiPiName(const RingAndCenter& rci, const RingAndCenter& rcj) const;
    IWString MakeHPiName(atom_number_t h, const RingAndCenter& rcj) const;
    IWString MakeHAccName(atom_number_t h, atom_number_t acc) const;
    IWString MakeAccDonName(atom_number_t acc, atom_number_t don) const;

    int WriteLigand(Molecule& m, IWString fname);

    int Write(Molecule& m1,  // not const because the distance matrix is used.
               const Set_of_Atoms& s1,
               const int* ring_system1,
               const Molecule& m2,
               const Set_of_Atoms& s2,
               const int* ring_system2,
               const IWString& name);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Process(IWString_and_File_Descriptor& output);

    int verbose() const {
      return _verbose;
    }

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;

  _protein = nullptr;
  _atoms_in_protein = 0;
  _protein_ring_system = nullptr;

  _protein_scratch = nullptr;

  _ligand = nullptr;
  _atoms_in_ligand = 0;
  _ligand_ring_system = nullptr;

  _distance = nullptr;
  _neighbours = nullptr;

  _rings_in_protein = 0;
  _protein_rings = nullptr;

  _rings_in_ligand = 0;
  _ligand_rings = nullptr;

  _acceptor_donor_found = 0;
  _donor_acceptor_found = 0;
  _pipi_found = 0;
  _donor_to_pi_found = 0;
}

Options::~Options() {
  if (_protein != nullptr) {
    delete _protein;
  }
  if (_ligand != nullptr) {
    delete _ligand;
  }
  if (_ligand_ring_system != nullptr) {
    delete [] _ligand_ring_system;
  }
  if (_protein_ring_system != nullptr) {
    delete [] _protein_ring_system;
  }
  if (_distance != nullptr) {
    delete [] _distance;
  }
  if (_neighbours != nullptr) {
    delete [] _neighbours;
  }
  if (_protein_rings != nullptr) {
    delete [] _protein_rings;
  }
  if (_protein_scratch != nullptr) {
    delete [] _protein_scratch;
  }
}

float
Options::DistanceBetweenLPAtoms(int ligand_atom, int protein_atom) const {
  const uint32_t ndx = ligand_atom * _atoms_in_protein + protein_atom;
  return _distance[ndx];
}

IdentifyConnectedAtoms::IdentifyConnectedAtoms() {
  _m = nullptr;
  _visited = nullptr;
  _max_radius = 3;
}

int
IdentifyConnectedAtoms::Initialise(Molecule* m) {
  _m = m;
  _visited = new int[m->natoms()];
  return 1;
}

IdentifyConnectedAtoms::~IdentifyConnectedAtoms() {
  delete [] _visited;
}

Set_of_Atoms
IdentifyConnectedAtoms::ConnectedAtoms(atom_number_t zatom, uint32_t max_to_find) {
  Set_of_Atoms result;
  std::fill_n(_visited, _m->natoms(), 0);
  result.resize_keep_storage(0);
  ConnectedAtoms(zatom, 0, max_to_find, result);
  return result;
}

void
IdentifyConnectedAtoms::ConnectedAtoms(atom_number_t zatom,
                                       int radius,
                                       uint32_t max_to_find,
                                       Set_of_Atoms& result) {
  _visited[zatom] = 1;
  result << zatom;
  if (result.size() >= max_to_find) {
    return;
  }
  if (radius == _max_radius) {
    return;
  }

  // Do a breadth first search.
  Set_of_Atoms next_shell;
  const Atom* a = _m->atomi(zatom);
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (_visited[j]) {
      continue;
    }
    next_shell << j;
  }

  if (next_shell.empty()) {
    return;
  }

  for (atom_number_t a : next_shell) {
    // If there are rings present...
    if (_visited[zatom]) {
      return;
    }

    ConnectedAtoms(a, radius + 1, max_to_find, result);
    if (result.size() > max_to_find) {
      return;
    }
  }
}

int
Options::InitialiseDistanceMatrix() {
  _distance = new float[_atoms_in_ligand * _atoms_in_protein];
  _neighbours = new resizable_array<uint32_t>[_atoms_in_ligand];

  Accumulator<double> acc;

  for (uint32_t i = 0; i < _atoms_in_ligand; ++i) {
    const Atom* ligand_atom = _ligand->atomi(i);
    for (uint32_t j = 0; j < _atoms_in_protein; ++j) {
      const Atom* protein_atom = _protein->atomi(j);
      // Do we want disconnected atoms or not?
      if (protein_atom->ncon() == 0) {
        continue;
      }

      const float d = ligand_atom->distance(*protein_atom);
      acc.extra(d);
      _distance[i * _atoms_in_protein + j] = d;

      // cerr << "Protein atom in " << _protein->atoms_in_fragment_containing_atom(j) << " fragment\n";
      if (_protein->atoms_in_fragment_containing_atom(j) < 5) {
        continue;
      }

      if (d < _config.include_all_within()) {
        _neighbours[i] << j;
      }
    }
  }

  if (_verbose) {
    cerr << "Ligand protein distances " << acc.n() << " btw " << acc.minval() << " and " << acc.maxval() << " ave " << acc.average() << '\n';
    Accumulator_Int<uint32_t> acc_nbrs;
    int zero_neighbour_count = 0;
    for (uint32_t i = 0; i < _atoms_in_ligand; ++i) {
      if (_neighbours[i].empty()) {
        ++zero_neighbour_count;
      } else {
        acc_nbrs.extra(_neighbours[i].size());
      }
    }
    cerr << zero_neighbour_count << " ligand atoms have no nearby protein atoms, dist " << _config.include_all_within() << '\n';
    cerr << "Ligand atoms have between " << acc_nbrs.minval() << " and " << acc_nbrs.maxval() << " ave " << acc_nbrs.average() << " ave neighbours in protein\n";
  }

  return 1;
}

int
Options::InitialiseRingsInProtein() {
  _rings_in_protein = _protein->nrings();
  _protein_rings = new RingAndCenter[_rings_in_protein];

  for (int i = 0; i < _rings_in_protein; ++i) {
    _protein_rings[i].Initialise(*_protein, _protein->ringi(i));
  }

  return 1;
}

int
Options::InitialiseRingsInLigand() {
  _rings_in_ligand = _ligand->nrings();
  _ligand_rings = new RingAndCenter[_rings_in_ligand];

  for (int i = 0; i < _rings_in_ligand; ++i) {
    _ligand_rings[i].Initialise(*_ligand, _ligand->ringi(i));
    if (_ligand->ringi(i)->is_aromatic()) {
      _ligand_rings[i].aromatic = true;
    }
  }

  return 1;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  set_store_pdb_atom_information(1);

  if (! cl.option_present('C')) {
    cerr << "Options::Initialise:must specify config proto via the -C option\n";
    return 0;
  }

  if (! cl.option_present('S')) {
    cerr << "Options::Initialise:must specify file name stem for output file(s) (-S)\n";
    return 0;
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "Output files created with name step '" << _stem << "'\n";
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    std::optional<LigandProtein::InteractionConfig> maybe_proto =
           iwmisc::ReadTextProto<LigandProtein::InteractionConfig>(fname);
    if (! maybe_proto) {
      cerr << "Options::Initialise:cannot read proto config '" << fname << "'\n";
      return 0;
    }
    _config = *maybe_proto;
    cerr << _config.ShortDebugString() << '\n';
  }

  if (! cl.option_present('P')) {
    cerr << "Must specify protein structure via the -P option\n";
    return 0;
  }

  set_ignore_self_bonds(1);

  if (cl.option_present('P')) {
    IWString fname = cl.string_value('P');
    if (! ReadProtein(fname.null_terminated_chars())) {
      cerr << "Options::Initialise:cannot read protein '" << fname << "'\n";
      return 0;
    }

    // We could wait to reconnect Hydrogens, only reconnecting those that are
    // later identified as being of interest. But this is fast.
    // int h_single_bonds = ReconnectHydrogens(*_protein);
    int h_single_bonds = MakeBondsFromDistances(*_protein);
    if (_verbose) {
      //cerr << "Reconnected " << h_single_bonds << " disconnected H atoms in protein\n";
      cerr << "Created " << h_single_bonds << " bonds protein\n";
    }
  }

  _input_type = FILE_TYPE_SDF;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SDF;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }


  if (! ReadLigand(cl[0])) {
    cerr << "Cannot read ligand '" << cl[0] << "'\n";
    return 0;
  }

  if (cl.option_present('I')) {
    IWString fname = cl.string_value('I');
    WriteLigand(*_ligand, fname);
  }

  return 1;
}

int
Options::WriteLigand(Molecule& m,
                     IWString fname) {
  Molecule_Output_Object output;
  output.add_output_type(FILE_TYPE_SDF);
  output.add_output_type(FILE_TYPE_SMI);
  if (! output.new_stem(fname)) {
    cerr << "Options::WriteLigand:cannot open '" << fname << "'\n";
    return 0;
  }

  Molecule mcopy(m);
  const int matoms = m.natoms();
  for (int i = 1; i < matoms; ++i) {
    mcopy.set_isotope(i, i);
  }

  return output.write(mcopy);
}

int
Options::ReadProtein(const char * fname) {
  data_source_and_type<Molecule> input(FILE_TYPE_PDB, fname);
  if (! input.good()) {
    cerr << "Options::ReadProtein:cannot open '" << fname << "'\n";
    return 0;
  }

  _protein = input.next_molecule();
  if (_protein == nullptr) {
    cerr << "Options::ReadProtein:cannot read protein from '" << fname << "'\n";
    return 0;
  }

  _atoms_in_protein = _protein->natoms();
  _protein_scratch = new int[_atoms_in_protein];

  if (_verbose) {
    cerr << "Protein contains " << _atoms_in_protein << " atoms\n";
  }

  _id_connected_atoms.Initialise(_protein);

  return 1;
}

int
Options::ReadLigand(const char * fname) {
  data_source_and_type<Molecule> input(_input_type, fname);
  if (! input.good()) {
    cerr << "Options::ReadLigand:cannot open '" << fname << "'\n";
    return 0;
  }

  _ligand = input.next_molecule();
  if (_ligand == nullptr) {
    cerr << "Options::ReadLigand:cannot read ligand from '" << fname << "'\n";
    return 0;
  }

  _atoms_in_ligand = _ligand->natoms();

  if (_verbose) {
    cerr << "Ligand contains " << _atoms_in_ligand << " atoms\n";
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Found " << setw(3) << _acceptor_donor_found << " acceptor in ligand -> donor    in protein\n";
  output << "Found " << setw(3) << _donor_acceptor_found << " donor    in ligand -> acceptor in protein\n";
  output << "Found " << setw(3) << _donor_to_pi_found << " donor    in ligand -> pi       in protein\n";
  output << "Found " << setw(3) << _pipi_found << " pi-pi interactions\n";
  output << "Wrote " << _file_name_index << " files\n";
  return 1;
}

int
Options::Process(IWString_and_File_Descriptor& output) {

  _ligand_ring_system = new_int(_atoms_in_ligand);
  if (_ligand->nrings()) {
    _ligand->compute_aromaticity_if_needed();
    int rs = IdentifyRingSystems(*_ligand, _ligand_ring_system);
    if (_verbose) {
      cerr << "Ligand contains " << rs << " ring systems\n";
    }
  }

  _protein_ring_system = new_int(_atoms_in_protein);
  int rs = IdentifyRingSystems(*_protein, _protein_ring_system);
  if (_verbose) {
    cerr << "Protein contains " << rs << " ring systems\n";
  }

  InitialiseDistanceMatrix();

  InitialiseRingsInLigand();
  InitialiseRingsInProtein();

  IdentifyHBonds();

  return 1;
}

// Just copied from numbered_smiles. See if this tool works or not...

int
NumberbyRingSystem(Molecule& m, int* ring_system) {
  const int nrings = m.nrings();
  if (nrings == 0) {
    return 1;
  }

  std::unique_ptr<int[]> ring_already_done(new_int(nrings + 1));

  int next_to_assign = 1;
  for (int i = 0; i < nrings; ++i, ++next_to_assign) {
    if (ring_already_done[i]) {
      continue;
    }
    const Ring* ri = m.ringi(i);
    ri->set_vector(ring_system, next_to_assign);
    for (int j = i + 1; j < nrings; ++j) {
      if (ring_already_done[j]) {
        continue;
      }
      const Ring* rj = m.ringi(j);
      if (rj->fused_system_identifier() != ri->fused_system_identifier()) {
        continue;
      }
      rj->set_vector(ring_system, next_to_assign);
      ring_already_done[j] = 1;
    }
  }

  // Add doubly bonded exocyclic bonds
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (ring_system[i] == 0) {
      continue;
    }
    const Atom* a = m.atomi(i);
    // This next condition fails for the protein where we do not have bonding info...
    if (a->ncon() == 2 || a->nbonds() == a->ncon()) {
      continue;
    }
    for (const Bond* b : *a) {
      if (! b->is_double_bond()) {
        continue;
      }
      const atom_number_t j = b->other(i);
      if (ring_system[j] == ring_system[i] || m.atomic_number(j) == 6) {
        continue;
      }
      // We have an exocyclic double bond.
      // If we want to restrict to just =O or =S, uncomment this.
      // if (m.ncon(j) != 1) {
      //   continue;
      // } 
      ring_system[j] = ring_system[i];
    }
  }

  return next_to_assign;
}

int
Options::IdentifyRingSystems(Molecule& m, int* ring_system) {
  int number_ring_systems = NumberbyRingSystem(m, ring_system);

  return number_ring_systems;
}

int
Options::OkPiPiDistance(const Space_Vector<float>& p1,
                        const Space_Vector<float>& p2) const {
  const float dist = p1.distance(p2);
  if (dist < _config.pipi_spec().min_distance()) {
    return 0;
  }
  if (dist > _config.pipi_spec().max_distance()) {
    return 0;
  }

  return 1;
}

int
Options::OkPiPiAngle(const Molecule& m1,
                     const RingAndCenter& rc1,
                     const Molecule& m2,
                     const RingAndCenter& rc2) const {
  const Atom* atom = m1.atomi(rc1.atoms[0]);
  Space_Vector<float> v10(*atom);
  atom = m1.atomi(rc1.atoms[1]);
  Space_Vector<float> v11(*atom);
  atom = m1.atomi(rc1.atoms[2]);
  Space_Vector<float> v12(*atom);

  v10 -= v11;
  v12.cross_product(v10);

  atom = m2.atomi(rc2.atoms[0]);
  Space_Vector<float> v20(*atom);
  atom = m2.atomi(rc2.atoms[1]);
  Space_Vector<float> v21(*atom);
  atom = m2.atomi(rc2.atoms[2]);
  Space_Vector<float> v22(*atom);

  v20 -= v21;
  v22.cross_product(v20);

  const float angle = v12.angle_between(v22);
  if (angle < _config.pipi_spec().min_angle()) {
    return 0;
  }
  if (angle > _config.pipi_spec().max_angle()) {
    return 0;
  }

  return 1;
}

int
Options::IdentifyHBonds() {
  for (uint32_t i = 0; i < _atoms_in_ligand; ++i) {
    const Atom* ligand_atom = _ligand->atomi(i);
    const atomic_number_t z = ligand_atom->atomic_number();
    // Skip carbon.
    if (z == 6) {
      continue;
    }
    // Skip halogens.
    if (z == 9 || z == 17 || z == 35 || z == 53) {
      continue;
    }

    if (z == 1) {
      IdentifyAcceptorsInProtein(i);
    } else {
      IdentifyDonorsInProtein(i);
    }
  }

  IdentifyPiPi();

  return 1;
}

IWString
Options::MakePiPiName(const RingAndCenter& rci,
                      const RingAndCenter& rcj) const {
  IWString result;
  result << "_PIPI_";
  for (uint32_t i = 0; i < rci.atoms.size(); ++i) {
    if (i > 0) {
      result << '-';
    }
    result << rci.atoms[i];
  }
  result << ResidueName(rcj.atoms[0]);

  return result;
}

// Return true if `ring_and_center` is likely aromatic.
// Update ring_and_center.aromatic with the result we
// determine.
int
IsAromatic(Molecule& m,
           RingAndCenter& ring_and_center) {
  if (ring_and_center.aromatic >= 0) {
    return ring_and_center.aromatic;
  }

  uint32_t sp2_count = 0;
  for (atom_number_t a : ring_and_center.atoms) {
    if (m.GeometryIsSp2(a)) {
      ++sp2_count;
    }
  }

  if (sp2_count == ring_and_center.atoms.size()) {
    ring_and_center.aromatic = 1;
    return 1;
  }

  if (sp2_count + 1 == ring_and_center.atoms.size()) {
    ring_and_center.aromatic = 1;
    return 1;
  }

  ring_and_center.aromatic = 0;
  return 0;
}

int
Options::IdentifyPiPi() {
  for (int i = 0; i < _rings_in_ligand; ++i) {
    const RingAndCenter& rci = _ligand_rings[i];
    if (! rci.aromatic) {
      continue;
    }
    for (int j = 0; j < _rings_in_protein; ++j) {
      RingAndCenter& rcj = _protein_rings[j];
      if (! OkPiPiDistance(rci.center, rcj.center)) {
        continue;
      }
      if (! OkPiPiAngle(*_ligand, rci, *_protein, rcj)) {
        continue;
      }
      if (!IsAromatic(*_protein, rcj)) {
        continue;
      }

      const IWString name = MakePiPiName(rci, rcj);
      if (_verbose) {
        cerr << name << '\n';
      }
      Write(*_ligand, rci.atoms, _ligand_ring_system,
            *_protein, rcj.atoms, _protein_ring_system, name);
      ++_pipi_found;
    }
  }

  return _pipi_found;
}

Set_of_Atoms
AsSet(atom_number_t a) {
  Set_of_Atoms result;
  result << a;
  return result;
}

// Make a name for H bond donor `h` interacting with a pi ring `rcj`.
IWString
Options::MakeHPiName(atom_number_t h,
                     const RingAndCenter& rcj) const {
  IWString result;
  result << "HPi_" << h << '_';
  result <<  ResidueName(rcj.atoms.first());

  return result;
}

IWString
Options::MakeHAccName(atom_number_t h,
                atom_number_t acc) const {
  IWString result;
  result << "HA_" << h << "_" << acc << '_';
  result << ResidueName(acc);

  return result;
}

int
Options::IdentifyAcceptorsInProtein(atom_number_t h) {
  Set_of_Atoms h_as_set;
  h_as_set << h;

  const Atom* ligand_atom = _ligand->atomi(h);
  const atom_number_t anchor = ligand_atom->other(h, 0);
  const Atom* anchor_atom = _ligand->atomi(anchor);

  for (uint32_t nbr : _neighbours[h]) {
    const Atom* protein_atom = _protein->atomi(nbr);
    if (! OkHBondDistance(DistanceBetweenLPAtoms(h, nbr))) {
      continue;
    }
    const atomic_number_t z = protein_atom->atomic_number();
    // the two most common cases first.
    if (z == 1 || z == 6) {
      continue;
    }

    // We want O or N atoms in the protein.
    if (z == 7) {
    } else if (z == 8) {
    } else {
      continue;
    }

    if (! OkHBondAngle(BondAngle(anchor_atom, ligand_atom, protein_atom) * RAD2DEG)) {
      continue;
    }

    const IWString name = MakeHAccName(h, nbr);
    if (_verbose) {
      cerr << name << '\n';
    }
    Write(*_ligand, h_as_set, _ligand_ring_system,
          *_protein, AsSet(nbr), _protein_ring_system, name);
    ++_donor_acceptor_found;
  }

  // Now look for `h` interacting with pi systems.
  // cerr << "Checking " << h << " against " << _rings_in_protein << " rings in the protein\n";
  for (int i = 0; i < _rings_in_protein; ++i) {
    if (! OkHPiBond(*_ligand, h, _protein_rings[i])) {
      continue;
    }
    if (! OkHPiAngle(*_ligand, h, *_protein, _protein_rings[i])) {
      continue;
    }
    if (! IsAromatic(*_protein, _protein_rings[i])) {
      continue;
    }

    const IWString name = MakeHPiName(h, _protein_rings[i]);
    if (_verbose) {
      cerr << name << '\n';
    }
    Write(*_ligand, h_as_set, _ligand_ring_system,
          *_protein, _protein_rings[i].atoms, _protein_ring_system, name);
    ++_donor_to_pi_found;
  }

  return 1;
}

int
Options::OkHPiBond(const Molecule& m,
                   atom_number_t hydrogen,
                   const RingAndCenter& rc) {
  
  const float distance = m.atomi(hydrogen)->distance(rc.center);
  if (distance < _config.hpi_spec().min_distance()) {
    return 0;
  }
  if (distance > _config.hpi_spec().max_distance()) {
    return 0;
  }

  return 1;
}

// Atom `zatom` in molecule `m1` is thinking about making a HBond-Pi
// interaction with `rc` which is in molecule `m2`. Check the angle.
int
Options::OkHPiAngle(const Molecule& m1,
                   atom_number_t zatom,
                   const Molecule& m2,
                   const RingAndCenter& rc) const {
  Space_Vector<float> v1(*m1.atomi(zatom));
  Space_Vector<float> v2(*m2.atomi(rc.atoms.front()));

  v1 -= rc.center;
  v2 -= rc.center;
  const float angle = v1.angle_between(v2) * RAD2DEG;
  if (angle < _config.hpi_spec().min_angle()) {
    return 0;
  }
  if (angle > _config.hpi_spec().max_angle()) {
    return 0;
  }
  return 1;
}

IWString
Options::ResidueName(atom_number_t zatom) const {
  const resizable_array_p<PDB_Stored_Atom_Information>& pdb_info= stored_pdb_atom_information_last_molecule_read();
  const PDB_Stored_Atom_Information* info = pdb_info[zatom];
  IWString result = info->residue_name();
  result << '.' << info->residue_number();
  return result;
}

int
OkDonorAtomType(const Atom* a) {
  const atomic_number_t z = a->atomic_number();
  if (z == 6) {
    return 0;
  }

  if (z == 7 || z == 8) {
    return 1;
  }

  return 0;
}

// Make a name for a found interaction between `acc` in the ligand and
// `don` in the protein.
IWString
Options::MakeAccDonName(atom_number_t acc,
                        atom_number_t don) const {
  IWString result;
  result << "AccDon_" << acc << '-' << don << '_';
  result << ResidueName(don);

  return result;
}

// `zatom` is a non-hydrogen atom in the ligand.
int
Options::IdentifyDonorsInProtein(atom_number_t zatom) {
  Set_of_Atoms zatom_as_set;
  zatom_as_set << zatom;

  const Atom* ligand_atom = _ligand->atomi(zatom);
  if (ligand_atom->atomic_number() == 7) {
  } else if (ligand_atom->atomic_number() == 8) {
  } else {
    return 0;
  }

  // Some atoms can both donate and accept, change if needed...
  if (_ligand->hcount(zatom)) {
    return 0;
  }

  for (uint32_t nbr : _neighbours[zatom]) {
    const Atom* protein_atom = _protein->atomi(nbr);
    if (protein_atom->atomic_number() != 1) {
      continue;
    }
    if (! OkHBondDistance(DistanceBetweenLPAtoms(zatom, nbr))) {
      continue;
    }
    if (protein_atom->ncon() == 0) {
      cerr << "Options::IdentifyDonorsInProtein:disconnected H in protein, ignoring\n";
      continue;
    }
    const atom_number_t anchor = protein_atom->other(nbr, 0);
    if (! OkDonorAtomType(_protein->atomi(anchor))) {
      continue;
    }
    if (! OkHBondAngle(BondAngle(ligand_atom, protein_atom, _protein->atomi(anchor)) * RAD2DEG)) {
      continue;
    }

    const IWString name = MakeAccDonName(zatom, nbr);
    if (_verbose) {
      cerr << name << '\n';
    }
    Write(*_ligand, zatom_as_set, _ligand_ring_system,
          *_protein, AsSet(nbr), _protein_ring_system, name);
    ++_acceptor_donor_found;
    if (_verbose) {
      const IWString residue_name = ResidueName(nbr);
      cerr << "Acceptor " << zatom << " " << _ligand->smarts_equivalent_for_atom(zatom) << " in ligand " << nbr << " in protein " << residue_name << ", d " << DistanceBetweenLPAtoms(zatom, nbr) << '\n';
    }
  }

  return 1;
}

// `zatom` is an N or O. At this stage, we do not know if it is
// an acceptor or donor. Figure out the attached hydrogen status
// and branch accordingly.
// Might want some more complex logic about what can, and cannot, be
// a donor/acceptor.
int
Options::IdentifyHBonds(atom_number_t zatom) {
  const Atom* a = _ligand->atomi(zatom);
  Set_of_Atoms hydrogens_attached;

  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    const Atom* aj = _ligand->atomi(j);
    if (aj->atomic_number() != 1) {
      continue;
    }
    hydrogens_attached << j;
  }

  if (hydrogens_attached.empty()) {
    return IdentifyHBondsAcceptorInLigand(zatom);
  } else {
    return IdentifyHBondsDonorInLigand(zatom, hydrogens_attached);
  }
}

int
Options::OkHBondDistance(distance_t distance) const {
  if (distance < _config.hbond_spec().min_distance()) {
    return 0;
  }
  if (distance > _config.hbond_spec().max_distance()) {
    return 0;
  }

  return 1;
}

int
Options::OkHBondAngle(angle_t angle) const {
#ifdef DEBUG_OKHBONDANGLE
  cerr << "OkHBondAngle:checking " << angle << " vs " << _config.hbond_spec().ShortDebugString() << '\n';
#endif

  if (angle < _config.hbond_spec().min_angle()) {
    return 0;
  }
  if (angle > _config.hbond_spec().max_angle()) {
    return 0;
  }

  return 1;
}

int
Options::IdentifyHBondsAcceptorInLigand(atom_number_t zatom) {
  const Atom* ligand_atom = _ligand->atomi(zatom);

  for (auto nbr : _neighbours[zatom]) {
    const Atom* nbr_atom = _protein->atomi(nbr);
    if (nbr_atom->atomic_number() != 1) {
      continue;
    }
    if (! OkHBondDistance(DistanceBetweenAtoms(ligand_atom, nbr_atom))) {
      continue;
    }
    const atom_number_t j = nbr_atom->other(nbr, 0);
    const angle_t angle = BondAngle(ligand_atom, nbr_atom, _protein->atomi(j)) * RAD2DEG;
    if (! OkHBondAngle(angle)) {
      continue;
    }
  }

  return 1;
}

int
Options::IdentifyHBondsDonorInLigand(atom_number_t zatom,
                                     const Set_of_Atoms& hydrogens_attached) {
  return 1;
}

void
AddAtom(const Atom* from, Molecule& m) {
  Atom* a = new Atom(from->element());
  a->setxyz(from->x(), from->y(), from->z());

  m.add(a);
}

// `destination` is a partially built subset of `from`.
// For each atom already in `destination` look for other
// atoms that share the same `ring_system` value and add
// those atoms to `destination`, updating xref as atoms are added.
int
AddRingSystem(const Molecule& from,
              const int* ring_system,
              int* xref,
              Molecule& destination) {
  const int matoms = from.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    // Skip if not in the subset.
    if (xref[i] < 0) {
      continue;
    }
    // Skip if not part of a ring.
    if (ring_system[i] <= 0) {
      continue;
    }
    // Add all atoms in the same ring system
    for (int j = 0; j < matoms; ++j) {
      // Skip if already part of the subset.
      if (xref[i] >= 0) {
        continue;
      }
      // Skip unless same ring system.
      if (ring_system[i] != ring_system[j]) {
        continue;
      }
      xref[i] = destination.natoms();
      AddAtom(from.atomi(j), destination);
      ++rc;
    }
  }

  return rc;
}

// Atoms have been added to `destination` from `from`, with
// `xref` specifying which atoms have been added. Add any bonds
// implied.
int
AddBonds(const Molecule& from,
         const int* xref,
         Molecule& destination) {
  const int nedges = from.nedges();
  int rc = 0;
  for (int i = 0; i < nedges; ++i) {
    const Bond* b = from.bondi(i);
    const atom_number_t a1 = b->a1();
    if (xref[a1] < 0) {
      continue;
    }
    const atom_number_t a2 = b->a2();
    if (xref[a2] < 0) {
      continue;
    }
    if (destination.are_bonded(xref[a1], xref[a2])) {
      continue;
    }
    destination.add_bond(xref[a1], xref[a2], b->btype());
    ++rc;
  }

  return rc;
}

// We have found an interaction involving the atoms `s1` in `m1` and
//the atoms `s2` in `m2`.
// Variables `ring_system1,2` describe ring system membership in each
// molecule, and are used to extend the subset written to nearby atoms and
// ring systems.
int
Options::Write(Molecule& m1,
               const Set_of_Atoms& s1,
               const int* ring_system1,
               const Molecule& m2,
               const Set_of_Atoms& s2,
               const int* ring_system2,
               const IWString& name) {
  const int m1natoms = m1.natoms();
  const int m2natoms = m2.natoms();
  std::unique_ptr<int[]> xref1(new_int(m1natoms, -1));
  std::unique_ptr<int[]> xref2(new_int(m2natoms, -1));

  Molecule to_write;
  to_write.set_name(name);

  for (atom_number_t a : s1) {
    if (xref1[a] >= 0) {
      continue;
    }
    xref1[a] = to_write.natoms();
    AddAtom(m1.atomi(a), to_write);
    if (_config.bond_extent() > 0) {
      for (int i = 0; i < m1natoms; ++i) {
        if (xref1[i] >= 0) {
          continue;
        }
        if (m1.bonds_between(a, i) <= _config.bond_extent()) {
          xref1[i] = to_write.natoms();
          AddAtom(m1.atomi(i), to_write);
        }
      }
    }
  }

  for (atom_number_t a : s2) {
    if (xref2[a] >= 0) {
      continue;
    }
    xref2[a] = to_write.natoms();
    AddAtom(m2.atomi(a), to_write);
    if (_config.bond_extent() > 0) {
      AddProximalAtoms(m2, xref2.get(), a, to_write);
    }
  }

  AddRingSystem(m1, ring_system1, xref1.get(), to_write);
  AddRingSystem(m2, ring_system2, xref2.get(), to_write);

  AddBonds(m1, xref1.get(), to_write);
  AddBonds(m2, xref2.get(), to_write);

  //to_write.make_implicit_hydrogens_explicit();

  IWString fname;
  fname << _stem << _file_name_index << ".sdf";
  ++_file_name_index;
  Molecule_Output_Object output;
  output.add_output_type(FILE_TYPE_SDF);
  if (! output.new_stem(fname)) {
    cerr << "Options::Write:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Writing " << fname << " with " << to_write.natoms() << " atoms\n";
  }

  return output.write(to_write);
}

// Adding proximal atoms in the protein cannot use the distance matrix, so
// we have an iterative expansion.
int
Options::AddProximalAtoms(const Molecule& m,
                          int* xref,
                          atom_number_t zatom,
                          Molecule& destination) {
  static constexpr int kMagicNumber = std::numeric_limits<int>::max();

  const Atom* a = m.atomi(zatom);

  Set_of_Atoms frontier;

  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (xref[j] >= 0) {
      continue;
    }
    frontier << j;
    xref[j] = kMagicNumber;
  }

  Set_of_Atoms next_frontier;
  for (int r = 0; r < _config.bond_extent() && ! frontier.empty(); ++r) {
    // cerr << "r = " << r << " frontier " << frontier << '\n';
    for (int frontier_atom : frontier) {
      const Atom* a = m.atomi(frontier_atom);
      for (const Bond* b : *a) {
        int j = b->other(frontier_atom);
        // Either already in the set, or being added here.
        if (xref[j] >= 0) {
          continue;
        }
        next_frontier << j;
        xref[j] = kMagicNumber;
      }
    }
    frontier = next_frontier;
    next_frontier.resize_keep_storage(0);
  }

  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (xref[i] != kMagicNumber) {
      continue;
    }
    xref[i] = destination.natoms();
    AddAtom(m.atomi(i), destination);
    ++rc;
  }

  AddBonds(m, xref, destination);

  return rc;
}

#ifdef NO_LONGER_USED
int
Options::WriteHBond(atom_number_t ligand_atom, atom_number_t protein_atom) {
  std::unique_ptr<int[]> xref_ligand(new_int(_atoms_in_ligand, -1));
  std::unique_ptr<int[]> xref_protein(new_int(_atoms_in_protein, -1));

  Molecule to_write;

  for (uint32_t i = 0; i < _atoms_in_ligand; ++i) {
    const float d = DistanceBetweenLPAtoms(i, protein_atom);
    if (d > _config.include_all_within()) {
      continue;
    }
    xref_ligand[i] = to_write.natoms();
    AddAtom(_ligand->atomi(i), to_write);
  }

  for (uint32_t i = 0; i < _atoms_in_protein; ++i) {
    const float d = DistanceBetweenLPAtoms(ligand_atom, i);
    if (d > _config.include_all_within()) {
      continue;
    }
    xref_protein[i] = to_write.natoms();
    AddAtom(_protein->atomi(i), to_write);
  }

  AddRingSystem(*_ligand, _ligand_ring_system, xref_ligand.get(), to_write);
  AddRingSystem(*_protein, _protein_ring_system, xref_protein.get(), to_write);

  AddBonds(*_ligand, xref_ligand.get(), to_write);
  AddBonds(*_protein, xref_protein.get(), to_write);

  to_write.make_implicit_hydrogens_explicit();

  return 1;
}
#endif

int
Options::ReconnectHydrogen(Molecule& m,
                  atom_number_t zatom) const {
  float cutoff;
  if (m.atomic_number(zatom) == 1) {
    cutoff = 1.3;
  } else {
    cutoff = _config.reconnect_within();
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == zatom) {
      continue;
    }
    const float d = m.distance_between_atoms(i, zatom);
    if (d > cutoff) {
      continue;
    }
    m.add_bond(i, zatom, SINGLE_BOND);
    return 1;
  }

  return 0;
}


int
Options::ReconnectHydrogens(Molecule& m) const {
  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Atom* a = m.atomi(i);
    if (a->ncon() > 0) {
      continue;
    }
    rc += ReconnectHydrogen(m, i);
  }
  return rc;
}

int
ImplicitHydrogens(const Atom* a) {
  int ih;
  if (! const_cast<Atom*>(a)->compute_implicit_hydrogens(ih)) {
    return 0;
  }

  return ih;
}

int
Options::MakeBondsFromDistances(Molecule& m) const {
  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Atom* a = m.atomi(i);
    int implicit_hydrogens = ImplicitHydrogens(a);
    if (implicit_hydrogens == 0) {
      continue;
    }
    float cutoff;
    if (a->atomic_number() == 1) {
      cutoff = 1.4;
    } else {
      cutoff = _config.reconnect_within();
    }
    for (int j = i + 1; j < matoms; ++j) {
      if (a->is_bonded_to(j)) {
        continue;
      }
      if (DistanceBetweenAtoms(a, m.atomi(j)) > cutoff) {
        continue;
      }
      m.add_bond(i, j, SINGLE_BOND);
      ++rc;
      // If atom I cannot take any more bonds, we are done with it.
      --implicit_hydrogens;
      if (implicit_hydrogens == 0) {
        break;
      }
    }
  }

  return rc;
}

int
LigandProteinInteractions(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return 1;
}

int
LigandProteinInteractions(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:P:C:S:I:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  if (! options.Process(output)) {
    cerr << "Fatal error\n";
    return 1;
  }

  options.Report(cerr);

  return 0;
}
}  // namespace ligand_protein

int
main(int argc, char ** argv) {

  int rc = ligand_protein::LigandProteinInteractions(argc, argv);

  return rc;
}
