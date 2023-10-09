// Implementation of xlogp
// J. Chem. Inf. Comput. Sci. 1997, 37, 3, 615–621
// Wang Fu Lai

#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/xlogp.h"

namespace xlogp {

using std::cerr;

// Values from Table 1 of the paper.

constexpr int kCH3NoPi = 1;
constexpr int kCH3Pi = 2;
constexpr int kCh3X = 3;
constexpr int kCh2R2NoPi = 4;
constexpr int kCh2R2Pi = 5;
constexpr int kCh2RX = 6;
constexpr int kCh2X2 = 7;
constexpr int kChR3NoPi = 8;
constexpr int kChR3Pi = 9;
constexpr int kChR2X = 10;
constexpr int kChRX2 = 11;

constexpr int kCR4NoPi = 12;
constexpr int kCR4Pi = 13;
constexpr int kCR3X = 14;
constexpr int kCR2X2 = 15;
constexpr int kCRX3 = 16;
constexpr int kCX4 = 17;

constexpr int kCarbonSp2RCH2 = 18;
constexpr int kCarbonSp2RdblCHR = 19;
constexpr int kCarbonSp2RdblCHX = 20;
constexpr int kCarbonSp2XdblCHR = 21;
constexpr int kCarbonSp2RdblCR2 = 22;
constexpr int kCarbonSp2RdbkCRX = 23;
constexpr int kCarbonSp2XdblCR2 = 24;
constexpr int kCarbonSp2XdblCX2 = 25;

constexpr int kCaromRCHR = 26;
constexpr int kCaromRCHX = 27;
constexpr int kCaromXCHX = 28;
constexpr int kCaromRCRR = 29;
constexpr int kCaromRCXR = 30;
constexpr int kCaromRCRX = 31;
constexpr int kCaromRCXX = 32;
constexpr int kCaromXCAX = 33;
constexpr int kCaromACAA = 34;

constexpr int kCarbonSpRCH = 35;
constexpr int kCarbonSpRCR = 36;

constexpr int kHydrogen = 37;

constexpr int kOxygenROHnopi = 38;
constexpr int kOxygenROHpi = 39;
constexpr int kOxygenXOH = 40;
constexpr int kOxygenROR = 41;
constexpr int kOxygenROX = 42;

constexpr int kOxygenAromatic = 43;
constexpr int kOxygenOdblR = 44;
constexpr int kOxygenOdblX = 45;

constexpr int kNitrogenSp3NH2Rnopi = 46;
constexpr int kNitrogenSp3NH2Rpi = 47;
constexpr int kNitrogenSp3X = 48;
constexpr int kNitrogenNS2RNHR = 49;
constexpr int kNitrogenNSp3RNH = 50;
constexpr int kNitrogenNSp3NR3 = 51;
constexpr int kNitrogenNSp3NR2X = 52;

constexpr int kNitrogenSp2RdblNH = 53;
constexpr int kNitrogenSp2RdblNX = 54;
constexpr int kNitrogenSp2XdblNR = 55;
constexpr int kNitrogenSp2XdblNX = 56;

constexpr int kNitrogenAromatic = 57;
constexpr int kNitrogenTrigonalPlanarRNHR = 58;
constexpr int kNitrogenTrigonalPlanarRNHX = 59;
constexpr int kNitrogenTrigonalPlanarAromatic = 60;
constexpr int kNitrogenTrigonalPlanarNA3nonring = 61;
constexpr int kNitrogenTrigonalPlanarNA3ring = 62;

constexpr int kAmideNitrogenNH2 = 63;
constexpr int kAmideNitrogenNH = 64;
constexpr int kAmideNitrogen = 65;

constexpr int kSulphurSH = 66;
constexpr int kSulphurD2 = 67;
constexpr int kSulphurArom = 68;
constexpr int kSulphurSdbleR = 69;
constexpr int kSulphurSulfoxide = 70;
constexpr int kSulphurSulfone = 71;

constexpr int kFluorine = 72;
constexpr int kChlorine = 73;
constexpr int kBromine = 74;
constexpr int kIodine = 75;
constexpr int kPhosphorus = 76;
constexpr int kCyano = 77;
constexpr int kNCS = 78;
constexpr int kNO = 79;
constexpr int kNitro = 80;
// Added IAW
constexpr int kQuatNitrogen = 81;
constexpr int kPositiveAromaticNitrogen = 82;
constexpr int kPositiveNitrogenJoinedN = 83;

constexpr int kAmidine = 84;

// Always 1 larger than the highest defined index.
constexpr int kMaxArrayIndex = 85;

double type_to_score[] = {
   0.0,    // 0
   0.484,  // 1
   0.168,  // 2
  -0.181,  // 3
   0.358,  // 4
   0.009,  // 5
  -0.344,  // 6
  -0.439,  // 7
   0.051,  // 8
  -0.138,  // 9
  -0.417,  // 10
  -0.454,  // 11
  -0.378,  // 12
   0.223,  // 13
  -0.598,  // 14
  -0.396,  // 15
  -0.699,  // 16
  -0.362,  // 17

   0.395,  // 18
   0.236,  // 19
  -0.166,  // 20
   1.726,  // 21
   0.098,  // 22
  -0.108,  // 23
   1.637,  // 24
   1.774,  // 25

   0.281,  // 26
   0.142,  // 27
   0.715,  // 28
   0.302,  // 29
  -0.064,  // 30
   0.079,  // 31
   0.200,  // 32
   0.869,  // 33
   0.316,  // 34

   0.054,  // 35
   0.347,  // 36

   0.046,  // 37

  -0.399,  // 38
  -0.029,  // 39
  -0.330,  // 40
   0.397,  // 41
   0.068,  // 42
   0.327,  // 43

  -2.057,  // 44
   0.218,  // 45

  -0.582,  // 46
  -0.449,  // 47
  -0.774,  // 48
   0.040,  // 49
  -0.381,  // 50
   0.443,  // 51
  -0.117,  // 52

  -2.052,  // 53
  -1.716,  // 54
   0.321,  // 55
  -0.921,  // 56

  -0.704,  // 57

   0.119,  // 58
   1.192,  // 59
   0.434,  // 60
   0.587,  // 61
   0.668,  // 62

  -0.791,  // 63
  -0.212,  // 64
   0.016,  // 65

   0.752,  // 66
   1.071,  // 67
   0.964,  // 68

  -1.817,  // 69
  -1.214,  // 70
  -0.778,  // 71
   0.493,  // 72
   1.010,  // 73
   1.187,  // 74
   1.489,  // 75
  -0.802,  // 76

  -0.256,  // 77
   1.626,  // 78
   0.077,  // 79
   0.264,  // 80

   // Added IAW
   -3.00,  // 81
   -4.50,  // 82
   -1.70,  // 83

   // Amidine
   -1.00   // 84
};

constexpr int kFailed = -1;

bool display_unclassified_atom_messages = 1;

int
ProcessNewFragmentParameter(const XLogP::XlogpParameter& proto) {
  if (! proto.has_index() || ! proto.has_value()) {
    cerr << "ProcessNewFragmentParameter:must have both index and value\n";
    return 0;
  }

  if (proto.index() == 0 || proto.index() > kMaxArrayIndex) {
    cerr << "ProcessNewFragmentParameter::invalid index\n";
    return 0;
  }

  type_to_score[proto.index()] = proto.value();

  return 1;
}

int
ReadNewFragmentParameters(const XLogP::XlogpParameters& proto) {
  for (const auto& contribution : proto.contribution()) {
    if (! ProcessNewFragmentParameter(contribution)) {
      cerr << "ReadNewFragmentParameters:cannot process " << contribution.ShortDebugString() << '\n';
      return 0;
    }
  }

  return 1;
}

int
ReadNewFragmentParameters(IWString& fname) {
  std::optional<XLogP::XlogpParameters> maybe_proto =
        iwmisc::ReadTextProto<XLogP::XlogpParameters>(fname);
  if (! maybe_proto) {
    cerr << "xlogp::ReadNewFragmentParameters:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadNewFragmentParameters(*maybe_proto);
}

// By default we apply the corrections, but for testing it may be helpful
// to be able to turn them off.
int apply_corrections = 1;

// for debugging it is helpful to get a detailed dump of what got assigned
// to each atom.
int display_assignments = 0;

void
SetDisplayAtomAssignments(int s) {
  display_assignments = s;
}

void
ForTestingSetApplyCorrections(int s) {
  apply_corrections = s;
}

struct PerMoleculeData {
  public:
    const int natoms;

    atomic_number_t* atomic_number;
    int* ncon;
    int* hcount;
    int* unsaturation;
    int* aromatic;
    int* pi_electrons;
    int* attached_carbons;
    int* ring_bond_count;

    int halogen_count;
    int positively_charged_atoms;

  // private:
    int DearomatizePyrroles(Molecule& m);

  public:
    PerMoleculeData(Molecule& m);
    ~PerMoleculeData();

    // return true if any of the neighbours of `zatom` have pi electrons.
    bool NbrsHavePiElectrions(const Molecule& m, atom_number_t zatom) const;

    // the number of neighbours that have pi electons.
    int NbrsWithPiElectrons(const Molecule& m,
                        atom_number_t zatom) const;
};

bool
PerMoleculeData::NbrsHavePiElectrions(const Molecule& m,
                        atom_number_t zatom) const {
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (pi_electrons[o]) {
      return true;
    }
  }

  return false;
}

int
PerMoleculeData::NbrsWithPiElectrons(const Molecule& m,
                        atom_number_t zatom) const {
  int rc = 0;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (pi_electrons[o]) {
      ++rc;
    }
  }

  return rc;
}

PerMoleculeData::PerMoleculeData(Molecule& m) : natoms(m.natoms()) {
  atomic_number = new atomic_number_t[natoms];
  ncon = new int[natoms];
  aromatic = new int[natoms];
  pi_electrons = new int[natoms];
  hcount = new int[natoms];
  ring_bond_count = new int[natoms];
  halogen_count = 0;
  positively_charged_atoms = 0;
  for (int i = 0; i < natoms; ++i) {
    const Atom& a = m[i];

    atomic_number[i] = a.atomic_number();
    ncon[i] = a.ncon();
    aromatic[i] = m.is_aromatic(i);
    hcount[i] = m.hcount(i);
    (void) m.pi_electrons(i, pi_electrons[i]);
    ring_bond_count[i] = m.ring_bond_count(i);
    if (m.formal_charge(i) > 0) {
      ++positively_charged_atoms;
    }
  }

  for (int i = 0; i < natoms; ++i) {
    if (atomic_number[i] == 6 ||
        atomic_number[i] == 7 ||
        atomic_number[i] == 8) {
      continue;
    }

    // Maybe we will need to separate these sometime.
    if (atomic_number[i] == 9) {
      ++halogen_count;
    } else if (atomic_number[i] == 17 ||
               atomic_number[i] == 35 ||
               atomic_number[i] == 53) {
      ++halogen_count;
    }
  }

  unsaturation = new_int(natoms);
  attached_carbons = new_int(natoms);
  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (atomic_number[a1] == 6) {
      ++attached_carbons[a2];
    }
    if (atomic_number[a2] == 6) {
      ++attached_carbons[a1];
    }

    if (b->is_single_bond()) {
      continue;
    }

    if (b->is_double_bond()) {
      ++unsaturation[a1];
      ++unsaturation[a2];
    } else {
      unsaturation[a1] += 2;
      unsaturation[a2] += 2;
    }
  }
}

PerMoleculeData::~PerMoleculeData() {
  delete [] atomic_number;
  delete [] ncon;
  delete [] pi_electrons;
  delete [] unsaturation;
  delete [] aromatic;
  delete [] hcount;
  delete [] attached_carbons;
  delete [] ring_bond_count;
}

const int hydrophobic_carbon[] {
  0,  // 0
  1,  // 1
  1,  // 2
  0,  // 3
  1,  // 4
  1,  // 5
  0,  // 6
  0,  // 7
  1,  // 8
  1,  // 9
  0,  // 10
  0,  // 11
  1,  // 12
  1,  // 13
  0,  // 14
  0,  // 15
  0,  // 16
  0,  // 17
  1,  // 18
  1,  // 19
  0,  // 20
  0,  // 21
  1,  // 22
  0,  // 23
  0,  // 24
  0,  // 25
  0,  // 26
  0,  // 27
  0,  // 28
  0,  // 29
  0,  // 30
  0,  // 31
  0,  // 32
  0,  // 33
  0,  // 34
  0,  // 35
  0,  // 36
  0,  // 37
  0,  // 38
  0,  // 39
  0,  // 40
  0,  // 41
  0,  // 42
  0,  // 43
  0,  // 44
  0,  // 45
  0,  // 46
  0,  // 47
  0,  // 48
  0,  // 49
  0,  // 50
  0,  // 51
  0,  // 52
  0,  // 53
  0,  // 54
  0,  // 55
  0,  // 56
  0,  // 57
  0,  // 58
  0,  // 59
  0,  // 60
  0,  // 61
  0,  // 62
  0,  // 63
  0,  // 64
  0,  // 65
  0,  // 66
  0,  // 67
  0,  // 68
  0,  // 69
  0,  // 70
  0,  // 71
  0,  // 72
  0,  // 73
  0,  // 74
  0,  // 75
  0,  // 76
  0,  // 77
  0,  // 78
  0,  // 79
  0   // 80
};

bool
ContainsHeteroatoms(const Molecule& m,
                    const PerMoleculeData& per_molecule_data) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (per_molecule_data.atomic_number[i] != 6) {
      return true;
    }
  }

  return false;
}


double
HydropphobicCarbonCorrection(Molecule& m,
                const PerMoleculeData& per_molecule_data,
                const int* status) {
  const int matoms = m.natoms();
  if (ContainsHeteroatoms(m, per_molecule_data)) {
    // cerr << "ContainsHeteroatoms, no HydropphobicCarbonCorrection\n";
    return 0.0;
  }

  int hydrophobic_atoms = 0;
  for (int i = 0; i < matoms; ++i) {
    if (hydrophobic_carbon[status[i]]) {
      hydrophobic_atoms += 1;
    }
  }

  return hydrophobic_atoms * 0.19;
}

// `oxygen` and `carbon` are part of a carboxyllic acid.
// `oxygen` is [OD1H]
// Is there an amine alpha to the carbon.
bool
AmineNearby(Molecule& m, atom_number_t oxygen, atom_number_t carbon,
            const PerMoleculeData& per_molecule_data,
            atom_number_t& nitrogen) {
  for (const Bond* b : m[carbon]) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t c = b->other(carbon);
    if (c == oxygen) {
      continue;
    }

    if (per_molecule_data.atomic_number[c] != 6) {
      return false;
    }
    if (per_molecule_data.unsaturation[c] > 0) {
      return false;
    }
    if (per_molecule_data.hcount[c] == 0) {
      return false;
    }

    for (const Bond* b2 : m[c]) {
      atom_number_t n = b2->other(c);
      if (per_molecule_data.atomic_number[n] != 7) {
        continue;
      }
      if (per_molecule_data.NbrsWithPiElectrons(m, n)) {
        continue;
      }
      if (per_molecule_data.hcount[n] > 0) {
        nitrogen = n;
        return true;
      }
    }
  }

  return false;
}

bool
IsCarboxyllicAcid(Molecule& m,
                  atom_number_t o,
                  const PerMoleculeData& per_molecule_data,
                  atom_number_t& carbon) {
  carbon = m.other(o, 0);
  if (per_molecule_data.atomic_number[carbon] != 6) {
    return false;
  }
  if (per_molecule_data.ncon[carbon] != 3) {
    return false;
  }
  if (per_molecule_data.unsaturation[carbon] != 1) {
    return false;
  }

  for (const Bond * b : m[carbon]) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o2 = b->other(carbon);
    if (per_molecule_data.atomic_number[o2] == 8) {
      return true;
    }
  }

  return false;
}

double
AminoAcidCorrection(Molecule& m,
                    const PerMoleculeData& per_molecule_data,
                    const int* status) {
  const int matoms = m.natoms();
  int namino_acids = 0;

  std::unique_ptr<int[]> aa(new_int(matoms));
  for (int i = 0; i < matoms; ++i) {
    if (aa[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 8 ||
        per_molecule_data.ncon[i] != 1 ||
        per_molecule_data.unsaturation[i]) {
      continue;
    }
    atom_number_t carbon;
    if (! IsCarboxyllicAcid(m, i, per_molecule_data, carbon)) {
      continue;
    }
    // cerr << "Ot acid, check nitrogen amine\n";
    atom_number_t nitrogen;
    if (! AmineNearby(m, i, carbon, per_molecule_data, nitrogen)) {
      continue;
    }
    if (aa[nitrogen]) {
      continue;
    }

    aa[i] = 1;
    aa[carbon] = 1;
    aa[nitrogen] = 1;
    ++namino_acids;
  }

  // cerr << "Found " << namino_acids << " AMINO ACIDS\n";
  return namino_acids * -2.27;
}

enum class GemHalogenStatus {
  kNot = 0,
  kHasFluorine,
  kNoFluorine
};

// `h1` is a halogen. Examine the atom to which it is connected
// and determine whether there are any gem halogen interactions on
// that atom.
GemHalogenStatus
GetGemHalogenStatus(const Molecule& m,
                    atom_number_t h1,
                    const PerMoleculeData& per_molecule_data,
                    int* processed,
                    int& fluorine_count,
                    int& other_halogen_count) {
  processed[h1] = 1;

  atom_number_t c = m.other(h1, 0);
  if (per_molecule_data.ncon[c] < 2) {
    return GemHalogenStatus::kNot;
  }

  fluorine_count = 0;
  other_halogen_count = 0;
  for (const Bond* b : m[c]) {
    atom_number_t o = b->other(c);

    if (per_molecule_data.atomic_number[o] == 9) {
      ++fluorine_count;
      processed[o] = 1;
     } else if (per_molecule_data.atomic_number[o] == 17 ||
                per_molecule_data.atomic_number[o] == 35 ||
                per_molecule_data.atomic_number[o] == 53) {
      ++other_halogen_count;
      processed[o] = 1;
    }
  }

  // cerr << "FOund fluorine_count " << fluorine_count << " other_halogen_count " << other_halogen_count << '\n';
  if (fluorine_count + other_halogen_count < 2) {
    return GemHalogenStatus::kNot;
  }

  if (fluorine_count) {
    return GemHalogenStatus::kHasFluorine;
  } else {
    return GemHalogenStatus::kNoFluorine;
  }
}

double
HalogenCorrection(Molecule& m,
            PerMoleculeData& per_molecule_data,
            const int* status) {
  // cerr << "Molecule contains " << per_molecule_data.halogen_count << " halogen atoms\n";
  if (per_molecule_data.halogen_count < 2) {
    return 0.0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<int[]> processed(new_int(matoms));

  double rc = 0.0;

  for (int i = 0; i < matoms; ++i) {
    if (processed[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] == 6 ||
        per_molecule_data.atomic_number[i] == 7 ||
        per_molecule_data.atomic_number[i] == 8) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] == 9 ||
        per_molecule_data.atomic_number[i] == 17 ||
        per_molecule_data.atomic_number[i] == 35 ||
        per_molecule_data.atomic_number[i] == 53) {
      int fluorine = 0;
      int other_halogen = 0;
      auto gem_status = GetGemHalogenStatus(m, i, per_molecule_data,
                        processed.get(), fluorine, other_halogen);
      if (gem_status == GemHalogenStatus::kHasFluorine) {
        int total = fluorine + other_halogen;
        if (total > 3) {
          total = 3;
        }
        total = 1;
        rc += 0.08 * total;
      } else if (gem_status == GemHalogenStatus::kNoFluorine) {
        int total = fluorine + other_halogen;
        if (total > 3) {
          total = 3;
        }
        total = 1;
        rc += total * -0.26;
      } else {
        // cerr << "not halogen type\n";
      }
    }
  }

  // cerr << "Apply halogen correction " << rc << '\n';
  return rc;
}

// `a1` and `a2` are a possible intramolecular hydrogen bonding pair.
// find ring bonds encountered in moving from a1 to a2 and if there
// are more than 1, return 1.
int
MoreThanOneRingBond(Molecule& m, atom_number_t a1, atom_number_t a2,
                    int& ring_bonds_found) {
  // We need to move one bond closer to `a2`.
  const int d = m.bonds_between(a1, a2) - 1;
  for (const Bond * b : m[a1]) {
    atom_number_t o = b->other(a1);
    if (o == a2) {
      return 0;
    }

    if (m.bonds_between(o, a2) != d) {
      continue;
    }

    if (b->nrings()) {
      ++ring_bonds_found;
      if (ring_bonds_found > 1) {
        return 1;
      }
    }

    if (MoreThanOneRingBond(m, o, a2, ring_bonds_found)) {
      return 1;
    }
  }

  return 0;
}

double
IntraMolecularHBondCorrection(Molecule& m,
            PerMoleculeData& per_molecule_data,
            const int* status) {
  const int matoms = m.natoms();
  int rc = 0;

  m.recompute_distance_matrix();
  m.ring_membership();

  std::unique_ptr<int[]> already_used(new_int(matoms));

  // First identify donors, then a corresponding acceptor.
  for (int i = 0; i < matoms; ++i) {
    if (per_molecule_data.hcount[i] == 0) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] == 6) {
      continue;
    }
    // We have a heteroatoms with at least 1 Hydrogen atoms.
    for (int j = 0; j < matoms; ++j) {
      if (per_molecule_data.atomic_number[j] == 6) {
        continue;
      }
      if (m.bonds_between(i, j) != 4) {
        continue;
      }
      // This is not in the paper, but seems to be how it is working.
      if (per_molecule_data.atomic_number[i] == per_molecule_data.atomic_number[j]) {
        continue;
      }
      if (per_molecule_data.atomic_number[j] == 9) {
      } else if (per_molecule_data.atomic_number[j] == 8 &&
                 per_molecule_data.unsaturation[j] == 1) {
      } else if (per_molecule_data.atomic_number[j] != 6 &&
                 per_molecule_data.hcount[j]) {
      } else {
        continue;
      }

      Set_of_Atoms between;
      m.atoms_between(i, j, between);
      if (between.any_members_set_in_array(already_used.get(), 1)) {
        continue;
      }
      int ring_bonds_found = 0;
      if (MoreThanOneRingBond(m, i, j, ring_bonds_found)) {
        continue;
      }

      between.set_vector(already_used.get(), 1);

      ++rc;
      // We only allow each atom to be in 1 internal intramolecular hbond.
      // Not great, but it would be hard to fix this...
      break;
    }
  }

  // cerr << "Detect " << rc << " intramolecular hbonds\n";
  return rc * 0.60;
}

// Note that this deliberately matches guanidine forms.
// `n1` is the =N in an amidine or ganidine.
bool
IsAmidine(const Molecule& m,
          atom_number_t n1,
          const PerMoleculeData& per_molecule_data) {
  const atom_number_t carbon = m.other(n1, 0);

  if (per_molecule_data.atomic_number[carbon] != 6) {
    return false;
  }
  if (per_molecule_data.ncon[carbon] != 3) {
    return false;
  }
  if (per_molecule_data.aromatic[carbon]) {
    return false;
  }

  atom_number_t singly_bonded_nitrogen = INVALID_ATOM_NUMBER;
  for (const Bond* b : m[carbon]) {
    if (b->is_double_bond()) {
      continue;
    }
    atom_number_t j = b->other(carbon);
    if (per_molecule_data.atomic_number[j] == 6) {
      continue;
    } else if (per_molecule_data.atomic_number[j] != 7){
      return false;
    } else if (per_molecule_data.ncon[j] == 1) {
      singly_bonded_nitrogen = j;
    }
  }

  return singly_bonded_nitrogen >= 0;
}

// this will also hit guanidines
double
AmidineCorrection(Molecule& m,
            PerMoleculeData& per_molecule_data,
            const int* status) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (per_molecule_data.atomic_number[i] != 7 ||
        per_molecule_data.ncon[i] != 1 ||
        per_molecule_data.unsaturation[i] != 1) {
      continue;
    }

    if (IsAmidine(m, i, per_molecule_data)) {
      ++rc;
    }
  }

  return rc * type_to_score[kAmidine];  // 84
}

double
Corrections(Molecule& m,
            PerMoleculeData& per_molecule_data,
            const int* status) {
  if (! apply_corrections) {
    return 0.0;
  }
  // cerr << "Corrections being applied\n";

  double rc = HydropphobicCarbonCorrection(m, per_molecule_data, status);
  rc += AminoAcidCorrection(m, per_molecule_data, status);
  rc += IntraMolecularHBondCorrection(m, per_molecule_data, status);
  rc += HalogenCorrection(m, per_molecule_data, status);
  rc += AmidineCorrection(m, per_molecule_data, status);

  return rc;
}

double
IdentifyNO(Molecule& m,
            PerMoleculeData& per_molecule_data,
            int* status) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int o1 = 0; o1 < matoms; ++o1) {
    if (status[o1]) {
      continue;
    }

    if (per_molecule_data.atomic_number[o1] != 8 ||
        per_molecule_data.ncon[o1] != 1 ||
        per_molecule_data.unsaturation[o1] != 1) {
      continue;
    }

    const Bond* b = m[o1][0];
    if (! b->is_double_bond()) {
      continue;
    }

    const atom_number_t n = b->other(o1);
    if (per_molecule_data.atomic_number[n] != 7 ||
        per_molecule_data.ncon[n] != 2 ||
        per_molecule_data.unsaturation[n] != 1) {
      continue;
    }

    status[o1] = kNO;  // 79
    status[n] = kNO;  // 79
    ++rc;
  }

  return rc * type_to_score[kNO];
}

double
IdentifyNO2(Molecule& m,
            PerMoleculeData& per_molecule_data,
            int* status) {
  int rc = 0;
  const int matoms = m.natoms();
  for (int o1 = 0; o1 < matoms; ++o1) {
    if (status[o1]) {
      continue;
    }

    if (per_molecule_data.atomic_number[o1] != 8 ||
        per_molecule_data.ncon[o1] != 1 ||
        per_molecule_data.unsaturation[o1] != 1) {
      continue;
    }

    const Bond* b = m[o1][0];
    if (! b->is_double_bond()) {
      continue;
    }

    const atom_number_t n = b->other(o1);
    if (per_molecule_data.atomic_number[n] != 7 ||
        per_molecule_data.ncon[n] != 3 ||
        per_molecule_data.unsaturation[n] != 2) {
      continue;
    }

    for (const Bond*b : m[n]) {
      if (! b->is_double_bond()) {
        continue;
      }
      atom_number_t o2 = b->other(n);
      if (o2 == o1) {
        continue;
      }
      if (per_molecule_data.atomic_number[o2] != 8 ||
          per_molecule_data.ncon[o2] != 1) {
        continue;
      }

      status[o1] = kNitro;  // 80
      status[n] = kNitro;  // 80
      status[o2] = kNitro;  // 80
      ++rc;
    }
  }

  return rc * type_to_score[kNitro];
}


double
IdentifyNCS(Molecule& m,
            PerMoleculeData& per_molecule_data,
            int* status) {
  int rc = 0;

  const int matoms = m.natoms();
  for (int s = 0; s < matoms; ++s) {
    if (status[s]) {
      continue;
    }

    if (per_molecule_data.atomic_number[s] != 16) {
      continue;
    }
    if (per_molecule_data.ncon[s] != 1) {
      continue;
    }

    const Bond* b = m[s][0];
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t c = b->other(s);
    if (per_molecule_data.atomic_number[c] != 6) {
      continue;
    }
    if (per_molecule_data.ncon[c] != 2) {
      continue;
    }

    for (const Bond* b : m[c]) {
      if (! b->is_double_bond()) {
        continue;
      }
      atom_number_t n = b->other(c);
      if (n == s) {
        continue;
      }
      if (per_molecule_data.atomic_number[n] != 7) {
        continue;
      }

      status[s] = kNCS;
      status[c] = kNCS;
      status[n] = kNCS;
      ++rc;
      break;
    }
  }

  return rc * type_to_score[kNCS];
}

double
IdentifyCyano(Molecule& m,
              PerMoleculeData& per_molecule_data,
              int* status) {
  int rc = 0;
  for (const Bond * b : m.bond_list()) {
    if(! b->is_triple_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    atomic_number_t z1 = m.atomic_number(b->a1());
    atomic_number_t z2 = m.atomic_number(b->a2());
    if (status[a1] > 0 || status[a2] > 0) {
      continue;
    }
    if (z1 == 6 && z2 == 7) {
    } else if (z1 == 7 && z2 == 6) {
    } else {
      continue;
    }
    ++rc;
    status[a1] = kCyano;
    status[a2] = kCyano;
  }

  return rc * type_to_score[kCyano];
}

bool
IsAmide(const Molecule& m,
        atom_number_t n,
        const PerMoleculeData& per_molecule_data) {
  for (const Bond* b : m[n]) {
    atom_number_t cs = b->other(n);  // Carbon or Sulphur
    if (per_molecule_data.unsaturation[cs] == 0) {
      continue;
    }

    if (per_molecule_data.atomic_number[cs] == 6 &&
        per_molecule_data.ncon[cs] == 3) {
    } else if (per_molecule_data.atomic_number[cs] == 16 &&
        per_molecule_data.ncon[cs] > 2) {
    } else {
      continue;
    }

    for (const Bond* b2 : m[cs]) {
      if (! b2->is_double_bond()) {
        continue;
      }
      atom_number_t o = b2->other(cs);
      if (per_molecule_data.atomic_number[o] == 8) {
        return true;
      }
      if (per_molecule_data.atomic_number[o] == 16) {
        return true;
      }
    }
  }

  return false;
}

double
IdentifyAmideNitrogens(const Molecule& m,
                       const PerMoleculeData& per_molecule_data,
                       int* status) {
  double rc = 0.0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 7 ||
        per_molecule_data.aromatic[i]) {
      continue;
    }
    if (per_molecule_data.unsaturation[i]) {
      continue;
    }

    if (! IsAmide(m, i, per_molecule_data)) {
      continue;
    }

    if (per_molecule_data.hcount[i] == 2) {
      status[i] = kAmideNitrogenNH2;
      rc += type_to_score[kAmideNitrogenNH2];  // 63
    } else if (per_molecule_data.hcount[i] == 1) {
      status[i] = kAmideNitrogenNH;
      rc += type_to_score[kAmideNitrogenNH];  // 64
    } else {
      status[i] = kAmideNitrogen;
      rc += type_to_score[kAmideNitrogen];  // 65
    }
  }

  return rc;
}

std::optional<double>
ClassifyNSP2(Molecule& m,
             atom_number_t zatom,
             const PerMoleculeData& per_molecule_data,
             int * status) {
  assert(m.atomic_number(zatom) == 7);

  if (per_molecule_data.ncon[zatom] == 1 ||
      per_molecule_data.attached_carbons[zatom] == 2) {
    status[zatom] = kNitrogenSp2RdblNH;
    return type_to_score[kNitrogenSp2RdblNH];  // 53
  }

  bool heteratom_at_end_of_double_bond = false;
  for (const Bond* b : m[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    heteratom_at_end_of_double_bond = (per_molecule_data.atomic_number[o] != 6);
    break;
  }

  if (! heteratom_at_end_of_double_bond) {
    status[zatom] = kNitrogenSp2RdblNX;
    return type_to_score[kNitrogenSp2RdblNX];  // 54
  }

  if (per_molecule_data.attached_carbons[zatom] == 1) {
    status[zatom] = kNitrogenSp2XdblNR;
    return type_to_score[kNitrogenSp2XdblNR];  // 55
  }

  status[zatom] = kNitrogenSp2XdblNX;
  return type_to_score[kNitrogenSp2XdblNX];  // 56
}

double
IdentifyNitrogenSP2(Molecule& m,
             const PerMoleculeData& per_molecule_data,
             int * status) {
  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 7){
      continue;
    }
    if (per_molecule_data.unsaturation[i] != 1) {
      continue;
    }
    if (per_molecule_data.aromatic[i]) {
      continue;
    }

    std::optional<double> x = ClassifyNSP2(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      status[i] = kFailed;
    }
  }

  return rc;
}

std::optional<double>
ClassifyNSP3(const Molecule& m,
             atom_number_t zatom,
             const PerMoleculeData& per_molecule_data,
             int * status) {
  // cerr << per_molecule_data.hcount[zatom] << " hcount, " << per_molecule_data.attached_carbons[zatom] << " carbon\n";
  if (per_molecule_data.hcount[zatom] == 2) {
    if (per_molecule_data.attached_carbons[zatom] == 0) {
      status[zatom] = kNitrogenSp3X;
      return type_to_score[kNitrogenSp3X];   // 48
    }

    if (per_molecule_data.NbrsHavePiElectrions(m, zatom)) {
      status[zatom] = kNitrogenSp3NH2Rpi;
      return type_to_score[kNitrogenSp3NH2Rpi];  // 47
    } else {
      status[zatom] = kNitrogenSp3NH2Rnopi;
      return type_to_score[kNitrogenSp3NH2Rnopi];  // 46
    }
  }

  if (per_molecule_data.attached_carbons[zatom] == 3) {
    status[zatom] = kNitrogenNSp3NR3;
    return type_to_score[kNitrogenNSp3NR3];   // 51
  }

  if (per_molecule_data.hcount[zatom] == 1) {
    if (per_molecule_data.attached_carbons[zatom] == 2) {
      status[zatom] = kNitrogenNS2RNHR;
      return type_to_score[kNitrogenNS2RNHR];  // 49
    }

    status[zatom] = kNitrogenNSp3RNH;
    return type_to_score[kNitrogenNSp3RNH];  // 50
  }

  if (per_molecule_data.ncon[zatom] == 3) {
    status[zatom] = kNitrogenNSp3NR2X;
    return type_to_score[kNitrogenNSp3NR2X];  // 51
  }

  if (display_unclassified_atom_messages) {
    cerr << "xlogp::ClassifyNSP3:Unrecognised atom " << m.name() << '\n';
  }

  return std::nullopt;
}

double
IdentifyNitrogenSP3(const Molecule& m,
            const PerMoleculeData& per_molecule_data,
            int * status) {
  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 7) {
      continue;
    }
    if (per_molecule_data.aromatic[i]) {
      continue;
    }
    if (per_molecule_data.unsaturation[i]) {
      continue;
    }

    std::optional<double> x = ClassifyNSP3(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      status[i] = kFailed;
    }
  }

  return rc;
}

std::optional<double>
ClassifyTrigonalPlanarNitrogen(Molecule& m,
            atom_number_t zatom,
            const PerMoleculeData& per_molecule_data,
            int * status) {
  if (per_molecule_data.hcount[zatom] == 1 &&
      per_molecule_data.ncon[zatom] == 2 &&
      per_molecule_data.attached_carbons[zatom] == 2) {
    status[zatom] = kNitrogenTrigonalPlanarRNHR;
    return type_to_score[kNitrogenTrigonalPlanarRNHR];  // 58
  }

  if (per_molecule_data.ncon[zatom] == 2 &&
      per_molecule_data.hcount[zatom] == 1) {
    status[zatom] = kNitrogenTrigonalPlanarRNHX;
    return type_to_score[kNitrogenTrigonalPlanarRNHX];  // 59
  }

  // The paper says only 5 membered rings, but things do not match with
  // what they report unless we include 6 membered rings - but not larger.
  if (per_molecule_data.ncon[zatom] == 3 &&
      per_molecule_data.ring_bond_count[zatom] &&
      (m.in_ring_of_given_size(zatom, 5) ||
       m.in_ring_of_given_size(zatom, 6)) ) {
    status[zatom] = kNitrogenTrigonalPlanarNA3ring;
    return type_to_score[kNitrogenTrigonalPlanarNA3ring];  // 62
  }

  if (per_molecule_data.ncon[zatom] == 3 &&
      per_molecule_data.ring_bond_count[zatom] == 0) {
    status[zatom] = kNitrogenTrigonalPlanarNA3nonring;
    return type_to_score[kNitrogenTrigonalPlanarNA3nonring];  // 61
  }

  // This does not make much sense. But the paper says that the ring version
  // is for 5 membered rings. We see above that we need to also include 6
  // membered rings to replicate what they have. So I guess larger rings
  // But the computed value would be closer to what they report if we used
  // the ring version...
  // are classified as non ring???
  if (per_molecule_data.ncon[zatom] == 3 &&
      per_molecule_data.ring_bond_count[zatom]) {
    status[zatom] = kNitrogenTrigonalPlanarNA3nonring;
    return type_to_score[kNitrogenTrigonalPlanarNA3nonring];  // 61
  }

  if (m.formal_charge(zatom) == 1) {
    status[zatom] = kPositiveNitrogenJoinedN;  // 83
    return type_to_score[kPositiveNitrogenJoinedN];
  }

  if (display_unclassified_atom_messages) {
    cerr << "xlogp::ClassifyTrigonalPlanarNitrogen:unrecognised form '" << m.name() << "'\n";
  }

  return std::nullopt;
}

double
IdentifyTrigonalPlanarNitrogen(Molecule& m,
            const PerMoleculeData& per_molecule_data,
            int* status) {
  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 7) {
      continue;
    }
    if (per_molecule_data.aromatic[i]) {
      continue;
    }
    if (per_molecule_data.unsaturation[i]) {
      continue;
    }
    if (per_molecule_data.NbrsWithPiElectrons(m, i) < 2) {
      continue;
    }

    std::optional<double> x = ClassifyTrigonalPlanarNitrogen(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      status[i] = kFailed;
    }
  }

  return rc;
}

// Beware different requirements if aromaticity definitions are changed
double
IdentifyTrigonalPlanarNitrogenAromatic(Molecule& m,
                const PerMoleculeData& per_molecule_data,
                int* status) {
  double rc = 0.0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 7) {
      continue;
    }
    if (per_molecule_data.ring_bond_count[i] != 2) {
      continue;
    }
    if (per_molecule_data.hcount[i] == 0) {
      continue;
    }
#ifdef DAYLIGHT_AROMATICITY
    if (! per_molecule_data.aromatic[i]) {
      continue;
    }
#endif
    if (! m.in_ring_of_given_size(i, 5)) {
      continue;
    }
    status[i] = kNitrogenTrigonalPlanarAromatic;
    rc += type_to_score[kNitrogenTrigonalPlanarAromatic];  // 60
  }

  return rc;
}

double
IdentifyNitrogenAromatic(Molecule& m,
                const PerMoleculeData& per_molecule_data,
                int* status) {
  double rc = 0.0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 7) {
      continue;
    }
    if (! per_molecule_data.aromatic[i]) {
      continue;
    }
    status[i] = kNitrogenAromatic;
    rc += type_to_score[kNitrogenAromatic];  // 57
  }

  return rc;
}

double
IdentifyChargedNitrogen(const Molecule& m,
                        const PerMoleculeData& per_molecule_data,
                        int * status) {
  if (per_molecule_data.positively_charged_atoms == 0) {
    return 0;
  }

  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    const atomic_number_t z = per_molecule_data.atomic_number[i];
    if (z != 7) {
      continue;
    }
    if (m.formal_charge(i) != 1) {
      continue;
    }

    if (m.saturated(i)) {
      if (per_molecule_data.NbrsHavePiElectrions(m, i)) {
      } else {
      }
      status[i] = kQuatNitrogen;  // 81
      rc += type_to_score[kQuatNitrogen];
    } else if (per_molecule_data.ncon[i] == per_molecule_data.attached_carbons[i]) {
      status[i] = kPositiveAromaticNitrogen;  // 82
      rc += type_to_score[kPositiveAromaticNitrogen];
    } else {
      status[i] = kPositiveNitrogenJoinedN;  // 83
      rc += type_to_score[kPositiveNitrogenJoinedN];
    }
  }

  return rc;
}

double
IdentifyHalogens(const Molecule& m,
            const PerMoleculeData& per_molecule_data,
            int * status) {
  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    const atomic_number_t z = per_molecule_data.atomic_number[i];
    if (z == 6 || z == 7 || z == 8) {
      continue;
    }

    if (z == 9) {
      status[i] = kFluorine;
      rc += type_to_score[kFluorine];
      continue;
    }
    if (z == 17) {
      status[i] = kChlorine;
      rc += type_to_score[kChlorine];
      continue;
    }
    if (z == 35) {
      status[i] = kBromine;
      rc += type_to_score[kBromine];
      continue;
    }
    if (z == 53) {
      status[i] = kIodine;
      rc += type_to_score[kIodine];
      continue;
    }
  }

  return rc;
}

std::optional<double>
ClassifySulphur(Molecule& m,
                atom_number_t zatom,
                const PerMoleculeData& per_molecule_data,
                int * status) {
  if (per_molecule_data.ncon[zatom] == 1 &&
      per_molecule_data.hcount[zatom] == 1) {
    status[zatom] = kSulphurSH;
    return type_to_score[kSulphurSH];  // 66
  }

  if (per_molecule_data.ncon[zatom] == 1 &&
      per_molecule_data.unsaturation[zatom] == 1) {
    status[zatom] = kSulphurSdbleR;
    return type_to_score[kSulphurSdbleR];  // 69
  }

  if (per_molecule_data.aromatic[zatom]) {
    status[zatom] = kSulphurArom;
    return type_to_score[kSulphurArom];  // 68
  }

  if (per_molecule_data.ncon[zatom] == 2 &&
      per_molecule_data.NbrsWithPiElectrons(m, zatom) == 2 &&
      m.in_ring_of_given_size(zatom, 5)) {
    status[zatom] = kSulphurArom;
    return type_to_score[kSulphurArom];  // 68
  }

  if (per_molecule_data.ncon[zatom] == 2) {
    status[zatom] = kSulphurD2;
    return type_to_score[kSulphurD2];  // 67
  }

  int doubly_bonded_oxygens = 0;
  for (const Bond* b : m[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(zatom);
    if (per_molecule_data.atomic_number[o] == 8) {
      ++doubly_bonded_oxygens;
    }
  }

  if (doubly_bonded_oxygens == 1) {
    status[zatom] = kSulphurSulfoxide;
    return type_to_score[kSulphurSulfoxide];  // 70
  }

  if (doubly_bonded_oxygens == 2) {
    status[zatom] = kSulphurSulfone;
    return type_to_score[kSulphurSulfone];  // 71
  }

  if (display_unclassified_atom_messages) {
    cerr << "xlogp::ClassifySulphur:unrecognised form " << m.name() << '\n';
  }

  return std::nullopt;
}

double
IdentifySulphur(Molecule& m,
            const PerMoleculeData& per_molecule_data,
            int * status) {
  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 16) {
      continue;
    }

    std::optional<double> x = ClassifySulphur(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      if (display_unclassified_atom_messages) {
        cerr << "xlogp::IdentifySulphur:unrecognised Sulphur form '" << m.name() << "'\n";
      }
      status[i] = kFailed;
    }
  }

  return rc;
}

std::optional<double>
ClassifyPhosphorus(const Molecule& m,
                atom_number_t zatom,
                const PerMoleculeData& per_molecule_data,
                int * status) {
  for (const Bond* b : m[zatom]) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t o = b->other(zatom);
    if (per_molecule_data.atomic_number[o] == 8) {
      status[zatom] = kPhosphorus;
      return  type_to_score[kPhosphorus];
    }
  }

  return std::nullopt;
}

double
IdentifyPhosphorus(const Molecule& m,
            const PerMoleculeData& per_molecule_data,
            int * status) {
  const int matoms = m.natoms();
  double rc = 0.0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 15) {
      continue;
    }

    std::optional<double> x = ClassifyPhosphorus(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      status[i] = kFailed;
    }
  }

  return rc;
}

std::optional<double>
ClassifyAromaticCarbon(Molecule& m,
                       atom_number_t zatom,
                       PerMoleculeData& per_molecule_data,
                       int * status) {

  if (per_molecule_data.attached_carbons[zatom] == 2 &&
     per_molecule_data.hcount[zatom] == 1) {
    status[zatom] = kCaromRCHR;
    return type_to_score[kCaromRCHR];  // 26
  } 
  
  if (per_molecule_data.attached_carbons[zatom] == 1 &&
      per_molecule_data.hcount[zatom] == 1) {
    status[zatom] = kCaromRCHX;
    return type_to_score[kCaromRCHX];  // 27
  } 
  
  if (per_molecule_data.attached_carbons[zatom] == 0 &&
             per_molecule_data.hcount[zatom] == 1) {
    status[zatom] = kCaromXCHX;
    return type_to_score[kCaromXCHX];  // 28
  } 
  
  if (per_molecule_data.ncon[zatom] == 3 &&
      per_molecule_data.ring_bond_count[zatom] == 2 &&
      per_molecule_data.attached_carbons[zatom] == 3) {
    status[zatom] = kCaromRCRR;
    return type_to_score[kCaromRCRR];  // 29
  }

  if (m.ncon(zatom) != 3) {
    cerr << "xlogp::ClassifyAromaticCarbon: invalid aromatic carbon " << m.name() << '\n';
    return std::nullopt;
  }

  m.compute_aromaticity_if_needed();

  int aromatic_bonds = 0;
  int heteroatoms_in_ring = 0;
  bool heteroatom_outside_ring = false;

  for (const Bond * b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (b->is_aromatic()) {
      ++aromatic_bonds;
    }

    if (b->nrings() && b->is_aromatic()) {
      if (per_molecule_data.atomic_number[o] != 6) {
        ++heteroatoms_in_ring;
      }
    } else {
      if (per_molecule_data.atomic_number[o] != 6) {
        heteroatom_outside_ring = true;
      }
    }
  }

  if (aromatic_bonds == 3) {
    status[zatom] = kCaromACAA;
    return type_to_score[kCaromACAA];   // 34
  }

  if (heteroatoms_in_ring == 0 && heteroatom_outside_ring) {
    status[zatom] = kCaromRCXR;  // 30
    return type_to_score[kCaromRCXR];
  }

  if ( heteroatoms_in_ring == 1 && ! heteroatom_outside_ring) {
    status[zatom] = kCaromRCRX;
    return type_to_score[kCaromRCRX];   // 31
  }

  if ( heteroatoms_in_ring == 1 && heteroatom_outside_ring) {
    status[zatom] = kCaromRCXX;
    return type_to_score[kCaromRCXX];   // 32
  }

  if (heteroatoms_in_ring == 2) {
    status[zatom] = kCaromXCAX;
    return type_to_score[kCaromXCAX];  // 33
  }

  status[zatom] = kCaromRCRR;
  return type_to_score[kCaromRCRR];
}

std::optional<double>
ClassifyAromaticNitrogen(Molecule& m,
                        atom_number_t zatom,
                        PerMoleculeData& per_molecule_data,
                        int * status) {
  double rc = 0;
  return rc;
}

double
ClassifySinglyConnectedOxygen(const Molecule& m,
                        atom_number_t zatom,
                        const PerMoleculeData& per_molecule_data,
                        int* status) {
  assert(m.atomic_number(zatom) == 8 && m.ncon(zatom) == 1);
  if (per_molecule_data.unsaturation[zatom] == 1) {
    if (per_molecule_data.attached_carbons[zatom] == 1) {
      status[zatom] = kOxygenOdblR;
      return type_to_score[kOxygenOdblR];  // 44
    } else {
      status[zatom] = kOxygenOdblX;
      return type_to_score[kOxygenOdblX];  // 45
    }
  }

  if (per_molecule_data.attached_carbons[zatom] == 1) {
    if (per_molecule_data.NbrsHavePiElectrions(m, zatom)) {
      status[zatom] = kOxygenROHpi;
      return type_to_score[kOxygenROHpi];  // 39
    } else {
      status[zatom] = kOxygenROHnopi;
      return type_to_score[kOxygenROHnopi];  // 38
    }
  } 

  status[zatom] = kOxygenXOH;
  return type_to_score[kOxygenXOH];  // 40
}

double
ClassifyDoublyConnectedOxygen(const Molecule& m,
                        atom_number_t zatom,
                        const PerMoleculeData& per_molecule_data,
                        int* status) {
  assert(m.atomic_number(zatom) == 8 && m.ncon(zatom) == 2);

  // Should we also check for 5 membered ring only?
  if (per_molecule_data.ring_bond_count[zatom] == 2 &&
      per_molecule_data.NbrsWithPiElectrons(m, zatom) == 2) {
    status[zatom] = kOxygenAromatic;
    return type_to_score[kOxygenAromatic];  // 43
  }

  if (per_molecule_data.attached_carbons[zatom] == 2) {
    status[zatom] = kOxygenROR;
    return type_to_score[kOxygenROR];  // 41
  }

  status[zatom] = kOxygenROX;
  return 0.068; // 42
}

double
IdentifyOxygen(Molecule& m,
               PerMoleculeData& per_molecule_data,
               int * status) {
  double rc = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 8) {
      continue;
    }

    if (per_molecule_data.ncon[i] == 1) {
      rc += ClassifySinglyConnectedOxygen(m, i, per_molecule_data, status);
    } else if (per_molecule_data.ncon[i] == 2) {
      rc += ClassifyDoublyConnectedOxygen(m, i, per_molecule_data, status);
    } else {
      status[i] = kFailed;
      if (display_unclassified_atom_messages) {
        cerr << "xlogp::IdentifyOxygen:unrecognised oxygen " << m.name() << '\n';
      }
    }
  }

  return rc;
}


// This is only defined for furan types, but we do not check that.
std::optional<double>
ClassifyAromaticOxygen(Molecule& m,
                        atom_number_t zatom,
                        PerMoleculeData& per_molecule_data,
                        int * status) {
  status[zatom] = kOxygenAromatic;
  return 0.327;  // 43
}

double
IdentifyAromaticAtoms(Molecule& m,
                PerMoleculeData& per_molecule_data,
                int* status) {
  const int matoms = m.natoms();
  double rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (! per_molecule_data.aromatic[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] == 6) {
      std::optional<double> x = ClassifyAromaticCarbon(m, i, per_molecule_data, status);
      if (x) {
        rc += *x;
      } else {
        status[i] = kFailed;
      }
    } else if (per_molecule_data.atomic_number[i] == 7) {
      std::optional<double> x = ClassifyAromaticNitrogen(m, i, per_molecule_data, status);
      if (x) {
        rc += *x;
      } else {
        status[i] = kFailed;
      }
    } else if (per_molecule_data.atomic_number[i] == 8) {
      std::optional<double> x = ClassifyAromaticOxygen(m, i, per_molecule_data, status);
      if (x) {
        rc += *x;
      } else {
        status[i] = kFailed;
      }
    } else if (per_molecule_data.atomic_number[i] == 16) {
      status[i] = kSulphurArom;
      rc += 0.964;  // 68
    } else {
      if (display_unclassified_atom_messages) {
        cerr << "xlogp::IdentifyAromaticAtoms:unrecognised aromatic type " << m.name() << 
             m.smarts_equivalent_for_atom(i) << '\n';
        status[i] = kFailed;
      }
    }
  }

  return rc;
}

double
IdentifyCH3(Molecule& m,
                PerMoleculeData& per_molecule_data,
                int* status) {
  const int matoms = m.natoms();
  double rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 6 ||
        per_molecule_data.ncon[i] != 1 ||
        per_molecule_data.unsaturation[i] != 0) {
      continue;
    }

    if (per_molecule_data.attached_carbons[i] == 0) {
      status[i] = kCh3X;
      rc += -0.181;
      continue;
    }

    const atom_number_t o = m.other(i, 0);
    if (per_molecule_data.pi_electrons[o]) {
      status[i] = kCH3Pi;
      rc += 0.168;
    } else {
      status[i] = kCH3NoPi;
      rc += 0.484;
    }
  }

  return rc;
}

double
IdentifyCH2(Molecule& m,
                PerMoleculeData& per_molecule_data,
                int* status) {
  const int matoms = m.natoms();
  double rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 6 ||
        per_molecule_data.ncon[i] != 2 ||
        per_molecule_data.hcount[i] != 2) {
      continue;
    }


    if (per_molecule_data.attached_carbons[i] == 2) {
      if (per_molecule_data.NbrsHavePiElectrions(m, i)) {
        status[i] = kCh2R2Pi;
        rc += 0.009;  // 5
      } else {
        status[i] = kCh2R2NoPi;
        rc += 0.359; //  4
      }
    } else if (per_molecule_data.attached_carbons[i] == 1) {
      status[i] = kCh2RX;
      rc += -0.344;
      continue;
    } else {  // 0 attached carbons
      status[i] = kCh2X2;
      rc += -0.439;
    }
  }

  return rc;
}

double
IdentifyCH(Molecule& m,
           PerMoleculeData& per_molecule_data,
           int* status) {
  const int matoms = m.natoms();
  double rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 6 ||
        per_molecule_data.ncon[i] != 3 ||
        per_molecule_data.hcount[i] != 1) {
      continue;
    }

    if (per_molecule_data.attached_carbons[i] == 3) {
      if (per_molecule_data.NbrsHavePiElectrions(m, i)) {
        status[i] = kChR3Pi;
        rc += -0.138;  // 9
      } else {
        status[i] = kChR3NoPi;
        rc += 0.051;  // 8
      }
    } else if (per_molecule_data.attached_carbons[i] == 2) {
      status[i] = kChR2X;
      rc += -0.417;
    } else {  // 0 or 1 attached_carbons
      status[i] = kChRX2;
      rc += -0.454;
    }
  }

  return rc;
}

double
IdentifyCD4(Molecule& m,
           PerMoleculeData& per_molecule_data,
           int* status) {
  const int matoms = m.natoms();
  double rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }
    if (per_molecule_data.atomic_number[i] != 6 ||
        per_molecule_data.ncon[i] != 4 ||
        per_molecule_data.hcount[i] != 0) {
      continue;
    }

    if (per_molecule_data.attached_carbons[i] == 4) {
      if (per_molecule_data.NbrsHavePiElectrions(m, i)) {
        status[i] = kCR4Pi;
        rc += 0.223;  // 13
      } else {
        status[i] = kCR4NoPi;
        rc += -0.378;  // 12
      }
    } else if (per_molecule_data.attached_carbons[i] == 3) {
      status[i] = kCR3X;
      rc += -0.598;  // 14
    } else if (per_molecule_data.attached_carbons[i] == 2) {
      status[i] = kCR2X2;
      rc += -0.396;  // 15 seems poorly parameterized, only 4 molecules in WFL
    } else if (per_molecule_data.attached_carbons[i] == 1) {
      status[i] = kCRX3;
      rc += -0.699;  // 16
    } else {  // 0 attached carbons
      status[i] = kCX4;
      rc += -0.362;  // 17
    }
  }

  return rc;
}

std::optional<double>
ClassifyCarbonSp(Molecule& m,
                 atom_number_t zatom,
                 PerMoleculeData& per_molecule_data,
                 int* status) {
  atom_number_t end_of_triple_bond = INVALID_ATOM_NUMBER;
  for (const Bond* b : m[zatom]) {
    if (! b->is_triple_bond()) {
      continue;
    }
    end_of_triple_bond = b->other(zatom);
    break;
  }

  if (end_of_triple_bond != INVALID_ATOM_NUMBER &&
      per_molecule_data.ncon[zatom] == 1 &&
      per_molecule_data.atomic_number[end_of_triple_bond] == 6) {
    status[zatom] = kCarbonSpRCH;
    return type_to_score[kCarbonSpRCH];  // 35
  } 

  status[zatom] = kCarbonSpRCR;
  return type_to_score[kCarbonSpRCR];  // 36
}

double
IdentifyCsp(Molecule& m,
            PerMoleculeData& per_molecule_data,
            int* status) {
  double rc = 0.0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }


    if (per_molecule_data.atomic_number[i] != 6 ||
        per_molecule_data.unsaturation[i] != 2) {
      continue;
    }

    std::optional<double> x = ClassifyCarbonSp(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      status[i] = kFailed;
    }
  }

  return rc;
}

std::optional<double>
ClassifyCarbonSp2(Molecule& m,
                  atom_number_t zatom,
             PerMoleculeData& per_molecule_data,
             int* status) {
  if (per_molecule_data.ncon[zatom] == 1 &&
      per_molecule_data.hcount[zatom] == 2 &&
      per_molecule_data.attached_carbons[zatom] == 1) {
    status[zatom] = kCarbonSp2RCH2;
    return 0.395;  // 18
  }

  if (per_molecule_data.ncon[zatom] == 2 &&
      per_molecule_data.hcount[zatom] == 1 &&
      per_molecule_data.attached_carbons[zatom] == 2) {
    status[zatom] = kCarbonSp2RdblCHR;
    return 0.235;  // 19
  }

  if (per_molecule_data.attached_carbons[zatom] == 3) {
    status[zatom] = kCarbonSp2RdblCR2;
    return 0.098;  // 22
  }

  if (per_molecule_data.ncon[zatom] == 3 &&
      per_molecule_data.attached_carbons[zatom] == 0) {
    status[zatom] = kCarbonSp2XdblCX2;
    return 1.774; // 25
  }

  // For all the other types, we need to know where the heteratoms are.

  atom_number_t dbl = INVALID_ATOM_NUMBER;
  atom_number_t single1 = INVALID_ATOM_NUMBER;
  atom_number_t single2 = INVALID_ATOM_NUMBER;
  for (const Bond* b : m[zatom]) {
    atom_number_t o = b->other(zatom);
    if (b->is_double_bond()) {
      dbl = o;
    } else if (single1 == INVALID_ATOM_NUMBER) {
      single1 = o;
    } else {
      single2 = o;
    }
  }

  (void) single2;  // we don't really need it.

  // 20 and 23
  if (per_molecule_data.atomic_number[dbl] == 6) {
    if (per_molecule_data.attached_carbons[zatom] == 1 &&
        per_molecule_data.hcount[zatom] == 1) {
      status[zatom] = kCarbonSp2RdblCHX;
      return type_to_score[kCarbonSp2RdblCHX];  // 20
    } else {
      status[zatom] = kCarbonSp2RdbkCRX;
      return type_to_score[kCarbonSp2RdbkCRX];  // 23
    }
  } else {
    if (per_molecule_data.hcount[zatom] == 1) {
      status[zatom] = kCarbonSp2XdblCHR;
      return type_to_score[kCarbonSp2XdblCHR];  // 21
    } else {
      status[zatom] = kCarbonSp2XdblCR2;
      return type_to_score[kCarbonSp2XdblCR2];  // 24
    }
  }

  if (display_unclassified_atom_messages) {
    cerr << "xlogp::ClassifyCarbonSp2:unrecognised type " << m.name() << '\n';
  }

  return std::nullopt;
}

double
IdentifyCSp2(Molecule& m,
             PerMoleculeData& per_molecule_data,
             int* status) {
  double rc = 0.0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (status[i]) {
      continue;
    }

    if (per_molecule_data.atomic_number[i] != 6 ||
        per_molecule_data.aromatic[i]) {
      continue;
    }
    if (per_molecule_data.unsaturation[i] != 1) {
      continue;
    }

    std::optional<double> x = ClassifyCarbonSp2(m, i, per_molecule_data, status);
    if (x) {
      rc += *x;
    } else {
      status[i] = kFailed;
    }
  }

  return rc;
}

double
IdentifyTerminalGroups(Molecule& m,
                PerMoleculeData& per_molecule_data,
                int* status) {
  double rc = 0.0;
  // Terminal forms first.
  rc += IdentifyCyano(m, per_molecule_data, status);
  rc += IdentifyNCS(m, per_molecule_data, status);
  rc += IdentifyNO(m, per_molecule_data, status);
  rc += IdentifyNO2(m, per_molecule_data, status);

  return rc;
}

// Runs under the assumption that WFL aromaticity has been set.
std::optional<double>
XLogPWFL(Molecule& m,
      int* status) {
  const int matoms = m.natoms();

  std::fill_n(status, matoms, 0);

  PerMoleculeData per_molecule_data(m);

  double rc = IdentifyTerminalGroups(m, per_molecule_data, status);

  rc += IdentifyAromaticAtoms(m, per_molecule_data, status);

  rc += IdentifyHalogens(m, per_molecule_data, status);
  rc += IdentifyOxygen(m, per_molecule_data, status);
  rc += IdentifySulphur(m, per_molecule_data, status);
  rc += IdentifyPhosphorus(m, per_molecule_data, status);

  rc += IdentifyCH3(m, per_molecule_data, status);
  rc += IdentifyCH2(m, per_molecule_data, status);
  rc += IdentifyCH(m, per_molecule_data, status);
  rc += IdentifyCD4(m, per_molecule_data, status);

  // Do rare carbon forms first. Inefficient, but helpful
  rc += IdentifyCsp(m, per_molecule_data, status);
  rc += IdentifyCSp2(m, per_molecule_data, status);

  // Identify certain very specific nitrogen forms first.
  rc += IdentifyAmideNitrogens(m, per_molecule_data, status);
  rc += IdentifyTrigonalPlanarNitrogenAromatic(m, per_molecule_data, status);
  rc += IdentifyTrigonalPlanarNitrogen(m, per_molecule_data, status);
  rc += IdentifyNitrogenAromatic(m, per_molecule_data, status);
  rc += IdentifyChargedNitrogen(m, per_molecule_data, status);
  rc += IdentifyNitrogenSP3(m, per_molecule_data, status);
  rc += IdentifyNitrogenSP2(m, per_molecule_data, status);

  rc += m.implicit_hydrogens() * type_to_score[kHydrogen];  // 37

  rc += Corrections(m, per_molecule_data, status);

  std::unordered_map<int, int> count;

  int classified = 0;
  int unclassified = 0;
  int failed = 0;
  Molecule mcopy(m);
  for (int i = 0; i < matoms; ++i) {
    if (status[i] > 0) {
      ++classified;
      // cerr << m.smarts_equivalent_for_atom(i) << " assigned " << status[i] << '\n';
      mcopy.set_isotope(i, status[i]);
      ++count[status[i]];
    } else if (status[i] == 0) {
      ++unclassified;
    } else if (status[i] == kFailed) {
      ++failed;
    }
  }

  if (! display_assignments) {
    return rc;
  }

  cerr << m.name() << ' ' << classified << " classified " << unclassified << " unclassified " << failed << " failed " << rc << '\n';
  cerr << mcopy.aromatic_smiles() << ' ' << m.name() << '\n';
  double sum = 0.0;
  for (const auto& [k, v] : count) {
    cerr << " type " << k << " count " << v << ' ' << type_to_score[k] << ' ' << (v * type_to_score[k]) << '\n';
    sum += (v * type_to_score[k]);
  }

  const int ih = m.implicit_hydrogens();
  cerr << ih << " imph " << type_to_score[kHydrogen] << ' ' << (ih * type_to_score[kHydrogen]) << '\n';
  sum += (ih * type_to_score[kHydrogen]);
  cerr << sum << '\n';
  for (int i = 0; i < matoms; ++i) {
    if (i > 0) {
      cerr << ", ";
    }
    cerr << status[i];
  }
  cerr << "\n";


  return rc;
}

std::optional<double>
XLogP(Molecule& m,
      int* status) {

  const auto aromsave = global_aromaticity_type();
  // Aromaticity definitions seem quite problematic, need more test
  // cases to figure out what is closest to their implementation.
  set_global_aromaticity_type(WangFuLai);
  m.compute_aromaticity();

  auto rc = XLogPWFL(m, status);

  set_global_aromaticity_type(aromsave);
  // Potentially wasted computation, need to do something about this...
  m.compute_aromaticity();

  return rc;
}

std::optional<double>
XLogP(Molecule& m) {
  std::unique_ptr<int[]> status = std::make_unique<int[]>(m.natoms());

  return XLogP(m, status.get());
}

} // namespace xlogp
