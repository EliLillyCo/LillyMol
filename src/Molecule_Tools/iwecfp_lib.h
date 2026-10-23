#ifndef MOLECULE_TOOLS_IWECFP_LIB_H_
#define MOLECULE_TOOLS_IWECFP_LIB_H_

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <unordered_map>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

namespace iwecfp {

// Return codes from Iwecfp::Fingerprint()
enum class FingerprintResult {
  kOk,
  kNoStartAtoms,
  kFatal
};

using atype_t = uint32_t;

class Iwecfp {
 private:
  int _verbose = 0;

  int _min_shell_radius = 0;
  int _max_shell_radius = 0;

  int _additive = 1;

  int _each_shell_gets_own_fingerprint = 0;
  int _all_bonds_same_type = 0;

  std::unordered_map<uint32_t, uint32_t> _bits_to_investigate;
  int _looking_for_bit_meanings = 0;
  int _write_smiles_with_bit_meanings = 0;
  uint32_t _bits_found = 0;

  IWString_and_File_Descriptor _stream_for_bit_meanings;
  IWString_and_File_Descriptor _stream_for_all_bits;

  // Experimental idea. Keeps track of how often each atom
  // participates in bit formation. At the end, if certain atoms
  // have been "ignored", generate some more bits centred on 
  // those atoms. Not sure this is useful.
  int _equalise_atom_coverage = 0;
  int _label_by_visited = 0;

  // If set, shells are only started at atoms that match a query.
  resizable_array_p<Substructure_Query> _start_atom_query;

  // If the centre atom of each shell has an isotopic label
  // use that isotopic value to form a chiral aware shell.
  int _central_atom_possible_chiral = 0;

  // Used by the methods that write bit meanings.
  IWDigits _iwdigits_center;
  IWDigits _iwdigits;

  // Working data for the molecule being processed.
  // This object is not thread safe.
  Molecule* _current_molecule = nullptr;
  atom_number_t _centre_of_shell = kInvalidAtomNumber;
  IWString _smarts_for_centre_of_shell;
  isotope_t _centre_atom_isotope = 0;

  // Experimental idea. As shells are formed, add bonds to the
  // previous shell as a "tail". Not sure the implementation is
  // correct. Needs work, do not use for production work.
  int _add_tails = 0;

 public:
  Iwecfp();
  ~Iwecfp();

  int Initialise(Command_Line& cl);

  FingerprintResult Fingerprint(Molecule& m, const atype_t* atom_constant,
                  Sparse_Fingerprint_Creator* sfc);

  // Generate a fingerprint for a strict subset of `m`. Shell centres are
  // restricted to atoms where `include_atom[i]` is non-zero, and shell
  // expansion never crosses to atoms where `include_atom[i]` is zero.
  // Atomic properties are taken from `atom_constant` exactly as supplied.
  FingerprintResult Fingerprint(Molecule& m, const atype_t* atom_constant,
                  const int* include_atom, Sparse_Fingerprint_Creator* sfc);

  int Finalise();
  int Report(std::ostream& output) const;

  int verbose() const { return _verbose; }
  int max_radius() const { return _max_shell_radius;}
  void set_max_radius(int s) { _max_shell_radius = s;}
  int each_shell_gets_own_fingerprint() const { return _each_shell_gets_own_fingerprint;}

 private:
  int ReadBitsToInvestigate(iwstring_data_source& input);
  int ReadBitsToInvestigate(const char* fname);

  int CheckAgainstList(Molecule& m, const IWString& smarts,
                       atom_number_t centre_of_shell, unsigned int sum_so_far,
                       int radius);
  void WriteLabelledSmiles(const Molecule& m, int centre_of_shell, int radius,
                           IWString_and_File_Descriptor& output);
  void WriteBit(int centre_of_shell,
                const IWString& smarts_for_centre_of_shell,
                int radius, unsigned int bit, Molecule& m,
                IWString_and_File_Descriptor& output);

  int BondConstant(const Bond* bond) const;
  void Increment(unsigned int& sum_so_far, int bc, atype_t atom_constant) const;

  int GenerateShells(int matoms, int radius, int max_radius,
                     const Atom* const* atoms, const atype_t* atom_constant,
                     int* processing_status, unsigned int sum_so_far,
                     Molecule& m, Sparse_Fingerprint_Creator* sfc);
  int GenerateShellsSubset(int matoms, int radius, int max_radius,
                     const Atom* const* atoms, const atype_t* atom_constant,
                     const int* include_atom, int* processing_status,
                     unsigned int sum_so_far, Molecule& m,
                     Sparse_Fingerprint_Creator* sfc);

  int AddFingerprint(const Sparse_Fingerprint_Creator& from,
                     Sparse_Fingerprint_Creator& to);
  void FormBit(Molecule& m, const atype_t* atom_constant,
               const Atom* const* atoms, atom_number_t zatom, int max_r,
               int* processing_status, Sparse_Fingerprint_Creator* sfc);
  void FormBitSubset(Molecule& m, const atype_t* atom_constant,
               const Atom* const* atoms, const int* include_atom,
               atom_number_t zatom, int max_r, int* processing_status,
               Sparse_Fingerprint_Creator* sfc);

  void IdentifyAtomsWithinRange(Molecule& m, Set_of_Atoms* atoms_within_range);
  int DoEqualiseAtomCoverage(Molecule& m, const Atom* const* atoms,
                             const atype_t* atom_constant,
                             int* processing_status,
                             Sparse_Fingerprint_Creator* sfc);

  int IdentifyStartAtoms(Molecule& m, int* processing_status, int matched_flag);
};

void Usage(int rc);
void DisplayDashGOptions(std::ostream& output);
void DisplayDashYOptions(std::ostream& output);

}  // namespace iwecfp

#endif  // MOLECULE_TOOLS_IWECFP_LIB_H_
