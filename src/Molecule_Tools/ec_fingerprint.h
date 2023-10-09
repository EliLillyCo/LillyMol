#ifndef EC_FINGERPRINT_H_
#define EC_FINGERPRINT_H_

#include <algorithm>
#include <cstdint>
#include <tuple>
#include <unordered_set>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/fingerprint_writer.h"

namespace ec_fingerprint {

#define NOT_PROCESSED 0
#define EXCLUDED 1
#define PROCESSING_COMPLETE 2
#define CURRENT_SHELL 3

using atom_type_t = uint32_t;

using count_type_t = uint32_t;

// This will be filled by ec_fingerprint_main.
struct JobParameters {
  JobParameters();

  fingerprint_writer::FingerprintWriter fp_writer;

  IWString smiles_tag;
  IWString identifier_tag;

  // The atom type information is written to the precedent data file.
  // There should be automatic enforcement of consistency when it is used.
  IWString atom_type_string;

  // Has the user requested output or not.
  bool produce_output;

  // Reading fingerprints from stdin.
  bool function_as_tdt_filter;

  // If positive, the width of the fingerprint to produce.
  int write_fixed_width_fingerprint;

  bool write_counted_sparse_fingerprint;
};

// Class holds transient results of a shell expansion. It largely exists
// to cut down on the number of arguments getting passed between functions.
class ShellInfo {

  public:
    // Reference to molecule being processed and its atom count.
    const Molecule& _m;
    const int _matoms;

    // The atom types, from the caller.
    const atom_type_t * _atom_type;

    // Atom at centre of each shell
    atom_number_t _a0;

    // Processing status of atoms during bit formation.
    int * _status;

    // The atoms that comprise the frontier shell.
    Set_of_Atoms _next_shell;

  public:
    ShellInfo(const Molecule& m, const int* include_atom, const atom_type_t* atom_type);
    ~ShellInfo();

    const Molecule& m() const {
      return _m;
    }

    const IWString& name() const { 
      return _m.name();
    }

    atom_number_t a0() const {
      return _a0;
    }

    int status(const atom_number_t a) const {
      return _status[a];
    }

    const atom_type_t* atom_type() const { return _atom_type;}
    atom_type_t atom_type(const atom_number_t a) const {
      return _atom_type[a];
    }

    void SetCentreAtom(const atom_number_t s);

    void AddToNextShell(atom_number_t s) {
      _next_shell.add(s);
      _status[s] = CURRENT_SHELL;
    }

    const Set_of_Atoms& NextShell() const { return _next_shell;}
    int NextShellSize() const { return _next_shell.number_elements();}

    // The atoms in the next shell are marked in the status array
    // and the _next_shell vector is cleared, waiting for the next
    // layer to be found.
    void GrabNextShell(Set_of_Atoms& next_shell);
};

// The ECFingerprint generator class takes a template function, that must support
// certain methods. This is the virtual base class.
class ECFunction
{
  protected:
    // Called just before the object begins processing a new molecule.
    virtual int PrepareToProcess(Molecule& m) = 0;

    // As new bits are discovered, they are provided via the Bit method.
    virtual void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius) = 0;

    // Called once fingerprint generation is complete.
    virtual int FingerprintingComplete(Molecule& m) = 0;

    // If the object does any per molecule output.
    virtual int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                            IWString_and_File_Descriptor& output) = 0;
};

// A number of the classes that get passed to the templatized EC fingerprint
// generator need an output file. A common base class that provides that, as
// well as an Open method.
class ECBaseWithOutput
{
  protected:
    IWString_and_File_Descriptor _output;

  public:
    int Open(IWString& fname);
};

// Class that produces a fingerprint. This is the most common use case.
class ProduceFingerprint : public ECFunction
{
  private:
    Sparse_Fingerprint_Creator _sfc;

    // private functions
    int WriteFixedWidthFingerprint(const JobParameters& job_parameters, IWString_and_File_Descriptor& output) const;

  public: 
    ProduceFingerprint() {
    }

    int PrepareToProcess(Molecule& m) override;

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius) override {
      _sfc.hit_bit(running_sum);
    }

    int FingerprintingComplete(Molecule& m) override {
      return 1;
    }
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters, 
                    IWString_and_File_Descriptor& output) override;

    const Sparse_Fingerprint_Creator& sfc() const { return _sfc;}
};

// Puts the atom coverage on the atom map number of each atom.
class AtomMapCoverage : public ECBaseWithOutput, public ECFunction
{
  private:
    // Over time, the _coverage array will resize to the largest number
    // of atoms encountered. But that will avoid reallocations.
    int _matoms;
    int * _coverage;

  public:
    AtomMapCoverage();
    ~AtomMapCoverage();

    int PrepareToProcess(Molecule& m);

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius) {
      shell_info.NextShell().increment_vector(_coverage);
    }

    int FingerprintingComplete(Molecule& m);
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output) {
      return 1;
    }
};

// Class WriteAllBits writes all bits encountered to a stream.
class WriteAllBits : public ECBaseWithOutput, public ECFunction
{
  private:
  public:
    WriteAllBits() {
    }

    int PrepareToProcess(Molecule& m) {
      return 1;
    }

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius);

    int FingerprintingComplete(Molecule& m) {
      return 1;
    }
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output) {
      return 1;
    }
};

// Look for a specific set of bits and write an explanation.
// Not implemented...
class ECBitMeanings : public ECBaseWithOutput, public ECFunction
{
  private:
    std::unordered_set<atom_type_t> _bits_to_find;

    int _bits_found;
    int _bits_not_found;

    IWString _buffer_current_molecule;

  // private functions

    int _ReadBitsToFind(iwstring_data_source& input);

  public:
    ECBitMeanings();
    ~ECBitMeanings();

    int ReadBitsToFind(IWString& fname);

    // 
    int  ParseOptionalSettings(const IWString & args);

    int PrepareToProcess(Molecule& m);

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius);

    int FingerprintingComplete(Molecule& m) {
      return 1;
    }
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output) ;
};

// When adjusting magic values in the ECFingerprint object, studying collisions
// can be useful.
class ECCheckCollisions : public ECBaseWithOutput
{
  private:
// If checking for collisions, there needs to be a data structure that holds a
// mapping from bit number to the first molecule to report this bit number, the
// atom type of the centre atom, the radius, and the number of instances.
    std::unordered_map<atom_type_t, std::tuple<IWString, atom_type_t, int>> _bit_to_description;

    // Across all molecules encountered.
    int _collisions_found;

  public:
    ECCheckCollisions() {
      _collisions_found = 0;
    }
    ~ECCheckCollisions() {
    }

    int PrepareToProcess(Molecule& m) {
      return 1;
    }

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius);

    int FingerprintingComplete(Molecule& m) {
      return 1;
    }
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output) {
      return 1;
    }
};

// Precedent data is built up by studying a set of existing molecules.
class ECBuildPrecedent
{
  private:
    // The data that is accumulated during builds.
    struct BitCount {
      // The smiles of the first exemplar structure.
      IWString smiles;
      // Atom type of the center atom.
      atom_type_t center_atom_type;
      // The radius
      int radius;
      // The number of instances found.
      count_type_t count;

      BitCount(const IWString& s, atom_type_t a, int r, count_type_t c) {
        smiles = s;
        center_atom_type = a;
        radius = r;
        count = c;
      }
    };

    // If recording the provenance of of bits formed, we need to ma the bit number to
    // something recording the first molecule, atom_type of centre atom, radius and count.

    std::unordered_map<atom_type_t, BitCount> _bit_count;

    // During construction is can be useful to keep track of the number of bit
    // collisions avoided.

    int _bit_collisions_avoided;
    // The radii at which collisions are observed.
    extending_resizable_array<int> _collision_at_radius;

  public:
    ECBuildPrecedent();
    ~ECBuildPrecedent();

    void PrepareToProcess(Molecule& m) {
      return;
    }

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius);

    int FingerprintingComplete(Molecule& m) {
      return 1;
    }
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output) {
      return 1;
    }

    // Once data is accumulated, write to file.
    int WritePrecedentData(const char sep, const JobParameters& job_parameters, IWString& fname) const;

    int Report(std::ostream& output) const;
};

// Once a set of precedent data has been built, class ECUsePrecedent can use that
// data to make precedent assessments of new molecules.
class ECUsePrecedent : public ECFunction
{
  private:
    // For precedent calculations, keep the number of occurrences and radius of each bit
    // read from the previously generated precedent data.
    struct BitCount {
      // The number of occurrences of this bit in the knowledge base.
      count_type_t count;
      // The radius.
      int radius;
      BitCount(count_type_t c, int r) {
        count = c;
        radius = r;
      }
    };

    std::unordered_map<atom_type_t, BitCount> _precedent;

    // The longest radius used in the current calculation.
    // Note that this might be different from what was used to generate
    // the precedent data. Clearly if this is longer, there will be
    // erroneous results. Need some means of automatically enforcing this.
    const int _max_radius;

    // Various arrays used during processing each molecule.
    // No this is not thread safe.
    // At each radius, the smallest precedent value
    count_type_t * _count;
    // At each radius, the centre atom of the rare shell.
    atom_number_t* _atom;

    // A count of the number of molecules that have a missing shell at
    // each radius.
    count_type_t* _missing_at_radius;

  // private functions
    
    // Read record generated by ECBuildPrecedent.
    int _ParsePrecedentRecord(const const_IWSubstring& line);

  public:
    ECUsePrecedent(const int max_radius);
    ~ECUsePrecedent();

    // Read a previously generated precedent data set - from ECBuildPrecedent.
    int ReadPrecedentData(IWString& fname);

    int PrepareToProcess(Molecule& m);

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius);

    // Once all bits have been generated, and precedent data has been collected, 
    int FingerprintingComplete(Molecule& m);
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output);
    int Report(std::ostream& output) const;
};

// A class that examines each molecule and if any of the bits in `_bits_to_find`
// are generated, that molecule is not written.
class ECFilterByBits : public ECFunction
{
  private:
    std::unordered_set<uint32_t> _bits_to_find;

    bool _write_current_molecule;

  public:
    ECFilterByBits();
    ~ECFilterByBits() {}

    // Read a previously generated file of bits to look for.
    int ReadBitsToFilter(IWString& fname);

    int PrepareToProcess(Molecule& m) {
      _write_current_molecule = true;
      return 1;
    }

    void Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius);

    // Once all bits have been generated, and precedent data has been collected, 
    int FingerprintingComplete(Molecule& m) {
      return 1;
    }
    int DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                    IWString_and_File_Descriptor& output);
    int Report(std::ostream& output) const;
};

// Class for generating EC fingerprints.
// Apart from various public set_* methods, the only public function is Fingerprint.
// The template function needed must have the same function signature as class ECFunction.
class ECFingerprint {
  private:
  private:
    int _min_radius;
    int _max_radius;

    // Bits can be incremented either via an additive, or multiplicative
    // means. If multiplicative, then information about where a feature
    // occurs in the shells is lost.
    bool _additive_shell_formation;

    // When the shell expands, and we form a new bit, it is optional
    // as to whether that bit includes the info about the atom to
    // which the new atom is attached. If _precise is true, that 
    // information will be included.
    bool _precise;

    // Settable parameters that can change the bit numbers of what is produced.
    atom_type_t _magic1;
    atom_type_t _magic2;
    atom_type_t _magic3;
    atom_type_t _magic4;
    atom_type_t _magic5;

    atom_type_t _bond_magic1;
    atom_type_t _bond_magic2;
    atom_type_t _bond_magic3;
    atom_type_t _bond_magic4;

//  Private functions.

    // A mirror of the public method, but the argument is const.
    template <typename T>
    int _Fingerprint(const Molecule& m,
        ShellInfo& shell_info,
        T& accumulate_bits);

    
    template <typename T>
    void _ExpandShell(ShellInfo& shell_info, atom_type_t running_sum, int radius, T& accumulate_bits);

    atom_type_t _BondConstant(const Bond& b) const;

    void _AddToRunningSum(ShellInfo& shell_info,
                 const atom_number_t a1,
                 const atom_type_t bond_constant,
                 const atom_number_t a2,
                 atom_type_t& running_sum) const;

    // If 'radius' is in range, set the bit. Call _ExamineBit if needed.
    template <typename T>
    void _MaybeRecordBit(ShellInfo& shell_info, const atom_type_t running_sum, const int radius, T& accumulate_bits);

    // If _max_radius is zero, process separately.
    template <typename T>
    int _DoZeroRadiusOnly(const Molecule& m, ShellInfo& shell_info, T& accumulate_bits);

  public:
    ECFingerprint();

    void set_min_radius(int s) { _min_radius = s;}
    void set_max_radius(int s) { _max_radius = s;}

    void set_precise(bool s) { _precise = s;}
    void set_additive(bool s) { _additive_shell_formation = s; }

    // T must be an object that supports a method
    // void T::Bit(const ShellInfo& shell_info, const atom_number_t running_sum, const int radius)
    // As each bit is formed, it is passed to the T object, which can do whatever it wants
    // with that information.
    // There must be two other methods PostMoleculeProcessing(Molecule& m);
    template <typename T>
    int Fingerprint(Molecule& m,
        const int* include_atom,
        const atom_type_t* atom_type,
        T& accumulate_bits);
};

//#define DEBUG_EC_FINGERPRINT

template <typename T>
int
ECFingerprint::Fingerprint(Molecule& m,
        const int* include_atom,
        const atom_type_t* atom_type,
        T& accumulate_bits) {
  const int matoms = m.natoms();
  if (0 == matoms)
    return 0;

#ifdef DEBUG_EC_FINGERPRINT
  for (int i = 0; i < matoms; ++i) {
    cerr << "Atom " << i << " atype " << atom_type[i] << endl;
  }
#endif

  m.compute_aromaticity_if_needed();

  accumulate_bits.PrepareToProcess(m);

  ShellInfo shell_info(m, include_atom, atom_type);

  return _Fingerprint(m, shell_info, accumulate_bits);
}

template <typename T>
int
ECFingerprint::_Fingerprint(const Molecule& m,
        ShellInfo& shell_info,
        T& accumulate_bits)
{
  // Process separately so we do not need to check inside the loop.
  if (0 == _max_radius)
    return _DoZeroRadiusOnly(m, shell_info, accumulate_bits);

  const int matoms = m.natoms();

  const atom_type_t* atom_type = shell_info.atom_type();

  for (int i = 0; i < matoms; ++i) {
    if (EXCLUDED == shell_info.status(i))
      continue;

    shell_info.SetCentreAtom(i);

    atom_type_t running_sum = atom_type[i] + _magic1;
    _MaybeRecordBit(shell_info, running_sum, 0, accumulate_bits);

    const Atom * a = m.atomi(i);

    for (const Bond * b : *a)
    {
      const atom_number_t k = b->other(i);
      if (EXCLUDED == shell_info.status(k))
        continue;

      if (CURRENT_SHELL == shell_info.status(k))  // Loop.
        continue;

      _AddToRunningSum(shell_info, i, _BondConstant(*b), k, running_sum);
#ifdef DEBUG_EC_FINGERPRINT
      cerr << "From atom " << i << " add atom " << k << endl;
#endif
      shell_info.AddToNextShell(k);
    }

    if (shell_info.NextShell().number_elements())
      _MaybeRecordBit(shell_info, running_sum, 1, accumulate_bits);

    if (_max_radius >= 2)
      _ExpandShell(shell_info, running_sum, 1, accumulate_bits);
  }

  return 1;
}

template <typename T>
void
ECFingerprint::_MaybeRecordBit(ShellInfo& shell_info,
                               const atom_type_t running_sum,
                               const int radius,
                               T& accumulate_bits)
{
  if (radius < _min_radius)
    return;

#ifdef DEBUG_EC_FINGERPRINT
    cerr << "Recording bit " << running_sum << " at radius " << radius << " centre " << shell_info.a0() << endl;
#endif

  accumulate_bits.Bit(shell_info, running_sum, radius);

  return;
}

template <typename T>
void
ECFingerprint::_ExpandShell(ShellInfo& shell_info,
                            atom_type_t running_sum,
                            const int radius,
                            T& accumulate_bits)
{
  const Set_of_Atoms& current_shell = shell_info.NextShell();
  if (current_shell.empty())
    return;

#ifdef DEBUG_EC_FINGERPRINT
  cerr << "_ExpandShell from " << shell_info.a0() << " radius " << radius << " atoms " << shell_info.NextShell() << " sum " << running_sum << endl;
#endif

  if (_additive_shell_formation)
    running_sum = _magic3 * running_sum + _magic2;

#ifdef DEBUG_EC_FINGERPRINT
  cerr << "running_sum updated to " << running_sum << endl;
#endif

  Set_of_Atoms next_shell;
  next_shell.resize(12);

  for (const atom_number_t i : current_shell) {
    assert (CURRENT_SHELL == shell_info.status(i));
    const Atom * a = shell_info.m().atomi(i);

    for (const Bond* b : *a) {
      const atom_number_t k = b->other(i);
      const int kstatus = shell_info.status(k);
      if (NOT_PROCESSED == kstatus || CURRENT_SHELL == kstatus) {
#ifdef DEBUG_EC_FINGERPRINT
        cerr << "_ExpandShell going to atom " << k << " from " << shell_info.a0() << " radius " << radius << endl;
#endif
        _AddToRunningSum(shell_info, i, _BondConstant(*b), k, running_sum);
        if (CURRENT_SHELL != kstatus)
          next_shell.add_if_not_already_present(k);
      }
    }
  }

  if (next_shell.empty())
    return;

  _MaybeRecordBit(shell_info, running_sum, radius + 1, accumulate_bits);

  if (radius + 1 >= _max_radius)  // We have just processed (radius+1).
    return;

  shell_info.GrabNextShell(next_shell);

  _ExpandShell(shell_info, running_sum, radius + 1, accumulate_bits);
}

template <typename T>
int
ECFingerprint::_DoZeroRadiusOnly(const Molecule& m,
                        ShellInfo& shell_info,
                        T& accumulate_bits)
{
  const int matoms = m.natoms();

  const atom_type_t* atom_type = shell_info.atom_type();

  for (int i = 0; i < matoms; ++i) {
    if (EXCLUDED == shell_info.status(i))
      continue;

    shell_info.SetCentreAtom(i);

    atom_type_t running_sum = atom_type[i] + _magic1;
    _MaybeRecordBit(shell_info, running_sum, 0, accumulate_bits);
  }

  return 1;
}

}  // namespace ec_fingerprint

#endif  // EC_FINGERPRINT_H_
