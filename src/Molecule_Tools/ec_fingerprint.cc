#include "ec_fingerprint.h"

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/data_source/iwstring_data_source.h"

namespace ec_fingerprint {

using std::cerr;

ShellInfo::ShellInfo(const Molecule& m, const int* include_atom,
                     const atom_type_t* atom_type)
    : _m(m), _matoms(m.natoms()), _atom_type(atom_type) {
  _status = new int[_matoms];
  if (nullptr == include_atom) {
    std::fill_n(_status, _matoms, NOT_PROCESSED);
  } else {
    for (int i = 0; i < _matoms; ++i) {
      if (include_atom[i]) {
        _status[i] = NOT_PROCESSED;
      } else {
        _status[i] = EXCLUDED;
      }
    }
  }

  _next_shell.make_room_for_extra_items(_matoms);

  return;
}

ShellInfo::~ShellInfo() {
  delete[] _status;
}

void
ShellInfo::SetCentreAtom(const atom_number_t a) {
  _a0 = a;

  for (int i = 0; i < _matoms; ++i) {
    if (EXCLUDED != _status[i]) {
      _status[i] = NOT_PROCESSED;
    }
  }

  _status[a] = PROCESSING_COMPLETE;

  _next_shell.resize_keep_storage(0);

  return;
}

void
ShellInfo::GrabNextShell(Set_of_Atoms& next_shell) {
  //  What used to be the next shall has been processed.
  _next_shell.set_vector(_status, PROCESSING_COMPLETE);

  next_shell.set_vector(_status, CURRENT_SHELL);

  _next_shell = std::move(next_shell);

  return;
}

//  #define MAGIC_FROM_SHELL_VARIABLES

#ifdef MAGIC_FROM_SHELL_VARIABLES
//  Fetch the value of numeric shell variable `vname` into `value`.
//  Useful for trying different combinations of numbers from outside the program.
template <typename T>
int
FromShellVariable(const char* vname, T& value) {
  const char* from_env = getenv(vname);
  if (from_env == nullptr) {
    return 0;
  }

  const const_IWSubstring s(from_env);
  return s.numeric_value(value);
}
#endif

ECFingerprint::ECFingerprint() {
  _min_radius = 0;
  _max_radius = 3;

  _bond_magic1 = 12541;
  _bond_magic2 = 14434;
  _bond_magic3 = 22690;
  _bond_magic4 = 1496;

  _magic1 = 12724;
  _magic2 = 20070;
  _magic3 = 3998;
  _magic4 = 30703;
  _magic5 = 3595;

#ifdef MAGIC_FROM_SHELL_VARIABLES
  FromShellVariable("MAGIC1", _magic1);
  FromShellVariable("MAGIC2", _magic2);
  FromShellVariable("MAGIC3", _magic3);
  FromShellVariable("MAGIC4", _magic4);
  FromShellVariable("MAGIC5", _magic5);
  FromShellVariable("BOND_MAGIC1", _bond_magic1);
  FromShellVariable("BOND_MAGIC2", _bond_magic2);
  FromShellVariable("BOND_MAGIC3", _bond_magic3);
  FromShellVariable("BOND_MAGIC4", _bond_magic4);
  cerr << _magic1 << ' ' << _magic2 << ' ' << _magic3 << ' ' << _magic4 << ' ' << _magic5
       << ' ' << _bond_magic1 << ' ' << _bond_magic2 << ' ' << _bond_magic3 << ' '
       << _bond_magic4 << '\n';
#endif

  _precise = true;

  _additive_shell_formation = true;
}

atom_type_t
ECFingerprint::_BondConstant(const Bond& b) const {
#ifdef DEBUG_EC_FINGERPRINT
  cerr << "_BondConstant from " << b.a1() << " to " << b.a2() << " aromatic "
       << b.is_aromatic() << " single " << b.is_single_bond() << '\n';
#endif
  if (b.is_aromatic()) {
    return _bond_magic4;
  }

  if (b.is_single_bond()) {
    return _bond_magic1;
  }

  if (b.is_double_bond()) {
    return _bond_magic2;
  }

  if (b.is_triple_bond()) {
    return _bond_magic3;
  }

  cerr << "EC_Fingerprint_:_BondConstant:what kind of bond ";
  b.debug_print(cerr);

  return 1;
}

void
ECFingerprint::_AddToRunningSum(ShellInfo& shell_info, const atom_number_t a1,
                                const atom_type_t bond_constant, const atom_number_t a2,
                                atom_type_t& running_sum) const {
#ifdef DEBUG_EC_FINGERPRINT
  cerr << "_AddToRunningSum from atom " << a1 << " type " << shell_info.atom_type(a1)
       << " to " << a2 << " type " << shell_info.atom_type(a2) << '\n';
#endif

  atom_type_t b;
  if (_precise) {
    b = shell_info.atom_type(a1) * _magic3 +
        _magic4 * bond_constant * (_magic5 + shell_info.atom_type(a2));
    //  cerr << "  _precise bond from " << a1 << " to " << a2 << " bond_constant " <<
    //  bond_constant
    //  << " bit " << b << '\n';
  } else {
    b = _magic4 * bond_constant + _magic5 * shell_info.atom_type(a2);
  }

  if (_additive_shell_formation) {
    running_sum += b;
  } else {
    running_sum *= b;
  }

#ifdef DEBUG_EC_FINGERPRINT
  cerr << "  bond from " << a1 << " to " << a2 << " bond_constant " << bond_constant
       << " bit " << b << " sum " << running_sum << '\n';
#endif
}

int
ECBaseWithOutput::Open(IWString& fname) {
  assert(!_output.is_open());

  if (!_output.open(fname)) {
    cerr << "ECBaseWithOutput::Open:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

JobParameters::JobParameters() {
  produce_output = true;

  function_as_tdt_filter = false;

  write_fixed_width_fingerprint = 0;

  write_counted_sparse_fingerprint = true;
}

int
ProduceFingerprint::PrepareToProcess(Molecule& m) {
  _sfc.clear();

  return 1;
}

int
ProduceFingerprint::DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                                IWString_and_File_Descriptor& output) {
  if (!job_parameters.produce_output) {
    return 1;
  }

  if (job_parameters.fp_writer.IsWritingDescriptors()) {
    return job_parameters.fp_writer.WriteFingerprint(m.name(), _sfc, output);
  }

  if (!job_parameters.function_as_tdt_filter) {
    output << job_parameters.smiles_tag << m.smiles() << ">\n";
    output << job_parameters.identifier_tag << m.name() << ">\n";
  }

  job_parameters.fp_writer.WriteFingerprint(m.name(), _sfc, output);

  return 1;
}

AtomMapCoverage::AtomMapCoverage() {
  _matoms = 0;
  _coverage = nullptr;
}

AtomMapCoverage::~AtomMapCoverage() {
  if (nullptr != _coverage) {
    delete[] _coverage;
  }
}

int
AtomMapCoverage::PrepareToProcess(Molecule& m) {
  if (m.natoms() > _matoms) {
    if (nullptr != _coverage) {
      delete[] _coverage;
    }
    _matoms = m.natoms();
    _coverage = new int[_matoms];
  }

  std::fill_n(_coverage, m.natoms(), 0);

  return 1;
}

int
AtomMapCoverage::FingerprintingComplete(Molecule& m) {
  assert(nullptr != _coverage);

  const int matoms = m.natoms();

  int min_coverage = _coverage[0];
  int max_coverage = _coverage[0];
  m.set_atom_map_number(0, _coverage[0]);

  for (int i = 1; i < matoms; ++i) {
    const int c = _coverage[i];
    m.set_atom_map_number(i, c);

    if (c < min_coverage) {
      min_coverage = c;
    } else if (c > max_coverage) {
      max_coverage = c;
    }
  }

  IWString tmp = m.name();
  tmp << " min " << min_coverage << " max " << max_coverage << " range "
      << (max_coverage - min_coverage);
  m.set_name(tmp);

  return 1;
}

void
WriteAllBits::Bit(const ShellInfo& shell_info, const atom_type_t running_sum,
                  const int radius) {
  const atom_number_t a0 = shell_info.a0();

  Molecule mcopy(shell_info.m());
  mcopy.set_isotope(a0, radius);

  _output << mcopy.smiles() << ' ' << shell_info.name() << " bit " << running_sum
          << " atom type" << shell_info.atom_type(a0) << " radius " << radius << "\n";

  _output.write_if_buffer_holds_more_than(8192);

  return;
}

void
ECCheckCollisions::Bit(const ShellInfo& shell_info, const atom_type_t running_sum,
                       const int radius) {
  const atom_number_t a0 = shell_info.a0();

  //  Would be more efficient to use the Molecule inside shell_info, but this is safer.
  //  This slows things down a LOT, even if we never look in the hash.
  Molecule mcopy(shell_info.m());
  mcopy.set_atom_map_number(a0, radius + 1);  //  Because radius 0 will not get marked.

  std::tuple<IWString, atom_type_t, int> tmp{mcopy.smiles(), shell_info.atom_type(a0),
                                             radius};

  auto f = _bit_to_description.find(running_sum);
  if (f == _bit_to_description.end()) {
    _bit_to_description.emplace(running_sum, std::move(tmp));
    return;
  }

  //  If same centre atom type and same radius, good enough.
  if (std::get<1>(tmp) == std::get<1>(f->second) &&
      std::get<2>(tmp) == std::get<2>(f->second)) {
    return;
  }

  _collisions_found++;

  _output << "Collision on bit " << running_sum << "\n";
  _output << mcopy.smiles() << ' ' << shell_info.name() << " atom type "
          << shell_info.atom_type(a0) << " radius " << _output << "centre atom types "
          << std::get<1>(tmp) << ' ' << std::get<1>(f->second) << '\n';
  _output << "radii             " << std::get<2>(tmp) << ' ' << std::get<2>(f->second)
          << '\n';

  _output.write_if_buffer_holds_more_than(8192);

  return;
}

ECBuildPrecedent::ECBuildPrecedent() {
  _bit_collisions_avoided = 0;
}

ECBuildPrecedent::~ECBuildPrecedent() {
}

int
ECBuildPrecedent::Report(std::ostream& output) const {
  output << "ECBuildPrecedent:avoided " << _bit_collisions_avoided << " bit collisions\n";
  for (int i = 0; i < _collision_at_radius.number_elements(); ++i) {
    if (_collision_at_radius[i]) {
      output << _collision_at_radius[i] << " at radius " << i << "\n";
    }
  }
  return output.good();
}

void
ECBuildPrecedent::Bit(const ShellInfo& shell_info, const atom_type_t running_sum,
                      const int radius) {
  auto f = _bit_count.find(running_sum);

  const atom_number_t a0 = shell_info.a0();

  if (f == _bit_count.end())  //  Newly discovered bit.
  {
    //  cerr << "got " << running_sum << '\n';
    Molecule mcopy(shell_info.m());
    mcopy.set_atom_map_number(a0, radius + 1);
    BitCount tmp(mcopy.smiles(), shell_info.atom_type(a0), radius,
                 static_cast<count_type_t>(1));

    _bit_count.emplace(running_sum, std::move(tmp));
  } else if (f->second.center_atom_type != shell_info.atom_type(a0)) {
    _bit_collisions_avoided++;
    //  cerr << "Atom type bit collision avoided " << running_sum << ' ' <<
    //  f->second.center_atom_type << " vs " << shell_info.atom_type(a0) << '\n';
    _collision_at_radius[radius]++;
  } else if (f->second.radius != radius) {
    _bit_collisions_avoided++;
    _collision_at_radius[radius]++;
    //  cerr << "Radius bit collision avoided " << running_sum << ' ' << f->second.radius
    //  << " vs "
    //  << radius << '\n';
  } else {
    f->second.count++;
  }

  return;
}

int
ECBuildPrecedent::WritePrecedentData(const char sep, const JobParameters& job_parameters,
                                     IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (!output.open(fname)) {
    cerr << "ECBuildPrecedent::WritePrecedentData:cannot open '" << fname << "'\n";
    return 0;
  }

  output << "#Atom type: \"" << job_parameters.atom_type_string << "\"\n";
  output << "#Bit Count Smiles Atype Radius\n";

  for (const auto& iter : _bit_count) {
    //  cerr << "writ " << iter.first << '\n';
    output << iter.first << sep <<  //  The bit
        iter.second.count << sep << iter.second.smiles << sep
           << iter.second.center_atom_type << sep << iter.second.radius << "\n";
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

ECUsePrecedent::ECUsePrecedent(const int max_radius) : _max_radius(max_radius) {
  _count = new count_type_t[max_radius + 1];
  _atom = new atom_number_t[max_radius + 1];

  _missing_at_radius = new count_type_t[max_radius + 1];
  std::fill_n(_missing_at_radius, max_radius + 1, 0);

  return;
}

ECUsePrecedent::~ECUsePrecedent() {
  delete[] _count;
  delete[] _atom;
  delete[] _missing_at_radius;
}

int
ECUsePrecedent::PrepareToProcess(Molecule& m) {
  std::fill_n(_count, _max_radius + 1, std::numeric_limits<count_type_t>::max());
  std::fill_n(_atom, _max_radius + 1, INVALID_ATOM_NUMBER);

  return 1;
}

//  #define DEBUG_USE_PRECEDENT

void
ECUsePrecedent::Bit(const ShellInfo& shell_info, const atom_type_t running_sum,
                    const int radius) {
  const auto f = _precedent.find(running_sum);
  if (f == _precedent.end()) {
#ifdef DEBUG_USE_PRECEDENT
    cerr << "Bit:no match for " << running_sum << ", radius " << radius << '\n';
#endif
    _count[radius] = 0;
    _atom[radius] = shell_info.a0();
    _missing_at_radius[radius]++;
    return;
  }

#ifdef DEBUG_USE_PRECEDENT
  cerr << "Bit: we have info on " << running_sum << " radius " << radius << " db "
       << f->second.radius << '\n';
#endif
  if (f->second.radius != radius) {  //  collision
    return;
  }

  const count_type_t db_count = f->second.count;

#ifdef DEBUG_USE_PRECEDENT
  cerr << "No collision, count is " << db_count << '\n';
#endif

  if (db_count < _count[radius]) {
    _count[radius] = static_cast<int>(db_count);
    _atom[radius] = shell_info.a0();
  }

  return;
}

int
ECUsePrecedent::FingerprintingComplete(Molecule& m) {
  IWString result;
  result << ' ';
  result.resize(20);  //  guess.

  atom_number_t rarest_atom = INVALID_ATOM_NUMBER;
  count_type_t lowest_count = std::numeric_limits<count_type_t>::max();

  //  cerr << "ECUsePrecedent::FingerprintingComplete:radius " << _max_radius << '\n';

  count_type_t previous_count = std::numeric_limits<count_type_t>::max();

  for (int r = 0; r <= _max_radius; ++r) {
    if (r > 0) {
      result << ',';
    }

#ifdef DEBUG_USE_PRECEDENT
    cerr << "Radius " << r << " count " << _count[r] << " ationm " << _atom[r] << '\n';
#endif
    result << r << ',';
    if (std::numeric_limits<count_type_t>::max() == _count[r]) {
      result << '0';
    } else {
      if (_count[r] < lowest_count) {
        lowest_count = _count[r];
        rarest_atom = _atom[r];
      }
      result << _count[r];
      if (_count[r] > previous_count) {
#ifdef NOTIFY_COUNT_INCREASE
        cerr << "ECUsePrecedent::FingerprintingComplete:unexpected count increase, rad "
             << (r - 1) << " count " << previous_count << " radius " << r << " count "
             << _count[r] << ' ' << m.name() << '\n';
#endif
      }
      previous_count = _count[r];
    }
  }

#ifdef DEBUG_USE_PRECEDENT
  cerr << "lowest_count " << lowest_count << " rarest_atom " << rarest_atom << '\n';
#endif

  if (rarest_atom != INVALID_ATOM_NUMBER) {
    m.set_atom_map_number(rarest_atom, lowest_count + 1);
  } else {
    cerr << "No rarest atom in " << m.smiles() << '\n';
  }

  m.append_to_name(result);

  return 1;
}

int
ECUsePrecedent::DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                            IWString_and_File_Descriptor& output) {
  if (!job_parameters.produce_output)  //  No fingerprint output, must output smiles
  {
    output << m.smiles() << ' ' << m.name() << '\n';
  } else {
    output << job_parameters.smiles_tag << m.smiles() << ">\n";
    output << job_parameters.identifier_tag << m.name() << ">\n";
  }

  output.write_if_buffer_holds_more_than(8192);

  return output.good();
}

int
ECUsePrecedent::Report(std::ostream& output) const {
  for (int i = 0; i <= _max_radius; ++i) {
    output << " radius " << i << " missing " << _missing_at_radius[i] << "\n";
  }

  return output.good();
}

int
ECUsePrecedent::ReadPrecedentData(IWString& fname) {
  cerr << "ReadPrecedentData from '" << fname << "'\n";
  iwstring_data_source input;
  if (!input.open(fname)) {
    cerr << "ECFingerprint::ReadPrecedentData:cannot open '" << fname << "'\n";
    return 0;
  }

  const_IWSubstring line;
  while (input.next_record(line)) {
    //  cerr << "Examining '" << line << "'\n";

    if (line.starts_with('#')) {
      continue;
    }

    if (!_ParsePrecedentRecord(line)) {
      cerr << "ECFingerprint::ReadPrecedentData:invalid input '" << line << "'\n";
      return 0;
    }
  }

  return _precedent.size() > 0;
}

int
ECUsePrecedent::_ParsePrecedentRecord(const const_IWSubstring& line) {
  int i = 0;
  const_IWSubstring token;

  if (!line.nextword(token, i)) {
    cerr << "ECFingerprint::_ParsePrecedentRecord:cannot extract first token\n";
    return 0;
  }

  atom_type_t b;
  if (!token.numeric_value(b)) {
    cerr << "efficient::_ParsePrecedentRecord:invalid bit number\n";
    return 0;
  }

  if (!line.nextword(token, i)) {
    cerr << "ECFingerprint::_ParsePrecedentRecord:cannot extract second token\n";
    return 0;
  }

  count_type_t count;
  if (!token.numeric_value(count)) {
    cerr << "ECFingerprint::_ParsePrecedentRecord:invalid count\n";
    return 0;
  }

  //  Skip over smiles and central atom type to get radius token
  if (!line.nextword(token, i) || !line.nextword(token, i) || !line.nextword(token, i)) {
    cerr << "ECFingerprint::_ParsePrecedentRecord:cannot extract tokens 3 and 4\n";
    return 0;
  }

  int radius;
  if (!token.numeric_value(radius)) {
    cerr << "ECFingerprint::_ParsePrecedentRecord:Invalid radius\n";
    return 0;
  }

  //  std::tuple<count_type_t, int> tmp{count, radius};
  BitCount tmp(count, radius);

  _precedent.emplace(b, std::move(tmp));

  return 1;
}

int
hash_from_file(iwstring_data_source& input,
               std::unordered_set<atom_type_t>& destination) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with('#') || buffer.empty()) {
      continue;
    }
    buffer.truncate_at_first(' ');
    atom_type_t b;
    if (!buffer.numeric_value(b)) {
      cerr << "hash_from_file::_ReadBitsToFind:invalid bit " << buffer << "'\n";
      return 0;
    }
    destination.insert(b);
  }

  return destination.size();
}

ECBitMeanings::ECBitMeanings() {
  _bits_found = 0;
  _bits_not_found = 0;
}

ECBitMeanings::~ECBitMeanings() {
}

int
ECBitMeanings::ReadBitsToFind(IWString& fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "ECBitMeanings::ReadBitsToFind:cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_strip_leading_blanks(1);
  input.set_strip_trailing_blanks(1);

  return hash_from_file(input, _bits_to_find);
}

int
ECBitMeanings::PrepareToProcess(Molecule& m) {
  _buffer_current_molecule.resize_keep_storage(0);
  return 1;
}

void
ECBitMeanings::Bit(const ShellInfo& shell_info, const atom_type_t running_sum,
                   const int radius) {
  if (const auto iter = _bits_to_find.find(running_sum); iter == _bits_to_find.end()) {
    _bits_not_found++;
    return;
  }

  constexpr char sep = ' ';

  _bits_found++;

  Molecule mcopy(shell_info.m());
  mcopy.set_atom_map_number(shell_info.a0(), radius);
  _buffer_current_molecule << mcopy.smiles() << sep << shell_info._m.name() << sep
                           << running_sum << '\n';
}

int
ECBitMeanings::DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                           IWString_and_File_Descriptor& output) {
  if (_buffer_current_molecule.empty()) {
    return 1;
  }

  output << _buffer_current_molecule;
  output.write_if_buffer_holds_more_than(8192);
  return 1;
}

ECFilterByBits::ECFilterByBits() {
  _write_current_molecule = true;
}

int
ECFilterByBits::ReadBitsToFilter(IWString& fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "ECFilterByBits:ReadBitsToFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return hash_from_file(input, _bits_to_find);
}

void
ECFilterByBits::Bit(const ShellInfo& shell_info, const atom_type_t running_sum,
                    const int radius) {
  //  If already discarded, do not check again.
  if (!_write_current_molecule) {
    return;
  }
  if (const auto iter = _bits_to_find.find(running_sum); iter != _bits_to_find.end()) {
    _write_current_molecule = false;
  }
}

int
ECFilterByBits::DoAnyOutput(Molecule& m, JobParameters& job_parameters,
                            IWString_and_File_Descriptor& output) {
  if (!_write_current_molecule) {
    return 0;
  }

  output << m.smiles() << ' ' << m.name() << '\n';
  output.write_if_buffer_holds_more_than(8192);
  return 1;
}

}  //  namespace ec_fingerprint
