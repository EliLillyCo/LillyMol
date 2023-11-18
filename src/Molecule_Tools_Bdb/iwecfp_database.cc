#include "iwecfp_database.h"

#include <stdlib.h>

#include <algorithm>
#include <limits>
#include <memory>

#include "db_cxx.h"

namespace iwecfp_database {

using std::cerr;

// Data is packed into a binary form in the most efficient way possible.
// The Pubchem version of this database is 3.4GB, corporate db 134MB.
void
form_key(uint64_t b, int radius, unsigned int atom_constant_centre_atom, DBKey& dbkey)
{
  // cerr << "RAD " << radius << " Acca " << atom_constant_centre_atom << '\n';

  static constexpr uint64_t maxu32 = static_cast<uint64_t>(std::numeric_limits<uint32_t>::max());

  if (b > maxu32) {
    dbkey._bit = b % maxu32;
  } else {
    dbkey._bit = b;
  }
  dbkey._acca = static_cast<unsigned int>(atom_constant_centre_atom);
  dbkey._radius = static_cast<unsigned char>(radius);
  dbkey._nu1 = static_cast<unsigned char>(0);
  dbkey._nu2 = static_cast<unsigned char>(0);
  dbkey._nu3 = static_cast<unsigned char>(0);

  return;
}

/*std::ostream &
operator << (const DBKey & k, std::ostream & os)
{
  os << " b " << k._bit << " acca " << k._acca << " rad " << k._radius;

  return os;
}*/

std::ostream&
operator<<(std::ostream& os, const DBKey& k)
{
  os << " b " << k._bit << " acca " << k._acca << " rad " << static_cast<int>(k._radius);

  return os;
}

int
debug_print_key_components(const DBKey& dbkey, std::ostream& os)
{
  os << "bit " << dbkey._bit << " R " << static_cast<int>(dbkey._radius) << " ACCA "
     << dbkey._acca << '\n';

  return 1;
}

int
operator==(const DBKey& k1, const DBKey& k2)
{
  if (k1._bit != k2._bit) {
    return 0;
  }

  if (k1._acca != k2._acca) {
    return 0;
  }

  return k1._radius == k2._radius;
}

/*
  We just return the first size_t of bytes - we are dealing with 64 bit
  items, so we should be OK
*/

size_t
IWdbkeyHash::operator()(DBKey const& s) const
{
  return s._bit * s._acca + s._radius;
}

Bit_Produced::Bit_Produced(unsigned int b, unsigned int r, int c, int acc)
    : _bit(b), _radius(r), _centre_atom(c), _atom_constant_centre_atom(acc)
{
  _count = 0;

  return;
}

Fingerprint_Characteristics::Fingerprint_Characteristics()
{
  _min_shell_radius = 0;
  _max_shell_radius = std::numeric_limits<int>::max();
  _additive = 1;
  _atype = (IWATTYPE_Z | DIFFERENTIATE_RINGS);

  _intra_molecular_bit_collisions = 0;
}

int
Fingerprint_Characteristics::build(const Command_Line& cl,
                                   const int mr,  // max radius from database
                                   const int verbose)
{
  if (cl.option_present('r')) {
    if (!cl.value('r', _min_shell_radius) || _min_shell_radius < 0) {
      cerr << "The min shell radius (-r) must be a whole +ve number\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will only fingerprint paths larger than " << _min_shell_radius
           << " bonds\n";
    }
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', _max_shell_radius) || _max_shell_radius < 0) {
      cerr << "The max shell radius (-R) must be a whole +ve number\n";
      return 0;
    }

    if (verbose) {
      cerr << "Max radius " << _max_shell_radius << '\n';
    }
  } else if (mr > 0) {
    _max_shell_radius = mr;
    if (verbose) {
      cerr << "Fingerprint_Characteristics::build:using database default radius "
           << _max_shell_radius << '\n';
    }
  } else {
    cerr << "Must specify maximum radius for synthetic example fragments\n";
    return 0;
  }

  if (_max_shell_radius < _min_shell_radius) {
    cerr << "Inconsistent min " << _min_shell_radius << " and max " << _max_shell_radius
         << " shell radius\n";
    return 0;
  }

  if (cl.option_present('m')) {
    _additive = 0;

    if (verbose) {
      cerr << "Fingerprints formed with multiplication operations\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');

    IWString tag("NC");  // not used
    _atype = determine_atom_type_and_set_tag(p, tag, verbose);

    if (_atype <= 0) {
      cerr << "Fingerprint_Characteristics::build:cannot determine atom typying to use '"
           << p << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Fingerprint_Characteristics::string_atom_type(IWString& s) const
{
  if (iwattype_convert_to_string_form(_atype, s)) {
    return 1;
  }

  cerr << "Fingerprint_Characteristics::string_atype:Unrecognised atom type !!! "
       << _atype << '\n';
  return 0;
}

void
Fingerprint_Characteristics::Increment(unsigned int processing_status,
                unsigned int bond_constant,
                unsigned int atom_constant,
                uint64_t& sum_so_far) const {
  if (_additive) {
    sum_so_far += (processing_status % 73) * bond_constant * atom_constant;
  } else {
    sum_so_far *= (processing_status % 73) * bond_constant * atom_constant;
  }
}

}  // namespace iwecfp_database
