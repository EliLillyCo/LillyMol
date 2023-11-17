#ifndef MOLECULE_TOOLS_IWECFP_DATABASE_H
#define MOLECULE_TOOLS_IWECFP_DATABASE_H

#include <iostream>
#include <unordered_map>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/atom_typing.h"

namespace iwecfp_database {

/*
  They key is a 10 byte entity.
  First four bytes are the bit number.
  Second four bytes are the atom constant of centre atom
  Then follows a byte is the radius
  The last byte is not used
*/

struct DBKey
{
  unsigned int   _bit;    // bit number
  unsigned int   _acca;   // atom constant centre atom
  unsigned char  _radius; // shell radius

// Seems compiler pads things to 12 bytes

  unsigned char  _nu1;
  unsigned char  _nu2;
  unsigned char  _nu3;
};

//extern std::ostream & operator << (const DBKey & k, std::ostream & os);
extern std::ostream & operator << (std::ostream & os, const DBKey & k);


typedef struct DBKey DBKey;

class IWdbkeyHash
{
  private:
  public:

#if defined (IW_INTEL_COMPILER)

    static const size_t bucket_size = 4;
    static const size_t min_buckets = 8;
    bool  operator () (const DBKey &, const DBKey &) const;

#endif

    size_t operator () (const DBKey &) const;
};

extern int operator == (const DBKey &, const DBKey &);

extern void form_key(uint64_t b,
          int radius,
          unsigned int atom_constant_centre_atom,
          DBKey & dkey);

extern int debug_print_key_components (const DBKey & dbkey, std::ostream & os);

struct Count_Radius
{
  int _count;
  int _radius;
};

typedef struct Count_Radius Count_Radius;

#define ATYPE_KEY "_ATYPE"
#define RADIUS_KEY "_RADIUS"


class Bit_Produced
{
  protected:
    const unsigned int _bit;
    const unsigned int _radius;

    const atom_number_t _centre_atom;

//  To help keep bit collisions down, we keep track of the
//  atom constant of the centre atom

    const unsigned int _atom_constant_centre_atom;

    unsigned int _count;

  public:
    Bit_Produced (unsigned int, unsigned int, atom_number_t, int);

    unsigned int bit () const { return _bit;}

    void extra() { _count++;}

    int radius () const { return _radius;}

    int count () const { return _count;}

    atom_number_t centre_atom() const { return _centre_atom;}

    unsigned int atom_constant_centre_atom () const { return _atom_constant_centre_atom;}
};

//typedef unordered_map<DBKey, Bit_Produced *, IWdbkeyHash> Set_of_Bits;
template <typename T>
class Set_of_Bits : public std::unordered_map<DBKey, T *, IWdbkeyHash>
{
  private:

// private functions

    void _delete_items();

  public:
    ~Set_of_Bits();

    void clear ();

    template <typename O> int debug_print(O & output) const;
};

template <typename T>
Set_of_Bits<T>::~Set_of_Bits()
{
//cerr << "Set_of_Bits::destructor, destroying " << this->size() << " items\n";
  for (auto i : *this)
  {
    delete i.second;
  }

  return;
}

template<typename T> template <typename O>
int
Set_of_Bits<T>::debug_print(O & output) const
{
  output << "Set_of_Bits::debug_print: " << this->size() << " bits\n";
  for (auto i : *this)
  {
    output << "b " << i.first._bit << " acca " << i.first._acca << " r " << static_cast<int>(i.first._radius) << '\n';
  }
  return 1;
}

template <typename T>
void
Set_of_Bits<T>::_delete_items()
{
  for (auto i : *this)
  {
    delete i.second;
  }

  return;
}

template <typename T>
void
Set_of_Bits<T>::clear()
{
  _delete_items();
  std::unordered_map<DBKey, T*, IWdbkeyHash>::clear();

  return;
}

class Fingerprint_Characteristics
{
  private:
    int _min_shell_radius;
    int _max_shell_radius;
    int _additive;
    int _atype;
    int _intra_molecular_bit_collisions;

//  private functions

  public:
    Fingerprint_Characteristics();

    void usage (const int rc) const;    // only exits if given non zero rc

    void set_min_shell_radius (int s) { _min_shell_radius = s;}
    void set_max_shell_radius (int s) { _max_shell_radius = s;}
    void set_additive (int s) { _additive = s;}
    void set_atype (int s) { _atype = s;}

    void set_only_set_bits_for_max_radius_shell() { _min_shell_radius = _max_shell_radius;}

    int min_shell_radius() const { return _min_shell_radius;}
    int max_shell_radius() const { return _max_shell_radius;}
    int additive() const { return _additive;}
    int atype() const { return _atype;}
    
    int intra_molecular_bit_collisions () const { return _intra_molecular_bit_collisions;}

    void intra_molecular_bit_collision() { _intra_molecular_bit_collisions++;}

    int string_atom_type (IWString & s) const;

    int build (const Command_Line & cl, const int mr, const int verbose);

    // update `sum_so_far` for a new member in the next shell.
    void Increment(unsigned int processing_status,
                unsigned int bond_constant,
                unsigned int atom_constant,
                uint64_t& sum_so_far) const;
};

template <typename T>
int
handle_newly_found_bit (atom_number_t centre_atom,
                        const unsigned int atom_constant_centre_atom,
                        const int radius,
                        const uint64_t b,
                        Set_of_Bits<T> &  sob,
                        Fingerprint_Characteristics & fc)
{
#ifdef DEBUG_BITS_FORMED
  std::cerr << "From atom " << centre_atom << " radius " << radius << " bit " << b << '\n';
#endif

  DBKey dbkey;

  form_key(b, radius, atom_constant_centre_atom, dbkey);

  auto f = sob.find(dbkey);

  if (f != sob.end()) {  // already encountered this bit
    T * bp = (*f).second;
    if (atom_constant_centre_atom != bp->atom_constant_centre_atom()) {
      fc.intra_molecular_bit_collision();
      return 0;
    }

    bp->extra();
    return 1;
  }

  T * bp = new T(b, radius, centre_atom, atom_constant_centre_atom);

  sob[dbkey] = bp;

  return 1;
}

#ifdef NOW_MEMBER_FUNCTION
static void
increment (uint64_t & sum_so_far,
           unsigned int processing_status,
           unsigned int bc,
           unsigned int atom_constant,
           const Fingerprint_Characteristics & fc)
{
  if (fc.additive())
    sum_so_far += (processing_status % 73) * bc * atom_constant;
  else
    sum_so_far *= (processing_status % 73) * bc * atom_constant;

  return;
}
#endif

static unsigned int 
bond_constant(const Bond * bondi)
{
  if (bondi->is_aromatic())
    return 11;
  if (bondi->is_triple_bond())
    return 7;
  if (bondi->is_double_bond())
    return 5;
  
  return 3;
}

#define PROCESSING_FINISHED 0

// This is an arbitrary number, but guaranteed that no 32 bit number, negative,
// will collide.
static constexpr int64_t kNotYetSeen = std::numeric_limits<int64_t>::min();

static constexpr int64_t kNextShell = std::numeric_limits<int64_t>::min() + 1;
static constexpr int64_t kWillBeNextShell = std::numeric_limits<int64_t>::min() + 2;
static constexpr int64_t kProcessingFinished = std::numeric_limits<int64_t>::min() + 3;

#ifdef OLD_VERSION_THAT_DOES_NOT_WORK
template <typename T>
int
generate_shells (int matoms,
                int radius,
                const Atom * const * a,
                atom_number_t centre_atom,
                const unsigned int * atom_constant,
                int64_t* processing_status,
                uint64_t sum_so_far,
                Set_of_Bits<T> & sob,
                Fingerprint_Characteristics & fc)
{
  radius++;

  if (fc.additive())
    sum_so_far *= 7879;   // an arbitrary prime number

//#define DEBUG_ECFP_BITS 1
#ifdef DEBUG_ECFP_BITS
  std::cerr << "On entry sum_so_far " << sum_so_far << '\n';
#endif

  int continue_processing = 0;

  for (int i = 0; i < matoms; i++) {
    if (processing_status[i] <= 0) {
      continue;
    }

    for (const Bond* b : *a[i]) {
      atom_number_t k = b->other(i);
      // std::cerr << "From atom " << i << " to " << k << '\n';

      // expanding the shell, either outward or closing a ring
      if (PROCESSING_FINISHED == processing_status[k]) {
        unsigned int bc = bond_constant(b);
        // std::cerr << "B4 inc " << sum_so_far << " bc " << bc << '\n';
        fc.Increment(processing_status[i], bc, atom_constant[i], sum_so_far);
        // std::cerr << "from " << centre_atom << " at rad " << radius << " atom " << i << " sum " << sum_so_far << " PC " << processing_status[i] << " AC " << atom_constant[i] << '\n';
      }
      else if (processing_status[k] > 0) {  // will be processed later in this loop
      } else if (kNotYetSeen == processing_status[k]) {   // mark for next time
        processing_status[k] = - static_cast<int64_t>(atom_constant[i]);
        continue_processing = 1;
      } else {
        // std::cerr << radius << "else processing_status " << processing_status[k] << '\n';
//      processing_status[k] -= atom_constant[i];
//      continue_processing = 1;
      }
    }
  }

  if (radius >= fc.min_shell_radius()) {
#ifdef DEBUG_ECFP_BITS
    std::cerr << "From " << centre_atom << " radius " << radius << " hit bit " << sum_so_far << '\n';
#endif
    handle_newly_found_bit(centre_atom, atom_constant[centre_atom], radius, sum_so_far, sob, fc);
  }

  if (! continue_processing) {
    return 1;
  }

  if (radius >= fc.max_shell_radius()) {
    return 1;
  }

#ifdef DEBUG_ECFP_BITS
  std::cerr << "from " << centre_atom << " at radius " << radius << " atoms";
#endif

  for (int i = 0; i < matoms; i++) {
    if (processing_status[i] > 0)
      processing_status[i] = PROCESSING_FINISHED;
    else if (kNotYetSeen == processing_status[i])
      ;
    else if (processing_status[i] < 0)
    {
      processing_status[i] = - processing_status[i];
#ifdef DEBUG_ECFP_BITS
      std::cerr << ' ' << i;
#endif
    }
  }

#ifdef DEBUG_ECFP_BITS
  std::cerr << '\n';
#endif

  return generate_shells(matoms, radius, a, centre_atom, atom_constant, processing_status, sum_so_far, sob, fc);
}

template <typename T>
int
compute_fingerprints (Molecule & m,
                      Fingerprint_Characteristics & fc,
                      const unsigned int * atom_constant,
                      Set_of_Bits<T> & sob)
{
  sob.clear();

#ifdef DEBUG_ECFP_BITS
  std::cerr << "Computing fingerprints on '" << m.unique_smiles() << " radius " << fc.min_shell_radius() << " to " << fc.max_shell_radius() << '\n';

  for (int i = 0; i < m.natoms(); i++)
  {
    std::cerr << "Atom " << i << " type " << atom_constant[i] << " " << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

//cerr << "Computing fingerprints type " << iwecfp_atom_type << '\n';

  const int matoms = m.natoms();

  Atom ** atoms = new Atom * [matoms]; std::unique_ptr<Atom *[]> free_atoms(atoms);

  m.atoms((const Atom **) atoms);   // disregard of const OK

  // Use signed int64 for this, so we can store the uint32_2 atom types,
  // as well as their negative values.
  std::unique_ptr<int64_t[]> processing_status = std::make_unique<int64_t[]>(matoms);

  for (int i = 0; i < matoms; i++) {
#ifdef DEBUG_ECFP_BITS
    if (0 == fc.min_shell_radius())
      std::cerr << "Starting with atom " << i << " bit " << atom_constant[i] << '\n';
#endif
    if (0 == fc.min_shell_radius()) {
      handle_newly_found_bit(i, atom_constant[i], 0, atom_constant[i], sob, fc);
      if (0 == fc.max_shell_radius())
        continue;
    }

    std::fill_n(processing_status.get(), matoms, kNotYetSeen);

    processing_status[i] = PROCESSING_FINISHED;

    const Atom * ai = atoms[i];
    for (const Bond* b : *ai) {
      atom_number_t k = b->other(i);

      processing_status[k] = atom_constant[i];
    }

    generate_shells(matoms, 0, atoms, i, atom_constant, processing_status.get(), 
                    static_cast<uint64_t>(atom_constant[i]), sob, fc);
  }

#ifdef DEBUG_ECFP_BITS
  std::cerr << "Fingerprint contqainsq " << sob.size() << " bits\n";
  sob.debug_print(std::cerr);
#endif

  return 1;
}

template <typename T>
int
compute_fingerprints(Molecule & m,
                     Fingerprint_Characteristics & fc,
                     Set_of_Bits<T> & sob)
{
  unsigned int * atype = new unsigned int[m.natoms()]; std::unique_ptr<unsigned int[]> free_atype(atype);

  if (! assign_atom_types(m, fc.atype(), atype)) {
    std::cerr << "compute_fingerprints::cannot assign atom types for " << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  return compute_fingerprints(m, fc, atype, sob);
}
#endif  // OLD_VERSION_THAT_DOES_NOT_WORK

template <typename T>
int
generate_shells(int matoms,
                int radius,
                const Atom * const * a,
                atom_number_t centre_atom,
                const unsigned int * atom_constant,
                int64_t* processing_status,
                uint64_t sum_so_far,
                Set_of_Bits<T> & sob,
                Fingerprint_Characteristics & fc)
{
  radius++;

  if (fc.additive())
    sum_so_far *= 7879;   // an arbitrary prime number

// #define DEBUG_ECFP_BITS 1
#ifdef DEBUG_ECFP_BITS
  std::cerr << "On entry sum_so_far " << sum_so_far << '\n';
#endif

  int continue_processing = 0;

  for (int i = 0; i < matoms; i++) {
    if (processing_status[i] != kNextShell) {
      continue;
    }

    for (const Bond* b : *a[i]) {
      atom_number_t k = b->other(i);
      // std::cerr << "From atom " << i << " to " << k << '\n';

      // Either going back to previous atom or closing a ring
      if (kProcessingFinished == processing_status[k]) {
        unsigned int bc = bond_constant(b);
        // std::cerr << "B4 inc " << sum_so_far << " bc " << bc << '\n';
        fc.Increment(atom_constant[i], bc, atom_constant[k], sum_so_far);
        // std::cerr << "from " << centre_atom << " at rad " << radius << " atom " << i << " sum " << sum_so_far << " PC " << processing_status[i] << " AC " << atom_constant[i] << '\n';
      } else if (processing_status[k] == kNotYetSeen) {
        processing_status[k] = kWillBeNextShell;
        continue_processing = 1;
      }
    }
  }

  if (radius >= fc.min_shell_radius()) {
#ifdef DEBUG_ECFP_BITS
    std::cerr << "From " << centre_atom << " radius " << radius << " hit bit " << sum_so_far << '\n';
#endif
    handle_newly_found_bit(centre_atom, atom_constant[centre_atom], radius, sum_so_far, sob, fc);
  }

  if (! continue_processing) {
    return 1;
  }

  if (radius >= fc.max_shell_radius()) {
    return 1;
  }

#ifdef DEBUG_ECFP_BITS
  std::cerr << "from " << centre_atom << " at radius " << radius << " atoms";
#endif

  for (int i = 0; i < matoms; i++) {
    if (processing_status[i] == kNextShell) {
      processing_status[i] = kProcessingFinished;
    } else if (processing_status[i] == kWillBeNextShell) {
      processing_status[i] = kNextShell;
    }
  }

#ifdef DEBUG_ECFP_BITS
  std::cerr << '\n';
#endif

  return generate_shells(matoms, radius, a, centre_atom, atom_constant, processing_status, sum_so_far, sob, fc);
}

template <typename T>
int
compute_fingerprints(Molecule & m,
                      Fingerprint_Characteristics & fc,
                      const unsigned int * atom_constant,
                      Set_of_Bits<T> & sob)
{
  sob.clear();

#ifdef DEBUG_ECFP_BITS
  std::cerr << "Computing fingerprints on '" << m.unique_smiles() << " radius " << fc.min_shell_radius() << " to " << fc.max_shell_radius() << '\n';

  for (int i = 0; i < m.natoms(); i++)
  {
    std::cerr << "Atom " << i << " type " << atom_constant[i] << " " << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

//cerr << "Computing fingerprints type " << iwecfp_atom_type << '\n';

  const int matoms = m.natoms();

  Atom ** atoms = new Atom * [matoms]; std::unique_ptr<Atom *[]> free_atoms(atoms);

  m.atoms((const Atom **) atoms);   // disregard of const OK

  // Use signed int64 for this, so we can store the uint32_2 atom types,
  // as well as their negative values.
  std::unique_ptr<int64_t[]> processing_status = std::make_unique<int64_t[]>(matoms);

  for (int i = 0; i < matoms; i++) {
#ifdef DEBUG_ECFP_BITS
    if (0 == fc.min_shell_radius())
      std::cerr << "Starting with atom " << i << " bit " << atom_constant[i] << '\n';
#endif
    if (0 == fc.min_shell_radius()) {
      handle_newly_found_bit(i, atom_constant[i], 0, atom_constant[i], sob, fc);
      if (0 == fc.max_shell_radius())
        continue;
    }

    std::fill_n(processing_status.get(), matoms, kNotYetSeen);

    processing_status[i] = kProcessingFinished;

    const Atom * ai = atoms[i];
    for (const Bond* b : *ai) {
      atom_number_t k = b->other(i);

      processing_status[k] = kNextShell;
    }

    generate_shells(matoms, 0, atoms, i, atom_constant, processing_status.get(), 
                    static_cast<uint64_t>(atom_constant[i]), sob, fc);
  }

#ifdef DEBUG_ECFP_BITS
  std::cerr << "Fingerprint contqainsq " << sob.size() << " bits\n";
  sob.debug_print(std::cerr);
#endif

  return 1;
}


template <typename T>
int
compute_fingerprints(Molecule & m,
                     Fingerprint_Characteristics & fc,
                     Set_of_Bits<T> & sob)
{
  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(m.natoms());

  if (! assign_atom_types(m, fc.atype(), atype.get())) {
    std::cerr << "compute_fingerprints::cannot assign atom types for " << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  return compute_fingerprints(m, fc, atype.get(), sob);
}

}  // namespace iwecfp_database
#endif  // MOLECULE_TOOLS_IWECFP_DATABASE_H
