/*
  Randomly change molecules
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "google/protobuf/descriptor.h"
#include "google/protobuf/message.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Utilities/GFP_Tools/gfp_standard.h"
#include "Utilities/GFP_Tools/sparse_collection.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/numass.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/random_molecular_permutations.pb.h"
#include "Molecule_Tools/set_of_target_molecules.h"

using std::cerr;

static int verbose = 0;

static int molecules_read = 0;
static int molecules_written = 0;

static int molecules_made = 0;

static Accumulator_Int<int> atoms_in_starting_molecules;
static Accumulator_Int<int> atoms_in_molecules_written;

static int write_starting_molecule_once = 0;
static int write_starting_molecule_before_each_copy = 0;
static int write_starting_molecule_before_every_copy = 0;

static int allow_fused_rings_to_swap_atoms = 0;

static double probability_add_double_bond = 0.5;
static double probability_remove_atom = 1.0;
static double probability_lower_bond_order = 0.9;
static double probability_break_aliphatic_ring = 0.4;
static double probability_make_ring = 0.5;
static double probability_add_carbon = 0.2;
static double probability_move_fragment = 1.0;
static double probability_swap_adjacent_atoms = 1.0;
static double probability_change_carbon_to_nitrogen = 0.5;
static double probability_change_nitrogen_to_carbon = 0.6;
static double probability_change_carbon_to_oxygen = 0.4;
static double probability_change_oxygen_to_something = 0.6;
static double probability_add_from_single_attachment_library = 0.6;
static double probability_add_from_aromatic_attachment_library = 0.4;
static double probability_add_from_double_attachment_library = 0.5;
static double probability_change_ring_substitution = 1.0;
static double probability_remove_terminal_functional_group = 0.8;
static double probability_remove_embedded_functional_group = 0.7;
static double probability_remove_terminal_ring = 0.1;
static double probability_switch_aromatic_ring = 0.1;
static double probability_switch_aliphatic_ring = 0.1;
static double probability_breed = 1.0;
static double probability_breed_close_to_target = 1.0;
static double probability_split_fused_ring = 0.03;
static double probability_create_fused_ring = 0.15;
static double probability_expand_ring = 0.07;
static double probability_shrink_ring = 0.06;

static double probability_inter_permutation=0.5;

static rmp::Config config;
using rmp::Config;

static Number_Assigner number_assigner;

static int can_form_strongly_fused_ring_systems = 0;

static int lower_atom_count_cutoff = 0;
static int upper_atom_count_cutoff = std::numeric_limits<int>::max();

static std::random_device rd;

/*
  Various ways of controlling the execution
*/

static IWString stop_file;

static int run_for = 0;

static time_t tzero;    // not initialised

static int stop_after_producing = std::numeric_limits<int>::max();

static int warn_invalid_valence = 1;

/*
  Destroying an aromatic ring is done within lower_bond_order
*/

static double probability_destroy_aromatic_ring = 0.2;

static int append_permutation_number_to_name = 0;

class RandomState {
  private:
    std::default_random_engine _rng;
    std::uniform_real_distribution<float> _fraction;

    // By having a number of pre-existing uniform_int_distribution
    // we can easily generate random integers in those ranges.
    std::vector<std::uniform_int_distribution<int>> _uniform;

  public:
    RandomState();

    float operator()() { 
      return _fraction(_rng);
    }

    float Fraction() {
      return _fraction(_rng);
    }

    template <typename T>
    int RandomInt(std::uniform_int_distribution<T>& dist) {
      return dist(_rng);
    }

    int OkProbability(const float threshold) {
      return threshold > _fraction(_rng);
    }

    // Return a random number in the range [0, max).
    int RandomInt(int max);

    // Return a random number in the range [min, max).
    int RandomInt(int min, int max);

    template <typename T>
    T * RandomItem(const resizable_array_p<T>& items);
    template <typename T>
    T RandomItem(const resizable_array<T>& items);
};

RandomState::RandomState() : _fraction(0.0f, 1.0f) {
  _uniform.reserve(50);
  for (int i = 0; i < 50; ++i) {
    _uniform.emplace_back(std::uniform_int_distribution<int>(0, i));
  }
}

int
RandomState::RandomInt(int max) {
  if (max < static_cast<int>(_uniform.size())) {
    return _uniform[max](_rng);
  }

  std::uniform_int_distribution<int> u(0, max);
  return u(_rng);
}

int
RandomState::RandomInt(int min, int max) {
  std::uniform_int_distribution<int> u(min, max);
  return u(_rng);
}

template <typename T>
T*
RandomState::RandomItem(const resizable_array_p<T>& items) {
  if (items.empty()) {
    return nullptr;
  }
  if (items.size() == 1) {
    return items[0];
  }
  int ndx = RandomInt(items.size() - 1);
  return items[ndx];
}

template <typename T>
T
RandomState::RandomItem(const resizable_array<T>& items) {
  if (items.empty()) {
    return T{};
  }
  if (items.size() == 1) {
    return items[0];
  }
  int ndx = RandomInt(items.size() - 1);
  return items[ndx];
}

// THis whole concept is broken.
// The idea was to identify some atoms that cannot change, and store them
// in this class.
// but as the molecule changes, atoms are removed and added, this array will
// quickly become wrong, either indices wrong, or too short for the actual
// molecule.
class Can_Change
{
  private:
    const int _matoms;
    int * _can_change;

  public:
    Can_Change(const int matoms);
    ~Can_Change();

    int operator[](const int s) const;

    const int * rawdata() const { return _can_change;}
    int * rawdata() { return _can_change;}

    int matoms() const { return _matoms;}

    int AnyMembersSet(const Set_of_Atoms& atom_numbers, int value) const;
};

Can_Change::Can_Change(const int m) : _matoms(m)
{
  _can_change = new int[m];
  std::fill_n(_can_change, m, 1);

  return;
}

Can_Change::~Can_Change()
{
  delete [] _can_change;
}

int
Can_Change::operator[](const int s) const
{
  if (s < _matoms)
    return _can_change[s];

  return 0;
}

int
Can_Change::AnyMembersSet(const Set_of_Atoms& atom_numbers, int value) const {
  for (int atom_number : atom_numbers) {
    if (operator[](atom_number) == value) {
      return 1;
    }
  }

  return 0;
}

/*
  designed to replicate what a for loop does, so it does NOT return the max value
*/

class Iterate_Range
{
  private:
    int _mymax;
    int _ndx;
  public:
    Iterate_Range();
    Iterate_Range(const int s);

    void set_range(const int s);

    void reset() { _ndx = -1;}
    void reset(const int s) { _mymax = s; _ndx = -1;}

    int next(int &);
};

Iterate_Range::Iterate_Range()
{
  _mymax = -1;
  _ndx = -1;
}

Iterate_Range::Iterate_Range(const int s) : _mymax(s)
{
  _ndx = -1;
}

void 
Iterate_Range::set_range(const int s)
{
  _mymax = s;
  _ndx = -1;

  return;
}

int
Iterate_Range::next(int & s)
{
  if (_ndx == _mymax)
    return 0;

  _ndx++;

  if (_ndx == _mymax)
    return 0;

  return _ndx;
}

class Iterate_Atom_Pairs
{
  private:
    const int _matoms;
    int _a1;
    int _a2;

  public:
    Iterate_Atom_Pairs(const int s);

    int next(int &, int &);
};

Iterate_Atom_Pairs::Iterate_Atom_Pairs(const int s) :_matoms(s)
{
  _a1 = -1;
  _a2 = -1;
}

std::optional<Config>
ReadConfiguration(IWString& fname) {
  return iwmisc::ReadTextProtoCommentsOK<Config>(fname);
}

/*
  We can run in a mode where many of the sidechains must present themselves
  with the first atom being the attachment point - and having an implicit hydrogen.
  Or we can choose which atom to use as the join point every time the fragment is
  added
*/

static int choose_attachment_points_randomly = 0;


class Molecular_Transformation
{
  protected:
    const IWString _name;

    int _molecules_examined;
    int _molecules_changed;

    double _probability;

    RandomState _random_state;

//  protected functions

    int _check_probability();

    virtual int _do_transformation(Molecule & m, const Can_Change & can_change) = 0;

  public:
    Molecular_Transformation (const char *);
    virtual ~Molecular_Transformation ();

    int determine_probability(const const_IWSubstring &, double);

    int report(std::ostream &) const;

    int do_transformation(Molecule & m, const Can_Change & can_change);

    virtual void reset_iterator() = 0;
};

Molecular_Transformation::Molecular_Transformation(const char * s) : _name (s)
{
  _molecules_examined = 0;
  _molecules_changed = 0;
  _probability = 0.0;

  return;
}

Molecular_Transformation::~Molecular_Transformation ()
{
  _molecules_examined = -9;

  return;
}

int
Molecular_Transformation::determine_probability(const const_IWSubstring & token,
                                  double p)
{
  if (token == _name) {
    return 0;
  }

  assert (p >= 0.0 && p <= 1.0);

  _probability = p;

  return 1;
}

int
Molecular_Transformation::_check_probability()
{
  return _random_state.OkProbability(_probability);
}

class Molecular_Transformation_Add_Double_Bond : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Add_Double_Bond (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation(Molecule & m, const Can_Change &);
    int next_transformation(Molecule & m, const Can_Change &);
};

class Molecular_Transformation_Remove_Atom : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Remove_Atom (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Lower_Bond_Order : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Lower_Bond_Order (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};


class Molecular_Transformation_Break_Aliphatic_Ring : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Break_Aliphatic_Ring (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Make_Ring : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Make_Ring (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Add_Carbon : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Add_Carbon (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Move_Fragment : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Move_Fragment (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Swap_Adjacent_Atoms : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Swap_Adjacent_Atoms (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Change_Carbon_to_Nitrogen : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Change_Carbon_to_Nitrogen (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Change_Nitrogen_to_Carbon : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Change_Nitrogen_to_Carbon (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Change_Carbon_to_Oxygen : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Change_Carbon_to_Oxygen (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Change_Oxygen_to_Something : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Change_Oxygen_to_Something (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};


class Molecular_Transformation_Add_From_Single_Attachment_Library : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Add_From_Single_Attachment_Library (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Add_From_Aromatic_Attachment_Library : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Add_From_Aromatic_Attachment_Library (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Add_From_Double_Attachment_Library : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Add_From_Double_Attachment_Library (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Change_Ring_Substitution : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Change_Ring_Substitution (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Remove_Terminal_Functional_Group : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Remove_Terminal_Functional_Group (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Remove_Embedded_Functional_Group : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Remove_Embedded_Functional_Group (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Remove_Terminal_Ring : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Remove_Terminal_Ring (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Switch_Aromatic_Ring : public Molecular_Transformation
{
  private:

    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Switch_Aromatic_Ring (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation (Molecule &, const Can_Change &);
    int next_transformation (Molecule &, const Can_Change &);
};

class Molecular_Transformation_Switch_Aliphatic_Ring : public Molecular_Transformation
{
  private:
    Iterate_Range _iter;

//  private functions

    int _do_transformation (Molecule &, const Can_Change &);

  public:
    Molecular_Transformation_Switch_Aliphatic_Ring (const char *);

    void reset_iterator() { _iter.reset();}

    int random_transformation(Molecule &, const Can_Change &);
    int next_transformation(Molecule &, const Can_Change &);
};

Molecular_Transformation_Add_Double_Bond::Molecular_Transformation_Add_Double_Bond (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_add_double_bond;

  return;
}
Molecular_Transformation_Remove_Atom::Molecular_Transformation_Remove_Atom (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_remove_atom;

  return;
}

Molecular_Transformation_Lower_Bond_Order::Molecular_Transformation_Lower_Bond_Order (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_lower_bond_order;

  return;
}

Molecular_Transformation_Break_Aliphatic_Ring::Molecular_Transformation_Break_Aliphatic_Ring (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_break_aliphatic_ring;

  return;
}

Molecular_Transformation_Make_Ring::Molecular_Transformation_Make_Ring (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_make_ring;

  return;
}

Molecular_Transformation_Add_Carbon::Molecular_Transformation_Add_Carbon (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_add_carbon;

  return;
}

Molecular_Transformation_Move_Fragment::Molecular_Transformation_Move_Fragment (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_move_fragment;

  return;
}

Molecular_Transformation_Swap_Adjacent_Atoms::Molecular_Transformation_Swap_Adjacent_Atoms (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_swap_adjacent_atoms;

  return;
}

Molecular_Transformation_Change_Carbon_to_Nitrogen::Molecular_Transformation_Change_Carbon_to_Nitrogen (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_change_carbon_to_nitrogen;

  return;
}

Molecular_Transformation_Change_Nitrogen_to_Carbon::Molecular_Transformation_Change_Nitrogen_to_Carbon (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_change_nitrogen_to_carbon;

  return;
}

Molecular_Transformation_Change_Carbon_to_Oxygen::Molecular_Transformation_Change_Carbon_to_Oxygen (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_change_carbon_to_oxygen;

  return;
}

Molecular_Transformation_Change_Oxygen_to_Something::Molecular_Transformation_Change_Oxygen_to_Something (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_change_oxygen_to_something;

  return;
}

Molecular_Transformation_Add_From_Single_Attachment_Library::Molecular_Transformation_Add_From_Single_Attachment_Library (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_add_from_single_attachment_library;

  return;
}

Molecular_Transformation_Add_From_Aromatic_Attachment_Library::Molecular_Transformation_Add_From_Aromatic_Attachment_Library (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_add_from_aromatic_attachment_library;

  return;
}

Molecular_Transformation_Add_From_Double_Attachment_Library::Molecular_Transformation_Add_From_Double_Attachment_Library (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_add_from_double_attachment_library;

  return;
}

Molecular_Transformation_Change_Ring_Substitution::Molecular_Transformation_Change_Ring_Substitution (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_change_ring_substitution;

  return;
}

Molecular_Transformation_Remove_Terminal_Functional_Group::Molecular_Transformation_Remove_Terminal_Functional_Group (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_remove_terminal_functional_group;

  return;
}

Molecular_Transformation_Remove_Embedded_Functional_Group::Molecular_Transformation_Remove_Embedded_Functional_Group (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_remove_embedded_functional_group;

  return;
}

Molecular_Transformation_Remove_Terminal_Ring::Molecular_Transformation_Remove_Terminal_Ring (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_remove_terminal_ring;

  return;
}

Molecular_Transformation_Switch_Aromatic_Ring::Molecular_Transformation_Switch_Aromatic_Ring (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_switch_aromatic_ring;

  return;
}

Molecular_Transformation_Switch_Aliphatic_Ring::Molecular_Transformation_Switch_Aliphatic_Ring (const char * s) : Molecular_Transformation (s)
{
  _probability = probability_switch_aliphatic_ring;

  return;
}
















int
Molecular_Transformation_Add_Double_Bond::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Remove_Atom::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Lower_Bond_Order::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Break_Aliphatic_Ring::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Make_Ring::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Add_Carbon::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Move_Fragment::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Swap_Adjacent_Atoms::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Change_Carbon_to_Nitrogen::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Change_Nitrogen_to_Carbon::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Change_Carbon_to_Oxygen::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Change_Oxygen_to_Something::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Add_From_Single_Attachment_Library::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Add_From_Aromatic_Attachment_Library::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Add_From_Double_Attachment_Library::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Change_Ring_Substitution::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Remove_Terminal_Functional_Group::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Remove_Embedded_Functional_Group::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Remove_Terminal_Ring::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Switch_Aromatic_Ring::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

int
Molecular_Transformation_Switch_Aliphatic_Ring::_do_transformation (Molecule & m,
                                          const Can_Change & can_change)
{
  _molecules_examined++;

  return 1;
}

// Return a map from name to probability for all the
// entries in `config.probability`.
std::unordered_map<std::string, float>
Probabilities(const Config& config) {
  std::unordered_map<std::string, float> result;

  // const /*Descriptor*/auto* descriptor = config.probability().GetDescriptor();
  const /*Reflection*/auto* reflection = config.probability().GetReflection();
  using google::protobuf::FieldDescriptor;
  std::vector<const FieldDescriptor*> fields;
  reflection->ListFields(config.probability(), &fields);

  for (const FieldDescriptor* f : fields) {
    cerr << f->name() << '\n';
    float tmp = reflection->GetFloat(config.probability(), f);
    result[f->name()] = tmp;
  }

  return result;
}

class Transformations {
  private:
    int _verbose;

    // Read from the proto via the -N option.
    Config _config;

    std::vector<Molecular_Transformation> _transformation;

  public:
    Transformations();

    int Initialise(Command_Line& cl);
};

Transformations::Transformations() {
  _verbose = 0;
}

float
SumValues(const std::unordered_map<std::string, float>& hash) {
  float result= 0.0;
  for (auto& [_, value] : hash) {
    result += value;
  }
  return result;
}

int
Transformations::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('N')) {
    cerr << "Transformations::Initialise:must specify config file via the -N option\n";
    return 0;
  }

  if (cl.option_present('N')) {
    IWString fname = cl.string_value('N');
    std::optional<Config> tmp = ReadConfiguration(fname);
    if (! tmp) {
      cerr << "Transformations::Initialise:cannot read config '" << fname << "'\n";
      return 0;
    }
    _config = *tmp;
  }

  std::unordered_map<std::string, float> probabilities = 
        Probabilities(_config);

  float sum_probabilties = SumValues(probabilities);
  // Normalize the probabilities to the range 1000
  for (auto& [_, value] : probabilities) {
    value /= sum_probabilties * 1000.0f;
  }
  return 1;
}

static resizable_array_p<Molecule> input_molecules;

static int max_attempts = 10;

static resizable_array_p<Molecule> single_attachment_fragment_library_join_at_first_atom;
static resizable_array_p<Molecule> single_attachment_fragment_library_join_from_any_atom; 
static resizable_array_p<Molecule> aromatic_attachment_fragment_library;
static resizable_array_p<Molecule> double_attachment_fragment_library;
static resizable_array_p<Molecule> aromatic_rings;
static resizable_array_p<Molecule> aliphatic_rings;
static resizable_array_p<Molecule> aromatic_ring_fusions;

/*
  We can set a global maximum number of rings, or we can set the maximum
  change in rings
*/

static int min_rings_in_each_permutation = 0;
static int max_rings_in_each_permutation = 0;
static int molecules_with_too_few_rings = 0;
static int molecules_with_too_many_rings = 0;

static int max_rings_to_add = -1;
static int max_rings_to_remove = -1;

static resizable_array_p<Substructure_Hit_Statistics> queries_to_avoid;
static resizable_array_p<Substructure_Hit_Statistics> queries_to_match;

/*
  We can preserve a core by specifying a query
*/

static resizable_array_p<Substructure_Hit_Statistics> never_change_query;

/*
  We can prevent rings of certain sizes from forming
*/

static int minimum_ring_size_allowed = 0;
static int maximum_ring_size_allowed = 0;

/*
  During debugging, it is interesting to see the discarded molecules
*/

static Molecule_Output_Object stream_for_discarded_molecules;

/*
  We can keep track of which transformations are applied. This should really
  be passed as an argument, but too cumbersome
*/

static IWString changed_by;

static int changed_by_contains_full_history = 1;

#ifdef NOT_USED
static int
initialise_never_change_atoms (Molecule & m,
                               Can_Change & can_change)
{
  Molecule_to_Match target(&m);

  for (int i = 0; i < never_change_query.number_elements(); i++)
  {
    Substructure_Hit_Statistics * q = never_change_query[i];

    Substructure_Results sresults;
    int nhits = q->substructure_search(target, sresults);

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      e->set_vector(can_change.rawdata(), 0);
    }
  }

  return 1;
}
#endif

/*
  The convention is that each member of the double attachment library
  is attached via the first atom, and then via an asterisk atom.
  When the library is read in, we identify the asterisk atom, remove it
  and then store in this array the atom that was attached to the asterisk
*/

static resizable_array<int> index_of_star_atom_attachment;

static int
permutation_add_double_bond(Molecule & m,
                            const int bond_number,
                            const Can_Change & can_change)
{
  const Bond * b = m.bondi(bond_number);

  if (! b->is_single_bond())
    return 0;

  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  if (! can_change[a1] || ! can_change[a2])
    return 0;

  if (m.is_aromatic(a1) || m.is_aromatic(a2))
    return 0;

  if (0 == m.hcount(a1) || 0 == m.hcount(a2))
    return 0;

  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  if (aa1->nbonds() > aa1->ncon())
    return 0;
  if (aa2->nbonds() > aa2->ncon())
    return 0;

  m.set_bond_type_between_atoms(a1, a2, DOUBLE_BOND);

  if (verbose > 1)
    changed_by << " Double_bond";
  if (verbose > 2)
    cerr << "Double bond added between atoms " << a1 << " and " << a2 << '\n';

#ifdef FREQUENTLY_CHECK_VALENCES
  assert( m.valence_ok());
#endif

  return 1;
}

static int
permutation_add_double_bond(Molecule & m,
                            const Can_Change & can_change,
                            RandomState& random_state)
{
  if (! random_state.OkProbability(probability_add_double_bond)) {
    return 0;
  }

  const int nedges = m.nedges();

  for (int i = 0; i < max_attempts; i++)
  {
    const int bond_number = random_state.RandomInt(nedges - 1);

    if (permutation_add_double_bond(m, bond_number, can_change))
      return 1;
  }

  return 0;
}

static void
set_all_transformation_probabilities(const double p)
{
  probability_add_double_bond = p;
  probability_remove_atom = p;
  probability_lower_bond_order = p;
  probability_break_aliphatic_ring = p;
  probability_make_ring = p;
  probability_add_carbon = p;
  probability_move_fragment = p;
  probability_swap_adjacent_atoms = p;
  probability_change_carbon_to_nitrogen = p;
  probability_change_nitrogen_to_carbon = p;
  probability_change_carbon_to_oxygen = p;
  probability_change_oxygen_to_something = p;
  probability_add_from_single_attachment_library = p;
  probability_add_from_aromatic_attachment_library = p;
  probability_add_from_double_attachment_library = p;
  probability_change_ring_substitution = p;
  probability_remove_terminal_functional_group = p;
  probability_remove_embedded_functional_group = p;
  probability_remove_terminal_ring = p;
  probability_switch_aromatic_ring = p;
  probability_switch_aliphatic_ring = p;
  probability_breed = p;
  probability_breed_close_to_target = p;
  probability_split_fused_ring = p;
  probability_create_fused_ring = p;
  probability_expand_ring = p;
  probability_shrink_ring = p;

  return;
}

/*
  We are contemplating joining atoms A1 and A2
*/

static int
any_undesirable_adjacencies(Molecule & m,
                             atom_number_t a1,
                             atom_number_t a2)
{
  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  atomic_number_t z1 = aa1->atomic_number();
  atomic_number_t z2 = aa2->atomic_number();

  if (6 != z1 && 6 != z2)
    return 1;    // never join heteroatoms

// Avoid forming alkyl halides

  if (aa1->is_halogen())
  {
    if (! m.is_aromatic(a2))
      return 1;
  }
  else if (aa2->is_halogen())
  {
    if (! m.is_halogen(a1))
      return 1;
  }

// Never make adjacent carbonyls or such...

  if (m.multiple_bond_to_heteroatom(a1) && m.multiple_bond_to_heteroatom(a2))
    return 1;

  return 0;
}

static int
permutation_remove_atom(Molecule & m,
                        const atom_number_t j,
                        const Can_Change & can_change)
{
  const Atom * aj = m.atomi(j);

  if (2 != aj->ncon())
    return 0;

  if (2 != aj->nbonds())
    return 0;

  if (m.is_aromatic(j))
    return 0;

  if (m.in_ring_of_given_size(j, minimum_ring_size_allowed + 1))    // smallest possible ring
    return 0;

  atom_number_t a1 = aj->other(j, 0);
  atom_number_t a2 = aj->other(j, 1);

  if (m.are_bonded(a1, a2))    // existing 3 membered ring
    return 0;

  if (any_undesirable_adjacencies(m, a1, a2))
    return 0;

  m.remove_bond_between_atoms(a1, j);
  m.remove_bond_between_atoms(a2, j);
  if (a1 > j)
    a1--;
  if (a2 > j)
    a2--;

  m.remove_atom(j);
  m.add_bond(a1, a2, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Rm_atom";
  if (verbose > 2)
    cerr << "Removed atom " << j << '\n';

  return 1;
}

static int
permutation_remove_atom(Molecule & m,
                        const Can_Change & can_change,
                        RandomState& random_state)
{
  if (! random_state.OkProbability(probability_remove_atom)) {
    return 0;
  }
  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; ++i)
  {
    atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_remove_atom(m, j, can_change))
      return 1;
  }

  return 0;
}


/*
  We are contemplating changing the multiple bond between atoms A1 and
  A2 to a single bond.  If that bond is part of some desirable group,
  like a carboxyllic acid, don't do it
*/

static int
would_break_up_nice_group(const Molecule & m,
                          atom_number_t a1,
                          atom_number_t a2)
{
  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  atomic_number_t z1 = aa1->atomic_number();
  atomic_number_t z2 = aa2->atomic_number();

  if (6 == z1 && 6 == z2)
    return 0;

  if ((8 == z1 && 16 == z2) || (16 == z1 && 8 == z2))    // some kind of S=O
    return 1;

  if ((8 == z1 && 15 == z2) || (15 == z1 && 8 == z2))    // some kind of P=O
    return 1;

  if ((8 == z1 && 7 == z2) || (7 == z1 && 8 == z2))    // some kind of N=O
    return 1;

// Quat's can't be changed

  if (7 == z1 && 1 == aa1->formal_charge())
    return 1;

  if (7 == z2 && 1 == aa2->formal_charge())
    return 1;

  if (aa1->nbonds() > 4 || aa2->nbonds() > 4)    // don't touch these
    return 1;

// If we have an oxygen bonded to something that has another heteroatom attached, keep the double bond

  if (8 == z1 && 1 != m.attached_heteroatom_count(a1))
    return 1;

  if (8 == z2 && 1 != m.attached_heteroatom_count(a2))
    return 1;

  if (0 != aa1->formal_charge() || 0 != aa2->formal_charge())    // who knows!
    return 1;

  return 0;
}

static int
bonded_to_halogen(const Molecule & m,
                  atom_number_t zatom)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t x = b->other(zatom);

    if (m.is_halogen(x))
      return 1;
  }

  return 0;
}

static int
permutation_lower_bond_order(Molecule & m,
                             const int bnumber,
                             const Can_Change & can_change,
                             RandomState& random_state)
{
#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  const Bond * b = m.bondi(bnumber);

  if (b->is_single_bond()) {
    return 0;
  }

  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  if (! can_change[a1] || ! can_change[a2]) {
    return 0;
  }

  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  if (0 != aa1->formal_charge() || 0 != aa2->formal_charge())
    return 0;

  if (7 == aa1->atomic_number() && 5 == aa1->nbonds())
    return 0;

  if (7 == aa2->atomic_number() && 5 == aa2->nbonds())
    return 0;

  if (b->is_aromatic())
  {

    if (! random_state.OkProbability(probability_destroy_aromatic_ring)) {
      return 0;
    }
    if (bonded_to_halogen(m, a1) || bonded_to_halogen(m, a2))   // don't create alkyl halides
      return 0;
  }
  else if (would_break_up_nice_group(m, a1, a2)) {
    return 0;
  }

  m.set_bond_type_between_atoms(a1, a2, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Lower_bond";
  if (verbose > 2)
    cerr << "Bond between atoms " << a1 << " and " << a2 << " set to single\n";

#ifdef FREQUENTLY_CHECK_VALENCES
  if (! m.valence_ok())
    cerr << "Lowering bond order bad v " << m.smiles() << '\n';
  assert (m.valence_ok());
#endif

  return 1;
}

static int
permutation_lower_bond_order(Molecule & m,
                             const Can_Change & can_change,
                             RandomState& random_state)
{
#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  if (! random_state.OkProbability(probability_lower_bond_order)) {
    return 0;
  }

  const int nedges = m.nedges();

  for (int i = 0; i < max_attempts; i++)
  {
    atom_number_t j = random_state.RandomInt(nedges - 1);

    const Bond * b = m.bondi(j);

    if (! can_change[b->a1()] || ! can_change[b->a2()])
      continue;

    if (permutation_lower_bond_order(m, j, can_change, random_state))
      return 1;
  }

  return 0;
}

static int
permutation_break_aliphatic_ring(Molecule & m,
                                 const int bond_number,
                                 const Can_Change & can_change)
{
  const Bond * b = m.bondi(bond_number);

  if (0 == b->nrings())
    return 0;

  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  if (! can_change[a1] || ! can_change[a2])
    return 0;

  if (m.is_aromatic(a1) || m.is_aromatic(a2))
    return 0;

  m.remove_bond_between_atoms(a1, a2);

  if (verbose > 1)
    changed_by << " Break_Aliph";
  if (verbose > 2)
    cerr << "Broke aliphatic ring between atoms " << a1 << " and " << a2 << '\n';

  return 1;
}

static int
permutation_break_aliphatic_ring(Molecule & m,
                                 const Can_Change & can_change,
                                 RandomState& random_state)
{
  if (! random_state.OkProbability(probability_break_aliphatic_ring)) {
    return 0;
  }

  (void) m.ring_membership();    // force ring computation

  const int nedges = m.nedges();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(nedges - 1);

    if (permutation_break_aliphatic_ring(m, j, can_change))
      return 1;
  }

  return 0;
}

static int
count_aromatic_bonds(Molecule & m,
                     atom_number_t astart,
                     atom_number_t destination,
                     int encountered_so_far)
{
  const Atom * a = m.atomi(astart);

  int acon = a->ncon();

  int d = m.bonds_between(astart, destination);

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    atom_number_t j = b->other(astart);

    if (j == destination)
      return b->is_aromatic() + encountered_so_far;

    int dj = m.bonds_between(j, destination);
    if (dj >= d)     // not moving closer
      continue;

    int ab = count_aromatic_bonds(m, j, destination, encountered_so_far + b->is_aromatic());
    if (ab > 1)
      return ab;
  }

  return encountered_so_far;
}

/*
  We are thinking about making a bond between atoms A1 and A2. If they are
  both part of the same fused system, we don't want to do that
*/

static int
looks_like_different_parts_of_fused_system(Molecule & m,
                                           const atom_number_t a1,
                                           const atom_number_t a2)
{
  if (0 == m.nrings(a1))
    return 0;

  if (0 == m.nrings(a2))
    return 0;

  const Ring * r1 = m.ring_containing_atom(a1);

  if (! r1->is_fused())
    return 0;

  const Ring * r2 = m.ring_containing_atom(a2);

  if (! r2->is_fused())
    return 0;

  if (r1->fused_system_identifier() == r2->fused_system_identifier())   // in same system, bad...
    return 1;

  return 0;
}

/*
  We are contemplating forming a ring. We don't want to form something that has
  more than 2 consecutive existing aromatic atoms. Note that this will preclude
  forming an aliphatic ring that spans multiple aromatic rings
*/

static int
would_form_strained_aromatic(Molecule & m,
                             const atom_number_t a1,
                             const atom_number_t a2)
{
  m.longest_intra_molecular_distance();    // force distance matrix computation
  m.compute_aromaticity_if_needed();

  int aromatic_bonds = count_aromatic_bonds(m, a1, a2, 0);

  if (aromatic_bonds > 1)
    return 1;

  return 0;
}

static int
would_form_three_membered_ring_with_heteroatom(Molecule & m,
                                               atom_number_t a1,
                                               atom_number_t a2)
{
  if (2 != m.bonds_between(a1, a2))
    return 0;

  if (6 != m.atomic_number(a1))
    return 1;

  if (6 != m.atomic_number(a2))
    return 1;

// Find the atom between a1 and a2

  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  int acon = aa1->ncon();

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = aa1->other(a1, i);

    if (! aa2->is_bonded_to(j))
      continue;

    return (6 != m.atomic_number(j));
  }

  cerr << "would_form_three_membered_ring_with_heteroatom:should not come to here\n";
  return 1;
}

static int
permutation_make_ring(Molecule & m,
                      const atom_number_t a1,
                      const atom_number_t a2,
                      const int initial_strongly_fused_ring_count,
                      const Can_Change & can_change)
{
  if (0 == m.hcount(a1))
    return 0;

  if (m.nrings(a1) > 1)
    return 0;

  if (a2 == a1 || m.are_bonded(a2, a1) || 0 == m.hcount(a2) || m.nrings(a2) > 1)
    return 0;

  if (! can_change[a2])
    return 0;

  if (any_undesirable_adjacencies(m, a1, a2))
    return 0;

  if (m.bonds_between(a1, a2) <= minimum_ring_size_allowed - 1)
    return 0;
  else if (maximum_ring_size_allowed > 0 && m.bonds_between(a1, a2) >= maximum_ring_size_allowed - 1)
    return 0;
  else if (m.in_same_aromatic_ring(a1, a2))
    return 0;

  if (would_form_strained_aromatic(m, a1, a2))
    return 0;

  if (looks_like_different_parts_of_fused_system(m, a1, a2))
    return 0;

  if (would_form_three_membered_ring_with_heteroatom(m, a1, a2))
    return 0;

  if (m.in_same_ring(a1, a2))    // generate too many wierd things otherwise
    return 0;

  m.add_bond(a1, a2, SINGLE_BOND);

  if (can_form_strongly_fused_ring_systems)
    ;
  else if (m.rings_with_strongly_fused_ring_neighbours() > initial_strongly_fused_ring_count)
  {
    m.remove_bond_between_atoms(a1, a2);
    return 0;
  }

  if (verbose > 1)
    changed_by << " Make_ring";
  if (verbose > 2)
    cerr << "made ring by joining atoms " << a1 << " and " << a2 << '\n';

  return 1;
}

static int
permutation_make_ring(Molecule & m,
                      const Can_Change & can_change,
                      RandomState& random_state)
{
  if (! random_state.OkProbability(probability_make_ring)) {
    return 0;
  }

  if (max_rings_in_each_permutation > 0 && m.nrings() + 1 >= max_rings_in_each_permutation)
    return 0;

  if (min_rings_in_each_permutation > 0 && m.nrings() - 1 <= min_rings_in_each_permutation)
    return 0;

  int initial_strongly_fused_ring_count = 0;
  if (! can_form_strongly_fused_ring_systems)
    initial_strongly_fused_ring_count = m.rings_with_strongly_fused_ring_neighbours();

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    const atom_number_t k = random_state.RandomInt(matoms - 1);
    if (j == k)
      continue;

    if (! can_change[k])
      continue;

    if (permutation_make_ring(m, j, k, initial_strongly_fused_ring_count, can_change))
      return 1;
  }

  return 0;
}

static const Element * element_carbon = nullptr;
static const Element * element_oxygen = nullptr;
static const Element * element_nitrogen = nullptr;
static const Element * element_sulphur = nullptr;

static int
permutation_change_nitrogen_to_carbon(Molecule & m,
                                      const atom_number_t j,
                                      const Can_Change & can_change)
{
  if (7 != m.atomic_number(j))
    return 0;

  if (0 != m.formal_charge(j))
    return 0;

  if (m.nbonds(j) > 4)
    return 0;

  m.set_element(j, element_carbon);

  if (verbose > 1)
    changed_by << " N2C";
  if (verbose > 2)
    cerr << "Nitrogen changed to Carbon\n";

  return 1;
}


static int
permutation_change_nitrogen_to_carbon(Molecule & m,
                                      const Can_Change & can_change,
                                      RandomState& random_state)
{
  if (! random_state.OkProbability(probability_change_nitrogen_to_carbon)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_change_nitrogen_to_carbon(m, j, can_change))
      return 1;
  }

  return 0;
}

static int
singly_bonded_to_oxygen_or_nitrogen(const Molecule & m,
                                    const atom_number_t a)
{
  const Atom * zatom = m.atomi(a);

  for (int i = 0; i < zatom->ncon(); i++)
  {
    atom_number_t j = zatom->other(a, i);

    if (7 == m.atomic_number(j))
      return 1;
    if (8 == m.atomic_number(j))
      return 1;
  }

  return 0;
}

static int
alpha_to_oxygen_or_nitrogen(Molecule & m, 
                            const atom_number_t astart)
{
  const Atom * a = m.atomi(astart);

  const int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const atom_number_t j = a->other(astart, i);

    const Atom * aj = m.atomi(j);

    const atomic_number_t zj = aj->atomic_number();

    if (6 == zj)
      continue;

    if (m.is_aromatic(j))
      continue;

    if (7 == zj || 8 == zj || 15 == zj || 16 == zj)
    {
      if (aj->ncon() == aj->nbonds())    // fully saturated
        return 1;
    }
  }

  return 0;
}

/*
  We are contemplating changing a carbon to an oxygen.  The carbon has
  no attached heteroatoms, but we want to avoid Nitrogen and Oxygen
  and Sulphur atoms in the beta position
*/

static int
beta_to_oxygen_or_nitrogen(Molecule & m, 
                           const atom_number_t astart)
{
  const Atom * a = m.atomi(astart);

  assert (6 == a->atomic_number());

  const int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(astart, i);

    if (0 == m.attached_heteroatom_count(j))
      continue;

    if (alpha_to_oxygen_or_nitrogen(m, j))
      return 1;
  }

  return 0;
}

static int
permutation_change_carbon_to_nitrogen(Molecule & m,
                                      const atom_number_t j,
                                      const Can_Change & can_change)
{
  if (6 != m.atomic_number(j))
    return 0;

  if (4 == m.ncon(j) || 4 == m.nbonds(j))
    return 0;

  if (singly_bonded_to_oxygen_or_nitrogen(m, j))
    return 0;

  if (beta_to_oxygen_or_nitrogen(m, j))
    return 0;

  if (m.in_ring_of_given_size(j, 3))
    return 0;

  m.set_element(j, element_nitrogen);

  if (verbose > 1)
    changed_by << " C2N";
  if (verbose > 2)
    cerr << "Carbon changed to Nitrogen\n";

  return 1;
}

static int
permutation_change_carbon_to_nitrogen(Molecule & m,
                                      const Can_Change & can_change,
                                      RandomState& random_state)
{
  if (! random_state.OkProbability(probability_change_carbon_to_nitrogen)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_change_carbon_to_nitrogen(m, j, can_change))
      return 1;
  }

  return 0;
}

/*
  Check whether or not atom J is bonded to any element of type Z
*/

/*static int
is_bonded_to_element (const Molecule & m,
                      atom_number_t j,
                      atomic_number_t z)
{
  const Atom * aj = m.atomi(j);

  for (int i = 0; i < aj->ncon(); i++)
  {
    atom_number_t o = aj->other(j, i);

    if (z == m.atomic_number(o))
      return 1;
  }

  return 0;
}*/

static int
permutation_change_carbon_to_oxygen(Molecule & m,
                                    const atom_number_t j,
                                    const Can_Change & can_change)
{
  if (6 != m.atomic_number(j))
    return 0;

  if (m.ncon(j) > 2 || m.nbonds(j) > 2)
    return 0;

  if (0 != m.attached_heteroatom_count(j))
    return 0;

  if (m.in_ring_of_given_size(j, 3))
    return 0;

  if (! can_change[j])
    return 0;

  if (beta_to_oxygen_or_nitrogen(m, j))
    return 0;

  m.set_element(j, element_oxygen);

  if (verbose > 1)
    changed_by << " C2O";
  if (verbose > 2)
    cerr << "Carbon changed to oxygen\n";

  return 1;
}

static int
permutation_change_carbon_to_oxygen(Molecule & m,
                                    const Can_Change & can_change,
                                    RandomState& random_state)
{
  if (! random_state.OkProbability(probability_change_carbon_to_oxygen)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_change_carbon_to_oxygen(m, j, can_change))
      return 1;
  }

  return 0;
}

static int
permutation_change_oxygen_to_something(Molecule & m,
                                       const atom_number_t j,
                                       const Can_Change & can_change,
                                       RandomState& random_state)
{
  const Atom * aj = m.atomi(j);

  if (8 != aj->atomic_number())
    return 0;

  if (m.multiple_bond_to_heteroatom(j))
    return 0;

  if (1 == aj->ncon() && 2 == aj->nbonds())   // don't change carbonyls
    return 0;

  const float rand = random_state();

  if (rand < 0.2)
    m.set_element(j, element_nitrogen);
  else if (2 == aj->ncon() && rand < 0.4)   // change ether-type to sulphur
    m.set_element(j, element_sulphur);
  else
    m.set_element(j, element_carbon);

  if (verbose > 1)
    changed_by << " O2*";
  if (verbose > 2)
    cerr << "Oxygen atom switched\n";

  return 1;
}

static int
permutation_change_oxygen_to_something(Molecule & m,
                                      const Can_Change & can_change,
                                      RandomState& random_state)
{
  if (! random_state.OkProbability(probability_change_oxygen_to_something)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_change_oxygen_to_something(m, j, can_change, random_state))
      return 1;
  }

  return 0;
}

static int
permutation_add_carbon(Molecule & m,
                       const atom_number_t j,
                       const Can_Change & can_change)
{
  if (0 == m.hcount (j))
    return 0;

  m.add(element_carbon);
  m.add_bond(j, m.natoms() - 1, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Add_C";
  if (verbose > 2)
    cerr << "Added carbon atom\n";

  return 1;
}

static int
permutation_add_carbon(Molecule & m,
                       const Can_Change & can_change,
                       RandomState& random_state)
{
  if (! random_state.OkProbability(probability_add_carbon)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_add_carbon(m, j, can_change))
      return 1;
  }

  return 0;
}

static int
identify_chain_single_bond(Molecule & m,
                           atom_number_t & a1,
                           atom_number_t & a2,
                           RandomState& random_state)
{
  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const int j = random_state.RandomInt(matoms - 1);

    const Atom * aj = m.atomi(j);

    const int jcon = aj->ncon();

    for (int k = 0; k < jcon; k++)
    {
      const Bond * bk = aj->item (k);
      if (! bk->is_single_bond ())
        continue;

      if (bk->nrings())
        continue;

      atom_number_t l = bk->other(j);

//    don't make any checks, allow aldehydes to form, hopefully they'll disappear later

//    if (m.is_aromatic(l))    // ok to break these
//      ;
//    else if (m.nbonds(l) > m.ncon(l))    // unsaturated group, don't break off from it - creates aldehydes and such
//      continue;

      a1 = j;
      a2 = l;

      return 1;
    }
  }

  return 0;
}

static int
identify_join_atom_in_other_fragment(Molecule & m, 
                                     const atom_number_t & a3,
                                     atom_number_t & a4,
                                     RandomState& random_state)
{
  const int matoms = m.natoms();

  const int f = m.fragment_membership(a3);

  for (int i = 0; i < max_attempts; i++)
  {
    const int j = random_state.RandomInt(matoms - 1);

    if (f == m.fragment_membership(j))
      continue;

    if (0 == m.hcount(j))
      continue;

    if (3 == m.ncon(j))
      continue;

    if (m.nrings(j) > 1)
      continue;

    a4 = j;
    return 1;
  }

// Try an exhaustive search

  for (int i = 0; i < matoms; i++)
  {
    if (f == m.fragment_membership(i))
      continue;

    if (0 == m.hcount(i))
      continue;

    if (3 == m.ncon(i))
      continue;

    if (m.nrings(i) > 1)
      continue;

    a4 = i;
    return 1;
  }

  return 0;
}

/*
  Break a single chain bond and re-attach one of the ends somewhere else
*/

static int
permutation_move_fragment(Molecule & m,
                          const atom_number_t a1,
                          const atom_number_t a2,
                          const Can_Change & can_change,
                          RandomState& random_state)
{
  if (! can_change[a1] || ! can_change[a2]) {
    return 0;
  }

  m.remove_bond_between_atoms(a1, a2);

  atom_number_t a3;
  if (random_state() < 0.5)
    a3 = a1;
  else
    a3 = a2;
  
  atom_number_t a4 = INVALID_ATOM_NUMBER;
  if (! identify_join_atom_in_other_fragment(m, a3, a4, random_state)) {
    m.add_bond(a1, a2, SINGLE_BOND);
    return 0;
  }

  if (any_undesirable_adjacencies(m, a3, a4)) {
    m.add_bond(a1, a2, SINGLE_BOND);
    return 0;
  }

  m.add_bond(a3, a4, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Move_Frag";
  if (verbose > 2)
    cerr << "moved fragment " << a1 << ',' << a2 << ',' << a3 << ',' << a4 << '\n';

  return 1;
}


static int
permutation_move_fragment(Molecule & m,
                          const Can_Change & can_change,
                          RandomState& random_state)
{
  if (! random_state.OkProbability(probability_move_fragment)) {
    return 0;
  }

  m.ring_membership();

  for (int i = 0; i < max_attempts; i++) {
    atom_number_t a1, a2;

    if (! identify_chain_single_bond(m, a1, a2, random_state))
      continue;

    if (permutation_move_fragment(m, a1, a2, can_change, random_state))
      return 1;
  }

  return 0;
}

/*
  Common code for swapping adjacent atoms. It can be called by the do and undo transformation
     A1-J-A2-A3
  becomes
     A1-A2-J-A3
*/

static void
swap_adjacent_atoms(Molecule & m,
                    const atom_number_t a1,
                    const atom_number_t j,
                    const atom_number_t a2,
                    const atom_number_t a3)
{
  m.remove_bond_between_atoms(a1, j);
  m.remove_bond_between_atoms(j, a2);
  m.remove_bond_between_atoms(a2, a3);

  m.add_bond(a1, a2, SINGLE_BOND);
  m.add_bond(a2, j, SINGLE_BOND);
  m.add_bond(j, a3, SINGLE_BOND);

  return;
}

/*
  We have swapped adjacent atoms, but this may lead to a violation of ring size
  constraints. Too hard to figure out ahead of time.
  Put the molecule back to its original state if we've violated a constraint
*/

static int
violates_ring_size_constraints(Molecule & m,
                               const atom_number_t a1,
                               const atom_number_t j,
                               const atom_number_t a2,
                               const atom_number_t a3)
{
  if (0 == minimum_ring_size_allowed && 0 == maximum_ring_size_allowed)
    return 0;     // obviously no violation

  const int nr = m.nrings();

  if (0 == nr)
    return 0;    // no violation here either

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    const int ring_size = ri->number_elements();

    if (ring_size <= minimum_ring_size_allowed)
    {
      swap_adjacent_atoms(m, a1, a2, j, a3);
      return 1;
    }

    if (maximum_ring_size_allowed > 0 && ring_size >= maximum_ring_size_allowed)
    {
      swap_adjacent_atoms(m, a1, a2, j, a3);
      return 1;
    }
  }
  
  return 0;      // no violations found
}

/*
  We are just about to swap atoms A1-J-A2-A3
   to A1-A2-J-A3
*/

static int
any_undesirable_adjacencies(Molecule & m,
                            const atom_number_t a1,
                            const atom_number_t j,
                            const atom_number_t a2,
                            const atom_number_t a3)
{
  if (any_undesirable_adjacencies(m, a1, a2))
    return 1;

  if (any_undesirable_adjacencies(m, j, a3))
    return 1;

  return 0;
}

static int
is_part_of_ring_system(Molecule & m,
                       const atom_number_t zatom)
{
  if (m.nrings(zatom) > 1)
    return 1;

  const Ring * r = m.ring_containing_atom(zatom);

  assert (nullptr != r);

  if (r->fused_system_identifier() >= 0)
    return 1;

  return 0;
}

/*
  Find some atoms either side of J
*/

#ifdef BY_ATOMS_QWEQWEQWE
static int
permutation_swap_adjacent_atoms(Molecule & m,
                                const atom_number_t j,
                                const Can_Change & can_change)
{
  if (m.is_aromatic(j))
    return 0;

  if (allow_fused_rings_to_swap_atoms)
    ;
  else if (0 == m.nrings(j))
    ;
  else if (is_part_of_ring_system(m, j))
    return 0;

  const Atom * aj = m.atomi(j);

  int jcon = aj->ncon();

  if (1 == jcon)
    return 0;

// We need to identify a sequence A1 - J - A2 ---

  atom_number_t a1 = INVALID_ATOM_NUMBER;
  atom_number_t a2 = INVALID_ATOM_NUMBER;

  for (int k = 0; k < jcon; k++)
  {
    const Bond * bk = aj->item(k);

    if (! bk->is_single_bond())
      continue;

    const atom_number_t l = bk->other(j);

    const Atom * al = m.atomi(l);

    if (al->atomic_number() == aj->atomic_number() && al->ncon() == aj->ncon() && al->nbonds() == aj->nbonds())     // the look to be the same
      continue;

    if (INVALID_ATOM_NUMBER == a1)
      a1 = l;
    else 
    {
      a2 = l;
      break;
    }
  }

  if (INVALID_ATOM_NUMBER == a1 || INVALID_ATOM_NUMBER == a2)
    return 0;

  if (! can_change[a1] || ! can_change[a2])
    return 0;

// Now we need to identify A3 :  A1 - J - A2 - A3

  atom_number_t a3 = INVALID_ATOM_NUMBER;

  const Atom * am = m.atomi(a2);

  int mcon = am->ncon();

  for (int k = 0; k < mcon; k++)
  {
    const Bond * b = am->item (k);
    if (! b->is_single_bond ())
      continue;

    atom_number_t l = b->other(a2);
    if (j == l)
      continue;

    a3 = l;
    break;
  }

  if(INVALID_ATOM_NUMBER == a3)
    return 0;

  if (m.are_bonded(a1, a2) || m.are_bonded(j, a3) || m.are_bonded(a1, a3))    // possible existing small ring
    return 0;

  if (any_undesirable_adjacencies(m, a1, j, a2, a3))
    return 0;

//  Great we got it. Change A1-J-A2-A3 to A1-A2-J-A3

  swap_adjacent_atoms(m, a1, j, a2, a3);

  if (violates_ring_size_constraints(m, a1, j, a2, a3))    // will restore it if needed
    return 0;

  if (verbose > 1)
    changed_by << " Swap_Adj";
  if (verbose > 2)
    cerr << "Swapping adjacent atoms " << a1 << ',' << j << ',' << a2 << ',' << a3 << '\n';

  return 1;
}
#endif

static atom_number_t
get_connected_atom(const Atom * a,
                   const atom_number_t zatom,
                   const atom_number_t exclude)
{
  const int acon = a->ncon();

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);

    if (j == exclude)
      continue;

    return j;
  }

  return INVALID_ATOM_NUMBER;
}

static int
permutation_swap_adjacent_atoms(Molecule & m,
                                const int bond_number,
                                const Can_Change & can_change)
{
  const Bond * b = m.bondi(bond_number);

  if (b->is_aromatic())
    return 0;

  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  if (! can_change[a1] || ! can_change[a2])
    return 0;

  if (allow_fused_rings_to_swap_atoms)
    ;
  else if (0 == b->nrings())
    ;
  else if (is_part_of_ring_system(m, a1))
    return 0;

  const Atom * aa1 = m.atomi(a1);
  const Atom * aa2 = m.atomi(a2);

  if (1 == aa1->ncon() || 1 == aa2->ncon())
    return 0;

  if ((aa1->atomic_number() == aa2->atomic_number()) &&    // nothing would change if we swapped them
     (aa1->ncon() == aa2->ncon()) &&
     (aa1->nbonds() == aa2->nbonds()))
    return 0;

// Get some atoms outside the bond

  const atom_number_t a0 = get_connected_atom(aa1, a1, a2);
  const atom_number_t a3 = get_connected_atom(aa2, a2, a1);

  assert (INVALID_ATOM_NUMBER != a0);
  assert (INVALID_ATOM_NUMBER != a3);

  if (! can_change[a0] || ! can_change[a3])
    return 0;

  if (m.are_bonded(a0, a2) || m.are_bonded(a1, a3))    // possible existing small ring
    return 0;

  if (any_undesirable_adjacencies(m, a0, a1, a2, a3))
    return 0;

//  Great we got it. Change A1-J-A2-A3 to A1-A2-J-A3

  swap_adjacent_atoms(m, a0, a1, a2, a3);

  if (violates_ring_size_constraints(m, a0, a1, a2, a3))    // will restore it if needed
    return 0;

  if (verbose > 1)
    changed_by << " Swap_Adj";
  if (verbose > 2)
    cerr << "Swapping adjacent atoms " << a0 << ',' << a1 << ',' << a2 << ',' << a3 << '\n';

  return 1;
}

static int
permutation_swap_adjacent_atoms(Molecule & m,
                                const Can_Change & can_change,
                                RandomState& random_state)
{
  if (! random_state.OkProbability(probability_swap_adjacent_atoms)) {
    return 0;
  }

  m.compute_aromaticity_if_needed();

  const int nedges = m.nedges();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(nedges - 1);

    if (permutation_swap_adjacent_atoms(m, j, can_change))
      return 1;
  }

  return 0;
}

/*
  When adding things from the fragment libraries, we need to avoid certain groupings
*/

static int
any_undesirable_adjacencies(const Molecule & m1,
                             atom_number_t a1,
                             const Molecule & m2,
                             atom_number_t a2)
{
  const Atom * aa1 = m1.atomi(a1);
  const Atom * aa2 = m2.atomi(a2);

  const atomic_number_t z1 = aa1->atomic_number();
  const atomic_number_t z2 = aa2->atomic_number();

  if (6 != z1 && 6 != z2)    // never join heteroatoms
    return 1;

// Avoid very crowded things

  if (aa1->ncon() + aa2->ncon() > 4)
    return 1;

  return 0;
}

static atom_number_t 
choose_random_attachment_point(Molecule & m,
                               RandomState& random_state)
{
  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    atom_number_t rc = random_state.RandomInt(matoms - 1);

    int h = m.hcount(rc);

    if (0 == h)
      continue;

    atomic_number_t z = m.atomic_number(rc);

    if (6 == z || 7 == z || 8 == z)
      ;
    else if (9 == z || 17 == z || 35 == z || 53 == z)   // halogens
      continue;
    else if (random_state() > 0.2)   // allow some possibility for joining
      continue;

//  dis-favour forming 4 connected, fully saturated carbons

    if (1 == h && 3 == m.ncon(rc))
    {
      if (m.nrings(rc))   // never make these
        continue;
      
      if (random_state() > 0.2)
        continue;
    }

    if (m.nrings(rc) > 1)
      continue;

    return rc;
  }

  return INVALID_ATOM_NUMBER;   
}

static int
permutation_add_from_single_attachment_library(Molecule & m,
                                               const atom_number_t j1,
                                               const Molecule * f,
                                               const atom_number_t j2)
{
  if (any_undesirable_adjacencies(m, j1, *f, j2))
    return 0;

  const int initial_matoms = m.natoms();

  m.add_molecule (f);
  m.add_bond(j1, initial_matoms + j2, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Sngle_Attch";
  if (verbose > 2)
    cerr << "Added single attachment fragment library member '" << f->name () << "'\n";

// if (random_atom_in_fragment)
//   cerr << "random_atom_in_fragment produced " << m.smiles() << "' from '" << f->smiles() << "'\n";

  return 1;
}

static int
permutation_add_from_single_attachment_library(Molecule & m,
                                        const Can_Change & can_change,
                                        const resizable_array_p<Molecule> & library,
                                        int random_atom_in_fragment,
                                        RandomState& random_state)
{
  if (library.empty())
    return 0;

  for (int i = 0; i < max_attempts; i++)
  {
    atom_number_t j1 = choose_random_attachment_point(m, random_state);

    if (INVALID_ATOM_NUMBER == j1)
      continue;

    Molecule * f = random_state.RandomItem(library);   // we do not change it, but choose_random_attachment_point uses some non const methods

    atom_number_t j2;

    if (random_atom_in_fragment)
    {
      j2 = choose_random_attachment_point(*f, random_state);
      if (INVALID_ATOM_NUMBER == j2)
        continue;
    }
    else
      j2 = 0;

    if (permutation_add_from_single_attachment_library(m, j1, f, j2))
      return 1;
  }

  return 0;
}

static int
permutation_add_from_single_attachment_library(Molecule & m,
                                               const Can_Change & can_change,
                                               RandomState& random_state)
{
  if (! random_state.OkProbability(probability_add_from_single_attachment_library)) {
    return 0;
  }

  return permutation_add_from_single_attachment_library(m, can_change,
              single_attachment_fragment_library_join_at_first_atom,
              choose_attachment_points_randomly,
              random_state);  // 0 means take first atom in fragment
}

static int
permutation_add_from_single_attachment_library_random_join(Molecule & m,
                                             const Can_Change & can_change,
                                             RandomState& random_state)
{
  if (! random_state.OkProbability(probability_add_from_single_attachment_library)) {
    return 0;
  }

  return permutation_add_from_single_attachment_library(m, can_change,
                        single_attachment_fragment_library_join_from_any_atom,
                        1,
                        random_state);   // 1 means use random atom in fragment
}

static int
permutation_add_from_aromatic_attachment_library(Molecule & m,
                                                 const atom_number_t j,
                                                 Molecule * f,
                                                 const Can_Change & can_change,
                                                 RandomState& random_state)
{
  if (0 == m.hcount (j))
    return 0;

// if (m.is_aromatic (j))    // don't form biphenyls
//   continue;

  atom_number_t a;
  if (choose_attachment_points_randomly)
    a = choose_random_attachment_point(*f, random_state);
  else
    a = 0;

  if (any_undesirable_adjacencies(m, j, *f, a))
    return 0;

  int initial_matoms = m.natoms();

  m.add_molecule (f);
  m.add_bond(j, initial_matoms + a, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Arom_Attch";
  if (verbose > 2)
    cerr << "Added aromatic attachment fragment library member '" << f->name () << "'\n";

  return 1;
}

static int
permutation_add_from_aromatic_attachment_library(Molecule & m,
                                                 const Can_Change & can_change,
                                                 RandomState& random_state)
{
  if (! random_state.OkProbability(probability_add_from_aromatic_attachment_library)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    atom_number_t j = random_state.RandomInt(matoms - 1);

    if (permutation_add_from_aromatic_attachment_library(m, j, 
             random_state.RandomItem(aromatic_attachment_fragment_library),
             can_change, random_state))
      return 1;
  }

  return 0;
}

static int
permutation_add_from_double_attachment_library(Molecule & m,
                                               const atom_number_t j,
                                               const Molecule * f,
                                               const atom_number_t star_atom,
                                               const Can_Change & can_change)
{
  const Atom * aj = m.atomi(j);
  if (6 != aj->atomic_number())
    return 0;

  const int jcon = aj->ncon();

  for (int k = 0; k < jcon; k++)
  {
    const Bond * b = aj->item (k);

    atom_number_t l = b->other(j);

    if (m.in_same_aromatic_ring(j, l))
      continue;

    if (b->is_single_bond ())
      ;
    else if (b->is_double_bond () && (m.is_aromatic (j) || m.is_aromatic (l)))
      continue;

    if (any_undesirable_adjacencies(m, j, *f, 0) || any_undesirable_adjacencies (m, l, *f, star_atom))
      continue;

    int initial_matoms = m.natoms();
    m.add_molecule(f);
    m.remove_bond_between_atoms(j, l);
    m.add_bond(j, initial_matoms, SINGLE_BOND);
    m.add_bond(l, initial_matoms + star_atom, SINGLE_BOND);

    if (verbose > 1)
      changed_by << " Add_Dbl_Attch";
    if (verbose > 2)
      cerr << "Added doubly attached fragment library member '" << f->name () << "'\n";

    return 1;
  }

  return 0;
}

static int
permutation_add_from_double_attachment_library(Molecule & m,
                                               const Can_Change & can_change,
                                               RandomState& random_state)
{
  if (! random_state.OkProbability(probability_add_from_double_attachment_library))
    return 0;

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; ++i)
  {
    const int n = random_state.RandomInt(double_attachment_fragment_library.number_elements() - 1);

    const Molecule * f = double_attachment_fragment_library[n];

    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (permutation_add_from_double_attachment_library(m, j, f, index_of_star_atom_attachment[n], can_change))
      return 1;
  }

  return 0;
}


static int
permutation_change_ring_substitution(Molecule & m,
                                     const int r,
                                     const Can_Change & can_change)
{
  atom_number_t has_hydrogen = INVALID_ATOM_NUMBER;
  atom_number_t non_ring_bond1 = INVALID_ATOM_NUMBER;
  atom_number_t non_ring_bond2 = INVALID_ATOM_NUMBER;

  const Ring * ri = m.ringi(r);

  for (int j = 0; j < ri->number_elements(); j++)
  {
    atom_number_t k = ri->item(j);

    if (m.hcount(k)) {
      has_hydrogen = k;
      continue;
    }

    const Atom * ak = m.atomi(k);

    if (3 != ak->ncon()) {
      continue;
    }

    for (const Bond* b : *ak) {
      if (b->nrings()) {
        continue;
      }

      non_ring_bond1 = k;
      non_ring_bond2 = b->other(k);
      break;
    }
  }

  if (INVALID_ATOM_NUMBER == has_hydrogen ||
      INVALID_ATOM_NUMBER == non_ring_bond1) {
    return 0;
  }

  if (! can_change[has_hydrogen] && ! can_change[non_ring_bond1] && ! can_change[non_ring_bond2]) {
    return 0;
  }

  m.remove_bond_between_atoms(non_ring_bond1, non_ring_bond2);
  m.add_bond(has_hydrogen, non_ring_bond2, SINGLE_BOND);

  if (verbose > 1)
    changed_by << " Chg_Ring_Sub";
  if (verbose > 2)
    cerr << "Changed ring substitution, atoms " << '\n';

#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  return 1;
}

static int
permutation_change_ring_substitution(Molecule & m,
                                     const Can_Change & can_change,
                                     RandomState& random_state)
{
  if (1.0 == probability_change_ring_substitution)
    ;
  else if (! random_state.OkProbability(probability_change_ring_substitution)) {
    return 0;
  }

  const int nr = m.nrings();
  if (nr == 0) {
    return 0;
  }

  for (int i = 0; i < max_attempts; i++)
  {
    const int r = random_state.RandomInt(nr - 1);

    if (permutation_change_ring_substitution(m, r, can_change))
      return 1;
  }

  return 0;
}

static int
identify_terminal_functional_group(Molecule & m,
                                   const atom_number_t xstart,
                                   Set_of_Atoms & s)
{
  s.add(xstart);

  atom_number_t a = m.other(xstart, 0);
  s.add(a);

  const Atom * a1 = m.atomi(a);

  int single_bond_to_carbon_found = 0;

  for (const Bond * b : *a1) {
    const atom_number_t j = b->other(a);

    if (j == xstart) {
      continue;
    }

    const int jcon = m.ncon(j);

    if (1 == jcon) {
      s.add(j);
      continue;
    }

    if (b->is_single_bond () && 6 == m.atomic_number(j)) {
      single_bond_to_carbon_found++;
      continue;
    }

    return 0;          // some other type of connection, cannot remove this group of atoms
  }

  return 1 == single_bond_to_carbon_found;
}

/*
  Find a singly connected atom that is multiply bonded to something and look for a functional
  group to get rid of. Should trim things like cyano, acids, etc...
*/

static int
permutation_remove_terminal_functional_group(Molecule & m,
                                             const atom_number_t j,
                                             const Can_Change & can_change)
{
  const Atom * aj = m.atomi(j);

  if (1 != aj->ncon())
    return 0;

  if (1 == aj->nbonds())    // we start with unsaturated atoms
    return 0;

  Set_of_Atoms s;

  if (! identify_terminal_functional_group(m, j, s)) {
    return 0;
  }

  m.remove_atoms(s);

  if (verbose > 1)
    changed_by << " Mv_Terminal";
  if (verbose > 2)
    cerr << "Removed terminal functional group with " << s.number_elements() << " atoms\n";

  return 1;
}

static int
permutation_remove_terminal_functional_group(Molecule & m,
                                             const Can_Change & can_change,
                                             RandomState& random_state)
{
  if (! random_state.OkProbability(probability_remove_terminal_functional_group)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (! can_change[j])
      continue;

    if (permutation_remove_terminal_functional_group(m, j, can_change))
      return 1;
  }

  return 0;
}

static int
remove_embedded_functional_group(Molecule & m,
                                 atom_number_t zstart)
{
  Set_of_Atoms s;    // the atoms that will be removed

  s.add(zstart);

  atom_number_t middle = m.other(zstart, 0);

  s.add(middle);

  const Atom * am = m.atomi(middle);

// Since we are excising the group, we need the attachment points on either side

//#define DEBUG_EXCISE_FUNCTIONAL_GROUP
#ifdef DEBUG_EXCISE_FUNCTIONAL_GROUP
  cerr << "Starting with atom " << zstart << '\n';
#endif

  atom_number_t singly_connected_carbon1 = INVALID_ATOM_NUMBER;
  atom_number_t singly_connected_carbon2 = INVALID_ATOM_NUMBER;

  for (const Bond*b : *am) {
    atom_number_t j = b->other(middle);

    if (j == zstart) {
      continue;
    }

    const Atom * aj = m.atomi(j);

    if (1 == aj->ncon()) {
      s.add(j);
      continue;
    }

    if (b->is_single_bond () && 6 == aj->atomic_number()) {
      if (INVALID_ATOM_NUMBER == singly_connected_carbon1)
        singly_connected_carbon1 = j;
      else if (INVALID_ATOM_NUMBER == singly_connected_carbon2)
        singly_connected_carbon2 = j;
      else
        return 0;          // more than 3 singly bonded carbon connections
    }

    if (INVALID_ATOM_NUMBER != singly_connected_carbon2)
      return 0;

//  think of the Nitrogen in an amide. We want to remove the Nitrogen and identify the atom next in the chain from it

    if (6 != aj->atomic_number() && 2 == aj->ncon() && 2 == aj->nbonds()) {
      atom_number_t atom_on_other_side_of_j = aj->other(j, 0);
      if (middle == atom_on_other_side_of_j)
        atom_on_other_side_of_j = aj->other(j, 1);

      if (INVALID_ATOM_NUMBER == singly_connected_carbon1)
        singly_connected_carbon1 = atom_on_other_side_of_j;
      else if (INVALID_ATOM_NUMBER == singly_connected_carbon2)
        singly_connected_carbon2 = atom_on_other_side_of_j;

      s.add(j);
      continue;
    }

    return 0;      // not sure what this must be
  }

#ifdef DEBUG_EXCISE_FUNCTIONAL_GROUP
  cerr << "Group is " << s << '\n';
#endif

  // we need to have two connections
  if (INVALID_ATOM_NUMBER == singly_connected_carbon1 ||
      INVALID_ATOM_NUMBER == singly_connected_carbon2) {
    return 0;
  }

// If we are in a 5 membered ring, this can happen

  if (m.are_bonded(singly_connected_carbon1, singly_connected_carbon2)) {
    return 0;
  }

  m.add_bond(singly_connected_carbon1, singly_connected_carbon2, SINGLE_BOND);   // temporary invalid valence

  m.remove_atoms (s);

  if (verbose > 1)
    changed_by << " Rm_Embedded";
  if (verbose > 2)
    cerr << "Removed embedded functional group with " << s.number_elements() << " atoms\n";

  return 1;
}

/*
  Embedded functional groups are things like an amide, or a -S(=O)(=O)-
*/

static int
permutation_remove_embedded_functional_group(Molecule & m,
                                             const atom_number_t j,
                                             const Can_Change & can_change)
{
  if (! can_change[j]) {
    return 0;
  }

  const Atom * aj = m.atomi(j);

  if (1 != aj->ncon()) {
    return 0;
  }

  if (1 == aj->nbonds()) {    // we start with unsaturated atoms
    return 0;
  }

  return remove_embedded_functional_group(m, j);
}

static int
permutation_remove_embedded_functional_group(Molecule & m,
                                             const Can_Change & can_change,
                                             RandomState& random_state)
{
  if (! random_state.OkProbability(probability_remove_embedded_functional_group)) {
    return 0;
  }

  const int matoms = m.natoms();

  for (int i = 0; i < max_attempts; i++)
  {
    const atom_number_t j = random_state.RandomInt(matoms - 1);

    if (permutation_remove_embedded_functional_group(m, j, can_change))
      return 1;
  }

  return 0;
}

static int
all_terminal_atoms(Molecule & m,
                   const atom_number_t afrom,
                   const atom_number_t astart,
                   Set_of_Atoms & to_be_removed)
{
  const Atom * a = m.atomi(astart);

  for (int i = 0; i < a->ncon(); i++)
  {
    const atom_number_t j = a->other(astart, i);

    if (afrom == j)
      continue;

    if (m.ncon(j) > 1)
      return 0;
  }

// Great, they are all single connections. Mark them to be deleted

  for (int i = 0; i < a->ncon(); i++)
  {
    const atom_number_t j = a->other(astart, i);

    if (afrom != j)
      to_be_removed.add(j);
  }

  return 1;
}

/*
  We have an atom outside a ring. Is it singly connected or something like a CF3, NO2
  that can come off as well
*/

static int
can_also_be_removed(Molecule & m,
                    const Ring & r,
                    const atom_number_t astart,
                    Set_of_Atoms & to_be_removed)
{
  const Atom * a = m.atomi(astart);

  for (int i = 0; i < a->ncon(); i++)
  {
    const atom_number_t j = a->other(astart, i);

    if (r.contains(j)) {
      continue;
    }

    if (1 == m.ncon(j)) {
      to_be_removed.add(j);
      continue;
    }

    if (all_terminal_atoms(m, astart, j, to_be_removed)) {
      to_be_removed.add(j);
      continue;
    }

    return 0;
  }

  return 1;
}

static int
permutation_remove_terminal_ring(Molecule & m,
                     const int ring_number,
                     const Can_Change & can_change)
{
  const Ring & r = *m.ringi(ring_number);

  if (r.is_fused()) {
    return 0;
  }

#ifdef CAN_CHANGE_BROKEN
  if (r.count_members_set_in_array(can_change.rawdata(), 0))
    return 0;
#endif

  int ring_size = r.number_elements();

  Set_of_Atoms to_be_removed;

  for (int i = 0; i < ring_size; i++) {
    to_be_removed.add(r[i]);
  }

  int attachment_points_found = 0;

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];

    const Atom * aj = m.atomi(j);

    if (2 == aj->ncon()) {
      continue;
    }

    if (! can_also_be_removed(m, r, j, to_be_removed)) {
      if (4 == aj->ncon()) {
        return 0;
      }

      attachment_points_found++;
      if (attachment_points_found > 1) {
        return 0;
      }
    }
  }

  if (1 != attachment_points_found)
    return 0;

  m.remove_atoms(to_be_removed);

  if (verbose > 1)
    changed_by << " Rm_Terminal_Ring";
  if (verbose > 2)
    cerr << "Removed " << to_be_removed.number_elements() << " atoms associated with terminal ring with " << ring_size << " atoms\n";

  return 1;
}

static int
permutation_remove_terminal_ring(Molecule & m,
                                 const Can_Change & can_change,
                                 RandomState& random_state) {
  if (! random_state.OkProbability(probability_remove_terminal_ring)) {
    return 0;
  }
  int nr = m.nrings();

  if (0 == nr)
    return 0;

  if (min_rings_in_each_permutation > 0 && (nr - 1) <= min_rings_in_each_permutation)
    return 0;

  int non_fused_rings_found = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);
    if (! ri->is_fused())
    {
      non_fused_rings_found++;
      break;
    }
  }

  if (! non_fused_rings_found)    // napthalene for example
    return 0;

  for (int i = 0; i < max_attempts; i++)
  {
    const int j = random_state.RandomInt(nr - 1);

    if (permutation_remove_terminal_ring(m, j, can_change))
      return 1;
  }

  return 0;
}

static void
randomly_permute(Set_of_Atoms & s)
{
  std::random_shuffle(s.begin(), s.end());
  return;
}

/*
  Used for keeping some index within a range
*/

static int
next_after_wrap(int i,
                int imax)
{
  i++;

  if (i < imax)   // note < and not <=
    return i;

  return 0;
}

#ifdef NOT_USED
static void
set_atoms_to_remove(const Ring & r, 
                    atom_number_t exclude1,
                    atom_number_t exclude2,
                    Set_of_Atoms & to_remove)
{
  to_remove = r;

  to_remove.remove_first(exclude1);
  to_remove.remove_first(exclude2);

  return;
}
#endif

static int
identify_doubly_bonded_atoms_in_ring(const Molecule & m,
                                     const Ring & r,
                                     atom_number_t & a1,
                                     atom_number_t & a2,
                                     RandomState& random_state)
{
  const int n = r.number_elements();

  int istart = random_state.RandomInt(n - 1);

  atom_number_t prev;

  if (0 == istart)
    prev = r.last_item();
  else
    prev = r[istart - 1];

  for (int i = 0; i < n; istart = next_after_wrap(istart, n), i++)
  {
    const atom_number_t a = r[istart];

    const Bond * b = m.bond_between_atoms(prev, a);

    if (2 != m.ncon(prev))    // should be some way to avoid checking atoms twice...
      ;
    else if (2 != m.ncon(a))
      ;
    else if (6 != m.atomic_number(prev))
      ;
    else if (6 != m.atomic_number(a))
      ;
    else if (b->is_double_bond())
    {
      a1 = prev;
      a2 = a;
      return 1;
    }

    prev = a;
  }

  return 0;    // very strange, an aromatic ring with no double bonds
}

static int
permutation_create_fused_ring(Molecule & m,
                              const Ring * ri,
                              const Molecule * parial_ring,
                              const Can_Change & can_change,
                              RandomState& random_state)
{
  m.compute_aromaticity_if_needed();

//#define DEBUG_FORM_FUSED_AROMATIC
#ifdef DEBUG_FORM_FUSED_AROMATIC
  IWString initial_smiles = m.smiles();
#endif

  if (! ri->is_aromatic())
    return 0;

  if (ri->is_fused())   // avoid building large systems
    return 0;

  atom_number_t a1, a2;
  if (! identify_doubly_bonded_atoms_in_ring(m, *ri, a1, a2, random_state))
    return 0;

#ifdef CAN_CHANGE_BROKEN
  if (0 == can_change[a1] || 0 == can_change[a2])
    return 0;
#endif

  int initial_natoms = m.natoms();

  m.add_molecule(parial_ring);

  if (parial_ring->name().contains("DB"))
  {
    m.set_bond_type_between_atoms(a1, a2, SINGLE_BOND);
    m.add_bond(a1, initial_natoms, DOUBLE_BOND);
    m.add_bond(a2, m.natoms() - 1, DOUBLE_BOND);
  }
  else
  {
    m.add_bond(a1, initial_natoms, SINGLE_BOND);
    m.add_bond(a2, m.natoms() - 1, SINGLE_BOND);
  }

#ifdef DEBUG_FORM_FUSED_AROMATIC
  cerr << "MAKE FUSED AROMATIC initial " << initial_smiles << '\n';
  cerr << "                    created " << m.smiles() << '\n';
#endif

#ifdef FREQUENTLY_CHECK_VALENCES
  if (! m.valence_ok())
  {
    cerr << "GACK, formed bad valence '" << m.smiles() << '\n';
    abort();
  }
#endif

  return 1;
}

static int
permutation_create_fused_ring(Molecule & m,
                              const Can_Change & can_change,
                              RandomState& random_state)
{
  if (! random_state.OkProbability(probability_create_fused_ring)) {
    return 0;
  }

  const int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  if (max_rings_in_each_permutation > 0 && nr >= max_rings_in_each_permutation) {
    return 0;
  }

  const Molecule * parial_ring = random_state.RandomItem(aromatic_ring_fusions);

  const int r = random_state.RandomInt(nr - 1);

  return permutation_create_fused_ring(m, m.ringi(r), parial_ring, can_change, random_state);
}

static int
permutation_split_apart_fused_ring(Molecule & m,
                                   const int ring_number,
                                   const Can_Change & can_change,
                                   RandomState& random_state)
{
  const Ring * ri = m.ringi(ring_number);

  if (! ri->is_fused())
    return 0;

  if (ri->largest_number_of_bonds_shared_with_another_ring() > 1)  // too hard
    return 0;

  atom_number_t a1, a2;    // ring fusion atoms

  const Ring * r2 = ri->fused_neighbour(0);

  if (! find_adjacent_atoms_in_common_between_two_rings(*ri, *r2, a1, a2))
    return 0;

  if (m.formal_charge(a1) || m.formal_charge(a2))
    return 0;

  Set_of_Atoms to_remove;  // other atoms in the ring being destroyed

  if (ri->is_aromatic() || r2->is_aromatic())
  {
    const Bond * b = m.bond_between_atoms(a1, a2);

    if (! b->is_double_bond())
      return 0;
  }

  if (random_state() < 0.5)
    to_remove = *ri;
  else
    to_remove = *r2;

  atom_number_t prev = to_remove.last_item();

  for (int i = 0; i < to_remove.number_elements(); i++)
  {
    atom_number_t a = to_remove[i];

    if (a == a1 && prev == a2)   // ring fusion bond
      ;
    else if (a == a2 && prev == a1)   // ring fusion bond
      ;
#ifdef CAN_CHANGE_BROKEN
    else if (0 == can_change[a] || 0 == can_change[prev])
      ;
#endif
    else
    {
      if (a == a1)   // break bond adjacent to ring
        m.remove_bond_between_atoms(a, prev);
      else
        m.set_bond_type_between_atoms(a, prev, SINGLE_BOND);
    }

    prev = a;
  }

  if (! m.valence_ok())
  {
    if (warn_invalid_valence) {
      cerr << "Splitting aromatic ring yields bad valence " << m.smiles() << '\n';
    }
  }

#ifdef DEBUG_SPLIT_FUSED_RINGS
  cerr << "Split apart fused rings, initial '" << initial_smiles << "' now '" << m.smiles() << "'\n";
#endif

  return 1;
}

static int
permutation_split_apart_fused_ring(Molecule & m,
                                   const Can_Change & can_change,
                                   RandomState& random_state)
{
  if (! random_state.OkProbability(probability_split_fused_ring)) {
    return 0;
  }

  int nr = m.nrings();

  if (nr < 2)
    return 0;

  if (min_rings_in_each_permutation > 0 && nr <= min_rings_in_each_permutation)
    return 0;

#ifdef DEBUG_SPLIT_FUSED_RINGS
  IWString initial_smiles = m.smiles();
#endif

  m.compute_aromaticity_if_needed();

  int istart = random_state.RandomInt(nr - 1);

  for (int i = 0; i < nr; istart = next_after_wrap(istart, nr), i++)
  {
    if (permutation_split_apart_fused_ring(m, istart, can_change, random_state))
      return 1;
  }

  return 0;
}

static int
switch_ring(Molecule & m,
            const Set_of_Atoms & ring_to_remove_arg,
            Molecule & ring_to_add)
{
  Set_of_Atoms ring_to_remove = ring_to_remove_arg;    // we need a copy because add_molecule destroys existing rings

  Set_of_Atoms attachment_points_for_old_ring;

  for (atom_number_t j : ring_to_remove) {
    const int jcon = m.ncon(j);
    
    if (4 == jcon) {     // can't do these right now
      return 0;
    }

    if (2 != jcon) {
      attachment_points_for_old_ring.add(j);
    }
  }

  Set_of_Atoms attachment_points_open;   // must be in a ring and have an implicit Hydrogen

  int ratoms = ring_to_add.natoms();
  for (int i = 0; i < ratoms; i++)
  {
    if (! ring_to_add.is_ring_atom(i)) {
      continue;
    }

    if (ring_to_add.hcount(i)) {
      attachment_points_open.add(i);
    }
  }

  if (attachment_points_open.size() < attachment_points_for_old_ring.size()) {
    return 0;
  }

  randomly_permute(attachment_points_open);

// Now the tricky part. We want to break all the bonds to the old ring, and re-make them to atoms
// in the new ring. Note that any previous double bond will be converted to a single bond

  const int previous_natoms = m.natoms();

  m.add_molecule(&ring_to_add);

// We need an index into the attachment_points_open array

  int apo = 0;

  for (atom_number_t j : ring_to_remove) {
    const Atom * aj = m.atomi(j);

    for (int k = 0; k < aj->ncon(); k++)
    {
      atom_number_t l = aj->other(j, k);
      if (ring_to_remove.contains(l))
        continue;

      m.remove_bond_between_atoms(j, l);
      m.add_bond(l, attachment_points_open[apo] + previous_natoms, SINGLE_BOND);
      apo++;
    }
  }

  m.remove_atoms(ring_to_remove);

  return 1;
}

static int
permutation_switch_ring(Molecule & m,
                        int aromatic_rings_only,
                        Molecule & ring_to_add,
                        RandomState& random_state)
{
  int nr = m.nrings();

  if (0 == nr) {
    return 0;
  }

  m.compute_aromaticity_if_needed();

  resizable_array<const Ring *> rings_to_try;

  for (int i = 0; i < nr; i++) {
    const Ring * ri = m.ringi(i);
    if (ri->is_fused()) {
      continue;
    }
    
    int arom = ri->is_aromatic();

    if ((arom && aromatic_rings_only) || (! arom && ! aromatic_rings_only))
    {
      rings_to_try.add(ri);
    }
  }

//cerr << "There are " << rings_to_try.number_elements() << " rings to try\n";

  if (rings_to_try.empty()) {
    return 0;
  }

  for (int i = 0; i < max_attempts; i++)
  {
    const Ring * rj = random_state.RandomItem(rings_to_try);

    if (switch_ring(m, *rj, ring_to_add))
      return 1;
  }

  return 0;
}

static int
permutation_switch_aromatic_ring(Molecule & m,
                                 const Can_Change & can_change,
                                 RandomState& random_state)
{
//cerr << "Can we switch aromatic rings " << probability_switch_aromatic_ring << '\n';
  if (! random_state.OkProbability(probability_switch_aromatic_ring)) {
    return 0;
  }

  Molecule * f = random_state.RandomItem(aromatic_rings);

  return permutation_switch_ring(m, 1, *f, random_state);
}

static int
permutation_switch_aliphatic_ring(Molecule & m,
                                  const Can_Change & can_change,
                                  RandomState& random_state)
{
//cerr << "Can we switch aliphatic rings " << probability_switch_aliphatic_ring << '\n';
  if (! random_state.OkProbability(probability_switch_aliphatic_ring)) {
    return 0;
  }

  Molecule * f = random_state.RandomItem(aliphatic_rings);

  return permutation_switch_ring(m, 0, *f, random_state);
}

static int
find_random_single_bond(Molecule & m,
                        const Can_Change & can_change,
                        RandomState& random_state)
{
  const int n = m.nedges();
  const int i = random_state.RandomInt(n - 1);

  for(int j=0;j<n;j++)
  {
    const int k = next_after_wrap(i, n);

    const Bond * b=m.bondi(k);

    if (! b->is_single_bond())
      continue;

    if (b->nrings())
      continue;

#ifdef CAN_CHANGE_BROKEN
    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();

    if (0 == can_change[a1] || 0 == can_change[a2])
      continue;
#endif

    return k;
  }

  return -1;
}

int
Molecular_Transformation_Add_Double_Bond::random_transformation(Molecule & m, const Can_Change & can_change)
{
  return permutation_add_double_bond(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Atom::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_atom(m, can_change, _random_state);
}

int
Molecular_Transformation_Lower_Bond_Order::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_lower_bond_order(m, can_change, _random_state);
}

int
Molecular_Transformation_Break_Aliphatic_Ring::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_break_aliphatic_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Make_Ring::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_make_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_Carbon::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_carbon(m, can_change, _random_state);
}

int
Molecular_Transformation_Move_Fragment::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_move_fragment(m, can_change, _random_state);
}

int
Molecular_Transformation_Swap_Adjacent_Atoms::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_swap_adjacent_atoms(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Carbon_to_Nitrogen::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_carbon_to_nitrogen(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Nitrogen_to_Carbon::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_nitrogen_to_carbon(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Carbon_to_Oxygen::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_carbon_to_oxygen(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Oxygen_to_Something::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_oxygen_to_something(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_From_Single_Attachment_Library::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_from_single_attachment_library(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_From_Aromatic_Attachment_Library::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_from_aromatic_attachment_library(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_From_Double_Attachment_Library::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_from_double_attachment_library(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Ring_Substitution::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_ring_substitution(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Terminal_Functional_Group::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_terminal_functional_group(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Embedded_Functional_Group::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_embedded_functional_group(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Terminal_Ring::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_terminal_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Switch_Aromatic_Ring::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_switch_aromatic_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Switch_Aliphatic_Ring::random_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_switch_aliphatic_ring(m, can_change, _random_state);
}



int
Molecular_Transformation_Add_Double_Bond::next_transformation(Molecule & m, const Can_Change & can_change)
{
  return permutation_add_double_bond(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Atom::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_atom(m, can_change, _random_state);
}

int
Molecular_Transformation_Lower_Bond_Order::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_lower_bond_order(m, can_change, _random_state);
}

int
Molecular_Transformation_Break_Aliphatic_Ring::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_break_aliphatic_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Make_Ring::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_make_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_Carbon::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_carbon(m, can_change, _random_state);
}

int
Molecular_Transformation_Move_Fragment::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_move_fragment(m, can_change, _random_state);
}

int
Molecular_Transformation_Swap_Adjacent_Atoms::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_swap_adjacent_atoms(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Carbon_to_Nitrogen::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_carbon_to_nitrogen(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Nitrogen_to_Carbon::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_nitrogen_to_carbon(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Carbon_to_Oxygen::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_carbon_to_oxygen(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Oxygen_to_Something::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_oxygen_to_something(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_From_Single_Attachment_Library::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_from_single_attachment_library(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_From_Aromatic_Attachment_Library::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_from_aromatic_attachment_library(m, can_change, _random_state);
}

int
Molecular_Transformation_Add_From_Double_Attachment_Library::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_add_from_double_attachment_library(m, can_change, _random_state);
}

int
Molecular_Transformation_Change_Ring_Substitution::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_change_ring_substitution(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Terminal_Functional_Group::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_terminal_functional_group(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Embedded_Functional_Group::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_embedded_functional_group(m, can_change, _random_state);
}

int
Molecular_Transformation_Remove_Terminal_Ring::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_remove_terminal_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Switch_Aromatic_Ring::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_switch_aromatic_ring(m, can_change, _random_state);
}

int
Molecular_Transformation_Switch_Aliphatic_Ring::next_transformation (Molecule & m, const Can_Change & can_change)
{
  return permutation_switch_aliphatic_ring(m, can_change, _random_state);
}





static int
breed(Molecule & m,
      const Can_Change & can_change,
      Molecule & f,
      RandomState& random_state)
{
#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
  assert (f.valence_ok());
#endif

  m.ring_membership();

  const int i=find_random_single_bond(m, can_change, random_state);
  if (i < 0)
    return 0;

  int j=find_random_single_bond(f, can_change, random_state);
  if (j < 0)
    return 0;

  int initial_matoms = m.natoms();

  IWString smiles1 = m.smiles();
  IWString smiles2 = f.smiles();

  j+=m.nedges();
  m.add_molecule(&f);
  const Bond *bi=m.bondi(i);
  const Bond *bj=m.bondi(j);
  int a1=bi->a1();
  int a2=bi->a2();
  int a3=bj->a1();
  int a4=bj->a2();
  m.remove_bond(j); //firstly remove higher bond index j, so i will not change
  m.remove_bond(i);

  if (1 == m.number_fragments())
  {
    cerr << "Huh, one fragment after removing bonds\n";
    m.resize(initial_matoms);
    m.add_bond(a1, a2, SINGLE_BOND);
    cerr << m.smiles() << " and " << f.smiles() << '\n';
    return 0;
  }

  int n1=m.atoms_in_fragment(m.fragment_membership(a1));
  int n2=m.atoms_in_fragment(m.fragment_membership(a2));
  int n3=m.atoms_in_fragment(m.fragment_membership(a3));
  int n4=m.atoms_in_fragment(m.fragment_membership(a4));
  if((n1+n3)>=lower_atom_count_cutoff && (n1+n3)<=upper_atom_count_cutoff
          && (6==m.atomic_number(a1) || 6==m.atomic_number(a3)))
  {
    m.add_bond(a1,a3,SINGLE_BOND);
    m.add_bond(a2,a4,SINGLE_BOND);
    m.remove_fragment_containing_atom(a2);
  } 
  else if((n1+n4)>=lower_atom_count_cutoff && (n1+n4)<=upper_atom_count_cutoff
          && (6==m.atomic_number(a1) || 6==m.atomic_number(a4)))
  {
    m.add_bond(a1,a4,SINGLE_BOND);
    m.add_bond(a2,a3,SINGLE_BOND);
    m.remove_fragment_containing_atom(a2);
  }
  else if((n2+n4)>=lower_atom_count_cutoff && (n2+n4)<=upper_atom_count_cutoff
          && (6==m.atomic_number(a2) || 6==m.atomic_number(a4)))
  {
    m.add_bond(a2,a4,SINGLE_BOND);
    m.add_bond(a1,a3,SINGLE_BOND);
    m.remove_fragment_containing_atom(a1);
  }
  else if((n2+n3)>=lower_atom_count_cutoff && (n2+n3)<=upper_atom_count_cutoff
          && (6==m.atomic_number(a2) || 6==m.atomic_number(a3)))
  {
    m.add_bond(a2,a3,SINGLE_BOND);
    m.add_bond(a1,a4,SINGLE_BOND);
    m.remove_fragment_containing_atom(a1);
  }
  else
  {
    m.resize(initial_matoms);
    m.add_bond(a1, a2, SINGLE_BOND);
    return 0;
  }

#ifdef FREQUENTLY_CHECK_VALENCES
  if (! m.valence_ok())
  {
    cerr << "Invalid valence in breed '" << m.smiles() << "'\n";
    cerr << smiles1 << '\n';
    cerr << smiles2 << '\n';
    abort();
  }
#endif

  return 1;
}

static int
permutation_breed(Molecule & m,
                  const Can_Change & can_change,
                  RandomState& random_state)
{
  if (1.0 == probability_breed)
    ;
  else if (! random_state.OkProbability(probability_breed)) 
    return 0;

  // DO not select the last molecule.
  int j = random_state.RandomInt(input_molecules.number_elements() - 1);

  Molecule *f =input_molecules[j];

  return breed(m, can_change, *f, random_state);
}


/*
  We have an array of transformation objects.

  It is very important that the inter molecular permutation(s) be last in the array
  Right now, we have just one inter molecular permutation
*/

#define NUMBER_KNOWN_INTRA_PERMUTATIONS 24

#define NUMBER_KNOWN_INTER_PERMUTATIONS 2

#define NUMBER_KNOWN_PERMUTATIONS (NUMBER_KNOWN_INTRA_PERMUTATIONS + NUMBER_KNOWN_INTER_PERMUTATIONS)

static int number_known_permutations = NUMBER_KNOWN_PERMUTATIONS;

/*
  We have an array of functions that all have the same signature. We
  call one at random
*/

static int (*known_permutation[NUMBER_KNOWN_PERMUTATIONS]) (Molecule &, const Can_Change &, RandomState& random_state);

/*
  We make NSTART copies of the molecule and make NUMBER_PERTURBATIONS to each
  of those copies
*/

static int nstart = 1000;

static int number_perturbations = 5;

static extending_resizable_array<int> molecules_made_per_copy;

/*
  We don't know if we will be keeping track of duplicates just for each
  starting molecule, or across all molecules
*/

static int discard_duplicates_within_each_molecule = 0;
static int discard_duplicates_across_all_molecules = 0;

static IW_STL_Hash_Map_int usmi;

static Set_of_Target_Molecules<GFP_Standard> set_of_target_molecules;

/*
  We want to create some molecules that appear close to our targets.
*/

class Molecule_and_Distance : public Molecule
{
  private:
    float _dist;

  public:
    Molecule_and_Distance(const Molecule & m, const float s);

    void set_dist(const float s) { _dist = s;}
    float dist() const { return _dist;}
};

Molecule_and_Distance::Molecule_and_Distance(const Molecule & m, const float d) : Molecule(m), _dist(d)
{
  return;
}

class Molecule_and_Distance_Comparator
{
  private:
  public:
    bool operator() (const float d, const Molecule_and_Distance * m1) const 
//  bool operator() (const Molecule_and_Distance * m1, const float d) const 
    {
//    return (*m1)->dist() < (*m2)->dist();
      return m1->dist() > d;
    };
};

class Set_of_Closer_Molecules
{
  private:
    // Seems this should be a resizable_array_p<Molecule_and_Distance>
    Molecule_and_Distance ** _m;
    int _n;
    int _capacity;

    float _max_distance;

    Molecule_and_Distance_Comparator _cmp;

    IW_STL_Hash_Set _smiles;

//  private functions

    int _check_sorted() const;

  public:
    Set_of_Closer_Molecules();
    ~Set_of_Closer_Molecules();

    int initialise(const int s);

    int number_molecules() const { return _n;}

    Molecule * moleculei(const int s) { return _m[s];}

    void set_max_distance(const float s) { _max_distance = s;}

    Molecule_and_Distance * a_random_molecule(RandomState& random_state) const;

    int would_retain(const float s) const;

    int consider(Molecule & m, const float d);
};

Set_of_Closer_Molecules::Set_of_Closer_Molecules()
{
  _m = nullptr;
  _n = 0;
  _capacity = 0;
  _max_distance = 1.0f;
}

Set_of_Closer_Molecules::~Set_of_Closer_Molecules()
{
  if (nullptr != _m)
    delete [] _m;
}

int
Set_of_Closer_Molecules::initialise(const int s)
{
  _capacity = s;

  _n = 0;
  _m = new Molecule_and_Distance *[_capacity];

  return 1;
}

int
Set_of_Closer_Molecules::would_retain(const float s) const
{
  if (0 == _n) {
    return 0;
  }

  return s < _m[_n-1]->dist();
}

int
Set_of_Closer_Molecules::consider(Molecule & m, const float d)
{
  if (d > _max_distance) {
    return 0;
  }

  if (_n < _capacity)   // have space, will store it
    ;
  else if (d >= _m[_n-1]->dist())    // longer than our farthest distance, do not want
    return 0;

  if (_smiles.contains(m.unique_smiles())) {
    return 0;
  }

  _smiles.emplace(m.unique_smiles());

  Molecule_and_Distance * mad = new Molecule_and_Distance(m, d);
  mad->set_name(m.name());

  if (0 == _n)
  {
    _m[0] = mad;
    _n = 1;
    return 1;
  }

  const auto f = std::upper_bound(_m, _m + _n, d, _cmp);

//cerr << "Looking for " << d << " in " << _n << " items, iterator " << (f - _m) << " short " << _m[0]->dist() << " long " << _m[_n-1]->dist() << '\n';
  if (f == _m + _n)   // nothing in the list is shorter
  {
    _m[_n] = mad;
    _n++;
    return 1;
  }

// list is full, we must insert somewhere

  if (_n == _capacity) {
    delete _m[_n-1];
  }

  const auto ndx = f - _m;

  for (int i = _n; i > ndx; --i)
  {
    _m[i] = _m[i-1];
  }

  _m[ndx] = mad;

  _check_sorted();

  if (_n < _capacity) {
    _n++;
  }

  return 1;
}

int
Set_of_Closer_Molecules::_check_sorted() const
{
  for (int i = 1; i < _n; ++i) {
    if (_m[i-1]->dist() > _m[i]->dist()) {
      cerr << "Set_of_Closer_Molecules::_check_sorted:out of order, i = " << i << " prev " << _m[i-1]->dist() << " curr " << _m[i]->dist() << '\n';
      return 0;
    }
  }

  return 1;
}

Molecule_and_Distance *
Set_of_Closer_Molecules::a_random_molecule(RandomState& random_state) const
{
  if (_n == 0) {
    return nullptr;
  }

  int j = random_state.RandomInt(_n - 1);
  return _m[j];
}

static Set_of_Closer_Molecules closer_molecules;

static IWString_and_File_Descriptor stream_for_more_proximal_molecules;

/*
  Since the target molecules set may add molecules to the input set, we impose a maximum
*/

static Accumulator<double> acc_dist_to_target;

/*
  Once we have read our input molecules, we can process any newly created molecules
*/

static int iterations_through_new_molecules = 0;

static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << "  -x <number>    number of copies of the starting molecule to use (" << nstart << ")\n";
  cerr << "  -p <perm>      number of permutations to make on each molecule (default " << number_perturbations << ")\n";
  cerr << "  -F query       never change fragments labelled with isotope <isotope>\n";
  cerr << "  -f <isotope>   use <isotope> to label the atoms not to be changed\n";
  cerr << "  -P <file>      file containing transformation probabilities\n";
  cerr << "  -r <size>      never form a ring with <size> atoms or smaller\n";
  cerr << "  -R <size>      never form a ring with <size> atoms or larger\n";
  cerr << "  -c <number>    exclude molecules with atom count below <number>\n";
  cerr << "  -C <number>    exclude molecules with atom count above <number>\n";
  cerr << "  -y <number>    minimum number of rings in each molecule\n";
  cerr << "  -Y <number>    maximum number of rings in each molecule\n";
  cerr << "  -G +<number>   maximum number of rings that can be added to each molecule\n";
  cerr << "  -G -<number>   maximum number of rings that can be removed from each molecule\n";
  cerr << "  -L single=<file>  library file of singly attached fragments - join at first atom\n";
  cerr << "  -L singleR=<file> library of singly attached fragments - join from any atom\n";
  cerr << "  -L double=<file>  library file of doubly attached fragments\n";
  cerr << "  -L arom=<file>    library file of fragments that are added only to aromatic rings\n";
  cerr << "  -L aromR=<file>   library file of aromatic rings\n";
  cerr << "  -L aliphR=<file>  library file of aliphatic rings\n";
  cerr << "  -L fsarom=<file>  library of aromatic rings to fuse (partial rings input)\n";
  cerr << "  -q <queries>   queries which must match - non matches are discarded\n";
  cerr << "  -Q <queries>   queries which must be avoided - matches are discarded\n";
  cerr << "  -K once        write the initial molecule to the output\n";
  cerr << "  -K copy        write the initial molecule to the output before each copy (-s)\n";
  cerr << "  -K every       write the initial molecule to the output before every perturbation\n";
  cerr << "  -J <name>      write discarded molecules to <name>\n";
  cerr << "  -U each        discard duplicate molecules - check molecules separately\n";
  cerr << "  -U all         global checking for duplicate molecules\n";
  cerr << "  -M ...         miscellaneous fine control options\n";
  cerr << "  -T <fname>     file of target molecules to try to make, enter -T help for info\n";
  //(void) display_standard_aromaticity_options(cerr);
  //cerr << "  -E <symbol>    create an element with symbol <symbol>\n";
  //cerr << "  -E autocreate  automatically create new elements when encountered\n";
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int
permutation_breed_with_close_to_target(Molecule & m,
                                       const Can_Change & can_change,
                                       RandomState& random_state)
{
  if (! random_state.OkProbability(probability_breed_close_to_target)) {
    return 0;
  }

  Molecule_and_Distance * mad = closer_molecules.a_random_molecule(random_state);

  if (nullptr == mad)
    return 0;

  return breed(m, can_change, *mad, random_state);
}

static void
display_miscellaneous_options(const char c,
                              std::ostream & output)
{
  output << "  -" << c << " afras      allow fused ring atoms to swap atoms - can generate strange ring systems\n";
  output << "  -" << c << " cblast     the changed by string contains only the last permuation\n";
  output << "  -" << c << " okcage     will allow formation of strongly fused ring systems\n";
  output << "  -" << c << " prd=<n>    produce <n> molecules then stop\n";
  output << "  -" << c << " run=<secs> run for <secs> seconds and terminate\n";
  output << "  -" << c << " stopfile=<fname> terminate processing if <fname> present\n";
  output << "  -" << c << " scrand     all sidechains  will be attached via random atoms\n";
  output << "  -" << c << " ssintra    perform non random sequential scan of all intra molecular transformations\n";
  output << "  -" << c << " all        assign unit probability to every transformation\n";

  exit(1);
}

/*
  when doing exhaustive enumeration we need a class to hold the iterator and the function
*/

template <typename T, typename I>
class Transformation_and_Iterator
{
  protected:
    T _transform;
    I _iterator;

  public:
    Transformation_and_Iterator(T & t);

    virtual int process(Molecule & m, const Can_Change & can_change) = 0;

    virtual void reset() { _iterator.reset();}
};

template <typename T>
class Transformation_and_Atom_Iterator : public Transformation_and_Iterator<T, Iterate_Range>
{
  private:
    using Transformation_and_Iterator<T, Iterate_Range>::_transform;
  public:
    Transformation_and_Atom_Iterator(const Molecule & m, T & t);

    int process(Molecule & m, const Can_Change & can_change);
    void process();
};

template <typename T>
Transformation_and_Atom_Iterator<T>::Transformation_and_Atom_Iterator(const Molecule & m, T & t) :
                                                Transformation_and_Iterator<T, Iterate_Range>(t), 
                                                Transformation_and_Iterator<T, Iterate_Range>::_iterator(m.natoms())
{
  return;
}

template <typename T>
int
Transformation_and_Atom_Iterator<T>::process(Molecule & m, const Can_Change & can_change)
{
  int atom_number;

  if (! Transformation_and_Iterator<T, Iterate_Range>::_iterator.next(atom_number))
    return 0;

  return _transform->process(m, atom_number, can_change);
}

template <typename T>
class Transformation_and_Bond_Iterator : public Transformation_and_Iterator<T, Iterate_Range>
{
  private:
    using Transformation_and_Iterator<T, Iterate_Range>::_transform;
  public:
    Transformation_and_Bond_Iterator(const Molecule & m, T & t);

    int process(Molecule & m, const Can_Change & can_change);
};

template <typename T>
Transformation_and_Bond_Iterator<T>::Transformation_and_Bond_Iterator(const Molecule & m, T & t) :
                                                Transformation_and_Iterator<T, Iterate_Range>(t), 
                                                Transformation_and_Iterator<T, Iterate_Range>::_iterator(m.nedges())
{
  return;
}

template <typename T>
int
Transformation_and_Bond_Iterator<T>::process(Molecule & m, const Can_Change & can_change)
{
  int bond_number;

  if (! Transformation_and_Iterator<T, Iterate_Range>::_iterator.next(bond_number))
    return 0;

  return _transform->process(m, bond_number, can_change);
}

template <typename T>
class Transformation_and_Ring_Iterator : public Transformation_and_Iterator<T, Iterate_Range>
{
  private:
    using Transformation_and_Iterator<T, Iterate_Range>::_transform;

  public:
    Transformation_and_Ring_Iterator(Molecule & m, T & t);

    int process(Molecule & m, const Can_Change & can_change);
};

template <typename T>
Transformation_and_Ring_Iterator<T>::Transformation_and_Ring_Iterator(Molecule & m, T & t) :
                                                Transformation_and_Iterator<T, Iterate_Range>(t), 
                                                Transformation_and_Iterator<T, Iterate_Range>::_iterator(m.nrings())
{
  return;
}

template <typename T>
int
Transformation_and_Ring_Iterator<T>::process(Molecule & m, const Can_Change & can_change)
{
  int ring_number;

  if (! Transformation_and_Iterator<T, Iterate_Range>::_iterator.next(ring_number))
    return 0;

  return _transform->process(m, ring_number, can_change);
}

template <typename T>
class
Transformation_Make_Ring : public Transformation_and_Iterator<T, Iterate_Range>
{
  private:
    using Transformation_and_Iterator<T, Iterate_Range>::_transform;
    using Transformation_and_Iterator<T, Iterate_Range>::_iterator;

    Iterate_Range _a2;

  public:
    Transformation_Make_Ring(T & t, Molecule & m);

    int process(Molecule & m, const Can_Change & can_change);

    void reset();
};

template <typename T>
Transformation_Make_Ring<T>::Transformation_Make_Ring(T & t, Molecule & m) : 
                                          Transformation_and_Iterator<T, Iterate_Range>::_transform(Transformation_and_Iterator<T, Iterate_Range>(t)),
                                          Transformation_and_Iterator<T, Iterate_Range>::_iterator(Transformation_and_Iterator<T, Iterate_Range>(m.natoms())),
                                          _a2(m.natoms())
{
  return;
}

template <typename T>
int
Transformation_Make_Ring<T>::process(Molecule & m, const Can_Change & can_change)
{
  int a1;
  while (_transform.next(a1))
  {
    int a2;
    while (_a2.next(a2))
    {
      if (_transform->process(m, a1, a2, can_change))
        return 1;
    }

    _a2.reset();
  }

  return 0;
}

template <typename T>
class
Transformation_Create_Fused_Ring : public Transformation_and_Iterator<T, Iterate_Range>
{
  private:
    using Transformation_and_Iterator<T, Iterate_Range>::_transform;
    using Transformation_and_Iterator<T, Iterate_Range>::_iterator;

    resizable_array_p<Molecule> & _aromatic_ring_fusions;
    Iterate_Range _lib_ndx;

  public:
    Transformation_Create_Fused_Ring(T & t, Molecule & m, resizable_array_p<Molecule> & aromatic_ring_fusions);

    int process(Molecule & m, const Can_Change & can_change);
};

template <typename T>
Transformation_Create_Fused_Ring<T>::Transformation_Create_Fused_Ring(T & t, Molecule & m,
                                        resizable_array_p<Molecule> & f) : 
                                          Transformation_and_Iterator<T, Iterate_Range>::_transform(Transformation_and_Iterator<T, Iterate_Range>(t)),
                                          Transformation_and_Iterator<T, Iterate_Range>::_iterator(Transformation_and_Iterator<T, Iterate_Range>(m.nrings())),
                                          _aromatic_ring_fusions(f)
{
}

template <typename T>
int
Transformation_Create_Fused_Ring<T>::process(Molecule & m,
                                       const Can_Change & can_change)
{
  int l;
  while (_lib_ndx.next(l))
  {
    int r;
    while (_iterator.next(r))
    {
      if (_transform->process(m, m.ringi(r), _aromatic_ring_fusions[l], can_change))
        return 1;
    }

    _iterator.reset();
  }
}

template <typename T>
class Transformation_Add_Lib_Member : public Transformation_and_Iterator<T, Iterate_Range>
{
  private:
    using Transformation_and_Iterator<T, Iterate_Range>::_transform;
    using Transformation_and_Iterator<T, Iterate_Range>::_iterator;

    resizable_array_p<Molecule> & _lib;
    Iterate_Range _iter_lib;

    const bond_type_t _bt;

  public:
    Transformation_Add_Lib_Member(T & t, Molecule &, resizable_array_p<Molecule> &, const bond_type_t);

    int process(Molecule & m, const Can_Change & can_change);

    void reset();
};

template <typename T>
Transformation_Add_Lib_Member<T>::Transformation_Add_Lib_Member(T & t, Molecule & m,
                                        resizable_array_p<Molecule> & lib,
                                        const bond_type_t bt) :
                                            Transformation_and_Iterator<T, Iterate_Range>::_transform(Transformation_and_Iterator<T, Iterate_Range>(t)),
                                            Transformation_and_Iterator<T, Iterate_Range>::_iterator(Transformation_and_Iterator<T, Iterate_Range>(m.nrings())),
                                            _lib(lib),
                                            _bt(bt)
{
}

template <typename T>
int
Transformation_Add_Lib_Member<T>::process(Molecule & m, const Can_Change & can_change)
{
  int l;
  while (_iter_lib.next(l))
  {
    int a;
    while (_iterator.next(a))
    {
      if (_transform(m, a, _lib[l], can_change))
        return 1;
    }

    _iterator.reset();
  }
}

template <typename T>
void
Transformation_Add_Lib_Member<T>::reset()
{
  _iter_lib.reset();

  Transformation_and_Iterator<T, Iterate_Range>::reset();

  return;
}

/*
  main entry point for exhaustive application of rules
  
  We have some number of permutations that are driven by an atom number.
  Some that are driven by a pair of atoms.
  Some that are driven by a ring number
*/

#ifdef TODO_IMPLEMENT_THIS_SOMETIME
static int
exhaustive_application(Molecule & m,
                       const Can_Change & can_change,
                       const int min_rings_this_molecule,
                       const int max_rings_this_molecule,
                       const int depth)
{

  for (int i = 0; i < NUMBER_KNOWN_PERMUTATIONS; ++i)
  {
  }

  return 1;
}
#endif

/*
  If the molecule looks closer to our target, we can add it to the global INPUT_MOLECULES array
*/

static int
check_proximity_to_target_set(Molecule & m,
                              const int perturbation_number,
                              Set_of_Target_Molecules<GFP_Standard>& set_of_target_molecules)
{
  Accumulator<float> acc;
  set_of_target_molecules.closest_distance(m, acc);

  const float d = acc.minval();

  acc_dist_to_target.extra(d);

  if (! closer_molecules.consider(m, d)) {
    return 0;
  }

  if (stream_for_more_proximal_molecules.is_open()) {
    stream_for_more_proximal_molecules << m.smiles() << ' ' << m.name() << ' ' << perturbation_number << ' ' << d << '\n';
    stream_for_more_proximal_molecules.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

/*
  We need two different kinds of number generators.
  One is a random one that returns random numbers in a range.
  the other is one that will cycle through the range.
*/

class Stream_of_Random_Numbers
{
  private:
    std::uniform_int_distribution<int> _u;
    std::mt19937_64 _rng;

  public:
    Stream_of_Random_Numbers(const int);
    Stream_of_Random_Numbers(const int, const unsigned int);

    int operator()();
};

Stream_of_Random_Numbers::Stream_of_Random_Numbers(const int imax) : _u(0, imax), _rng(rd())
{
}


// Only useful for test suite, do not use otherwise
Stream_of_Random_Numbers::Stream_of_Random_Numbers(const int imax, const unsigned int seed) : _u(0, imax), _rng(seed)
{
}

int
Stream_of_Random_Numbers::operator()()
{
  return _u(_rng);
}

class Stream_of_Sequential_Numbers
{
  private:
    const int _nmax;
    int _ndx;

  public:
    Stream_of_Sequential_Numbers(const int);

    int operator()();
};

Stream_of_Sequential_Numbers::Stream_of_Sequential_Numbers(const int imax) : _nmax(imax)
{
  _ndx = -1;

  return;
}

int
Stream_of_Sequential_Numbers::operator()()
{
  _ndx++;
  if (_ndx <= _nmax)
    return _ndx;

  _ndx = 0;

  return _ndx;
}

template <typename C>
int
do_a_single_random_permutation(Molecule & m,
                               const Can_Change & can_change,
                               C & transformation_chooser,
                               RandomState& random_state)
{
#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  if (m.empty()) {
    cerr << "Yipes, empty structure, cannot continue\n";
    return 0;
  }

  assert (m.ok());

  int attempts = 0;

  while (1)
  {
    int p =0;
    if (random_state.OkProbability(probability_inter_permutation)) {
      p=random_state.RandomInt(NUMBER_KNOWN_INTRA_PERMUTATIONS, number_known_permutations - 1);
    }
    else {     // an intra molecular transformation
      p = transformation_chooser();
    } 

    if (verbose > 2)
      cerr << "Trying permutation " << p << '\n';

    int (*f) (Molecule &, const Can_Change &, RandomState&) = known_permutation[p];

    if (f(m, can_change, random_state))
    {
#ifdef FREQUENTLY_CHECK_VALENCES
      assert (m.valence_ok());
#endif
      return 1;
    }

    attempts++;
#ifdef FREQUENTLY_CHECK_VALENCES
    assert (m.valence_ok());
#endif

    if (attempts > max_attempts) {
      return 0;
    }
  }

  return 1;
}

static int
do_the_output(Molecule & m,
              const int molecules_made,
              const int which_copy,
              const int which_permutation,
              Molecule_Output_Object & output)
{
  const int matoms = m.natoms();
  if (matoms >= lower_atom_count_cutoff && matoms <= upper_atom_count_cutoff)
    ;
  else
  {
    if (stream_for_discarded_molecules.active()) {
      return stream_for_discarded_molecules.write(m);
    }

    return 1;
  }

  atoms_in_molecules_written.extra(matoms);

  const IWString mname(m.name());

  if (number_assigner.active()) {
    number_assigner.process(m);

    IWString tmp(m.name());
    if (append_permutation_number_to_name)
      tmp << ' ' << which_copy << ',' << which_permutation;
    m.set_name(tmp);
    output.write(m);

    m.set_name(mname);

    return 1;
  }

  IWString tmp;

  tmp << mname << ' ' << molecules_made << " molecule " << molecules_read << " copy " << which_copy << " perm " << which_permutation;

  if (0 == changed_by.length())
    ;
  else if (changed_by_contains_full_history)
    tmp << changed_by;
  else
  {
    const_IWSubstring tmp2(changed_by);

    tmp2.remove_leading_words(changed_by.nwords() - 1);
    tmp << tmp2;
  }

  m.set_name(tmp);

  int rc = output.write(m);

  m.set_name(mname);

  ++molecules_written;

  return rc;
}

static int
violates_uniqueness(Molecule & m)
{
  if (0 == discard_duplicates_across_all_molecules &&
      0 == discard_duplicates_within_each_molecule)
    return 0;

  IWString smi = m.unique_smiles();

  if (! usmi.contains(smi)) {
    //usmi[smi] = 1;
    usmi.emplace(std::make_pair(smi, 1));
    return 0;
  }

  usmi[smi]++;

  return 1;
}

static int
_violates_query_constraints(Molecule & m,
                            const int min_rings_this_molecule,
                            const int max_rings_this_molecule)
{
  const int nr = m.nrings();
#ifdef DEBUG_VIOLATES_QUERY_CONSTRAINTS
  cerr << nr << " rings, min_rings_this_molecule " << min_rings_this_molecule << " max_rings_this_molecule " << max_rings_this_molecule << '\n';
#endif

  if (min_rings_this_molecule > 0 && nr <= min_rings_this_molecule)
  {
    molecules_with_too_few_rings++;
    return 1;
  }

  if (max_rings_this_molecule > 0 && nr >= max_rings_this_molecule)
  {
    molecules_with_too_many_rings++;
    return 1;
  }

  if (! m.valence_ok())
  {
    if (warn_invalid_valence) {
      cerr << "Warning invalid valence encountered " << molecules_made << '\n';
      cerr << m.smiles() << ' ' << m.name() << '\n';
    }
    return 1;
  }

// Need to drive out peroxides - remove when problem found

//#define KILL_PEROXIDES
#ifdef KILL_PEROXIDES
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    const Atom * ai = m.atomi(i);
    if (8 != ai->atomic_number())
      continue;

    for (int j = 0; j < ai->ncon(); j++)
    {
      atom_number_t k = ai->other(i, j);
      if (8 == m.atomic_number(k))
      {
        cerr << "Found peroxide\n";
        abort();
      }
    }
  }
#endif

  if (minimum_ring_size_allowed > 0 || maximum_ring_size_allowed > 0)
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = m.ringi(i);

      if (ri->number_elements() <= minimum_ring_size_allowed)
        return 1;

      if (maximum_ring_size_allowed > 0 && ri->number_elements() >= maximum_ring_size_allowed)
        return 1;
    }
  }

  if (queries_to_match.empty() && queries_to_avoid.empty()) {
    return 0;
  }

  Molecule_to_Match target(&m);

  // If multiple must match queries are present, all must be matched.
  // A non-match means we have violated the constraint.
  // That seems harsh, an OR match would seem better.
  for (Substructure_Hit_Statistics * q : queries_to_match) {
    if (0 == q->substructure_search(target)) {
      return 1;
    }
  }

  // If any avoid query is present, we are violating the constraint.
  for (Substructure_Hit_Statistics * q : queries_to_avoid) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;    // no violations
}

static int
violates_query_constraints(Molecule & m,
                           const int min_rings_this_molecule,
                           const int max_rings_this_molecule)
{
  if (m.empty()) {
    return 1;
  }

  const int matoms = m.natoms();

  if (lower_atom_count_cutoff > 0 && matoms < lower_atom_count_cutoff) {
    return 1;
  }

  if (matoms > upper_atom_count_cutoff) {
    return 1;
  }

  if (violates_uniqueness(m)) {
    return 1;
  }

  if (! _violates_query_constraints(m, min_rings_this_molecule, max_rings_this_molecule)) {
    return 0;
  }

  if (stream_for_discarded_molecules.active()) {
    stream_for_discarded_molecules.write(m);
  }

  return 1;
}


template <typename C>
int
random_molecular_permutations(Molecule & m,
                              const Can_Change & can_change,
                              int istart,
                              const int min_rings_this_molecule,
                              const int max_rings_this_molecule,
                              C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  int molecules_made_this_molecule = 0;

  for (int i = 0; i < number_perturbations; i++) {
    if (! do_a_single_random_permutation(m, can_change, transformation_chooser, random_state))  {
      continue;
    }

    molecules_made++;
    molecules_made_this_molecule++;

    if (violates_query_constraints(m, min_rings_this_molecule, max_rings_this_molecule)) {
      if (verbose > 2)
        cerr << "Permutation " << i << " violates constraints\n";
      break;
    }

    if (set_of_target_molecules.number_molecules()) {
      check_proximity_to_target_set(m, i, set_of_target_molecules);
    }

#ifdef FREQUENTLY_CHECK_VALENCES
    assert (m.valence_ok());
#endif

    do_the_output(m, molecules_made, istart, i, output);
  }

  molecules_made_per_copy[molecules_made_this_molecule]++;

  return 1;
}

template <typename C>
int
random_molecular_permutations(Molecule & m,
                              const Can_Change & can_change,
                              int min_rings_this_molecule,
                              int max_rings_this_molecule,
                              C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  if (discard_duplicates_within_each_molecule) {
    usmi.clear();
  }

  if (write_starting_molecule_once) {
    output.write(m);
  }

  if (write_starting_molecule_before_every_copy) {
    output.write(m);
  }

// We need to jump through some hoops to make sure we only write the
// starting molecule when there has been some output from a previous
// perturbation

  int molecules_previously_made = molecules_made;

  // For each copy of the starting molecule.
  for (int i = 0; i < nstart; i++) {
    if (i > 0) {
      if (write_starting_molecule_before_every_copy
          && molecules_previously_made < molecules_made) {
        output.write(m);
      }

      molecules_previously_made = molecules_made;
      changed_by.resize_keep_storage(0);
    }

    Molecule mcopy(m);
    mcopy.set_name(m.name());

    if (! random_molecular_permutations(mcopy, can_change, i, 
                min_rings_this_molecule, max_rings_this_molecule,
                transformation_chooser,
                random_state,
                output))
      return 0;
  }

  return 1;
}

static void
determine_min_and_max_rings(int nrings,
                            int & min_rings_this_molecule,
                            int & max_rings_this_molecule)
{
  if (min_rings_in_each_permutation > 0) {
    min_rings_this_molecule = min_rings_in_each_permutation;
  }

  if (max_rings_in_each_permutation > 0) {
    max_rings_this_molecule = max_rings_in_each_permutation;
  }

  if (max_rings_to_add >= 0) {
    int r = nrings + max_rings_to_add;
    if (max_rings_in_each_permutation > 0 && r > max_rings_in_each_permutation)
      r = max_rings_in_each_permutation;

    max_rings_this_molecule = r;
  }

  if (max_rings_to_remove >= 0) {
    int r = nrings - max_rings_to_remove;
    if (r < 0) {
      r = 0;
    }

    if (min_rings_in_each_permutation > 0 && r < min_rings_in_each_permutation) {
      r = min_rings_in_each_permutation;
    }

    min_rings_this_molecule = r;
  }

  return;
}

static void
preprocess(Molecule & m)
{
  m.reduce_to_largest_fragment();

  m.remove_all_chiral_centres();

  m.revert_all_directional_bonds_to_non_directional();

  arrange_kekule_forms_in_fused_rings(m);

  return;
}

static int
stop_for_any_reason() {
  if (run_for > 0)
  {
    time_t tnow = time(NULL);

    if (tnow - tzero > run_for) {
      if (verbose)
        cerr << "Computation complete after " << run_for << " seconds\n";
      return 1;
    }
  }

  if (molecules_made > stop_after_producing) {
    if (verbose)
      cerr << "Computation complete having produced " << molecules_made << " molecules\n";
    return 1;
  }

  if (stop_file.length() && dash_s(stop_file.null_terminated_chars())) {
    cerr << "Processing stopped by '" << stop_file << "'\n";
    return 1;
  }

  return 0;    // not stopped for any reason
}

template <typename C>
int
random_molecular_permutations(Molecule & m,
                              C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
#ifdef FREQUENTLY_CHECK_VALENCES
    assert (m.valence_ok());
#endif

  if (verbose > 1) {
    cerr << "Processing '" << m.name() << "'\n";
  }

  const int matoms = m.natoms();

  Can_Change can_change(matoms);

  if (never_change_query.number_elements()) {
    cerr << "The never change query option does not work now. See Ian\n";
#ifdef CAN_CHANGE_BROKEN
    initialise_never_change_atoms(m, can_change);
#endif
  }

  atoms_in_starting_molecules.extra(matoms);

  int max_rings_this_molecule = 0;
  int min_rings_this_molecule = 0;
  determine_min_and_max_rings(m.nrings(), min_rings_this_molecule, max_rings_this_molecule);

  if (verbose > 2)
    cerr << "For a molecule with " << m.nrings() << " rings, nrings must be between " << min_rings_this_molecule << " and " << max_rings_this_molecule << '\n';

#ifdef FREQUENTLY_CHECK_VALENCES
  assert (m.valence_ok());
#endif

  return random_molecular_permutations(m, can_change,
                                min_rings_this_molecule, max_rings_this_molecule,
                                transformation_chooser,
                                random_state,
                                output);
}


template <typename C>
int
random_molecular_permutations(C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
  for (int i = 0; i < input_molecules.number_elements(); i++)   // not fixed size since it may grow
  {
    if (stop_for_any_reason()) {
      return 1;
    }

    molecules_read++;
    
    Molecule * m=input_molecules[i];

    if (! random_molecular_permutations(*m, transformation_chooser,
                random_state, output)) {
      break;
    }

    if (run_for > 0 && (input_molecules.number_elements()-1) == i) {
      i = 0;
    }
  }

  for (int i = 0; i < iterations_through_new_molecules; ++i) {
    const int n = closer_molecules.number_molecules();

    for (int j = 0; j < n; ++j) {
      Molecule mcopy(*closer_molecules.moleculei(j));   // need to create a copy since we might replace ourselves in the closer set
      mcopy.set_name(closer_molecules.moleculei(j)->name());

      if (! random_molecular_permutations(mcopy, transformation_chooser, random_state, output)) {
        break;
      }
    }
  }

  return 1;
}

template <typename C>
int
random_molecular_permutations(data_source_and_type<Molecule> & input,
                              C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
  // Read all the molecules in `input` into `input_molecules`.
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    preprocess(*m);

    if (! m->valence_ok()) {
      cerr << "Invalid valence in input\n";
      cerr << m->smiles() << ' ' << m->name() << "\n";
      cerr << "Skipped\n";
      delete m;
      continue;
    }

    input_molecules.add(m);

    // Avoid regenerating starting molecules.
    if (discard_duplicates_across_all_molecules) {
      usmi[m->unique_smiles()] = 1;
    }
  }

  return random_molecular_permutations(transformation_chooser, random_state, output);
}

template <typename C>
int
random_molecular_permutations(const Command_Line& cl, 
                              FileType input_type,
                              C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! random_molecular_permutations(cl[i], input_type, transformation_chooser, random_state, output))
    {
      rc = i + 1;
      break;
    }
  }
  return rc;
}

template <typename C>
int
random_molecular_permutations(const char * fname, 
                              FileType input_type,
                              C & transformation_chooser,
                              RandomState& random_state,
                              Molecule_Output_Object & output)
{
  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return random_molecular_permutations(input, transformation_chooser, random_state, output);
}

/*
  Records look like

  directive probability

*/

static int
read_transformation_probability_file_record(const const_IWSubstring & buffer)
{
  const_IWSubstring directive, prob;

  if (! buffer.split(directive, ' ', prob))
  {
    cerr << "read_transformation_probability_file_record: cannot split record\n";
    return 0;
  }

  double p;

  if (! prob.numeric_value(p) || p < 0.0 || p > 1.0)
  {
    cerr << "Invalid probability '" << prob << "'\n";
    return 0;
  }

  if ("add_double_bond" == directive)
  {
    probability_add_double_bond = p;
    if (verbose)
      cerr << "Probability of adding double bond " << probability_add_double_bond << '\n';
  }
  else if ("remove_atom" == directive)
  {
    probability_remove_atom = p;
    if (verbose)
      cerr << "Probability of removing atom " << probability_remove_atom << '\n';
  }
  else if ("lower_bond_order" == directive)
  {
    probability_lower_bond_order = p;
    if (verbose)
      cerr << "Probability of lowering bond order " << probability_lower_bond_order << '\n';
  }
  else if ("break_aliphatic_ring" == directive)
  {
    probability_break_aliphatic_ring = p;
    if (verbose)
      cerr << "Probability of breaking aliphatic ring " << probability_break_aliphatic_ring << '\n';
  }
  else if ("make_ring" == directive)
  {
    probability_make_ring = p;
    if (verbose)
      cerr << "Probability of making ring " << probability_make_ring << '\n';
  }
  else if ("add_carbon" == directive)
  {
    probability_add_carbon = p;
    if (verbose)
      cerr << "Probability of adding carbon " << probability_add_carbon << '\n';
  } 
  else if ("probability_move_fragment" == directive)
  {
    probability_move_fragment = p;
    if (verbose)
      cerr << "Probability of moving fragment " << probability_move_fragment << '\n';
  }
  else if ("swap_adjacent_atoms" == directive)
  {
    probability_swap_adjacent_atoms = p;
    if (verbose)
      cerr << "Probability of swapping adjacent atoms " << probability_swap_adjacent_atoms << '\n';
  }
  else if ("destroy_aromatic_ring" == directive)
  {
    probability_destroy_aromatic_ring = p;
    if (verbose)
      cerr << "Probability of destroying an aromatic ring " << probability_destroy_aromatic_ring << '\n';
  }
  else if ("change_carbon_to_nitrogen" == directive)
  {
    probability_change_carbon_to_nitrogen = p;
    if (verbose)
      cerr << "Probability of changing carbon to nitrogen " << probability_change_carbon_to_nitrogen << '\n';
  }
  else if ("change_nitrogen_to_carbon" == directive)
  {
    probability_change_nitrogen_to_carbon = p;
    if (verbose)
      cerr << "Probability of changing nitrogen to carbon " << probability_change_nitrogen_to_carbon << '\n';
  }
  else if ("change_carbon_to_oxygen" == directive)
  {
    probability_change_carbon_to_oxygen = p;
    if (verbose)
      cerr << "Probability of changing carbon to oxygen " << probability_change_carbon_to_oxygen << '\n';
  }
  else if ("change_oxygen_to_something" == directive)
  {
    probability_change_oxygen_to_something = p;
    if (verbose)
      cerr << "Probability of changing Oxygen to something else " << probability_change_oxygen_to_something << '\n';
  }
  else if ("add_from_single_attachment_library" == directive)
  {
    probability_add_from_single_attachment_library = p;
    if (verbose)
      cerr << "Probability of adding from single attachment library " << probability_add_from_single_attachment_library << '\n';
  }
  else if ("add_from_aromatic_attachment_library" == directive)
  {
    probability_add_from_aromatic_attachment_library = p;
    if (verbose)
      cerr << "Probability of adding from aromatic attachment library " << probability_add_from_aromatic_attachment_library << '\n';
  }
  else if ("add_from_double_attachment_library" == directive)
  {
    probability_add_from_double_attachment_library = p;
    if (verbose)
      cerr << "Probability of adding from double attachment library " << probability_add_from_double_attachment_library << '\n';
  }
  else if ("change_ring_substitution" == directive)
  {
    probability_change_ring_substitution = p;
    if (verbose)
      cerr << "Probability of changing ring substitution " << probability_change_ring_substitution << '\n';
  }
  else if ("remove_terminal_functional_group" == directive)
  {
    probability_remove_terminal_functional_group = p;
    if (verbose)
      cerr << "Probability of removing a terminal functional group " << probability_remove_terminal_functional_group << '\n';
  }
  else if ("remove_embedded_functional_group" == directive)
  {
    probability_remove_embedded_functional_group = p;
    if (verbose)
      cerr << "Probability of removing an embedded functional group " << probability_remove_embedded_functional_group << '\n';
  }
  else if ("remove_terminal_ring" == directive)
  {
    probability_remove_terminal_ring = p;
    if (verbose)
      cerr << "Probability of removing a terminal ring " << probability_remove_terminal_ring << '\n';
  }
  else if ("switch_aromatic_ring" == directive)
  {
    probability_switch_aromatic_ring = p;
    if (verbose)
      cerr << "Probability of switching an aromatic ring " << probability_switch_aromatic_ring << '\n';
  }
  else if ("switch_aliphatic_ring" == directive)
  {
    probability_switch_aliphatic_ring = p;
    if (verbose)
      cerr << "Probability of switching an aliphatic ring " << probability_switch_aliphatic_ring << '\n';
  }
  else if("breed" == directive)
  {
    probability_breed=p;
    if (verbose)
      cerr << "Probability of breed " << probability_breed << '\n'; 
  }
  else if("inter_molecular" == directive)
  {
    probability_inter_permutation=p;
    if(verbose)
      cerr << "Probability of inter molecular permutation " << probability_inter_permutation << '\n'; 
  }
  else if ("split_fused_ring" == directive)
  {
    probability_split_fused_ring = p;
    if (verbose)
      cerr << "Probability of splitting a fused ring " << probability_split_fused_ring << '\n';
  }
  else if ("create_fused_ring" == directive)
  {
    probability_create_fused_ring = p;
    if (verbose)
      cerr << "Probability of creating a fused ring " << probability_create_fused_ring << '\n';
  }
  else
  {
    cerr << "What kind of transformation is this '" << directive << "'\n";
    return 0;
  }

  return 1;
}

static int
read_transformation_probability_file(iwstring_data_source & input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (0 == buffer.length())
      continue;

    if (buffer.starts_with('#'))
      continue;

    if (! read_transformation_probability_file_record(buffer))
    {
      cerr << "Cannot process transformation probability record '" << buffer << "'\n";
      return 0;
    }
  }
  return 1;
}

static int
read_transformation_probability_file(const IWString & fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_transformation_probability_file(input);
}

template <typename M>
int
read_fragments(data_source_and_type<M> & input,
               resizable_array_p<M> & molecules)
{
  M * m;

  while (nullptr != (m = input.next_molecule()))
  {
    molecules.add(m);
  }

  return molecules.number_elements();
}

template <typename M>
int
read_fragments(const IWString & fname,
               resizable_array_p<M> & m)
{
  data_source_and_type<M> input(FILE_TYPE_SMI, fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return read_fragments(input, m);
}

static int
read_double_attachment_library(const const_IWSubstring & fname,
                               resizable_array_p<Molecule> & lib)
{
  if (! read_fragments(fname, lib))
  {
    cerr << "Cannot read doubly attached fragment library '" << fname << "'\n";
    return 8;
  }

  for (int i = 0; i < lib.number_elements(); i++)
  {
    Molecule * f = lib[i];

    if (0 == f->hcount(0))
    {
      cerr << "No Hydrogens on first atom in fragment library '" << f->name() << "'\n";
      return 0;
    }

    int fatoms = f->natoms();

    atom_number_t asterisk_atom = INVALID_ATOM_NUMBER;
    for (int j = 0; j < fatoms; j++)
    {
      if (1 != f->ncon(j))
        continue;

      if (0 != f->atomic_number(j))
        continue;

      asterisk_atom = j;
      break;
    }

    if (INVALID_ATOM_NUMBER == asterisk_atom)
    {
      cerr << "No singly connected asterisk atom in '" << f->name() << "'\n";
      return 0;
    }

    atom_number_t attachment_point = f->other(asterisk_atom, 0);
    if (attachment_point > asterisk_atom)
      attachment_point--;

    f->remove_atom(asterisk_atom);
    index_of_star_atom_attachment.add(attachment_point);
  }

  return lib.number_elements();
}

static int
read_single_attachment_library(const const_IWSubstring & fname,
                               resizable_array_p<Molecule> & lib,
                               int first_atom_number_have_hydrogen)
{
  if (! read_fragments(fname, lib))
  {
    cerr << "Cannot read fragment library '" << fname << "'\n";
    return 0;
  }

  if (choose_attachment_points_randomly)
    return lib.number_elements();

  if (! first_atom_number_have_hydrogen)
    return lib.number_elements();

  for (int i = 0; i < lib.number_elements(); i++)
  {
    Molecule * f = lib[i];

    if (f->hcount(0))
      continue;

//  special case of the Nitro fragment

    const Atom * f0 = f->atomi(0);

    if (7 == f0->atomic_number() && 2 == f0->ncon() && 4 == f0->nbonds())
      continue;

    cerr << "No Hydrogens on first atom in fragment library '" << f->name() << "'\n";
    return 0;
  }

  return lib.number_elements();
}

static int
read_aromatic_ring_fusions(const_IWSubstring & fname,
                           resizable_array_p<Molecule> & lib)
{
  if (! read_fragments(fname, lib))
  {
    cerr << "Cannot read aromatic ring fusions '" << fname << "'\n";
    return 0;
  }

  return lib.number_elements();
}

static void
display_dash_T_options(std::ostream & output)
{
  output << "Options for controlling building towards a target population of molecules\n";
  output << " -T <fname>          file containing the desired molecules\n";
  output << " -T nmax=<n>         the number of proximal molecules to retain\n";
  output << " -T dmax=<d>         only retain proximal molecules that are closer than <d>\n";
  output << " -T fdbk=<n>         feed back newly identified proximal molecules <n> times\n";
  output << " -T tom=<d>          turn off any target molecule matched to within <d>\n";
  output << " -T write=<fname>    write all molecules added to the proximal collection to <fname>\n";

  exit(1);
}

static int
random_molecular_permutations(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vo:A:E:x:c:C:p:S:L:f:F:P:r:R:s:q:Q:K:y:Y:G:J:U:M:a:T:N:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (! process_standard_aromaticity_options(cl, verbose > 1))
    {
      cerr << "Cannot process -A option\n";
      usage(11);
    }
  }
  else {
    set_global_aromaticity_type(Daylight);
  }

  set_copy_name_in_molecule_copy_constructor(1);

  if (cl.option_present('x')) {
    if (! cl.value('x', nstart)  || nstart < 1) {
      cerr << "The number of copies of each molecule to start with (-x) must be a whole positive number\n";
      usage(15);
    }
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', number_perturbations) || number_perturbations < 1) {
      cerr << "The number of permutations of each copy (-p) must be a positive whole number\n";
      usage(15);
    }
  }

  if (verbose) {
    cerr << "Will make " << nstart << " copies of each input molecule and make "
         << number_perturbations << " permutations of that molecule\n";
  }

  if (cl.option_present('P')) {
    IWString p = cl.string_value('P');

    if (! read_transformation_probability_file(p)) {
      cerr << "Cannot open transformation probability file '" << p << "'\n";
      return 8;
    }
  }

  if (cl.option_present('U')) {
    const_IWSubstring u = cl.string_value('U');

    if ("each" == u) {
      discard_duplicates_within_each_molecule = 1;
      if (verbose)
        cerr << "Duplidates within each molecule will be discarded\n";
    }
    else if ("all" == u)
    {
      discard_duplicates_across_all_molecules = 1;
      if (verbose)
        cerr << "Duplicates derived from any molecule willb e discarded\n";
    }
    else
    {
      cerr << "Unrecognised -U qualifier '" << u << "'\n";
      usage(5);
    }
  }

  int sequential_scan_of_transformations = 0;

  if (cl.option_present('M'))
  {
    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++))
    {
      if ("help" == m)
        display_miscellaneous_options('M', cerr);
      else if ("afras" == m)
      {
        allow_fused_rings_to_swap_atoms = 1;
        if (verbose)
          cerr << "Fused ring atoms will be allowed to swap atoms\n";
      }
      else if ("cblast" == m)
      {
        changed_by_contains_full_history = 0;
        if (verbose)
          cerr << "The changed by information will contain only the last change\n";
      }
      else if ("okcage" == m)
      {
        can_form_strongly_fused_ring_systems = 1;
        if (verbose)
          cerr << "Will allow transformations that produce cage molecules\n";
      }
      else if (m.starts_with("stopfile="))
      {
        m.remove_leading_chars(9);
        stop_file = m;

        if (verbose)
          cerr << "Will terminate processing if stop file '" << stop_file << "' present\n";
      }
      else if (m.starts_with("run="))
      {
        m.remove_leading_chars(4);

        if (! m.numeric_value(run_for) || run_for < 1) {
          cerr << "The '-M run=nnn' option must specify a valid +ve number of seconds\n";
          return 3;
        }

        if (verbose)
          cerr << "Will run for " << run_for << " seconds\n";

        time(&tzero);
      }
      else if (m.starts_with("prd="))
      {
        m.remove_leading_chars(4);

        if (! m.numeric_value(stop_after_producing) || stop_after_producing < 1) {
          cerr << "Will stop after producing " << stop_after_producing << " valid molecules\n";
          return 3;
        }

        if (verbose)
          cerr << "Will produce only " << stop_after_producing << " valid molecules\n";
      }
      else if ("scrand" == m)
      {
        choose_attachment_points_randomly = 1;

        if (verbose)
          cerr << "Will choose sidechain attachment points randomly\n";
      }
      else if ("ssintra" == m)
      {
        sequential_scan_of_transformations = 1;
        if (verbose)
          cerr << "Will sequentially scan all " << NUMBER_KNOWN_INTRA_PERMUTATIONS << " intra molecular transformations\n";
      }
      else if ("all" == m)
      {
        set_all_transformation_probabilities(1.0);

        if (verbose)
          cerr << "All transformations will always happen\n";
      }
      else if (m == "nowarnvalence")
      {
        warn_invalid_valence = 0;
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_miscellaneous_options('M', cerr);
      }
    }
  }

  if (cl.option_present('r'))
  {
    if (! cl.value('r', minimum_ring_size_allowed) || minimum_ring_size_allowed < 3)
    {
      cerr << "The minimum ring size (-r) option must be a whole number > 2\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will never form a ring of size " << minimum_ring_size_allowed << " or smaller\n";
  }

  if (cl.option_present('R'))
  {
    if (! cl.value('R', maximum_ring_size_allowed) || maximum_ring_size_allowed < 3)
    {
      cerr << "The maximum ring size (-R) option must be a whole number > 2\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will never form a ring of size " << maximum_ring_size_allowed << " or larger\n";
  }

  if (cl.option_present('y'))
  {
    if (! cl.value('y', min_rings_in_each_permutation) || min_rings_in_each_permutation < 1)
    {
      cerr << "The minimum number of rings (-y) option must be a whole positive number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard any molecule having " << min_rings_in_each_permutation << " or fewer rings\n";
  }

  if (cl.option_present('Y'))
  {
    if (! cl.value('Y', max_rings_in_each_permutation) || max_rings_in_each_permutation < 1)
    {
      cerr << "The maximum number of rings (-Y) option must be a whole positive number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will discard any molecule having " << max_rings_in_each_permutation << " or more rings\n";
  }

  if (cl.option_present('a'))
  {
    if (! number_assigner.initialise(cl, 'a', verbose))
    {
      cerr << "Cannot initialise number assigner\n";
      return 4;
    }
  }

  if (cl.option_present('T'))
  {
    IWString fname;
    IWString t;
    int max_size_extra_molecules = 0;
    float set_turn_off_molecules_matched_to_within = -1.0f;

    for (int i = 0; cl.value('T', t, i); ++i)
    {
//    cerr << "TVALUE " << t << '\n';
      if (t.starts_with("nmax="))
      {
        t.remove_leading_chars(5);
        if (! t.numeric_value(max_size_extra_molecules) || max_size_extra_molecules < 1)
        {
          cerr << "Max size of input molecule array (-T max=) must be a whole +ve integer\n";
          usage(1);
        }
      }
      else if (t.starts_with("dmax="))
      {
        t.remove_leading_chars(5);
        float d;
        if (! t.numeric_value(d) || d < 0.0f || d > 1.0f) {
          cerr << "The dmax (max distance to consider) option must contain a valid distance\n";
          usage(1);
        }

        closer_molecules.set_max_distance(d);
        if (verbose)
          cerr << "Newly created molecules added to proximal list if closer than " << d << '\n';
      }
      else if (t.starts_with("write="))
      {
        t.remove_leading_chars(6);

        if (! t.ends_with(".smi")) {
          t << ".smi";
        }

        if (! stream_for_more_proximal_molecules.open(t.null_terminated_chars())) {
          cerr << "Cannot open stream for more proximal molecules '" << t << "'\n";
          return 1;
        }

        if (verbose)
          cerr << "More proximal molecules written to '" << t << "'\n";
      }
      else if (t.starts_with("fdbk="))
      {
        t.remove_leading_chars(5);
        if (! t.numeric_value(iterations_through_new_molecules) || iterations_through_new_molecules < 1) {
          cerr << "The number of feedback iterations (fdbk=) must be a whole +ve number\n";
          display_dash_T_options(cerr);
        }
        if (verbose)
          cerr << "Will iterate through accumulated proximal molecules " << iterations_through_new_molecules << " times\n";
      }
      else if (t.starts_with("tom="))
      {
        t.remove_leading_chars(4);
        if (! t.numeric_value(set_turn_off_molecules_matched_to_within) || set_turn_off_molecules_matched_to_within < 0.0f || set_turn_off_molecules_matched_to_within > 1.0f)
        {
          cerr << "The turn off if matched within option (zwd=) must be a valid distance\n";
          display_dash_T_options(cerr);
        }
        if (verbose)
          cerr << "Ay target molecule matched to within " << set_turn_off_molecules_matched_to_within << " will no longer be examined\n";
      }
      else if ("help" == t) {
        display_dash_T_options(cerr);
      }
      else if (0 == fname.length()) {
        fname = t;
      }
      else {
        cerr << "Unrecognised -T qualifier '" << t << "'\n";
        display_dash_T_options(cerr);
      }
    }

    if (fname.empty()) {
      cerr << "Must specify name of the target molecules file via the -T option\n";
      usage(1);
    }

    initialise_properties_ratios();

    if (! set_of_target_molecules.build(fname.null_terminated_chars())) {
      cerr << "Cannot build set of target molecules '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Read " << set_of_target_molecules.number_molecules() << " target molecules from '" << fname << "'\n";
    }

    if (0 == max_size_extra_molecules) {
      cerr << "New proximal molecules being saved, must specify max size of container\n";
      display_dash_T_options(cerr);
    }

    closer_molecules.initialise(max_size_extra_molecules);

    if (set_turn_off_molecules_matched_to_within > 0.0f) {
      set_of_target_molecules.set_turn_off_molecules_matched_to_within(set_turn_off_molecules_matched_to_within);
    }
  }

  if (cl.option_present('G'))
  {
    const_IWSubstring g;
    int i = 0;
    while (cl.value('G', g, i++))
    {
      if (g.starts_with('+'))
      {
        g.remove_leading_chars (1);
        if (! g.numeric_value(max_rings_to_add) || max_rings_to_add < 0)
        {
          cerr << "INvalid max rings to add '" << g << "'\n";
          usage(3);
        }

        if (verbose)
          cerr << "Will discard molecules that have more than " << max_rings_to_add << " extra rings over the starting molecule\n";

        max_rings_to_add++;
      }
      else if (g.starts_with('-'))
      {
        g.remove_leading_chars(1);
        if (! g.numeric_value(max_rings_to_remove) || max_rings_to_remove < 0)
        {
          cerr << "INvalid max rings to remove '" << g << "'\n";
          usage(3);
        }

        if (verbose)
          cerr << "Will discard molecules that have lost more than " << max_rings_to_remove << " rings from the starting molecule\n";

        max_rings_to_remove++;
      }
      else
      {
        cerr << "Unrecognised ring count change qualifier '" << g << "'\n";
        usage(17);
      }
    }
  }

  bool use_seed = false;
  uint32_t seed = 0;
  if (cl.option_present('s')) {
    if (! cl.value('s', seed)) {
      cerr << "Invalid random number seed value\n";
      return 1;
    }
    if (verbose)
      cerr << "Random number seed " << seed << '\n';
    use_seed = true;
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, queries_to_avoid, verbose > 1, 'q'))
    {
      cerr << "Cannot process queries to avoid (-q option)\n";
      usage(5);
    }

    if (verbose)
      cerr << "Defined " << queries_to_avoid.number_elements() << " queries to avoid\n";
  }

  if (cl.option_present('Q'))
  {
    if (! process_queries(cl, queries_to_match, verbose > 1, 'Q'))
    {
      cerr << "Cannot process queries to match (-Q option)\n";
      usage(5);
    }

    if (verbose)
      cerr << "Defined " << queries_to_match.number_elements() << " queries to match\n";
  }

  if (cl.option_present('K'))
  {
    const_IWSubstring k = cl.string_value('K');

    if ("once" == k)
    {
      write_starting_molecule_once = 1;
    }
    else if ("copy" == k)
    {
      write_starting_molecule_before_each_copy = 1;
    }
    else if ("every" == k)
    {
      write_starting_molecule_before_every_copy = 1;
    }
    else
    {
      cerr << "Unrecognised -K qualifier '" << k << "'\n";
      usage(5);
    }
  }

  element_carbon = get_element_from_atomic_number(6);
  element_nitrogen = get_element_from_atomic_number(7);
  element_oxygen = get_element_from_atomic_number(8);
  element_sulphur = get_element_from_atomic_number(16);

  if (cl.option_present('L'))
  {
    int i = 0;
    const_IWSubstring l;
    while (cl.value('L', l, i++))
    {
      if (l.starts_with("single="))
      {
        l.remove_up_to_first('=');
        if (! read_single_attachment_library(l, single_attachment_fragment_library_join_at_first_atom, 1))    // 1 means first atom must have hydrogen
        {
          cerr << "Cannot read singly bonded attachments from '" << l << "'\n";
          return 1;
        }

        if (verbose)
          cerr << "Read " << single_attachment_fragment_library_join_at_first_atom.number_elements() << " singly attached fragments from '" << l << "'\n";
      }
      else if (l.starts_with("singleR="))
      {
        l.remove_up_to_first('=');
        if (! read_single_attachment_library(l, single_attachment_fragment_library_join_from_any_atom, 0))   // 0 means first atom does not need to have a hydrogen
        {
          cerr << "Cannot read singly bonded attachments from '" << l << "'\n";
          return 1;
        }

        if (verbose)
          cerr << "Read " << single_attachment_fragment_library_join_from_any_atom.number_elements() << " singly attached (any atom) fragments from '" << l << "'\n";
      }
      else if (l.starts_with("double="))
      {
        l.remove_up_to_first('=');
        if (! read_double_attachment_library(l, double_attachment_fragment_library))
        {
          cerr << "Cannot read singly bonded attachments from '" << l << "'\n";
          return 1;
        }

        if (verbose)
          cerr << "Read " << double_attachment_fragment_library.number_elements() << " doubly connected attachment fragments from '" << l << "'\n";
      }
      else if (l.starts_with("arom="))
      {
        l.remove_up_to_first('=');
        if (! read_single_attachment_library(l, aromatic_attachment_fragment_library, 1))
        {
          cerr << "Cannot read aromatic attachments from '" << l << "'\n";
          return 3;
        }

        if (verbose)
          cerr << "Read " << aromatic_attachment_fragment_library.number_elements() << " aromatic attachment fragments from '" << l << "'\n";
      }
      else if (l.starts_with("aromR="))
      {
        l.remove_up_to_first('=');
        if (! read_single_attachment_library(l, aromatic_rings, 0))
        {
          cerr << "Cannot read aromatic rings from '" << l << "'\n";
          return 3;
        }

        if (verbose)
          cerr << "Read " << aromatic_rings.number_elements() << " aromatic rings from '" << l << "'\n";
      }
      else if (l.starts_with("aliphR="))
      {
        l.remove_up_to_first('=');
        if (! read_single_attachment_library(l, aliphatic_rings, 0))
        {
          cerr << "Cannot read aliphatic rings from '" << l << "'\n";
          return 3;
        }

        if (verbose)
          cerr << "Read " << aliphatic_rings.number_elements() << " aliphatic rings from '" << l << "'\n";
      }
      else if (l.starts_with("fsarom="))
      {
        l.remove_up_to_first('=');
        if (! read_aromatic_ring_fusions(l, aromatic_ring_fusions))
        {
          cerr << "Cannot read aromatic ring fusions from '" << l << "'\n";
          return 2;
        }

        if (verbose)
          cerr << "Read " << aromatic_ring_fusions.number_elements() << " aromatic ring fusions from '" << l << "'\n";
      }
      else
      {
        cerr << "Unrecognised library qualifier '" << l << "'\n";
        usage(18);
      }
    }
  }

// Note that this next case does not consider the join at any point library...

  if (single_attachment_fragment_library_join_at_first_atom.empty() && 0.0 != probability_add_from_single_attachment_library)
  {
    cerr << "No single attachment fragments available, probability set to 0.0\n";
    
    probability_add_from_single_attachment_library = 0.0;
  }

  if (double_attachment_fragment_library.empty() && 0.0 != probability_add_from_double_attachment_library)
  {
    cerr << "No double attachment fragments available, probability set to 0.0\n";
    
    probability_add_from_double_attachment_library = 0.0;
  }

  if (aromatic_attachment_fragment_library.empty() && 0.0 != probability_add_from_aromatic_attachment_library)
  {
    cerr << "No aromatic attachment fragments available, probability set to 0.0\n";
    
    probability_add_from_aromatic_attachment_library = 0.0;
  }

  if (aromatic_rings.empty() && 0.0 != probability_switch_aromatic_ring)
  {
    cerr << "No aromatic rings available, probability set to 0.0\n";
    
    probability_switch_aromatic_ring = 0.0;
  }

  if (aliphatic_rings.empty() && 0.0 != probability_switch_aliphatic_ring)
  {
    cerr << "No aliphatic rings available, probability set to 0.0\n";
    
    probability_switch_aliphatic_ring = 0.0;
  }

  if (aromatic_ring_fusions.empty())
  {
    cerr << "No aromatic rings to fuse, probability set to 0.0\n";

    probability_create_fused_ring = 0.0;
  }

  if (cl.option_present('F'))
  {
    if (! process_queries(cl, never_change_query, verbose, 'F'))
    {
      cerr << "Cannot process never change queries (-F)\n";
      return 6;
    }

    if (verbose)
      cerr << "Defined " << never_change_query.number_elements() << " queries for fixed atoms\n";
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  }
  else if (! all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.empty())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  number_known_permutations = NUMBER_KNOWN_PERMUTATIONS;
  known_permutation[0] = &permutation_add_double_bond;
  known_permutation[1] = &permutation_remove_atom;
  known_permutation[2] = &permutation_lower_bond_order;
  known_permutation[3] = &permutation_break_aliphatic_ring;
  known_permutation[4] = &permutation_make_ring;
  known_permutation[5] = &permutation_add_carbon;
  known_permutation[6] = &permutation_change_nitrogen_to_carbon;
  known_permutation[7] = &permutation_move_fragment;
  known_permutation[8] = &permutation_swap_adjacent_atoms;
  known_permutation[9] = &permutation_change_carbon_to_nitrogen;
  known_permutation[10] = &permutation_change_carbon_to_oxygen;
  known_permutation[11] = &permutation_change_oxygen_to_something;
  known_permutation[12] = &permutation_add_from_single_attachment_library;
  known_permutation[13] = &permutation_add_from_aromatic_attachment_library;
  known_permutation[14] = &permutation_add_from_double_attachment_library;
  known_permutation[15] = &permutation_change_ring_substitution;
  known_permutation[16] = &permutation_remove_terminal_functional_group;
  known_permutation[17] = &permutation_remove_embedded_functional_group;
  known_permutation[18] = &permutation_remove_terminal_ring;
  known_permutation[19] = &permutation_switch_aromatic_ring;
  known_permutation[20] = &permutation_switch_aliphatic_ring;
  known_permutation[21] = &permutation_split_apart_fused_ring;
  known_permutation[22] = &permutation_create_fused_ring;
  known_permutation[23] = &permutation_add_from_single_attachment_library_random_join;

  known_permutation[24] = &permutation_breed;
  known_permutation[25] = &permutation_breed_with_close_to_target;

  if (! cl.option_present('c'))
    ;
  else if (cl.value('c', lower_atom_count_cutoff) && lower_atom_count_cutoff > 0)
  {
    if (verbose)
      cerr << "Will exclude molecules with fewer than " << lower_atom_count_cutoff << " atoms\n";
  }
  else
  {
    cerr << "Cannot discern lower atom count cutoff from " << cl.option_value('c') << "'\n";
    usage(48);
  }

  if (!cl.option_present('C'))
	;
  else if (! cl.value('C',upper_atom_count_cutoff) || upper_atom_count_cutoff < lower_atom_count_cutoff)
  {
    cerr << "The upper atom count cutoff must be a whole positive number >= " << lower_atom_count_cutoff << '\n';
    usage(5);
  }
  else if (verbose)
  {
     cerr << "Will exclude molecules with more than " << upper_atom_count_cutoff << " atoms\n";
  }


  if (cl.option_present('J'))
  {
    const_IWSubstring j = cl.string_value('J');

    if (! cl.option_present('o'))
      stream_for_discarded_molecules.add_output_type(FILE_TYPE_SMI);
    else if (! stream_for_discarded_molecules.determine_output_types(cl, 'o'))
    {
      cerr << "Cannot determine output type(s) for rejection file\n";
      return 5;
    }

    if (stream_for_discarded_molecules.would_overwrite_input_files(cl, j))
    {
      cerr << "Rejection file cannot overwrite input file(s)\n";
      return 8;
    }

    if (! stream_for_discarded_molecules.new_stem(j))
    {
      cerr << "Cannot initialise stem for discarded molecules '" << j << "'\n";
      return 5;
    }

    if (verbose)
      cerr << "Discarded molecules written to '" << j << "'\n";
  }

  Molecule_Output_Object output;

  if (! cl.option_present('o'))
  {
    output.add_output_type(FILE_TYPE_SMI);

    if (verbose)
      cerr << "Default smiles output\n";
  }
  else if (! output.determine_output_types(cl, 'o'))
  {
    cerr << "Cannot determine output type(s)\n";
    usage(8);
  }

  if (cl.option_present('S'))
  {
    const_IWSubstring s = cl.string_value('S');

    if (output.would_overwrite_input_files(cl, s))
    {
      cerr << "Cannot overwrite input file(s)\n";
      return 6;
    }

    if (! output.new_stem(s))
    {
      cerr << "Cannot initialise output stem '" << s << "'\n";
      return 11;
    }
  }
  else if (! output.new_stem("-"))
  {
    cerr << "Cannot recirect output to stdout\n";
    return 4;
  }

  RandomState random_state;

  int rc = 0;
  if (sequential_scan_of_transformations)
  {
    Stream_of_Sequential_Numbers s(NUMBER_KNOWN_INTRA_PERMUTATIONS);

    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! random_molecular_permutations(cl[i], input_type, s, random_state, output))
      {
        rc = i + 1;
        break;
      }
    }
  }
  else
  {
    if (use_seed) {
      Stream_of_Random_Numbers srn(NUMBER_KNOWN_INTRA_PERMUTATIONS, seed);
      rc = random_molecular_permutations(cl, input_type, srn, random_state, output);
    }
    else
    {
      Stream_of_Random_Numbers srn(NUMBER_KNOWN_INTRA_PERMUTATIONS);
      rc = random_molecular_permutations(cl, input_type, srn, random_state, output);
    }
  }

  if (verbose)
  {
    cerr << "Processed " << molecules_read << " molecules\n";
    cerr << "Created " << molecules_made << " molecules\n";
    cerr << "Wrote " << molecules_written << " molecules\n";
    cerr << "Input molecules had between " << atoms_in_starting_molecules.minval() << " and " << atoms_in_starting_molecules.maxval();
    if (atoms_in_starting_molecules.n() > 1)
      cerr << " ave " << atoms_in_starting_molecules.average();
    cerr << " atoms\n";
    cerr << "Final molecules had between " << atoms_in_molecules_written.minval() << " and " << atoms_in_molecules_written.maxval();
    if (atoms_in_molecules_written.n() > 1)
      cerr << " ave " << atoms_in_molecules_written.average();
    cerr << " atoms\n";

    for (int i = 0; i < molecules_made_per_copy.number_elements(); i++)
    {
      if (molecules_made_per_copy[i]) {
        cerr << molecules_made_per_copy[i] << " copies made " << i << " permutations\n";
      }
    }

    if (verbose > 1)
    {
      if (molecules_with_too_few_rings)
        cerr << molecules_with_too_few_rings << " permutations which had too few rings\n";
      if (molecules_with_too_many_rings)
        cerr << molecules_with_too_many_rings << " permutations which had too many rings\n";
    }

    if (acc_dist_to_target.n())
    {
      cerr << acc_dist_to_target.n() << " shortest distances to target molecules, btw " << acc_dist_to_target.minval() << " and " << acc_dist_to_target.maxval() << " ave " << static_cast<float>(acc_dist_to_target.average()) << '\n';
    }
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  int rc = random_molecular_permutations(argc, argv);

  return rc;
}
