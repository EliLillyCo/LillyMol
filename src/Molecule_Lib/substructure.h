#ifndef MOLECULE_LIB_SUBSTRUCTURE_H_
#define MOLECULE_LIB_SUBSTRUCTURE_H_

#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#ifdef IW_USE_TBB_SCALABLE_ALLOCATOR
#include "tbb/scalable_allocator.h"
#endif

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/iwminmax.h"
#include "Foundational/iwmisc/logical_expression.h"
#include "Foundational/iwmisc/matcher.h"
#include "Foundational/iwmisc/minmaxspc.h"
#include "Foundational/iwmisc/msi_object.h"
#include "Foundational/iwmisc/set_or_unset.h"
#include "Foundational/iwstring/iwstring.h"

#include "atom_alias.h"
#include "iwmtypes.h"
#include "geometric_constraints.h"

class Molecule_to_Match;
class Target_Atom;
class MDL_Molecule;
class Atom_Typing_Specification;

#include "molecule.h"

#include "chiral_centre.h"
#include "set_of_atoms.h"
#include "Molecule_Lib/substructure.pb.h"

class MDL_File_Data;
class MDL_Atom_Data;
class Bond_and_Target_Atom;
class ISIS_Atom_List;

#define NAME_OF_QUERY_OBJECT "query"
#define NAME_OF_COMPOSITE_QUERY_OBJECT "composite_query"

#define SUBSTRUCTURE_NOT_SPECIFIED -637

#define NAME_OF_OPERATOR_ATTRIBUTE "operator"

class Substructure_Atom;

class Atomic_Smarts_Component;

class ThreeDots;

/*
  Whenever we get an embedding, we need to keep the identity of the
  Substructure_Atom objects which made the match.

  _preference_value will be derived from the preference values of the
  atoms which comprise the match
*/

class Query_Atoms_Matched : public resizable_array<Substructure_Atom *>
{
  private:
    int _preference_value;

  public:
    Query_Atoms_Matched ();
    ~Query_Atoms_Matched ();

    int preference_value () const { return _preference_value;}

    int increment_preference_value (int p) {_preference_value += p; return _preference_value;}
};


/*
  The chiral centre object is a little strange. During building, all
  we will know are atom numbers, so we use a regular Chiral_Centre
  object. Once the whole query has been built, then we can convert
  these numeric values into pointers to atoms
*/

class Substructure_Chiral_Centre
{
  private:
    Chiral_Centre * _numeric;

    const Substructure_Atom * _centre;
    const Substructure_Atom * _top_front;
    const Substructure_Atom * _top_back;
    const Substructure_Atom * _left_down;
    const Substructure_Atom * _right_down;

//  We can have a chiral centre with either 4 of 3 explicit connections

    int  _number_explicit_connections;

  public:
    Substructure_Chiral_Centre();
    ~Substructure_Chiral_Centre();

    int debug_print (std::ostream &) const;

    void set_centre(const Substructure_Atom * a) { _centre = a;}
    void set_top_front(const Substructure_Atom * a);
    void set_top_back(const Substructure_Atom * a);
    void set_left_down(const Substructure_Atom * a);
    void set_right_down(const Substructure_Atom * a);

    const Chiral_Centre * numeric_chiral_centre() const { return _numeric;}

    const Substructure_Atom * centre () const { return _centre;}

    int invert();

    int copy_chirality(const Chiral_Centre &);   // fills the _numeric object
  
    int extract_connections(const Chiral_Centre & c);

    int construct_from_msi_attribute (const msi_attribute &);
    int write_msi (std::ostream &, const IWString &, const char *) const;
    int BuildProto(SubstructureSearch::SubstructureChiralCenter& proto) const;

    int is_matched (const Molecule * m) const;
};

// Class describing the ... directive from a smarts.
// That can be either
//  [smarts]...[smarts] which just means no matched atoms
//                      between those matched atoms, with no
//                      constraints on what is between.
// The ... directive can be qualified with a {} modifier
//   [smarts]...{}[smarts]
// That modifier can contain:
//   <digit>           exact distance between atoms.
//   <digit>,<digit>   set of distances.
//   <digit>-<digit>   range of distances.
//  e.g.  OC(=O)...C(=O)-O   to find esters "facing away" from each other.
//  e.g.  OC(=O)...{2}C(=O)-O  exactly two unmatched atoms between the carbons.
//  e.g.  OC(=O)...{3,5}C(=O)-O  3 or 5 unmatched atoms between the carbons.
//  e.g.  OC(=O)...{3-5}C(=O)-O  3 - 5 unmatched atoms between the carbons.
//  e.g.  OC(=O)...{>3}C(=O)-O  4 or more unmatched atoms between the carbons.
// We can also specify what kinds of atoms are (not) found in the 
// ... atoms. This is specified with single atom smarts.
// {[c]}   every one of the atoms must be an aromatic carbon.
// {0[c]}  none of the atoms can be an aromatic carbon.
// {>0[c]} there must be at least one atom that is an aromatic carbon.
// {1[c]}  there must be exactly one atom that is an aromatic carbon.
// {2[c]}  there must be exactly two atoms that are aromatic carbon.
// {>2[c]}  there must more than two atoms that are aromatic carbon.
// {1[c]&0[#7]} Exactly 1 c, and no Nitrogens.
// {1[c];0[#7]} Exactly 1 c, and no Nitrogens.
// {4;0[#7]}  Four bonds separating, but no Nitrogens.
// {4&0[#7]}  Four bonds separating, but no Nitrogens.
// {0[#7];4}  Four bonds separating, but no Nitrogens.
// {0[#7]&4}  Four bonds separating, but no Nitrogens.
// {>3;0[c]}  More than 3 bonds, and no aromatic Carbons.
// {3-5;[c]}  Between 3 and 5 bonds, all atoms must be aromatic carbon.
// {3;2[C]}   Three bonds separating, two Aliphatic Carbon - same as {3;[C]}.
// {3,4}      3 or 4 bonds between.
// {>0[$(C(=O)N)]}  an aimde carbon atom in ...
// So a fully specified query might look like
// N-c...{>4[c]}c-N which would be two nitrogen substituents, separated
// by more than 4 aromatic carbons (not counting the two shown here).

// A class used for parsing these directives.
// Each token in the {} directive can be parsed into a NMAB_Token.
class NMAB_Token 
{
  private:
    // The two matched atoms that define the region.
    const int _a1;
    const int _a2;

    // One of IW_LOGEXP_* logical operators.
    int _op;
    // An optional < (-1) or > (+1)
    int _relational ;
    // Any numbers perceived. If a range, it will be enumerated.
    resizable_array<int> _number;

    // The smarts. Must be enclosed in [].
    IWString _smarts;

    // Private functions.

    bool IsComplete() const;
    // Look for either a single number, or a range starting at s[i].
    // s[i] must be a digit.
    // Returns 0 if failure.
    int GetNumberOrRange(const const_IWSubstring& s, int & i);

  public:
    NMAB_Token(int a1, int a2);

    int op() const { return _op;}
    int relational() const { return _relational;}
    const resizable_array<int>& numbers() const { return _number;}
    int IsRange() const {
      return _number.size() >= 1;
    }

    // If this is a single number, and 'rhs' is also a single number
    // with an or operator, consolidate with 'this'.
    int CombineAsRange(const NMAB_Token& rhs);

    // Parse starting at s[i].
    // Returns 1 if successful, zero otherwise.
    // Updates i if successful.
    int Parse(const const_IWSubstring& s, int& i);
    int Parse(const const_IWSubstring& s);

    // Mostlu used for debugging.
    const IWString& smarts() const { return _smarts;}
};

// Once the NMAB directive is parsed into NMAB_Token's, those can be transformed
// into NMAB_Operator's.
// If _query is not set, then _number is interpreted as match atom count.
// Kind of a lkudge, I did think of having a base class and virtual...
class NMAB_Operator 
{
  private:
    // An optional < (-1) or > (+1)
    int _relational;

    // Any numbers perceived. If a range, it will be enumerated.
    resizable_array<int> _number;

    // If built from smarts.
    std::unique_ptr<Substructure_Atom> _query;

    // Private functions.

    // Does 's' match the numeric specifications.
    int _relational_matches(const int s) const;

  public:
    NMAB_Operator();

    int Build(const NMAB_Token& token);

    int Matches(Molecule_to_Match& target, const Set_of_Atoms& unmatched_atoms);
};
  
// A testable function for parsing the {} directive.
// Note argument is a local copy.
extern int
TokeniseNMABSpecification(const_IWSubstring s, resizable_array_p<NMAB_Token>& tokens);

// A resizable_array_p<NMAB_Token> can be turned into one of these, that is actually
// used for matching. 

class No_Matched_Atoms_Between {
  private:
    // The matched atoms that define the ends of the no matched atoms between segment.
    const int _a1;
    const int _a2;

    IW_Logical_Expression _logexp;

    resizable_array_p<NMAB_Operator> _specs;

  // Private functions

    int _matches_exhaustive(Molecule_to_Match& target,
                                              const Query_Atoms_Matched& matched_atoms,
                                              const int * matched,
                                              int * in_path,
                                              const atom_number_t my_atom,
                                              const atom_number_t destination,
                                              const int current_distance);
    int _is_a_match(Molecule_to_Match& target,
                                const Query_Atoms_Matched& matched_atoms,
                                const int * in_path);
  public:
    No_Matched_Atoms_Between(const int a1, const int a2);

    // Parse any {...} directive given.
    int Initialise(const ThreeDots& three_dots);

    // Given that 'matched_atoms' is an embedding, and 'matched'
    // contains the atoms matched by 'matched_atoms', do we
    // match this molecule.
    // 'tmp' will be a scratch array sized to the number of atoms
    // in the molecule being searched.
    int Matches(Molecule_to_Match& target,
                        Query_Atoms_Matched & matched_atoms,
                        const int * matched,
                        int * tmp);
    // The Matches method looks only at the first shortest path
    // between the atoms found. Sometimes it will be important
    // to examine all paths found.
    int MatchesExhaustive(Molecule_to_Match& target,
                        Query_Atoms_Matched & matched_atoms,
                        const int * matched,
                        int * tmp);
};

namespace down_the_bond {

// Looking down from one matched atom to another, how many atoms are
// there - matched or otherwise. Includes matched atom _a2.
// the -{} directive is parsed into one of these.
// TODO: expand to other concepts.
class DownTheBond {
  private:
    int _a1;
    int _a2;

    Min_Max_Specifier<int> _natoms;

  public:
    DownTheBond();
    DownTheBond(int a1);

    int ConstructFromProto(const SubstructureSearch::DownTheBond& proto);
    int BuildProto(SubstructureSearch::DownTheBond& proto) const;

    // Called when parsing a smarts. The opening and closing {}
    // should have been removed.
    int Build(const const_IWSubstring& buffer);

    // When parsing smarts the identity of the atom at the other end
    // is not yet known, so it must be settable.
    void set_a2(int s) {
      _a2 = s;
    }

    int Complete() const {
      return _a2 >= 0;
    }

    int Matches(Molecule_to_Match& target,
                Query_Atoms_Matched & matched_atoms,
                int * visited);
};

}  // namespace down_the_bond

/*
  A bond can be one of several kinds.
    It can match a single bond type (1, 2, 3)
    It can match an aromatic bond
    It can match a ring bond
*/

class Substructure_Bond_Specifier_Base
{
  protected:

    Substructure_Bond_Specifier_Base * _next;

  public:
    Substructure_Bond_Specifier_Base ();
    virtual ~Substructure_Bond_Specifier_Base ();

    virtual int ok () const = 0;
    virtual int debug_print (std::ostream &, const IWString & = "") const = 0;

    Substructure_Bond_Specifier_Base * next () const { return _next;}

    void add_to_chain (Substructure_Bond_Specifier_Base *);

    virtual int involves_aromatic_bond_specification (int &) const = 0;

    int bond_type_as_string (IWString & btype) const;

//  virtual int construct_from_msi_object (const msi_object *) = 0;
//  virtual int write_as_msi_attribute (std::ostream & os, const IWString & ind, atom_number_t other_atom) const = 0;

//  virtual int write_msi (std::ostream &, int &, int = 0) = 0;

    // Write smarts to a stream.
    virtual int smarts (std::ostream &) const = 0;
    // Append smarts to a string.
    virtual int Smarts(IWString& smt) const = 0;

    virtual int matches (Bond_and_Target_Atom & b) const = 0;

    int copy (const Bond * b, int);

/*  virtual int construct_from_smarts (const char * smarts,
                             int chars_to_process,
                             int & characters_processed) = 0;*/
};

/*
  Bond specifier object used for single, double or triple bonds
*/

class Substructure_Bond_Specifier_Type : public Substructure_Bond_Specifier_Base
{
  private:
    bond_type_t _bond_type;

  public:
    Substructure_Bond_Specifier_Type (bond_type_t);

    int ok () const;
    int debug_print (std::ostream &, const IWString &) const;

    int smarts (std::ostream &) const;
    int Smarts(IWString& smt) const;

    int involves_aromatic_bond_specification (int & r) const;

    int construct_from_msi_object (const msi_object &);

    int matches (Bond_and_Target_Atom & b) const;
};

/*
  There are two kinds of ring bond specifiers. One which matches a given number
  of rings, and the other which matches ring or non-ring
*/

class Substructure_Bond_Specifier_Ring : public Substructure_Bond_Specifier_Base
{
  private:
    boolean _ring;

  public:
    Substructure_Bond_Specifier_Ring (int);

    int ok () const;
    int debug_print (std::ostream &, const IWString &) const;

    int smarts (std::ostream &) const;
    int Smarts(IWString& smt) const;

    int involves_aromatic_bond_specification (int & r) const;

    int construct_from_msi_object (const msi_object &);

    int matches (Bond_and_Target_Atom & b) const;
};

class Substructure_Bond_Specifier_NRings : public Substructure_Bond_Specifier_Base
{
  private:

    int _nrings;

  public:
    Substructure_Bond_Specifier_NRings (int);

    int ok () const;
    int debug_print (std::ostream &, const IWString &) const;

    int smarts (std::ostream &) const;
    int Smarts(IWString& smt) const;

    int involves_aromatic_bond_specification (int & r) const;

    int construct_from_msi_object (const msi_object &);

    int matches (Bond_and_Target_Atom & b) const;
};

class Substructure_Bond_Specifier_Aromatic : public Substructure_Bond_Specifier_Base
{
  private:
    boolean _aromatic;

  public:
    Substructure_Bond_Specifier_Aromatic (boolean);

    int ok () const;
    int debug_print (std::ostream &, const IWString &) const;

    int smarts (std::ostream &) const;
    int Smarts(IWString& smt) const;

    int involves_aromatic_bond_specification (int & r) const;

    int construct_from_msi_object (const msi_object &);

    int matches (Bond_and_Target_Atom & b) const;
};

/*
  A substructure bond is matched either by a value of _bond_types or a
  logical expression of Substructure_Bond_Specifier_Base objects linked
  by a logical_expression

  The bits set on _bond_types are treated as OR matches, so that's how
  we do the smiles default, of single or aromatic
*/

class Substructure_Bond
{
  private:

    bond_type_t _bond_types;

    Substructure_Atom * _a1;

//  We have a logical expression involving the bond smarts

    IW_Logical_Expression _logexp;

//  and the beginning of a linked list of bond specifications

    Substructure_Bond_Specifier_Base * _b;

// private functions

    void _default_values ();
    int  __write_as_smarts (std::ostream &) const;
    int  _write_as_smarts (std::ostream &) const;
    int _construct_from_smarts (const char * smarts, int chars_to_process, int & characters_processed);

  public:
    Substructure_Bond ();
    ~Substructure_Bond ();

    int ok () const;
    int debug_print (std::ostream &, const IWString &) const;

    int construct_from_smarts (const char * smarts, int chars_to_process, int & characters_processed);
    int construct_from_msi_attribute (const msi_object *, extending_resizable_array<Substructure_Atom *> &);
    int construct_from_msi_attribute (const msi_attribute *);
    int write_as_msi_attribute (std::ostream & os, int indentation) const;
    // Since the bond is represented either as a SubstructureBond or 
    // as a bond smarts, building the proto requires access to the parent.
    int BuildProto(SubstructureSearch::SubstructureAtom& proto) const;

    int BuildProtoNoBond(SubstructureSearch::SubstructureBond& proto) const;

    int make_single_or_aromatic ();

    int set_must_be_in_a_ring (int);

    Substructure_Atom * a () const { return _a1;}
    void  set_atom (Substructure_Atom * a);

    void set_type (bond_type_t b) { _bond_types = b;}
    void set_match_any ();

    bond_type_t types_matched () const { return _bond_types;}

    int bond_type_as_string (IWString &) const;

    int involves_aromatic_bond_specification (int & r) const;

    int can_be_written_as_attribute () const;

    int copy (const Bond *, int);

    int matches (Bond_and_Target_Atom & b);

    void set_bond_type (bond_type_t b) { _bond_types = b;}

    // Represent the bond as a smarts component.
    IWString Smarts() const;
};

/*
  An atom can be specified in many ways:
    By atomic number: for now, no isotopic information is allowed.
    By ncon
    By nbonds
    By ncon2     (ncon of this atom + (ncon - 1) of its connections)
    By formal charge
    By nrings    (the number of rings the atom must be in)

  _ring_size can be used to specify the size of the ring the atom is in.
  Note that I'm not entirely happy with this, and may again change it to
  be a resizable_array<int> of possible ring sizes.
*/

class Substructure_Atom_Specifier
{
  protected:

//  This next group of things can either be specified for the Substructure_Atom,
//  or with individual Substructure_Atom_Specifiers.

    // _element is no longer used for matching, _element_unique_id is now used.
    // The _element array is retained so elemental properties can be retrieved.
    resizable_array<const Element *> _element;
    resizable_array<int>   _element_unique_id;
    extending_resizable_array<int> _element_uid;

    iwmatcher::Matcher<int> _ncon;
    Min_Max_Specifier<int> _ncon2;
    iwmatcher::Matcher<int> _nbonds;
    iwmatcher::Matcher<int> _valence;  // From the 'v' directive in smarts.
    Min_Max_Specifier<int> _formal_charge;
    iwmatcher::Matcher<int> _nrings;
    iwmatcher::Matcher<int> _ring_bond_count;
    iwmatcher::Matcher<int> _ring_size;
    iwmatcher::Matcher<int> _hcount;
    aromaticity_type_t     _aromaticity;
    int                    _chirality;
    iwmatcher::Matcher<int> _aromatic_ring_size;
    iwmatcher::Matcher<int> _aliphatic_ring_size;
    iwmatcher::Matcher<int> _attached_heteroatom_count;
    iwmatcher::Matcher<int> _lone_pair_count;
    iwmatcher::Matcher<int> _unsaturation;
    iwmatcher::Matcher<int> _daylight_x;
    iwmatcher::Matcher<uint32_t> _isotope;
    int                    _userAtomType = 0;
    iwmatcher::Matcher<int> _aryl;
    iwmatcher::Matcher<int> _fused_system_size;
    iwmatcher::Matcher<int> _vinyl;
    int                    _all_rings_kekule;
    // initialise to negative to indicate no specification.
    // 0 means NOT spiro.
    // 1 means IS spiro.
    int                    _spiro;

//  int                    _carbocycle;

    iwmatcher::Matcher<int> _heteroatoms_in_ring;

    CahnIngoldPrelog _cip_chirality;

//  We can specify that an atom only match an atom in the molecular
//  spinach or not. Note that this is not an atom property, but can
//  only be discerned by looking at all the matched atoms after the
//  match is done, so we don't slow down our matching by doing this.

//  Value will be -1 (unspecified), 0 (not in spinach) or 1(must be in spinach)

    int _match_spinach_only;

//  We want to be able to look for an atom in a "terminal" ring.

    iwmatcher::Matcher<int> _scaffold_bonds_attached_to_ring;

    // If the query contains atom typing.
    atom_type_t _atom_type;

    // Matching non organic atoms is done as a special case. This is because
    // when a query is defined, not all non-organic elements may be known;
    int _match_non_organic;

//  Every atom specifier has a preference value. When a Substructure_Atom
//  is matched, it will determine an aggregate preference value for that
//  atom based on which Substructure_Atom_Specifier was matched.

    int _preference_value;

//  There are two symmetry concepts.
//   a. the number of atoms in a symmetry group (F in CF3 is 3) [DEGREE]
//   b. relationship to other matched atoms   [GROUP]

    iwmatcher::Matcher<int> _symmetry_degree;
    int _symmetry_group;

//  During matching is it desirable to keep track of how many attributes
//  have been processed. Once all attributes have been checked, we are
//  done

    int _attributes_specified;

//  private functions

    void _default_values ();

    int _adjust_nrings (int);
    int _adjust_ring_sizes (const List_of_Ring_Sizes &);

    int _set_implicit_hydrogens (Molecule & m, atom_number_t a) const;
    int _fill_min_ncon (Molecule & m, atom_number_t a) const;

    int _matches (Target_Atom &);

    int _create_environment_from_smarts (const_IWSubstring &);
    int _match_scaffold_bonds_attached_to_ring (Target_Atom & target);
    int _match_scaffold_bonds_attached_to_ring (Target_Atom & target_atom, const Ring & r) const;
    int _match_symmetry_degree(Target_Atom & target_atom) const;

    int _get_atomic_number_or_symbol (const char * smarts, const int characters_to_process);

    int _fetch_symmetry_group(const msi_object & msi);

    int _add_element (const atomic_number_t z);
    int _add_element (const Element * e);

    int _fetch_aromaticity(const SubstructureSearch::SubstructureAtomSpecifier & proto);
    void AddOrganicElements();
    void AddNonOrganicElements();

  public:
    Substructure_Atom_Specifier ();
    Substructure_Atom_Specifier (const Element *);

    ~Substructure_Atom_Specifier ();

    int ok () const;
    int debug_print (std::ostream &, const IWString & = "") const;
    int terse_details (std::ostream &, const IWString &) const;

    int BuildProto(SubstructureSearch::SubstructureAtomSpecifier& proto) const;

    int  preference_value () const { return _preference_value;}
    void set_preference_value (int s) { _preference_value = s;}

    int construct_from_msi_object (const msi_object &);

    int construct_from_smarts_token (const const_IWSubstring &);
    int construct_from_smarts_token (const char *, int);
    int construct_from_smiles_token (const const_IWSubstring &);
    int construct_from_smiles_token (const char *, int);

    int construct_from_proto(const SubstructureSearch::SubstructureAtomSpecifier&);

    int smarts (IWString &) const;

    // Non const method that will count everthing and set _attributes_specified.
    int count_attributes_specified ();
    int attributes_specified() const {
      return _attributes_specified;
    }

    int spinach_match_specified () const { return _match_spinach_only >= 0;}
    int spinach_match_value () const { return _match_spinach_only;}

    int check_internal_consistency (int) const;

    int chirality () const { return _chirality;}

    Atom * create_atom () const;

    void assign_unique_atom_numbers (int &);

    int symmetry_group() const { return _symmetry_group;}

    int set_element (const Element *);
    int atomic_number (atomic_number_t &) const;

//  Functions min_ncon and min_nbonds are somewhat unusual, because they are used
//  by create_molecule only

    int min_ncon () const;
    int min_nbonds () const;

    int set_ncon (int);
    int set_min_ncon (int);
    int set_max_ncon (int);

    int set_ncon2 (int);
    int set_min_ncon2 (int);
    int set_max_ncon2 (int);

//  This next group of functions are used to tell whether or not an
//  attribute has been specified for an atom. Used mostly by the
//  fingerprint creation functions, and create_molecule.

    int ncon_specification (int &) const;    // used by fingerprint routines.
    int hcount_specification (int &) const;
    int formal_charge_specification (formal_charge_t &) const;

    aromaticity_type_t aromaticity() const { return _aromaticity;}

    int set_nbonds (int);
    int set_min_nbonds (int);
    int set_max_nbonds (int);

    int set_valence (int);
    int set_min_valence (int);
    int set_max_valence (int);

    int set_nrings (int);
    int set_min_nrings (int i);
    int set_max_nrings (int i);

    int set_ring_membership (Molecule &);

    int formal_charge (formal_charge_t &) const;
    int set_formal_charge (formal_charge_t);

    int set_aromaticity (aromaticity_type_t);
    int update_aromaticity (aromaticity_type_t);

    int set_attached_heteroatom_count (int);
    
    int set_chirality (int);

    int involves_rings () const;
    int nrings () const;     // the smallest number of rings for any component
    int determine_ring_or_non_ring (int &) const;  // nrings value if specified
    int ring_sizes_specified (resizable_array<int> &) const;
    int involves_ring_specifications_for_bonds () const;

    int adjust_ring_specifications (int, const List_of_Ring_Sizes &);

//  The extra arguments are for the unary NOT operator, and any binary operator

    int write_msi (std::ostream &, int, int = 1, int = IW_LOGEXP_UNDEFINED);

    int matches (Target_Atom &);

    int create_molecule (Molecule &, int, int) const;

    int reconcile_and_conditions (const Substructure_Atom_Specifier * s);

    const resizable_array<const Element *> & element () const { return _element;};
    const Min_Max_Specifier<int> & formal_charge () const { return _formal_charge;}
    const iwmatcher::Matcher<uint32_t> & isotope () const { return _isotope;}
    const int & userAtomType () const { return _userAtomType;} 

    int spiro() const {
      return _spiro;
    }
};

class Substructure_Atom_Environment : public resizable_array_p<Substructure_Atom>
{
  private:
    IW_Logical_Expression _operator;

//  If all our components are just simple atoms, we don't have to allocate
//  a second already_matched array

    int _all_components_single_atoms;

//  private functions

    int _matches (Target_Atom & target_atom, int atoms_in_target, const int * already_matched, int * tmp);
    int _match_component (Target_Atom & target_atom,
                          int which_component, const int * already_matched_by_query, int * already_matched);
    int _perform_environment_search (Substructure_Atom * root_atom, int * already_matched);

  public:
    Substructure_Atom_Environment ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int number_elements () const { return _number_elements;}
    int active () const { return _number_elements;}

    int BuildProto(SubstructureSearch::SubstructureAtom& proto) const;
    int BuildFromProto(const SubstructureSearch::SubstructureAtomEnvironment& proto);

    int create_from_smarts (const Atomic_Smarts_Component &);
    int create_from_msi_object (msi_object &);

    int write_msi (std::ostream &, int &, int) const;

    int involves_aromatic_bond_specifications (int &) const;

    int matches (Target_Atom & target_atom, int, const int *);
};

class Smiles_Ring_Status;

class Molecule_to_Query_Specifications;

/*
  A Substructure_Atom object is a resizable_array of Substructure_Atom_Specifier
  objects, which are used for matching. In addition, there is a resizable array
  of Substructure_Atom_Specifier objects which are interpreted as exclusions.
  Exclusions are applied before matches, so, for example,
    if the exclusion list contained Carbon, but one of the 
    Substructure_Atom_Specifier's to be matched was also Carbon,
    that match would never be found, as it would be rejected by
    the exclusion first.

  There are really two types of Substructure_Atom objects. Root_Substructure_Atoms
  have no notion of a parent. At one stage I separated the Substructure_Atom
  object into two types, but this lead to too many other problems. These
  problems could have been solved with virtual functions, but I was unwilling
  to incur the associated run-time inefficiencies.

  TODO:ianwatson Not a good decision to make this inherit from a Substructure_Atom_Specifier.
  Instead, rely on the _components array.
*/

class Parse_Smarts_Tmp;

class Substructure_Atom : public Substructure_Atom_Specifier
{
  protected:

//  Each substructure atom has a unique id. This is re-assigned dynamically
//  during execution. If the user wants to respect the original numbering
//  of the atoms, we also need to retain the initial number
//  Note. These 3 items, _unique_id, _initial_atom_number and _atom_map_number are not used with
//  great clarity, their roles should be more clearly separated. Bad initial design, working on it...

    int _unique_id;

    int _initial_atom_number;

    int _atom_map_number;

    // Global matching conditions, like ring systems, can set global match numbers
    // for each atom. If requested, we can only match with an atom that has a given
    // global match id.
    // Note that in smarts, these are restricted to single digits.
    int _global_match_id;

    resizable_array_p<Substructure_Atom> _children;

//  For alternate matches, we use the or id

    int _or_id;

//  During matching, the atom with which we are currently matched

    Target_Atom * _current_hold_atom;

//  The Daylight atom environment 

    Substructure_Atom_Environment _environment;

//  The ring_id variable is used to force atoms to be in the same ring

    int _ring_id;

//  The _fused_system_id is used to force atoms to be in the same ring system

    int _fused_system_id;

//  We are really a set of Substructure_Atom_Specifier's

    IW_Logical_Expression _operator;

//  The components are evaluated under the control of _operator;

    resizable_array_p<Substructure_Atom_Specifier> _components;

//  Many times we need a unique means of identifying a given atom in more
//  than one substructure query - for example when the atom of interest may
//  be different query atom numbers in two parts of a composite query.
//  Each query_atom has a text identifier. Warning, we don't check for uniqueness
//  among these.

    IWString _text_identifier;

//  For many applications, we may want to have a numeric value associated with
//  each atom

    Set_or_Unset<double> _numeric_value;

//  When a match is found, we can optionally exclude a query atom from
//  the resulting Set_of_Atoms embedding objects.

    int _include_in_embedding;

//  When we match this atom, we compute a preference value based on
//  which of our Substructure_Atom_Specifier's match

    int _preference_value_current_match;

//  In addition to preference values derived from the Substructure_Atom_Specifiers,
//  a Substructure_Atom object may also have some preference objects.
//  Preference values from these are added to those which come from the
//  specifiers

    resizable_array_p<Substructure_Atom_Specifier> _preferences;

//  Normally when summing preference values, we find the first preference
//  which matches, and then break. Optionally we can sum _all_ preferences
//  which match

    int _sum_all_preference_hits;

    // The number of unmatched atoms attached.
    // Note this is matched once all atoms have been matched.
    Min_Max_Specifier<int> _unmatched_atoms_attached;

    Substructure_Atom * _parent;
    Substructure_Bond * _bond_to_parent;

//  for component level grouping we can specify an arbitrary fragment id.
//  all root atoms with the same fragment ID must be in the same fragment.
//  for efficiency, the fragment ID is propagated to the root atom

    int _fragment_id;

    // Atoms can be assigned an atom type group. All matched atoms with the
    // same atom type group must have the same atom type.
    // atom_type_group values must be positive.
    unsigned int _atom_type_group;

//  connections to other atoms. First is same as _bond_to_parent. Here
//  is where ring closing bonds in the query will be found.
//  Note that only the second atom knows about a ring closure

    resizable_array_p<Substructure_Bond> _bonds;

    int _match_as_match_or_rejection;

//  Here begins all the stuff associated with the matching process.
//  All the variables below here could be separated out into a separate
//  class from root_substructure_atom.

//  _current_hold_atom  is also the _icon'th connection to atom _anchor.
//  This match is via Substructure_Bond b;

    Target_Atom * _anchor;      // should be _parent->current_hold_atom ();
    int           _con;
    int           _anchor_ncon;

    // Feb 2023. Keep track of whether or not we are on the queue of atoms awaiting
    // matching. This avoids a search.
    int _on_queue_for_matching;

//  End of non Root_Substructure_Atom variables

//  private functions

    void _default_values ();

    int _add_component  (const msi_object & msi);
    int _add_child      (const msi_object * msi, extending_resizable_array<Substructure_Atom *> & completed);

    int _adjust_nrings (int);
    int _adjust_ring_sizes (const List_of_Ring_Sizes &);
    int _set_match_any_aromatic (int convert_all_aromatic_atoms_to_generic_aromatic, 
                                            Molecule & m,
                                            atom_number_t my_atom_number);
    int _add_bond_to_molecule (Molecule & m,
                               const Substructure_Atom * other_atom,
                               const Substructure_Bond * b) const;
    int _add_bond_to_molecule (Molecule & m,
                               const Substructure_Bond * b) const;

    int _set_implicit_hydrogens (Molecule & m, atom_number_t a) const;
    int _fill_min_ncon (Molecule & m, atom_number_t a) const;

    int _create_preference_from_msi_object (const msi_object & msi);
    int _create_preference_from_proto (const SubstructureSearch::SubstructureAtomSpecifier& proto);
    int _create_environment_from_msi_object (const msi_object & ,
                                             extending_resizable_array<Substructure_Atom *> &);
    int _create_environment_from_proto(const SubstructureSearch::SubstructureAtomEnvironment& proto,
                                       extending_resizable_array<Substructure_Atom *> & completed);
    int _construct_from_msi_object (const msi_object &, extending_resizable_array<Substructure_Atom *> &);

    int _parse_smarts_environment (const Atomic_Smarts_Component &);

    int _extract_initial_atom_number (const_IWSubstring & mysmarts);

    int _parse_smiles_specifier (Molecule & m, const MDL_File_Data &, extending_resizable_array<Substructure_Atom *> & completed);
    int _parse_smarts_specifier (const const_IWSubstring & ,
                                 Parse_Smarts_Tmp & pst,
                                 int atoms_in_previous_disconnected_sections,
                                 Smiles_Ring_Status &);

    int _build_atomic_number_matches_from_atom_list(const ISIS_Atom_List & ali);

    int _write_msi (std::ostream &, int &, int = 0);
    int _write_children_msi (std::ostream & os, int & object_id, int indentation);
    int _match_all_ring_ids (Target_Atom *, const Query_Atoms_Matched &);

    void _release_hold ();

//  int _check_ring_id (Target_Atom *);
    int _matches (Target_Atom &, const int *);

    int _evaluate_environment (int * already_matched);
    int _evaluate_environment (Query_Atoms_Matched &, int * already_matched);
    int _evaluate_environments (Target_Atom *, int * already_matched, int atoms_in_target);
    int _evaluate_environments (Target_Atom *);

//  These next functions relate to non-root substructure atoms

    int _ring_closure_bonds_match (Target_Atom * a) const;

    int _add_bond(const msi_object * msi, extending_resizable_array<Substructure_Atom *> & completed);
    int _add_bond(Substructure_Bond * b);
    int _add_bonds(const google::protobuf::RepeatedField<uint32_t>& atoms,
                   const bond_type_t btype,
                   extending_resizable_array<Substructure_Atom *> & completed);
    int _process_smarts_bond(const IWString& input,
                             extending_resizable_array<Substructure_Atom *> & completed);
    int _process_attribute_bond(const msi_attribute *, bond_type_t, extending_resizable_array<Substructure_Atom *> &);
    int _process_attribute_bonds(const SubstructureSearch::SubstructureAtom& proto,
                                 extending_resizable_array<Substructure_Atom *> & completed);
    int _process_attribute_bonds(const msi_object & msi, 
                                 extending_resizable_array<Substructure_Atom *> & completed);
    int _process_attribute_smarts_bond (const msi_attribute * att,
                                        extending_resizable_array<Substructure_Atom *> & completed);
    int _process_substructure_bond(const SubstructureSearch::SubstructureBond& bond,
                                        extending_resizable_array<Substructure_Atom *> & completed);
    int _construct_from_proto(const SubstructureSearch::SubstructureAtom& proto,
                              extending_resizable_array<Substructure_Atom *> & completed);

    int _add_ncon_query_atom_preferences(const Molecule & m, atom_number_t my_atom_number);

    const Substructure_Bond * _bond_between_atoms(int a1, int a2) const;
    const Substructure_Bond * _bond_between_atom_map_numbers(int a1, int a2) const;

  public:
    Substructure_Atom ();
    Substructure_Atom (const Element *);

    ~Substructure_Atom ();

    int ok () const;
    int ok_recursive () const;     // recursively checks all children
    int debug_print (std::ostream &, const IWString &) const;
    int recursive_debug_print (std::ostream &, const IWString & = "") const;
    int terse_details (std::ostream &, const IWString &) const;
    int recursive_terse_details (std::ostream &, const IWString &) const;
    int print_connectivity_graph(std::ostream &) const;
    // Transfer attributes to `proto`.
    int BuildProto(SubstructureSearch::SubstructureAtom& proto) const;

    void set_match_as_match_or_rejection (int s) { _match_as_match_or_rejection = s;}

    int  add_ncon_preference_object (int n, int p);

    int atom_map_number() const { return _atom_map_number;}
    int global_match_id() const { return _global_match_id;}

    int number_descendants () const;

    int preferences_present () const;     // does descendants too
    int set_preference_value (int);

    int fragment_ids_present () const;

    int symmetry_groups_present() const;
    int first_symmetry_group() const;

    int atom_type_groups_present() const;
    unsigned int atom_type_group() const {
      return _atom_type_group;
    }

    int spinach_match_specified () const;

    void adjust_initial_atom_numbers (const int * xref);     // subset of molecule has been used, need to adjust initial atom numbers

    void assign_unique_atom_numbers (int &);
    // Non const method that counts values and sets _attributes_specified.
    int  count_attributes_specified ();
    int  attributes_specified ();
    int  unique_id () const { return _unique_id;}
    int  initial_atom_number () const { return _initial_atom_number;}

    int ring_ids_present() const;
    int fused_system_ids_present() const;

//  when building from a smarts, we can explicitly set the initial atom number
   
    void set_initial_atom_number (int n) { _initial_atom_number = n;}

    Substructure_Atom * query_atom_with_initial_atom_number (atom_number_t);
    Substructure_Atom * query_atom_with_atom_map_number (atom_number_t);

    Substructure_Atom_Specifier * component(const int s) const { return _components[s];}
    int ncomponents() const { return _components.number_elements();}

    int  highest_initial_atom_number () const;
    int  highest_atom_map_number () const;
    void identify_atom_numbers (extending_resizable_array<int> &) const;
    void identify_atom_map_numbers (extending_resizable_array<int> &) const;
    void assign_unique_id_from_atom_number_if_set (extending_resizable_array<int> & numbers_in_use);
    int  assign_atom_map_numbers(int & amap);

    const IWString & text_identifier () const { return _text_identifier;}

    int numeric_value (double &) const;
    void set_numeric_value (double n) { _numeric_value = n;}

//  int set_do_not_match_these_elements (const resizable_array<const Element *> & atn);    // used with ISIS atom lists

    int fragment_id (int &) const;
    void set_fragment_id (int f) { _fragment_id = f;}

    int preference_value_current_match () const { return _preference_value_current_match;}

    int parse_smarts_specifier (const msi_attribute * msi);
    int parse_smarts_specifier (const const_IWSubstring &, Parse_Smarts_Tmp &, int);
    int parse_smarts_specifier (const const_IWSubstring &);

    int parse_smiles_specifier (const msi_attribute * msi);
    int parse_smiles_specifier (const IWString & smiles);

    int construct_from_msi_object (const msi_object &, extending_resizable_array<Substructure_Atom *> &);
    int construct_from_proto (const SubstructureSearch::SubstructureAtom & proto, extending_resizable_array<Substructure_Atom *> & completed);
    int FormBonds(const SubstructureSearch::SubstructureAtom& proto, extending_resizable_array<Substructure_Atom *> & completed);
    int create_from_molecule (Molecule & m, 
                              const MDL_File_Data & mdlfd,
                              atom_number_t my_atom_number,
                              Substructure_Atom * my_parent,
                              const Molecule_to_Query_Specifications &,
                              extending_resizable_array<Substructure_Atom *> & completed,
                              const int * include_these_atoms = nullptr);

    int construct_from_smarts_token (const const_IWSubstring & smarts);
    int construct_from_smarts_token (const char * smarts, int nchars);
    int smarts (IWString &) const;
    int add_bond (const msi_object *, int, int, int, extending_resizable_array<Substructure_Atom *> &);

    int collect_all_atoms (extending_resizable_array<Substructure_Atom *> & completed);

    int check_internal_consistency (int) const;

    int or_id () const { return _or_id;}

    int ring_id () const { return _ring_id;}
    int fused_system_id () const { return _fused_system_id;}

    Atom * create_atom () const;

    int print_current_hold (std::ostream & os) const;

    int min_atoms_for_match () const;

    int atomic_number (atomic_number_t &) const;

//  During a substructure search, we gain efficiencies by only searching
//  those parts of the molecule where a match is possible. Given a
//  Molecule_to_Match object, we can examine our _atomic_number values
//  and return starting and stopping points

    int determine_start_stop (const Molecule_to_Match &, int & istart, int & istop) const;

//  Functions min_ncon and min_nbonds are somewhat unusual, because they are used
//  by create_molecule only

    int min_ncon () const;
    int min_nbonds () const;

//  This next group of functions are used to tell whether or not an
//  attribute has been specified for an atom. Used mostly by the
//  fingerprint creation functions, and create_molecule.

    int ncon_specification (int &) const;    // used by fingerprint routines.
    int hcount_specification (int &) const;
    int formal_charge_specification (formal_charge_t &) const;

    int set_ring_membership (Molecule &);

    int formal_charge (formal_charge_t &) const;
//  int update_aromaticity (aromaticity_type_t);

    int set_hold (Target_Atom *, int *);
    int release_hold (int *);
    int recursive_release_hold ();

    Target_Atom * current_hold_atom () const { return _current_hold_atom;}
    atom_number_t atom_number_matched () const;
    int is_matched () const { return nullptr != _current_hold_atom;}

    int unmatched_connections (const int *) const;

    int include_in_embedding () const { return _include_in_embedding;}
    void set_include_in_embedding(int s) {
      _include_in_embedding = s;
    }

    void set_unique_id (int id);

    int unique_id_from_initial_atom_number ();

    int  add_ring_closure_bond (Substructure_Bond *);
    void notify_extra_child (Substructure_Atom * a);

    int involves_rings () const;
    int nrings () const;     // the smallest number of rings for any component
    int determine_ring_or_non_ring (int &) const;  // nrings value if specified
    int ring_sizes_specified (resizable_array<int> &) const;
    int involves_ring_specifications_for_bonds () const;
    int involves_aromatic_bond_specifications (int &) const;
    int unmatched_atom_counts_present() const;
    int unmatched_atoms_attached_specified() const;

    int adjust_ring_specifications (int, const List_of_Ring_Sizes &);

    int write_msi (std::ostream &, int &, const const_IWSubstring &, int = 0);
    int write_proto (std::ostream &, int &, const const_IWSubstring &, int = 0);

    int matches (Target_Atom &, const int *);

//  Once the atom is embedded, we can check whether or not any chirality is OK

//  int chirality_matched () const;

//  These next functions is used by Substructure_Atom's which are part
//  of an environment

    int matches (int * previously_matched_atoms);
    int environment_search(Query_Atoms_Matched & matched_query_atoms,
                           int * previously_matched_atoms);

    int add_your_children (Query_Atoms_Matched &);
    int remove_your_children (Query_Atoms_Matched &, int *);

    int create_molecule (Molecule &, int, int) const;

//  Next come the functions which are specific to non-root atoms

    Target_Atom * anchor () const { return _anchor;}
    int  set_anchor (Target_Atom *);

    int  prepare_for_matching (Target_Atom *);

    int move_to_next_match_from_current_anchor (int *, const Query_Atoms_Matched &);

    int nbonds () const { return _bonds.number_elements ();}
    const resizable_array_p<Substructure_Bond> & bonds() const { return _bonds;}

    int no_longer_anchored ();

    Substructure_Atom * parent () const { return _parent;}
    const Substructure_Bond * bond_to_parent () const { return _bond_to_parent;}

    int number_children () const { return _children.number_elements ();}
    Substructure_Atom * child (int i) const { return _children[i];}
    void add_child (Substructure_Atom * s) { _children.add(s);}

    int number_ring_closure_bonds() const { return _bonds.number_elements();}
    const Substructure_Bond * ring_closure_bond(int i) const { return _bonds[i];}

    const Substructure_Bond * bond_between_atoms (int, int) const;
    const Substructure_Bond * bond_between_atom_map_numbers(int a1, int a2) const;

    int has_ring_closure_bond_to (const Substructure_Atom *) const;

//  this function is used in the environment

    void set_parent (Substructure_Atom * a, Substructure_Bond * b, int add_bond_to_bonds_array = 0);
    void no_parent () { _parent = nullptr; _bond_to_parent = nullptr;};

//  Used during building from smarts

    int first_chirality_value_in_any_component() const;
    int complete_chiral_centre_via_ring_closure (Substructure_Chiral_Centre & c, int & ndx) const;

    int is_bonded_to (atom_number_t) const;

//  smirks parsing reads the products as a query, and we need these attributes

    const Element * first_specified_element () const;
    int first_specified_isotope () const;
    int first_specified_formal_charge () const;

    int query_atom_with_isotope (int) const;   // first query atom that specifies a given isotope

    template <typename T> int any_query_atom(T) const;

    // If there is a specification of unmatched query atoms, does it match?
    int unmatched_atoms_attached_matches(const Query_Atoms_Matched& matched_query_atoms,
                                         const int * already_matched,
                                         const Molecule_to_Match& target);

    // Used when parsing environments. Note that we take ownership of `spec`.
    int add_component(Substructure_Atom_Specifier* spec);

    int on_queue_for_matching() const {
      return _on_queue_for_matching;
    }
    void set_on_queue_for_matching(int s) {
      _on_queue_for_matching = s;
    }
};

extern std::ostream & operator << (std::ostream &, const Substructure_Atom &);


/*
  An environment owns the attachment points, and the number of hits
  needed.

  The environment consists of any number of possibilities, which are
  applied as an OR operator. These are the components of the resizable_array
*/

class Substructure_Environment : public resizable_array_p<Substructure_Atom>
{
  protected:
    int _unique_id;

//  Our attachment points

    resizable_array<Substructure_Atom *> _possible_parents;

//  and the type of bond by which we might be joined

    Substructure_Bond _bond;

    int _or_id;
    int _and_id;

//  We can also specify the number of times the environment must match
//  Note that this is only useful in cases where multiple attachment points
//  are specified.

    Min_Max_Specifier<int> _hits_needed;

//  In the case of multiple possible attachment points, we often want to
//  be able to say that either the environment is attached at a point,
//  or nothing is there. Say we were looking for chlorobenzenes, but
//  did not want to hit nitrobenzene. Use no_other_substituents_allowed
//  to do that.

    int _no_other_substituents_allowed;

//  With something like O1C(C)(C)CCCC1 PBCHM12395962, should the environment
//  [CH3] match once or twice?

    int _max_environment_matches_per_attachment_point;

//  When multiple environments are present, do they each have their own
//  attchment point, or can different environments share attachment points

    int _environments_can_share_attachment_points;

//  We can limit the number of matches identified. this is mostly useful when
//  there are multiple environments, and they cannot share attachment points

    int _max_matches_to_find;

//  This variable will be 0 if matching the ENTIRE environment should
//  result in a mis-match rather than a match. Note the interesting
//  interplay between this flag and Substructure_Atom::_match_as_match_or_rejection;
//  This flag is tested after the whole environment tree has been matched.

  protected:
    int _query_environment_match_as_rejection;

  private:

//  If is useful to keep track of the number of times each environment
//  component is matched.

    extending_resizable_array<int> _matches;

//  Feb 2016. Need some way of saying H is a valid environment - without having to use explicit H

    int _hydrogen_ok_as_environment_match;

// private functions

    int _print_common_info (std::ostream & os, const IWString & indentation) const;

    int _process_attribute_bond (const msi_attribute * att,
                                            bond_type_t bond_type,
                                            extending_resizable_array<Substructure_Atom *> & completed);
    int _process_attribute_bonds (const msi_object & msi,
                                             extending_resizable_array<Substructure_Atom *> & completed);
    int _process_how_connected(const SubstructureSearch::SubstructureEnvironment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed);
    int _process_attachment_bonds(const SubstructureSearch::EnvironmentAttachment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed);
    int _process_attachment_via_substructure_bond(const SubstructureSearch::EnvironmentAttachment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed);
    int _add_possible_parent (atom_number_t possible_parent,
                                                bond_type_t possible_parent_bond_type,
                                                extending_resizable_array<Substructure_Atom *> & completed);

  public:
    Substructure_Environment ();

    int ok () const;
    int ok_recursive () const;     // recursively checks all children
    int debug_print (std::ostream &, const IWString &) const;
    int recursive_debug_print (std::ostream &, const IWString & = "") const;
    int terse_details (std::ostream &, const IWString &) const;
    int recursive_terse_details (std::ostream &, const IWString &) const;

    int and_id () const { return _and_id;}
    int or_id  () const { return _or_id;}

    void assign_unique_atom_numbers(int &);
    int  count_attributes_specified();

    void set_environments_can_share_attachment_points(int s) { _environments_can_share_attachment_points = s;}

    int involves_aromatic_bond_specifications(int &) const;

    int no_other_substituents_allowed() const { return _no_other_substituents_allowed;}

    int construct_from_msi_object(const msi_object & msi,
                               extending_resizable_array<Substructure_Atom *> & completed,
                               atom_number_t possible_parent = INVALID_ATOM_NUMBER,
                               bond_type_t possible_parent_bond_type = INVALID_BOND_TYPE);
    int construct_from_proto(const SubstructureSearch::SubstructureEnvironment& proto,
                               extending_resizable_array<Substructure_Atom *> & completed,
                               atom_number_t possible_parent = INVALID_ATOM_NUMBER,
                               bond_type_t possible_parent_bond_type = INVALID_BOND_TYPE);

    int write_msi(std::ostream & os, int & object_id, int indentation);
    int BuildProto(SubstructureSearch::SubstructureEnvironment& proto) const;

    int matches(int * previously_matched_atoms,
                const std::unique_ptr<int[]>& matched_by_global_specs,
                int * anchor_atom_available);

    int print_environment_matches(std::ostream &) const;
};

class Single_Substructure_Query;

class MatchedAtomMatch
{
  private:
    int _unique_id;

    // The matched atom numbers that we must match
    resizable_array<int> _atoms;

    // Queries derived from the 'smarts' directive in the proto are
    // divided into positive and netagive matches, depending on whether
    // the smarts starts with ! or not.

    // If any of the _positive_matches queries matches, that is a match.
    resizable_array_p<Single_Substructure_Query> _positive_matches;
    // If any of the _negative_matches queries matches, that is a non match.
    resizable_array_p<Single_Substructure_Query> _negative_matches;

  public:
    MatchedAtomMatch();

    void set_id(int s) {
      _unique_id = s;
    }
    int id() const {
      return _unique_id;
    }

    int ConstructFromProto(const SubstructureSearch::MatchedAtomMatch& proto);
    int BuildProto(SubstructureSearch::MatchedAtomMatch& proto) const;

    int Matches(Query_Atoms_Matched& matched_query_atoms,
                Molecule_to_Match& target,
                int* already_matched);
};


/*
  We can specify how many of a given atom type must match

  Should make this inherit from the Element_Matcher object
*/

class Elements_Needed
{
  private:
    resizable_array<const Element*> _elements;

    enum class Op {
      NOT_SET,
      AND,
      OR,
      XOR
    };
    Op _operator;

    Min_Max_Specifier<int> _hits_needed;

    // private functions;
    int _matches_single_element(const Query_Atoms_Matched& qam) const;
    int _determine_match(const int * found, const int nhits) const;
    int _matches_single_element(Molecule_to_Match & target_molecule) const;

  public:
    Elements_Needed();

    int add_atomic_number(int atomic_number);

    int ok() const;
    int debug_print(std::ostream &, const IWString &) const;

    void set_min(int s) {
      _hits_needed.set_min(s);
    }

    int construct_from_msi_object(const msi_object &);
    int ConstructFromProto(const SubstructureSearch::ElementsNeeded& proto);
    int write_msi(std::ostream & os, const char *, int & object_id, int indentation) const;
    int BuildProto(SubstructureSearch::ElementsNeeded& proto) const;

    int matches(Query_Atoms_Matched & qam) const;
    int matches(Molecule_to_Match & target) const;
};

// We can sometimes speed up searches by doing an initial scan of the
// bond list.
// Note that we deliberately only support atomic numbers as ints, rather than
// using the atomic_symbol_hash_value. While it is an int today, maybe it becomes
// unsigned later. And that need is too obscure.
class RequiredBond {
  private:
    // element and number of connections of first atom.
    int _atomic_number_1;
    int _ncon1;

    bond_type_t _btype;

    // element and number of connections of second atom.
    int _atomic_number_2;
    int _ncon2;

    // How many instances needed for a match, defaults to 1.
    // Note that we do not support zero, could be changed.
    int _min_count;

  // private functions.
    int MatchesSingle(const Molecule& m) const;
    int MatchesDouble(const Molecule& m) const;
    int MatchesTriple(const Molecule& m) const;

  public:
    RequiredBond();

    int ConstructFromProto(const SubstructureSearch::RequiredBond& proto);
    int BuildProto(SubstructureSearch::RequiredBond& proto) const;

    int Matches(const Molecule& m) const;
};

/*
  These next two classes were never implemented
*/

class Substructure_Environment_Match : public Substructure_Environment
{
  private:
  public:
    Substructure_Environment_Match ();
};

class Substructure_Environment_Rejection : public Substructure_Environment
{
  private:
  public:
    Substructure_Environment_Rejection ();
};

class Substructure_Environment_Fuzzy : public Substructure_Environment
{
  private:
    Min_Max_Specifier<int> _distance;    // how many bonds from root atom

    int _path_includes_matched_atoms;

  public:
    int debug_print (std::ostream &, const const_IWSubstring &) const;

    int construct_from_msi_object (const msi_object &);
    int write_msi (std::ostream & os, int & object_id, int indentation) const;

};

/*
  Class for the features common to the ring environments
*/

class Substructure_Ring_Environment : public Substructure_Atom
{
  private:
    int _active;

    Min_Max_Specifier<int> _hits_needed;

  public:
    Substructure_Ring_Environment ();

    int active () const { return _active;}

    int construct_from_msi_object (const msi_object &);
    int write_msi (std::ostream &, int & object_id, const const_IWSubstring &) const;

    int matches (Molecule_to_Match &, const int *);
};

class Substructure_Query;
// A ring or ring system can have any number of Substituent's.
class Substituent {
  private:
    // How many of these must there be attached to a ring.
    Min_Max_Specifier<int> _hits_needed;
    // the atoms in the substituent.
    Min_Max_Specifier<int> _natoms;
    // the number of rings in the substituent.
    Min_Max_Specifier<int> _nrings;
    // The max bond separation from attachment point.
    Min_Max_Specifier<int> _length;

    // Substituent's can set a global_id.
    int _set_global_id;

    // Queries that are applied to the atoms in the substituent.
    resizable_array_p<Substructure_Query> _required;
    resizable_array_p<Substructure_Query> _disqualifying;

    // In order to enable writing a proto that is the same as the
    // input, keep a copy of the smarts.
    resizable_array_p<IWString> _required_smarts;
    resizable_array_p<IWString> _disqualifying_smarts;

    // private functions.
    int OkNrings(Molecule& m, const int* storage, int flag) const;
    int OkLength(Molecule& m, const int* storage, atom_number_t anchor, int flag) const;
    int RunQueries(Molecule_to_Match& target, const int * storage, int flag,
                       int& got_required_match,
                       int& got_rejected_match);
    int MatchesInner(Molecule_to_Match& target, const int * ring_atoms, int * storage);
    void SetTargetGlobalIds(int * storage, int flag, Molecule_to_Match& target) const;
    void FillMatchedAtomsArray(std::unique_ptr<int[]>& matched_by_global_specs,
                        const int matoms,
                        const int * storage,
                        int flag) const;

  public:
    Substituent();

    int ConstructFromProto(const SubstructureSearch::Substituent& proto);
    int BuildProto(SubstructureSearch::Substituent& proto) const;

    // Given a set of matched atoms (which may, or may not) describe
    // a ring, or ring system, does this specification match?
    // `storage` is an array of m.natoms() that each Substituent uses to
    // hold the identities of the atoms it is processing.
    int Matches(Molecule_to_Match& target, const int * ring_atoms, int * storage,
                std::unique_ptr<int[]>& matched_by_global_specs);
};

// In processing inter-ring regions, there are a number of attributes
// that are computed. Keeping them all in a struct avoids having
// too many function arguments.

struct InterRingRegionData {
  // The molecule being examined.
  Molecule* m;
  // The value that will be set for entries in `region`.
  int region_number;
  int * region = nullptr;

  // Counters that are updated as atoms in the region are identified.
  int atoms_in_region = 0;
  // The ring atoms that define this region.
  Set_of_Atoms ring_connections;

  public:
    InterRingRegionData(Molecule* mol, int * r) : m(mol), region(r) {
    }

    // At the start of identifying a new region, set the region number
    // and initialise per-region counters.
    void BeginNewRegion(int number);
};

class InterRingAtoms {
  private:
    // The number of atoms in an inter-ring region.
    Min_Max_Specifier<int> _natoms;

    // The number of ring atoms connected. If 2, then a traditional
    // linker group.
    Min_Max_Specifier<int> _ring_connections;

    // The number of instances of this needed for a match.
    Min_Max_Specifier<int> _hits_needed;

    // The longest bonds_between between the ring atoms defining an inter-ring
    // region. If there are > 2 ring atoms connected, then all inter-ring
    // distances must satisfy this constraint.
    Min_Max_Specifier<int> _length;

    // If specified, all lengths specified here must be matched. This enables
    // specification of AND conditions.
    resizable_array<int> _required_length;

    // An inter-ring region can be defined by global_id values set by ring
    // or ring system specifications.
    resizable_array<int> _defined_by_ring_global_id;

    // InterRingAtoms can set a global_id.
    int _set_global_id;

    // Queries that are applied to the atoms in the region.
    resizable_array_p<Substructure_Query> _required;
    resizable_array_p<Substructure_Query> _disqualifying;

    // In order to enable writing a proto that is the same as the
    // input, keep a copy of the smarts.
    resizable_array_p<IWString> _required_smarts;
    resizable_array_p<IWString> _disqualifying_smarts;

  // private functions.

    // Returns true of the ring_connections in `data` satisfy _length.
    int OkLengths(Molecule& m, InterRingRegionData& data, int region_number);
    int OkRequiredLengths(Molecule& m, InterRingRegionData& data, int region_number);
    int RunQueries(Molecule_to_Match& target, InterRingRegionData& data, int region_number);
    int RunQueriesInner(Molecule_to_Match& subset_target, InterRingRegionData& data, int region_number);

  public:
    InterRingAtoms();

    int global_id() const {
      return _set_global_id;
    }

    int ConstructFromProto(const SubstructureSearch::InterRingAtoms& proto);
    int BuildProto(SubstructureSearch::InterRingAtoms& proto) const;

    // Given a set of matched atoms, does this specification match?
    // The atoms where data.region[i] == region_number are examined.
    int Matches(Molecule_to_Match& target,
                int region_number,
                InterRingRegionData& data,
                std::unique_ptr<int[]>& matched_by_global_specs);
    // After all InterRingAtoms are checked against all possible inter ring
    // regions, we need to check whether the number of matches found is 
    // consistent with any specification of _hits_needed.
    int MatchesNhits(int nhits) const;
};

/*
  We can specify details of either a ring, or a ring system. Those
  specifications have certain elements in common
*/

class Substructure_Ring_Base
{
  protected:

    int _match_as_match_or_rejection;

    Min_Max_Specifier<int> _hits_needed;

    Min_Max_Specifier<int> _attached_heteroatom_count;

    Min_Max_Specifier<int> _heteroatom_count;

//  Mar 2004. Allow run-time definition of heteroatoms

    int * _is_heteroatom;

    Min_Max_Specifier<int> _ncon;

//  When specifying hits_needed, it often makes sense to require
//  that all the hits be in one fragment

    int _all_hits_in_same_fragment;

    int _only_keep_matches_in_largest_fragment;

//  unsaturation is defined as being within the ring

    Min_Max_Specifier<int> _within_ring_unsaturation;

    Min_Max_Specifier<int> _atoms_with_pi_electrons;
    Min_Max_Specifier<int> _fully_saturated_atoms;

    Min_Max_Specifier<int> _largest_number_of_bonds_shared_with_another_ring;

    Min_Max_Specifier<int> _strongly_fused_ring_neighbours;

//  We can also specify a query to hang off the ring. Note that the
//  root Substructure_Atom must match one of the ring atoms

//  We do this as a logical expression: 1a-[OH]||<3c-F
//  Therefore we need an array of minmaxspc objects, an array of Substructure_Atom and a logical expression

    resizable_array_p<Min_Max_Specifier<int> > _environment_numerical_requirement;
    resizable_array_p<Substructure_Atom> _environment_atom;
    IW_Logical_Expression _environment_logexp;

    int _environment_can_match_in_ring_atoms;     // Aug 2014. Need to be able to describe the ring itself.

    // Jan 2022. Optionally set a global match id value for each member of the ring/system.
    int _set_global_id;

    // May 2022. Extend atoms covered by a global id to attached carbonyl and similar.
    bool _ring_extends_to_carbonyl;

    // June 2022. Implement substituent idea.
    resizable_array_p<Substituent> _substituent;

    IWString _comment;

//  private functions

    int _construct_environment(const const_IWSubstring &);

  protected:
    int ExtendToCarbonyl(const Molecule& m, int* in_system) const;

  public:
    Substructure_Ring_Base();
    ~Substructure_Ring_Base();

    int ok() const;

    int set_global_id() const {
      return _set_global_id;
    }

  protected:
    int debug_print(std::ostream &) const;

    int write_msi_attributes(std::ostream & os, int & object_id, const const_IWSubstring & ind) const;
    int construct_from_msi_object(const msi_object &, int &);
    int ConstructFromProto(const SubstructureSearch::SubstructureRingBase& proto);
    int BuildProto(SubstructureSearch::SubstructureRingBase& proto) const;

    int _environment_matches(Molecule_to_Match &, int *);
    int _environment_matches(Molecule_to_Match & target,
                             int * ring,
                             int * already_matched);
    int _environment_matches(Molecule_to_Match &, const int *, Substructure_Atom &);
    int _environment_matches(Molecule_to_Match & target,
                             const int * ring,
                             Substructure_Atom & root_atom,
                             int * already_matched);
    int _environment_matches(Molecule_to_Match & target,
                             Substructure_Atom & root_atom,
                             const int * ring,
                             int * already_matched);
    int _environment_matches(Molecule_to_Match & target,
                             Query_Atoms_Matched & matched_query_atoms,
                             int * previously_matched_atoms);

};

/*
  There can be any number of specifications about a ring.

  These must all be true unless an or operator is present

  Note however that saying that we want a 6 memberd aromatic
  carbocycle means that somewhere in the molecule is one of
  these rings. It does NOT mean that every ring must match
  that condition
*/

class Substructure_Ring_Specification : public Substructure_Ring_Base
{
  private:
    Min_Max_Specifier<int> _ring_size;

    int _aromatic;

    Min_Max_Specifier<int> _fused;

    Min_Max_Specifier<int> _fused_aromatic_neighbours;

    Min_Max_Specifier<int> _fused_non_aromatic_neighbours;

    int _spiro_fusion_count;

//  private functions

    int _environment_matches(Molecule_to_Match &, const Ring &);

  // private functions.
    int OkSpiroFusionCount(Molecule& m, const Ring& ring) const;

  public:
    Substructure_Ring_Specification();
    ~Substructure_Ring_Specification();

    int ok() const;
    int debug_print(std::ostream &, const IWString &) const;
    int terse_details(std::ostream &, const IWString &) const;

    int construct_from_msi_object(const msi_object &);
    int ConstructFromProto(const SubstructureSearch::SubstructureRingSpecification& proto);
    int BuildProto(SubstructureSearch::SubstructureRingSpecification& proto) const;
    int write_msi(std::ostream & os, int & object_id, int indentation) const;

    int matches(Molecule_to_Match &, std::unique_ptr<int[]>& matched_by_global_specs);
};

/*
  Within a Substructure_Ring_System_Specification we can specify how many
  of individual ring sizes are needed.
  The ring size can be one or more values as well as a range.
  But the overall RingSizeRequirement matches, with no requirement that any
  individual values specified will be matched. If you need to make sure
  there are two six membered rings and one five membered ring, that will
  take two RingSizeRequirement.
*/

class RingSizeRequirement
{
  private:
    // The size of the ring.
    Min_Max_Specifier<int> _ring_size;

    // For any/all ring size(s) that satisifies _ring_size, the number needed.
    Min_Max_Specifier<int> _count;

  public:
    RingSizeRequirement();

    //int ring_size() const { return _ring_size;}

    int ConstructFromProto(const SubstructureSearch::RingSizeRequirement&);
    int BuildProto(SubstructureSearch::RingSizeRequirement& proto) const;

    // Given a set of ring sizes that are in a ring system, return whether
    // or not we match that state.
    int Matches(const int * ring_sizes_encountered,
                int max_ring_size) const;

    // Does this specification apply to rings of size `ring_size`.
    int MatchRingSize(int ring_size) const {
      return _ring_size.matches(ring_size);
    }

    // Do we match `count` occurrences. Presumably the caller will first
    // call MatchRingSize to determine if this specifiation applies to
    // the ring size under consideration.
    int MatchesCount(const int count) const {
      return _count.matches(count);
    }
};

class Substructure_Ring_System_Specification : public Substructure_Ring_Base
{
  private:
//  We could for example specify that we are looking for a ring
//  system with between 3 and 5 rings (_rings_in_system), all of
//  which are size 5 or 6 (_ring_sizes)

    Min_Max_Specifier<int> _rings_in_system;
    Min_Max_Specifier<int> _ring_sizes;

//  Dec 2003. As currently implemented, every ring in a ring system must
//  match one of the ring_sizes specified in _ring_sizes. Make that flexible

    Min_Max_Specifier<int> _rings_that_must_match_ring_sizes;

    // Jul 2019. Allow specifying a numerical requirement for any ring size.

    resizable_array_p<RingSizeRequirement> _ring_size_requirement;

    Min_Max_Specifier<int> _aromatic_ring_count;
    Min_Max_Specifier<int> _non_aromatic_ring_count;

//  As we examing a ring system we look for the ring with the largest number
//  of fused neighbours. That is the "degree of fusion";

    Min_Max_Specifier<int> _degree_of_fusion;

    Min_Max_Specifier<int> _atoms_in_system;

//  We may want to impose conditions on the "spinach" hanging off the ring. 

    Min_Max_Specifier<int> _number_spinach_groups;

//  We can identity "terminal" rings as those with just one non-spinach connection

    Min_Max_Specifier<int> _number_non_spinach_groups;

//  We match if any one of the "spinach" groups match the specification

    Min_Max_Specifier<int> _atoms_in_spinach_group;

    // If set, every spinach group must match `_atoms_in_spinach_group`.
    // Defaults to false since that is the historical behaviour.
    bool _every_group_matches_atoms_in_spinach_group;

//  Or the longest distance of a spinach atom to the attachment point

    Min_Max_Specifier<int> _length_of_spinach_group;

    // If set, every spinach group must match `_length_of_spinach_group`.
    // Defaults to false since that is the historical behaviour.
    int _every_group_matches_length_of_spinach_group;

//  Once we have spinach identified, atoms that aren't in a ring, and aren't spinach
//  are atoms between ring systems

    Min_Max_Specifier<int> _distance_to_another_ring;

//  May 2004. Realised that even with all the various things I had about
//  fusion and bonds shared with other rings, there was nothing to just
//  count the number of rings involved in a strongly fused relationship

    Min_Max_Specifier<int> _strongly_fused_ring_count;

//  If we have certain conditions set, we may need a per-atom array duing searches

    int _need_per_atom_array;

    // Optionally, ring systems can span spiro fusions.
    bool _ring_systems_extend_across_spiro_fusions;

//  private functions

    int _matches(Molecule_to_Match &, atom_number_t *, std::unique_ptr<int[]>& matched_by_global_specs);

    int _spinach_matches(Molecule_to_Match & target, const atom_number_t *) const;
    int _check_length_of_spinach(const Molecule & m, const int * in_system, atom_number_t atom_in_ring, atom_number_t first_spinach_atom) const;
    int _match_distance_to_another_ring(Molecule_to_Match & target, const int * in_ring_system) const;

    int _check_heteroatoms(const int * tmp,
                    Molecule_to_Match & target) const;

    int _check_ncon(const int * tmp,
                    Molecule_to_Match & target) const;

    int _ring_size_requirements_matched(const int * ring_sizes_encountered,
                int largest_ring_size);

    int GetRingSize(const SubstructureSearch::SubstructureRingSystemSpecification& proto);

    int IdentifyRingsInSystem(
        Molecule_to_Match& target,
        int* ring_already_done,
        resizable_array<int>& rings_in_system,
        int * atoms_in_system) const;

  public:
    Substructure_Ring_System_Specification();
    ~Substructure_Ring_System_Specification();

    int ok() const;
    int debug_print(std::ostream &, const IWString &) const;
    int terse_details(std::ostream &, const IWString &) const;

    int construct_from_msi_object(const msi_object &);
    int ConstructFromProto(const SubstructureSearch::SubstructureRingSystemSpecification& proto);
    int BuildProto(SubstructureSearch::SubstructureRingSystemSpecification& proto) const;
    int write_msi(std::ostream & os, int & object_id, int indentation) const;

    int matches(Molecule_to_Match &, std::unique_ptr<int[]>& matched_by_global_specs);
};

/*
  A common task is to have a link atom enumerate possibilities. For that
  it must know the atom in the target molecule to attach a1 and the
  atom to which to attach a2. It also needs to know how many atoms it
  has placed
*/

class Link_Atom;

class Link_Atom_Current_State
{
  private:
    atom_number_t _lhs;
    atom_number_t _rhs;
    Set_of_Atoms  _placed;

//  _bt will be a single bonding type used for constructing molecules,
//  it will not be a query type

    bond_type_t _bt;

  public:
    Link_Atom_Current_State();
    Link_Atom_Current_State(atom_number_t a1, atom_number_t a2) : _lhs(a1), _rhs(a2) {};

    int initialise(const Link_Atom *);

    atom_number_t lhs() const { return _lhs;}
    atom_number_t rhs() const { return _rhs;}

    bond_type_t btype_for_molecule() const { return _bt;}

    void reset() {_placed.resize_keep_storage(0);}

    atom_number_t last_atom_placed() const { return _placed.last_item();}

    int number_placed() const { return _placed.number_elements();}

    void add(atom_number_t a) { _placed.add(a);} 
};

/*
  A Link atom is just like no-matched_atoms_between, except that it has
  limits on the number of bonds involved.
  Also, there can be a specification of what atom(s) can be within the linked
  region
  The bond involved must be fully specified - different bonding types,
  bond topology, etc...
*/

class Link_Atom
{
  protected:
    int _a1;
    bond_type_t _bt;
    int _a2;

    int _bond_topology;

    Min_Max_Specifier<int> _distance;

//  the symbol could be almost anything. C, N, O, [OD2], *, A, Q, etc

    IWString _symbol;

//  Quite often, this will be a simple atomic symbol, so if we can convert it to
//  an element, we do that

    const Element * _e;

//  If we came from an MDL file, we may have some MDL atom data

    MDL_Atom_Data * _mdl_atom_data;

  public:
    Link_Atom();
    Link_Atom(const Link_Atom &);
    ~Link_Atom();

    int debug_print(std::ostream &) const;

    void set_a1(int a) { _a1 = a;}
    void set_a2(int a) { _a2 = a;}

    int a1() const { return _a1;}
    int a2() const { return _a2;}

    bond_type_t btype() const { return _bt;}
    void set_bond_type(bond_type_t b) { _bt = b;}

    int set_symbol(const const_IWSubstring & s);

    int set_mdl_atom_data(const MDL_Atom_Data * s);

    void set_bond_topology(int s) { _bond_topology = s;}

    int max_distance() const;

    int construct_from_msi_attribute(const msi_attribute *);
    int ConstructFromProto(const SubstructureSearch::LinkAtoms & proto);
    int write_msi(std::ostream &, const IWString &, const char *) const;
    int initialise_from_mdl_record(const const_IWSubstring & buffer, int matoms, int & a);
    int write_M_LIN(atom_number_t, std::ostream &) const;
    int BuildProto(SubstructureSearch::LinkAtoms& proto) const;

    int initialise_from_smarts(const const_IWSubstring &);

    int condition_satisfied(Molecule &, atom_number_t, atom_number_t) const;

//  void set_no_matched_atoms_between(int s) { _no_matched_atoms_between = s;}
//  int no_matched_atoms_between() const { return _no_matched_atoms_between;}

    int satisfies_separation(int d) const { return _distance.matches(d);}

    int swap_atoms(atom_number_t, atom_number_t);

    int atom_is_being_removed(atom_number_t);

    int create_next_variant(MDL_Molecule &, Link_Atom_Current_State &) const;
};

extern std::ostream & operator << (std::ostream &, const Link_Atom &);

/*
  As read in, an ISIS atom list has four pieces of information. Three
  of them are common with the Link_Atom object. The extra one is the
  actual link atom, which will be not appear explicitly in the query
  constructed
*/

class ISIS_Link_Atom : public Link_Atom
{
  private:
    int _a;

// private functions

    int _adjust_atom_number(const int * xref, atom_number_t & a);

    void _default_values();

  public:
    ISIS_Link_Atom();
    ISIS_Link_Atom(const ISIS_Link_Atom &);

    int debug_print(std::ostream &) const;

    int write_M_LIN(std::ostream &) const;

    int central_atom() const { return _a;}

    int construct_from_M_ISIS_record(const IWString &);

    int adjust_atom_numbers(const int *);

    int swap_atoms(atom_number_t, atom_number_t);
};

extern std::ostream & operator << (std::ostream &, const ISIS_Link_Atom &);

class SetOfEmbeddings : public resizable_array_p<Set_of_Atoms> {
  private:
  public:
};
/*
  We need an object to hold the results of a substructure query, as well as
  the intermediate working data used by a substructure query. This is for
  thread safety and const'ness.
  Unfortunately the object holds a mixture of things that are relevent
  while matches are being found mixed in with the actual results.
*/

class Substructure_Results
{
  private:
    uint64_t _hits_found;

    int _atoms_in_target_molecule;

//  When checking for unique embeddings, we need a temporary array

    int * _just_matched;

//  When checking if embeddings overlap, we accumulate the atoms hit

    int * _already_matched;

//  We copy this from the substructure query which is using us

    int _save_matched_atoms;

//  During a search, each query keeps a list of the embeddings it has made. 

    SetOfEmbeddings _embedding;

//  Some programmes (fragmenter for example) need to know the identity
//  of the Substructure_Atom matched to each atom in a molecule, so we
//  keep track of embeddings that way.

    resizable_array_p<Query_Atoms_Matched> _query_atoms_matched;

//  of course, saving the query atoms matched is optional

    int _save_query_atoms_matched;

//  things are easier in the substructure searching code if we can have a flag
//  telling if matching is done - got the required number of hits for example

    int _complete;

//  If we are not perceiving symmetrically equivalent atoms, we need
//  to know the symmetry class of each atom in the target molecule

    const int * _symmetry_class;

    int _embeddings_violating_distance_constraints;

//  Some queries are constrained by needing the matches to be in the same
//  or different fragments

    extending_resizable_array<int> _hits_per_fragment;

//  When we have atoms that are excluded from embeddings, we can get embeddings that
//  contain INVALID_ATOM_NUMBER. We can optionally remove those as we get the embeddings

    int _remove_invalid_atom_numbers_from_new_embeddings;

    int _max_query_atoms_matched_in_search;

//  private functions

    int _are_symmetry_related(const Set_of_Atoms & e1, const Set_of_Atoms & e2) const;

    int _remove_hits_violating_distance_check_all_atoms(Molecule & m,
                              const Min_Max_Specifier<int> & distance_between_hits, int ncheck);

    int _remove_hits_not_in_largest_fragment_multiple_largest(Molecule & m);

  public:
    Substructure_Results();
    ~Substructure_Results();

    int ok() const;

    void reset();            // clears out previous results
    void initialise(int);    // called before a search

    int copy_embeddings(const Substructure_Results & rhs);   // pretty close to operator=
    int add_embeddings(const Substructure_Results & rhs);    // just does the embeddings, nothing else changed

    int set_save_matched_atoms(int s);
    int set_save_query_atoms_matched(int s);

    int save_matched_atoms() const { return _save_matched_atoms;}
    int save_query_atoms_matched() const { return _save_query_atoms_matched;}

    void set_remove_invalid_atom_numbers_from_new_embeddings(int s) { _remove_invalid_atom_numbers_from_new_embeddings = s;}
    unsigned int hits_found() const { return _hits_found;}

    int return_code() const;

    int matching_complete() const { return _complete;}
    void set_complete() { _complete = 1;}

    void matched_this_many_query_atoms(int m) { if (m > _max_query_atoms_matched_in_search) _max_query_atoms_matched_in_search = m;}
    int  max_query_atoms_matched_in_search() const { return _max_query_atoms_matched_in_search;}

    void set_symmetry_class(const int * s) { _symmetry_class = s;}

    int embedding_is_unique(const Set_of_Atoms & new_embedding);
    int embedding_overlaps_previous_embeddings(const Set_of_Atoms & new_embedding);

    int embedding_is_symmetry_related(const Set_of_Atoms & embedding) const;

//  When the matched atoms are not being saved, getting a hit is just a matter
//  of incrementing the counter.

    void got_embedding() { _hits_found++;}

    int add_embedding(Set_of_Atoms * embedding, const Query_Atoms_Matched & qam);

    int sort_by_preference_value();
    int sort_matches(Molecule_to_Match & target, int sort_specification);

    int remove_hits_violating_distance(Molecule_to_Match & target,
              const Min_Max_Specifier<int> & distance_between_hits, int ncheck);

    int overlaps_with_existing_embedding(const Set_of_Atoms & s);     // used by embeddings_do_not_overlap

    int remove_hits_not_in_largest_fragment(Molecule_to_Match & target);

    int remove_low_preference_hits();

    int remove_embedding(int);

    const Set_of_Atoms * embedding(int i) const { return _embedding[i];}
    const Query_Atoms_Matched * query_atoms_matching(int i) const { return _query_atoms_matched[i];}

    uint32_t number_embeddings() const { return _embedding.size();};
    int print_embeddings(std::ostream &, int = 0) const;
    int print_embeddings(std::ostream &, const Molecule *) const;

    void size_hits_per_fragment_array(int s);
    void got_hit_in_fragment(int f) { _hits_per_fragment[f]++;}
    const extending_resizable_array<int> & hits_per_fragment() const { return _hits_per_fragment;}

//  A common operation is to mark each atom hit by a substructure query

    void each_embedding_set_vector(int * v, int s) const;
    // A version that is used for python bindings.
    std::vector<int> EachEmbeddingSetVector(int matoms, int value) const;

//  We may need to re-order the embeddings - trxn. The strategy would be for the caller to
//  get all the embeddings, then pass them back 

    int set_embeddings(const Set_of_Atoms **, const int n);

    const SetOfEmbeddings& embeddings() const {
      return _embedding;
    }

    // If multiple embeddings, compress them to just one.
    int CompressToSingleEmbedding();

    // If per-atom global id values have been set, remove embeddings whose matched query
    // atoms do not match those global match conditions.
    int RemoveEmbeddingsNotSatisfyingGlobalId(const int * matched_by_global_conditions);
};

// Specify the bond separation between pairs of matched atoms.
// This is likely to be most useful in the case of a query that consists
// of multiple fragments, but can be applied to any query.
class SeparatedAtoms {
  private:
    // A pair of matched atoms - indices within the embedding.
    int _a1;
    int _a2;
    // The bonds_between value that must be met.
    Min_Max_Specifier<int> _separation;

  public:
    SeparatedAtoms();

    int Build(const SubstructureSearch::SeparatedAtoms& proto);
    int BuildProto(SubstructureSearch::SeparatedAtoms& proto) const;

    // Return true if atoms molecule[embedding[separated_atoms.a1]]
    //                  and molecule[embedding[aseparated_atoms.2]]
    // are compatible with the bond separation constraint in
    // separated_atoms.separation.
    int Matches(Molecule& m,
                const Set_of_Atoms& embedding) const;
};

// See the proto for description of what this is.
class Region {
  private:
    // The two or more matched atoms that define the region.
    resizable_array<int> _matched_atom;

    // The number of atoms that are present in the region.
    // This includes the matched atoms that define it.
    Min_Max_Specifier<int> _natoms;

    // The number of rings fully contained within the region
    Min_Max_Specifier<int> _nrings;

    Min_Max_Specifier<int> _atoms_not_on_shortest_path;

  public:
    int ConstructFromProto(const SubstructureSearch::Region& proto);
    int BuildProto(SubstructureSearch::Region& proto) const;

    // Given an embedding, return true if we match.
    int Matches(Molecule_to_Match& target,
                Query_Atoms_Matched & matched_atoms,
                const int* already_matched) const;
};

// Before atom matching is attempted, we can quickly check a number
// whole molecule properties.
// For example if the query specifies atoms in two different aromatic
// rings, we can restrict searches to just  molecules having at least
// two aromatic rings.
class RequiredMolecularProperties {
  private:
    // When doing matches, if nothing is specified, quick exit.
    int _attributes_specified;

    // The number of atoms in the molecule.
    iwmatcher::Matcher<int> _natoms;

    // The number of rings in the molecule.
    iwmatcher::Matcher<int> _nrings;

    // The number of heteroatoms in the molecule.
    iwmatcher::Matcher<int> _heteroatoms_in_molecule;

    // The number of fused rings in the molecule.
    // Note that this can never be 1.
    iwmatcher::Matcher<int> _fused_rings;

    // The number of strongly fused rings in the molecule.
    // Note that this can never be 1.
    iwmatcher::Matcher<int> _strongly_fused_rings;

    // The number of isolated rings in the molecule.
    iwmatcher::Matcher<int> _isolated_rings;

    // The number of ring systems in the molecule.
    // Note that by default, ring fusions do not span spiro fusions.
    iwmatcher::Matcher<int> _ring_systems;

    // The number of aromatic rings in the molecule.
    iwmatcher::Matcher<int> _aromatic_rings;

    // The number of non aromatic rings in the molecule.
    iwmatcher::Matcher<int> _non_aromatic_rings;

    // The number of aromatic atoms in the molecule.
    iwmatcher::Matcher<int> _aromatic_atoms;

    // The number of isotopic atoms in the molecule.
    iwmatcher::Matcher<int> _number_isotopic_atoms;

    // The number of fragments in the molecule.
    iwmatcher::Matcher<int> _number_fragments;

    // The number of atoms in the spinach (non scaffold) in the molecule.
    iwmatcher::Matcher<int> _atoms_in_spinach;

    // The number of atoms in between rings in the molecule.
    iwmatcher::Matcher<int> _inter_ring_atoms;

    // The net formal charge of the molecule.
    iwmatcher::Matcher<int> _net_formal_charge;

    // Requirements on counts of certain elements.
    resizable_array_p<Elements_Needed> _elements_needed;

    // A set of bond specifications that are checked before any
    // atom matching is done.
    resizable_array_p<RequiredBond> _required_bonds;

    // The most common use for net formal charge is to look for molecules
    // that have a net formal charge, either positive or negative, so make
    // that easy to specify.
    int _any_net_formal_charge;

    // If further conditions are added they must be accounted for in at 
    // least two places.
    //    Matches
    //    DiscernAttributesSpecified

  // private functions
    int ParseElementsNeeded(const SubstructureSearch::RequiredMolecularProperties& proto);

    int MatchElementsNeeded(Molecule_to_Match & target_molecule) const;
  
    int DiscernAttributesSpecified();

    int SpinachAtomsMatch(Molecule_to_Match & target) const;

  public:
    RequiredMolecularProperties();

    int ConstructFromProto(const SubstructureSearch::RequiredMolecularProperties& proto);
    int BuildProto(SubstructureSearch::RequiredMolecularProperties& proto) const;

    int attributes_specified() const {
      return _attributes_specified;
    }

    int Matches(Molecule_to_Match & target_molecule);
};

// Specify that within some distance of a matched atom, there must be
// another query matched.
class NearbyAtoms {
  private:
    resizable_array_p<Substructure_Query> _query;

    // can the atoms matched here overlap with the matoms matched by the main query.
    boolean _can_overlap_matched_atoms;

    // How many instances of this do we require.
    Min_Max_Specifier<int> _hits_needed;

    // By default, we measure ourselves relative to any matched atom, but
    // if this is set, we measure relative to only some matched atoms.
    resizable_array<uint32_t> _matched_atom;

    // The constraint on the number of bonds between our matched atom(s)
    // and those in the query.
    Min_Max_Specifier<int> _bonds_between;

    // By default, we are a positive requirement for a match. We can instead
    // be a rejection.
    int _rejection;

  // private functions.
    int AnyMatchedAtomInRange(Molecule_to_Match& target,
                const int* already_matched, const Set_of_Atoms& embedding) const;
    int SpecifiedmatchedAtomInRange(Molecule_to_Match& target, const Query_Atoms_Matched& qam,
                const Set_of_Atoms& embedding) const;

  public:
    NearbyAtoms();

    int ConstructFromProto(const SubstructureSearch::NearbyAtoms& proto);
    int BuildProto(SubstructureSearch::NearbyAtoms& proto) const;

    int Matches(Molecule_to_Match& target, Query_Atoms_Matched& qam,
                const int* already_matched);
};


/*
  The process of verifying the internal consistency of a substructure
  query is somewhat expensive. To avoid recomputing this, we keep
  a record of the "atom" count of the query when consistency was 
  last done. If the count is different, we need to compute consistency.
  Note that this can be defeated by someone directly manipulating
  any of the Substructure_ objects which comprise the query.
*/

class Substructure_Query;

class Single_Substructure_Query
{
  friend
    std::ostream &
      operator << (std::ostream &, const Single_Substructure_Query &);

  private:

    magic_number_t _magic;

//  In order to allow multiple queries within the same molecule, we allow
//  for multiple root atoms. Note that in this case, the queries cannot
//  overlap, whereas OR queries in a composite query can overlap

    resizable_array_p<Substructure_Atom> _root_atoms;

//  When doing matches we need to keep track of which root atom is being matched
    
    int _iroot;

//  We can specify the number of bonds attached to the matched atoms
//  Note that this refers only to atom matches, NOT to global conditions

    Min_Max_Specifier<int> _ncon;

//  The environment of the match

    resizable_array_p<Substructure_Environment> _environment;
    int _environment_must_match_unmatched_atoms;

//  Environments which cause a non-match are handled similarly

    resizable_array_p<Substructure_Environment> _environment_rejections;

//  It is useful to know how many matches were obtained before the
//  various environment things happened.

    int _matches_before_checking_environment;

//  the number of matches rejected by no match to environment

    int _no_match_to_environment;

    int _match_to_environemt_rejection;

//  Various whole molecule conditions for a match

    RequiredMolecularProperties _required_molecular_properties;

    resizable_array_p<Substructure_Ring_Specification> _ring_specification;

//  Apr 03. Why not link _ring_specification's via a logical expression

    IW_Logical_Expression _ring_specification_logexp;

    resizable_array_p<Substructure_Ring_System_Specification> _ring_system_specification;

    IW_Logical_Expression _ring_system_specification_logexp;

//  when testing _ncon, we can optionally ignore singly connected
//  substituents. This was first introduced for the bis-napthylene
//  query, where I wanted only one point of attachment, but also 
//  wanted to allow the possibility of a CH3 or F substituent elsewhere
//  on the ring

    int _ncon_ignore_singly_connected;

//  When specifying hits_needed, it often makes sense to require
//  that all the hits be in one fragment

    int _all_hits_in_same_fragment;

//  When specifying an attached heteroatom count, one can also specify
//  the atoms which will be considered "heteroatoms"

    Min_Max_Specifier<int> _attached_heteroatom_count;
    resizable_array<atomic_number_t> _heteroatoms;

//  Sometimes it is useful to put constraints on the number of ring
//  atoms which are matched

    Min_Max_Specifier<int> _ring_atoms_matched;

// We can also specify that there must be a distance between any hits
// We first find all the embeddings, and then reject those which are
// outside the constraint relative to any higher earlier match

    Min_Max_Specifier<int> _distance_between_hits;

//  Jan 2003. By default, every atom in every embedding is checked. If we
//  know about the query, maybe we just need to check some limited
//  number of atoms

    int _matched_atoms_to_check_for_hits_too_close;

//  Sometimes we just want to reject molecules where the matches are too close

    int _fail_if_embeddings_too_close;

    // START OF WHAT IS NOW IN REQUIRED_MOLECULAR_PROPERTIES
//  We can only check molecules that have a given number of fragments
    Min_Max_Specifier<int> _number_fragments;

//  We can set constraints on the number of spinach atoms
    Min_Max_Specifier<int> _atoms_in_spinach;

//  The atoms between rings
    Min_Max_Specifier<int> _inter_ring_atoms;

//  We can also put constraints on the number of heteroatoms in the molecule.
//  These will depend on the definition of heteroatoms stored in _heteroatoms

    // Only attempt matching if this may heteroatoms in the molecule.
    Min_Max_Specifier<int> _heteroatoms_in_molecule;

    // Only attempt matching if the net formal charge of the molecule matches.
    Min_Max_Specifier<int> _net_formal_charge;

    Min_Max_Specifier<int> _natoms;
    Min_Max_Specifier<int> _nrings;

//  Simple counts of the number of each type of ring

    Min_Max_Specifier<int> _aromatic_rings;
    Min_Max_Specifier<int> _aromatic_atoms;
    Min_Max_Specifier<int> _non_aromatic_rings;
    Min_Max_Specifier<int> _fused_rings;
    Min_Max_Specifier<int> _strongly_fused_rings;
    Min_Max_Specifier<int> _isolated_rings;

    Min_Max_Specifier<int> _number_isotopic_atoms;

//  Nov 98. Needed something to get molecules with two isolated rings or
//  ring systems.

    Min_Max_Specifier<int> _isolated_ring_objects;

    // Requirements on counts of certain elements.
    resizable_array_p<Elements_Needed> _elements_needed;

    // END OF WHAT IS NOW IN REQUIRED_MOLECULAR_PROPERTIES

//  Jan 2005. We may have constraints on the number of unmatched atoms

    Min_Max_Specifier<int> _unmatched_atoms;

    // Within each embedding, what fraction of the atoms in the
    // molecule are matched.
    float _min_fraction_atoms_matched;
    float _max_fraction_atoms_matched;

    // Oct 2022. Across all embeddings, what is the fraction of
    // atoms in the molecule that are hit by one or more embeddings.
    float _min_all_matches_fraction_atoms_matched;
    float _max_all_matches_fraction_atoms_matched;

    // If atom types are specified this will be active.
    // Incorporate as pointer to lessen the dependency on this.
    // Not sure it is a good idea yet.
    Atom_Typing_Specification* _atom_typing;

    // May 2022.
    // Given 2 matched atoms, we can count the number of atoms in between them.
    // Not just the shortest path, but all atoms that are reachable
    // starting at one and arriving at the other.
    // Note this will be confounded if atoms are in the same ring system.

    Min_Max_Specifier<int> _atoms_in_region;

    protected:

    int _embeddings_violating_distance_constraints;

//  Dec 2008. I want to be able to specify whether a matched atom is in the
//  molecular spinach or not. We query our atoms at the beginning of the
//  search if spinach requirements are present, we check the matched atoms
//  whenever we get a match

    int _spinach_match_requirements_present;

    private:

//  When I implemented _distance_between_hits, I quickly ran into the
//  case where I needed to prioritise the hits. The default behaviour
//  is to remove all embeddings having less than the highest preference
//  value. When _distance_between_hits is active, it may be preferable
//  to just sort the embeddings by their preference value, and then
//  scan for lower preference value embeddings which violate a distance
//  constraint relative to a higher priority embedding

    int _sort_by_preference_value;

//  A search can be restricted in several ways.

    int _find_one_embedding_per_start_atom;
    int _find_unique_embeddings_only;

//  We can insist on a range of hits in order to get a match

    Min_Max_Specifier<int> _hits_needed;

//  We can also specify that multiple embeddings have no atoms in
//  common. Note that this will imply _find_one_embedding_per_start_atom

    int _embeddings_do_not_overlap;

//  When doing a query in which min_hits_needed is specified, we can
//  optionally modify our return code to reflect any minimum.

    int _normalise_rc_per_hits_needed;

//  We can also subtract an offset from the number of hits found.
//  This was first used for the cyclohexane demerit. We wanted to, in effect,
//  determine the number of cyclohexanes beyond 2

    int _subtract_from_rc;

//  During matching the order may or may not be the same as what the user
//  entered. We can respect the initial numbering if needed

    int _respect_initial_atom_numbering;

//  and when we are respecting the initial atom numbering, we need to know
//  what is the largest atom number present

    int _highest_initial_atom_number;

//  We can specify conditions on the number of heteroatoms matched

    Min_Max_Specifier<int> _heteroatoms_matched;

//  We can specify arbitrary conditions on the number of given
//  atom types which must match

    resizable_array_p<Elements_Needed> _element_hits_needed;


    // A set of bond specifications that are checked before any
    // atom matching is done.
    resizable_array_p<RequiredBond> _required_bonds;

//  We can also stop looking once we have a given number of matches.

    int _max_matches_to_find;

//  We may or may not want to save the atom numbers of the matches,
//  sometimes we just want the number of hits. Note that if we want
//  unique embeddings only, then we must record the matched atoms

    int _save_matched_atoms;

    int _rings_in_query;
    int _consistency_count;    // atom count when last consistency check done

//  Since a query may contains OR conditions, there is not a fixed number
//  of atoms in a query. This number is the smallest possible number of
//  atoms needed for a full match.

    int _min_atoms_in_query;
    int _max_atoms_in_query;

    IWString _comment;

//  sometimes it is convenient to have a numeric value associated with
//  one of the queries.

    resizable_array<double> _numeric_value;

//  If any of our bonds specify anything about rings, or aromaticity
//  we will need to ask the molecule in advance to compute ring membership

    int _need_to_compute_ring_membership;

//  Similarly, if any of the bonds specify aromaticity, we need to ask
//  the molecule first to do the computation - Sept 2000, not really needed
//  any more, I now compute aromaticity all the time

    int _need_to_compute_aromaticity;

//  If any of out Substructure_Atom's have an atom type group, those must
//  be reconciled at the end of matching.
    int _need_to_check_atom_type_groups;

//  We also keep track of the largest number of atoms matched during a substructure search

    iwmax<int> _max_query_atoms_matched;

//  Does the caller want us to perceive symmetry

    int _do_not_perceive_symmetry_equivalent_matches;

//  If any of our Substructure_Atoms have preference values, we will
//  need to sort and classify the hits whenever we have more than one
//  embedding.

    int _preferences_present;

//  for efficiency we store whether or not fragment_id's are present

    int _fragment_ids_present;

//  We can impose constraints on the distance between root atoms

    Min_Max_Specifier<int> _distance_between_root_atoms;

//  In cases where we have multiple root atoms we can specify that there be no
//  matched atoms on any path between matched root atoms. We need an object to
//  describe two matched atoms. For laziness we use the Bond object, just because
//  it already has two numbers which cannot be identical. We never use the bond
//  type or other attribute of a bond.
//  NOTE that the current implementation does not pay attention to the length of
//  the path between the atoms. The match succeeds if ANY path of unmatched atoms
//  can be found between the pairs of atoms. Only a problem with matched atoms in
//  ring systems

    resizable_array_p<Bond> _no_matched_atoms_between;

    resizable_array_p<Link_Atom> _link_atom;

    resizable_array_p<down_the_bond::DownTheBond> _down_the_bond;

    // Substituents that are not anchored to a ring.
    resizable_array_p<Substituent> _substituent;
    // If any of these match, the whole query match fails.
    resizable_array_p<Substituent> _substituent_no_match;

    // The no matched atoms between matching can be done on either the
    // first path between the ends of the No_Matched_Atoms_Between, or
    // an exhaustive enumeration of all possible paths.
    // Consider the molecule (same moleule, just different smiles)
    // "Oc1cnc(O)cc1"
    // "Oc1ccc(O)nc1"
    // with the smarts
    // "[OH]-c...{>0[n]}c-[OH]". Without exhaustive enumeration of the
    // shortest paths, the first path between the matched atoms in the
    // first smiles will not include the nitrogen atom. Exhaustive
    // enumeration will check all shortest paths found, until a match
    // is found.
    int _no_matched_atoms_between_exhaustive;
    resizable_array_p<No_Matched_Atoms_Between> _nmab;

//  We can also impose distance constraints among matched atoms. Again, we use
//  the bond object just for convenience

//  implement this some time! Need a complex object with a min_max_specifier and
//  two matched atom numbers. Need an array of them

//  September 2000. Bob Coner wanted to restrict queries from forming rings.
//  We can generalise the idea to either exclude implicit rings or requiring one
//  Note that the current implementation only works reliably with queries that are a single chain

    int _implicit_ring_condition;

//  We can also invert the meaning of a query.

    int _rejection;

//  Apr 2004. Restrict matches to largest fragment

    int _only_keep_matches_in_largest_fragment;

    IW_Bits_Base * _fingerprint;

//  Aug 2005. Implement chirality

    resizable_array_p<Substructure_Chiral_Centre> _chirality;

//  Aug 2007. My initial implementation of ring_id's tried to
//  check ring_id's during matching, but that was not correct.
//  Move checking ring_ids to the final match part. For efficiency
//  we store whether or not ring_ids are present in any of the query atoms

    int _ring_ids_present;

    int _fused_system_ids_present;

//  Sept 2011. I want to be able to differentiate "internal" from "external" matches.
//  do that by being able to sort the matches by various criteria

#define SORT_MATCHES_BY_INCREASING_NCON 1
#define SORT_MATCHES_BY_DECREASING_NCON 2
#define SORT_MATCHES_BY_INCREASING_MAXD 4
#define SORT_MATCHES_BY_DECREASING_MAXD 8

    int _sort_matches_by;

//  Mar 2016. Can we gain some efficiencies by keeping track of whether any of the global
//  conditions are set. Also useful for detecting completely unspecified queries

    int _global_conditions_present;

//  Aug 2017. Do any of our Substructure_Atom components have symmetry group info

    int _first_root_atom_with_symmetry_group;

//  One can set a query to report 0 matches unless there are at least
//  a specified number of matches actually present. Seems redundant, not
//  sure why this exists, just use min_hits_needed.

    int _no_match_unless_nhits;

    // If set, all hits will be merged into one. Yields a single
    // embedding with all atoms matched by any of the embeddings.
    int _compress_embeddings;

//  Dec 2019. If per-atom unmatched atom counts are present, store that info
//  in a class variable.

    int _unmatched_atoms_attached_specified;

//  If any of our component Substructure_Atom's have atom type groups,
//  those will need to be checked after matching.
    int _atom_type_groups_present;

    // Jul 2021. Introduce geometric constraints.

    resizable_array_p<geometric_constraints::SetOfGeometricConstraints> _geometric_constraints;

    // Aug 2021. Bond separation between individual matched atoms.
    resizable_array_p<SeparatedAtoms> _separated_atoms;

    // May 2022. MatchedAtomMatch, as an alternative to recursive smarts.
    resizable_array_p<MatchedAtomMatch> _matched_atom_match;
    IW_Logical_Expression _matched_atom_match_operator;

    // Dec 2021. If this query was constructed from a smarts,
    // store that.
    IWString _originating_smarts;

    // June 2022. Implement inter-ring regions
    resizable_array_p<InterRingAtoms> _inter_ring_region;

    // May 2023. Regions
    resizable_array_p<Region> _region;

    // Oct 2023.
    // Groups that are to be found near the matched atoms of the query.
    resizable_array_p<NearbyAtoms> _nearby_atoms;

//  private functions

    void _default_values();

    int _compute_attribute_counts();
    int _locate_chiral_atoms();

    void _collect_all_atoms(extending_resizable_array<Substructure_Atom *> &) const;

    int _add_environment_according_to_matched_atoms(Molecule_to_Query_Specifications & mqs);

//  When creating from smiles or smarts, we need a means of turning off all
//  the global conditions

    int _relax_all_global_conditions();

    int _parse_and_consume_optional_leading_numeric_qualifiers(const_IWSubstring & smarts);

//  int _initialise_sub_array(Molecule & m, Substructure_Query & substitutions_only_at,
//                             Molecule_to_Query_Specifications & mqs);

    int _embedding_is_symmetry_related(const Set_of_Atoms * new_embedding) const;
    int _embedding_is_unique(const Set_of_Atoms *, int *) const;
    int _embedding_is_unique(int, const Set_of_Atoms *) const;

    int _add_embedding(Query_Atoms_Matched & matched_atoms,
                        Set_of_Atoms *,
                        Substructure_Results & results);
    int _find_next_root_atom_embedding(Query_Atoms_Matched & matched_atoms,
                                    Molecule_to_Match & target_molecule,
                                    int * already_matched,
                                    std::unique_ptr<int[]>& matched_by_global_specs,
                                    Substructure_Results & results);
    int _got_embedding(Query_Atoms_Matched & matched_atoms,
                            Molecule_to_Match & target_molecule,
                            int * already_matched,
                            std::unique_ptr<int[]>& matched_by_global_specs,
                            Substructure_Results & results);
    std::unique_ptr<Set_of_Atoms> _make_new_embedding(const Query_Atoms_Matched& matched_atoms) const;

    int _set_start_atom(Target_Atom *);

    int _ring_sizes_specified(resizable_array<int> &) const;

    int _examine_bond_specifications();

    int _atom_numbers_in_geometric_constraints_ok(const geometric_constraints::SetOfGeometricConstraints* constraints) const;

    void _determine_if_ring_ids_are_present();
    void _determine_if_spinach_specifications_are_present();
    void _determine_if_fused_system_ids_are_present();
    void _determine_if_symmetry_groups_present();
    void _determine_if_unmatched_atom_counts_are_present();
    void _determine_if_atom_type_groups_are_present();
    int _one_time_initialisations();

    int _match_elements_needed(Molecule_to_Match & target_molecule) const;

    int RequiredBondsMatch(const Molecule& m);

    int _aromatic_atoms_matches(Molecule_to_Match& target_molecule) const;

//  Function to handle _ring_specification

    int _match_ring_specifications(Molecule_to_Match & target_molecule, std::unique_ptr<int[]>& matched_by_global_specs);

//  Function to handle _ring_system_specification

    int _match_ring_system_specifications(Molecule_to_Match & target_molecule, std::unique_ptr<int[]>& matched_by_global_specs);

    int _discern_global_conditions_present();

//  Function to process the fused_ring_sizes query specifier

    int _match_ring_type_specifications(Molecule_to_Match & target_molecule);
    int _match_nrings_specifications   (Molecule_to_Match & target_molecule);
    int _match_global_specifications   (Molecule_to_Match & target_molecule, std::unique_ptr<int[]>& matched_by_global_specs);
    int _spinach_atoms_match(Molecule_to_Match & target) const;

//  Function to process the heteroatoms specifier

    int _match_heteroatom_specifications(Molecule_to_Match & target_molecule);

//  We pre-compute various properties of the molecule.

    void _allocate_rings_of_sizes(resizable_array<int> &, Molecule *);
    int  _adjust_for_internal_consistency();
    int  _set_ring_membership(Molecule &);

    int _substructure_atoms_known(const Substructure_Atom *, const Substructure_Atom *) const;

//  Find am embedding, with a given atom number matched to query atom 0

    int _find_embedding(Molecule_to_Match & target_molecule, Target_Atom & a,
                                     Query_Atoms_Matched & matched_atoms,
                                     int * already_matched,
                                     std::unique_ptr<int[]>& matched_by_global_specs,
                                     Substructure_Atom * root_atom,
                                     Substructure_Results & results);
    int _find_embedding(Molecule_to_Match & target_molecule,
                        Target_Atom & a, int *,
                        std::unique_ptr<int[]>& matched_by_global_specs,
                        Substructure_Atom *root_atom, Substructure_Results & sresults);
    int _find_embedding(Molecule_to_Match &, Query_Atoms_Matched &,
                        int * already_matched,
                        std::unique_ptr<int[]>& matched_by_global_specs,
                        Substructure_Atom * root_atom, 
                        Substructure_Results & sresults);
                                                 
    int _remove_atoms_with_same_or(Query_Atoms_Matched & matched_atoms,
                                    int atom_to_process, int);
    uint32_t _substructure_search(Molecule_to_Match & target, int * already_matched, std::unique_ptr<int[]>& matched_by_global_specs, Substructure_Results & sresults);
    uint32_t _substructure_search(Molecule_to_Match &, Substructure_Results &);

    int _parse_ring_specifier_object(const msi_object & msi);
    int _parse_ring_specifier(const SubstructureSearch::SingleSubstructureQuery& proto);
    int _parse_ring_system_specifier_object(const msi_object & msi);
    int _parse_ring_system_specifier(const SubstructureSearch::SingleSubstructureQuery& proto);
    int _parse_element_hits_needed_object(const msi_object & msi);
    int _parse_element_hits_needed(const SubstructureSearch::SingleSubstructureQuery& proto);
    int _parse_elements_needed_object(const msi_object & msi);
    int _construct_from_msi_object(const msi_object &, int *);
    int _construct_from_proto(const SubstructureSearch::SingleSubstructureQuery & proto);
    int _construct_environment_from_msi_object(const msi_object &, 
                               extending_resizable_array<Substructure_Atom *> &,
                               resizable_array_p<Substructure_Environment> & env);
    int _construct_environment_from_proto(
    const google::protobuf::RepeatedPtrField<SubstructureSearch::SubstructureEnvironment>& env,
    extending_resizable_array<Substructure_Atom *> & completed,
    resizable_array_p<Substructure_Environment> & destination);
    int _construct_matched_atoms_match(
        const google::protobuf::RepeatedPtrField<SubstructureSearch::MatchedAtomMatch>& matched_atoms_match);
    int _parse_smarts_specifier(const const_IWSubstring & smarts);
    int ParseRange(const_IWSubstring& smarts);
    int _create_from_molecule(MDL_Molecule & m, Molecule_to_Query_Specifications & mqs, const int * include_these_atoms = nullptr);
    int _build_element_hits_needed(const MDL_Molecule & m, const Molecule_to_Query_Specifications & mqs, const int * include_these_atoms);

    int _add_chiral_centre(const Molecule & m, const Chiral_Centre & c);
    int _build_chirality_specification_from_msi_attribute(const IWString & s);
    int _build_chirality_component(const SubstructureSearch::SubstructureChiralCenter& proto,
                                                      void (Substructure_Chiral_Centre::* ptrmfn)(const Substructure_Atom *),
                                                      const SubstructureSearch::AtomNumberHydrogenLonePair& atom_or,
                                                      Substructure_Chiral_Centre& c,
                                                      resizable_array<int>& numbers_encountered);
    int _build_chirality_specification_from_proto(const SubstructureSearch::SubstructureChiralCenter & proto);
    int _build_chirality_specifications_from_atom_chiral_info();
    int _build_chirality_specifications_from_atom_chiral_info(const Substructure_Atom * r);
    int _build_chirality_specification_from_atom_chiral_info(const Substructure_Atom * a);
    int _initialise_chirality_info();
    int _initialise_chirality_info(Substructure_Chiral_Centre * c);
  
    int _build_no_matched_atoms_between(const resizable_array<ThreeDots*>& three_dots);

//  Check global query conditions (attached heteroatom count, etc)

    int  _has_implicit_rings(Query_Atoms_Matched & matched_query_atoms, const int * already_matched) const;
    int  _no_matched_atoms_between_satisfied(Query_Atoms_Matched & matched_atoms) const;
    int  _nmab_satisfied(Molecule_to_Match& target, Query_Atoms_Matched & matched_atoms) const;
    int  _down_the_bond_satisfied(Molecule_to_Match& target_molecule, Query_Atoms_Matched& matched_query_atoms) const;
    int  _substituent_satisfied(Molecule_to_Match& target_molecule, Query_Atoms_Matched& matched_query_atoms) const;
    int _substituent_no_match_matches(Molecule_to_Match& target_molecule, Query_Atoms_Matched& matched_query_atoms) const;
    int  _link_atoms_satisfied(Query_Atoms_Matched & matched_atoms) const;
    int  _link_atom_satisfied(const Link_Atom & l,
                               Query_Atoms_Matched & matched_atoms) const;
//  int  _no_matched_atoms_between_satisfied(Molecule * m, const int * matched, atom_number_t a1, atom_number_t a2) const;
    int  _heteroatoms_matched_satisfied(Query_Atoms_Matched & matched_atoms) const;
    int  _ring_atoms_matched_satisfied(Query_Atoms_Matched & matched_atoms) const;
    int  _distance_between_root_atoms_satisfied(Query_Atoms_Matched & matched_atoms) const;
    int  _fragment_id_conditions_satisfied(Query_Atoms_Matched & matched_atoms) const;
    int  _ring_id_conditions_satisfied(Query_Atoms_Matched & matched_query_atoms) const;
    int  _fused_system_id_conditions_satisfied(Query_Atoms_Matched & matched_query_atoms) const;
    int  _spinach_match_requirements_satisfied(Query_Atoms_Matched & matched_query_atoms, Molecule_to_Match & target_molecule) const;
    int  _attached_heteroatom_count_satisfied(Query_Atoms_Matched & matched_atoms,
                           const int * already_matched) const;
    int  _embedding_ncon_satisfied(Query_Atoms_Matched & matched_atoms,
                           const int * already_matched) const;
    int  _global_query_conditions_also_matched(Query_Atoms_Matched & matched_atoms,
                           const int * already_matched,
                           Molecule_to_Match & target_molecule) const;
    int _chiral_atoms_matched(Query_Atoms_Matched & matched_atoms,
                               Molecule_to_Match & target_molecule) const;
    int AllMatchesFractionAtomsMatched(const Substructure_Results & sresults,
                Molecule_to_Match& target_molecule) const;
    int _symmetry_group_specifications_matches(const Query_Atoms_Matched & query_atoms_matched,
                                        Molecule_to_Match & target) const;
    int _unmatched_atoms_attached_matched(const Query_Atoms_Matched& matched_query_atoms,
                                          const int * already_matched,
                                          Molecule_to_Match& target) const;

    int _atom_type_groupings_matched(const Query_Atoms_Matched& matched_query_atoms) const;

    int InterRingRegionsMatch(Molecule_to_Match& target,
                std::unique_ptr<int[]>& matched_by_global_specs);

    int RegionsMatch(Molecule_to_Match& target, const Query_Atoms_Matched& matched_atoms,
                        const int* already_matched);

//  Various functions for query environments

    int _match_member_of_query_environment(Query_Atoms_Matched & matched_query_atoms,
                        int * previously_matched_atoms,
                        Substructure_Environment * e);
    int _query_environment_xor_group_matched(int * previously_matched_atoms,
                        const std::unique_ptr<int[]>& matched_by_global_specs,
                        const resizable_array<Substructure_Environment *> & xor_group,
                        int * anchor_atom_available);
    int _query_environment_and_group_matched(int * previously_matched_atoms,
                        const std::unique_ptr<int[]>& matched_by_global_specs,
                        const resizable_array<Substructure_Environment *> & or_group,
                        int * anchor_atom_available);
    int _query_environment_or_group_matched(int * previously_matched_atoms,
                        const std::unique_ptr<int[]>& matched_by_global_specs,
                        const resizable_array<Substructure_Environment *> & or_group,
                        int * anchor_atom_available);
    int _environment_rejections_matched(const int atoms_in_target_molecule, int * previously_matched_atoms,
                        const std::unique_ptr<int[]>& matched_by_global_specs,
                        int * env_already_done);
    int _query_environment_also_matched(const int atoms_in_target_molecule, int * previously_matched_atoms,
                        const std::unique_ptr<int[]>& matched_by_global_specs,
                        int * env_already_done);
    int _matched_atoms_match_also_matched(Query_Atoms_Matched & matched_query_atoms, 
                Molecule_to_Match& target);
    int _query_environment_also_matched(Query_Atoms_Matched &, int atoms_in_target_molecule,
                        const std::unique_ptr<int[]>& matched_by_global_specs);

  public:
    Single_Substructure_Query();
    Single_Substructure_Query(const char *);
    Single_Substructure_Query(const IWString &);
    Single_Substructure_Query(const const_IWSubstring &);
    Single_Substructure_Query(iwstring_data_source &);

    ~Single_Substructure_Query();

    int valid() const;
    int ok() const;
    int debug_print(std::ostream &, const IWString & = "");    // not const for this class
    int terse_details(std::ostream &, const IWString &);  // not const for this class
    int print_connectivity_graph(std::ostream &) const;

    int construct_from_msi_object(const msi_object &);
    int ConstructFromProto(const SubstructureSearch::SingleSubstructureQuery& proto);

//  We can create a molecule from a query. This will necessarily be very crude.

    int create_molecule(Molecule &, int = 0, int = 0);   // used for ring perception and fragment identification

//  And we can create a query from a molecule. Again, just approximate

    int create_from_molecule(const Molecule &, Molecule_to_Query_Specifications &, const int * include_these_atoms = nullptr);
    int create_from_molecule(MDL_Molecule &, Molecule_to_Query_Specifications &, const int * include_these_atoms = nullptr);

    void set_comment(const IWString & string) { _comment = string;}
    const IWString & comment() const { return _comment;}

    int numeric_value(double & d, int = 0) const;     // fetch value of I'th numeric value
    void add_numeric_value(double d) { _numeric_value.add(d);}
    void set_numeric_value(double d, int ndx) { _numeric_value[ndx] = d;}
    void discard_all_numeric_values() { _numeric_value.resize(0);}
    void set_respect_initial_atom_numbering(int s) { _respect_initial_atom_numbering = s;}

    void set_embeddings_do_not_overlap(int s) { _embeddings_do_not_overlap = s;}

    void set_only_keep_matches_in_largest_fragment(int s) { _only_keep_matches_in_largest_fragment = s;}
    int  only_keep_matches_in_largest_fragment() const { return _only_keep_matches_in_largest_fragment;}

    int set_min_atoms_to_match(int);
    int set_max_atoms_to_match(int);

    int assign_unique_numbers();
    void assign_unique_id_from_atom_number_if_set(extending_resizable_array<int> & numbers_in_use);

    int unique_numbers_from_initial_atom_numbers();

    int find_one_embedding_per_start_atom() const { return _find_one_embedding_per_start_atom;}

    Substructure_Atom * query_atom_with_initial_atom_number(atom_number_t) const;
    Substructure_Atom * query_atom_with_atom_map_number(int) const;

    const Substructure_Bond * bond_between_atoms(int, int) const;
    const Substructure_Bond * bond_between_atom_map_numbers(int a1, int a2) const;

    int  highest_initial_atom_number() const;    // always computed, even though the variable _highest_initial_atom_number is in the object. Problem of when determined....
    void identify_atom_numbers(extending_resizable_array<int> &) const;
    void identify_atom_map_numbers(extending_resizable_array<int> &) const;
    int  assign_atom_map_numbers(int & amap);
    int  highest_atom_map_number() const;

    int set_find_one_embedding_per_atom(int);
    int set_find_unique_embeddings_only(int);
    int set_min_matches_to_find(int);
    int set_max_matches_to_find(int);
    int set_save_matched_atoms(int);
    void set_ncon(const Min_Max_Specifier<int> & n) {_ncon = n;}
    void set_ncon(int s);
    void set_distance_between_hits(const Min_Max_Specifier<int> & n) { _distance_between_hits = n;}
    void set_no_matched_atoms_between_exhaustive(int s) { _no_matched_atoms_between_exhaustive = s;}

    int add_link_atom(const Link_Atom &);
    int add_no_matched_atoms_between_initial_atom_numbers(int, int);

    int set_normalise_rc_per_hits_needed(int);

    int set_do_not_perceive_symmetry_equivalent_matches(int s);
    int set_perceive_symmetry_equivalent_matches(int s) {
      return set_do_not_perceive_symmetry_equivalent_matches(!s);
    }

    int involves_rings() const;
    int involves_ring_specifications_for_bonds() const;

    int max_atoms_in_query();
    int min_atoms_in_query();

//  Search all atoms in the molecule for a match to the query.

    uint32_t substructure_search(Molecule *, Substructure_Results &);
    uint32_t substructure_search(Molecule_to_Match &, Substructure_Results &);

//  Once we allow decisions into a query, there can be a different number
//  of query atoms which are matched for each match.

    int query_atoms_in_match(int) const;

//  The set of query atoms which made a given match

    const Query_Atoms_Matched * query_atoms_matching(int) const;

    const Substructure_Atom *   query_atom_matching(int, int) const;

    int   root_atoms() const { return _root_atoms.number_elements();}
    const Substructure_Atom * root_atom(int i) const { return _root_atoms[i];}
    int   add_root_atom(Substructure_Atom * r) { return _root_atoms.add(r);}

//  Does a particular atom match the query?

//  int match_atom(Molecule *, atom_number_t);

//  The list of all atoms which match this query

    int matching_atoms(Molecule *, resizable_array<atom_number_t> &);

    int max_query_atoms_matched_in_search() const { return _max_query_atoms_matched.maxval();}

    int print_environment_matches(std::ostream &) const;

    int matches_before_checking_environment() const { return _matches_before_checking_environment;}
    int no_match_to_environment() const { return _no_match_to_environment;}
    int match_to_environemt_rejection() const { return _match_to_environemt_rejection;}
    int environment_present() const { return _environment.number_elements() || _environment_rejections.number_elements();}

//  We may want to ask how many unmatched atoms are bonded to the atoms
//  of the I'th embedding

    int  connections_to_embedding(const Molecule &, int) const;

    int write(const char *);
    int WriteProto(const char *);
    int write_msi(std::ostream &);
    int WriteProto(std::ostream &);
    int write_msi(std::ostream &, int &, const const_IWSubstring &, int = 0);

    int read(iwstring_data_source &);
    int read(const char *);
    int read(const IWString &);
    int ReadProto(const char *);
    int ReadProto(const IWString &);
    int ReadProto(iwstring_data_source& input);
    int next(iwstring_data_source &);

    // Fill in 'proto' with all settings available.
    int BuildProto(SubstructureSearch::SingleSubstructureQuery& proto) const;

    int create_from_smiles(const IWString & smiles);
    int create_from_smarts(const IWString & smarts);

    void set_implicit_ring_condition(int i) { _implicit_ring_condition = i;}

    int query_atom_with_isotope(int) const;   // first query atom that specifies a given isotope

    template <typename T> int any_query_atom(T) const;

    // It can be convenient to forget about any smarts used to build a query.
    // Returns 1 if there was an _originating_smarts to clear.
    int ForgetOriginatingSmarts();

    // THis is a public interface to a private function.
    // If you have a Substructure_Atom and you want to start
    // a substructure search with `root_atom` matched to 'a',
    // this can do that.
    // Implemented this to support the MatchedAtomMatch functionality
    // where we have Substructure_Atom's as the query form.
    int find_embedding(Molecule_to_Match & target_molecule,
                                     Target_Atom & a,
                                     Query_Atoms_Matched & matched_atoms,
                                     int * already_matched,
                                     Substructure_Atom * root_atom,
                                     Substructure_Results & results);
};

/*
  A substructure_query object consists of a set of single_substructure_queries
  linked by logical operators
*/

class Substructure_Query : public resizable_array_p<Single_Substructure_Query>
{
  protected:
    IWString _comment;

  private:

    // Describes the logical operator relationship between the
    // component Single_Substructure_Query.
    IW_Logical_Expression _operator;

//  If set, the Substructure_Results returned by substructure_search will be the
//  aggregate of all atoms matched by all the component queries that match.

    int _each_component_search;

//  private functions

    void _default_values();
    int _parse_smarts_components(const IWString & smarts, char separator);
    int _single_query_construct_from_msi_object(const msi_object & msi);
    int _composite_query_construct_from_msi_object(const msi_object & msi);
    int _add_component_from_smarts(const const_IWSubstring & smarts, int op);
    int _create_query_and_add(MDL_Molecule & m, Molecule_to_Query_Specifications & mqs, const int * include_these_atoms);
    int _enumerate_list_atom_possibilities(const MDL_Molecule & m,
                                                        Molecule_to_Query_Specifications & mqs,
                                                        const resizable_array_p<Link_Atom> & link_atoms,
                                                        Link_Atom_Current_State * lacc,
                                                        int nlink,
                                                        int ndx,
                                                        const int * include_these_atoms);

    int _substructure_search(Molecule_to_Match &, Substructure_Results &);

  public:
    Substructure_Query();
    Substructure_Query(const char *);
    Substructure_Query(const const_IWSubstring &);
    Substructure_Query(const IWString &);
    ~Substructure_Query();

    int ok() const;
    int debug_print(std::ostream &) const;
    int terse_details(std::ostream &) const;
    int active() const { return _number_elements;}
    int print_connectivity_graph(std::ostream &) const;

    const IWString & comment() const { return _comment;}
    void set_comment(const IWString & cm) { _comment = cm;}

    int highest_initial_atom_number() const;
    int highest_atom_map_number() const;
    void identify_atom_numbers(extending_resizable_array<int> &) const;
    void identify_atom_map_numbers(extending_resizable_array<int> &) const;
    void assign_unique_id_from_atom_number_if_set(extending_resizable_array<int> & numbers_in_use);
    int  assign_atom_map_numbers(int & amap);

    const Substructure_Bond * bond_between_atoms(int, int) const;
    const Substructure_Bond * bond_between_atom_map_numbers(int a1, int a2) const;

    void set_each_component_search(int s) { _each_component_search = s;}

    void set_respect_initial_atom_numbering(int s);

    void set_only_keep_matches_in_largest_fragment(int s);
    int  only_keep_matches_in_largest_fragment() const;

    void set_embeddings_do_not_overlap(int s);

    int add(Single_Substructure_Query *, int = IW_LOGEXP_OR);    // add a component and include the operator

    int read(iwstring_data_source &);
    int read(const char *);
    int read(const IWString &);
    int read(const const_IWSubstring &);

    int ReadProto(iwstring_data_source& input);
    int ReadProto(const IWString& fname);

    int construct_from_msi_object(const msi_object &);
    int ConstructFromProto(const SubstructureSearch::SubstructureQuery& proto);

    int create_from_smiles(const IWString & smiles);
    int create_from_smarts(const IWString & smarts);
    int CreateFromSmarts(const std::string & smarts);

    int create_from_molecule(const Molecule &, Molecule_to_Query_Specifications &, const int * include_these_atoms = nullptr);
    int create_from_molecule(MDL_Molecule &, Molecule_to_Query_Specifications &, const int * include_these_atoms = nullptr);
//  int create_from_molecule(Molecule_to_Query_Specifications &,
//                            Molecule &,
//                            const MDL_File_Data &);
    int unique_numbers_from_initial_atom_numbers();

    Substructure_Atom * query_atom_with_initial_atom_number(int) const;
    Substructure_Atom * query_atom_with_atom_map_number(int) const;

    int write_msi(IWString& fname);
    int write_msi(std::ostream& fname);
    int write_msi(std::ostream& fname, int &, int = 0);

    // Write proto representation to 'os'.
    int WriteProto(std::ostream& os) const;
    int WriteProto(IWString& fname) const;
    int WriteProto(IWString_and_File_Descriptor& output) const;

    // Return a proto containing the contents of the Substructure_Query;
    SubstructureSearch::SubstructureQuery BuildProto() const;

    void set_ncon(const Min_Max_Specifier<int> &);
    void set_distance_between_hits(const Min_Max_Specifier<int> &);

    // Substructure searching with a composite query is complex.
    // If there is only one component, the return value will be what that
    // single component found, and Substructure_Results will contain
    // the matched atoms for that query.
    // If there are multiple components, there can only be a yes/no return
    // value. And in that case, Substructure_Results will be undefined.
    // On a negative match, it might contain some matched atoms (the 
    // expression was an XOR).
    // On a positive match, it may (or may not) contain some matched atoms,
    // because the last component matched might have been a non-match in
    // an XOR.
    uint32_t substructure_search(Molecule *);
    uint32_t substructure_search(Molecule *, Substructure_Results &);
    uint32_t substructure_search(Molecule &);
    uint32_t substructure_search(Molecule &, Substructure_Results &);
    uint32_t substructure_search(Molecule_to_Match &);
    uint32_t substructure_search(Molecule_to_Match &, Substructure_Results &);

//  Sometimes we want to search each component of a query regardless of what operators
//  might be present.

    uint32_t substructure_search_do_each_component(Molecule_to_Match & target, Substructure_Results & sresults);

    int set_find_one_embedding_per_atom(int);
    int set_find_unique_embeddings_only(int);
    int set_min_matches_to_find(int);
    int set_max_matches_to_find(int);
    int set_do_not_perceive_symmetry_equivalent_matches(int);
    int set_perceive_symmetry_equivalent_matches(int s) {
      return set_do_not_perceive_symmetry_equivalent_matches(!s);
    }
    int set_save_matched_atoms(int);
    int set_min_atoms_to_match(int);   // set the _natoms attribute
    int set_max_atoms_to_match(int);   // set the _natoms attribute

    int max_query_atoms_matched_in_search() const;
    int max_atoms_in_query();

    int print_environment_matches(std::ostream &) const;

    int numeric_value(double & d, int = 0) const;
    void set_numeric_value(double d, int ndx);
    void discard_all_numeric_values();
    void add_numeric_value(double d);

    template <typename T> int any_query_atom(T) const;
    template <typename T> void each_query_atom(T) const;

    // Applied to each component Single_Substructure_Query
    int ForgetOriginatingSmarts();
};


//extern Substructure_Query * create_query(Molecule *, int = 0, int = 0, int = 0);
//extern int create_query(Molecule *, Substructure_Query &, int = 0, int = 0);

//extern Single_Substructure_Query * new_substructure_query(const char *, int);
//extern Single_Substructure_Query * next_substructure_query(iwstring_data_source &, int);

extern int atom_matches_queries(Molecule *, atom_number_t, resizable_array_p<Single_Substructure_Query> &);

extern void set_respect_aliphatic_smarts(int);

extern void set_report_multiple_hits_threshold(unsigned int s);

extern int use_fingerprints_for_screening_substructure_searches();
extern void set_use_fingerprints_for_screening_substructure_searches(int);

/*
  0 means initial behaviour - if there are multiple largest fragments with the same atom
  atom count, the first one will be taken as the largest. 
  1 means if there are multiple largest fragments, then matches in all of them will be
  preserved.
*/

extern void set_remove_hits_not_in_largest_fragment_behaviour(int s);

/*
  by default, H in a smarts means exactly one Hydrogen.
  We can change that behaviour so it means a minimum of one
*/

extern void set_h_means_exactly_one_hydrogen(int s);

extern void set_atom_environment_only_matches_unmatched_atoms(int);

extern void set_ignore_chirality_in_smarts_input(int);
extern int  ignore_chirality_in_smarts_input();
extern void set_query_environment_must_match_unmatched_atoms(const int s);

// This is an internal thing used during matches. Not for external use

extern int remove_atoms_with_same_or(Query_Atoms_Matched & matched_atoms,
                           int istart, int orid);

template <typename T>
int
process_queries(Command_Line & cl, resizable_array_p<T> & queries,
                 int verbose = 0, char option = 'q');

template <typename T>
int
process_files_of_queries(Command_Line & cl, resizable_array_p<T> & queries,
                 int inherit_directory_path, int verbose = 0, char option = 'Q');

template <typename T>
int
queries_from_file(const const_IWSubstring & fname, resizable_array_p<T> & queries,
                   int inherit_directory_path, int verbose = 0);
template <typename T>
int
queries_from_file(const IWString & fname, resizable_array_p<T> & queries,
                   int inherit_directory_path, int verbose = 0);
template <typename T>
int
queries_from_ISIS_query_file(const const_IWSubstring & fname,
                              resizable_array_p<T> & queries,
                              int verbose);
template <typename T>
int
query_from_ISIS_query_file(MDL_Molecule & m,
                            Molecule_to_Query_Specifications & mqs,
                            resizable_array_p<T> & queries, int verbose);
template <typename T>
int
queries_from_ISIS_query_file(const const_IWSubstring & fname,
                              resizable_array_p<T> & queries, int verbose);
template <typename T>
int
queries_from_file_of_molecules(const const_IWSubstring & fname,
                                resizable_array_p<T> & queries,
                                int verbose);
template <typename T>
int
process_cmdline_token(const char, const const_IWSubstring & token,
                       resizable_array_p<T> & queries, int verbose);
/*
  When parsing the form '-q M:fname' we can add directives in the form
  '-q M:fname%onlysub=smarts' 
  Here is the definition of the separator
*/

#define DIRECTIVE_SEPARATOR_TOKEN '%'

template <typename M, typename Q>
int
first_query_to_match(M & m, resizable_array_p<Q> & q)
{
  for (auto i = 0; i < q.number_elements(); ++i)
  {
    if (q[i]->substructure_search(m))
      return i;
  }

  return -1;
}

template <typename T>
int
Substructure_Query::any_query_atom(T todo) const
{
  for (int i = 0; i < _number_elements; ++i)
  {
    if (_things[i]->any_query_atom(todo))
      return 1;
  }

  return 0;
}
template <typename T> int 
Single_Substructure_Query::any_query_atom(T todo) const
{
  for (int i = 0; i < _root_atoms.number_elements(); ++i)
  {
    if (_root_atoms[i]->any_query_atom(todo))
      return 1;
  }

  return 0;
}
template <typename T> int 
Substructure_Atom::any_query_atom(T todo) const
{
  for (int i = 0; i < _children.number_elements(); ++i)
  {
    if (todo(_children[i]))
      return 1;
  }

  return 0;
}

namespace substructure_spec {

// Exposed only for testing
int SmartsNumericQualifier(const char * input,
                       int max_chars,
                       Min_Max_Specifier<int>& result);
int SmartsFetchNumeric(const char * string, int nchars, int & value, 
                     int & qualifier);

int
SmartsNumericQualifier(const char * input,
                 int max_chars,
                 Min_Max_Specifier<int>& result);

}  // namespace substructure_spec

namespace ss_ring {
void
FillMatchedAtomsArray(const Ring& r, 
                      int value,
                      int matoms,
                      std::unique_ptr<int[]>& matched_by_global_specs);
void
SetTargetGlobalIds(const Set_of_Atoms& r, 
                   int value,
                   Molecule_to_Match& target);
}  // namespace ss_ring

namespace iwsubstructure {
// Return the contents of the next
// query {
// ...
// }
// grouping found in `input`.
std::optional<std::string> GetNextQueryTextProto(iwstring_data_source& input);
}

#endif  // MOLECULE_LIB_SUBSTRUCTURE_H_
