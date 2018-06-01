#ifndef IWREACTION_H
#define IWREACTION_H

#include <random>

#include "cmdline.h"
#include "msi_object.h"

#include "substructure.h"
#include "istream_and_type.h"

#include "reaction_match_conditions.h"
#include "toggle_kekule_form.h"

/*
  Names used in msi_object .rxn files
*/

#define NAME_OF_REACTION_OBJECT "reaction"
#define NAME_OF_SCAFFOLD_OBJECT "scaffold"
#define NAME_OF_SIDECHAIN_OBJECT "reagent"
#define OLD_NAME_OF_SIDECHAIN_OBJECT "sidechain"
#define NAME_OF_REACTION_COMMENT "comment"
#define NAME_OF_NO_REACTION_OBJECT "no_reaction"
#define NAME_OF_MATCH_VIA_ATOM_MAP "match_via_atom_map"
#define NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE "max_matches"
#define NAME_OF_BOND_TO_BREAK_ATTRIBUTE "break_bond"
#define NAME_OF_SINGLE_BOND_ATTRIBUTE "single_bond"
#define NAME_OF_DOUBLE_BOND_ATTRIBUTE "double_bond"
#define NAME_OF_TRIPLE_BOND_ATTRIBUTE "triple_bond"
#define NAME_OF_WEDGE_BOND_ATTRIBUTE "wedge_bond"
#define NAME_OF_AROMATIC_BOND_ATTRIBUTE "aromatic_bond"
#define NAME_OF_REMOVE_ATOM_ATTRIBUTE "remove_atom"
#define NAME_OF_REMOVE_FRAGMENT_ATTRIBUTE "remove_fragment"
#define NAME_OF_KEEP_FRAGMENT_ATTRIBUTE "keep_fragment"
#define NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE "join"
#define NAME_OF_INTER_PARTICLE_SUBSTITUTE_ATOM_ATTRIBUTE "inter_particle_substitute_atom"
#define NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE "substitute"
#define NAME_OF_APPLY_TO_GROWING_MOLECULES_ATTRIBUTE "apply_to_growing_molecule"
#define NAME_OF_SIDECHAIN_MOLECULE "molecule"
#define NAME_OF_SIDECHAIN_SMILES "smiles"
#define NAME_OF_SIDECHAIN_FILE "file"
#define NAME_OF_CHANGE_ELEMENT_ATRRIBUTE "change_element"
#define NAME_OF_FORMAL_CHARGE_ATTRIBUTE "formal_charge"
#define NAME_OF_CHANGE_FORMAL_CHARGE_ATTRIBUTE "change_formal_charge"
#define NAME_OF_ISOTOPE_ATTRIBUTE "isotope"
#define NAME_OF_INCREMENT_ISOTOPE_ATTRIBUTE "increment_isotope"
#define NAME_OF_INVERT_ISOTOPE_ATTRIBUTE "invert_isotope"
#define NAME_OF_DIHEDRAL_ANGLE_ATTRIBUTE "dihedral_angle"
#define NAME_OF_BOND_LENGTH_ATTRIBUTE "bond_length"
#define NAME_OF_BOND_ANGLE_ATTRIBUTE "bond_angle"
//#define NAME_OF_ROTATE_FRAGMENT_ATTRIBUTE "rotate_fragment"
#define NAME_OF_3D_REPLACE_ATTRIBUTE "3d_replace"
#define NAME_OF_QUERY_SMARTS_ATTRIBUTE "smarts"
#define NAME_OF_INACTIVE_ATTRIBUTE "inactive"
#define NAME_OF_STEREO_CENTRE_ATTRIBUTE "stereo"
#define NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE "invert_chirality"
#define NAME_OF_OPTIONAL_STEREO_CENTRE_ATTRIBUTE "optional_stereo"
#define NAME_OF_RECOMPUTE_IMPLICIT_HYDROGENS_ATTRIBUTE "recompute_implicit_hydrogens"
#define NAME_OF_FIND_KEKULE_FORM_FOR_INVALID_VALENCES_ATTRIBUTE "find_kekule"
#define NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE "one_embedding_per_start_atom"
#define NAME_OF_UNIQUE_EMBEDDINGS_ONLY_ATTRIBUTE "unique_embeddings_only"
#define NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE "embeddings_do_not_overlap"
#define NAME_OF_IGNORE_SYMMETRY_RELATED_MATCHES_ATTRIBUTE "ignore_symmetry_related_embeddings"
#define NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE "match_each_component"
#define NAME_OF_TOGGLE_KEKULE_FORM_ATTRIBUTE "toggle_kekule_form"
#define NAME_OF_REMOVE_CHIRAL_CENTRE_ATTRIBUTE "remove_chiral_centre"
#define NAME_OF_MAKE_IMPLICIT_HYDROGEN_EXPLICIT_ATTRIBUTE "hexplicit"
#define NAME_OF_UNLINK_UNMATCHED_ATOMS "unlink_unmatched_atoms"
#define NAME_OF_APPEND_TO_NAME "append_to_name"
#define NAME_OF_NOOP_REACTION "noop_reaction"

/*
  We have the capability to read the query from a separate file rather
  than embedding it in the reaction file
*/

#define NAME_OF_QUERY_FILE_ATTRIBUTE "query_file"

/*
  When specifying inter-particle bonds and substitution bonds we
  cannot use Bond objects, because the two atom numbers may be the
  same.

  The Substitute_Atom operation is designed for cases where the
  stereochemistry about an atom must be preserved. Specifically,
  imagine the case of a 4 connected stereo centre in the scaffold.
  One of its attachments is removed and an atom from the sidechain
  joined. Using an Inter_Particle_Bond operation for this doesn't
  preserve the stereochemistry, because of how the reaction is
  performed. First the existing bond is broken, and then the new
  atom is attached, but when the existing bond is broken, the
  stereo centre is destroyed.

  The difference between an Inter_Particle_Bond and a Substitute_Atom
  directive is that the Substitute_Atom operation will preserve a
  stereo centre
*/

class Pair_of_Atoms
{
  protected:
    int _a1;
    int _a2;

  public:
    Pair_of_Atoms (atom_number_t q1, atom_number_t q2) { _a1 = q1; _a2 = q2;}

    int a1    () const { return _a1;}
    int a2    () const { return _a2;}

    void set_a1 (atom_number_t a) { _a1 = a;}     // no checking - you could set _a1 == _a2!!
    void set_a2 (atom_number_t a) { _a2 = a;}     // no checking
};

/*
  For dealing with multi-component reactions, we have matched atoms in any component.
  By convention, component -1 is the scaffold.
  Often _component will remain set to -2, which means in the current component
*/

#define MAIC_CURRENT_COMPONENT -2

class Matched_Atom_in_Component
{
  protected:
    int _matched_atom;
    int _component;

  public:
    Matched_Atom_in_Component ();

    int debug_print(std::ostream & os) const;

    int construct (const const_IWSubstring &, int = -1);

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref);

    int in_scaffold  () const { return -1 == _component;}
    void set_in_scaffold () { _component = -1;}
    int in_component () const { return _component;}
    int matched_atom () const { return _matched_atom;}
    int atom () const { return _matched_atom;}
    int in_current_component () const { return MAIC_CURRENT_COMPONENT == _component;}

    void set_matched_atom (int m) { _matched_atom = m;}
    void set_in_component (int c) { _component = c;}
};

extern std::ostream &
operator << (std::ostream &, const Matched_Atom_in_Component &);


class Inter_Particle_Bond
{
  private:
    Matched_Atom_in_Component _a1;
    Matched_Atom_in_Component _a2;
    bond_type_t   _bt;

  public:
    Inter_Particle_Bond ();
    Inter_Particle_Bond (int f1, int q1, int q2, bond_type_t bt);

    int ok () const { return OK_BOND_TYPE (_bt);}
    int debug_print (std::ostream &) const;

    const Matched_Atom_in_Component & a1 () const { return _a1;}
    const Matched_Atom_in_Component & a2 () const { return _a2;}

    bond_type_t   btype () const { return _bt;}

    int build (const IWString & s);
    int write_msi (std::ostream & os, const IWString & ind) const;

    int construct_from_msi_attribute (const msi_attribute &);

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> &);

    int is_part_of_component (int);
};

std::ostream & operator << (std::ostream &, const Inter_Particle_Bond &);

class Reaction_Wedge_Bond : public Pair_of_Atoms
{
  private:
    int _direction;

  public:
    Reaction_Wedge_Bond ();

    int construct_from_msi_attribute (const msi_attribute *);

    int write_msi (std::ostream &, const IWString &, const const_IWSubstring &) const;

    int direction () const { return _direction;}

    int process (Molecule &, const Set_of_Atoms &, int) const;
};

/*
  A stereo preserving substitution can occur on an implicit hydrogen
*/

/*class Substitute_Atom : public Pair_of_Atoms
{
  private:
    int _anchor;

  public:
    Substitute_Atom (int c, int q1, int q2) : Pair_of_Atoms (q1, q2) , _anchor (c) {}

    int write_msi (std::ostream &, const IWString &, const const_IWSubstring &) const;

    int anchor () const { return _anchor; }
};*/

/*
  In order to specify element transformations, we need an object to hold
  those specifications
*/

class Reaction_Change_Element
{
  private:
    int _atom;
    const Element * _element;

  public:
    Reaction_Change_Element ();

    int process (Molecule &, const Set_of_Atoms &, int);

    void set_atom (int a) { _atom = a;}
    void set_element (const Element * e) { _element = e;}

    int atom () const { return _atom;}

    int construct_from_msi_attribute (const msi_attribute *);

    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;
};

/*
  We may want to assign or change some formal charges
*/

class Reaction_Formal_Charge
{
  private:
    int _atom;
    formal_charge_t _fc;

  public:
    Reaction_Formal_Charge ();

    void set_atom (int a) { _atom = a;}
    void set_formal_charge (formal_charge_t f) { _fc = f;}

    int atom () const { return _atom;}

    int process (Molecule &, const Set_of_Atoms &, int);

    int construct_from_msi_attribute (const msi_attribute *);

    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;
};

/*
  We may want to increment or decrement a formal charge
*/

class Reaction_Change_Formal_Charge
{
  private:
    int _atom;
    int _delta;

  public:
    Reaction_Change_Formal_Charge ();

    int atom () const { return _atom;}

    int process (Molecule &, const Set_of_Atoms &, int);

    int construct_from_msi_attribute (const msi_attribute *);

    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;
};

class Reaction_Place_Isotope
{
  protected:
    int _atom;
    int _isotope;

  public:
    Reaction_Place_Isotope ();
    Reaction_Place_Isotope (int a, int iso) : _atom (a), _isotope (iso) {}

    int atom () const { return _atom;}
    int isotope () const { return _isotope;}

    void set_atom (int s) { _atom = s;}
    void set_isotope (int s) { _isotope = s;}

    int process (Molecule &, const Set_of_Atoms &, int);

    int construct_from_msi_attribute (const msi_attribute *);

    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;
};

/*
  Actually the Reaction_Place_Isotope and Reaction_Increment_Isotope should inherit from
  a common base class, but I had already written Reaction_Place_Isotope.
*/

class Reaction_Increment_Isotope : public Reaction_Place_Isotope
{
  private:
  public:
    Reaction_Increment_Isotope () : Reaction_Place_Isotope () {}
    Reaction_Increment_Isotope (int a, int iso) : Reaction_Place_Isotope (a, iso) {}

    int process (Molecule &, const Set_of_Atoms &, int);
};

/*
  Needed to be able to change all labelled atoms to zero isotope
  and set all zero isotopes to some number
*/

class Reaction_Invert_Isotope : public Reaction_Place_Isotope
{
  private:
  public:
    Reaction_Invert_Isotope () : Reaction_Place_Isotope () {}
    Reaction_Invert_Isotope (int a, int iso) : Reaction_Place_Isotope (a, iso) {}

    int process (Molecule &, const Set_of_Atoms &, int);
};

/*
  The sidechain objects make extensive use of a molecule and 
  a set of atoms within that molecule

  At first I tried to declare this as 'public Molecule, public Set_of_Atoms'
  but it seemed a little too complicated. For example, I could not figure
  out the syntax to access the Set_of_Atoms part alone?

*/

class Sidechain_Reaction_Site;

class Molecule_and_Embedding : public Molecule
{
  private:
    Set_of_Atoms _embedding;

//  Each Molecule_and_Embedding reacts under the guidance of a sidechain.

    Sidechain_Reaction_Site * _owning_sidechain;

  public:
    Molecule_and_Embedding ();

    const Set_of_Atoms * embedding () const { return &_embedding;}

//  A substructure query has just matched. We need to copy the
//  appropriate embedding to ourself.

    int collect_matched_atoms (const Substructure_Results &, int = 0);   // 0 means the 0'th embedding

    Sidechain_Reaction_Site * owning_sidechain () const { return _owning_sidechain;}
    void set_owning_sidechain (Sidechain_Reaction_Site * s) { _owning_sidechain = s;}

    int set_embedding (const Set_of_Atoms &);

    int do_toggle_kekule_form (Toggle_Kekule_Form &);
};

/*
  Mar 98. Ran into cases where under certain circumstances the
  reaction does not occur. For example, in the scaffold we have
  OH-CH2 or Cl-CH2, and in the sidechain we have OH-a or Br-a
  The OH in the sidechain can displace both the OH and Cl in the
  scaffold, but the Br cannot displace the Cl.
  Seems the most general way to implement this is as substructure
  queries, but, we don't have a way of ensuring that the atoms
  matched are the same as those matched by the underlying query.
  Perhaps this is OK!

  Either or both of these queries can be active
*/

class No_Reaction
{
  private:
    IWString _comment;

    Substructure_Query _scaffold_no_reaction;
    Substructure_Query _sidechain_no_reaction;

  public:
    int ok () const;
    int debug_print (std::ostream &) const;

    int construct_from_msi_object (const msi_object &);
    int write_msi (std::ostream &, int &, int);
};

/*
  Replace_Atom is the same as Substitute_Atom, but without the requirement of
  an anchor atom. Also, the atoms can be in either the scaffold or the sidechain
*/

class Replace_Atom
{
  private:
    Matched_Atom_in_Component _a1;
    Matched_Atom_in_Component _a2;

  public:
    int construct_from_msi_object (const msi_attribute *, int);

    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;

    const Matched_Atom_in_Component & a1 () const { return _a1;}
    const Matched_Atom_in_Component & a2 () const { return _a2;}

    Matched_Atom_in_Component & a1 () { return _a1;}
    Matched_Atom_in_Component & a2 () { return _a2;}

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> &);

    int process (Molecule &, const Set_of_Atoms *, const Set_of_Atoms *, int) const;
};

extern std::ostream & operator << (std::ostream &, const Replace_Atom &);

/*
  There are a bunch of temporary items that are needed during an enumeration.
  What I found was that argument lists were getting long, so we put these
  temporary things into an object
*/

class Enumeration_Temporaries
{
  private:
    Set_of_Atoms _atoms_to_be_removed;

    int * _atoms_in_growing_molecule;

    resizable_array<const Atom *> _fragments_to_be_removed;

    const Molecule_and_Embedding ** _reagent;

    resizable_array_p<Chiral_Centre> _saved_chiral_centre;

  public:
    Enumeration_Temporaries (int);
    ~Enumeration_Temporaries ();

//  int initialise (int);

    Set_of_Atoms & atoms_to_be_removed () { return _atoms_to_be_removed;}
    const Set_of_Atoms & atoms_to_be_removed () const { return _atoms_to_be_removed;}

    int * atoms_in_growing_molecule () { return _atoms_in_growing_molecule;}
    const int * atoms_in_growing_molecule () const { return _atoms_in_growing_molecule;}

    resizable_array<const Atom *> & fragments_to_be_removed () { return _fragments_to_be_removed;}
    const resizable_array<const Atom *> & fragments_to_be_removed () const { return _fragments_to_be_removed;}

    const Molecule_and_Embedding ** reagent () const { return _reagent;}
};


/*
  In order to do things like the Diels-Alder reaction, we need to be
  able to specify stereo centres to be constructed

  Sometimes we don't know if there will be a matched atom or an
  implicit hydrogen present. The value SCC_WHATEVER is assigned in
  such cases
*/

#define SCC_WHATEVER -12

class Stereo_Centre_Component : public Matched_Atom_in_Component
{
  private:
    int _implicit_hydrogen;

  public:
    Stereo_Centre_Component ();

    int ok () const;
    int active () const;

    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    void set_matched_atom (int m);

    void set_whatever () {_matched_atom = SCC_WHATEVER;}
    int  is_whatever () const { return _matched_atom == SCC_WHATEVER;}

    int  implicit_hydrogen () const { return _implicit_hydrogen;}
    void set_implicit_hydrogen (int);
};

extern std::ostream & operator << (std::ostream &, const Stereo_Centre_Component &);

class Reaction_Stereo_Centre
{
  private:
    Stereo_Centre_Component _ssc[5];

//  We can have optional stereo centres - form if you can, but if
//  there aren't enough atoms, don't complain (too much)

    int _optional;

//  private functions

    int _find_missing_connection (Molecule & m, atom_number_t zatom, Chiral_Centre * c) const;

  public:
    Reaction_Stereo_Centre ();

    int debug_print (std::ostream &) const;

//  When we have a stereo centre that is just in the scaffold, we don't
//  want to have to specify all atoms as "S.1", "S.2", ...

    void all_atoms_in_scaffold ();

    int construct_from_msi_attribute (const msi_attribute *);
    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;

    void set_optional (int s) { _optional = s;}
    int  optional () const { return _optional;}

    int process (Molecule & result,
                                 const Set_of_Atoms * scaffold_embedding,
                                 const Enumeration_Temporaries & etmp) const;

//  int make_stereo_centre (Molecule & result,
//                           const Set_of_Atoms * scaffold_embedding,
//                           const Set_of_Atoms * sidechain_embedding,
//                           const Set_of_Atoms & atoms_to_be_removed,
//                           int atoms_already_present);

//  When creating one of these out of an ISIS reaction file, it is handy to be able to grab individual components

    Stereo_Centre_Component & centre     () { return _ssc[0];}
    Stereo_Centre_Component & top_front  () { return _ssc[1];}
    Stereo_Centre_Component & top_back   () { return _ssc[2];}
    Stereo_Centre_Component & left_down  () { return _ssc[3];}
    Stereo_Centre_Component & right_down () { return _ssc[4];}

    const Stereo_Centre_Component & centre     () const { return _ssc[0];}

    int atom () const { return _ssc[0].matched_atom ();}
};

/*
  Feb 2001, Minmin needed to do geometrically correct reactions
*/

class Reaction_Dihedral_Angle
{
  private:
    Matched_Atom_in_Component _atom[4];

    angle_t _desired_angle;

//  private functions

    int _twist_around_bond (Molecule & m, atom_number_t a2, atom_number_t a3, int * to_move) const;

  public:
    Reaction_Dihedral_Angle ();

    void all_atoms_in_scaffold ();

    int construct_from_msi_attribute (const msi_attribute *);
    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> &);

    int process (Molecule &, const Set_of_Atoms *, const Enumeration_Temporaries &) const;
};

class Reaction_Bond_Angle
{
  private:
    Matched_Atom_in_Component _atom[3];

    angle_t _desired_angle;

  public:
    Reaction_Bond_Angle ();

    void all_atoms_in_scaffold ();

    int construct_from_msi_attribute (const msi_attribute *);
    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> &);

    int process (Molecule &, const Set_of_Atoms *, const Enumeration_Temporaries &) const;
};

extern std::ostream &
operator << (std::ostream &, const Reaction_Dihedral_Angle &);

class Reaction_Bond_Length {
  private:
    Matched_Atom_in_Component _atom[2];

    distance_t _desired_length;

  public:
    Reaction_Bond_Length ();

    const Matched_Atom_in_Component & a1 () const { return _atom[0];}
    const Matched_Atom_in_Component & a2 () const { return _atom[1];}

    void all_atoms_in_scaffold ();

    int construct_from_msi_attribute (const msi_attribute *);
    int write_msi (std::ostream &, const const_IWSubstring &, const const_IWSubstring &) const;

    int process (Molecule &, const Set_of_Atoms *, const int *, const Molecule_and_Embedding **) const;

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> &);

    distance_t desired_length () const { return _desired_length;}
};

class Reaction_Rotate_Fragment
{
  private:
    Matched_Atom_in_Component _atom[3];

    angle_t _desired_angle;

  public:
    Reaction_Rotate_Fragment();

    void all_atoms_in_scaffold ();

    int construct_from_msi_attribute (const msi_attribute * att);

    int write_msi (std::ostream & os, const const_IWSubstring & ind,
                   const const_IWSubstring & attribute_name) const;

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref);

    int process (Molecule & m, const Set_of_Atoms * scaffold_embedding,
                 const Enumeration_Temporaries & etmp) const;
};

class Reaction_3D_Replace
{
  private:
    int _n;

    Matched_Atom_in_Component * _a1;
    Matched_Atom_in_Component * _a2;

    double * _weight;

  public:
    Reaction_3D_Replace();
    ~Reaction_3D_Replace();

    int construct_from_msi_attribute (const msi_attribute * att);

    int write_msi (std::ostream & os, const const_IWSubstring & ind,
                   const const_IWSubstring & attribute_name) const;

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref);

    int process (Molecule & m, const Set_of_Atoms * scaffold_embedding,
                 const Enumeration_Temporaries & etmp) const;
};


/*
  Jul 99. Finally need to do something about "libraries" which consist
  of different types of reactions.
  Each of these "libraries"
*/

/*
  The sidechain and the scaffold are described by a substructure query
  We can also specify bonds to be made and removed within the scaffold
  or sidechain. Also, as a result of breaking bonds, or removing atoms,
  we may end up with a fragment to be removed

  We have a generic object which applies to both, and then scaffold
  and sidechain types are further specialised in their own classes

  Ran into several cases where multiple query matches cause problems - the most
  glaring being one that looks for O=C-C-C*C where the last bond is double or
  aromatic. If there is an aromatic ring, this will always match twice, which
  is silly. So, the _matched_atom_changed array is used to keep track of
  which members of the query matches are changed. This will break if we get
  embeddings with different numbers of atoms in them
*/

class Reaction_Site : public Substructure_Query
{

  protected:

//  When reading reaction files, each component has a unique identifier

    int _unique_id;

    IWString _comment;

    resizable_array_p<Bond> _bonds_to_be_made;
    resizable_array_p<Bond> _bonds_to_be_broken;
    resizable_array<int>    _atoms_to_be_removed;
    resizable_array<int>    _fragments_to_be_removed;
    resizable_array<int>    _fragments_to_be_kept;

    resizable_array_p<Reaction_Change_Element>  _elements_to_change;
    resizable_array_p<Reaction_Formal_Charge>   _formal_charges_to_assign;
    resizable_array_p<Reaction_Change_Formal_Charge> _formal_charges_to_change;
    resizable_array_p<Reaction_Place_Isotope>   _isotopes_to_assign;
    resizable_array_p<Reaction_Increment_Isotope> _isotopes_to_increment;
    resizable_array_p<Reaction_Invert_Isotope>  _isotopes_to_invert;

    resizable_array_p<Reaction_Dihedral_Angle> _reaction_dihedral_angle;
    resizable_array_p<Reaction_Bond_Length> _reaction_bond_length;
    resizable_array_p<Reaction_Bond_Angle> _reaction_bond_angle;
//  resizable_array_p<Reaction_Rotate_Fragment> _reaction_rotate_fragment;
    resizable_array_p<Reaction_3D_Replace> _reaction_3d_replace;

//  Do we need to worry about coordinates when doing stereo preserving substitutions

    int _3d;

//  Jan 2000. Michal Vieth needed something to get rid of explicit hydrogens

    resizable_array<int> _unfix_implicit_hydrogens;

//  Jan 2001 wedge bonds

    resizable_array_p<Reaction_Wedge_Bond> _wedge_bonds_to_place;

//  May 2002 new way of doing atom replacements

    resizable_array_p<Replace_Atom> _replace_atom;

//  May 2002 Inactive features

    resizable_array_p<Substructure_Query> _inactive;

//  May 2002 invert a chiral centre.

    resizable_array<int> _stereo_centres_to_invert;

//  Jan 2003. Optionally allow for automatic toggling of Kekule forms

    Toggle_Kekule_Form _toggle_kekule_form;

//  Jan 2003. We can remove an existing chiral centre

    resizable_array<int> _chiral_centres_to_remove;

//  which matched atoms are changed by

    int * _matched_atom_changed;

//  Jan 2003. Various ways of dealing with the problem of multiple hits

    int _ignore_multiple_matches_involving_atoms_not_changing;
    int _ignore_multiple_matches_involving_changing_atoms;

// Aug 2003. Chiral centres are difficult. During do_makes_breaks we identify any
// chiral centres that will be destroyed. We extract them from the molecule and
// save them in the Enumeration_Temporaries object

    resizable_array_p<Chiral_Centre> _saved_chiral_centre;

// Jun 2009. Get rid of all the unmatched atoms on an atom. Current
// implementation not only breaks the bonds, but also removes the fragments
  
    resizable_array<int> _unlink_unmatched_atoms;

//  Jan 2017. Need the concept of a reaction that does nothing - just appends the name

    int _noop_reaction;

//  Sep 2017

    int _match_via_atom_map;

//  private functions

    int _write_msi (std::ostream & os, int &, const IWString & ind);
    int _parse_recompute_implicit_hydrogens (const msi_attribute *);
    int _parse_wedge_bonds (const msi_attribute *);

    int _read_query_file_from_molecule(IWString & fname);
    int _read_query_file(const IWString & fname, const int query_files_in_current_directory, const IWString & reaction_directory);

    int _read_inactive_query (const msi_attribute * att);
    int _read_inactive_query (const IWString & smarts);
    int _read_inactive_from_query_file (const IWString & fname);

    int _atom_slated_for_removal (const Matched_Atom_in_Component & a) const;
    int _atom_slated_for_removal (int a) const;

    int _do_invert_stereo_centres (Molecule & result,
                                          const Set_of_Atoms & embedding,
                                          int offset) const;
    int _do_remove_stereo_centres (Molecule & result,
                                          const Set_of_Atoms & embedding,
                                          int offset) const;

    int _remove_multiple_hits_that_do_not_involve_changing_atoms (Molecule &, Substructure_Results & sresults) const;

    int _common_matched_atoms_do_not_change (const Set_of_Atoms & e1,
                                             const Set_of_Atoms & e2,
                                             const int * symm) const;

    int _remove_multiple_hits_that_hit_atoms_being_changed (int, Substructure_Results & sresults) const;

    int _common_matched_atoms_change (const int * changed_by_other_embedding, const Set_of_Atoms & e2) const;

    int _discern_chiral_centres_to_be_saved_around_removed_atoms (Molecule & result,
                                          const Set_of_Atoms & embedding,
                                          int offset,
                                          int * chiral_centre_available);

    int _do_restore_saved_chiral_centre (Molecule & result,
                                            const Set_of_Atoms & scaffold_embedding,
                                            Enumeration_Temporaries & etmp,
                                            Chiral_Centre * c);

  protected:

//  int _do_toggle_kekule_form (Molecule & m, const Substructure_Results & sresults);

    int _determine_matched_atoms_checking_inactives (Molecule & m, Substructure_Results & sresults);


  public:
    Reaction_Site ();
    ~Reaction_Site ();

    int ok () const;
    int debug_print (std::ostream &, const const_IWSubstring &) const;

    int check_internal_consistency ();

    int unique_id () const { return _unique_id;}
    void set_unique_id (int s) { _unique_id = s;}

    int will_remove_atoms () const { return _atoms_to_be_removed.number_elements() + _fragments_to_be_removed.number_elements();}

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref);

    int construct_from_msi_object (const msi_object &, int, const IWString &);

    int determine_which_matched_atoms_are_changed ();
    int another_reagent_changes_your_matched_atom (int);

    void set_ignore_multiple_matches_involving_atoms_not_changing (int s) { _ignore_multiple_matches_involving_atoms_not_changing = s;}
    void set_ignore_multiple_matches_involving_changing_atoms (int s) { _ignore_multiple_matches_involving_changing_atoms = s;}

    void set_match_via_atom_map(const int s) { _match_via_atom_map = s;}

    void add_atom_to_be_removed (int a) { _atoms_to_be_removed.add (a);}    // should check for A being a valid number
    int add_bond_to_be_broken (int, int);
    int add_bond_to_be_made (int, int, bond_type_t);
    int add_fragment_to_be_removed (int);
    int add_chiral_centre_to_be_inverted (int a) { _stereo_centres_to_invert.add (a); return 1;}
    int add_replace_atom (Replace_Atom * r) { _replace_atom.add (r); return 1;}
    int add_element_to_be_changed (Reaction_Change_Element * e) { _elements_to_change.add (e); return 1;}
    int add_formal_charge_to_assign (int, formal_charge_t);
    int add_isotope_to_be_placed (int a, int iso);
    int add_isotope_to_be_incremented (int a, int iso);

    int has_saved_chiral_centres () const { return _saved_chiral_centre.number_elements ();}

    int do_makes_breaks (Molecule &, const Set_of_Atoms &, int, 
                         Enumeration_Temporaries & etmp);

    int set_dihedral_angles (Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Set_of_Atoms * sidechain_embedding,
                             const Set_of_Atoms & atoms_to_be_removed,
                             int atoms_already_present) const;
    int set_bond_lengths    (Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Set_of_Atoms * sidechain_embedding,
                             const Set_of_Atoms & atoms_to_be_removed,
                             int atoms_already_present) const;

    int do_unfix_implicit_hydrogens (Molecule & result, const Set_of_Atoms & embedding, int aap) const;

    int remove_and_invert_stereo_centres (Molecule & result, const Set_of_Atoms & embedding, int offset) const;

    int do_restore_saved_chiral_centres (Molecule & result, const Set_of_Atoms & embedding, Enumeration_Temporaries & etmp);

//  const Reaction_Stereo_Centre * stereo_centre (int i) const { return _reaction_stereo_centre[i];}

    int number_bond_length_specifications () const { return _reaction_bond_length.number_elements ();}
    const Reaction_Bond_Length * bond_length_specification (int i) const { return _reaction_bond_length[i];}

    int number_bond_angle_specifications () const { return _reaction_bond_angle.number_elements ();}
    const Reaction_Bond_Angle * bond_angle_specification (int i) const { return _reaction_bond_angle[i];}

    int number_dihedral_angle_specifications () const { return _reaction_dihedral_angle.number_elements ();}
    const Reaction_Dihedral_Angle * dihedral_angle_specification (int i) const { return _reaction_dihedral_angle[i];}

//  int number_rotate_fragment_specifications () const { return _reaction_rotate_fragment.number_elements ();}
//  const Reaction_Rotate_Fragment * rotate_fragment_specification (int i) const { return _reaction_rotate_fragment[i];}

    int number_3d_replace_specifications () const { return _reaction_3d_replace.number_elements ();}
    const Reaction_3D_Replace * r3d_replace_specification (int i) const { return _reaction_3d_replace[i];}

    int number_atom_replacements () const { return _replace_atom.number_elements ();}
    const Replace_Atom * replace_atom (int i) const { return _replace_atom[i];}

    int add_toggle_kekule_form (int, int, bond_type_t);
};

/*
  We differentiate between scaffold and sidechain reaction sites
*/

class Scaffold_Reaction_Site : public Reaction_Site
{
  protected:

  protected:

// protected functions

    int _do_intra_particle_replacements (Molecule & result,
                                          const Set_of_Atoms * scaffold_embedding) const;

// private functions

  public:
    Scaffold_Reaction_Site ();
    ~Scaffold_Reaction_Site ();

    int construct_from_msi_object (const msi_object & msi, int, const IWString &);
    int write_msi (std::ostream &, int &, int);
};

/*
  Sidechains are different from scaffolds for a number of reasons

    1. Sidechains own the points of attachment - that way we can have multiple
       sidechains join at once

    2. A sidechain_reaction_site may also have a SINGLE molecule which it
       attaches.
*/

class Sidechain_Reaction_Site : public Reaction_Site
{
  private:
    int _sidechain_number;

    resizable_array_p<Molecule_and_Embedding> _reagents;

    resizable_array_p<No_Reaction> _no_reaction;

//  how this is joined onto the reaction

    resizable_array_p<Inter_Particle_Bond> _inter_particle_bonds;

    int _make_implicit_hydrogens_explicit;

    Sidechain_Match_Conditions _match_conditions;

//  private functions

    int _copy_match_conditions_to_query ();

    int _check_inter_particle_atoms (int atoms_in_scaffold_query,
                            int scaffold_atom,
                            const resizable_array<int> & scaffold_atoms_to_be_removed,
                            const resizable_array<int> & scaffold_fragments_to_be_removed,
                            int sidechain_atom) const;

    int _add_reagents_from_smiles (const msi_attribute & att);
    int _add_reagents_from_file (const msi_attribute & att, const Sidechain_Match_Conditions & smc);

  public:
    Sidechain_Reaction_Site ();
    ~Sidechain_Reaction_Site ();

    int ok () const;
    int debug_print (std::ostream &, const const_IWSubstring &) const;

    int sidechain_number () const { return _sidechain_number;}
    void set_sidechain_number (int i) { _sidechain_number = i;}

    int adjust_matched_atoms_in_component (const extending_resizable_array<int> & xref);

    void set_make_implicit_hydrogens_explicit (int s) { _make_implicit_hydrogens_explicit = s;}

    int number_inter_particle_bonds () const { return _inter_particle_bonds.number_elements ();}
    const Inter_Particle_Bond * inter_particle_bond (int i) const { return _inter_particle_bonds[i];}

    int construct_from_msi_object (const msi_object & msi, int, const IWString &,
                                   const Sidechain_Match_Conditions & smc);

    int number_reagents () const { return _reagents.number_elements ();}

    int single_reagent_only () const { return 1 == _reagents.number_elements ();}

    int add_inter_particle_bond (int, int, int, bond_type_t);

    int remove_first_reagent_no_delete();

    int add_reagent  (Molecule_and_Embedding * m, const Sidechain_Match_Conditions & smc);
    int add_reagents (const char *, int, const Sidechain_Match_Conditions &);
    int add_reagents (data_source_and_type<Molecule_and_Embedding> &, const Sidechain_Match_Conditions &);
    int add_reagents_with_stringbuffer (const char * fname, int input_type,
                                       const Sidechain_Match_Conditions & smc, const char* stringbuffer, int stringbuffer_size, const char ** options, int optionsCount);

    int add_reagent_embedding_identified  (Molecule_and_Embedding * m, const Sidechain_Match_Conditions & smc);
    Molecule_and_Embedding * reagent (int r) const { return _reagents[r];}

    int empty_reagents_array ();

    int check_internal_consistency (int, int, const resizable_array<int> &, const resizable_array<int> &);

    int notify_scaffold_of_atoms_involved_in_changes (Reaction_Site & r) const;
    int notify_sidechain_of_atoms_involved_in_changes (Sidechain_Reaction_Site * r) const;

    int write_msi (std::ostream &, int &, int);

    int set_single_reagent (Molecule &);

    int do_makes_breaks (Molecule & result,
                         const Set_of_Atoms * embedding,
                         int offset,
                         Enumeration_Temporaries & etmp);

    int make_inter_particle_substitutions (Molecule & result,
                                          const Set_of_Atoms * scaffold_embedding,
                                          int atoms_in_scaffold,
                                          const Set_of_Atoms * sidechain_embedding);
    int make_inter_particle_substitutions (Molecule & result,
                                   const Set_of_Atoms * scaffold_embedding,
                                   int atoms_in_scaffold);
    int make_inter_particle_bonds (Molecule & result,
                                   const Set_of_Atoms * scaffold_embedding,
                                   int atoms_in_scaffold,
                                   const Set_of_Atoms * sidechain_embedding);
    int make_inter_particle_bonds (Molecule & result,
                                   const Set_of_Atoms * scaffold_embedding,
                                   int atoms_in_scaffold);
//  int make_stereo_centres (Molecule & result,
//                           const Set_of_Atoms * scaffold_embedding,
//                           const Set_of_Atoms * sidechain_embedding,
//                           const Set_of_Atoms & atoms_to_be_removed,
//                           int atoms_already_present) const;
    int do_inter_particle_replacements (Molecule & result,
                                          const Set_of_Atoms * scaffold_embedding,
                                          int atoms_in_scaffold,
                                          const Set_of_Atoms * sidechain_embedding) const;
};

/*
  A Reaction_Iterator describes just which reagent from each sidechain is currently
  being used
*/

class IWReaction;
class Molecule_Output_Object;

class Reaction_Iterator
{
  private:

    int _number_sidechains;

    int * _reagents_in_sidechain;

    int * _reagent;

    int _active;

  public:
    Reaction_Iterator ();
    ~Reaction_Iterator ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int initialise (const IWReaction &);

    int active () const { return _active;}

    void operator ++ (int);

    int reagent (int i) const { assert (i >= 0 && i < _number_sidechains); return _reagent[i];}

#if (GCC_VERSION >= 40900)
    void randomise (std::mt19937_64 & rng);
#endif
};

extern std::ostream & operator << (std::ostream &, const Reaction_Iterator &);

/*
  A reaction is derived (ultimately) from the query which defines
  how the scaffold grows
*/

class IWReaction : public Scaffold_Reaction_Site
{
  protected:
    IWString _comment;

    resizable_array_p<Sidechain_Reaction_Site> _sidechains;  // all the various things which get added

//  Although reagents are owned by various sidechains, it is often convenient
//  to have all the reagents in a given order - this was originally driven by
//  the parallel_synthesis application, but it is a good idea generally

    resizable_array<Molecule_and_Embedding *> _reagent;

// Sept 2003.  Move the stereo centre stuff to the IWReaction object. 
// This is because the atoms involved in the chiral centre may involve
// atoms in the scaffold and/or any number of sidechains

    resizable_array_p<Reaction_Stereo_Centre> _reaction_stereo_centre;

//  Do we want names appended or not

    int _append_names;

//  We can also have a specific string to append to the name

    IWString _append_to_name;

//  When looking for a query file, do we look in the current directory or
//  in the same directory as the reaction file

    int _query_files_in_current_directory;

//  In that case, we will need the directory

    IWString _reaction_directory;

//  May 2002. We may have a Kekule form problem. We can try to automatically fix the problem

    int _find_kekule_forms_for_bad_valence;

//  Passed to all the sidechains

    int _make_implicit_hydrogens_explicit;

    Scaffold_Match_Conditions _match_conditions;

    void * _user_specified_void_ptr;

//  When reading a reaction, are the transformation numbers based on matched atom
//  indices or atom map numbers

    int _match_via_atom_map;

//  private functions

    int _construct_from_msi_object (const msi_object &,
                                    const Sidechain_Match_Conditions & sidechain_match_conditions);
    int _build_from_smirks ();
    int _from_smirks_add_inter_particle_bond_involves_scaffold (int a1, int a2, bond_type_t bt);
    int _from_smirks_add_inter_particle_bond(int a1, 
                                             Reaction_Site * a1site,
                                             int a2,
                                             Reaction_Site * a2site,
                                             bond_type_t bt);

    int _add_molecule (Molecule & result, const Molecule & added);

    int _perform_scaffold_transformations (Molecule & result, const Set_of_Atoms * scaffold_embedding, const Enumeration_Temporaries & etmp) const;

    int _do_atom_removals (Molecule & result, Enumeration_Temporaries & etmp) const;

    int _perform_reaction_recheck_scaffold_reactivity (const Molecule * scaffold,
                                                           const Substructure_Results & sresults,
                                                           Molecule & result,
                                                           Enumeration_Temporaries & etmp);
    int _perform_reaction (Molecule &, const Sidechain_Reaction_Site *,
                           const Set_of_Atoms *,
                           const Molecule *, const Set_of_Atoms *,
                           Set_of_Atoms &, resizable_array<const Atom *> &);
    int _perform_reaction (Molecule &, const Set_of_Atoms *,
                           const Molecule *, const Set_of_Atoms *);

    int _perform_reaction (const Molecule * scaffold,
                               const Set_of_Atoms * e, 
                               Sidechain_Reaction_Site * s,
                               Molecule_and_Embedding * reagent,
                               Molecule & result);

    int _examine_possible_reagent (Molecule_and_Embedding *, int &, int verbose);

    void _copy_iterator_data_to_etmp (const Reaction_Iterator & iter,
                                         Enumeration_Temporaries & etmp) const;

    int _take_first_reagent_from_each_sidechain (Enumeration_Temporaries & etmp) const;

    int _in_place_transformation (Molecule & m,
                                  const Set_of_Atoms * scaffold_embedding,
                                  Enumeration_Temporaries & etmp);

    int _perform_reaction (Molecule & result,
                               const Set_of_Atoms * scaffold_embedding,
                               Enumeration_Temporaries & etmp);

//  int _assemble_reagents (const Reaction_Iterator & riter,
//                          Molecule_and_Embedding ** reagent) const;

    int _do_replacement (Molecule & result,
                             const Set_of_Atoms * scaffold_embedding,
                             const Replace_Atom & r,
                             const Enumeration_Temporaries & etmp) const;
    int _do_replacements (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              Enumeration_Temporaries & etmp) const;
    int _make_inter_particle_bond (Molecule & result,
                                       const Set_of_Atoms * scaffold_embedding,
                                       const Inter_Particle_Bond & b,
                                       const Enumeration_Temporaries & etmp) const;
    int _make_inter_particle_bonds (Molecule & result,
                                        const Set_of_Atoms * scaffold_embedding,
                                        const Enumeration_Temporaries & etmp) const;
    int _make_stereo_centres (Molecule & result,
                                  const Set_of_Atoms * scaffold_embedding,
                                  const Enumeration_Temporaries & etmp) const;
    int _set_bond_length (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Reaction_Bond_Length & b,
                              const Enumeration_Temporaries & etmp) const;
    int _set_bond_lengths (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Reaction_Bond_Length & b,
                              const Enumeration_Temporaries & etmp,
                              const Reaction_Site &) const;
    int _set_bond_lengths (Molecule & result,
                               const Set_of_Atoms * scaffold_embedding,
                               const Enumeration_Temporaries & etmp,
                               const Reaction_Site & r) const;
    int _set_bond_lengths (Molecule & result,
                               const Set_of_Atoms * scaffold_embedding,
                               const Enumeration_Temporaries & etmp) const;
    int _set_dihedral_angles (Molecule & result,
                                  const Set_of_Atoms * scaffold_embedding,
                                  const Enumeration_Temporaries & etmp) const;
   int _set_dihedral_angles (Molecule & result,
                                  const Set_of_Atoms * scaffold_embedding,
                                  const Enumeration_Temporaries & etmp,
                                  const Reaction_Site & r) const;
    int _set_bond_angles (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp) const;
    int _set_bond_angles (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp,
                              const Reaction_Site & r) const;

    int _do_3d_replacements (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp,
                              const Reaction_Site & r) const;
    int _do_3d_replacements (Molecule & result,
                              const Set_of_Atoms * scaffold_embedding,
                              const Enumeration_Temporaries & etmp) const;

    int _do_find_kekule_forms_for_bad_valence (Molecule & result) const;

    int _break_bonds_to_atoms_and_fragments_to_be_removed (Molecule & result, const Enumeration_Temporaries & etmp) const;

    int _do_restore_saved_chiral_centres (Molecule & result,
                                          const Set_of_Atoms & scaffold_embedding,
                                          Enumeration_Temporaries & etmp);

    int _determine_which_matched_atoms_are_changed ();
    int _convert_from_msi_object_numbers_to_internal_indices ();

    int _identify_changes_from_smirks (const resizable_array<int> & a);
    int _identify_changes_from_smirks (const extending_resizable_array<int> & a, const int istop, const Molecule & orphan_molecule, const Substructure_Query &);
    int _discern_atomic_changes(Reaction_Site &, const Substructure_Atom & , const Substructure_Atom &);
    int _discern_atomic_changes_specifier(Reaction_Site & r,
                                    const atom_number_t,
                                    const Substructure_Atom_Specifier & q1,
                                    const Substructure_Atom_Specifier & q2);

    int _create_orphan_molecule(const resizable_array<int> & orphan_atoms,
                                Molecule & orphan_molecule,
                                    const Substructure_Query & product_molecule,
                                    Sidechain_Reaction_Site & sc);

    int _construct_from_smirks(const const_IWSubstring & smirks);

    const Substructure_Atom * query_atom_with_initial_atom_number (atom_number_t) const;
    Reaction_Site * _reaction_site_with_initial_atom_number (atom_number_t z);
    const Substructure_Bond * bond_between_atoms (atom_number_t, atom_number_t) const;
    Reaction_Site * _component_with_bond (const atom_number_t a1, const atom_number_t a2);
    Reaction_Site * _reaction_site_with_atom_map_number(const int amap);
    Reaction_Site * _component_with_bond_between_mapped_atoms (const int a1, const int a2);
    int _sidechain_with_mapped_atom(const int x) const;

  public:
    IWReaction ();
    ~IWReaction ();

    int ok () const;
    int debug_print (std::ostream &) const;

    const IWString & comment () const { return _comment;}
    void set_comment (const IWString & c) { _comment = c;}
    const IWString & name () const { return _comment;}

    int number_reagents () const;

    int number_products_per_scaffold_embedding () const;

    int number_sidechains () const { return _sidechains.number_elements ();}

    int set_number_sidechains (int s);

    Sidechain_Reaction_Site * sidechain (int s) const { return _sidechains[s];}

    int number_sidechains_with_reagents () const;

    int all_sidechains_have_reagents () const;

    int construct_from_msi_object (const msi_object &);

    int construct_from_msi_object (const msi_object &, const Sidechain_Match_Conditions & sidechain_match_conditions);

    int construct_from_smirks (const const_IWSubstring & smirks);
    int construct_from_smirks (const const_IWSubstring & reagents, const const_IWSubstring & products);

    int write_msi (const char * fname);
    int write_msi (std::ostream &);

    int check_internal_consistency ();

    int setup_to_skip_multiple_embeddings_involving_non_changing_atoms ();

    void set_append_names (int s) { _append_names = s;}

    Scaffold_Match_Conditions & scaffold_match_conditions () { return _match_conditions;}
    const Scaffold_Match_Conditions & scaffold_match_conditions () const { return _match_conditions;}

    int do_read (iwstring_data_source &, const Sidechain_Match_Conditions &);
    int do_read (const const_IWSubstring &, const Sidechain_Match_Conditions &);
    int do_read (const char *, const Sidechain_Match_Conditions &);
    int do_read (const IWString &, const Sidechain_Match_Conditions &);
    int do_read_with_stringbuffer (const const_IWSubstring & fname, const Sidechain_Match_Conditions & sidechain_match_conditions, const char* stringbuffer, int stringbuffer_size, const char ** options, int optionsCount);

    void set_one_embedding_per_start_atom(const int s);

    void set_query_files_in_current_directory (int q) { _query_files_in_current_directory = q;}

    void set_find_kekule_forms_for_bad_valence (int f) { _find_kekule_forms_for_bad_valence = f;}
    void set_make_implicit_hydrogens_explicit (int s) { _make_implicit_hydrogens_explicit = s;}

    int add_reaction_stereo_centre (Reaction_Stereo_Centre * r) { _reaction_stereo_centre.add (r); return 1;}
    int number_stereo_centres () const { return _reaction_stereo_centre.number_elements ();}

//  When determining matches, we need to consider not only the underlying Substructure_Query object,
//  but also any _inactive queries

    int determine_matched_atoms (Molecule &, Substructure_Results &);

    Molecule_and_Embedding * reagent (const Reaction_Iterator &) const;

//  Sometimes we are just making changes to a single molecule and don't
//  need to preserve the "scaffold" for the next reagent.

    int in_place_transformation  (Molecule & m, const Set_of_Atoms * scaffold_embedding);
    int in_place_transformation  (Molecule & m, const Set_of_Atoms * scaffold_embedding, const Reaction_Iterator &);
    int in_place_transformations (Molecule & m, const Substructure_Results & sresults);
    int in_place_transformations (Molecule & m);

    int perform_reaction (Molecule & m, const Substructure_Results & sresults);

    int perform_reaction (const Molecule * scaffold, const Substructure_Results & sresults, Molecule & result);
    int perform_reaction (const Molecule * scaffold, const Set_of_Atoms * embedding, Molecule & result);

    int perform_reaction (const Molecule * scaffold, const Substructure_Results & sresults, const Reaction_Iterator &, Molecule & result);
    int perform_reaction (const Molecule * scaffold, const Set_of_Atoms * embedding, const Reaction_Iterator &, Molecule & result);
    int perform_reaction_recheck_scaffold_reactivity (const Molecule * scaffold, const Substructure_Results & sresults, const Reaction_Iterator & iterator, Molecule & result);
    int perform_reaction_recheck_scaffold_reactivity (const Molecule * scaffold, const Substructure_Results & sresults, Molecule & result);

//  Used by parallel synthesis to access the I'th reagent - regardless of which sidechain owns it

    int perform_reaction (const Molecule * scaffold, const Set_of_Atoms * embedding, int, Molecule & result);

    int perform_reaction (const Molecule *, const Set_of_Atoms *,
                          const Molecule *, const Set_of_Atoms *,
                          Molecule &);

    int perform_reaction (Molecule &, resizable_array_p<Molecule> &);
//  int perform_reaction (Molecule & m, Molecule_Output_Object &);
//  int perform_reaction (Molecule & m, const Set_of_Atoms * embedding, Molecule_Output_Object &);

    void   set_user_specified_void_ptr(void * v) { _user_specified_void_ptr = v;}
    void * user_specified_void_ptr() const { return _user_specified_void_ptr;}
};

extern void set_strip_reagents_to_largest_fragment (int);
extern void set_component_separator (const const_IWSubstring &);
extern void set_reaction_transfer_text_info_to_products (int);
extern int set_stream_for_sidechains_not_matching_query (const const_IWSubstring & fname);

extern int
determine_atom_number (const Set_of_Atoms & scaffold_embedding,
                       const Matched_Atom_in_Component & ma,
                       const Enumeration_Temporaries & etmp,
                       const char * caller,
                       atom_number_t & zresult);

int
read_reactions (const Command_Line & cl,
                resizable_array_p<IWReaction> & rxn,
                Sidechain_Match_Conditions & smc,
                char flag);

extern int
count_number_of_reactions (const Command_Line & cl,
                           char flag);

/*
  When things go wrong, it can be handy to know which atoms
  are problematic
*/

extern int write_isotopically_labelled_smiles(Molecule &, std::ostream &);

extern void set_iwreaction_display_no_atoms_in_query_message (int s);

extern void set_smirks_lost_atom_means_remove_frgment(int s);

extern void set_iwreaction_display_take_first_reagent_from_each_sidechain_warning_message(const int s);

#endif

/* arch-tag: faea2f4a-364f-4841-b65b-158d34ea3a17

*/
