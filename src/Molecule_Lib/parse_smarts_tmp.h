#ifndef MOLECULE_LIB_PARSE_SMARTS_TMP_H_
#define MOLECULE_LIB_PARSE_SMARTS_TMP_H_

/*
  When parsing a smarts, we need a lot of things passed around,
  so rather than having functions with way too many arguments
  we stick things into an object
*/

// During parsing, three dots specifications consist of
// two query atoms and a possible qualification, like number
// of atoms between, or smarts.

#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"

// A base class for three dots and down the bond

class AtomsAndQualifier
{
  protected:
    int _matched_atom_1;
    int _matched_atom_2;
    IWString _qualifier;

  public:
    AtomsAndQualifier() {
      _matched_atom_1 = -1;
      _matched_atom_2 = -1;
    }
    AtomsAndQualifier(const int a1, const int a2) : _matched_atom_1(a1), _matched_atom_2(a2)
    {
    }

    int a1() const { return _matched_atom_1;}
    int a2() const { return _matched_atom_2;}

    void set_qualifier(const const_IWSubstring& s) {
      _qualifier = s;
    }

    const IWString& qualifier() const { return _qualifier;}
};

class ThreeDots : public AtomsAndQualifier {
  private:
  public:
    ThreeDots(const int a1, const int a2);

    int BuildFromProto(const SubstructureSearch::NoMatchedAtomsBetween& proto);
};

class DownTheBondInfo : public AtomsAndQualifier {
  private:
  public:
    DownTheBondInfo(int previous_atom);

    void set_a2(int s) {
      _matched_atom_2 = s;
    }

    int Build(const const_IWSubstring& buffer);
};

class Parse_Smarts_Tmp
{
  private:
    int _last_query_atom_created;

    resizable_array<Substructure_Atom *> _root;
    resizable_array<Bond *> _no_matched_atoms_between;
    resizable_array<Link_Atom *> _link_atom;
    resizable_array<ThreeDots *> _three_dots;
    resizable_array<down_the_bond::DownTheBond *> _down_the_bond;

    extending_resizable_array<Substructure_Atom *> _completed;

  public:
    Parse_Smarts_Tmp ();
    ~Parse_Smarts_Tmp ();

    int set_natoms (int);

    void add_root_atom (Substructure_Atom * r) { _root.add (r);}
    void add_no_matched_atoms_between (Bond * r) { _no_matched_atoms_between.add (r);}
    void add_link_atom (Link_Atom * r) { _link_atom.add (r);}
    void add_three_dots(ThreeDots* three_dots) {
      _three_dots.add(three_dots);
    }
    void add_down_the_bond(down_the_bond::DownTheBond* dtb) {
      _down_the_bond.add(dtb);
    }

    extending_resizable_array<Substructure_Atom *> & completed () { return _completed;}

    void set_last_query_atom_created (int s) { _last_query_atom_created = s;}
    int last_query_atom_created () const { return _last_query_atom_created;}

    const resizable_array<Substructure_Atom *> & root_atoms () const { return _root;}
    const resizable_array<Bond *> & no_matched_atoms_between () const { return _no_matched_atoms_between;}
    const resizable_array<Link_Atom *> & link_atoms () const { return _link_atom;}
    const resizable_array<ThreeDots *> & three_dots () const { return _three_dots;}
    const resizable_array<down_the_bond::DownTheBond *> & down_the_bond () const { return _down_the_bond;}

    resizable_array<Substructure_Atom *> & root_atoms () { return _root;}
    resizable_array<Bond *> & no_matched_atoms_between () { return _no_matched_atoms_between;}
    resizable_array<Link_Atom *> & link_atoms () { return _link_atom;}
    resizable_array<ThreeDots *> & three_dots () { return _three_dots;}
    resizable_array<down_the_bond::DownTheBond *> & down_the_bond () { return _down_the_bond;}

    // An error has occurred, and the smarts cannot be parsed. This object holds a bunch
    // of things that have been allocated and not yet transferred to the
    // Substructure_Query object. so they need to be manually deleted.
    void DeletePointers();

    // Check that each DownTheBond object is Complete()
    int DownTheBondSpecificationsComplete() const;
};

#endif  // MOLECULE_LIB_PARSE_SMARTS_TMP_H_
