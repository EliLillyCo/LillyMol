#ifndef MOLECULAR_ABS_SPECS_H
#define MOLECULAR_ABS_SPECS_H

#include <iostream>

#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/molecule.h"

class Molecule_With_Info_About_Parent : public Molecule
{
  private:
    int _parent_natoms;

  public:
    Molecule_With_Info_About_Parent() : _parent_natoms(0) {};

    void compute_parent_natoms() { _parent_natoms = Molecule::natoms();}

    int parent_natoms () const { return _parent_natoms;}
};

/*
  We need a language for defining molecular abstractions. That
  way people have complete flexibility in what they can do.

  The language looks like

  token1(args).token2(args)

  args are optional directives for the preceeding token.

  We can have branching

  token1(args)(token2(args))token3(args)

  In this case, we get a tree rooted at token1 with two children
*/

// We need string and integer versions of these

#define MAD_SCAFFOLD "scaffold"
#define MAD_TYPE_SCAFFOLD 1

#define MAD_REPLACE_LINKER "rplink"
#define MAD_TYPE_REPLACE_LINKER 2

#define MAD_RMSCAFFOLD "rmscaffold"
#define MAD_TYPE_RMSCAFFOLD 3

#define MAD_RINGS "rings"
#define MAD_TYPE_RINGS 4

#define MAD_BIGRING "bigring"
#define MAD_TYPE_BIGRING 5

#define MAD_TRANSLATE "translate"
#define MAD_TYPE_TRANSLATE 6

#define MAD_RMSPINACH "rmspinach"
#define MAD_TYPE_RMSPINACH 7

#define MAD_REMOVE_ATOM "rmat"
#define MAD_TYPE_REMOVE_ATOM 8

#define MAD_CHANGE_BOND_TYPE "cbt"
#define MAD_TYPE_CHANGE_BOND_TYPE 9

#define MAD_ABSTRACT_RING_FORM "arf"
#define MAD_TYPE_ABSTRACT_RING_FORM 10

#define MAD_FRAGMENT "frag"
#define MAD_TYPE_FRAGMENT 11

#define MAD_REMOVE_ATOMS "rmatoms"
#define MAD_TYPE_REMOVE_ATOMS 12

#define MAD_PLACE_ISOTOPE "isotope"
#define MAD_TYPE_PLACE_ISOTOPE 13

#define MAD_PLACE_CHARGE "charge"
#define MAD_TYPE_PLACE_CHARGE 14

#define MAD_ALL_ATOMS_TRANSFORM "allatoms"
#define MAD_TYPE_ALL_ATOMS_TRANSFORM 15

#define MAD_ALL_BONDS_TRANSFORM "allbonds"
#define MAD_TYPE_ALL_BONDS_TRANSFORM 16

#define MAD_COMPRESS_CONSECUTIVE "comprconsec"
#define MAD_TYPE_COMPRESS_CONSECUTIVE 17

#define MAD_RINGSYS "ringsys"
#define MAD_TYPE_RINGSYS 18

#define MAD_RMBOND "rmbond"
#define MAD_TYPE_RMBOND 19

#define MAD_RMRD2 "rmrd2"
#define MAD_TYPE_RMRD2 20

#define MAD_INVSCAF "invscaf"
#define MAD_TYPE_INVSCAF 21

#define MAD_SSS "sss"
#define MAD_TYPE_SSS 22

#define MAD_SPINACH "spinach"
#define MAD_TYPE_SPINACH 23

extern void DisplayAbstractionNames(std::ostream& output);
extern void DisplayUsageExamples(std::ostream& output);

class Molecular_Abstraction_Directives_Node
{
  private:
    int _type;

    IWString _directive;
    IWString _args;

    Molecular_Abstraction_Directives_Node * _next;

// private functions

    int _build (const const_IWSubstring &);
    int _build_children (const const_IWSubstring &);
    int _extract_directive_and_args(const const_IWSubstring & s);

  public:
    Molecular_Abstraction_Directives_Node();
    ~Molecular_Abstraction_Directives_Node();

    int build(const const_IWSubstring &);

    const IWString & args() const { return _args;}

    int debug_print(std::ostream &) const;

    int ztype() const { return _type;}

    int directive_recognised();

    int number_abstractions() const;

    const Molecular_Abstraction_Directives_Node * next() const { return _next;}
};

#endif
