#ifndef MOLECULAR_ABSTRACTION_FUNCTIONS_H
#define MOLECULAR_ABSTRACTION_FUNCTIONS_H

#include "substructure.h"

#include "molecular_abstraction_specifications.h"

class Molecular_Abstraction_Base_Class
{
  protected:
    int _molecules_processed;
    int _molecules_changed;

    IWString _smiles_tag;
    IWString _fingerprint_tag;
    IWString _write_tag;

    int _nbits;

    int _isotope;

    int _append_count_to_tag;

    int _write_only_if_changes;

    int _min_atoms_needed_for_write;
    int _max_atoms_allowed_for_write;

    float _min_atom_ratio_needed_for_write;
    float _max_atom_ratio_allowed_for_write;

    Atom_Typing_Specification _atom_typing_specification;

//  protected functions

    int _process(const const_IWSubstring &, const char *, int &);

    int _do_any_writing_needed(Molecule_With_Info_About_Parent &, int, IWString_and_File_Descriptor &);
    int _handle_no_match_to_query(Molecule_With_Info_About_Parent & m, IWString_and_File_Descriptor & output);

    int _identify_scaffold (Molecule_With_Info_About_Parent &, int *, int) const;
    int _is_spinach(Molecule_With_Info_About_Parent & m,
                     int * in_scaffold,
                     atom_number_t aprev,
                     atom_number_t zatom) const;


  private:
    int _parse_write_directive(const const_IWSubstring & token, const char *);
    int _parse_fp_directive(const const_IWSubstring & token, const char *);
    int _parse_at_directive (const const_IWSubstring & token, const char * caller);
    int _parse_isotope_directive(const const_IWSubstring & token, const char * caller);

  public:
    Molecular_Abstraction_Base_Class();
    virtual ~Molecular_Abstraction_Base_Class();

    int ok() const;

    int report(ostream &) const;

    const IWString & fingerprint_tag() const { return _fingerprint_tag;}
    const IWString & write_tag() const { return _write_tag;}

    virtual int build (const Molecular_Abstraction_Directives_Node &) = 0;

    virtual int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &) = 0;
};

class Molecular_Abstraction_Transform : public Molecular_Abstraction_Base_Class
{
  private:
    resizable_array_p<Single_Substructure_Query> _smarts;
    resizable_array<const Element *> _eto;

  public:
    Molecular_Abstraction_Transform();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_All_Transform : public Molecular_Abstraction_Base_Class
{
  private:
    const Element * _eto;

//  private functions

  public:
    Molecular_Abstraction_All_Transform();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Remove_Atom : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;
    int _rejoin_all;
    extending_resizable_array<int> _rejoin;

//  private functions

    int _do_removals_with_rejoins(Molecule_With_Info_About_Parent & m, int * to_remove) const;

  public:
    Molecular_Abstraction_Remove_Atom();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

/*
  Proper name should be ...Remove_Atoms, but that name is too close to 
  the one that just works on a single atom, so we create a different
  sounding name
*/

class Molecular_Abstraction_Delete_Atoms : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;
    int _rejoin_all;
    extending_resizable_array<int> _rejoin;

//  private functions

    int _do_removals_with_rejoins(Molecule_With_Info_About_Parent & m, int * to_remove) const;

  public:
    Molecular_Abstraction_Delete_Atoms();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Largest_Ring_System : public Molecular_Abstraction_Base_Class
{
  private:
    int _spiro;   // do we include spiro fusions or not?

  public:
    Molecular_Abstraction_Largest_Ring_System();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Rings : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;
    const Element * _eto;

  public:
    Molecular_Abstraction_Rings();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Scaffold : public Molecular_Abstraction_Base_Class
{
  int _keep_first_ring_attachment;

  private:

//  private function

  public:
    Molecular_Abstraction_Scaffold();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Change_Bond_Type : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;
    bond_type_t _bt;

  public:
    Molecular_Abstraction_Change_Bond_Type();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Change_All_Bonds : public Molecular_Abstraction_Base_Class
{
  private:
    bond_type_t _bt;

  public:
    Molecular_Abstraction_Change_All_Bonds();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};


class Molecular_Abstraction_Replace_Linker : public Molecular_Abstraction_Base_Class
{
  private:
    const Element * _linker_atom;

  public:
    Molecular_Abstraction_Replace_Linker();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Joins_to_Abstract_Ring;

class Molecular_Abstraction_Abstract_Ring_Form : public Molecular_Abstraction_Base_Class
{
  private:
    const Element * _arom_ele;
    const Element * _aliph_ele;
    int _label_by_ring_size;
    int _label_by_heteroatom_count;

//  private functions

    int _spread_apart_any_spiro_rings(Molecule_With_Info_About_Parent & m, int * tmp) const;

    int _place_abstract_rings (Molecule_With_Info_About_Parent & m,
                                  const int * ring_membership,
                                  Joins_to_Abstract_Ring * jar);
    void _do_apply_isotopic_label (Molecule_With_Info_About_Parent & m,
                                        atom_number_t zatom,
                                        const Joins_to_Abstract_Ring & jar) const;

  public:
    Molecular_Abstraction_Abstract_Ring_Form();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Fragment : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _fragment_to_keep;
    Single_Substructure_Query _fragment_to_remove;

//  private functions

    int _identify_fragments_hit (Molecule_With_Info_About_Parent & m, Single_Substructure_Query & q,
                                        int * frag);
  public:
    Molecular_Abstraction_Fragment();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};


class Molecular_Abstraction_Place_Isotope : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;

//  ideally this would be a hash_map<int, int>, but I just didn't want
//  to include hashes for such a simple thing.

    resizable_array<int> _matched_atom;
    resizable_array<int> _isotope;

//  private functions

    void _initialise_isotope_on_first_matched_atom();

  public:
    Molecular_Abstraction_Place_Isotope();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Place_Charge : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;

//  ideally this would be a hash_map<int, int>, but I just didn't want
//  to include hashes for such a simple thing.

    resizable_array<int> _matched_atom;
    resizable_array<formal_charge_t> _charge;

//  private functions

    void _initialise_charge_on_first_matched_atom();

  public:
    Molecular_Abstraction_Place_Charge();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};


class Molecular_Abstraction_Compress_Consecutive : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;

//  private functions

    int _compress_group(Molecule_With_Info_About_Parent & m, int * to_remove,
               int flag, resizable_array<const Bond *> & join_points) const;

  public:
    Molecular_Abstraction_Compress_Consecutive();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};


class Molecular_Abstraction_Ring_Systems : public Molecular_Abstraction_Base_Class
{
  private:

    const Element * _ele;

//  private functions

    int _replace_ring_system (Molecule_With_Info_About_Parent & m, 
                     int matoms,
                     int flag,
                     const int * ring_systems) const;

  public:
    Molecular_Abstraction_Ring_Systems();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Remove_Bond : public Molecular_Abstraction_Base_Class
{
  private:
    Single_Substructure_Query _smarts;

    int _remove_fragment;

//  private functions

  public:
    Molecular_Abstraction_Remove_Bond();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Remove_Ring_CH2 : public Molecular_Abstraction_Base_Class
{
  private:
    const Element * _element;

// private functions

    int _identify_ch2_to_remove (Molecule_With_Info_About_Parent & m, atom_number_t zatom,
                        Set_of_Atoms & lhs,
                        Set_of_Atoms & rhs) const;
    int _identify_sequence (Molecule_With_Info_About_Parent & m, atom_number_t avoid,
                                                atom_number_t zatom,
                                                Set_of_Atoms & s) const;

  public:
    Molecular_Abstraction_Remove_Ring_CH2();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Inverse_Scaffold : public Molecular_Abstraction_Base_Class
{
  private:
    const Element * _scaffold_chain_element;

// private functions

  public:
    Molecular_Abstraction_Inverse_Scaffold();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Substructure_Search : public Molecular_Abstraction_Base_Class
{
  private:
    resizable_array_p<Substructure_Query> _queries;

    int _must_hit_at_least_one_query;

// private functions

  public:
    Molecular_Abstraction_Substructure_Search();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};

class Molecular_Abstraction_Spinach : public Molecular_Abstraction_Base_Class
{
  private:

    int _remove_doubly_bonded_atoms_in_spinach;

    const Element * _aromatic_element_replacement;
    const Element * _aliphatic_element_replacement;
    const Element * _chain_element_replacement;

// private functions

  public:
    Molecular_Abstraction_Spinach();

    int debug_print(ostream &) const;

    int build(const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent &, IWString_and_File_Descriptor &);
};


class Set_of_Molecular_Abstractions
{
  private:
    int _n;
    Molecular_Abstraction_Base_Class ** _a;

  public:
    Set_of_Molecular_Abstractions();
    ~Set_of_Molecular_Abstractions();

    int report(ostream &) const;

    int build (const Molecular_Abstraction_Directives_Node &);

    int process(Molecule_With_Info_About_Parent & m, IWString_and_File_Descriptor & output);

    int what_is_being_written(int &, int &) const;   // smiles or fingerprints
};

extern void set_write_only_if_changes (int s);
extern void set_append_count_to_tag (int s);
extern void set_remove_invalid_chiral_centres_before_writing (int s);
extern void set_write_empty_molecule_on_no_match(const int s);

#endif
