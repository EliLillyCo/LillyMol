#ifndef MDL_FUNCTIONS_H
#define MDL_FUNCTIONS_H

#include "iwstring.h"
#include "iwmtypes.h"

#include "atom.h"
#include "atom_alias.h"
#include "set_of_atoms.h"

/*
  There are a number of attributes and optional behaviours associated with
  reading and writing MDL files. We have a single object that holds all
  those descriptions, and performs some of the work involved
*/

class MDL_File_Supporting_Material
{
  private:
    int _isis_standard_records;
    int _ignore_unrecognised_m_records;
    int _report_unrecognised_records;
    int _die_on_erroneous_m_input;
    int _write_mdl_dollar_signs;
    int _write_mdl_m_end_record;
    int _write_v30_mdl_files;
    int _ignore_self_bonds;
    int _write_mdl_chiral_flags;
    int _include_chiral_info_in_mdl_outputs;
    int _mdl_read_h_correct_chiral_centres;
    int _mdl_write_h_correct_chiral_centres;
    int _extract_isis_extregno;
    int _fetch_all_sdf_identifiers;
    int _take_first_tag_as_name;
    int _prepend_sdfid;
    int _discard_sdf_molecule_name;
    int _multi_record_tag_data_present;
    int _mdl_write_aromatic_bonds;
    int _mdl_write_aromatic_atoms;
    int _display_non_organic_chirality_messages;
    int _mdl_display_invalid_chiral_connectivity;
    int _truncate_long_symbols;
    int _discern_chirality_from_wedge_bonds;
    int _write_isotopes_as_numbers_rather_than_differences_from_normal;
    int _write_M_isotopes_as_numbers_rather_than_differences_from_normal;
    int _write_fixed_width_m_iso_fields;
    int _read_isotopes_as_numbers_rather_than_differences_from_normal;
    int _read_M_isotopes_as_numbers_rather_than_differences_from_normal;
    int _allow_deuterium;
    int _allow_tritium;
    int _mdl_g_records_hold_atom_symbols;
    int _write_mdl_charges_as_m_chg;
    int _set_elements_based_on_atom_aliases;
    int _write_Rn_groups_as_element;
    int _convert_single_atom_sgroup_to_element;
    IWString _change_long_symbols_to;
    IWString _insert_between_sdf_name_tokens;
    IW_Regular_Expression _sdf_identifier;
    IWString _name_in_m_tag;
    IWString _replace_first_sdf_tag;
    int * _input_bond_type_translation_table;

    IWString * _digits2;
    IWString * _digits3;

//  Stuff for RDfiles

    resizable_array_p<IWString> _rdfile_identifiers;
    IWString _rdfile_start_of_record;

    char _gsub_mdl_file_data;      // drive out spaces

//  I initially implemented a separate object to hold these items, but it just made the code
//  messier, so we just put them here. Not ideal, but the whole idea of silently accumulating
//  various things somewhere is horrible. So this is not that great either...

    int _accumulate_mdl_chirality_features;

    Set_of_Atoms _unspecified_chiral_atoms_last_molecule;
    Set_of_Atoms _down_bonds_last_molecule;
    Set_of_Atoms _up_bonds_last_molecule;
    Set_of_Atoms _squiggle_bonds_last_molecule;
    Set_of_Atoms _unspecified_double_bond_atoms;

    int _a_records_found;
    int _g_records_found;
    resizable_array_p<Atom_Alias> _aliases;

// private functions

    void _default_values ();
    void _initialise_digits ();

  public:
    MDL_File_Supporting_Material ();
    ~MDL_File_Supporting_Material ();

    void reset_for_next_molecule ();

    void set_write_isis_standard (int s) { _isis_standard_records = s;}
    void set_ignore_unrecognised_mdl_m_records (int s) { _ignore_unrecognised_m_records = s;}
    void set_report_unrecognised_records (int s) { _report_unrecognised_records = s;}
    void set_die_on_erroneous_m_input (int s) { _die_on_erroneous_m_input = s;}
    void set_write_mdl_dollars (int s) { _write_mdl_dollar_signs = s;}
    void set_write_mdl_m_end_record (int s) { _write_mdl_m_end_record = s;}
    void set_write_v30_mdl_files (int s) { _write_v30_mdl_files = s;}
    void set_ignore_self_bonds (int s) { _ignore_self_bonds = s;}
    void set_write_mdl_chiral_flags (int s) { _write_mdl_chiral_flags = s;}
    void set_include_chiral_info_in_mdl_outputs (int s) { _include_chiral_info_in_mdl_outputs = s;}
    void set_mdl_read_h_correct_chiral_centres (int s) { _mdl_read_h_correct_chiral_centres = s;}
    void set_mdl_write_h_correct_chiral_centres (int s) { _mdl_write_h_correct_chiral_centres = s;}
    void set_mdl_insert_between_sdf_name_tokens (const const_IWSubstring & s) { _insert_between_sdf_name_tokens = s;}
    int  set_sdf_identifier (const const_IWSubstring & s);
    void set_extract_isis_extregno (int s) { _extract_isis_extregno = s;}
    void set_fetch_all_sdf_identifiers (int s) { _fetch_all_sdf_identifiers = s;}
    void set_take_first_tag_as_name (int s) { _take_first_tag_as_name = s;}
    void set_prepend_sdfid (int s) { _prepend_sdfid = s;}
    void set_discard_sdf_molecule_name (int s) { _discard_sdf_molecule_name = s;}
    void set_multi_record_tag_data_present (int s) { _multi_record_tag_data_present = s;}
    void set_mdl_write_aromatic_bonds (int s) { _mdl_write_aromatic_bonds = s;}
    void set_mdl_write_aromatic_atoms (int s) { _mdl_write_aromatic_atoms = s;}
    void set_display_non_organic_chirality_messages (int s) { _display_non_organic_chirality_messages = s;}
    void set_mdl_display_invalid_chiral_connectivity (int s) { _mdl_display_invalid_chiral_connectivity = s;}
    void set_truncate_long_symbols (int s) { _truncate_long_symbols = s;}
    void set_mdl_change_long_symbols_to (const const_IWSubstring & s);
    void set_discern_chirality_from_wedge_bonds (int s) { _discern_chirality_from_wedge_bonds = s;}
    void set_write_isotopes_as_numbers_rather_than_differences_from_normal (int s) { _write_isotopes_as_numbers_rather_than_differences_from_normal = s;}
    void set_write_M_isotopes_as_numbers_rather_than_differences_from_normal (int s) { _write_M_isotopes_as_numbers_rather_than_differences_from_normal = s;}
    void set_read_isotopes_as_numbers_rather_than_differences_from_normal (int s) { _read_isotopes_as_numbers_rather_than_differences_from_normal = s;}
    void set_read_M_isotopes_as_numbers_rather_than_differences_from_normal (int s) { _read_M_isotopes_as_numbers_rather_than_differences_from_normal = s;}
    void set_replace_first_sdf_tag (const const_IWSubstring & s) { _replace_first_sdf_tag = s;}
    void set_allow_deuterium (int s) { _allow_deuterium = s;}
    void set_allow_tritium (int s) { _allow_tritium = s;}
    int  set_mdl_input_bond_type_translation (int zfrom, int zto);
    int  process_mdl_bond_translation(const const_IWSubstring &);
    void set_mdl_g_records_hold_atom_symbols (int s) { _mdl_g_records_hold_atom_symbols = s;}
    void set_write_mdl_charges_as_m_chg (int s) { _write_mdl_charges_as_m_chg = s;}
    void set_set_elements_based_on_atom_aliases (int s) { _set_elements_based_on_atom_aliases = s;}
    void set_write_Rn_groups_as_element (int s) { _write_Rn_groups_as_element = s;}
    void set_mdl_name_in_m_tag(const const_IWSubstring & s);
    void set_convert_single_atom_sgroup_to_element(const int s) { _convert_single_atom_sgroup_to_element = s;}
    void set_gsub_mdl_file_data(const char c) { _gsub_mdl_file_data = c;}

    int isis_standard_records () const { return _isis_standard_records;}
    int ignore_unrecognised_m_records () const { return _ignore_unrecognised_m_records;}
    int report_unrecognised_records () const { return _report_unrecognised_records;}
    int die_on_erroneous_m_input () const { return _die_on_erroneous_m_input;}
    int write_mdl_dollars () const { return _write_mdl_dollar_signs;}
    int write_mdl_m_end_record () const { return _write_mdl_m_end_record;}
    int write_v30_mdl_files () const { return _write_v30_mdl_files; }
    int ignore_self_bonds () const { return _ignore_self_bonds;}
    int write_mdl_chiral_flags () const { return _write_mdl_chiral_flags;}
    int include_chiral_info_in_mdl_outputs () const { return _include_chiral_info_in_mdl_outputs;}
    int mdl_read_h_correct_chiral_centres () const { return _mdl_read_h_correct_chiral_centres;}
    int mdl_write_h_correct_chiral_centres () const { return _mdl_write_h_correct_chiral_centres;}
    int extract_isis_extregno () const { return _extract_isis_extregno;}
    int fetch_all_sdf_identifiers () const { return _fetch_all_sdf_identifiers;}
    int take_first_tag_as_name () const { return _take_first_tag_as_name;}
    int prepend_sdfid () const { return _prepend_sdfid;}
    int discard_sdf_molecule_name () const { return _discard_sdf_molecule_name;}
    int multi_record_tag_data_present () const { return _multi_record_tag_data_present;}
    int mdl_write_aromatic_bonds () const { return _mdl_write_aromatic_bonds;}
    int mdl_write_aromatic_atoms () const { return _mdl_write_aromatic_atoms;}
    int display_non_organic_chirality_messages () const { return _display_non_organic_chirality_messages;}
    int mdl_display_invalid_chiral_connectivity () const { return _mdl_display_invalid_chiral_connectivity;}
    int truncate_long_symbols () const { return _truncate_long_symbols;}
    int discern_chirality_from_wedge_bonds () const { return _discern_chirality_from_wedge_bonds;}
    int write_isotopes_as_numbers_rather_than_differences_from_normal () const { return _write_isotopes_as_numbers_rather_than_differences_from_normal;}
    int write_M_isotopes_as_numbers_rather_than_differences_from_normal () const { return _write_M_isotopes_as_numbers_rather_than_differences_from_normal;}
    int read_isotopes_as_numbers_rather_than_differences_from_normal () const { return _read_isotopes_as_numbers_rather_than_differences_from_normal;}
    int read_M_isotopes_as_numbers_rather_than_differences_from_normal () const { return _read_M_isotopes_as_numbers_rather_than_differences_from_normal;}
    int allow_deuterium () const { return _allow_deuterium;}
    int allow_tritium () const { return _allow_tritium;}
    int mdl_g_records_hold_atom_symbols () const { return _mdl_g_records_hold_atom_symbols;}
    int write_mdl_charges_as_m_chg () const { return _write_mdl_charges_as_m_chg;}
    int set_elements_based_on_atom_aliases () const { return _set_elements_based_on_atom_aliases;}
    int write_Rn_groups_as_element () const { return _write_Rn_groups_as_element;}
    int convert_single_atom_sgroup_to_element() const { return _convert_single_atom_sgroup_to_element;}
    char gsub_mdl_file_data() const { return _gsub_mdl_file_data;}

    int translate_input_bond_type (int b) const;

    const IWString & insert_between_sdf_name_tokens () const { return _insert_between_sdf_name_tokens;}
    const IWString & name_in_m_tag () const { return _name_in_m_tag;}
    const IWString & replace_first_sdf_tag () { return _replace_first_sdf_tag;}
    int sdf_identifier_matches (const IWString &);
     
    template <typename T> int write_atoms_and_bonds (int, int, T &);

    Atom * create_mdl_atom (const const_IWSubstring & ss, int msdif, int chg, int is_radical) const;

// functions associated with last molecules read. Yes, this is bad...

    void set_accumulate_mdl_chirality_features (int s) { _accumulate_mdl_chirality_features = s;}
    int  accumulate_mdl_chirality_features () const { return _accumulate_mdl_chirality_features;}

    const Set_of_Atoms & mdl_unspecified_chiral_atoms() const { return _unspecified_chiral_atoms_last_molecule;}
    const Set_of_Atoms & down_bonds_last_molecule() const { return _down_bonds_last_molecule;}
    const Set_of_Atoms & up_bonds_last_molecule() const { return _up_bonds_last_molecule;}
    const Set_of_Atoms & squiggle_bonds_last_molecule() const { return _squiggle_bonds_last_molecule;}
    const Set_of_Atoms & unspecified_double_bond_atoms() const { return _unspecified_double_bond_atoms;}

    void add_alias (Atom_Alias * a) { _aliases.add(a); _a_records_found++;}
    int  number_aliases () const { return _aliases.number_elements();}
    const resizable_array_p<Atom_Alias> & atom_aliases () const { return _aliases;}

    void extra_g_record_found () { _g_records_found++;}

    int unspecified_chiral_atom_if_interested (atom_number_t s);

    int parse_bond_record (const_IWSubstring & buffer, int na, atom_number_t & a1, atom_number_t & a2,
                           int & bond_type_read_in, int & directionality);

    void append_isotope_information (IWString & output_buffer, int iso, int normal_isotope) const;

    const IWString & digits2 (int n);
    const IWString & digits3 (int n);

    int add_rdfile_identifier (const IWString & new_identifier);
    const resizable_array_p<IWString> & rdfile_identifiers() const { return _rdfile_identifiers;}

    int set_rdfile_start_of_record (const const_IWSubstring & s);
    const IWString & rdfile_start_of_record () const { return _rdfile_start_of_record;}
};

extern MDL_File_Supporting_Material * global_default_MDL_File_Supporting_Material ();

/*
  A couple of supporting functions used by the various mdl routines
*/

extern int int3d (const const_IWSubstring &, int &, int &, int * = NULL);

#define MDL_RADICAL -998

extern int convert_from_mdl_charge (int);

extern void add_to_text_info (resizable_array_p<IWString> & text_info, const const_IWSubstring & zextra);

extern int extract_sdf_identifier(const IWString & buffer, IWString & id);

extern int convert_from_mdl_number_to_bond_type (int int_rep, bond_type_t & bt);

/*
  When reading the M lines in MDL files, we have pairs of atom numbers and atom properties.
  This class describes such a grouping.
*/

struct Aprop
{
  int _atom_number;
  int _property;
};

typedef struct Aprop Aprop;

extern int write_v30_record (IWString & buffer, std::ostream & output);

/*
  According to the documentation, there can be a max of MAX_PAIRS of these pairs on a record

  March 2014. Some MDL writers will put arbitrary numbers of items on a record. Use a large
  number (never guaranteed to be large enough)
*/

#define MAX_PAIRS 100

extern int fill_atom_property_array (const IWString & buffer, int &, Aprop * atom_properties);

extern int parse_bond_record (const_IWSubstring & buffer,
                   int na,
                   atom_number_t & a1, atom_number_t & a2,
                   int & bond_type_read_in,
                   int & directionality);

class Atom;

#endif
