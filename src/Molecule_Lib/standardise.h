#ifndef MOLECULE_LIB_STANDARDISE_H_
#define MOLECULE_LIB_STANDARDISE_H_

// Chemical standardisation.
// Transformations that change one molecular representation to another.

#include <iostream>
#include <string>

class Molecule;
class Atom;
class Command_Line;

#include "molecule.h"
#include "set_of_atoms.h"
#include "substructure.h"
#include "Molecule_Lib/standardise.pb.h"

namespace standardisation {
// Transformations can be externally specified.
class ExternalTransformation {
  private:
    uint64_t _molecules_processed;
    uint64_t _molecules_changed;

    // The proto from which we are built.
    standardisation::Standardisation _standardisation;

    // The query built either from the smarts or the query file
    // in the proto.
    Substructure_Query _query;

    // This will be constructed from the `smiles` attribute of
    // the proto. It stores the bonds to be made. This all works
    // because the atom numbers in `_molecule` correspond to the
    // query atom numbers in `_query`.
    Molecule _molecule;

  // Private functions.

    int Process(Molecule& m, const Set_of_Atoms& embedding);

  public:
    ExternalTransformation();

    // This will copy `proto` to _standardisation.
    // `fname` is needed if there is a query file. The default
    // assumption is that the query file will be in the same
    // directory as the proto.
    int Build(const standardisation::Standardisation& proto,
                              const IWString& fname);

    // Returns the name() field in the proto.
    const std::string& name() const;

    uint64_t molecules_processed() const {
      return _molecules_processed;
    }
    uint64_t molecules_changed() const {
      return _molecules_changed;
    }

    int Process(Molecule& m);
};

}  // namespace standardisation

class Chemical_Transformation
{
  private:
    int _active;
    uint64_t _groups_changed;
    uint64_t _molecules_changed;

  public:
    Chemical_Transformation ();
    ~Chemical_Transformation ();

    int active () const { return _active;}
    void activate () { _active = 1;}
    void deactivate () { _active = 0;}

    uint64_t molecules_changed () const { return _molecules_changed;}
    uint64_t groups_changed ()    const { return _groups_changed;}

    void extra (int);

    int report (std::ostream & ) const;
};

class Possible_Lactim_Lactam
{
  private:
    const atom_number_t _oxygen;
    const atom_number_t _carbon;
    atom_number_t _nitrogen;

    atom_number_t _second_nitrogen;    // two connected, could be part of another lactam/lactim

    int _total_nitrogen_attachments;      // might be 2, but one might be 3 connected

    int _shared_nitrogen_group;

    atom_number_t _alpha_nitrogen;       // bonded to either _nitrogen or _second_nitrogen

    int _is_ring;
    int _fused_system_size;
    int _fused_system_identifier;
    int _lactims_in_fused_system;
    int _aromatic;
    int _ring_size;
    int _lactam_form;
    int _is_urea;

    int _made_change;    // keep track of result of to_lactim_form() call

//  private functions
    
    int _alpha_nitrogen_present (const Molecule & m, const atom_number_t n) const;

  public:
    Possible_Lactim_Lactam (atom_number_t o, atom_number_t c, atom_number_t n);

    atom_number_t oxygen () const { return _oxygen;}
    atom_number_t carbon () const { return _carbon;}
    atom_number_t nitrogen () const { return _nitrogen;}

    void set_is_ring (int s) { _is_ring = s;}
    void set_fused_system_size (int s) { _fused_system_size = s;}
    void set_fused_system_identifier (int s) { _fused_system_identifier = s;}
    void set_aromatic (int s) { _aromatic = s;}
    void set_ring_size (int s) { _ring_size = s;}
    void set_lactam_form (int s) { _lactam_form = s;}
    void set_second_nitrogen (atom_number_t s) { _second_nitrogen = s;}
    void set_total_nitrogen_attachments (atom_number_t s) { _total_nitrogen_attachments = s;}
    void set_shared_nitrogen_group (int s) { _shared_nitrogen_group = s;}
    void set_is_urea (int s) { _is_urea = s;}

    int is_ring () const { return _is_ring;}
    int fused_system_size () const { return _fused_system_size;}
    int fused_system_identifier () const { return _fused_system_identifier;}
    int ring_size () const { return _ring_size;}
    int aromatic () const { return _aromatic;}
    int lactam_form () const { return _lactam_form;}
    atom_number_t second_nitrogen () const { return _second_nitrogen;}
    atom_number_t total_nitrogen_attachments () const { return _total_nitrogen_attachments;}
    int shared_nitrogen_group () const { return _shared_nitrogen_group;}
    int is_urea () const { return _is_urea;}

    int shares_nitrogen_with (const Possible_Lactim_Lactam & rhs) const;

    int discern_alpha_nitrogen (Molecule & m);
    atom_number_t alpha_nitrogen() const { return _alpha_nitrogen;}

    void increment_lactims_in_fused_system () { _lactims_in_fused_system++;}
    int lactims_in_fused_system () const { return _lactims_in_fused_system;}

    int to_lactim_form (Molecule & m, int check_nitrogen_hcount);
    int to_lactam_form (Molecule & m, int check_nitrogen_hcount);
    int made_change () const { return _made_change;}

    int hcount (Molecule &) const;
    int could_change_to_lactim_with_current_bonding (const Molecule & m) const;

    // Returns non zero if `m` has _nitrogen=_carbon=_oxygen
    int TwoDoubleBonds(const Molecule& m) const;

    int add_unique_nitrogens (Set_of_Atoms & unique_nitrogens) const;
    int reperceive (Molecule & m);
};

/*
  for thread safety, the Chemical_Standardisation object needs a stack
  based object in which it can hold attributes of the current molecule 
  being processed
*/

class IWStandard_Current_Molecule
{
  private:
    int _matoms;
    int _nrings;
    atomic_number_t * _atomic_number;
    int * _ncon;
    int * _ring_membership;
    int * _ring_size;
    int * _ring_is_fused;
    int * _atom_is_aromatic;
    int * _ring_is_aromatic;
    int * _ring_bond_count;
    int * _fsid;
    const Atom ** _atom;

    int * _ring_nitrogen_count;

    int _sulphur;
    int _nneg;    // number of negative charges
    int _npos;    // number of positive charges
    int _ominus;  // O-
    int _sminus;  // S-
    int _nplus;   // N+
    int _cminus;  // C-
    int _splus;   // S+

//  needed for reverse transformations

    int _nitrogens;
    int _oxygens;
    int _phosphorus;

    int _non_organic;
    int _isolated_metal;
    int _isolated_halogen;
    int _singly_connected_metal;
    Set_of_Atoms _possible_guanidine;
    int _explicit_hydrogen_count;

    int _possible_valence_errors;

    int _isotope;

    resizable_array_p<Possible_Lactim_Lactam> _possible_lactam;

// We save the rings because they may get recomputed during our changes

    resizable_array_p<Set_of_Atoms> _rings;

//  private functions

  public:
    IWStandard_Current_Molecule();
    ~IWStandard_Current_Molecule();

    int initialise (Molecule &);

    int processing_needed () const;

    const atomic_number_t * atomic_number () const { return _atomic_number;}
    const Atom * const * atoms () { return _atom;}
    const int *  ncon () const { return _ncon;}
    const int *  ring_membership () const { return _ring_membership;}
    const int *  ring_size () const { return _ring_size;}
    const int *  ring_is_fused () const { return _ring_is_fused;}
    const int *  ring_is_aromatic () const { return _ring_is_aromatic;}
    const int *  atom_is_aromatic () const { return _atom_is_aromatic;}
    const int *  ring_nitrogen_count () const { return _ring_nitrogen_count;}
    const int *  ring_bond_count() const { return _ring_bond_count;}

//  Some methods are non const

    int *  ncon () { return _ncon;}

    int nrings () const { return _nrings;}

    int aromatic_rings_with_multiple_nitrogens () const;

    void change_sulphur (int s) { _sulphur += s;}
    void change_nneg (int s) { _nneg += s;}
    void change_npos (int s) { _npos += s;}
    void change_ominus (int s) { _ominus += s;}
    void change_sminus (int s) { _sminus += s;}
    void change_splus (int s) { _splus += s;}
    void change_nplus (int s) { _nplus += s;}
    void change_cminus (int s) { _cminus += s;}
    void change_nitrogens (int s) { _nitrogens += s;}
    void change_oxygens (int s) { _oxygens += s;}
    void change_isolated_metal (int s) { _isolated_metal += s;}
    void change_isolated_halogen (int s) { _isolated_halogen += s;}
    void change_singly_connected_metal (int s) { _singly_connected_metal += s;}
//  void change_possible_guanidine (int s) { _possible_guanidine += s;}
    void change_phosphorus (int s) { _phosphorus += s;}
    void change_explicit_hydrogen_count (int s) { _explicit_hydrogen_count += s;}
    void change_possible_valence_errors (int s) { _possible_valence_errors += s;}

    void set_sulphur (int s) { _sulphur = s;}
    void set_nneg (int s) { _nneg = s;}
    void set_npos (int s) { _npos = s;}
    void set_ominus (int s) { _ominus = s;}
    void set_sminus (int s) { _sminus = s;}
    void set_nplus (int s) { _nplus = s;}
    void set_cminus (int s) { _cminus = s;}
    void set_nitrogens (int s) { _nitrogens = s;}
    void set_oxygens (int s) { _oxygens = s;}
    void set_isolated_metal (int s) { _isolated_metal = s;}
    void set_isolated_halogen (int s) { _isolated_halogen = s;}
    void set_singly_connected_metal (int s) { _singly_connected_metal = s;}
//  void set_possible_guanidine (int s) { _possible_guanidine = s;}
    void set_phosphorus (int s) { _phosphorus = s;}
    void set_explicit_hydrogen_count (int s) { _explicit_hydrogen_count = s;}
    void set_possible_valence_errors (int s) { _possible_valence_errors = s;}

    int sulphur () const { return  _sulphur;}
    int nneg () const { return  _nneg;}
    int npos () const { return  _npos;}
    int ominus () const { return  _ominus;}
    int sminus () const { return  _sminus;}
    int splus () const { return  _splus;}
    int nplus () const { return  _nplus;}
    int cminus () const { return  _cminus;}
    int nitrogens () const { return  _nitrogens;}
    int oxygens () const { return  _oxygens;}
    int isolated_metal () const { return  _isolated_metal;}
    int isolated_halogen () const { return  _isolated_halogen;}
    int singly_connected_metal () const { return  _singly_connected_metal;}
    int non_organic() const { return _non_organic;}
    const Set_of_Atoms & possible_guanidine () const { return  _possible_guanidine;}
    int phosphorus () const { return  _phosphorus;}
    int explicit_hydrogen_count () const { return  _explicit_hydrogen_count;}
    int possible_valence_errors () const { return  _possible_valence_errors;}
    const resizable_array_p<Possible_Lactim_Lactam> & possible_lactam() const { return _possible_lactam;}
    int isotope() const {
      return _isotope;
    }

    int remove_possible_guanidine (const atom_number_t);

    int ring_number_containing_atom (const atom_number_t s) const;
    const Set_of_Atoms * ring_containing_atom (const atom_number_t) const;
    const Set_of_Atoms * ringi (const int) const;

    int fused_system_identifier (const atom_number_t s) const { return _fsid[s];}
    const resizable_array_p<Set_of_Atoms>& rings() const {
      return _rings;
    }
};

// Names of chemical standardisations. Pass to Activate() to turn on
// individual transformations. Should transition to an enum.
#define CS_NITRO "nitro"
#define CS_NpOm  "n+o-"
#define CS_NpNm  "n+n-"
#define CS_SpCm  "s+c-"
#define CS_ALLpm "all+-"
#define CS_XH    "xh"
#define CS_NpH3  "n+h3"
#define CS_AMINE "amine"
#define CS_Om    "o-"
#define CS_Nm    "n-"
#define CS_Cm    "c-"
#define CS_ALL   "all"
#define CS_NRMCH "nrmch"
#define CS_COVM  "covm"
#define CS_ISOLC "isolc"
#define CS_GUAND "guan"
#define CS_GUANDR "Rguan"
#define CS_SPOM  "s+o-"
#define CS_ACID  "acid"
#define CS_EHLST "ehlast"
#define CS_FMRK  "fmrk"
#define CS_AZID  "azid"
#define CS_MSDUR "msdur"
#define CS_MSDSA "msdsa"
#define CS_FCRN  "fcor"
//#define CS_RNPNM "Rn+n-"
#define CS_FWIH  "fwih"
#define CS_IMIDAZOLE  "imidazole"
#define CS_CHARGED_IMIDAZOLE  "charged_imidazole"
#define CS_IMIDAZOLE_EXONH "imidazoleNH"
#define CS_PYRAZOLE  "pyrazole"
#define CS_TRIAZOLE  "triazole"
#define CS_TETRAZOLE  "tetrazole"
#define CS_PYRIMIDINE "pyrimidine"
#define CS_LACTIM_LACTAM "ltlt"
#define CS_LACTIM_LACTAM_RING "ltltr"
#define CS_REVERSE_NITRO "rvnitro"
#define CS_REVERSE_NV5 "rvnv5"
#define CS_ISOXAZOLE  "isoxazole"
#define CS_ARGUAN "arguan"
#define CS_PYRAZOLONE "pyrazolone"
#define CS_AMINO_THIAZOLE "aminothazole"
#define CS_KETO_ENOL "keto_enol"
#define CS_4_PYRIDONE "4-pyridone"
#define CS_SULFONYL_UREA "surea"
#define CS_124TRIAZINE "124-triazine"
#define CS_ENOL_FUSED "enol-fused"
#define CS_FCNO "fcno"
#define CS_ISOTOPE "isotope"


namespace standardise {
// Mar 2022. Some transformations can result in different molecules depending on
// the ordering of the connection table. Optionally switch to unique smiles order
// to stop that happening.
enum class Canonicalise {
  kNone,
  kReorderAtoms,
  kReinterpretSmiles
};
}  // namespace standardise

/*
  As you add standardisations, make sure you update the code around the "all" directive
*/

class Chemical_Standardisation
{
  private:
    int _ok;
    int _verbose;
    int _active;

    uint64_t _molecules_processed;
    uint64_t _molecules_changed;

//  Our default set of conventions

    Chemical_Transformation _transform_amines;
    Chemical_Transformation _transform_nitro;
    Chemical_Transformation _transform_nplus_ominus;
    Chemical_Transformation _transform_plus_minus;
    Chemical_Transformation _transform_n_charge_sep;
    Chemical_Transformation _protonate_no;
    Chemical_Transformation _remove_hydrogens;
    Chemical_Transformation _protonate_carboxyllic_acids;
    Chemical_Transformation _protonate_sulfonic_acids;
    Chemical_Transformation _protonate_sulfinic_acids;
    Chemical_Transformation _transform_splus_cminus;
    Chemical_Transformation _transform_cminus;
    Chemical_Transformation _transform_ominus;
    Chemical_Transformation _transform_nminus;
    Chemical_Transformation _transform_covalent_metals;
    Chemical_Transformation _transform_single_atom_ions;
    Chemical_Transformation _transform_guanidine;
    Chemical_Transformation _transform_guanidine_ring;
    Chemical_Transformation _protonate_sulfur_acids;
    Chemical_Transformation _protonate_phosphorous_acids;
    Chemical_Transformation _from_mrk_standardisations;
    Chemical_Transformation _explicit_hydrogens_last;
    Chemical_Transformation _transform_tetrazole;
    Chemical_Transformation _transform_azid;
    Chemical_Transformation _transform_isoxazole;
    Chemical_Transformation _transform_misdrawn_sulfonamide;
    Chemical_Transformation _transform_misdrawn_urea;
    Chemical_Transformation _transform_imidazole;
    Chemical_Transformation _transform_charged_imidazole;
    Chemical_Transformation _transform_imidazole_exocyclic_nh;
    Chemical_Transformation _transform_pyrazole;
    Chemical_Transformation _transform_triazole;
    Chemical_Transformation _transform_pyrimidine;
    Chemical_Transformation _transform_lactim_lactam;
    Chemical_Transformation _transform_lactim_lactam_ring;
    Chemical_Transformation _transform_aromatic_guanidine_ring;
    Chemical_Transformation _transform_pyrazolone;
    Chemical_Transformation _transform_amino_thiazole;
    Chemical_Transformation _transform_enol_to_keto;
    Chemical_Transformation _transform_to_4_pyridone;
    Chemical_Transformation _transform_sulfonyl_urea;
    Chemical_Transformation _transform_124_triazine;
    Chemical_Transformation _transform_enol_fused;
    Chemical_Transformation _transform_charged_non_organic;
    Chemical_Transformation _transform_isotopes;

//  Various reverse direction transformations

    Chemical_Transformation _transform_nitro_reverse;
    Chemical_Transformation _transform_back_to_nplus_nminus;
    Chemical_Transformation _transform_nv5_to_charge_separated;
    Chemical_Transformation _transform_to_charge_separated_azid;

    Chemical_Transformation _transform_obvious_implicit_hydrogen_errors;

    int _remove_hydrogens_attached_to_chiral_centres;

//  Sometimes it is useful to flag changes. Either a single text appended to
//  all changed molecules

    IWString _append_to_changed_molecules;

//  or a string that tells exactly what was done

    int _append_string_depending_on_what_changed;

    int _check_valence_before_and_after;

    standardise::Canonicalise _convert_to_canonical_order;

    // Any number of externally specified transformations.
    resizable_array_p<standardisation::ExternalTransformation> _external_transformation;

//  Some other possibilities


// private functions

    void _default_values ();

    int ConvertToCanonicalOrder(Molecule& m);

    void _do_transform_plus_minus_pair (Molecule & m, atom_number_t a1, atom_number_t a2, IWStandard_Current_Molecule & current_molecule_data);

    int  _do_unset_isotopes(Molecule& m);
    int  _do_transform_amines   (Molecule &, Set_of_Atoms &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_nitro    (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_nplus_ominus (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_nv5_to_charge_separated(Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int _do_nv5_to_charge_separated(Molecule& m, atom_number_t zatom, const IWStandard_Current_Molecule& current_molecule_data);
    int  _do_transform_plus_minus   (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_n_charge_sep (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_azid_to_charge_separated (Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_protonate_no (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_remove_hydrogens (Molecule &);
    int  _do_protonate_carboxyllic_acids (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_protonate_sulfonic_acids (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_protonate_sulfinic_acids (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_splus_cminus (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_cminus       (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_ominus       (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_nminus       (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_covalent_metals (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_covalent_metals (Molecule & m);
    int  _do_transform_guanidine (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_ring_guanidine (Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_aromatic_ring_guanidine (Molecule & m, const Set_of_Atoms & r, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_aromatic_ring_guanidine (Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_single_atom_ions (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_protonate_sulfur_acids (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_protonate_phosphorous_acids (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_from_mrk_standardisations (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_tetrazole (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_tetrazole (Molecule &, const Set_of_Atoms &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_charged_imidazole(Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_charged_imidazole(Molecule &, const int ring_number, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_imidazole (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_imidazole (Molecule &, const int, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_imidazole_exocyclic_nh(Molecule& m, IWStandard_Current_Molecule& current_molecule_data);
    int  _do_imidazole_exocyclic_nh(Molecule & m,
                        IWStandard_Current_Molecule& current_molecule_data,
                        atom_number_t nh, atom_number_t carbon);
    int  _swap_imidazole (Molecule & m, atom_number_t n1, atom_number_t c, atom_number_t n2) const;
    int ResolveImidazole(Molecule& m,
                const Set_of_Atoms& r,
                atom_number_t c1,
                atom_number_t nh0,
                atom_number_t c2,
                atom_number_t c3,
                atom_number_t nh1);
    int  _do_pyrazole (Molecule &, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_pyrazole (Molecule &, int * atom_already_changed, const int, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_pyrimidine(Molecule& m, IWStandard_Current_Molecule& current_molecule_data);
    int  _do_pyrimidine(Molecule& m, IWStandard_Current_Molecule& current_molecule_data,
                const Set_of_Atoms& ring);
    int  _swap_charged_pyrazole(Molecule& m, const int ring_number, IWStandard_Current_Molecule & current_molecule_data, const atom_number_t n1, const atom_number_t n2);
    int  _do_triazole (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_triazole (Molecule &, const Set_of_Atoms &, const Set_of_Atoms * is_fused, IWStandard_Current_Molecule & current_molecule_data);
//  int  _do_123_triazole (Molecule & m, const Ring & r, int n1_index_in_ring, int n2_index_in_ring, int n3_index_in_ring, const atomic_number_t * z, const int * ncon, Atom ** atoms);
    int  _do_123_triazole (Molecule & m, const Set_of_Atoms & r, int n1_index_in_ring, int n2_index_in_ring, int n3_index_in_ring, int c4_index_in_ring, int c5_index_in_ring, IWStandard_Current_Molecule &, const Set_of_Atoms * is_fused) const;
    int  _do_134_triazole (Molecule & m, const Set_of_Atoms & r, int n1_index_in_ring, int c2_index_in_ring, int n3_index_in_ring, int n4_index_in_ring, int c5_index_in_ring, IWStandard_Current_Molecule &) const;
//  int  _do_134_triazole (Molecule & m, const Ring & r, int n1_index_in_ring, int n2_index_in_ring, int n3_index_in_ring, const atomic_number_t * z, const int * ncon, Atom ** atoms);
    int  _do_isoxazole (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_isoxazole (Molecule &, const Set_of_Atoms &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_azid  (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_misdrawn_urea (Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_misdrawn_sulfonamide (Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_misdrawn_sulfonamide (Molecule & m, const atom_number_t s, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_back_to_nplus_nminus  (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_explicit_hydrogens_last (Molecule &);
    int  _do_amino_thiazole (Molecule &, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_amino_thiazole (Molecule & m, const Set_of_Atoms & r, const int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_enol_to_keto(Molecule & m, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_diketo_to_enol(Molecule & m,
                        IWStandard_Current_Molecule& current_molecule_data,
                        int * diketone13);
    int _identify_current_keto_and_enol_forms(Molecule & m, IWStandard_Current_Molecule& current_molecule_data,
                        Set_of_Atoms& current_keto_form,
                        Set_of_Atoms& current_enol_form) const;

    int _lactim_lactam_process_if_all_groups_non_overlapping (Molecule & m, const resizable_array<int> & in_system,
                                int * already_done, IWStandard_Current_Molecule & current_molecule_data);
    int _do_lactam_lactim_pyrazole_triazole (Molecule & m, const Possible_Lactim_Lactam & p, const Set_of_Atoms & r, IWStandard_Current_Molecule & current_molecule_data);
    int _do_lactam_lactim_pyrazole_triazole (Molecule & m, const resizable_array_p<Set_of_Atoms> &, int * atom_already_changed, int * already_done, IWStandard_Current_Molecule & current_molecule_data);
    int _do_transform_ring_lactim(Molecule &, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int __do_transform_ring_lactim (Molecule & m, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int _do_transform_non_ring_lactim(Molecule &, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int _toggle_kekule_forms_to_lactim_form (Molecule & m, IWStandard_Current_Molecule & current_molecule_data) const;
    int _change_molecule_kekule_form_for_lactam_canonical (Molecule & m, Possible_Lactim_Lactam & p);
    int _handle_one_lactim_lactam_in_ring_system (Molecule & m, Possible_Lactim_Lactam & p,
                                                IWStandard_Current_Molecule & current_molecule_data) const;
    int _handle_two_nitrogens_sharing_an_oxygen (Molecule & m, Possible_Lactim_Lactam & p1,
                                                        Possible_Lactim_Lactam & p2,
                                                        IWStandard_Current_Molecule & current_molecule_data) const;
    int _identify_and_handle_two_nitrogens_sharing_an_oxygen (Molecule & m, resizable_array<Possible_Lactim_Lactam *> & s,
                                        IWStandard_Current_Molecule & current_molecule_data) const;
    int _process_lactim_in_isolated_aromatic_ring (Molecule & m, Possible_Lactim_Lactam & p,
                                                IWStandard_Current_Molecule & current_molecule_data);

    int  _do_transform_reverse_nitro (Molecule &, IWStandard_Current_Molecule & current_molecule_data);
//    int  _do_transform_reverse_azid  (Molecule &, IWStandard_Current_Molecule & current_molecule_data);

    int  _do_transform_pyrazolone (Molecule & m, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);
    int  _do_transform_pyrazolone (Molecule & m, const Set_of_Atoms & ri, const int ring_is_fused, int * atom_already_changed, IWStandard_Current_Molecule & current_molecule_data);

    int _do_transform_to_4_pyridone(Molecule& m, IWStandard_Current_Molecule& current_molecule_data);
    int _do_transform_to_4_pyridone(Molecule& m,
                        IWStandard_Current_Molecule& current_molecule_data,
                        const Set_of_Atoms& ring) const;
    int _do_transform_to_4_pyridone(Molecule& m, const Set_of_Atoms& r, int n_index, int oh_index, atom_number_t oh) const;

    int _do_transform_sulfonyl_urea(Molecule& m, IWStandard_Current_Molecule& current_molecule_data);
    int _do_transform_124_triazine(Molecule& m, IWStandard_Current_Molecule& current_molecule_data);
    int _do_transform_124_triazine(Molecule& m,
                        const Set_of_Atoms& ring,
                        IWStandard_Current_Molecule& current_molecule_data);
    int _do_transform_124_triazine(Molecule& m, atom_number_t o1, atom_number_t c1,
                atom_number_t n1,
                atom_number_t c2, atom_number_t o2,
                atom_number_t n2, atom_number_t n3, atom_number_t c3) const;

    int _do_transform_enol_fused(Molecule& m, IWStandard_Current_Molecule& current_molecule_data);

    int _do_transform_charged_non_organics(Molecule& m,
                IWStandard_Current_Molecule& current_molecule_data);
    int  _do_transform_implicit_hydrogen_known_errors (Molecule & m, IWStandard_Current_Molecule & current_molecule_data);

    int _activate_nohmove_transformations();

    int  _process (Molecule &);
    int  _process (Molecule &, IWStandard_Current_Molecule & current_molecule_data);

    int _processing_needed (const IWStandard_Current_Molecule & current_molecule_data) const;

    int AddExternalSpecification(const const_IWSubstring& directive);
    int ReadFileOfExternalProtos(iwstring_data_source& input, const const_IWSubstring& fname);
    int ReadFileOfExternalProtos(const const_IWSubstring& fname);
    int DoExternalTransformations(Molecule& m);

  public:
    Chemical_Standardisation ();
    ~Chemical_Standardisation ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int active () const { return _active;}

    void set_verbose (int v) { _verbose = v;}

    void set_convert_to_canonical_order(standardise::Canonicalise s) {
      _convert_to_canonical_order = s;
    }
    void set_append_string_depending_on_what_changed(int s);

    int construct_from_command_line (Command_Line &, int = 0, char = 'g');

    // Turn on individual chemical transformations based on the name.
    int Activate(const IWString& directive, const int verbose);

    int process (Molecule &);

    int report (std::ostream &) const;

    void activate_all ();

    int activate_all_except_hydrogen_removal ();

    void deactivate_lactim_lactam ();

    int activate_from_corina_transformations();

    void deactivate () { _active = 0;}

    int deactivate (const const_IWSubstring &);

    // `buffer` must contain a TextProto representation of a substitution::Substitution
    // proto. `fname` is needed in the case of there being a query_file specification
    // in the proto.
    int ReadExternalProto(const const_IWSubstring& buffer, const const_IWSubstring& fname);
};

extern int display_standard_chemical_standardisation_options (std::ostream &, char);

extern void set_update_chemical_standardisation_accumulators (int s);

#endif  // MOLECULE_LIB_STANDARDISE_H_
