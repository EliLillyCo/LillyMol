/*
  Chemical standardisation routines.
*/

#include <memory>
#include <algorithm>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "cmdline.h"
#include "iwqsort.h"
#include "misc.h"

#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "path.h"
#include "toggle_kekule_form.h"
#include "smiles.h"
#include "misc2.h"

class Molecule_Data_for_Standardisation
{
  private:
    int _natoms;
    int _nrings;
    int * _ncon;
    int * _atom_is_aromatic;
    int * _ring_is_aromatic;
    atomic_number_t _atomic_number;

  public:
    Molecule_Data_for_Standardisation();
    ~Molecule_Data_for_Standardisation();

    int initialise (Molecule & m);
};

Possible_Lactim_Lactam::Possible_Lactim_Lactam (atom_number_t o,
                                                atom_number_t c,
                                                atom_number_t n) : _oxygen(o), _carbon(c), _nitrogen(n)
{
  _second_nitrogen = INVALID_ATOM_NUMBER;
  _total_nitrogen_attachments = 0;
  _shared_nitrogen_group = 0;
  _alpha_nitrogen = INVALID_ATOM_NUMBER;

  _is_ring = 0;
  _fused_system_size = 0;
  _fused_system_identifier = 0;
  _lactims_in_fused_system = 1;          // by default, we are in the system (will be ignored for non ring 
  _ring_size = 0;
  _aromatic = 0;
  _lactam_form = 0;
  _is_urea = 0;
  _made_change = 0;

  return;
}

int
Possible_Lactim_Lactam::shares_nitrogen_with (const Possible_Lactim_Lactam & rhs) const
{
  if (_nitrogen == rhs._nitrogen)
    return 1;

  if (_nitrogen == rhs._second_nitrogen)
    return 1;

  if (INVALID_ATOM_NUMBER == _second_nitrogen)
    return 0;

  if (_second_nitrogen == rhs._nitrogen)
    return 1;

  if (_second_nitrogen == rhs._second_nitrogen)
    return 1;

  return 0;
}

int
Possible_Lactim_Lactam::to_lactim_form (Molecule & m, int check_nitrogen_hcount)
{
  if (! m.bond_between_atoms(_oxygen, _carbon)->is_double_bond())
    return 0;

  if (! m.bond_between_atoms(_carbon, _nitrogen)->is_single_bond())
    return 0;

  if (! check_nitrogen_hcount)
    ;
  else if (0 == m.hcount(_nitrogen))
    return 0;

  m.set_bond_type_between_atoms(_oxygen, _carbon, SINGLE_BOND);
  m.set_bond_type_between_atoms(_carbon, _nitrogen, DOUBLE_BOND);

  _made_change = 1;

  return 1;
}

int
Possible_Lactim_Lactam::to_lactam_form (Molecule & m, int check_nitrogen_hcount)
{
  if (! m.bond_between_atoms(_oxygen, _carbon)->is_single_bond())
    return 0;

  if (! m.bond_between_atoms(_carbon, _nitrogen)->is_double_bond())
    return 0;

  if (! check_nitrogen_hcount)
    ;
  else if (m.hcount(_nitrogen))
    return 0;

  if (-1 == m.formal_charge(_oxygen))
    m.set_formal_charge(_oxygen, 0);

  m.set_bond_type_between_atoms(_oxygen, _carbon, DOUBLE_BOND);
  m.set_bond_type_between_atoms(_carbon, _nitrogen, SINGLE_BOND);
  m.set_implicit_hydrogens_known(_nitrogen, 0);

  _made_change = 1;

  return 1;
}

int
Possible_Lactim_Lactam::hcount (Molecule & m) const
{
  int rc = m.hcount(_nitrogen);

  if (INVALID_ATOM_NUMBER != _second_nitrogen)
    rc += m.hcount(_second_nitrogen);

  return rc;
}

int
Possible_Lactim_Lactam::could_change_to_lactim_with_current_bonding (const Molecule & m) const
{
  if (! m.bond_between_atoms(_oxygen, _carbon)->is_single_bond())
    return 0;

  if (! m.bond_between_atoms(_carbon, _nitrogen)->is_double_bond())
    return 0;

  return 1;
}

int
Possible_Lactim_Lactam::add_unique_nitrogens (Set_of_Atoms & unique_nitrogens) const
{
  int rc = 0;

  if (unique_nitrogens.add_if_not_already_present(_nitrogen))
    rc = 1;

  if (INVALID_ATOM_NUMBER == _second_nitrogen)
    return rc;

  if (unique_nitrogens.add_if_not_already_present(_second_nitrogen))
    rc++;

  return rc;
}

/*
  Someone may have switched a Kekule form
*/

int
Possible_Lactim_Lactam::reperceive (Molecule & m)
{
  if (INVALID_ATOM_NUMBER == _second_nitrogen)    // no change possible
    return 0;

  const Bond * b1 = m.bond_between_atoms(_carbon, _nitrogen);

  if (b1->is_double_bond())   // we are good
    return 1;

  const Bond * b2 = m.bond_between_atoms(_carbon, _second_nitrogen);

  if (b2->is_single_bond())   // we are good
    return 1;

//cerr << "reperceive, bond1 " << b1->is_double_bond() << " bond2 " << b2->is_single_bond() << endl;

  if (b1->is_single_bond() && b2->is_double_bond())
  {
    std::swap(_nitrogen, _second_nitrogen);
    return 1;
  }

  return 0;
}

/*
  In a multi-threaded environment, we cannot update accumulators
*/

static int update_accumulators = 1;

void
set_update_chemical_standardisation_accumulators (int s)
{
  update_accumulators = s;
}

Chemical_Transformation::Chemical_Transformation()
{
  _active = 0;
  _groups_changed = 0;
  _molecules_changed = 0;

  return;
}

Chemical_Transformation::~Chemical_Transformation()
{
  _active = 0;
  _groups_changed = -1;

  return;
}

/*
  A molecule has been processed, and N groups were changed in it.
*/

void
Chemical_Transformation::extra(int n)
{
  if (! update_accumulators)
    ;
  else if (n)
  {
    _groups_changed += n;
    _molecules_changed++;
  }

  return;
}

int
Chemical_Transformation::report(std::ostream & os) const
{
  os << "changed " << _groups_changed << " groups";
  if (_groups_changed)
    os << " in " << _molecules_changed << " molecules";

  os << endl;

  return os.good();
}

void
Chemical_Standardisation::_default_values()
{
  _ok = 1;
  _active = 0;

  _append_string_depending_on_what_changed = 0;

  _remove_hydrogens_attached_to_chiral_centres = 1;

  _check_valence_before_and_after = 0;
  
  _verbose = 0;

  return;
}

Chemical_Standardisation::Chemical_Standardisation()
{
  _default_values();

  return;
}

Chemical_Standardisation::~Chemical_Standardisation()
{
  _ok = 0;

  return;
}

int
Chemical_Standardisation::ok() const
{
  return _ok;
}

int
Chemical_Standardisation::debug_print (std::ostream & os) const
{
  os << "Chemical Standardisation object\n";

  return os.good();
}

int
display_standard_chemical_standardisation_options (std::ostream & os, char zoption)
{
  os << "  -" << zoption << " <qualifier> chemical standardisations, enter \"-" << zoption << " help\" for usage\n";

  return os.good();
}

void
Chemical_Standardisation::activate_all()
{
  _transform_amines.activate();
  _transform_nplus_ominus.activate();    // no need to do nitro, as they are a subset
  _transform_n_charge_sep.activate();
  _remove_hydrogens.activate();
  _protonate_carboxyllic_acids.activate();
  _protonate_no.activate();
  _protonate_sulfonic_acids.activate();
  _protonate_sulfinic_acids.activate();
  _transform_splus_cminus.activate();
  _transform_ominus.activate();
  _transform_nminus.activate();
  _transform_covalent_metals.activate();
  _transform_guanidine.activate();
  _transform_guanidine_ring.activate();
  _transform_tetrazole.activate();
  _transform_azid.activate();
  _transform_misdrawn_urea.activate();
  _transform_imidazole.activate();
  _transform_pyrazole.activate();
  _transform_lactim_lactam.activate();
  _transform_lactim_lactam_ring.activate();
  _transform_triazole.activate();
  _transform_isoxazole.activate();
  _transform_aromatic_guanidine_ring.activate();
  _transform_pyrazolone.activate();
  _transform_amino_thiazole.activate();

  _active = 1;

  return;
}

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
#define CS_FCRN  "fcor"
#define CS_RNPNM "Rn+n-"
#define CS_FWIH  "fwih"
#define CS_IMIDAZOLE  "imidazole"
#define CS_PYRAZOLE  "pyrazole"
#define CS_TRIAZOLE  "triazole"
#define CS_TETRAZOLE  "tetrazole"
#define CS_LACTIM_LACTAM "ltlt"
#define CS_LACTIM_LACTAM_RING "ltltr"
#define CS_REVERSE_NITRO "rvnitro"
#define CS_REVERSE_NV5 "rvnv5"
#define CS_ISOXAZOLE  "isoxazole"
#define CS_ARGUAN "arguan"
#define CS_PYRAZOLONE "pirazolone"
#define CS_AMINO_THIAZOLE "aminothazole"

int
display_all_chemical_standardisation_options (std::ostream & os, char zoption)
{
  os << "  -" << zoption << ' ' << CS_NITRO << "       transform nitro groups to N(=O)=O\n";
  os << "  -" << zoption << ' ' << CS_NpOm << "        transform charge separated [N+]-[O-] (includes nitro)\n";
  os << "  -" << zoption << ' ' << CS_NpNm << "        transform charge separated [N+]-[N-] to N=N\n";
  os << "  -" << zoption << ' ' << CS_SpCm << "        transform [S+]-[C-] to S=C\n";
  os << "  -" << zoption << ' ' << CS_ALLpm << "       transform all [X+]-[Y-] to X=Y\n";
  os << "  -" << zoption << ' ' << CS_XH << "          remove hydrogens\n";
  os << "  -" << zoption << ' ' << CS_AMINE << "       change all amines\n";
  os << "  -" << zoption << ' ' << CS_Om << "          protonate all O- groups\n";
  os << "  -" << zoption << ' ' << CS_Nm << "          protonate all N- groups\n";
  os << "  -" << zoption << ' ' << CS_NRMCH << "       remove all hydrogens except those to chiral centres\n";
  os << "  -" << zoption << ' ' << CS_COVM << "        break covalent bonds between Oxygen and Na,K\n";
  os << "  -" << zoption << ' ' << CS_ISOLC << "       assign formal charges to isolated Na, K, .. and Halogens\n";
  os << "  -" << zoption << ' ' << CS_GUAND << "        convert guanidines to -N-C(=N)-N form\n";
  os << "  -" << zoption << ' ' << CS_GUANDR << "        convert Ring type guanidines to -N-C(=N)-N form\n";
  os << "  -" << zoption << ' ' << CS_AZID << "        convert charge separated azids [N-]=[N+]=N to N#N=N\n";
  os << "  -" << zoption << ' ' << CS_MSDUR << "       convert misdrawn ureas, O-C(=N)-N to O=C(-N)-N\n";
  os << "  -" << zoption << ' ' << CS_FCRN << "        for converting back from corina mangled structures\n";
  os << "  -" << zoption << ' ' << CS_EHLST << "      move all explicit Hydrogen atoms to last in the connection table\n";
  os << "  -" << zoption << ' ' << CS_FMRK << "        reverse transformations applied to .mrk files\n";
  os << "  -" << zoption << ' ' << CS_RNPNM << "       transform 5 valent N=N to charge separated form\n";
  os << "  -" << zoption << ' ' << CS_FWIH << "        fix obviously wrong implicit hydrogen settings\n";
  os << "  -" << zoption << ' ' << CS_IMIDAZOLE << "   convert imidazoles to have nH near cD3\n";
  os << "  -" << zoption << ' ' << CS_PYRAZOLE << "    convert pyrazoles to have nH near electron withdrawing\n";
  os << "  -" << zoption << ' ' << CS_TRIAZOLE << "    convert triazoles to have nH near electron withdrawing\n";
  os << "  -" << zoption << ' ' << CS_TETRAZOLE << "    convert tetrazoles to have nH near attachment\n";
  os << "  -" << zoption << ' ' << CS_LACTIM_LACTAM << "        convert lactim to lactam form (non ring)\n";
  os << "  -" << zoption << ' ' << CS_LACTIM_LACTAM_RING << "       convert lactim to lactam form (ring)\n";
  os << "  -" << zoption << ' ' << CS_ISOXAZOLE << "       convert Hydroxy isoxazoles to O= forms\n";
  os << "  -" << zoption << ' ' << CS_ARGUAN << "       aromatic \"gauanidines\" - adjacent to =O, better name needed\n";
  os << "  -" << zoption << ' ' << CS_PYRAZOLONE << "       convert pyrazolone to keto form\n";
  os << "  -" << zoption << ' ' << CS_AMINO_THIAZOLE << "     convert -N=c1scc[nH]1 to -[NH]c1sccn1\n";
  os << "  -" << zoption << ' ' << CS_ALL << "         ALL the above standardistions\n";
  os << "  -" << zoption << ' ' << CS_REVERSE_NITRO << "     convert O=N=O nitro groups to charge separated\n";
  os << "  -" << zoption << ' ' << CS_REVERSE_NV5 << "       convert all 5 valent N atoms to charge separated\n";
  os << "  -" << zoption << ' ' << "APP=<xxx>   append 'xxx' to the name of changed molecules\n";
  os << "  -" << zoption << ' ' << "APP=EACH    append the reason for each change\n";

  return os.good();
}

#include "cmdline.h"

int
Chemical_Standardisation::construct_from_command_line (Command_Line & cl,
                                            int verbose,
                                            char flag)
{
  assert (ok());

  _verbose = verbose;

  if (verbose)
    set_update_chemical_standardisation_accumulators(1);

  int i = 0;
  IWString tmp;
  int rc = 0;
  while (cl.value(flag, tmp, i++))
  {
    bool negation = false;
    if ('-' == tmp[0])
    {
      negation = true;
      tmp.remove_leading_chars(1);
    }

    if ("CKV" == tmp)
    {
      _check_valence_before_and_after = 1;
    }
    else if (CS_ALL == tmp)    // do first on purpose
    {
      activate_all();
      if (verbose)
        cerr << "All chemical standardisation transformations enabled\n";

      rc++;
    }
    else if (CS_ALLpm == tmp)
    {
      _transform_plus_minus.activate();
      if (verbose)
        cerr << "CS: all charge separated [X+]-[Y-] will be transformed to X=Y\n";

      _active = 1;
      rc++;
    }
    else if (CS_NITRO == tmp)
    {
      _transform_nitro.activate();
      if (verbose)
        cerr << "CS: nitro groups will be transformed to N(=O)=O\n";

      _active = 1;
      rc++;
    }
    else if (CS_NpOm == tmp)
    {
      _transform_nplus_ominus.activate();
      if (verbose)
        cerr << "CS: charge separated [N+]-[O-] will be transformed to N=O\n";

      _active = 1;
      rc++;
    }
    else if (CS_NpNm == tmp)
    {
      _transform_n_charge_sep.activate();
      if (verbose)
        cerr << "CS: charge separated [N+]-[N-] will be transformed to N=N\n";

      _active = 1;
      rc++;
    }
    else if (CS_XH == tmp)
    {
      if (negation)
      {
        _remove_hydrogens.deactivate();
        if (verbose)
          cerr << "Explicit Hydrogens NOT removed\n";
      }
      else
      {
        _remove_hydrogens.activate();
        if (verbose)
          cerr << "Hydrogens will be removed\n";

        _active = 1;
        rc++;
      }
    }
    else if (CS_SpCm == tmp)
    {
      _transform_splus_cminus.activate();
      if (verbose)
        cerr << "[S+]-[C-] will be transformed to S=C\n";

      _active = 1;
      rc++;
    }
    else if (CS_Om == tmp)
    {
      _transform_ominus.activate();
      if (verbose)
        cerr << "All free [O-] groups will be protonated\n";

      _active = 1;
      rc++;
    }
    else if (CS_Nm == tmp)
    {
      _transform_nminus.activate();
      if (verbose)
        cerr << "All [N-] groups will be protonated\n";

      _active = 1;
      rc++;
    }
    else if (CS_AMINE == tmp)
    {
      _transform_amines.activate();
      if (verbose)
        cerr << "All charged amines will be deprotonated\n";

      _active = 1;
      rc++;
    }
    else if (CS_COVM == tmp)
    {
      _transform_covalent_metals.activate();
      if (verbose)
        cerr << "Will break bonds to covalently bonded Na and K\n";

      _active = 1;
      rc++;
    }
    else if (CS_ISOLC == tmp)
    {
      _transform_single_atom_ions.activate();

      if (verbose)
        cerr << "Isolated metals and halogens will be assigned charges\n";

      _active = 1;

      rc++;
    }
    else if (CS_GUAND == tmp)
    {
      _transform_guanidine.activate();

      if (verbose)
        cerr << "Guanidines will be transformed\n";

      _active = 1;

      rc++;
    }
    else if (CS_GUANDR == tmp)
    {
      _transform_guanidine_ring.activate();

      if (verbose)
        cerr << "Ring guanidines will be transformed\n";

      _active = 1;

      rc++;
    }
    else if (CS_NRMCH == tmp)
    {
      _remove_hydrogens.activate();
      _remove_hydrogens_attached_to_chiral_centres = 0;
      _active = 1;
      rc++;
    }
    else if (CS_ACID == tmp)
    {
      _protonate_carboxyllic_acids.activate();
      _protonate_sulfinic_acids.activate();
      _protonate_sulfonic_acids.activate();
      _protonate_sulfur_acids.activate();
      _protonate_phosphorous_acids.activate();
      _active = 1;
      rc++;

      if (verbose)
        cerr << "All acids will be protonated\n";
    }
    else if (CS_EHLST == tmp)
    {
      _explicit_hydrogens_last.activate();
      _active = 1;
      rc++;
    }
    else if (CS_FMRK == tmp)
    {
      _from_mrk_standardisations.activate();
      _active = 1;
      rc++;
    }
    else if (CS_RNPNM == tmp)
    {
      _transform_back_to_nplus_nminus.activate();
      _transform_to_charge_separated_azid.activate();
      _active = 1;
      rc++;
    }
    else if (CS_FWIH == tmp)
    {
      _transform_obvious_implicit_hydrogen_errors.activate();
      _active = 1;
      rc++;
    }
    else if (CS_AZID == tmp)
    {
      _transform_azid.activate();
      _active = 1;
      rc++;
    }
    else if (CS_MSDUR == tmp)
    {
      _transform_misdrawn_urea.activate();
      _active = 1;
      rc++;
    }
    else if (CS_FCRN == tmp)
    {
      _transform_nitro.activate();
      _transform_azid.activate();
      _transform_nplus_ominus.activate();
      _active = 1;
      rc++;
    }
    else if (CS_IMIDAZOLE == tmp)
    {
      _transform_imidazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_TETRAZOLE == tmp)
    {
      _transform_tetrazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_PYRAZOLE == tmp)
    {
      _transform_pyrazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_TRIAZOLE == tmp)
    {
      _transform_triazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_ISOXAZOLE == tmp)
    {
      _transform_isoxazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_ARGUAN == tmp)
    {
      _transform_aromatic_guanidine_ring.activate();
      _active = 1;
      rc++;
    }
    else if (CS_PYRAZOLONE == tmp)
    {
      _transform_pyrazolone.activate();
      _active = 1;
      rc++;
    }
    else if (CS_AMINO_THIAZOLE == tmp)
    {
      _transform_amino_thiazole.activate();
      _active = 1;
      rc++;
    }
    else if (CS_LACTIM_LACTAM == tmp)
    {
      _transform_lactim_lactam.activate();
      _active = 1;
      rc++;
    }
    else if (CS_LACTIM_LACTAM_RING == tmp)
    {
      _transform_lactim_lactam_ring.activate();
      _active = 1;
      rc++;
    }
    else if (CS_REVERSE_NITRO == tmp)
    {
      _transform_nitro_reverse.activate();
      _active = 1;
      rc++;
    }
    else if (CS_REVERSE_NV5 == tmp)
    {
      _transform_nv5_to_charge_separated.activate();
      _active = 1;
      rc++;
    }
    else if (tmp.starts_with("APP="))
    {
      tmp.remove_leading_chars(4);

      if ("EACH" == tmp)
      {
        _append_string_depending_on_what_changed = 1;

        if (verbose)
          cerr << "Will append standardisation details to molecule names\n";
      }
      else
      {
        _append_to_changed_molecules = tmp;

        if (_verbose)
          cerr << "Will append '" << _append_to_changed_molecules << "' to changed molecules\n";
      }
    }
    else if ("help" == tmp)
    {
      display_all_chemical_standardisation_options(cerr, flag);
      exit(0);    // note the very different behaviour in this case!!
    }
    else
    {
      cerr << "Chemical_Standardisation::construct...: unrecognised directive '" << tmp << "'\n";
      _ok = 0;
      _active = 0;
      return 0;
    }
  }

  return rc;
}

int
Chemical_Standardisation::_do_remove_hydrogens(Molecule & m)
{
  assert (ok());

  if (0 == m.natoms(1))     // quick check
    return 0;

  const int matoms = m.natoms();

  if (1 == matoms && 1 == m.atomic_number(0))   // Hydrogen molecule, do not "disappear" it...
    return 0;

  int rc = 0;

  for (int i = m.natoms() - 1; i >= 0; i--)
  {
    const Atom * a = m.atomi(i);

//  cerr << "Atom " << i << " is " << m.smarts_equivalent_for_atom(i) << endl;
    if (1 != a->atomic_number())
      continue;

//  cerr << "H atom " << i << " has ncon " << a->ncon() << endl;
    if (a->ncon() > 1)     // don't touch those bridging hydrogens
      continue;

    if (a->formal_charge() && a->ncon())     // formally charged Hydrogen attached to something
    {
      if (_verbose)
        cerr << "Chemical_Standardisation::_do_remove_hydrogens: formally charged, covalently bonded Hydrogen not removed\n";
      continue;
    }

    if (a->is_isotope())
      continue;

    if (a->ncon() > 0)     // need to check for chiral centre involvement
    {
      const atom_number_t j = a->other(i, 0);
//    cerr << "Bonded to " << j << " type " << m.smarts_equivalent_for_atom(j) << endl;

      if (NULL == m.chiral_centre_at_atom(j))   // nothing to worry about, can remove it
        ;
      else if (0 == _remove_hydrogens_attached_to_chiral_centres)  // leave these alone
        continue;
      else     // check for non-organic chirality
      {
        const Atom * aj = m.atomi(j);
//      cerr << "ORGANIC ? " << aj->element()->organic() << endl;
        if (! aj->element()->organic())
        {
//        m.set_implicit_hydrogens(j, 1, 1);
//        cerr << "After setting " << m.implicit_hydrogens(j) << endl;
//        m.set_implicit_hydrogens_known(j, 1);
//        cerr << "NON ORGANIC CHIRIALITY IH SET\n";
        }
      }
    }

//  const atom_number_t j = a->other(i, 0);
    m.remove_atom(i);
//  cerr << "AFTER removing H residual j IH " << m.implicit_hydrogens(j) << endl;
    rc++;
  }

//cerr << "After removing explicit Hydrogens\n";
//m.debug_print(cerr);

  if (rc)
  {
    _remove_hydrogens.extra(rc);

    if (_verbose) 
      cerr << "Removed " << rc << " hydrogen atoms\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:remove H";
  }

  return rc;
}

/*
  Transform [N+][O-]=O to N(=O)=O
*/

int
Chemical_Standardisation::_do_transform_nitro (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  int rc = 0;      // the number of nitro groups we change
  int natoms = m.natoms();

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  atom_number_t negative_oxygen_j = INVALID_ATOM_NUMBER;    // not really used yet...

  for (int nitrogen = 0; nitrogen < natoms; nitrogen++)
  {
    if (7 != z[nitrogen])    // central atom must be N
      continue;
    if (3 != ncon[nitrogen])            // must have three connections
      continue;

    const Atom * an = atoms[nitrogen];

    if (1 != an->formal_charge())
      continue;
    if (4 != an->nbonds())          // all cases have 4 connections
      continue;

//  Now we must examine the atoms connected.

    int oxygens_attached = 0;
    int singly_bonded_oxygen = 0;
    int singly_bonded_oxygen_with_neg_charge = 0;
    int doubly_bonded_oxygen = 0;
    atom_number_t negative_oxygen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[nitrogen]; j++)
    {
      const Bond * b = an->item(j);

      atom_number_t k = b->other(nitrogen);
      if (8 != z[k] || 1 != ncon[k])
        continue;

      oxygens_attached++;

      const Atom * ak = atoms[k];

      if (b->is_double_bond())
        doubly_bonded_oxygen++;

      else if (-1 == ak->formal_charge() && b->is_single_bond())
      {
        singly_bonded_oxygen_with_neg_charge++;
        negative_oxygen = k;
        negative_oxygen_j = j;
      }
      else if (0 == ak->formal_charge() && b->is_single_bond())
      {
        singly_bonded_oxygen++;
        negative_oxygen = k;
        negative_oxygen_j = j;
      }
    }

    if (2 != oxygens_attached || 1 != doubly_bonded_oxygen)
      continue; 

    if (singly_bonded_oxygen_with_neg_charge && singly_bonded_oxygen)   // cannot have both
      continue;

    if (singly_bonded_oxygen_with_neg_charge > 1 ||     // can have only one of either kind
        singly_bonded_oxygen > 1)
      continue;

    if (0 == singly_bonded_oxygen_with_neg_charge &&     // must have either kind
        0 == singly_bonded_oxygen)
      continue;

//  We have a nitro!

    m.set_formal_charge(nitrogen, 0);
    m.set_formal_charge(negative_oxygen, 0);
    m.set_bond_type_between_atoms(nitrogen, negative_oxygen, DOUBLE_BOND);

//  adjust the global counters - not the oxygen because maybe it was neutral..

    current_molecule_data.change_nplus(-1);
    current_molecule_data.change_npos(-1);
    rc++;
  }

  if (rc)
  {
    _transform_nitro.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " nitro groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Nitro";
  }

  return rc;
}

/*
  Transform [N+][O-] to N(=O)
*/

int
Chemical_Standardisation::_do_transform_nplus_ominus (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int nitrogen = 0; nitrogen < natoms; nitrogen++)
  {
    if (7 != z[nitrogen])    // central atom must be N
      continue;

    const Atom * an = atoms[nitrogen];

    if (1 != an->formal_charge())
      continue;

    if (3 == an->ncon() && 3 == an->nbonds())     // would introduce a valence error
      continue;

    for (int j = 0; j < ncon[nitrogen]; j++)
    {
      const Bond * b = an->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t oxygen = b->other(nitrogen);

      if (8 != z[oxygen] || 1 != ncon[oxygen])
        continue;

      if (-1 != atoms[oxygen]->formal_charge())
        continue;

//    We have an [N+]-[O-] pair

      m.set_formal_charge(nitrogen, 0);
      m.set_formal_charge(oxygen, 0);
      m.set_bond_type_between_atoms(nitrogen, oxygen, DOUBLE_BOND);
      rc++;

      current_molecule_data.change_nplus(-1);
      current_molecule_data.change_ominus(-1);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);


      break;     // this N is no longer positively charged! - Comgenex example of [O-][N+]([O-])=O for a nitro group
    }
  }

  if (rc)
  {
    _transform_nplus_ominus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [N+]-[O-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[N+]-[O-]";
  }

  return rc;
}

/*
  After all other transformations are done, there may be other [O-]
  atoms to transform

  Jul 97, also do S-
  Sept 98. Make sure these are singly connected - it was changing
  things like C[N+]1=CC(=NC)[O-]=N1 
  which destroyed the aromaticity
*/

int
Chemical_Standardisation::_do_transform_ominus (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();
  for (int i = 0; i < natoms; i++)
  {
    if (8 == z[i])
      ;
    else if (16 == z[i])
      ;
    else
      continue;

    if (-1 != atoms[i]->formal_charge())
      continue;

    if (2 == ncon[i])
      continue;

    m.set_formal_charge(i, 0);
    current_molecule_data.change_nneg(-1);
    if (8 == z[i])
      current_molecule_data.change_ominus(-1);
    else
      current_molecule_data.change_sminus(-1);

    rc++;
  }

  if (rc)
  {
    _transform_ominus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [O-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[O-]";
  }

  return rc;
}

static int
bonded_to_positively_charged_anything (const Molecule & m,
                 atom_number_t nitrogen)
{
  const Atom * a = m.atomi(nitrogen);

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(nitrogen, i);
    
    const Atom * aj = m.atomi(j);

    if (aj->formal_charge() <= 0)
      continue;

//  if (7 != aj->atomic_number())
//    continue;

    return 1;
  }

  return 0;
}

/*
  Change N- that aren't next to anything positive
*/

int
Chemical_Standardisation::_do_transform_nminus (Molecule & m,
                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();
  for (int nitrogen = 0; nitrogen < natoms; nitrogen++)
  {
    if (7 != z[nitrogen])
      continue;

    const Atom * n = atoms[nitrogen];

    if (-1 != n->formal_charge())
      continue;

    if (bonded_to_positively_charged_anything(m, nitrogen))
      continue;

    m.set_formal_charge(nitrogen, 0);
    rc++;
  }

  if (rc)
  {
    _transform_nminus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [N-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[N-]";
  }

  return rc;
}

/*
  Transform [N+]-[N-] to N=N
  Make sure that the N+ has 4 bonds, and the N- two

  Ran into a big problem with this in the case where the N+ and N-
  are in an aromatic ring, so we don't do the transformation in
  those circumstances.

  Oct 2001. Don't want to do things like

  C(=[N+]1CC(=O)[N-]1)(C)C
*/

static int
ok_to_change_charge_separated_pair (const Molecule & m,
                        atom_number_t a1,
                        atom_number_t a2)
{
  const Atom * aa1 = m.atomi(a1);
  if (7 != aa1->atomic_number())     // not a Nitrogen, OK
    return 1;

  const Atom * aa2 = m.atomi(a2);
  if (7 != aa2->atomic_number())     // not a Nitrogen, OK
    return 1;

// we have a charge separated pair of Nitrogens

  int a1con = m.ncon(a1);
  int a2con = m.ncon(a2);

  if (2 == a1con && 1 == a2con)     // An azide in the form =[N+]=[N-]
    return 1;

  return 0;     // all other cases, don't change
}

int
Chemical_Standardisation::_do_transform_n_charge_sep (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change

// Ran into problems with LLY 373212 where the transformation
// destroyed aromaticity. I then fixed the aromaticity rules, but
// to be careful, let's recompute the aromaticity if likely

  int need_to_recompute_aromaticity = 0;

  int natoms = m.natoms();
  for (int n1 = 0; n1 < natoms; n1++)
  {
    if (7 != z[n1])    // central atom must be N
      continue;

    Atom * a1 = const_cast<Atom *>(atoms[n1]);

    if (1 != a1->formal_charge())
      continue;
    if (4 != (a1->nbonds() + a1->implicit_hydrogens()))
      continue;

//  beware cases like C12=CC=CC=C1[N-][N+H2][N-]2 p10, don't change them

    int negative_nitrogens_attached = 0;   // multiple N- attached to one N+
    int negative_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[n1]; j++)
    {
      const Bond * b = a1->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t n2 = b->other(n1);

      if (7 != z[n2])
        continue;

      Atom * an2 = const_cast<Atom *>(atoms[n2]);

      if (-1 != an2->formal_charge())
        continue;

      if (2 != (an2->nbonds() + an2->implicit_hydrogens()))
        continue;

//    We have an [N+]-[N-] pair. Are they in an aromatic ring
     
//    if (m.in_same_aromatic_ring(n1, n2))    August 2002, decided to drive out as much charge separation as possible
//      continue;

      negative_nitrogens_attached++;
      if (1 == negative_nitrogens_attached)
        negative_nitrogen = n2;
    }

    if (1 == negative_nitrogens_attached)
    {
      m.set_formal_charge(n1, 0);
      m.set_formal_charge(negative_nitrogen, 0);
      m.set_bond_type_between_atoms(n1, negative_nitrogen, DOUBLE_BOND);
      current_molecule_data.change_nplus(-1);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);
      rc++;
    }
  }

  if (need_to_recompute_aromaticity)
    m.compute_aromaticity();

  if (rc)
  {
    _transform_n_charge_sep.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [N+]-[N-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[N+]-[N-]";
  }

  return rc;
}

static int
two_negatively_charged_connections (atom_number_t zatom,
                                    const Atom & a,
                                    const Atom * const * atoms,
                                    atom_number_t & n1,
                                    atom_number_t & n2)
{
  n1 = INVALID_ATOM_NUMBER;
  n2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < a.ncon(); i++)
  {
    const Bond * b = a.item(i);
    if (! b->is_single_bond())
      continue;

    atom_number_t n = b->other(zatom);

    if (-1 != atoms[n]->formal_charge())
      continue;

    if (INVALID_ATOM_NUMBER == n1)
      n1 = n;
    else
      n2 = n;
  }

  if (INVALID_ATOM_NUMBER == n1 || INVALID_ATOM_NUMBER == n2)
    return 0;

  return 1;
}

void
Chemical_Standardisation::_do_transform_plus_minus_pair (Molecule & m,
                                 atom_number_t a1,
                                 atom_number_t a2,
                                 IWStandard_Current_Molecule & current_molecule_data)
{
  m.set_formal_charge(a1, 0);
  m.set_formal_charge(a2, 0);
  m.set_bond_type_between_atoms(a1, a2, DOUBLE_BOND);

  current_molecule_data.change_npos(-1);
  current_molecule_data.change_nneg(-1);

  return;
}

/*
  Transform all [X+]-[Y-] to X=Y
*/

int
Chemical_Standardisation::_do_transform_plus_minus (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;      // the number of groups we change
  int natoms = m.natoms();
  for (int i = 0; i < natoms; i++)
  {
    const Atom * ai = atoms[i];

    if (1 != ai->formal_charge())
      continue;

    atom_number_t n1, n2;
    (void) two_negatively_charged_connections(i, *ai, atoms, n1, n2);

    if (INVALID_ATOM_NUMBER == n1)    // no negatively charged neighbours
      continue;

    if (INVALID_ATOM_NUMBER == n2)    // only N1 is a negatively charged neighbour
    {
      if (ok_to_change_charge_separated_pair(m, i, n1))
        continue;

      _do_transform_plus_minus_pair(m, i, n1, current_molecule_data);
      rc++;
      continue;
    }

//  We have two negatively charged neighbours. Which one should we use?

    int ok1 = ok_to_change_charge_separated_pair(m, i, n1);
    int ok2 = ok_to_change_charge_separated_pair(m, i, n2);
//  cerr << ok1 << " and " << ok2 << endl;
    if (ok1 && ! ok2)
      _do_transform_plus_minus_pair(m, i, n1, current_molecule_data);
    else if (ok2 && ! ok1)
      _do_transform_plus_minus_pair(m, i, n2, current_molecule_data);
    else if (atoms[n1]->nbonds() < atoms[n2]->nbonds())
      _do_transform_plus_minus_pair(m, i, n1, current_molecule_data);
    else if (atoms[n1]->nbonds() > atoms[n2]->nbonds())
      _do_transform_plus_minus_pair(m, i, n2, current_molecule_data);
    else
      _do_transform_plus_minus_pair(m, i, n2, current_molecule_data);

    rc++;
  }

  if (rc)
  {
    _transform_plus_minus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [X+]-[Y-] groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[X+]-[Y-]";
  }

  return rc;
}

/*
  Change C(=O)[O-] to C(=O)O
*/

int
Chemical_Standardisation::_do_protonate_carboxyllic_acids (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a carbon.

    const Bond * b = ai->item(0);

    if (! b->is_single_bond())
      continue;

    atom_number_t c = b->other(i);

    if (6 != z[c])
      continue;

    if (3 != ncon[c])
      continue;

    const Atom * ac = atoms[c];

    int found_doubly_bonded_oxygen = 0;
    for (int j = 0; j < 3; j++)
    {
      const Bond * b2 = ac->item(j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t o = b2->other(c);

      if (8 == z[o] && 1 == ncon[o])
      {
        found_doubly_bonded_oxygen = 1;
        break;
      }
    }

    if (! found_doubly_bonded_oxygen)
      continue;

    m.set_formal_charge(i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_ominus(-1);
  }

  if (rc)
  {
    _protonate_carboxyllic_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " carboxyllic acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[OH]-C=O";
  }

  return rc;
}

/*
  Change [O-,S-]-P=[O,S] to [O,S]-P=[O,S]
*/

int
Chemical_Standardisation::_do_protonate_phosphorous_acids (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (15 != z[i])
      continue;

    if (ncon[i] < 3)
      continue;

    const Atom * ai = atoms[i];

//  Find a negatively charged Oxygen or Sulphur, and a doubly bonded O or S - let's hope we don't see too many P~S bonds, ugly!

    atom_number_t negatively_charged_OS = INVALID_ATOM_NUMBER;
    atom_number_t doubly_bonded_OS = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item(j);

      atom_number_t k = b->other(i);

      const Atom * ak = atoms[k];

      if (8 == ak->atomic_number())
        ;
      else if (16 == ak->atomic_number())
        ;
      else
        continue;

      if (b->is_single_bond())
      {
        if (-1 == ak->formal_charge())
          negatively_charged_OS = k;
      }
      if (b->is_double_bond())
        doubly_bonded_OS = k;
      else if (b->is_single_bond() && -1 == ak->formal_charge())
        negatively_charged_OS = k;
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_OS || INVALID_ATOM_NUMBER == negatively_charged_OS)
      continue;

    m.set_formal_charge(negatively_charged_OS, 0);
    rc++;

    current_molecule_data.change_nneg(-1);

    if (8 == atoms[negatively_charged_OS]->atomic_number())
      current_molecule_data.change_ominus(-1);
    else
      current_molecule_data.change_sminus(-1);

    current_molecule_data.change_phosphorus(-1);
    if (0 == current_molecule_data.phosphorus())
      break;
  }

  if (rc)
  {
    _protonate_phosphorous_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " phosphorus acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[OH,SH]-P=[O,S]";
  }

  return rc;
}

/*
  Change [S-]-[C,P]=[O,S] to S-[C,P]=[O,S]
*/

int
Chemical_Standardisation::_do_protonate_sulfur_acids (Molecule & m,
                             IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (16 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

    const Bond * b = ai->item(0);

    if (! b->is_single_bond())
      continue;

    atom_number_t cps = b->other(i);

    if (6 == z[cps])
      ;
    else if (15 == z[cps])
      ;
    else if (16 == z[cps])
      ;
    else
      continue;

    if (ncon[cps] < 3)
      continue;

    const Atom * acps = atoms[cps];

    int found_doubly_bonded_OS = 0;
    for (int j = 0; j < ncon[cps]; j++)
    {
      const Bond * b2 = acps->item(j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t os = b2->other(cps);

      if (1 != ncon[os])
        continue;

      if (8 == z[os] || 16 == z[os])
      {
        found_doubly_bonded_OS = 1;
        break;
      }
    }

    if (! found_doubly_bonded_OS)
      continue;

    m.set_formal_charge(i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_sminus(-1);
  }
  if (rc)
  {
    _protonate_sulfur_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " Sulphur acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[SH]-[C,P,S]=[O,S]";
  }

  return rc;
}

/*
 April 2014. Expand this to handle

 O-[S+](=O)
*/

int
Chemical_Standardisation::_do_protonate_sulfonic_acids (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

//  if (-1 != ai->formal_charge())
//    continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a sulphur.

    const Bond * b1 = ai->item(0);

    if (! b1->is_single_bond())
      continue;

    atom_number_t s = b1->other(i);

    if (16 != z[s])
      continue;

    if (4 != ncon[s])
      continue;

    const Atom * as = atoms[s];

    int doubly_bonded_oxygen_count = 0;
    int ih_known = 0;
    for (int j = 0; j < 3; j++)
    {
      const Bond * b2 = as->item(j);
      if (! b2->is_double_bond())
        continue;

      const atom_number_t o = b2->other(s);

      if (8 != z[o])
        continue;

      if (1 != ncon[o])
        continue;

      doubly_bonded_oxygen_count++;

      if (atoms[o]->implicit_hydrogens_known())
        ih_known++;
    }

    if (0 == doubly_bonded_oxygen_count)
      continue;

    if (ih_known)     // [S+]([O])(=O)(OCCCNC(=O)C=C)ON PBCHM101160818  oxygen atom is radical???
      continue;

    if (2 == doubly_bonded_oxygen_count && 1 == ai->formal_charge())    
    {
      m.set_formal_charge(i, 0);
      current_molecule_data.change_nneg(-1);
      current_molecule_data.change_ominus(-1);
    }
    else if (1 == doubly_bonded_oxygen_count && 1 == as->formal_charge() && 0 == ai->formal_charge())
    {
      m.set_formal_charge(s, 0);
      m.set_bond_type_between_atoms(i, s, DOUBLE_BOND);
      current_molecule_data.change_splus(-1);
      current_molecule_data.change_npos(-1);
    }
    else
      continue;

    rc++;
  }

  if (rc)
  {
    _protonate_sulfonic_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " sulfonic acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Sulfonic";
  }

  return rc;
}

/*
  A sulfinic acid is S(=O)O

  April 2014. Change to also process -[S+](-O)- from NIH
*/

int
Chemical_Standardisation::_do_protonate_sulfinic_acids (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

//  if (-1 != ai->formal_charge())
//    continue;

//  look for the Sulphur

    const Bond * b = ai->item(0);
    if (! b->is_single_bond())
      continue;

    atom_number_t s = b->other(i);

    if (16 != z[s])
      continue;

    if (3 != ncon[s])
      continue;

    const Atom * as = atoms[s];

    if (1 == as->formal_charge() && 3 == as->nbonds() && 0 == ai->formal_charge())   // create a sulfinic acid
    {
      m.set_formal_charge(s, 0);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_splus(-1);
      m.set_bond_type_between_atoms(i, s, DOUBLE_BOND);
      rc++;
      continue;
    }

    if (3 == as->nbonds() || 0 == ai->formal_charge())      // no double bonds here, and/or no charge to neutralise
      continue;

    int doubly_bonded_oxygen_count = 0;
    for (int j = 0; j < 3; j++)
    {
      const Bond * b2 = as->item(j);
      if (! b2->is_double_bond())
        continue;

      atom_number_t o = b2->other(s);

      if (8 == z[o] && 1 == ncon[o])
        doubly_bonded_oxygen_count++;
    }

    if (1 != doubly_bonded_oxygen_count)
      continue;

    if (1 == ai->formal_charge() && 0 == as->formal_charge())
    {
      m.set_formal_charge(i, 0);
      current_molecule_data.change_nneg(-1);
      current_molecule_data.change_ominus(-1);
    }
    else
      continue;

    rc++;
  }

  if (rc)
  {
    _protonate_sulfinic_acids.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " sulfinic acids\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Sulfinic";
  }

  return rc;
}

int
Chemical_Standardisation::_do_protonate_no (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

//  At this stage, we have a singly bonded, negatively charged Oxygen. It should
//  be bonded to a Nitrogen.

    const Bond * b = ai->item(0);
    if (! b->is_single_bond())
      continue;

    atom_number_t n = b->other(i);;
    if (7 != z[n])
      continue;

    if (0 != m.formal_charge(n))   // hard to know what might be going on. The [N+]-[O-] rule should have got it
      continue;

    m.set_formal_charge(i, 0);
    rc++;

    current_molecule_data.change_nneg(-1);
    current_molecule_data.change_ominus(-1);
  }

  if (rc)
  {
    _protonate_no.extra(rc);
    if (_verbose)
      cerr << "Protonated " << rc << " NO- groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:N-[O-]";
  }

  return rc;
}

/*
  Transform [C-]-[S+] to C=S
  if the S+ has 5 bonds
  October 2000. Relax the restrictions on the environment of the S. Now allow either
    4 connections and 5 bonds
  or
    3 connections and 3 bonds

  Feb 2001. Extend so that we also transform [O-]-[S+]
*/

int
Chemical_Standardisation::_do_transform_splus_cminus (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (16 != z[i])
      continue;

    const Atom * ai = atoms[i];

    if (1 != ai->formal_charge())
      continue;

    if (4 == ncon[i] && 5 == ai->nbonds())     // yep, do this one
      ;
    else if (3 == ncon[i] && 3 == ai->nbonds())     // yep, do this one
      ;
    else
      continue;

//  At this stage, we have an appropriate S+. Is it bonded to a C- or O-

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item(j);
      if (! b->is_single_bond())
        continue;

      atom_number_t c = b->other(i);
      if (6 == z[c])
        ;
      else if (8 == z[c])
        ;
      else
        continue;

      if (-1 == atoms[c]->formal_charge())
      {
        current_molecule_data.change_nneg(-1);
        m.set_formal_charge(c, 0);
      }
      else if (8 == z[c] && 1 == ncon[c])    // OK to have a neutral oxygen here
        ;
      else
        continue;

      m.set_formal_charge(i, 0);
      m.set_bond_type_between_atoms(i, c, DOUBLE_BOND);
      rc++;

      if (NULL != m.chiral_centre_at_atom(c))
        m.remove_chiral_centre_at_atom(c);

      current_molecule_data.change_npos(-1);

      break;
    }
  }

  if (rc)
  {
    _transform_splus_cminus.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " [S+]-[C-]\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:[S+]-[C-]";
  }

  return rc;
}

/*
  This next one is greatly complicated by the possibility of explicit hydrogens
*/

int
Chemical_Standardisation::_do_transform_amines (Molecule & m,
                                        Set_of_Atoms & atoms_to_be_removed,
                                        IWStandard_Current_Molecule & current_molecule_data) 
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    Atom * a = const_cast<Atom *>(atoms[i]);

    if (1 != a->formal_charge())
      continue;

//  It must have a hydrogen.

    int hc = m.hcount(i);     // check both implicit and explicit

    if (0 == hc)     // must have at least 1 Hydrogen attachment
      continue;

    int implicit_hydrogens;
    if (hc)
      implicit_hydrogens = a->implicit_hydrogens();
    else
      implicit_hydrogens = 0;

//  If all hydrogens are explicit, we must find one to remove. Check also for
//  adjacent N- which do not get processed

    atom_number_t explicit_hydrogen_to_remove = INVALID_ATOM_NUMBER;
    int found_negative_nitrogen = 0;

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);
      if (1 == z[k])
        explicit_hydrogen_to_remove = k;
      else if (7 == z[k] && -1 == atoms[k]->formal_charge())
      {
        found_negative_nitrogen = 1;
        break;
      }
    }

    if (found_negative_nitrogen)
      continue;

    if (implicit_hydrogens)
      ;
    else if (INVALID_ATOM_NUMBER == explicit_hydrogen_to_remove)   // should not happen
      continue;
    else
      atoms_to_be_removed.add_if_not_already_present(explicit_hydrogen_to_remove);

    m.set_formal_charge(i, 0);

    rc++;

    current_molecule_data.change_nplus(-1);
    current_molecule_data.change_npos(-1);
  }

  if (rc)
  {
    _transform_amines.extra(rc);
    if (_verbose)
      cerr << "Transformed " << rc << " protonated amines\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:amines";
  }

  return rc;
}

int
Chemical_Standardisation::_do_transform_covalent_metals (Molecule & m)
{
  const int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    if (1 != a->ncon())       // we only process single atoms
      continue;

    const atomic_number_t z = a->atomic_number();

    if (11 == z || 19 == z)     // Na and K only
    {
      atom_number_t o = a->other(i, 0);

      const Atom * ao = m.atomi(o);

      if (8 == ao->atomic_number())
        ;
      else if (16 == ao->atomic_number())
        ;
      else
        continue;

      if (2 != ao->ncon())
        continue;

      m.remove_bond_between_atoms(i, o);
      m.set_formal_charge(i, 1);
      m.set_formal_charge(o, -1);

      rc++;
    }
  }

  if (rc)
  {
    _transform_covalent_metals.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " covalently bonded metals\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:cv_metal";
  }

  return rc;
}


#ifdef OLD_VERSION_USING_MOLECULE_CACHE
int
Chemical_Standardisation::_do_transform_covalent_metals (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  int matoms = m.natoms();

  int rc = 0;

  const atomic_number_t * z = current_molecule_data.atomic_number();
  int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int i = 0; i < matoms; i++)
  {
    if (1 != ncon[i])       // we only process single atoms
      continue;

    if (11 == z[i] || 19 == z[i])     // Na and K only
    {
      atom_number_t o = atoms[i]->other(i, 0);

      if (8 == z[o])
        ;
      else if (16 == z[o])
        ;
      else
        continue;

      if (2 != ncon[o])
        continue;

      m.remove_bond_between_atoms(i, o);
      m.set_formal_charge(i, 1);
      m.set_formal_charge(o, -1);
      ncon[i] = 0;
      ncon[o] = 1;
      current_molecule_data.change_nneg(1);
      if (8 == z[o])
        current_molecule_data.change_ominus(1);
      else
        current_molecule_data.change_sminus(1);

      rc++;

      current_molecule_data.change_singly_connected_metal(-1);

      if (0 == current_molecule_data.singly_connected_metal())
        break;
    }
  }

  if (rc)
  {
    _transform_covalent_metals.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " covalently bonded metals\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:cv_metal";
  }

  return rc;
}
#endif

static int
resolved (const int * score,
          atom_number_t & n0,
          atom_number_t & n1,
          atom_number_t & n2)
{
// If any two are the same, we are not resolved

  if (score[0] == score[1] || score[0] == score[2] || score[1] == score[2])
    return 0;

// The must be all different

  std::pair<int, int> tmp[3];
  for (int i = 0; i < 3; ++i)
  {
    tmp[i].second = score[i];
  }

  tmp[0].first = n0;
  tmp[1].first = n1;
  tmp[2].first = n2;

  std::sort(tmp, tmp + 3, [] (const std::pair<int, unsigned int> & p1, const std::pair<int, unsigned int> & p2)
                                        {
                                          return p1.second < p2.second;
                                        });
//cerr << "sorted " << tmp[0].second << ' ' << tmp[1].second << ' ' << tmp[2].second << endl;
  n0 = tmp[0].first;
  n1 = tmp[1].first;
  n2 = tmp[2].first;

  return 1;
}

static int
initialise (const Molecule & m,
            const atom_number_t zatom,
            const int * claimed,
            Set_of_Atoms & edge_atoms,
            const IWStandard_Current_Molecule & current_molecule_data)
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

  int rc = 0;

//cerr << "initialise from atom " << zatom <<endl;

  for (int i = 0; i < acon; ++i)
  {
    const atom_number_t j = a->other(zatom, i);

    if (claimed[j])
      continue;

    edge_atoms.add(j);

    rc += 5 * current_molecule_data.atomic_number()[j] + current_molecule_data.ncon()[j];

//  cerr << "Atom " << j << " attached, rc " << rc << endl;
  }

  return rc;
}

static int
advance (const Molecule & m,
         const Set_of_Atoms & edge_atoms,
         Set_of_Atoms & next_edge_atoms,
         const int * claimed,
         IWStandard_Current_Molecule & current_molecule_data)
{
  const Atom * const * atoms = current_molecule_data.atoms();
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();

  const int n = edge_atoms.number_elements();

  next_edge_atoms.resize_keep_storage(0);

  if (0 == n)
    return 0;

  int rc = 0;

  for (int i = 0; i <n; ++i)
  {
    const atom_number_t j = edge_atoms[i];

    const Atom * aj = atoms[j];

    const int jcon = aj->ncon();

    for (int k = 0; k < jcon; ++k)
    {
      const atom_number_t l = aj->other(j, k);
      if (claimed[l])
        continue;

      rc += 5 * z[l] + ncon[l];
      next_edge_atoms.add_if_not_already_present(l);
    }
  }

  return rc;
}

static void
write_numbered_smiles (const Molecule & m,
                       std::ostream & output)
{
  Molecule mcopy(m);
  for (int i = 0; i < m.natoms(); ++i)
  {
    mcopy.set_isotope(i, i);
  }
  output << mcopy.smiles();

  return;
}

static int
identify_the_one_still_active(const int * score,
                              atom_number_t & n0,
                              atom_number_t & n1,
                              atom_number_t & n2)
{
  for (int i = 0; i < 3; ++i)
  {
  }
//cerr << "identify_the_one_still_active scores " << score[0] << ' ' << score[1] << ' ' << score[2] << endl;
  if (score[0] == score[1])    // n2 is distinct
    std::swap(n0, n2);
  else if (score[0] == score[2])  // n1 is distinct
    std::swap(n0, n1);
  else
  {
// default case means 1 and 2 the same

    assert (score[1] == score[2]);
  }

  return 1;
}

/*
  Kind of like what we do in path_scoring, but we do not look at the bond types
*/

static int
resolve_by_shell_expansion(const Molecule & m,
                           atom_number_t & n1,
                           atom_number_t & n2,
                           atom_number_t & n3,
                           const atom_number_t carbon,
                           IWStandard_Current_Molecule & current_molecule_data)
{
  const int matoms = m.natoms();

  int * claimed = new_int(matoms); std::unique_ptr<int[]> free_claimed(claimed);
  claimed[carbon] = 1;
  claimed[n1] = 1;
  claimed[n2] = 1;
  claimed[n3] = 1;

  Set_of_Atoms n(3);
  n.add(n1);
  n.add(n2);
  n.add(n3);

//#define DEBUG_RESOLVE_BY_SHELL_EXPANSION
#ifdef DEBUG_RESOLVE_BY_SHELL_EXPANSION
  cerr << "resolve_by_shell_expansion atoms " << n1 << ' ' << m.smarts_equivalent_for_atom(n1) << " n2 " << n2 << ' ' << m.smarts_equivalent_for_atom(n2) << " n3 " << n3 << ' ' << m.smarts_equivalent_for_atom(n3) << endl;
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << " starting molecule\n";
#endif


  Set_of_Atoms e1[3], e2[3];   // edge atoms for each iteration. We keep two sets and alternate between them - using one as current and the other as next

  int score[3];
  for (int i = 0; i < 3; ++i)
  {
    score[i] = initialise(m, n[i], claimed, e1[i], current_molecule_data);
#ifdef DEBUG_RESOLVE_BY_SHELL_EXPANSION
    cerr << " i = " << i << " atom " << n[i] << " score " << score[i] << endl;
#endif
  }

  if (resolved(score, n1, n2, n3))
  {
#ifdef DEBUG_RESOLVE_BY_SHELL_EXPANSION
    cerr << "Resolved, n1 is " << n1 << ", n2 " << n2 << " n3 " << n3 << endl;
#endif
    return 1;
  }

// We alternate between the e1 and e2 arrays - one cycle the current edge atoms and the next edge atoms

  Set_of_Atoms * p1 = e1;
  Set_of_Atoms * p2 = e2;

  const int max_iterations = 10;

  int nactive = 3;

  for (int i = 1; i <= max_iterations; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if (score[j] > 0)   // still active, compute new score
      {
        score[j] = advance(m, p1[j], p2[j], claimed, current_molecule_data);
        if (0 == score[j])    // mark as exhausted
        {
          nactive--;
          score[j] = - i;    // mark as being exhausted
        }
      }
    }

#ifdef DEBUG_RESOLVE_BY_SHELL_EXPANSION
    cerr << " iteration " << i << " resolved? " << resolved(score, n1, n2, n3) << ", nactive " << nactive << " " << score[0] << ' ' << score[1] << ' ' << score[2] << endl;
#endif
    if (resolved(score, n1, n2, n3))
      return 1;

    if (nactive <= 1)
      return identify_the_one_still_active(score, n1, n2, n3);

    for (int j = 0; j < 3; ++j)
    {
      p1[j].set_vector(claimed, 1);
    }

#ifdef DEBUG_RESOLVE_BY_SHELL_EXPANSION
    for (int j = 0; j < 3; ++j)
    {
      cerr << " iteration " << i << " j = " << j << " score " << score[j] << " E1 " << e1[j] << " E2 " << e2[j] << endl;
    }
#endif

    std::swap(p1, p2);
  }

  return 0;
}

static int
collect_attached_nitrogen_atoms (Molecule & m,
                                 const atom_number_t c,
                                 atom_number_t & first_nh,
                                 atom_number_t & second_nh,
                                 atom_number_t & doubly_bonded_n,
                                 const IWStandard_Current_Molecule & current_molecule_data)
{
  first_nh = INVALID_ATOM_NUMBER;
  second_nh = INVALID_ATOM_NUMBER;
  doubly_bonded_n = INVALID_ATOM_NUMBER;

  const Atom * ac = m.atomi(c);

  const atomic_number_t * z = current_molecule_data.atomic_number();

  if (3 != ac->ncon())   // should not happen
    return 0;

  for (int i = 0; i < 3; ++i)
  {
    const Bond * b = ac->item(i);
    if (b->is_directional())
      return 0;

    const atom_number_t n = b->other(c);
    if (7 != z[n])    // should not happen
      return 0;

    if (b->is_single_bond() && m.hcount(n))
    {
      if (INVALID_ATOM_NUMBER == first_nh)
        first_nh = n;
      else
        second_nh = n;
    }
    else if (b->is_double_bond())
      doubly_bonded_n = n;
  }

  if (INVALID_ATOM_NUMBER == second_nh || INVALID_ATOM_NUMBER == doubly_bonded_n)
    return 0;
  
  return 1;
}

/*
  We convert -N=C(-NH2)-NH2  to -N-C(=NH)-NH2

  March 2016. Extend this to the generic case of an unsaturated carbon surrounded by three
  Nitrogen atoms.
  Note that it will fail for molecules like

  N(=C(N)N)C(=NC(C)C)N

  where it will random as to which one gets standardised first.
  We could mitigate this problem by first doing the "terminal" guanidines, then the interior
  ones, but these molecules are very rare...
*/

int
Chemical_Standardisation::_do_transform_guanidine (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * ring_membership = current_molecule_data.ring_membership();

  const Set_of_Atoms & g = current_molecule_data.possible_guanidine();

  const int ng = g.number_elements();

  if (0 == ng)
    return 0;

  int rc = 0;

  for (int i = ng - 1; i >= 0; --i)    // because we remove things from the array
  {
    const atom_number_t c = g[i];
    if (ring_membership[c])
      continue;

    assert (6 == z[c] && 3 == ncon[c]);

    atom_number_t first_nh, second_nh, doubly_bonded_n;
    if (! collect_attached_nitrogen_atoms(m, c, first_nh, second_nh, doubly_bonded_n, current_molecule_data))
      continue;

    if (1 == ncon[doubly_bonded_n])
      continue;

    if (1 == ncon[second_nh])
      ;
    else if (1 == ncon[first_nh])
      std::swap(first_nh, second_nh);
    else
    {
//    write_numbered_smiles(m, cerr);
//    cerr << " resolve by shell expansion\n";

      resolve_by_shell_expansion(m, second_nh, first_nh, doubly_bonded_n, c, current_molecule_data);
//    cerr << first_nh << ' ' << second_nh << ' ' << doubly_bonded_n << endl;
    }

    m.set_bond_type_between_atoms(c, doubly_bonded_n, SINGLE_BOND);
    m.set_bond_type_between_atoms(c, second_nh, DOUBLE_BOND);
    m.set_bond_type_between_atoms(c, first_nh, SINGLE_BOND);
    m.set_implicit_hydrogens_known(second_nh, 0);
    m.set_implicit_hydrogens_known(doubly_bonded_n, 0);

//  Lose any formal charge already applied

    if (1 == atoms[doubly_bonded_n]->formal_charge())
      m.set_formal_charge(doubly_bonded_n, 0);

    current_molecule_data.remove_possible_guanidine(c);

    rc++;
  }

  if (rc)
  {
    _transform_guanidine.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " guanidine groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:guanidine";
  }

  return rc;
}

static int
compute_guanidine_bond_acceptance_desirability (Molecule & m,
                                                atom_number_t n,
                                                const int * ring_membership,
                                                const int * ncon)
{
  if (1 == ncon[n])
    return 100;

  if (ring_membership[n])
    return 50 + 4 * ncon[n] + m.attached_heteroatom_count(n);

  return 4 * ncon[n] + m.attached_heteroatom_count(n);
}

#ifdef OLDER_VERSIONJ
static int
compute_guanidine_bond_acceptance_desirability (Molecule & m,
                                                atom_number_t n,
                                                const int * ring_membership,
                                                const int * ncon)
{
  if (0 == m.hcount(n))
    return 0;
  if (1 == ncon[n])
    return 10;
  if (0 == ring_membership[n])
    return 5 + m.attached_heteroatom_count(n);

   return m.attached_heteroatom_count(n) + ncon[n];
}
#endif

#define DEBUG_RING_GUANIDINE_STANDARDISE

/*
  We have something like


  N     N
   \   /
    \ /
     C
     ||
     N

  try to get the N=C bond into the ring
*/

int
Chemical_Standardisation::_do_transform_ring_guanidine(Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  if (0 == m.nrings())
    return 0;

  int rc = 0;

  const Atom * const * atoms = current_molecule_data.atoms();
  const int * ring_membership = current_molecule_data.ring_membership();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

//m.compute_aromaticity();    // must compute it. Molecule may have aromaticity definition computed with Pearlman rules

  const Set_of_Atoms & g = current_molecule_data.possible_guanidine();

  const int ng = g.number_elements();

  if (0 == ng)
    return 0;

  for (int i = ng - 1; i >= 0; --i)
  {
    const atom_number_t c = g[i];

    if (atom_is_aromatic[c])
      continue;

    if (0 == ring_membership[c])
      continue;

    atom_number_t n1, n2, doubly_bonded_n;
    if (! collect_attached_nitrogen_atoms(m, c, n1, n2, doubly_bonded_n, current_molecule_data))
      continue;

    if (atom_is_aromatic[n1] || atom_is_aromatic[n2])
      continue;

    if (ring_membership[doubly_bonded_n])    // too hard, probably taken care of elsewhere
      continue;

    const int h1 = m.hcount(n1);
    const int h2 = m.hcount(n2);

    if (0 == h1 && 0 == h2)     // neither of the other atoms can accept a double bond
      continue;

// too hard. Thought about doing via expansion, but run into problems with other groups
// in ring. Even this is buggy because other transformations may change the H count

    if (h1 && h2)
      continue;

    if (m.hcount(n1))
      m.set_bond_type_between_atoms(c, n1, DOUBLE_BOND);
    else if (m.hcount(n2))
      m.set_bond_type_between_atoms(c, n2, DOUBLE_BOND);
    else                                      // surprising
      continue;

    m.set_bond_type_between_atoms(c, doubly_bonded_n, SINGLE_BOND);
    if (1 == atoms[doubly_bonded_n]->formal_charge())
      m.set_formal_charge(doubly_bonded_n, 0);

    current_molecule_data.remove_possible_guanidine(c);

    rc++;
  }

  if (rc)
  {
    _transform_guanidine_ring.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " ring guanidine groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:Rguanidine";
  }

  return rc;
}

static int
singly_bonded_exocyclic_connection(Molecule & m,
                                   const Set_of_Atoms & r,
                                   const atom_number_t zatom,
                                   bond_type_t & bt,
                                   atom_number_t & exocyclic)
{
  const Atom * a = m.atomi(zatom);

  assert (3 == a->ncon());

  for (int i = 0; i < 3; ++i)
  {
    const Bond * b = a->item(i);

    const atom_number_t j = b->other(zatom);

    if (1 != m.ncon(j))
      continue;

    exocyclic = j;
    if (b->is_single_bond())
      bt = SINGLE_BOND;
    else
      bt = DOUBLE_BOND;

    return 1;
  }

  return 0;
}

static atom_number_t
atom_is_connected_to (const Molecule & m,
                      const atom_number_t zatom,
                      const atomic_number_t z, 
                      const bond_type_t bt)
{
  const Atom * a = m.atomi(zatom);

  const int acon = a->ncon();

//cerr << "Search from atom " << zatom << " " << m.smarts_equivalent_for_atom(zatom) << endl;

  for (int i = 0; i < acon; ++i)
  {
    const Bond * b = a->item(i);
//  cerr << " i = " << i << " type " << BOND_TYPE_ONLY(b->btype()) << " cmp " << bt << endl;

    if (bt != BOND_TYPE_ONLY(b->btype()))
      continue;

    const atom_number_t j = b->other(zatom);

//  cerr << " Atom " << j << " z " << m.atomic_number(j) << " cmp " << z << endl;

    if (m.atomic_number(j) == z)
      return j;
  }

  return INVALID_ATOM_NUMBER;
}

/*
  We are looking for a nitrogen that is adjacent to a carbonyl
*/

static int
attached_to_aromatic_carbonyl (const Molecule & m,
                               const atom_number_t n,
                               const atom_number_t c,
                               const IWStandard_Current_Molecule & current_molecule_data)
{
  const Atom * an = m.atomi(n);

  if (2 != an->ncon())
    return 0;

  atom_number_t carbonyl = an->other(n, 0);
  if (c == carbonyl)
    carbonyl = an->other(n,1);

  if (3 != current_molecule_data.ncon()[carbonyl])
    return 0;

  if (! current_molecule_data.atom_is_aromatic()[carbonyl])
    return 0;

  const atom_number_t oxygen = atom_is_connected_to(m, carbonyl, 8, DOUBLE_BOND);

  if (INVALID_ATOM_NUMBER == oxygen)
    return 0;

  return 1;
}

/*
  Make sure N1 is between C and a carbonyl (if present)
*/

static int
identify_nitrogen_adjacent_to_carbonyl(const Molecule & m,
                                       const atom_number_t c,
                                       atom_number_t & n1,
                                       atom_number_t & n2,
                                       const IWStandard_Current_Molecule & current_molecule_data)
{
  const int x1 = attached_to_aromatic_carbonyl(m, n1, c, current_molecule_data);
  const int x2 = attached_to_aromatic_carbonyl(m, n2, c, current_molecule_data);

//cerr << "x1 " << x1 << " x2 " << x2 << endl;

  if (x1 && x2)
    return 0;

  if (x1)
    return 1;
  else if (x2)
  {
    std::swap(n1, n2);
    return 1;
  }

  return 0;
}

int
Chemical_Standardisation::_do_transform_aromatic_ring_guanidine (Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ncon = current_molecule_data.ncon();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  const Set_of_Atoms & g = current_molecule_data.possible_guanidine();

  const int ng = g.number_elements();

//cerr << "Chemical_Standardisation::_do_transform_aromatic_ring_guanidine: " << ng << " possible guanidines\n";

  if (0 == ng)
    return 0;

  int rc = 0;

  for (int i = ng - 1; i >= 0; --i)
  {
    const atom_number_t c = g[i];

    if (! atom_is_aromatic[c])
      continue;

    atom_number_t n1, n2, doubly_bonded_n;

    if (! collect_attached_nitrogen_atoms(m, c, n1, n2, doubly_bonded_n, current_molecule_data))
      continue;

//  cerr << "Atoms: C " << c << " n1 " << n1 << " n2 " << n2 << " DOUBLY " << doubly_bonded_n << endl;

    if (! atom_is_aromatic[doubly_bonded_n])   // must have double bond in ring
      continue;

    if (2 == ncon[n1] && 1 == ncon[n2])    // make sure N1 is in the ring
      ;
    else if (1 == ncon[n1] && 2 == ncon[n2])
      std::swap(n1, n2);
    else
      continue;

//  cerr << "N1 " << n1 << " aromatic " << atom_is_aromatic[n1] << endl;
    if (! atom_is_aromatic[n1])
      continue;

    if (! identify_nitrogen_adjacent_to_carbonyl(m, c, n1, doubly_bonded_n, current_molecule_data))    // will swap these two so N1 is adjacent to carbonyl
      continue;

//  cerr << "Adjacent to carbonyl " << n1 << " hcount " << m.hcount(n1) << endl;

    if (m.hcount(n1))    // already good
      continue;

#ifdef DEBUG_AROMATIC_RING_GUANIDINE
    cerr << "Before setting bond " << m.smiles() << endl;
#endif
    m.set_bond_type_between_atoms(c, n1, SINGLE_BOND);
    m.set_bond_type_between_atoms(c, doubly_bonded_n, DOUBLE_BOND);    // remember, may have been switched above..
#ifdef DEBUG_AROMATIC_RING_GUANIDINE
    cerr << "After  setting bond " << m.smiles() << endl;
#endif

    current_molecule_data.remove_possible_guanidine(c);

    rc++;
  }

  if (rc)
  {
    _transform_aromatic_guanidine_ring.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " aromatic guanidines\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:arguan";
  }

  return rc;
}

/*
  For Research Records, we often want counterions to be properly charged

  We put + charges on
   Na, K
  Two + on
   Ca
  Negative charges on
   F, Cl, Br, I
*/

int
Chemical_Standardisation::_do_transform_single_atom_ions (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  int matoms = m.natoms();

  int rc = 0;

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 != ncon[i])       // we only process isolated atoms
      continue;

    if (atoms[i]->formal_charge())    // already set, we don't change it
      continue;

    if (atoms[i]->implicit_hydrogens_known() && const_cast<Atom *>(atoms[i])->implicit_hydrogens() > 0)    // don't change it
      continue;

    if (atoms[i]->element()->is_halogen())
    {
      m.set_formal_charge(i, -1);

      rc++;
    }
    else if (! atoms[i]->element()->is_metal())     // we only process metals here
      ;
    else if (3 == z[i] || 11 == z[i] || 19 == z[i])     // Li, Na, K
    {
      m.set_formal_charge(i, 1);

      rc++;
    }
    else if (12 == z[i] || 20 == z[i] || 30 == z[i] || 56 == z[i])   // Mg, Ca, Zn, Ba
    {
      m.set_formal_charge(i, 2);

      rc++;
    }
  }

  if (rc)
  {
    _transform_single_atom_ions.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " isolated metals/halogens\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:isolated";
  }

  return rc;
}

int
Chemical_Standardisation::_do_explicit_hydrogens_last (Molecule & m)
{
  int matoms = m.natoms();

// Since other standardisations may have removed Hydrogen atoms, we count
// the number still present

  atom_number_t last_non_hydrogen_atom = INVALID_ATOM_NUMBER;

  Set_of_Atoms eh;
  for (int i = 0; i < matoms; i++)
  {
    if (1 == m.atomic_number(i))
      eh.add(i);
    else 
      last_non_hydrogen_atom = i;
  }

  if (0 == eh.number_elements())    // none present, that's easy
    return 0;

  if (INVALID_ATOM_NUMBER == last_non_hydrogen_atom)    // Huh, a molecule with just Hydrogen atoms
    return 0;

// Any explicit Hydrogen that is already past the last non-hydrogen atom doesn't need to be removed

  for (int i = eh.number_elements() - 1; i >= 0; i--)
  {
    atom_number_t h = eh[i];

    if (h > last_non_hydrogen_atom)   // doesn't need to be moved around
      eh.remove_item(i);
  }

  int neh = eh.number_elements();

  if (0 == neh)    // order must have already been OK
    return 1;

  for (int i = neh - 1; i >= 0; i--)
  {
    atom_number_t h = eh[i];

    m.move_atom_to_end_of_atom_list(h);
  }

  _explicit_hydrogens_last.extra(neh);

  if (neh)
  {
    if (_verbose)
      cerr << "Moved " << neh << " explicit Hydrogens\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:ehlast";
  }

  return neh;
}


int
Chemical_Standardisation::_process (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  int rc = 0;

//cerr << m.smiles() << ' ' << __LINE__ << " rc " << rc << endl;
  int * atom_already_changed = new_int(m.natoms()); std::unique_ptr<int[]> free_atom_already_changed(atom_already_changed);

  if (_transform_nitro_reverse.active())
    rc += _do_transform_reverse_nitro(m, current_molecule_data);

  if (_transform_azid_reverse.active())
    rc += _do_transform_reverse_azid(m, current_molecule_data);

  if (_transform_nv5_to_charge_separated.active())
    rc += _do_nv5_to_charge_separated(m, current_molecule_data);

  if (_transform_misdrawn_urea.active())
    rc += _do_transform_misdrawn_urea(m, current_molecule_data);

  if (_transform_nitro.active() && current_molecule_data.nplus() && current_molecule_data.ominus())
    rc += _do_transform_nitro(m, current_molecule_data);

  if (_transform_single_atom_ions.active())
    rc += _do_transform_single_atom_ions(m, current_molecule_data);

  if (_transform_nplus_ominus.active() && current_molecule_data.nplus() && current_molecule_data.ominus())
    rc += _do_transform_nplus_ominus(m, current_molecule_data);

  if (_transform_plus_minus.active() && current_molecule_data.nneg() && current_molecule_data.npos())
    rc += _do_transform_plus_minus(m, current_molecule_data);

// I used to have this towards the end, but I think it should come to the front.

  if (_transform_nminus.active() && current_molecule_data.nneg())
    rc += _do_transform_nminus(m, current_molecule_data);

  if (_transform_azid.active() && current_molecule_data.nneg() && current_molecule_data.npos() && current_molecule_data.nitrogens() >= 2)
    rc += _do_transform_azid(m, current_molecule_data);

  if (_transform_n_charge_sep.active() && current_molecule_data.nneg() && current_molecule_data.nplus())
    rc += _do_transform_n_charge_sep(m, current_molecule_data);

  if (_protonate_no.active() && current_molecule_data.ominus())
    rc += _do_protonate_no(m, current_molecule_data);

  if (_protonate_carboxyllic_acids.active() &&(current_molecule_data.ominus() || current_molecule_data.sminus()))
    rc += _do_protonate_carboxyllic_acids(m, current_molecule_data);

// The order of the sulphur acid types is important.
// Make sure sulfonic are done before sulfinic

  if (_protonate_sulfonic_acids.active() &&(current_molecule_data.ominus() || current_molecule_data.splus()) && current_molecule_data.sulphur())
    rc += _do_protonate_sulfonic_acids(m, current_molecule_data);

  if (_protonate_sulfinic_acids.active() &&(current_molecule_data.ominus() || current_molecule_data.splus()) && current_molecule_data.sulphur())
    rc += _do_protonate_sulfinic_acids(m, current_molecule_data);

  if (_protonate_sulfur_acids.active() && current_molecule_data.sminus())
    rc += _do_protonate_sulfur_acids(m, current_molecule_data);

  if (_protonate_phosphorous_acids.active() && current_molecule_data.phosphorus() &&(current_molecule_data.ominus() || current_molecule_data.sminus()))
    rc += _do_protonate_phosphorous_acids(m, current_molecule_data);

  if (_transform_splus_cminus.active() && current_molecule_data.nneg() && current_molecule_data.npos() && current_molecule_data.sulphur() && current_molecule_data.nneg())
    rc += _do_transform_splus_cminus(m, current_molecule_data);

  if (_transform_guanidine.active() && current_molecule_data.possible_guanidine().number_elements())
    rc += _do_transform_guanidine(m, current_molecule_data);

//cerr << m.smiles() << ' ' << __LINE__ << " rc " << rc << endl;

  if (_transform_guanidine_ring.active() && current_molecule_data.possible_guanidine().number_elements())
    rc += _do_transform_ring_guanidine(m, current_molecule_data);

  if (_transform_tetrazole.active())
    rc += _do_tetrazole(m, current_molecule_data);

// Explicit Hydrogens may get removed from Amines

  Set_of_Atoms atoms_to_be_removed;

  if (_transform_amines.active() && current_molecule_data.nplus())
    rc += _do_transform_amines(m, atoms_to_be_removed, current_molecule_data);

// Because some molecules 
// have both a lactam/lactim AND an aromatic nitrogen issue, we must do the 
// lactam/lactim thing first

  if (_transform_pyrazolone.active())
    rc += _do_transform_pyrazolone(m, atom_already_changed, current_molecule_data);

//cerr << "after _do_transform_pyrazolone " << rc << endl;

  if (_transform_lactim_lactam.active() && current_molecule_data.possible_lactam().number_elements() > 0)
    rc += _do_transform_non_ring_lactim(m, atom_already_changed, current_molecule_data);

//cerr << " after _do_transform_non_ring_lactim " << rc << endl;

  if (_transform_lactim_lactam_ring.active() && current_molecule_data.possible_lactam().number_elements() > 0)
    rc += _do_transform_ring_lactim(m, atom_already_changed, current_molecule_data);

//cerr << " after _do_transform_ring_lactim " << rc << endl;

//#define DEBUG_LACTAM_LACTIM
#ifdef DEBUG_LACTAM_LACTIM
  cerr << "After lactam/lactim " << m.smiles() << endl;
#endif

  if (_transform_imidazole.active())
    rc += _do_imidazole(m, current_molecule_data);

//cerr << "Before _do_pyrazole rc " << rc << endl;
  if (_transform_pyrazole.active())
    rc += _do_pyrazole(m, atom_already_changed, current_molecule_data);
//cerr << "After _do_pyrazole rc " << rc << endl;

//cerr << "After pyrazole " << m.unique_smiles() << endl;

  if (_transform_triazole.active())
    rc += _do_triazole(m, current_molecule_data);

  if (_transform_isoxazole.active())
    rc += _do_isoxazole(m, current_molecule_data);

  if (_from_mrk_standardisations.active())
    rc += _do_from_mrk_standardisations(m, current_molecule_data);

// Make sure we do O- after everything else

  if (_transform_ominus.active() &&(current_molecule_data.ominus() || current_molecule_data.sminus()))
    rc += _do_transform_ominus(m, current_molecule_data);

  if (_transform_back_to_nplus_nminus.active() && current_molecule_data.nitrogens() > 1)
    rc += _do_transform_back_to_nplus_nminus(m, current_molecule_data);

  if (_transform_to_charge_separated_azid.active() && current_molecule_data.nitrogens() > 2)
    rc += _do_transform_azid_to_charge_separated(m, current_molecule_data);

  if (_transform_obvious_implicit_hydrogen_errors.active())
    rc += _do_transform_implicit_hydrogen_known_errors(m, current_molecule_data);

  if (_transform_aromatic_guanidine_ring.active() && current_molecule_data.possible_guanidine().number_elements())
    rc += _do_transform_aromatic_ring_guanidine(m, current_molecule_data);

  if (_transform_amino_thiazole.active())
    rc += _do_amino_thiazole(m, atom_already_changed, current_molecule_data);

  assert (current_molecule_data.npos() >= 0);
  assert (current_molecule_data.nneg() >= 0);
  assert (current_molecule_data.nplus() >= 0);
  assert (current_molecule_data.ominus() >= 0);
  assert (current_molecule_data.sminus() >= 0);

  if (atoms_to_be_removed.number_elements())
    m.remove_atoms(atoms_to_be_removed);

// Since this changes the order of the atoms we don't use any of the arrays

//cerr << m.smiles() << ' ' << __LINE__ << " rc " << rc << endl;

  if (current_molecule_data.explicit_hydrogen_count() && _explicit_hydrogens_last.active())
    rc += _do_explicit_hydrogens_last(m);

//cerr << m.smiles() << ' ' << __LINE__ << " rc " << rc << endl;
  return rc;
}

int
Chemical_Standardisation::process (Molecule & m) 
{
//cerr << "Chemical_Standardisation::process:active " << _active << endl;

  if (0 == _active)
    return 1;

  int asave = global_aromaticity_type();

  if (Daylight != asave)
  {
    set_global_aromaticity_type(Daylight);
    m.compute_aromaticity();
  }

  const int rc = _process(m);


  set_global_aromaticity_type(asave);

//cerr << "AFter Chemical_Standardisation::process\n";
//m.debug_print(cerr);

  return rc;
}

/*
  There are some fundamental difficulties here. Some transformations will collide with others.
  So, there is order dependency. Generally, this is not a problem, but there are amsll number of
  molecules for which things go wrong. For now, let's just ignore the small number of cases
  that fail.
*/

int
Chemical_Standardisation::_process(Molecule & m)
{
  _molecules_processed++;

  int valence_before_processing = 1;

  if (_check_valence_before_and_after)
    valence_before_processing = m.valence_ok();

// Removing atoms can mess up anything which comes after, so make 
// sure we do that before anything else

  if (_append_string_depending_on_what_changed)
    _append_to_changed_molecules.resize_keep_storage(0);

  int rc = 0;

  if (_remove_hydrogens.active())
    rc += _do_remove_hydrogens(m);

  if (0 == m.natoms())
    return 0;

  if (_transform_covalent_metals.active())
    _do_transform_covalent_metals(m);

  IWStandard_Current_Molecule current_molecule_data;

  if (! current_molecule_data.initialise(m))
  {
    cerr << "Chemical_Standardisation::process:could not initialise '" << m.name() << "'\n";
    return 0;
  }

  if (! _processing_needed(current_molecule_data))
    return 0;

  rc += _process(m, current_molecule_data);

  if (rc)
  {
    _molecules_changed++;
    m.recompute_implicit_hydrogens();

    if (_append_to_changed_molecules.length())
      m.append_to_name(_append_to_changed_molecules);

    if (_check_valence_before_and_after)
    {
      int valence_check_after = m.valence_ok();

      if (valence_check_after == valence_before_processing)
        ;
      else if (valence_before_processing && ! valence_check_after)
        cerr << "Chemical_Standardisation::process:invalid valence introduced '" << m.name() << "'\n";
    }
  }

  return rc;
}

int
Chemical_Standardisation::report(std::ostream & os) const
{
  os << "Report on chemical standardisation object\n";

  if (0 == _molecules_processed)
    return os.good();

  if (_transform_amines.active())
  {
    os << "  Amines ";
    _transform_amines.report(os);
  }
  if (_transform_nitro.active())
  {
    os << "  Nitro groups ";
    _transform_nitro.report(os);
  }

  if (_transform_nplus_ominus.active())
  {
    os << "  N+O- ";
    _transform_nplus_ominus.report(os);
  }

  if (_transform_nv5_to_charge_separated.active())
  {
    os << "  Nv5 ";
    _transform_nv5_to_charge_separated.report(os);
  }

  if (_transform_plus_minus.active())
  {
    os << "  Plus Minus ";
    _transform_plus_minus.report(os);
  }

  if (_transform_covalent_metals.active())
  {
    os << "  Covalent metals ";
    _transform_covalent_metals.report(os);
  }

  if (_transform_n_charge_sep.active())
  {
    os << "  Nitrogen charge separated ";
    _transform_n_charge_sep.report(os);
  }

  if (_remove_hydrogens.active())
  {
    os << "  Remove hydrogens ";
    _remove_hydrogens.report(os);
  }

  if (_protonate_carboxyllic_acids.active())
  {
    os << "  Protonate carboxyllic acids ";
    _protonate_carboxyllic_acids.report(os);
  }

  if (_protonate_no.active())
  {
    os << "  Protonate NO ";
    _protonate_no.report(os);
  }

  if (_protonate_sulfonic_acids.active())
  {
    os << "  Protonate sulfonic acids ";
    _protonate_sulfonic_acids.report(os);
  }

  if (_protonate_sulfinic_acids.active())
  {
    os << "  Protonate sulfinic acids ";
    _protonate_sulfinic_acids.report(os);
  }

  if (_transform_splus_cminus.active())
  {
    os << "  [S+]-[C-] ";
    _transform_splus_cminus.report(os);
  }

  if (_transform_ominus.active())
  {
    os << "  Protonate O- ";
    _transform_ominus.report(os);
  }

  if (_transform_nminus.active())
  {
    os << "  Protonate N- ";
    _transform_nminus.report(os);
  }

  if (_protonate_sulfinic_acids.active())
  {
    os << "  Protonate [S-]-[C,S,P]=[C,S] ";
    _protonate_sulfinic_acids.report(os);
  }

  if (_protonate_phosphorous_acids.active())
  {
    os << "  Protonate [O,S]-P=[O,S] ";
    _protonate_phosphorous_acids.report(os);
  }

  if (_explicit_hydrogens_last.active())
  {
    os << "  Explicit Hydrogens last ";
    _explicit_hydrogens_last.report(os);
  }

  if (_transform_tetrazole.active())
  {
    os << "  Tetrazole ";
    _transform_tetrazole.report(os);
  }

  if (_transform_triazole.active())
  {
    os << "  Triazole ";
    _transform_triazole.report(os);
  }

  if (_transform_isoxazole.active())
  {
    os << "  Isoxazole ";
    _transform_isoxazole.report(os);
  }

  if (_transform_imidazole.active())
  {
    os << "  Imidazole ";
    _transform_imidazole.report(os);
  }

  if (_transform_pyrazole.active())
  {
    os << "  Pyrazole ";
    _transform_pyrazole.report(os);
  }

  if (_transform_guanidine.active())
  {
    os << "  Guanidine ";
    _transform_guanidine.report(os);
  }

  if (_transform_guanidine_ring.active())
  {
    os << "  Guanidine Ring ";
    _transform_guanidine_ring.report(os);
  }

  if (_transform_lactim_lactam.active())
  {
    os << "  Lactim_lactam ";
    _transform_lactim_lactam.report(os);
  }

  if (_transform_lactim_lactam_ring.active())
  {
    os << "  Lactim_lactam ring ";
    _transform_lactim_lactam_ring.report(os);
  }

  if (_transform_azid.active())
  {
    os << "  Azid ";
    _transform_azid.report(os);
  }

  if (_transform_misdrawn_urea.active())
  {
    os << "  MSDUR ";
    _transform_misdrawn_urea.report(os);
  }

  if (_transform_back_to_nplus_nminus.active())
  {
    os << "  =N=N- to =[N+]-[N-] ";
    _transform_back_to_nplus_nminus.report(os);
  }

  if (_transform_amino_thiazole.active())
  {
    os << "  amino thiazole ";
    _transform_amino_thiazole.report(os);
  }

  return os.good();
}

/*int
Chemical_Standardisation::deactivate (const const_IWSubstring & d)
{
  if (CS_NITRO  == d)
    _transform_nitro.deactivate();
  else if (CS_NpOm == d)
    _transform_nplus_ominus.deactivate();
  else if (CS_NpNm == d)
      _transform_nplus_ominus.deactivate();
  else if (CS_SpCm == d)
    _transform_splus_cminus.deactivate();
  else if (CS_ALLpm == d)
    _transform_plus_minus.deactivate();
  else if (CS_XH == d)
    _remove_hydrogens.deactivate();
  else if (CS_AMINE == d)
    _transform_amines.deactivate();
  else if (CS_Om == d)
    _transform_ominus.deactivate();
  else if (CS_Nm == d)
    _transform_nminus.deactivate();
  else
  {
    cerr << "Chemical_Standardisation::deactivate: unrecognised directive '" << d << "'\n";
    return 0;
  }

  return 1;
}*/

int
Chemical_Standardisation::activate_all_except_hydrogen_removal()
{
  activate_all();

  _remove_hydrogens.deactivate();

  return 1;
}

void
Chemical_Standardisation::deactivate_lactim_lactam ()
{
  _transform_lactim_lactam_ring.deactivate();
  _transform_lactim_lactam.deactivate();

  return;
}

/*
  Change N#N=N to [N-]=[N+]=N
*/

int
Chemical_Standardisation::_do_transform_reverse_azid (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    const Bond * b = ai->item(0);

    if (! b->is_triple_bond())
      continue;

    atom_number_t n2 = b->other(i);

    if (7 != z[n2])
      continue;

    if (2 != ncon[n2])
      continue;

    const Atom * an2 = atoms[n2];

    if (5 != an2->nbonds())
      continue;

    for (int j = 0; j < 2; j++)
    {
      const Bond * b = an2->item(j);

      if (! b->is_double_bond())
        continue;

      atom_number_t n3 = b->other(n2);

      if (7 != z[n3])
        continue;

      m.set_bond_type_between_atoms(i, n2, DOUBLE_BOND);
      m.set_formal_charge(i, -1);
      m.set_formal_charge(n2, 1);

      rc++;
    }
  }

  if (rc)
  {
    _transform_azid_reverse.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " azid groups to charge separated\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:azid+-";
  }

  return rc;
}

/*
  Change [N-]=[N+]=N to N#N=N
*/

int
Chemical_Standardisation::_do_transform_azid(Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (-1 != ai->formal_charge())
      continue;

    const Bond * b = ai->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t n2 = b->other(i);

    if (7 != z[n2])
      continue;

    if (2 != ncon[n2])
      continue;

    const Atom * an2 = atoms[n2];

    if (1 != an2->formal_charge())
      continue;

    if (4 != an2->nbonds())
      continue;

    for (int j = 0; j < 2; j++)
    {
      const Bond * b = an2->item(j);

      if (! b->is_double_bond())    // how could that happen?
        continue;

      atom_number_t n3 = b->other(n2);
      if (i == n3)
        continue;

      if (7 == z[n3] || 6 == z[n3])   // not sure what else could be there...
        ;
      else
        continue;

      if (-1 == atoms[n3]->formal_charge())    // pathological case of [N+](=[N-])=[N-] 16606
        continue;

      m.set_bond_type_between_atoms(i, n2, TRIPLE_BOND);
      m.set_formal_charge(i, 0);
      m.set_formal_charge(n2, 0);
      current_molecule_data.change_nplus(-1);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_nneg(-1);

      rc++;

      break;
    }
  }

  if (rc)
  {
    _transform_azid.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " charge separated azid groups to neutral\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:azid";
  }

  return rc;
}

int
Chemical_Standardisation::_do_transform_reverse_nitro (Molecule & m,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * o1 = atoms[i];

    const Bond * b1 = o1->item(0);

    if (! b1->is_double_bond())
      continue;

    atom_number_t n = b1->other(i);

    if (7 != z[n])
      continue;
     
    if (3 != ncon[n])
      continue;

    const Atom * na = atoms[n];

    if (5 != na->nbonds())
      continue;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = na->item(j);

      if (! b->is_double_bond())
        continue;

      atom_number_t o2 = b->other(n);

      if (i == o2)
        continue;

      if (8 != z[o2])
        continue;

      m.set_bond_type_between_atoms(i, n, SINGLE_BOND);
      m.set_formal_charge(i, -1);
      m.set_formal_charge(n, +1);

      rc++;
    }
  }

  if (rc)
  {
    _transform_nitro_reverse.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " nitro groups to charge separated\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:NO2+-";
  }

  return rc;
}

int
Chemical_Standardisation::_do_nv5_to_charge_separated(Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (ncon[i] < 2)
      continue;

    const Atom * a = atoms[i];

    if (5 != a->nbonds())
      continue;

    if (0 != a->formal_charge())
      continue;

    atom_number_t doubly_bonded_singly_connected = INVALID_ATOM_NUMBER;
    atom_number_t triply_connected_n = INVALID_ATOM_NUMBER;
    bond_type_t bond_to_be_placed = SINGLE_BOND;
    atom_number_t double_bonded_2_connected_n = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = a->item(j);

      if (b->is_single_bond())
        continue;

      atom_number_t k = b->other(i);

      if (1 == ncon[k] && 8 == z[k] && b->is_double_bond())
        doubly_bonded_singly_connected = k;
      else if (1 == ncon[k] && 7 == z[k] && b->is_triple_bond())
      {
        triply_connected_n = k;
        bond_to_be_placed = DOUBLE_BOND;
      }
      else if (2 == ncon[k] && 7 == z[k] && b->is_double_bond())
        double_bonded_2_connected_n = k;
    }

    if (INVALID_ATOM_NUMBER != doubly_bonded_singly_connected)
      ;
    else if (INVALID_ATOM_NUMBER != triply_connected_n)
      doubly_bonded_singly_connected = triply_connected_n;
    else if (INVALID_ATOM_NUMBER != double_bonded_2_connected_n)
      doubly_bonded_singly_connected = double_bonded_2_connected_n;
    else
      continue;

    m.set_formal_charge(i, 1);
    m.set_formal_charge(doubly_bonded_singly_connected, -1);
    m.set_bond_type_between_atoms(i, doubly_bonded_singly_connected, bond_to_be_placed);
    rc++;
  }

  if (rc)
  {
    _transform_nv5_to_charge_separated.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " 5 valent Nitrogens to charge separated form\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:RNV5";
  }

  return rc;
}

int
Chemical_Standardisation::_do_from_mrk_standardisations (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int matoms = m.natoms();

  if (0 == current_molecule_data.nneg() && 0 == current_molecule_data.npos())
    return 0;

  int rc = 0;

  if (current_molecule_data.phosphorus() && current_molecule_data.nneg() > 1 && current_molecule_data.npos())   //  [O-]-[P+]-[O-]
  {
    for (int i = 0; i < matoms; i++)
    {
      if (15 != z[i])
        continue;

      if (4 != ncon[i])
        continue;

      const Atom * ai = atoms[i];

      if (ai->formal_charge() < 1)
        continue;

      atom_number_t n1;
      atom_number_t n2;

      if (! two_negatively_charged_connections(i, *ai, atoms, n1, n2))
        continue;

      m.set_formal_charge(i, 0);
      m.set_formal_charge(n1, 0);
      m.set_formal_charge(n2, 0);
      m.set_bond_type_between_atoms(i, n1, DOUBLE_BOND);
      m.set_bond_type_between_atoms(i, n2, DOUBLE_BOND);
      current_molecule_data.change_nneg(-2);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_phosphorus(-1);
      rc++;
    } 
  }

  if ((current_molecule_data.sulphur() || current_molecule_data.phosphorus()) && current_molecule_data.nneg() > 1 && current_molecule_data.npos())    // try to fix [*-]-[S,++]-[*-]
  {
    for (int i = 0; i < matoms; i++)
    {
      if (16 != z[i])
        continue;

      if (4 != ncon[i])
        continue;

      const Atom * ai = atoms[i];

      if (2 != ai->formal_charge())
        continue;

      atom_number_t n1;
      atom_number_t n2;

      if (! two_negatively_charged_connections(i, *ai, atoms, n1, n2))
        continue;

      m.set_formal_charge(i, 0);
      m.set_formal_charge(n1, 0);
      m.set_formal_charge(n2, 0);
      m.set_bond_type_between_atoms(i, n1, DOUBLE_BOND);
      m.set_bond_type_between_atoms(i, n2, DOUBLE_BOND);
      current_molecule_data.change_nneg(-2);
      current_molecule_data.change_npos(-1);
      current_molecule_data.change_sulphur(-1);
      rc++;
    }
  }

  if (rc)
  {
    _from_mrk_standardisations.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " moities that had been changed in ->MRK form\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:MRK";
  }

  return rc;
}

/*
  Make sure that the Nitrogen with the Hydrogen is adjacent to the carbon
*/

int
Chemical_Standardisation::_do_tetrazole (Molecule & m,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_size = current_molecule_data.ring_size();
  const int * ring_is_fused = current_molecule_data.ring_is_fused();

  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    if (4 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    if (5 != ring_size[i])
      continue;

    if (ring_is_fused[i])
      continue;

    const Set_of_Atoms * r = current_molecule_data.ringi(i);

    rc += _do_tetrazole(m, *r, current_molecule_data);
  }

  if (rc)
  {
    _transform_tetrazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " tetrazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:tetrazole";
  }

  return rc;
}

int
Chemical_Standardisation::_do_triazole (Molecule & m,
                                        const Set_of_Atoms & r,
                                        const Set_of_Atoms * is_fused,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();

//#define DEBUG_DO_TRIAZOLE
#ifdef DEBUG_DO_TRIAZOLE
  cerr << "Checking possible triazole " << r << endl;
#endif

  int n1_index_in_ring = -1;
  int n2_index_in_ring = -1;
  int n3_index_in_ring = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    if (7 != z[j])
      continue;

    if (2 != ncon[j])    // If there is a 3 connected Nitrogen, we cannot change the bonding in the ring
      return 0;

    if (n1_index_in_ring < 0)   // there must be exactly three nitrogens in the ring
      n1_index_in_ring = i;
    else if (n2_index_in_ring < 0)
      n2_index_in_ring = i;
    else if (n3_index_in_ring < 0)
      n3_index_in_ring = i;
    else
      return 0;
  }

#ifdef DEBUG_DO_TRIAZOLE
  cerr << " Indices N1 " << n1_index_in_ring << " (atom " << r[n1_index_in_ring] << ") N2 " << n2_index_in_ring << " (atom " << r[n2_index_in_ring] << ") and N3 " << n3_index_in_ring << " (atom " << r[n3_index_in_ring] << ")\n";
#endif

  if (n3_index_in_ring < 0)
    return 0;

// We must decide if we have a 1,2,4 or 1,2,3 triazole

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring && 2 == n3_index_in_ring)
    return _do_123_triazole(m, r, 0, 1, 2, 3, 4, current_molecule_data, is_fused);

  if (0 == n1_index_in_ring && 3 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_123_triazole(m, r, 3, 4, 0, 1, 2, current_molecule_data, is_fused);

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_123_triazole(m, r, 4, 0, 1, 2, 3, current_molecule_data, is_fused);

  if (1 == n1_index_in_ring && 2 == n2_index_in_ring && 3 == n3_index_in_ring)
    return _do_123_triazole(m, r, 1, 2, 3, 4, 0, current_molecule_data, is_fused);

  if (2 == n1_index_in_ring && 3 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_123_triazole(m, r, 2, 3, 4, 0, 1, current_molecule_data, is_fused);

// Several more 134 types

  if (0 == n1_index_in_ring && 2 == n2_index_in_ring && 3 == n3_index_in_ring)
    return _do_134_triazole(m, r, 0, 1, 2, 3, 4, current_molecule_data);

  if (1 == n1_index_in_ring && 3 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_134_triazole(m, r, 1, 2, 3, 4, 0, current_molecule_data);

  if (1 == n1_index_in_ring && 2 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_134_triazole(m, r, 4, 3, 2, 1, 0, current_molecule_data);

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring && 3 == n3_index_in_ring)
    return _do_134_triazole(m, r, 3, 2, 1, 0, 4, current_molecule_data);

  if (0 == n1_index_in_ring && 2 == n2_index_in_ring && 4 == n3_index_in_ring)
    return _do_134_triazole(m, r, 2, 3, 4, 0, 1, current_molecule_data);

  cerr << "Unrecognised triazole form n1 = " << n1_index_in_ring << " n2 = " << n2_index_in_ring << " n3 = " << n3_index_in_ring << endl;
  return 0;
}

static int
acyl_group_attached (atom_number_t zatom,
                     const Atom * const * atoms)
{
  const Atom * a = atoms[zatom];

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other(zatom);

    return 8 == atoms[j]->atomic_number();
  }

  return 0;
}

static int
is_cf3 (atom_number_t zatom,
        const Atom * const * atoms)
{
  const Atom * a = atoms[zatom];

  assert (4 == a->ncon());

  int nfluorine = 0;

  for (int i = 0; i < 4; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (9 != atoms[j]->atomic_number())
      continue;

    nfluorine++;

    if (3 == nfluorine)
      return 1;
  }

  return 0;
}

/*
  Sept 2007. Google search for "electron withdrawing groups" yielded.

  http://www.mhhe.com/physsci/chemistry/carey/student/olc/graphics/carey04oc/ref/ch12substituenteffects.html


*/

static double
compute_electron_donating_power (Molecule & m,
                                 atom_number_t astart,
                                 atom_number_t avoid,
                                 const Atom * const * atoms)
{
#ifdef DEBUG_COMPUTE_ELECTRON_DONATING_POWER
  cerr << "Computing electron donating tendency for atom " << astart << " avoid " << avoid << endl;
#endif

  const Atom * a = atoms[astart];

  atomic_number_t z = a->atomic_number();

  int acon = a->ncon();

  if (1 == acon)
  {
    if (6 == z)
      return 1.3;

    if (7 == z)
      return 3.3;

    if (8 == z)
      return 3.5;

    if (9 == z)
      return -1.4;
    else if (17 == z)
      return -1.3;
    else if (35 == z)
      return -1.2;
    else if (53 == z)
      return -1.1;

    if (16 == z)
      return 3.25;

    return 0.0;   // hmm, what is this?
  }

// We must differentiate NR2 and NO2

  if (7 == z && 3 == acon)    // NR2 or NO2
  {
    if (3 == a->nbonds())   // NR2
      return 3.4;
    else if (5 == a->nbonds())    // presumably NO2
      return -3.6;
  }

  if (7 == z && 4 == acon)
    return -3.5;

  if (16 == z)
  {
    if (4 == acon)   // most likely SO2
      return -3.3;
    else
      return 2.0;   // maybe a thioether?
  }

  if (6 == z)
  {
    if (2 == acon && 4 == a->nbonds())   // cyano or acetylene. Actually reference lists cyano but not acetylene...
      return -3.2;

    if (4 == acon)
    {
      if (is_cf3(astart, atoms))
        return -3.0;
      else
        return 1.3;
    }
  }

// Looking at what is left, it seems that if there is an =O group attached
// we are withdrawing.
// IAW makes up a bunch of other heuristics

  if (acon == a->nbonds())   // saturated
    return static_cast<double>(m.attached_heteroatom_count(astart)) * -0.20;   // IAW
  else if (acyl_group_attached(astart, atoms))
    return -2.0;

// the unsaturated case, which may include aromatic. Cannot compute aromaticity

  int ahc = m.attached_heteroatom_count(astart);

  if (m.multiple_bond_to_heteroatom(astart))   // most likely a Nitrogen
    return 0.57;

  if (ahc > 0)
    return -0.21 * static_cast<double>(ahc);

  return 0.11 * static_cast<double>(acon);

  return 0.0;
}

/*
  The caller must present things with this numbering

       N2
      /  \
     /    \
    N1    N3
    |      |
    |      |
    C5----C4

  We try to put the Hydrogen on N1
*/

static int
place_123_triazole_bonds (Molecule & m,
                          atom_number_t a1,
                          atom_number_t a2,
                          atom_number_t a3,
                          atom_number_t a4,
                          atom_number_t a5,
                          const Atom * const * atoms,
                          const Set_of_Atoms & r,
                          const Set_of_Atoms * is_fused)
{
  if (2 == atoms[a1]->nbonds())   // already set up properly
    return 0;

//#define DEBUG_PLACE_123_TRIAZOLE_BONDS
#ifdef DEBUG_PLACE_123_TRIAZOLE_BONDS
  cerr << "Setting triazole " << m.smiles() << ", fused " << is_fused << endl;
#endif

  m.set_bond_type_between_atoms(a1, a2, SINGLE_BOND);
  m.set_bond_type_between_atoms(a2, a3, DOUBLE_BOND);
  m.set_bond_type_between_atoms(a3, a4, SINGLE_BOND);
  m.set_bond_type_between_atoms(a4, a5, DOUBLE_BOND);
  m.set_bond_type_between_atoms(a5, a1, SINGLE_BOND);

  if (NULL == is_fused)
    return 1;

#ifdef DEBUG_PLACE_123_TRIAZOLE_BONDS
  cerr << "After transform " << m.smiles() << ", bonds " << m.nbonds(a4) << " and " << m.nbonds(a5) << endl;
#endif

// Now this gets messy. We do not know which of the atoms above are the fused carbon atoms

  atom_number_t c1 = INVALID_ATOM_NUMBER;
  atom_number_t c2 = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 5; ++i)
  {
    const atom_number_t j = r[i];
    if (6 != m.atomic_number(j))
      continue;

    if (INVALID_ATOM_NUMBER == c1)
      c1 = j;
    else
      c2 = j;
  }

  if (INVALID_ATOM_NUMBER == c2)   // should not happen
    return 0;

  int ndx1 = is_fused->index(c1);
  int ndx2 = is_fused->index(c2);

// We need to switch the bonds in the attached ring - single bonds attached to both a4 and a5

#ifdef DEBUG_PLACE_123_TRIAZOLE_BONDS
  cerr << "Carbon atoms are " << c1 << " index " << ndx1 << " atom " << c2 << " index " << ndx2 << endl;
  cerr << "Ring " << *is_fused <<endl;
#endif

  if (ndx1 < 0 || ndx2 < 0)   // should not happen
    return 0;

// We need to make sure we traverse the ring the long way around. Make ndx2 is one larger than ndx1

  if (ndx1 + 1 == ndx2)
    ;
  else if (ndx2 + 1 == ndx1)
    std::swap(ndx1, ndx2);
  else if (ndx1 < ndx2)    // change 0 5 to 5 0
    std::swap(ndx1, ndx2);

#ifdef DEBUG_PLACE_123_TRIAZOLE_BONDS
  cerr << "Before setting bonds " << m.smiles() << " ring indices " << ndx1 << " and " << ndx2 << endl;
#endif

  m.set_bond_type_between_atoms(is_fused->item(ndx2),       is_fused->item((ndx2+1)%6), SINGLE_BOND);
  m.set_bond_type_between_atoms(is_fused->item((ndx2+1)%6), is_fused->item((ndx2+2)%6), DOUBLE_BOND);
  m.set_bond_type_between_atoms(is_fused->item((ndx2+2)%6), is_fused->item((ndx2+3)%6), SINGLE_BOND);
  m.set_bond_type_between_atoms(is_fused->item((ndx2+3)%6), is_fused->item((ndx2+4)%6), DOUBLE_BOND);
  m.set_bond_type_between_atoms(is_fused->item((ndx2+4)%6), is_fused->item((ndx2+5)%6), SINGLE_BOND);

  return 1;
}

static atom_number_t
identify_extra_ring_atom (const Atom * a, 
                          atom_number_t zatom,
                          atom_number_t avoid1,
                          atom_number_t avoid2)
{
  if (a->ncon() < 3)
    return INVALID_ATOM_NUMBER;

  for (int i = 0; i < 3; i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (j == avoid1)
      continue;

    if (j == avoid2)
      continue;

    return j;
  }

  return INVALID_ATOM_NUMBER;   // not sure how this could happen
}

/*
  The caller must present things with this numbering

       N2
      /  \
     /    \
    N1    N3
    |      |
    |      |
    C5----C4

  We try to put the Hydrogen on N1 or N3.

  BUT, we run into problems with uniqueness. Therefore, if
  we cannot resolve N1 and N3, we put the Hydrogen on N2
*/

int
Chemical_Standardisation::_do_123_triazole (Molecule & m,
                                        const Set_of_Atoms & r,
                                        int n1_index_in_ring, 
                                        int n2_index_in_ring,
                                        int n3_index_in_ring,
                                        int c4_index_in_ring,
                                        int c5_index_in_ring,
                                        IWStandard_Current_Molecule & current_molecule_data,
                                        const Set_of_Atoms * is_fused) const
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  atom_number_t n1 = r[n1_index_in_ring];
  atom_number_t n2 = r[n2_index_in_ring];

//#define DEBUG_123_TRIAZOLE
#ifdef DEBUG_123_TRIAZOLE
  cerr << "Checking possible 123 triazole " << r << endl;
  cerr << "n1 " << n1 << " n2 " << n2 << " n3 " << n3 << " c4 " << c4 << " c5 " << c5 << endl;
  cerr << "nbonds " << atoms[n2]->nbonds() << " ncon " << ncon[n2] << endl;
#endif

  if (3 == ncon[n2])   // rare, cannot change anything
    return 0;

  atom_number_t n3 = r[n3_index_in_ring];

  if (3 == ncon[n1] || 3 == ncon[n3])   // cannot change anything
    return 0;

  atom_number_t c4 = r[c4_index_in_ring];
  atom_number_t c5 = r[c5_index_in_ring];

  if (6 != z[c4])
    return 0;

  if (6 != z[c5])
    return 0;

// do we have unsubstituted 123-triazole

  if (2 == ncon[c4] && 2 == ncon[c5])   // 123-triazole
    return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms, r, is_fused);

// Now we need to decide whether to put the H on N1 or N3.
// We know that at least one of C4 and C5 are substituted, 

  atom_number_t cc4 = identify_extra_ring_atom(atoms[c4], c4, n3, c5);
  atom_number_t cc5 = identify_extra_ring_atom(atoms[c5], c5, n1, c4);

  if (INVALID_ATOM_NUMBER == cc4)   // substituted at c5
  {
    double electron_donating_power = compute_electron_donating_power(m, cc5, c5, atoms);
    if (electron_donating_power > 0.0)   // put H on N1
      return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms, r, is_fused);
    else    // put H on N3
      return place_123_triazole_bonds(m, n3, c4, c5, n1, n2, atoms, r, is_fused);
  }

  if (INVALID_ATOM_NUMBER == cc5)   // sustituted at c4
  {
    double electron_donating_power = compute_electron_donating_power(m, cc4, c4, atoms);
    if (electron_donating_power > 0)   // put H on N3
      return place_123_triazole_bonds(m, n3, c4, c5, n1, n2, atoms, r, is_fused);
    else
      return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms, r, is_fused);
  }

// Substituted at both positions

  double ed4 = compute_electron_donating_power(m, cc4, c4, atoms);
  double ed5 = compute_electron_donating_power(m, cc5, c5, atoms);

#ifdef DEBUG_123_TRIAZOLE
  cerr << "Atoms C4 " << c4 << " and C5 " << c5 << endl;
  cerr << "ed4 " << ed4 << " ed5 " << ed5 << endl;
#endif

// Wierd compiler but with gcc-4.0.2 on Linux. The comparisons failed!!!!
// C1(=CC=C(F)C=C1)C1=C(N=NN1)C1=CC=NC=C1 PBCHM22725982
// Restructure code so equality is checked first. Bizzare stuff

  if (fabs(ed4 - ed5) < 1.0e-09)
    return place_123_triazole_bonds(m, n2, n3, c4, c5, n1, atoms, r, is_fused);
  else if (ed4 < ed5)   // put Hydrogen on N1
    return place_123_triazole_bonds(m, n1, n2, n3, c4, c5, atoms, r, is_fused);
  else
    return place_123_triazole_bonds(m, n3, c4, c5, n1, n2, atoms, r, is_fused);
}

/*
  The caller must present things with this numbering

       N1
      /  \
     /    \
    C5    C2
    |      |
    |      |
    N4----N3

  We put double bonds between C2=N3 and N4=C5
*/

int
Chemical_Standardisation::_do_134_triazole (Molecule & m,
                                        const Set_of_Atoms & r,
                                        int n1_index_in_ring, 
                                        int c2_index_in_ring,
                                        int n3_index_in_ring,
                                        int n4_index_in_ring,
                                        int c5_index_in_ring,
                                        IWStandard_Current_Molecule & current_molecule_data) const
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

//#define DEBUG_134_TRIAZOLE
#ifdef DEBUG_134_TRIAZOLE
  cerr << "Checking 134 triazole " << r[n1_index_in_ring] << ' ' << r[c2_index_in_ring] << ' ' << r[n3_index_in_ring] << ' ' <<  r[n4_index_in_ring] << ',' << r[c5_index_in_ring] << endl;
#endif

  atom_number_t n1 = r[n1_index_in_ring];

#ifdef DEBUG_134_TRIAZOLE
  cerr << "What is n1? " << n1 << " bonds " << atoms[n1]->nbonds() << " ncon " << ncon[n1] << endl;
#endif

  if (atoms[n1]->nbonds() == ncon[n1])   // already single bonds at N1
    return 0;

  atom_number_t c2 = r[c2_index_in_ring];
  atom_number_t n3 = r[n3_index_in_ring];
  atom_number_t n4 = r[n4_index_in_ring];
  atom_number_t c5 = r[c5_index_in_ring];

#ifdef DEBUG_134_TRIAZOLE
  cerr << "Carbons " << c2 << ' ' << z[c2] << " and " << c5 << ' ' << z[c5] << endl;
#endif

  if (6 != z[c2])
    return 0;

  if (6 != z[c5])
    return 0;

// Beware N1=C(N)NNC1=S PBCHM2723869

  if (2 == atoms[n3]->nbonds() && 2 == atoms[n4]->nbonds() && 2 == ncon[n1] && 3 == atoms[n1]->nbonds())
  {
    const Bond * b15 = m.bond_between_atoms(n1, c5);

    if (b15->is_double_bond())
    {
      m.set_bond_type_between_atoms(n1, c5, SINGLE_BOND);
      m.set_bond_type_between_atoms(c5, n4, DOUBLE_BOND);
    }
    else 
    {
      const Bond * b12 = m.bond_between_atoms(n1, c2);
      if (b12->is_double_bond())
      {
        m.set_bond_type_between_atoms(n1, c2, SINGLE_BOND);
        m.set_bond_type_between_atoms(c2, n3, DOUBLE_BOND);
      }
    }
  }
  else
  {
    m.set_bond_type_between_atoms(n1, c5, SINGLE_BOND);
    m.set_bond_type_between_atoms(n1, c2, SINGLE_BOND);
    m.set_bond_type_between_atoms(n3, n4, SINGLE_BOND);
    m.set_bond_type_between_atoms(c2, n3, DOUBLE_BOND);
    m.set_bond_type_between_atoms(n4, c5, DOUBLE_BOND);
  }

  return 1;
}

/*
  When standardising imidazoles, we can easily do fused rings, IF they are
  simple rings like benzene, pyridine, alternating single and double bonds

  We need to find the one ring that is fused to R
*/

static const Set_of_Atoms *
fused_to_kekule_variable_ring (Molecule & m,
                               const int ring_number,
                               IWStandard_Current_Molecule & current_molecule_data)
{
  const Set_of_Atoms & r1 = *(current_molecule_data.ringi(ring_number));

  const int nr = current_molecule_data.nrings();

// First make sure there is just a fused system size of 2

  const int fsid1 = current_molecule_data.fused_system_identifier(r1[0]);

  int fused_ring = -1;

  for (int i = 0; i < nr; ++i)
  {
    if (i == ring_number)
      continue;

    const Set_of_Atoms * ri = current_molecule_data.ringi(i);

    if (current_molecule_data.fused_system_identifier(ri->item(0)) != fsid1)
      continue;

    if (fused_ring >= 0)    // we have more than one extra ring in the fused system
      return NULL;

    fused_ring = i;
  }

  if (fused_ring < 0)   // should not happen
    return NULL;

  const Set_of_Atoms * r2 = current_molecule_data.ringi(fused_ring);

  if (6 != r2->number_elements())
    return NULL;

  if (! current_molecule_data.ring_is_aromatic()[fused_ring])
    return NULL;

// Six membered aromatic ring R2 is fused to R1. Is it strictly alternating single and double bonds

  atom_number_t p = r2->last_item();

  int previous_was = 0;

  for (int i = 0; i < 6; ++i)
  {
    const atom_number_t n = r2->item(i);

    const Bond * b = m.bond_between_atoms(p, n);

    const int sd = b->is_single_bond() ? 1 : 2;

    if (0 == previous_was)
      previous_was = sd;
    else if (sd == previous_was)
      return 0;
    else
      previous_was = sd;

    p = n;
  }

  return r2;
}

int
Chemical_Standardisation::_do_triazole (Molecule & m,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_size = current_molecule_data.ring_size();
  const int * ring_is_fused = current_molecule_data.ring_is_fused();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

  const int nr = m.nrings();

  int rc = 0;

#ifdef DEBUG_DO_TRIAZOLE
  cerr << "_do_triazole, nr " << nr << endl;
  write_isotopically_labelled_smiles(m, false, cerr);
  cerr << " startring molecule\n";
#endif

  const Set_of_Atoms * fused_triazole = NULL;

  for (int i = 0; i < nr; i++)
  {
    if (3 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    if (5 != ring_size[i])
      continue;

    const Set_of_Atoms * r = current_molecule_data.ringi(i);

    if (! ring_is_fused[i])
      ;
    else if (NULL != (fused_triazole = fused_to_kekule_variable_ring(m, i, current_molecule_data)))
      ;
    else
      continue;

    rc += _do_triazole(m, *r, fused_triazole, current_molecule_data);
  }

  if (rc)
  {
    _transform_triazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " triazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:triazole";
  }

  return rc;
}

#ifdef COUNT_ATOMIC_NUMBER_IN_SET
static int
count_atomic_number_in_set (const Molecule & m,
                            const Set_of_Atoms & e,
                            const atomic_number_t z)
{
  int rc = 0;

  const int s = e.number_elements();

  for (int i = 0; i < s; ++i)
  {
    if (m.atomic_number(e[i]) == z)
      rc++;
  }

  return rc;
}
#endif

int
Chemical_Standardisation::_do_isoxazole (Molecule & m,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_size = current_molecule_data.ring_size();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

  const int nr = m.nrings();

  int rc = 0;

//#define DEBUG_DO_ISOXAZOLE
#ifdef DEBUG_DO_ISOXAZOLE
  cerr << "_do_isoxazole, nr " << nr << endl;
  write_numbered_smiles(m, cerr);
  cerr << ' ' << m.name() << " starting molecule\n";
#endif

  for (int i = 0; i < nr; i++)
  {
//  cerr << "Ring " << i << " rnc " << ring_nitrogen_count[i] << " aromatic " << ring_is_aromatic[i] << " size " << ring_size[i] <<endl;
    if (1 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    if (5 != ring_size[i])
      continue;

    const Set_of_Atoms * r = current_molecule_data.ringi(i);

    rc += _do_isoxazole(m, *r, current_molecule_data);
  }

  if (rc)
  {
    _transform_isoxazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " isoxazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:isoxazole";
  }

  return rc;
}

int
Chemical_Standardisation::_do_isoxazole (Molecule & m,
                                         const Set_of_Atoms & r,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  int n = -1;
  int o = -1; 

  assert (5 == r.number_elements());

  for (int i = 0; i < 5; ++i)
  {
    const atom_number_t j = r[i];

    const Atom * aj = m.atomi(j);

    const atomic_number_t zj = aj->atomic_number();

    if (7 == zj)
    {
      if (n >= 0)
        return 0;
      n = i;
    }
    else if (8 == zj)
    {
      if (o >= 0)
        return 0;

      o = i;
    }
    else if (6 != zj)
      return 0;
  }

  if (n < 0 || o < 0)
    return 0;

  if (1 == m.implicit_hydrogens(r[n]))    // already in the right form
    return 0;

  if (! m.are_bonded(r[o], r[n]))
    return 0;

  const Atom * an = m.atomi(r[n]);

  if (2 != an->ncon())
    return 0;

// Need to find a 3 connected atom adjacent to the Nitrogen

  atom_number_t anchor = INVALID_ATOM_NUMBER;
  atom_number_t exocyclic = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 2; ++i)
  {
    const Bond * b = an->item(i);

    if (b->is_single_bond())   // to the oxygen
      continue;

    const atom_number_t x = b->other(r[n]);

    const Atom * ax = m.atomi(x);

    if (6 != ax->atomic_number() || 3 != ax->ncon())
      continue;

    for (int j = 0; j < 3; ++j)
    {
      const Bond * b = ax->item(j);    // i know we are shadowing something in the outer loop

      if (! b->is_single_bond())
        continue;
      
      const atom_number_t y = b->other(x);

      const Atom * ay = m.atomi(y);

      if (1 != ay->ncon())
        continue;

      if (8 == ay->atomic_number() || 16 == ay->atomic_number() || 7 == ay->atomic_number())
      {
        anchor = x;
        exocyclic = y;
        break;
      }
    }
  }

  if (INVALID_ATOM_NUMBER == anchor)
    return 0;

  m.set_bond_type_between_atoms(exocyclic, anchor, DOUBLE_BOND);
  m.set_bond_type_between_atoms(anchor, r[n], SINGLE_BOND);

  return 1;
}

int
Chemical_Standardisation::_do_imidazole (Molecule & m,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_size = current_molecule_data.ring_size();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

  int nr = m.nrings();

  int rc = 0;

  for (int i = 0; i < nr; i++)
  {
    if (2 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    if (5 != ring_size[i])
      continue;

    rc += _do_imidazole(m, i, current_molecule_data);
  }

  if (rc)
  {
    _transform_imidazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " imidazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:imidazole";
  }

  return rc;
}

static int
expand_shell (const Molecule & m,
              int * visited)
{
#define VISITED_HERE 9

  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (1 != visited[i])
      continue;

    const Atom * a = m.atomi(i);

    int acon = a->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = a->other(i, j);

      if (0 == visited[k])
      {
        visited[k] = VISITED_HERE;
        rc += 5 * m.atomic_number(k) + m.ncon(k);     // just some arbitrary thing
      }
    }

    visited[i] = -1;    // do not expand from here again
  }

  for (auto i = 0; i < matoms; ++i)
  {
    if (VISITED_HERE == visited[i])
      visited[i] = 1;
  }
  
  return rc;
}

int 
Chemical_Standardisation::_swap_imidazole (Molecule & m,
                                           atom_number_t n1,
                                           atom_number_t c,
                                           atom_number_t n2) const
{
  m.set_bond_type_between_atoms(n1, c, DOUBLE_BOND);
  m.set_bond_type_between_atoms(c, n2, SINGLE_BOND);
  m.set_implicit_hydrogens_known(n1, 0);
  m.set_implicit_hydrogens_known(n2, 0);

  return 1;
}

int
Chemical_Standardisation::_do_imidazole (Molecule & m,
                                         const int ring_number,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  const Set_of_Atoms & r = *(current_molecule_data.ringi(ring_number));

  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const Atom * const * atoms = current_molecule_data.atoms();

  int ndxh0 = -1;         // index of nitrogen atom with zero hydrogens
  int ndxh1 = -1;         // index of nitrogen atom with one  hydrogens
  int c31 = -1;         // indices of carbon atoms
  int c32 = -1;
  int c33 = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t a = r[i];

    atomic_number_t za = z[a];

    if (7 == za)
    {
      if (0 != m.formal_charge(a))
        continue;

      if (m.hcount(a))
      {
        if (ndxh1 >= 0)
          return 0;

        ndxh1 = i;
      }
      else
      {
        if (ndxh0 >= 0)
          return 0;

        ndxh0 = i;
      }
    }
    else if (6 == za)
    {
      if (c31 < 0)
        c31 = i;
      else if (c32 < 0)
        c32 = i;
      else if (c33 < 0)
        c33 = i;
      else              // Not an imidazole
        return 0;
    }
    else
      return 0;
  }

  if (ndxh0 < 0 || ndxh1 < 0)
    return 0;

  if (c31 < 0 || c32 < 0 || c33 < 0)
    return 0;

// Make sure this is an imidazole - nitrogens separated by one atom

  atom_number_t nh1, c1, nh0, c2, c3;
  if ((ndxh0 + 2) % 5 == ndxh1)
  {
    nh0 = r[ndxh0];
    c1  = r[(ndxh0 + 1) % 5];
    nh1 = r[(ndxh0 + 2) % 5];
    c3  = r[(ndxh0 + 3) % 5];
    c2  = r[(ndxh0 + 4) % 5];
  }
  else if ((ndxh1 + 2) % 5 == ndxh0)
  {
    nh1 = r[ndxh1];
    c1  = r[(ndxh1 + 1) % 5];
    nh0 = r[(ndxh1 + 2) % 5];
    c2  = r[(ndxh1 + 3) % 5];
    c3  = r[(ndxh1 + 4) % 5];
  }
  else         // nitrogens not separated by two bonds
    return 0;

//#define DEBUG_DO_IMIDAZOLE
#ifdef DEBUG_DO_IMIDAZOLE
  cerr << "Atoms " << nh1 << " " << c1 << " " << nh0 << " " << c2 << " " << c3 << endl;
#endif

  const Bond * b = atoms[c1]->bond_to_atom(nh0);
  if (! b->is_double_bond())
    return 0;

// At this stage we have an imidazole. See if we can resolve it by
// number of connections at c2 and c3

/*             c1
             /    \
            /      \
         nh1        nh0
           |        |
          c3 ------ c2
*/

  assert (m.are_bonded(c2, nh0));
  assert (m.are_bonded(c3, nh1));
  assert (m.are_bonded(c2, c3));
  assert (m.are_bonded(c1, nh0));
  assert (m.are_bonded(c1, nh1));

#ifdef DEBUG_DO_IMIDAZOLE
  cerr << "Ncon " << atoms[c2]->ncon() << " and " << atoms[c3]->ncon() << endl;
#endif

  if (atoms[c2]->ncon() < atoms[c3]->ncon())   // already as we want it
    return 0;

  if (atoms[c2]->ncon() > atoms[c3]->ncon())
    return _swap_imidazole(m, nh1, c1, nh0);

// Not differentiated by the connectivity, start expanding

  int matoms = m.natoms();

  int * tmp = new_int(matoms + matoms); std::unique_ptr<int[]> free_tmp(tmp);

  int * tmp1 = tmp;
  int * tmp2 = tmp + matoms;

// Mark all atoms in the ring as visited

  r.set_vector(tmp1, -1);
  r.set_vector(tmp2, -1);

// The carbon atoms are the starting points for expansion

  tmp1[c2] = 1;
  tmp2[c3] = 1;

#ifdef DEBUG_DO_IMIDAZOLE
  cerr << "Beginning out of ring imidazole detection\n";
#endif

  while (1)
  {
    int e1 = expand_shell(m, tmp1);
    int e2 = expand_shell(m, tmp2);

#ifdef DEBUG_DO_IMIDAZOLE
    cerr << " e1 " << e1 << " and " << e2 << endl;
#endif

    if (e1 > e2)    // correct as is
      return 0;
    else if (e1 < e2)
      return _swap_imidazole(m, nh1, c1, nh0);

    if (0 == e1)    // no more expansion, cannot be resolved
      return 0;
  }

  return 1;
}

/*
  Make sure that the Nitrogen with the Hydrogen is adjacent to the carbon
*/

/*int
Chemical_Standardisation::_do_imidazole (Molecule & m,
                                         const atomic_number_t * z,
                                         const int * ncon,
                                         Atom ** atoms,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  int nr = m.nrings();

  int rc = 0;

//cerr << "_possible_imidazole " << _possible_imidazole << ", nr = " << nr << endl;

  for (int i = current_molecule_data.possible_imidazole(); i < nr; i++)
  {
    const Ring * r = m.ringi(i);
    if (5 != r->number_elements())
      continue;

    if (r->is_fused())
      continue;

    rc += _do_imidazole(m, *r, z, ncon, atoms);
  }

  if (rc)
  {
    _transform_imidazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " imidazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:imidazole";
  }

  return rc;
}*/

int
Chemical_Standardisation::_do_pyrazole (Molecule & m,
                                        int * atom_already_done,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();
  const int * ring_size = current_molecule_data.ring_size();

  int nr = m.nrings();

  int rc = 0;

//#define DEBUG_DO_PYRAZOLE
#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Processing pyrazoles, nrings " << nr << endl;
#endif

  for (int i = 0; i < nr; i++)
  {
    if (2 != ring_nitrogen_count[i])
      continue;

    if (! ring_is_aromatic[i])
      continue;

    if (5 != ring_size[i])
      continue;

    rc += _do_pyrazole(m, atom_already_done, i, current_molecule_data);
  }

#ifdef DEBUG_DO_PYRAZOLE
  if (rc)
  {
    write_numbered_smiles(m, cerr);
    cerr << " pyrazole\n";
  }
  else
    cerr << "No Pyrazole change\n";
#endif

  if (rc)
  {
    _transform_pyrazole.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " pyrazole groups\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:pyrazole";
  }

  return rc;
}

int
Chemical_Standardisation::_do_tetrazole (Molecule & m,
                                         const Set_of_Atoms & r,
                                         IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int c = -1;     // index of the Carbon atom
  
  atom_number_t nitrogen_with_one_hydrogen = INVALID_ATOM_NUMBER;
  int nitrogens_with_no_implicit_hydrogens = 0;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    Atom * aj = const_cast<Atom *>(atoms[j]);    // the implicit_hydrogens method is possibly non-const

    if (aj->formal_charge() > 0)   // we may have a negatively charged N in the ring
      continue;

    if (7 == z[j] && 2 == ncon[j])
    {
      if (-1 == aj->formal_charge())
      {
        if (INVALID_ATOM_NUMBER != nitrogen_with_one_hydrogen)
          return 0;

        nitrogen_with_one_hydrogen = j;
        continue;
      }

      int hcount = aj->implicit_hydrogens();
      if (0 == hcount)
      {
        nitrogens_with_no_implicit_hydrogens++;
        continue;
      }

      if (INVALID_ATOM_NUMBER != nitrogen_with_one_hydrogen)
        return 0;

      nitrogen_with_one_hydrogen = j;
    }
    else if (6 == z[j] && 3 == ncon[j] && 4 == aj->nbonds())
    {
      if (c >= 0)     // can be only one Carbon atom in a tetrazole
        return 0;

      c = i;
    }
    else
      return 0;
  }

  if (c < 0 || nitrogen_with_one_hydrogen < 0 || 3 != nitrogens_with_no_implicit_hydrogens)
    return 0;

// Identify the four nitrogens

  atom_number_t n1 = r.next_after_wrap(c);

  if (nitrogen_with_one_hydrogen == n1)
    return 0;

  atom_number_t n2 = r.next_after_wrap(c);
  atom_number_t n3 = r.next_after_wrap(c);
  atom_number_t n4 = r.next_after_wrap(c);

  if (nitrogen_with_one_hydrogen == n4)
    return 0;

// Not in the right form. 

  if (2 == atoms[n3]->nbonds())
    m.set_bond_type_between_atoms(n1, n2, SINGLE_BOND);
  else
    m.set_bond_type_between_atoms(n3, n4, SINGLE_BOND);

  m.set_bond_type_between_atoms(n2, n3, DOUBLE_BOND);

  if (-1 == m.formal_charge(nitrogen_with_one_hydrogen))
    m.set_formal_charge(nitrogen_with_one_hydrogen, 0);

  m.set_implicit_hydrogens_known(nitrogen_with_one_hydrogen, 0);

  return 1;
}

/*
  Make sure the NH is adjacent to the carbon
*/

static int
switch_pyrazole (Molecule & m,
                  atom_number_t c1,
                  atom_number_t n1,
                  atom_number_t n2,
                  atom_number_t c2,
                  atom_number_t c_opposite)
{
  assert (0 == m.hcount(n1));

  m.set_bond_type_between_atoms(n1, c1, SINGLE_BOND);
  m.set_bond_type_between_atoms(c2, c_opposite, SINGLE_BOND);

  m.set_implicit_hydrogens_known(n2, 0);

  m.set_bond_type_between_atoms(c1, c_opposite, DOUBLE_BOND);
  m.set_bond_type_between_atoms(n2, c2, DOUBLE_BOND);

  return 1;
}

/*
  We have a pyrazole and need to know if is part of a fused
  aromatic system
*/

static int
determine_fused_aromatic_pyrazole (const Molecule & m,
                                   const int ring_number,
                                   const atom_number_t c1,
                                   const atom_number_t c_opposite,
                                   const atom_number_t c2,
                                   IWStandard_Current_Molecule & current_molecule_data)
{
  if (! current_molecule_data.ring_is_fused()[ring_number])
    return 0;

  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  if (3 != ncon[c_opposite])   // must be fused in some strange way
    return 0;

  for (int i = 0; i < 3; i++)
  {
    atom_number_t j = atoms[c_opposite]->other(c_opposite, i);

    if (c1 == j || c2 == j)
      continue;

    return atom_is_aromatic[j];
  }

  return 0;   // should never come here
}

/*
  Basically make sure there is a double bond between atoms A1 and A2
*/

static int
switch_fused_pyrazole (Molecule & m,
                       atom_number_t a1,
                       atom_number_t a2)
{
//#define DEBUG_SWITCH_FUSED_PYRAZOLE
#ifdef DEBUG_SWITCH_FUSED_PYRAZOLE
  cerr << "switch_fused_pyrazole setting double bond between " << a1 << " and " << a2 << " begin " << m.smiles() << endl;
#endif

  Toggle_Kekule_Form tkf;

  tkf.set_allow_pyrrole_to_change(1);

  tkf.set_display_error_messages(0);

  int changed;

  if (! tkf.process(m, a1, a2, DOUBLE_BOND, changed))
  {
#ifdef DEBUG_SWITCH_FUSED_PYRAZOLE
    cerr << "Could not be changed\n";
#endif
    return 0;
  }

#ifdef DEBUG_SWITCH_FUSED_PYRAZOLE
  cerr << " switch_fused_pyrazole " << m.smiles() << " changed " << changed << endl;
#endif

  return changed;
}

int
Chemical_Standardisation::_do_pyrazole (Molecule & m,
                                        int * atom_already_done,
                                        const int ring_number,
                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const Set_of_Atoms & r = *(current_molecule_data.ringi(ring_number));

  assert (5 == r.number_elements());

#ifdef DEBUG_DO_PYRAZOLE
  write_numbered_smiles(m, cerr);
  cerr << " pyrazole atoms " << r << endl;
#endif

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int n1_index_in_ring = -1;
  int n2_index_in_ring = -1;

  for (int i = 0; i < 5; i++)
  {
    atom_number_t j = r[i];

    if (7 != z[j])
      continue;

    if (2 != ncon[j])
      continue;

    if (n1_index_in_ring < 0)   // there must be exactly two nitrogens in the ring
      n1_index_in_ring = i;
    else if (n2_index_in_ring < 0)
      n2_index_in_ring = i;
    else
      return 0;
  }

//cerr << "_do_pyrazole: indices " << n1_index_in_ring << " and " << n2_index_in_ring << endl;

  if (n2_index_in_ring < 0)
    return 0;

  if (n1_index_in_ring > n2_index_in_ring)
  {
    int tmp = n1_index_in_ring;
    n1_index_in_ring = n2_index_in_ring;
    n2_index_in_ring = tmp;
  }
//  std::swap(n1_index_in_ring, n2_index_in_ring);

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Got two nitrogens, n1 " << n1_index_in_ring << " and n2 " << n2_index_in_ring << endl;
#endif

// Thet two nitrogens must be adjacent. Identify the adjoining carbons

  int c1_index = -1;;
  int c2_index = -1;;
  int c_opposite_index = -1;

  if (0 == n1_index_in_ring && 1 == n2_index_in_ring)
  {
    c1_index = 4;
    c2_index = 2;
    c_opposite_index = 3;
  }
  else if (1 == n1_index_in_ring && 2 == n2_index_in_ring)
  {
    c1_index = 0;
    c2_index = 3;
    c_opposite_index = 4;
  }
  else if (2 == n1_index_in_ring && 3 == n2_index_in_ring)
  {
    c1_index = 1;
    c2_index = 4;
    c_opposite_index = 0;
  }
  else if (3 == n1_index_in_ring && 4 == n2_index_in_ring)
  {
    c1_index = 2;
    c2_index = 0;
    c_opposite_index = 1;
  }
  else if (0 == n1_index_in_ring && 4 == n2_index_in_ring)
  {
    c1_index = 1;
    c2_index = 3;
    c_opposite_index = 2;
  }
  else
    return 0;

// Now we need to differentiate the two nitrogens...

  atom_number_t c1 = r[c1_index];
  atom_number_t n1 = r[n1_index_in_ring];
  atom_number_t n2 = r[n2_index_in_ring];
  atom_number_t c2 = r[c2_index];

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Pyrazole atoms c1 " << c1 << " n1 " << n1 << " n2 " << n2 << " c2 " << c2 << endl;
  cerr << m.smiles() << endl;
#endif

  if (6 != atoms[c1]->atomic_number())
    return 0;

  if (6 != atoms[c2]->atomic_number())
    return 0;

  if (2 == ncon[c1] && 2 == ncon[c2])    // cannot do anything, symmetric and isolated
    return 0;

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "First two carbons OK\n";
#endif

  atom_number_t c_opposite = r[c_opposite_index];

  if (6 != atoms[c_opposite]->atomic_number())
    return 0;

  int h1 = m.hcount(n1);
  int h2 = m.hcount(n2);

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Pyrazole atoms c1 " << c1 << " n1 " << n1 << " (" << h1 << " H) n2 " << n2 << " (" << h2 << " H) c2 " << c2 << " opposite " << c_opposite << endl;
#endif

  if (0 == h1 && 0 == h2)
    return 0;

  if (h1 && h2)
    return 0;

// Ensure that n1 is the atom with the hydrogen

  if (h1)
  {
    if (3 != atoms[n2]->nbonds())
      return 0;
  }
  else
  {
    if (3 != atoms[n1]->nbonds())
      return 0;

    std::swap(n1, n2);
    std::swap(c1, c2);
  }

  assert (m.hcount(n1) > 0);
  assert (m.bond_between_atoms(n1, c1)->is_single_bond());
  assert (m.bond_between_atoms(n2, c2)->is_double_bond());

// Many pyrazoles are fused, but to an aliphatic ring, so 
// they do not need to be treated by Toggle_Kekule_Form

  int fused_aromatic = determine_fused_aromatic_pyrazole(m, ring_number, c1, c_opposite, c2, current_molecule_data);

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "After possible swap, pyrazole atoms c1 " << c1 << " n1 " << n1 << " n2 " << n2 << " c2 " << c2 << " opposite " << c_opposite << " fused? " << fused_aromatic << " connections " << ncon[c1] << ' ' << ncon[c2] << " fused? " << fused_aromatic << endl;
  cerr << "Carbons " << c1 << " ncon " << ncon[c1] << ' ' << m.smarts_equivalent_for_atom(c1) << " and " << c2 << " ncon " << ncon[c2] << " " << m.smarts_equivalent_for_atom(c2) << endl;
#endif

// See if we can resolve things by connectivity

  if (ncon[c1] < ncon[c2])
  {
    if (fused_aromatic)
      return switch_fused_pyrazole(m, c1, n1);
    else
      return switch_pyrazole(m, c2, n2, n1, c1, c_opposite);
  }
  else if (ncon[c1] > ncon[c2])
    return 0;

// Both adjacent carbons have 3 connections. Resolve by shells

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Pyrazole being resolved by shell expansion\n";
  cerr << m.smiles() << endl;
#endif

  const int matoms = m.natoms();

  int * tmp = new_int(matoms + matoms); std::unique_ptr<int[]> free_tmp(tmp);
  int * tmp1 = tmp;
  int * tmp2 = tmp1 + matoms;

  r.set_vector(tmp1, -1);
  r.set_vector(tmp2, -1);
  tmp1[c1] = 1;
  tmp2[c2] = 1;

  while (1)
  {
    int e1 = expand_shell(m, tmp1);
    int e2 = expand_shell(m, tmp2);

//  cerr << "Expanded to " << e1 << " and " << e2 << endl;

    if (e1 == e2)
    {
      if (0 == e1)   // done, cannot be resolved
        return 0;

      continue;
    }

    if (e1 < e2)
      return 0;

    if (fused_aromatic)
      return switch_fused_pyrazole(m, c2, n2);
    else
      return switch_pyrazole(m, c2, n2, n1, c1, c_opposite);
  }

#ifdef DEBUG_DO_PYRAZOLE
  cerr << "Pyrazole not resolved by shell expansion\n";
#endif
  return 0;   // cannot figure out what to do
}


/*
  Take N=N to [N+]-[N-]

  We also do CN(C)(C)=N
*/

int
Chemical_Standardisation::_do_transform_back_to_nplus_nminus (Molecule & m,
                                          IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (ncon[i] < 3)
      continue;

    const Atom * ai = atoms[i];

    if (0 != ai->formal_charge())
      continue;

    if (5 != ai->nbonds())
      continue;

    atom_number_t doubly_bonded_nitrogen = INVALID_ATOM_NUMBER;
    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item(j);
      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other(i);
      if (7 != z[k])
        continue;

      if (0 != atoms[k]->formal_charge())
        continue;

      doubly_bonded_nitrogen = k;
      break;
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_nitrogen)
      continue;

    m.set_formal_charge(i, 1);
    m.set_formal_charge(doubly_bonded_nitrogen, -1);
    m.set_bond_type_between_atoms(i, doubly_bonded_nitrogen, SINGLE_BOND);
    rc++;
  }

  if (rc)
  {
    _transform_back_to_nplus_nminus.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " =N=N- to =[N+]-[N-]-\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:RN+N-";
  }

  return rc;
}

/*
  Convert N#N=N- to [N-]=[N+]=N-
*/

int 
Chemical_Standardisation::_do_transform_azid_to_charge_separated (Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
//const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  int matoms = m.natoms();

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * a = atoms[i];

    if (0 != a->formal_charge())
      continue;

    const Bond * b = a->item(0);

    if (! b->is_triple_bond())
      continue;

    atom_number_t n2 = b->other(i);

    if (7 != z[n2])
      continue;

    if (2 != ncon[n2])
      continue;

    const Atom * a2 = atoms[n2];

    if (5 != a2->nbonds())
      continue;

    atom_number_t n3;

    if (i == a2->other(n2, 0))
      n3 = a2->other(n2, 1);
    else
      n3 = a2->other(n2, 0);

    if (7 != z[n3])
      continue;

    if (0 != atoms[n3]->formal_charge())
      continue;

    m.set_formal_charge(i, -1);
    m.set_formal_charge(n2, 1);
    m.set_bond_type_between_atoms(i, n2, DOUBLE_BOND);
    rc++;
  }

  if (rc)
  {
    _transform_to_charge_separated_azid.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " N#N=N- to N=[N+]=[N-]-\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:RAZID";
  }

  return rc;
}

/*
  Change -N=C(-[OH])-N to -N-C(=O)-N
*/

int
Chemical_Standardisation::_do_transform_misdrawn_urea (Molecule & m,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();
  const int * atom_is_aromatic = current_molecule_data.atom_is_aromatic();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (1 != ncon[i])
      continue;

    const Atom * ai = atoms[i];

    if (1 != ai->nbonds())
      continue;

    atom_number_t c = ai->other(i, 0);

    if (3 != ncon[c])
      continue;

    if (atom_is_aromatic[c])
      continue;

    const Atom * ac = atoms[c];

    if (4 != ac->nbonds())
      continue;

    atom_number_t singly_bonded_nitrogen = INVALID_ATOM_NUMBER;
    atom_number_t doubly_bonded_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon[c]; j++)
    {
      const Bond * b = ac->item(j);

      atom_number_t k = b->other(c);

      if (7 != z[k])
        continue;

// don't change N1C(=O)C=CC2=C1C=CC=C2.C(=O)(O)/C=C/C(=O)O PBCHM6038

      if (b->is_double_bond())
      {
        if (! b->part_of_cis_trans_grouping())
          doubly_bonded_nitrogen = k;
      }
      else
        singly_bonded_nitrogen = k;
    }

    if (INVALID_ATOM_NUMBER == doubly_bonded_nitrogen || INVALID_ATOM_NUMBER == singly_bonded_nitrogen)
      continue;

    m.set_bond_type_between_atoms(i, c, DOUBLE_BOND);
    m.set_bond_type_between_atoms(c, doubly_bonded_nitrogen, SINGLE_BOND);
    rc++;
  }

  if (rc)
  {
    _transform_misdrawn_urea.extra(rc);

    if (_verbose)
      cerr << "Transformed " << rc << " misdrawn ureas\n";

    if (_append_string_depending_on_what_changed)
      _append_to_changed_molecules << " STD:MSDUR";
  }

  return rc;
}

int
Chemical_Standardisation::_do_transform_implicit_hydrogen_known_errors (Molecule & m,
                                                        IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ncon = current_molecule_data.ncon();
  const Atom * const * atoms = current_molecule_data.atoms();

  int rc = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    Atom * ai = const_cast<Atom *>(atoms[i]);

    if (! ai->implicit_hydrogens_known())
      continue;

    if (ai->valence_ok())
      continue;

//  we have an atom with implicit hydrogens known, and an invalid valence.

    int ih = ai->implicit_hydrogens();

    if (ncon[i] >= ai->element()->normal_valence() + ai->formal_charge())
      continue;

    int hdiff = ncon[i] - ai->element()->normal_valence() - ai->formal_charge();

//  cerr << " ih = " << ih << " hdiff " << hdiff << endl;

    ai->set_implicit_hydrogens(ih - hdiff, 1);
    rc++;
  }

  return rc;
}

static int
single_bond_to_oxygen (const Molecule & m,
                       atom_number_t n)
{
  const Atom * a = m.atomi(n);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    int o = b->other(n);

    if (8 == m.atomic_number(o))
      return 1;
  }

  return 0;
}

/*
  We have identified a possible lactim, but want to make sure the nitrogen is not shared
  between two lactim forms
*/

/*static int
single_bond_to_carbon_oxygen (const Molecule & m,
                              const atom_number_t nitrogen,
                              const atom_number_t carbon)
{
  const Atom * an = m.atomi(nitrogen);

  const int ncon = an->ncon();

  for (int i = 0; i < ncon; ++i)
  {
    atom_number_t c = an->other(nitrogen, i);

    if (c == carbon)
      continue;

    if (6 != m.atomic_number(c))
      continue;

    const Atom * ac = m.atomi(c);

    const int ccon = ac->ncon();

    for (auto j = 0; j < ccon; ++j)
    {
      const Bond * b = ac->item(j);

      if (! b->is_single_bond())
        continue;

      atom_number_t o = b->other(c);

      if (o == nitrogen)
        continue;

      if (8 != m.atomic_number(o))
        continue;

      if (1 == m.ncon(o))
        return 1;
    }
  }

  return 0;
}*/

//#define DEBUG_DO_TRANSFORM_LACTIM


#ifdef LATER_WORK_ENABLING_BACKOUT
static int
gather_existing_bond_types (Molecule & m,
                            resizable_array<bond_type_t> & existing_bond_order)
{
  const int nb = m.nedges();

  existing_bond_order.resize(nb);

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = m.bondi(i);

    if (b->is_single_bond())
      existing_bond_order.add(SINGLE_BOND);
    else if (b->is_double_bond())
      existing_bond_order.add(DOUBLE_BOND);
    else if (b->is_triple_bond())
      existing_bond_order.add(TRIPLE_BOND);
  }

  return 1;
}
#endif

/*static void
set_all_bonds_between_aromatic_atoms_to_single_bonds (Molecule & m,
                                                      const int * aromatic)
{

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == aromatic[i])
     continue;

    const Atom * a = m.atomi(i);

    if (7 == a->atomic_number())
      m.set_implicit_hydrogens_known(i, 0);

    int acon = a->ncon();

    for (int j = 0; j < acon; ++j)
    {
      const Bond * b = a->item(j);

      if (b->is_single_bond())
        continue;

      const atom_number_t k = b->other(i);

      if (aromatic[k])
        m.set_bond_type_between_atoms(i, k, SINGLE_BOND);
    }
  }

  return;
}*/

static int
is_isolated_lactim_lactam (const Molecule & m,
                           int ndx,
                           IWStandard_Current_Molecule & current_molecule_data)
{
  const resizable_array_p<Possible_Lactim_Lactam> & possible_lactam = current_molecule_data.possible_lactam();

  const int n = possible_lactam.number_elements();

  const Possible_Lactim_Lactam * p = possible_lactam[ndx];

  for (int i = 0; i < n; ++i)
  {
    if (i == ndx)
      continue;

    if (possible_lactam[i]->shares_nitrogen_with(*p))     // p not isolated
      return 0;
  }

  return 1;     // not shared with any of the others, isolated
}

/*
  These are horrendously complicated
*/

int
Chemical_Standardisation::_do_transform_ring_lactim (Molecule & m,
                                           int * atom_already_changed,
                                           IWStandard_Current_Molecule & current_molecule_data)
{
#ifdef DEBUG_DO_TRANSFORM_LACTIM
  cerr << "Into _do_transform_ring_lactim\n";
#endif

  int rc = __do_transform_ring_lactim(m, atom_already_changed, current_molecule_data);

  if (0 == rc)
    return 0;

  _transform_lactim_lactam_ring.extra(rc);

  if (_verbose)
    cerr << "Transformed " << rc << " lactim->lactam ring\n";

  if (_append_string_depending_on_what_changed)
    _append_to_changed_molecules << " STD:LTLTR";

  return rc;
}

int
Chemical_Standardisation::__do_transform_ring_lactim (Molecule & m,
                                                int * atom_already_changed,
                                                IWStandard_Current_Molecule & current_molecule_data)
{

  const resizable_array_p<Possible_Lactim_Lactam> & possible_lactam = current_molecule_data.possible_lactam();

  const int n = possible_lactam.number_elements();

#ifdef DEBUG_DO_TRANSFORM_LACTIM
  cerr << "Into __do_transform_ring_lactim, n = " << n << "\n";
#endif

// We run into problems if we make changes, and the molecule recomputes it's rings. So, we
// save them in the order we perceived them.

  resizable_array_p<Set_of_Atoms> rings;

  int * already_done = new_int(n); std::unique_ptr<int[]> free_already_done(already_done);

// We need to keep track of how many of the N possible lactams remain
// to be processed.  Some may be tried and fail, so we decrement
// remaining_to_be_processed even if the transform does not happen

  int remaining_to_be_processed = n;

// first process any non-aromatic forms

  int rc = 0;

  for (int i = 0; i < n; ++i)
  {
    Possible_Lactim_Lactam * p = possible_lactam[i];

#ifdef DEBUG_LACTAM_LACTIM
    cerr << "Possible lactam " << p->oxygen() << ' ' << p->carbon() << ' ' << p->nitrogen() << " ring " << p->is_ring() << " arom " << p->aromatic() << " fss " << p->fused_system_size() << endl;
#endif

    if (! p->is_ring())
    {
      already_done[i] = 1;
      remaining_to_be_processed--;
      continue;
    }

    if (p->lactam_form())    // already lactim form, OK
    {
      already_done[i] = 1;
      remaining_to_be_processed--;
      continue;
    }

    if (p->aromatic())
      continue;

    rc += p->to_lactam_form(m, 1);
    already_done[i] = 1;
    remaining_to_be_processed--;
  }

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "after processing aliphatic rings " << remaining_to_be_processed << " of " << n << " remaining\n";
#endif

  if (0 == remaining_to_be_processed)
    return rc;

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "Aromatic forms present, starting smiles " << m.smiles() << endl;
#endif

  int tmp =_do_lactam_lactim_pyrazole_triazole(m, rings, atom_already_changed, already_done, current_molecule_data);

  if (tmp)
  {
    rc += tmp;
    remaining_to_be_processed -= tmp;
    if (remaining_to_be_processed <= 0)
      return rc;
  }

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "After processing pyrazole and triazole forms " << remaining_to_be_processed << " of " << n << " remaining to be processed\n";
#endif

#ifdef LATER_WORK_ENABLING_BACKOUT
  int * aromatic = new int[matoms]; std::unique_ptr<int[]> free_aromatic(aromatic);
  copy_vector(aromatic, current_molecule_data.atom_is_aromatic(), matoms);

  const int initial_aromatic_atom_count = count_non_zero_occurrences_in_array(aromatic, matoms);

  resizable_array<bond_type_t> save_bond_types;
  gather_existing_bond_types(m, save_bond_types);
#endif

// for aromatics, process any that are the only ones in a ring/ring system, and
// where there is only one nitrogen

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];

#ifdef DEBUG_LACTAM_LACTIM
    cerr << "Processing " << p->oxygen() << ' ' << p->carbon() << ' ' << p->nitrogen() << " arom? " << p->aromatic() << " nitrogen attachments " << p->total_nitrogen_attachments() << endl;
#endif

    if (p->lactims_in_fused_system() > 1)
      continue;

    if (1 != p->total_nitrogen_attachments())
      continue;

    if (_process_lactim_in_isolated_aromatic_ring(m, *p, current_molecule_data))
      rc++;

    already_done[i] = 1;
    remaining_to_be_processed--;
  }

  if (0 == remaining_to_be_processed)
    return rc;

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "After processing things in isolated rings, remaining? " << remaining_to_be_processed << endl;
#endif

// Now try things where there is just one in a fused system

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];

    if (1 != p->total_nitrogen_attachments())
      continue;

    if (1 == p->fused_system_size())   // these were tried above
      continue;

    if (1 != p->lactims_in_fused_system())
      continue;

    if (_process_lactim_in_isolated_aromatic_ring(m, *p, current_molecule_data))
      rc++;

    remaining_to_be_processed--;
  }

  if (0 == remaining_to_be_processed)
    return rc;

// Within any ring system, if there are the same number of oxygen and nitrogen
// atoms, those should be independent C12=C(C(=NC(=N1)O)O)NC=N2 PBCHM1188

  resizable_array<int> fsid_completed;

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    const Possible_Lactim_Lactam * pi = possible_lactam[i];

    const int fsid = current_molecule_data.fused_system_identifier(pi->carbon());

    if (fsid_completed.contains(fsid))
      continue;

    fsid_completed.add(fsid);

    Set_of_Atoms nitrogens;
    pi->add_unique_nitrogens(nitrogens);
    int oxygen_count = 1;

    resizable_array<int> in_system;
    in_system.add(i);

    for (int j = i + 1; j < n; ++j)
    {
      if (already_done[j])
        continue;

      const Possible_Lactim_Lactam * pj = possible_lactam[j];

      if (fsid != current_molecule_data.fused_system_identifier(pj->carbon()))
        continue;

      pj->add_unique_nitrogens(nitrogens);
      oxygen_count++;
      in_system.add(j);
    }

#ifdef DEBUG_LACTAM_LACTIM
    cerr << "In fused system " << fsid << " found " << nitrogens.size() << " nitrogens and " << oxygen_count << " oxygens " << m.smiles() << endl;
#endif

    if (nitrogens.number_elements() != oxygen_count)
    {
      int tmp = _lactim_lactam_process_if_all_groups_non_overlapping(m, in_system, already_done, current_molecule_data);

      if (tmp)
        remaining_to_be_processed -= tmp;

      continue;
    }

//  Process in two loops. First do those groups that have one Nitrogen

    for (int j = 0; j < in_system.number_elements(); ++j)
    {
      int k = in_system[j];

      Possible_Lactim_Lactam * pk = possible_lactam[k];

      if (1 != pk->total_nitrogen_attachments())
        continue;

      if (INVALID_ATOM_NUMBER != pk->alpha_nitrogen())
        continue;

      if (_process_lactim_in_isolated_aromatic_ring(m, *pk, current_molecule_data))
        rc++;
      
      already_done[k] = 1;
      remaining_to_be_processed--;
    }

//  Now do those with two nitrogens

    for (int j = 0; j < in_system.number_elements(); ++j)
    {
      int k = in_system[j];

      if (already_done[k])
        continue;

      Possible_Lactim_Lactam * pk = possible_lactam[k];

      if (2 != pk->total_nitrogen_attachments())
        continue;

      if (INVALID_ATOM_NUMBER != pk->alpha_nitrogen())
        continue;

      pk->reperceive(m);

      if (_process_lactim_in_isolated_aromatic_ring(m, *pk, current_molecule_data))
        rc++;
      
      already_done[k] = 1;
      remaining_to_be_processed--;
    }
  }

  if (0 == remaining_to_be_processed)
    return rc;

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "After dealing with balanced systems, have " << remaining_to_be_processed << " of " << n << " to process\n";
#endif

// Deal with things like C1(=NN=CN1C)S 

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];

    if (INVALID_ATOM_NUMBER != p->second_nitrogen())
      continue;

    if (1 != p->lactims_in_fused_system())
      continue;

    if (p->to_lactam_form(m, 1))
      rc++;

    already_done[i] = 1;
    remaining_to_be_processed--;
  }

  if (0 == remaining_to_be_processed)
    return rc;

// Try processing things like SC1=NC2=CC=CC=C2N1 PBCHM72611734

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "After one per system need to process " << remaining_to_be_processed << " of " << n << " possible lactams\n";
#endif

// Dead code for now...

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];

//  const int fsid = current_molecule_data.fused_system_identifier(p->carbon());

    if (2 != p->total_nitrogen_attachments())
      continue;

    const atom_number_t n1 = p->nitrogen();
    const atom_number_t n2 = p->second_nitrogen();

    if (2 != m.ncon(n1) || 2 != m.ncon(n2))
      continue;

    if (m.hcount(n1) || m.hcount(n2))   // already got a hydrogen nearby
      continue;

#ifdef DEBUG_LACTAM_LACTIM
    cerr << "Trying to switch kekule form to get a hydrogen nearby, starting with " << m.smiles() << endl;
#endif

    Toggle_Kekule_Form tkf;
    tkf.set_display_error_messages(0);
    tkf.set_allow_pyrrole_to_change(1);

    int changed;
    tkf.process(m, p->carbon(), n1, SINGLE_BOND, changed);
#ifdef DEBUG_LACTAM_LACTIM
    cerr << "Changed " << changed << endl;
#endif

    if (changed)
      p->reperceive(m);
  }

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];

    if (INVALID_ATOM_NUMBER == p->second_nitrogen())
      continue;

    if (1 != m.hcount(p->second_nitrogen()))
      continue;

    if (1 == n)    // must be isolated
      ;
    else if (! is_isolated_lactim_lactam(m, i, current_molecule_data))
      continue;

    _change_molecule_kekule_form_for_lactam_canonical(m, *p);

    if (_process_lactim_in_isolated_aromatic_ring(m, *p, current_molecule_data))
      rc++;

    already_done[i] = 1;
    remaining_to_be_processed--;
  }

  if (0 == remaining_to_be_processed)
    return rc;

//cerr << m.smiles() << ' ' << m.name() << " multiple lactims on ring/ring system\n"; //    implement something one of these days...

  for (int i = 0; i < n; ++i)     // need to implement something here....
  {
    if (already_done[i])
      continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];
  }

  return rc;
}

int
Chemical_Standardisation::_process_lactim_in_isolated_aromatic_ring (Molecule & m,
                                                Possible_Lactim_Lactam & p,
                                                IWStandard_Current_Molecule & current_molecule_data)
{
  if (! p.could_change_to_lactim_with_current_bonding(m))
  {
    Toggle_Kekule_Form tkf;

    tkf.set_display_error_messages(0);

//  cerr << "_process_lactim_in_isolated_aromatic_ring::processing " << m.smiles() << endl;

    int changed;
    tkf.process(m, p.carbon(), p.nitrogen(), DOUBLE_BOND, changed);

//  cerr << "after toggling " << m.smiles() << " changed? " << changed << endl;

    if (! changed)
      return 0;

    if (! p.could_change_to_lactim_with_current_bonding(m))    // huh????
      return 0;
  }

  p.to_lactam_form(m, 1);

  return 1;
}

int
Chemical_Standardisation::_lactim_lactam_process_if_all_groups_non_overlapping (Molecule & m,
                                const resizable_array<int> & in_system,
                                int * already_done,
                                IWStandard_Current_Molecule & current_molecule_data)
{
  const resizable_array_p<Possible_Lactim_Lactam> & possible_lactam = current_molecule_data.possible_lactam();

  const int ns = in_system.number_elements();

  Set_of_Atoms n1;

  for (int i = 0; i < ns; ++i)
  {
    const int j = in_system[i];

    assert (! already_done[j]);

    Possible_Lactim_Lactam * pj = possible_lactam[j];

    if (! n1.add_if_not_already_present(pj->nitrogen()))
      return 0;
  }

  for (int i = 0; i < ns; ++i)
  {
    const int j = in_system[i];

    Possible_Lactim_Lactam * pj = possible_lactam[j];

    pj->to_lactam_form(m, 1);

    already_done[j] = 1;
  }

  return ns;
}

/*
  We have a lactam/lactim and need to find a canonical order for the sides
*/

int
Chemical_Standardisation::_change_molecule_kekule_form_for_lactam_canonical (Molecule & m,
                                                Possible_Lactim_Lactam & p)
{
  const int matoms = m.natoms();

  int * tmp = new_int(matoms + matoms); std::unique_ptr<int[]> free_tmp(tmp);
  int * tmp1 = tmp;
  int * tmp2 = tmp1 + matoms;

  tmp1[p.oxygen()] = -1;
  tmp1[p.carbon()] = -1;
  tmp1[p.nitrogen()] = 1;

  tmp2[p.oxygen()] = -1;
  tmp2[p.carbon()] = -1;
  tmp2[p.second_nitrogen()] = 1;

  while (1)
  {
    int e1 = expand_shell(m, tmp1);
    int e2 = expand_shell(m, tmp2);

//  cerr << "Expanded to " << e1 << " and " << e2 << endl;

    if (e1 == e2)
    {
      if (0 == e1)   // done, cannot be resolved
        return 0;

      continue;
    }

    if (e1 < e2)    // we are good
      return 0;

    Toggle_Kekule_Form tkf;

    tkf.set_display_error_messages(0);

//  cerr << "_process_lactim_in_isolated_aromatic_ring::processing " << m.smiles() << endl;

    int changed;
    tkf.process(m, p.carbon(), p.second_nitrogen(), DOUBLE_BOND, changed);

    if (changed)
      p.reperceive(m);

    return 1;
  }
}

//#define DEBUG_DO_LACTAM_LACTIM_PYRAZOLE_TRIAZOLE

int
Chemical_Standardisation::_do_lactam_lactim_pyrazole_triazole (Molecule & m,
                                const resizable_array_p<Set_of_Atoms> & rings,
                                int * atom_already_changed,
                                int * already_done,
                                IWStandard_Current_Molecule & current_molecule_data)
{
  const resizable_array_p<Possible_Lactim_Lactam> & possible_lactam = current_molecule_data.possible_lactam();
  const int * ring_size = current_molecule_data.ring_size();

  const int n = possible_lactam.number_elements();

  int rc = 0;

  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

#ifdef DEBUG_DO_LACTAM_LACTIM_PYRAZOLE_TRIAZOLE
  cerr << "_do_lactam_lactim_pyrazole_triazole:examining " << n << " possible lactim forms\n";
  write_numbered_smiles(m, cerr);
  cerr << endl;
#endif

  for (int i = 0; i < n; ++i)
  {
    if (already_done[i])
      continue;

    Possible_Lactim_Lactam * p = possible_lactam[i];

#ifdef DEBUG_DO_LACTAM_LACTIM_PYRAZOLE_TRIAZOLE
    cerr << "Is lactim at " << p->carbon() << " aromatic? " << p->aromatic() << endl;
#endif

    if (! p->aromatic())
      continue;

    const int ndx = current_molecule_data.ring_number_containing_atom(p->carbon());

    if (5 != ring_size[ndx])
      continue;

    if (ring_nitrogen_count[ndx] < 2)
      continue;

    const Set_of_Atoms * r = current_molecule_data.ringi(ndx);

    if (_do_lactam_lactim_pyrazole_triazole(m, *p, *r, current_molecule_data))
    {
      already_done[i] = 1;
      rc++;
    }
  }

#ifdef DEBUG_DO_LACTAM_LACTIM_PYRAZOLE_TRIAZOLE
  cerr << "after _do_lactam_lactim_pyrazole_triazole ";
  if (rc)
  {
    write_numbered_smiles(m, cerr);
    cerr << ' ' << m.name() << " RC " << rc << endl;
  }
  else
    cerr << "_do_lactam_lactim_pyrazole_triazole, no change\n";
#endif

  return rc;
}

// handle things that look like: OC1=CC=NN1 p39

int
Chemical_Standardisation::_do_lactam_lactim_pyrazole_triazole (Molecule & m,
                                                    const Possible_Lactim_Lactam & p,
                                                    const Set_of_Atoms & r,
                                                    IWStandard_Current_Molecule & current_molecule_data)
{
  return 0;
#ifdef REALLYBADIDEATHATDOESNOTWORK
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon          = current_molecule_data.ncon();

  const atom_number_t c = p.carbon();

  int ndx = r.index(c);

  assert (ndx >= 0);

  Set_of_Atoms myring;    // a version of R that starts with the carbon

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "Pyrazole: ring " << r << " carbon is index " << ndx << " (atom " << r[ndx] << ") z " << z[r[ndx]] << endl;
#endif

  ndx--;    // safe even if 0 == ndx

  for (int i = 0; i < 5; ++i)
  {
    myring.add(r.next_after_wrap(ndx, 1));
  }

#ifdef DEBUG_LACTAM_LACTIM
  cerr << "Ring is " << myring << " atomic numbers " << z[myring[1]] << ' ' << z[myring[2]] << ' ' << z[myring[3]] << endl;
#endif

  if (7 == z[myring[1]] && 7 == z[myring[2]] && 6 == z[myring[3]])
    ;
  else
  {
    std::swap(myring[1], myring[4]);
    std::swap(myring[2], myring[3]);
//  cerr << "Reversed " << myring << " z " << z[myring[1]] << ' ' << z[myring[2]] << ' ' << z[myring[3]] << endl;

    if (7 == z[myring[1]] && 7 == z[myring[2]] && 6 == z[myring[3]])
      ;
    else
      return 0;
  }

//cerr << "Possibly after swap " << myring << endl;

  const atom_number_t n1 = myring[1];
  const atom_number_t n2 = myring[2];
  const atom_number_t x3 = myring[3];
  const atom_number_t x4 = myring[4];

  if (2 != ncon[n1])
    return 0;

  if (2 != ncon[n2])
    return 0;

//cerr << "Checking hcount " << m.hcount(n1) << " " << m.hcount(n2) << endl;

  if (0 == m.hcount(n1) && 0 == m.hcount(n2))   // must be a hydrogen somewhere
    return 0;

// ring looks like: p.oxygen - c - n1 - n2 - x3 - x4 which joins c

  m.set_formal_charge(p.oxygen(), 0);
  m.set_bond_type_between_atoms(c,  x4, SINGLE_BOND);
  m.set_bond_type_between_atoms(n2, x3, SINGLE_BOND);
  m.unset_all_implicit_hydrogen_information(n1);
  m.unset_all_implicit_hydrogen_information(n2);

  m.set_bond_type_between_atoms(p.oxygen(), p.carbon(), DOUBLE_BOND);
  m.set_bond_type_between_atoms(x3, x4, DOUBLE_BOND);

  cerr << "Pyrazole/triazole result " << m.smiles() << ' ' << m.name() << endl;

  return 1;
#endif
}

int
Chemical_Standardisation::_do_transform_non_ring_lactim (Molecule & m,
                                    int * atom_already_changed,
                                    IWStandard_Current_Molecule & current_molecule_data)
{
#ifdef DEBUG_DO_TRANSFORM_LACTIM
  cerr << "Into _do_transform_non_ring_lactim " << m.smiles() << ' ' << m.name() << endl;
#endif

  const resizable_array_p<Possible_Lactim_Lactam> & possible_lactam = current_molecule_data.possible_lactam();

  const int n = possible_lactam.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    Possible_Lactim_Lactam * p = possible_lactam[i];

    if (p->is_ring())
      continue;

    if (INVALID_ATOM_NUMBER != p->second_nitrogen())
      continue;

    if (p->to_lactam_form(m, 1))
      rc++;
  }

  for (int i = 0; i < n; i++)
  {
    Possible_Lactim_Lactam * p = possible_lactam[i];

    if (p->is_ring())
      continue;

    if (INVALID_ATOM_NUMBER == p->second_nitrogen())
      continue;

    if (p->to_lactam_form(m, 1))
      rc++;
  }

  if (0 == rc)
    return 0;

  _transform_lactim_lactam.extra(rc);

  if (_verbose)
    cerr << "Transformed " << rc << " lactim->lactam\n";

  if (_append_string_depending_on_what_changed)
    _append_to_changed_molecules << " STD:LTLT";

  return rc;
}

int
Chemical_Standardisation::_do_transform_pyrazolone(Molecule & m, int * atom_already_changed,
                                                   IWStandard_Current_Molecule & current_molecule_data)
{
  const int * ring_is_aromatic = current_molecule_data.ring_is_aromatic();
  const int * ring_is_fused = current_molecule_data.ring_is_fused();
  const int * ring_size = current_molecule_data.ring_size();
  const int * ring_nitrogen_count = current_molecule_data.ring_nitrogen_count();

  int rc = 0;

  const int nr = current_molecule_data.nrings();

  if (0 == nr)
    return 0;

  for (int i = 0; i < nr; ++i)
  {
    if (! ring_is_aromatic[i])
      continue;

    if (5 != ring_size[i])
      continue;

    if (2 != ring_nitrogen_count[i])
      continue;

    const Set_of_Atoms * ri = current_molecule_data.ringi(i);

    rc += _do_transform_pyrazolone(m, *ri, ring_is_fused[i], atom_already_changed, current_molecule_data);
  }

  if (0 == rc)
    return 0;

  _transform_pyrazolone.extra(rc);

  if (_verbose)
    cerr << "Transformed " << rc << " pyrazolone\n";

  if (_append_string_depending_on_what_changed)
    _append_to_changed_molecules << " STD:pyrazolone";

  return rc;
}

static int
identity_exocyclic_singly_bonded_oxygen(const Molecule & m,
                                        const atom_number_t c,
                                        atom_number_t & oxygen,
                                        const IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon          = current_molecule_data.ncon();

  assert (3 == ncon[c]);
  assert (6 == z[c]);

  const Atom * a = m.atomi(c);

  for (int i = 0; i < 3; ++i)
  {
    const Bond * b = a->item(i);

    if (! b->is_single_bond())
      continue;

    const atom_number_t j = b->other(c);

    if (8 == z[j] && 1 == ncon[j])
    {
      oxygen = j;
      return 1;
    }
  }

  return 0;
}

/*
  Change C1(=CC(=NN1)O)N to O=C1NNC=C1
*/

int
Chemical_Standardisation::_do_transform_pyrazolone(Molecule & m,
                                Set_of_Atoms const & r, 
                                const int ring_is_fused,
                                int * atom_already_changed,
                                IWStandard_Current_Molecule & current_molecule_data)
{
  assert (5 == r.number_elements());

  const atomic_number_t * z = current_molecule_data.atomic_number();
  const int * ncon          = current_molecule_data.ncon();

  int ndx_exocyclic_singly_bonded_oxygen = -1;   // the carbon atom with the oxygen
  int ndx_n1 = -1;
  int ndx_n2 = -1;
  atom_number_t exocyclic_singly_bonded_oxygen = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 5; ++i)
  {
    const atom_number_t j = r[i];

    if (7 == z[j])
    {
      if (ndx_n1 < 0)
        ndx_n1 = i;
      else
        ndx_n2 = i;
    }
    else if (2 == ncon[j])
      ;
    else if (identity_exocyclic_singly_bonded_oxygen(m, j, exocyclic_singly_bonded_oxygen, current_molecule_data))
      ndx_exocyclic_singly_bonded_oxygen = i;
  }
   
  if (ndx_n2 < 0 || ndx_exocyclic_singly_bonded_oxygen < 0)
    return 0;

// Ensure that the Nitrogen atoms are adjacent

  if ((ndx_n1 + 1) % 5 == ndx_n2)
    ;
  else if ((ndx_n2 + 1) % 5 == ndx_n1)
    ;
  else
    return 0;

  if (7 != z[r[ndx_n1]] || 7 != z[r[ndx_n2]])
    return 0;

/* 
  Identify the atoms involved

        O
        |
       / \
      /   \
     c1    n1
     |     |
     |     |
     c2 -- n2
*/

  atom_number_t c1, c2, n1, n2;
  if (7 == z[r[(ndx_exocyclic_singly_bonded_oxygen+1)%5]])
  {
    n1 = r[(ndx_exocyclic_singly_bonded_oxygen+1)%5];
    n2 = r[(ndx_exocyclic_singly_bonded_oxygen+2)%5];
    c2 = r[(ndx_exocyclic_singly_bonded_oxygen+3)%5];
    c1 = r[(ndx_exocyclic_singly_bonded_oxygen+4)%5];
  }
  else if (6 == z[r[(ndx_exocyclic_singly_bonded_oxygen+1)%5]])
  {
    c1 = r[(ndx_exocyclic_singly_bonded_oxygen+1)%5];
    c2 = r[(ndx_exocyclic_singly_bonded_oxygen+2)%5];
    n2 = r[(ndx_exocyclic_singly_bonded_oxygen+3)%5];
    n1 = r[(ndx_exocyclic_singly_bonded_oxygen+4)%5];
  }
  else
    return 0;

  const atom_number_t c = r[ndx_exocyclic_singly_bonded_oxygen];
  assert (6 == z[c] && 3 == ncon[c]);

  if (7 != z[n1] || 7 != z[n2])
    return 0;

  const Bond * b = m.bond_between_atoms(c, n1);

  if (b->is_double_bond())    // the easy case OC1=NNC=C1
  {
    m.set_bond_type_between_atoms(c, exocyclic_singly_bonded_oxygen, DOUBLE_BOND);
    m.set_bond_type_between_atoms(c, n1, SINGLE_BOND);
    m.set_implicit_hydrogens_known(n1, 0);
  }
  else                  // The harder case of O-C1=CC=NN1
  {
    const Bond * b = m.bond_between_atoms(c, c1);
    if (! b->is_double_bond())
      return 0;
    b = m.bond_between_atoms(c2, n2);
    if (! b->is_double_bond())
      return 0;

    m.set_bond_type_between_atoms(c,  exocyclic_singly_bonded_oxygen, DOUBLE_BOND);
    m.set_bond_type_between_atoms(c,  n1, SINGLE_BOND);
    m.set_bond_type_between_atoms(n1, n2, SINGLE_BOND);
    m.set_bond_type_between_atoms(n2, c2, SINGLE_BOND);
    m.set_bond_type_between_atoms(c2, c1, DOUBLE_BOND);
    m.set_bond_type_between_atoms(c1, c, SINGLE_BOND);
    m.set_implicit_hydrogens_known(n1, 0);
    m.set_implicit_hydrogens_known(n2, 0);
  }

  atom_already_changed[n1] = 1;
  atom_already_changed[n2] = 1;
  atom_already_changed[c] = 1;
  atom_already_changed[exocyclic_singly_bonded_oxygen] = 1;

//write_numbered_smiles(m, cerr);
//cerr << " pyrazolone changed, n1 = " << n1 << " n2 " << n2 << " c2 " << c2 << " c1 " << c1 << " c " << c << " oxygen " << exocyclic_singly_bonded_oxygen << endl;

  return 1;
}

int
Chemical_Standardisation::activate_from_corina_transformations()
{
  _transform_nitro.activate();
  _transform_azid.activate();
  _transform_nplus_ominus.activate();
  _active = 1;

  return 1;
}

/*
  Corina puts a negative charge on -P(=O)(=O)-
*/

#ifdef PMINUS

int
Chemical_Standardisation::_do_p_minus (Molecule & m,
                                       IWStandard_Current_Molecule & current_molecule_data)
{
  const atomic_number_t * z = current_molecule_data.atomic_number();
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (15 != z[i])
      continue;
  }
}
#endif

IWStandard_Current_Molecule::IWStandard_Current_Molecule()
{
  _atomic_number = NULL;
  _ncon = NULL;
  _ring_membership  = NULL;
  _ring_size        = NULL;
  _ring_is_fused    = NULL;
  _atom_is_aromatic = NULL;
  _atom = NULL;
  _ring_nitrogen_count = NULL;
  _fsid = NULL;

  _npos = 0;
  _nneg = 0;
  _ominus = 0;
  _sminus = 0;
  _splus = 0;
  _nplus = 0;
  _cminus = 0;
  _phosphorus = 0;
  _sulphur = 0;
  _isolated_metal = 0;
  _isolated_halogen = 0;
  _singly_connected_metal = 0;
  _possible_guanidine = 0;
  _phosphorus = 0;
  _explicit_hydrogen_count = 0;
  _possible_valence_errors = 0;

  _nitrogens = 0;
  _oxygens = 0;

  return;
}

IWStandard_Current_Molecule::~IWStandard_Current_Molecule ()
{
  if (NULL != _atomic_number)
    delete [] _atomic_number;

  if (NULL != _ncon)
    delete [] _ncon;

  if (NULL != _ring_membership)
    delete [] _ring_membership;

  if (NULL != _ring_size)
    delete [] _ring_size;

  if (NULL != _ring_is_fused)
    delete [] _ring_is_fused;

  if (NULL != _atom_is_aromatic)
    delete [] _atom_is_aromatic;

  if (NULL != _atom)
    delete [] _atom;

  if (NULL != _ring_nitrogen_count)
    delete [] _ring_nitrogen_count;

  if (NULL != _ring_is_aromatic)
    delete [] _ring_is_aromatic;

  if (NULL != _fsid)
    delete [] _fsid;

  return;
}

static int
gather_bonds_to_nitrogen_atoms (const Molecule & m,
                                const atom_number_t carbon,
                                const atom_number_t oxygen,
                                resizable_array<const Bond *> & bonds_to_nitrogen)
{
  const Atom * ac = m.atomi(carbon);

  assert (3 == ac->ncon());

  for (int i = 0; i < 3; ++i)
  {
    const Bond * b = ac->item(i);

    const atom_number_t n = b->other(carbon);

    if (n == oxygen)
      continue;

    if (7 != m.atomic_number(n))
      continue;

    bonds_to_nitrogen.add(b);
  }

  if (0 == bonds_to_nitrogen.number_elements())
    return 0;

  if (1 == bonds_to_nitrogen.number_elements())
    return 1;

  if (bonds_to_nitrogen[0]->is_single_bond() && bonds_to_nitrogen[1]->is_double_bond())
    bonds_to_nitrogen.swap_elements(0, 1);

  return bonds_to_nitrogen.number_elements();
}

int
Possible_Lactim_Lactam::_alpha_nitrogen_present (const Molecule & m,
                                                 const atom_number_t n) const
{
  const Atom * a = m.atomi(n);

  if (2 != a->ncon())    // should never happen
    return 0;

  for (int i = 0; i < 2; ++i)
  {
    atom_number_t j = a->other(n, i);

    if (j == _carbon)
      continue;

    return 7 == m.atomic_number(j);
  }

  return 0;         // should never come here
}

int
Possible_Lactim_Lactam::discern_alpha_nitrogen (Molecule & m)
{
  if (! _aromatic)
    return 0;

  if (_alpha_nitrogen_present(m, _nitrogen))
  {
    _alpha_nitrogen = _nitrogen;
    return 1;
  }

  if (INVALID_ATOM_NUMBER == _second_nitrogen)
    return 0;

  if (_alpha_nitrogen_present(m, _second_nitrogen))
  {
    _alpha_nitrogen = _second_nitrogen;
    return 1;
  }

  return 0;
}

#ifdef NOT_BEING_USED
static int
remove_aromatic_lactims_with_alpha_nitrogen (Molecule & m,
                                             resizable_array_p<Possible_Lactim_Lactam> & possible_lactim)
{
  return 0;
  for (int i = possible_lactim.number_elements() - 1; i >= 0; --i)
  {
    const Possible_Lactim_Lactam * pi = possible_lactim[i];

    if (! pi->aromatic())
      continue;

    const Atom * an = m.atomi(pi->nitrogen());

    assert (2 == an->ncon());

    for (int j = 0; j < 2; ++j)
    {
      atom_number_t k = an->other(pi->nitrogen(), j);

      if (7 == m.atomic_number(k) && m.is_aromatic(k) && m.hcount(k))    // second check for aromaticity probably not necessary
      {
        possible_lactim.remove_item(i);
        break;
      }
    }
  }

  return 1;
}
#endif

static int
remove_two_lactims_sharing_a_single_nitrogen (resizable_array_p<Possible_Lactim_Lactam> & possible_lactim)
{
  const int n = possible_lactim.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const Possible_Lactim_Lactam * pi = possible_lactim[i];

    if (1 != pi->total_nitrogen_attachments())
      continue;

    for (int j = i + 1; j < n; ++j)
    {
      const Possible_Lactim_Lactam * pj = possible_lactim[j];

      if (1 != pj->total_nitrogen_attachments())
        continue;

      if (pi->nitrogen() != pj->nitrogen())
        continue;

      possible_lactim.remove_item(j);
      possible_lactim.remove_item(i);
      return remove_two_lactims_sharing_a_single_nitrogen(possible_lactim);    // just lazy, this will hardly evern happen
    }
  }

  return 1;
}

static int
carbonyl_present (const Molecule & m,
                  const atom_number_t carbon,
                  const atom_number_t nitrogen)
{
  const Atom * n = m.atomi(nitrogen);

  if (2 != n->ncon())
    return 0;

  for (int i = 0; i < 2; ++i)
  {
    const atom_number_t c2 = n->other(nitrogen, i);

    if (c2 == carbon)
      continue;

    const Atom * ac2 = m.atomi(c2);

    if (3 != ac2->ncon())
      return 0;

    for (int j = 0; j < 3; ++j)
    {
      const Bond * b = ac2->item(j);

      const atom_number_t o = b->other(c2);

      if (b->is_double_bond() && 8 == m.atomic_number(o))
        return 1;
    }
  }

  return 0;
}

static int
looks_like_urea (const Molecule & m,
                 const atom_number_t carbon,
                 const resizable_array<const Bond *> & bonds_to_nitrogen)
{
  assert (2 == bonds_to_nitrogen.number_elements());

  if (1 == m.ncon(bonds_to_nitrogen[0]->other(carbon)))
    return 1;

  if (1 == m.ncon(bonds_to_nitrogen[1]->other(carbon)))
    return 1;

  if (carbonyl_present(m, carbon, bonds_to_nitrogen[0]->other(carbon)) ||
      carbonyl_present(m, carbon, bonds_to_nitrogen[1]->other(carbon)))
    return 0;

  return 1;

}

static int
possible_lactam_comparator (const Possible_Lactim_Lactam * pll1, const Possible_Lactim_Lactam * pll2)
{
  if (! pll1->is_ring() && ! pll2->is_ring())    // not in a ring, the ordering really does not matter, do something arbitrary. Implement sometime...
  {
    const atom_number_t c1 = pll1->carbon();
    const atom_number_t c2 = pll2->carbon();

    if (c1 < c2)        // this is of course invalid, but does not matte right now...
      return -1;
    else
      return 1;
  }

  if (pll1->is_ring() && ! pll2->is_ring())
    return -1;

  if (! pll1->is_ring() && pll2->is_ring())
    return 1;

// both are in a ring.

  if (pll1->fused_system_identifier() < pll2->fused_system_identifier())
    return 1;

  if (pll1->fused_system_identifier() > pll2->fused_system_identifier())
    return -1;

  if (pll1->ring_size() < pll2->ring_size())
    return 1;

  if (pll1->ring_size() > pll2->ring_size())
    return -1;

  return 0;
}


int
IWStandard_Current_Molecule::initialise (Molecule & m)
{
  _matoms = m.natoms();

  if (0 == _matoms)
    return 0;

  _atom = new const Atom * [_matoms];

  m.atoms(_atom);

  _atomic_number = new atomic_number_t[_matoms];
  _ncon = new int[_matoms];

  _nrings = m.nrings();

  _atom_is_aromatic = new_int(_matoms);

  if (_nrings)
  {
    _ring_membership = new int[_matoms];
    m.ring_membership(_ring_membership);
    _ring_is_aromatic = new_int(_nrings);
    _ring_nitrogen_count = new_int(_nrings);
    _ring_size = new int[_nrings];
    _ring_is_fused = new int[_nrings];
    _fsid = new_int(_matoms, -1);

    m.compute_aromaticity_if_needed();

    for (int i = 0; i < _nrings; i++)
    {
      const Ring * ri = m.ringi(i);

      Set_of_Atoms * s = new Set_of_Atoms(*ri);

      _rings.add(s);

      _ring_size[i] = ri->number_elements();
      _ring_is_fused[i] = ri->is_fused();

      if (ri->is_aromatic())
      {
        _ring_is_aromatic[i] = 1;
        ri->set_vector(_atom_is_aromatic, 1);
      }

      ri->set_vector(_fsid, ri->fused_system_identifier());
    }
  }
  else
  {
    _ring_membership = new_int(_matoms);
    _ring_is_aromatic = NULL;
  }

  atom_number_t first_singly_connected_oxygen = INVALID_ATOM_NUMBER;
  atom_number_t first_singly_connected_sulphur = INVALID_ATOM_NUMBER;

  for (int i = 0; i < _matoms; i++)
  {
    Atom * ai = const_cast<Atom *>(_atom[i]);

    if (! ai->valence_ok())
      _possible_valence_errors++;

    formal_charge_t fc = ai->formal_charge();

    _ncon[i] = ai->ncon();

    atomic_number_t z = ai->atomic_number();

    _atomic_number[i] = z;

    if (6 == z)
    {
      if (fc < 0)
        _cminus++;
    }
    else if (7 == z)
    {
      _nitrogens++;
      if (fc > 0)
        _nplus++;
    }
    else if (8 == z)
    {
      _oxygens++;
      if (ai->formal_charge() < 0)
        _ominus++;
      if (INVALID_ATOM_NUMBER == first_singly_connected_oxygen && 1 == _ncon[i])
        first_singly_connected_oxygen = i;
    }
    else if (16 == z)
    {
      _sulphur++;
      if (0 == fc)
        ;
      else if (fc < 0)
        _sminus++;
      else if (fc > 0)
        _splus++;

      if (INVALID_ATOM_NUMBER == first_singly_connected_sulphur && 1 == _ncon[i])
        first_singly_connected_sulphur = i;
    }
    else if (15 == z)
      _phosphorus++;
    else if (ai->element()->is_halogen() && 0 == _ncon[i])
      _isolated_halogen++;
    else if (ai->element()->is_metal())
    {
      if (0 == _ncon[i])
        _isolated_metal++;
      else if (1 == _ncon[i])
        _singly_connected_metal++;
    }
    else if (1 == z)
      _explicit_hydrogen_count++;

    if (0 == fc)
      continue;

    if (fc < 0)
      _nneg++;
    else if (fc > 0)
      _npos++;
  }

  if (_nitrogens >= 3)
  {
    for (int i = 0; i < _matoms; ++i)
    {
      if (6 != _atomic_number[i] || 3 != _ncon[i])
        continue;

      const Atom * a = _atom[i];

      if (4 != a->nbonds())
        continue;

      int ncount = 0;

      for (int j = 0; j < 3; ++j)
      {
        const atom_number_t n = a->other(i,j);
        if (7 != _atomic_number[n])
          break;
        ncount++;
      }

      if (3 == ncount)
        _possible_guanidine.add(i);
    }
  }

  if (_nitrogens >= 1 && _nrings)
  {
    for (int i = 0; i < _nrings; i++)
    {
      const Ring * ri = m.ringi(i);

      if (5 != ri->number_elements())
        continue;

      int nitrogens = 0;
      for (int j = 0; j < 5; j++)
      {
        atom_number_t k = ri->item(j);
        if (7 == _atomic_number[k])
          nitrogens++;
      }

      _ring_nitrogen_count[i] = nitrogens;
    }
  }

// Only one remaining is lactam lactim

  if ((INVALID_ATOM_NUMBER == first_singly_connected_oxygen && 
       INVALID_ATOM_NUMBER == first_singly_connected_sulphur) || 0 == _nitrogens)
    return 1;

  int istart = first_singly_connected_oxygen;
  if (INVALID_ATOM_NUMBER == first_singly_connected_oxygen)
    istart = first_singly_connected_sulphur;
  else if (first_singly_connected_sulphur >= 0 && first_singly_connected_sulphur < first_singly_connected_oxygen)
    istart = first_singly_connected_sulphur;

  for (int i = istart; i < _matoms; ++i)
  {
    if (1 != _ncon[i])
      continue;

    if (8 == _atomic_number[i])
      ;
    else if (16 == _atomic_number[i])
      ;
    else
      continue;

    const Bond * boc = _atom[i]->item(0);

    if (! boc->is_single_bond())    // already in correct form
      continue;

    const atom_number_t c = boc->other(i);

    if (6 != _atomic_number[c])
      continue;

    if (3 != _ncon[c])
      continue;

    const Atom * ac = _atom[c];

    if (4 != ac->nbonds())    // must be unsaturated
      continue;

    resizable_array<const Bond *> bonds_to_nitrogen;
    const int attached_nitrogens = gather_bonds_to_nitrogen_atoms(m, c, i, bonds_to_nitrogen);

    if (0 == attached_nitrogens)
      continue;

    if (2 == attached_nitrogens && boc->is_double_bond() && looks_like_urea(m, c, bonds_to_nitrogen))
      continue;

    const Bond * bn1 = NULL;
    const Bond * bn2 = NULL;

    for (int j = 0; j < attached_nitrogens; ++j)
    {
      const Bond * bcn = bonds_to_nitrogen[j];

      if (bcn->is_double_bond() && bcn->part_of_cis_trans_grouping()) 
        continue;

      const atom_number_t n = bcn->other(c);

      if (2 != _ncon[n])
        continue;

      if (single_bond_to_oxygen(m, n))
        continue;

      if (NULL == bn1)
        bn1 = bcn;
      else
        bn2 = bcn;
    }

    if (NULL == bn1)
      continue;

    if (NULL != bn2 && bn1->is_single_bond())     // make sure bn1 is the double bond (if present)
      std::swap(bn1, bn2);

    const atom_number_t n1 = bn1->other(c);

    Possible_Lactim_Lactam * p = NULL;

    if (boc->is_double_bond())    // already correct - actually, will not happen
    {
      p = new Possible_Lactim_Lactam(i, c, n1);
      p->set_lactam_form(1);
    }
    else             // single bond from oxygen to carbon, may need changing
    {
      p = new Possible_Lactim_Lactam(i, c, n1);
    }

    if (NULL != bn2)
      p->set_second_nitrogen(bn2->other(c));

    p->set_total_nitrogen_attachments(attached_nitrogens);

    if (_ring_membership[c] > 0)
    {
      p->set_is_ring(1);
      p->set_fused_system_identifier(_fsid[c]);

      if (_atom_is_aromatic[c])
      {
        p->set_aromatic(1);
        p->set_fused_system_size(m.fused_system_size(c));
      }
    }

    _possible_lactam.add(p);
  }

  for (int i = 0; i < _possible_lactam.number_elements(); ++i)
  {
    _possible_lactam[i]->discern_alpha_nitrogen(m);
  }

//remove_aromatic_lactims_with_alpha_nitrogen(m, _possible_lactam);

  if (1 == _possible_lactam.number_elements())
    return 1;

  remove_two_lactims_sharing_a_single_nitrogen(_possible_lactam);

  const int n = _possible_lactam.number_elements();

  for (int i = 0; i < n; ++i)
  {
    Possible_Lactim_Lactam * pi = _possible_lactam[i];

    if (! pi->aromatic())
      continue;

    for (int j = i + 1; j < n; ++j)
    {
      Possible_Lactim_Lactam * pj = _possible_lactam[j];

      if (! pj->aromatic())
        continue;

      if (pi->fused_system_identifier() == pj->fused_system_identifier())
      {
        pi->increment_lactims_in_fused_system();
        pj->increment_lactims_in_fused_system();
      }
    }
  }

  _possible_lactam.iwqsort(possible_lactam_comparator);

  return 1;
}

int
Chemical_Standardisation::_processing_needed (const IWStandard_Current_Molecule & current_molecule_data) const
{
  if (current_molecule_data.nneg() || current_molecule_data.npos() || 
     (current_molecule_data.singly_connected_metal() && _transform_covalent_metals.active()) ||
     (current_molecule_data.isolated_metal() && _transform_single_atom_ions.active()) ||
     (current_molecule_data.isolated_halogen() && _transform_single_atom_ions.active()) ||
      current_molecule_data.possible_guanidine().number_elements() ||
      current_molecule_data.explicit_hydrogen_count() ||
      current_molecule_data.aromatic_rings_with_multiple_nitrogens() > 0 ||
      (current_molecule_data.nitrogens()> 1 && _transform_back_to_nplus_nminus.active()) ||
      (current_molecule_data.nitrogens() > 1 && current_molecule_data.oxygens() > 0 && _transform_misdrawn_urea.active()) ||
      current_molecule_data.possible_lactam().number_elements() > 0 ||
      current_molecule_data.possible_valence_errors() ||
      (current_molecule_data.nitrogens() && current_molecule_data.oxygens() && _transform_isoxazole.active()) ||
      (current_molecule_data.sulphur() && _transform_amino_thiazole.active()))

    return 1;

  return 0;
}

int
IWStandard_Current_Molecule::aromatic_rings_with_multiple_nitrogens () const
{
  if (0 == _nrings)
    return 0;

  int rc = 0;

  for (int i = 0; i < _nrings; i++)
  {
    if (! _ring_is_aromatic[i])
      continue;

    if (_ring_nitrogen_count[i] > 1)
      rc++;
  }

  return rc;
}

int
IWStandard_Current_Molecule::ring_number_containing_atom (const atom_number_t a) const
{
  const int nr = _rings.number_elements();

  for (int i = 0; i < nr; ++i)
  {
    if (_rings[i]->contains(a))
      return i;
  }

  return -1;
}


const Set_of_Atoms *
IWStandard_Current_Molecule::ring_containing_atom (const atom_number_t a) const
{
  const int nr = _rings.number_elements();

  for (int i = 0; i < nr; ++i)
  {
    if (_rings[i]->contains(a))
      return _rings[i];
  }

  return NULL;
}

const Set_of_Atoms *
IWStandard_Current_Molecule::ringi (const int s) const
{
  return _rings[s];
}

int
IWStandard_Current_Molecule::remove_possible_guanidine (const atom_number_t c)
{
  return _possible_guanidine.remove_first(c);
}

int
Chemical_Standardisation::_do_amino_thiazole(Molecule & m,
                                             int * atom_already_changed,
                                             IWStandard_Current_Molecule & current_molecule_data)
{
  if (0 == current_molecule_data.sulphur())
    return 0;

  const int nr = current_molecule_data.nrings();
  if (0 == nr)
    return 0;

  int rc = 0;

  for (int i = 0; i < nr; ++i)
  {
    const Set_of_Atoms * ri = current_molecule_data.ringi(i);

    if (ri->number_elements() < 5)
     continue;

   if (ri->number_elements() > 5)
     break;

   if (! current_molecule_data.ring_is_aromatic()[i])
     continue;

//  cerr << "Possible amino thiazole " << *ri << endl;

    rc += _do_amino_thiazole(m, *ri, atom_already_changed, current_molecule_data);
  }

  if (0 == rc)
    return 0;

  _transform_amino_thiazole.extra(rc);

  if (_verbose)
    cerr << "Transformed " << rc << " amino thiazole\n";

  if (_append_string_depending_on_what_changed)
    _append_to_changed_molecules << " STD:amino_thiazole";

  return rc;
}

int
Chemical_Standardisation::_do_amino_thiazole (Molecule & m,
                                              const Set_of_Atoms & r,
                                              const int * atom_already_changed,
                                              IWStandard_Current_Molecule & current_molecule_data)
{
  int s = -1;
  int n = -1;

  int c_exo_n = -1;     // the carbon with a =N outside the ring
  atom_number_t n_exo_c = INVALID_ATOM_NUMBER;    // the =N atom

  assert (5 == r.number_elements());

  for (int i = 0; i < 5; ++i)
  {
    const atom_number_t j = r[i];

    const atomic_number_t zj = current_molecule_data.atomic_number()[j];
    if (8 == zj || 16 == zj)     // try oxygen too
    {
      if (s >= 0)   // can only be one in the ring
        return 0;

      s = i;
      continue;
    }

    if (7 == zj)
    {
      if (n >= 0)   // can only be one in the ring
        return 0;

      if (2 != current_molecule_data.ncon()[j])
        return 0;

      n = i;
      continue;
    }

    if (6 != zj)    // ring can contain only S C and N
      return 0;

    if (2 == current_molecule_data.ncon()[j])   // if it is going to have an exocyclic Nitrogen, must have 3 connections
      continue;

    const Atom * aj = m.atomi(j);

    for (int k = 0; k < 3; ++k)
    {
      const Bond * b = aj->item(k);

      if (! b->is_double_bond())
        continue;

      const atom_number_t x = b->other(j);
      if (7 != current_molecule_data.atomic_number()[x])
        continue;

      if (INVALID_ATOM_NUMBER != n_exo_c)
        return 0;

      c_exo_n = i;
      n_exo_c = x;
    }
  }

  if (INVALID_ATOM_NUMBER == n_exo_c || s < 0 || n < 0)
    return 0;

// The Nitrogen and carbon=N must be adjacent in the ring

  if ((c_exo_n+1)%5 == n)
    ;
  else if ((n+1)%5 == c_exo_n)
    ;
  else
    return 0;

  if ((s+1)%5 == c_exo_n)
    ;
  else if ((c_exo_n+1)%5 == s)
    ;
  else
    return 0;

  m.set_bond_type_between_atoms(n_exo_c, r[c_exo_n], SINGLE_BOND);
  m.set_bond_type_between_atoms(r[c_exo_n], r[n], DOUBLE_BOND);
  m.recompute_implicit_hydrogens(n_exo_c);
  m.recompute_implicit_hydrogens(r[n]);

  return 1;
}

template resizable_array_base<Possible_Lactim_Lactam*>::~resizable_array_base();
template resizable_array<Possible_Lactim_Lactam*>::resizable_array();
