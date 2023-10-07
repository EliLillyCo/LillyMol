/*
  Main purpose is to generate fingerprint information around the raction core, defined by the
  changing atoms.
*/

#include <algorithm>
#include <cctype>
#include <iostream>
#include <memory>
#include <string>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/md5.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/rxn_file.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/extended_connectivity_fp.h"
#include "Molecule_Tools/linear_path_fingerprint.h"

using std::cerr;
using std::endl;

const char * prog_name = nullptr;

static int verbose = 0;

static int reactions_read = 0;

static Chemical_Standardisation chemical_standardisation;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString product_smiles_tag("PSMI<");

static int take_first_of_multiple_reagents = 0;
static int took_first_of_multiple_reagents = 0;
static int skip_reactions_with_multiple_ragents = 0;

static int reagents_with_multiple_reagents_skipped = 0;

static extending_resizable_array<int> acc_changing_atoms;

static int nr = 0;
static int * radius = nullptr;
static IWString * reagent_tag = nullptr;
static IWString * product_tag = nullptr;

static Atom_Typing_Specification atom_typing_specification;

static int default_nbits_dense = 2048;
static int default_nbits_sparse = 200000;

static int isotope = 0;

static int discard_existing_atom_map = 0;

static int expand_across_double_bonds = 0;

static int fingerprint_changing_product_atoms = 0;

static IWString difference_fingerprint_tag;

static IWString_and_File_Descriptor stream_for_changing_atoms;
static IWString_and_File_Descriptor stream_for_reagent_count;
static IWString_and_File_Descriptor stream_for_multiple_reagents;

static IWString_and_File_Descriptor stream_for_reaction_smiles;

static int write_agent = 1;

static std::ofstream stream_for_reactions;

static int write_failed_changing_atom_counts_to_wfile = 0;

static int append_changing_atoms_count_to_name = 0;

static IWString changing_atom_count_tag;

static int any_changing_bond_means_a_changing_atoms = 0;
static int only_consider_largest_reagent_fragment = 0;
static int include_changing_bonds_in_changing_atom_count = 0;

static resizable_array<int> changing_atoms_should_be;

static int reactions_with_unexpected_changing_atom_counts = 0;
static int reactions_with_corrected_atom_counts = 0;

/*
  No, we cannot do this. If we remove atoms from a reagent,
  that will mess up the atom map. We would need to build
  a remove atom method for the RXN_File object
*/

static int reduce_to_largest_fragment = 0;

static resizable_array_p<Substructure_Query> changing_atoms_must_match;

static int reactions_rejected_for_changing_atoms_must_match = 0;

static IW_STL_Hash_Map_int changing_atom_smarts[100];

static int discard_reaction_unless_all_atoms_mapped = 0;
static int reactions_discarded_for_incomplete_atom_map = 0;

static int remove_duplicate_reagents_atom_maps_scrambled = 0;

static extending_resizable_array<int> unmapped_atom_count;

static Accumulator<double> acc_unmapped_atom_fraction;

static int write_not_full_mapped_reactins_to_reaction_stream = 0;

static int reaction_name_is_file_name = 0;

static int strip_rxn_from_file_name = 0;

static int reaction_name_is_path_name = 0;

static int use_first_token_of_name = 0;

static int convert_products_to_agents_when_possible = 0;

static LFP::Linear_Fingerprint_Defaults lfpd;

static int ignore_bad_reactions = 0;

static int bad_reactions_ignored = 0;

static IWString_and_File_Descriptor stream_for_bad_reactions;

static int discard_reactions_with_isotopes = 0;
static int reactions_discarded_for_isotopic_atoms = 0;

static int discard_reactions_with_no_changing_atoms = 0;
static int reactions_with_no_changing_atoms_discarded = 0;

static int input_is_reaction_smiles = 0;

static int discern_changing_atoms_only_in_first_fragment = 0;

static int set_extra_bits_for_changing_atoms = 0;

static int reactions_with_no_reagent_atoms_in_products = 0;

static int reactions_containing_duplicate_atom_map_numbers = 0;

static int remove_cis_trans_bonding = 0;

static int remove_non_participating_fragments = 0;

static int only_create_query_from_first_reagent = 1;

static int consider_aromatic_bonds = 0;

class Columnar_Output_Specifications
{
  private:
    int _active;
    int _produce_NC_fingerprints;
    int _produce_md5_fingerprints;
    IWString _output_separator;
    IWString _bit_separator;

  public:
    Columnar_Output_Specifications();

    int build(const Command_Line & cl, const char flag, const int verbose);

    int display_help_message(const char flag, std::ostream & output) const;

    const IWString & output_separator() const { return _output_separator;}
    const IWString & bit_separator() const { return _bit_separator;}
    int active() const { return _active;}
    int produce_NC_fingerprints() const { return _produce_NC_fingerprints;}
    int produce_md5_fingerprints() const { return _produce_md5_fingerprints;}
};

Columnar_Output_Specifications::Columnar_Output_Specifications()
{
  _active = 0;
  _produce_NC_fingerprints = 0;
  _produce_md5_fingerprints = 0;
  _output_separator = ' ';
  _bit_separator = ',';

  return;
}

static void
set_separator(const const_IWSubstring & k,
              IWString & sep)
{
  if ("space" == k)
    sep = ' ';
  else if ("tab" == k)
    sep = '\t';
  else if ("comma" == k)
    sep = ',';
  else if ("semic" == k)
   sep = ';';
  else if ("colon" == k)
   sep = ':';
  else if ("squote" == k)
   sep = '\'';
  else if ("dquote" == k)
   sep = '"';
  else if ("vbar" == k)
    sep = '|';
  else
    sep = k;

  return;
}

int
Columnar_Output_Specifications::build(const Command_Line & cl,
                                      const char flag,
                                      const int verbose)
{
  const_IWSubstring k;
  cerr << "FLAG " << flag << endl;
  for (int i = 0; cl.value(flag, k, i); ++i)
  {
    cerr << "Processing '" << k << "'\n";
    if ("MD5" == k || "md5" == k)
    {
      _produce_md5_fingerprints = 1;
      if (verbose)
        cerr << "Will produce md5 hashes of fingerprints\n";
    }
    else if ("NC" == k || "nc" == k)
    {
      _produce_NC_fingerprints = 1;
      if (verbose)
        cerr << "Will produce NC fingerprints\n";
    }
    else if (k.starts_with("lin"))
    {
      _produce_NC_fingerprints = 0;

      if (verbose)
        cerr << "Will produce linear path fingerprints\n";
    }
    else if (k.starts_with("sep="))
    {
      k.remove_leading_chars(4);
      set_separator(k, _output_separator);
    }
    else if (k.starts_with("bsep="))
    {
      k.remove_leading_chars(5);
      set_separator(k, _bit_separator);
    }
    else if ("def" == k)
      ;
    else if ("help" == k)
    {
      display_help_message(flag, cerr);
    }
    else
    {
      cerr << "Columnar_Output_Specifications::build:unrecognised directive '" << k << "'\n";
      return 0;
    }
  }

  _active = 1;

  if (verbose)
    cerr << "Will produce columnar output\n";

  return 1;
}

int
Columnar_Output_Specifications::display_help_message(const char flag, std::ostream & output) const
{
  output << "Columnar Output via -" << flag << endl;
  output << " -" << flag << " MD5          produce MD5 of fingerprints\n";
  output << " -" << flag << " NC           produce non colliding fingerprints\n";
  output << " -" << flag << " sep='x'      output separator, default space\n";
  output << " -" << flag << " bsep='x'     separator between bits separator, default comma\n";
  output << " -" << flag << " def          default behaviour for all\n";

  return 1;
}

static class Columnar_Output_Specifications columnar_output_specifications;

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
  cerr << "  -i            input is reaction smiles\n";
  cerr << "  -r <rad>      radius from changing atoms to fingerprint\n";
  cerr << "  -J <tag>      tag for each radius (FP, NC), one per each -r option. NCEC means circular\n";
  cerr << "  -K <...>      produce columnar output\n";
  cerr << "  -m first      when multiple reagents present, process first\n";
  cerr << "  -m skip       when multiple reagents present, skip\n";
  cerr << "  -d            expand shell to include doubly bonded atoms\n";
  cerr << "  -h            an atom is considered changing if any new bond appears\n";
  cerr << "  -u            fingerprint changing product atoms also\n";
  cerr << "  -P <type>     atom typing specification\n";
  cerr << "  -x            discard existing atom mapping in input files\n";
  cerr << "  -j .          discard reaction if atom mapping does not map all atoms\n";
  cerr << "  -C <fname>    write changed atom counts to <fname> (csv)\n";
  cerr << "  -a            append changing atoms count to name\n";
  cerr << "  -c <n>        if the computed number of changing atoms is NOT <n>, recompute atom map\n";
  cerr << "  -I <iso>      mark changing atoms with <isotope>\n";
  cerr << "  -D <fname>    write reagent counts to <fname>\n";
  cerr << "  -W <fname>    write reactions to <fname>\n";
  cerr << "  -w            write reactions that cannot be re-mapped to -W file\n";
  cerr << "  -s <smarts>   specify what the changing atoms should be\n";
  cerr << "  -q <query>    specify what the changing atoms should be\n";
  cerr << "  -U p          use the reaction path name as the reaction name\n";
  cerr << "  -U f          use the reaction file name as the reaction name\n";
  cerr << "  -U nrxn       strip '.rxn' from any file name generated\n";
  cerr << "  -f            truncate reaction names to the first token\n";
  cerr << "  -B            remove duplicate reagents, even if atom maps scrambled\n";
  cerr << "  -y            do NOT include ring information in fingerprints\n";
  cerr << "  -X <fname>    write reaction smiles\n";
  cerr << "  -e            only look at reagent atoms in the largest fragment when discerning changing atoms\n";
  cerr << "  -M ...        miscellaneous options, enter '-M help' for info\n";
  cerr << "  -F <fname>    ignore otherwise bad reactions and write them to <fname>\n";
//cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess(Molecule & m)
{
  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment_carefully();

  return;
}

typedef unsigned int atype_t;

/*
  Avoid passing too many arguments between functions by putting a bunch of things in a class.
*/

class RXNFP_Temporary_Arrays
{
  private:
    int * _changed;
    atype_t * _atype;
    int * _tmp;

    atype_t * _atsave;

    const int _natoms;

  public:
    RXNFP_Temporary_Arrays(const int n, const int x);
    ~RXNFP_Temporary_Arrays();

    int * changed () const { return _changed;}
    atype_t * atype () const { return _atype;}
    int * tmp () const { return _tmp;}
    atype_t * atsave () const { return _atsave;}

    void save_and_change_atom_types_of_changed_atoms(const int offset);
    void restore_changed_atom_types();
    void dump(std::ostream& output) const;  // for debugging
};

RXNFP_Temporary_Arrays::RXNFP_Temporary_Arrays(const int n, const int x) : _natoms(n)
{
  _changed = new_int(_natoms);

  _atype = new atype_t[_natoms + _natoms];
  _atsave = _atype + _natoms;

  _tmp = new_int(x, 1);
}

RXNFP_Temporary_Arrays::~RXNFP_Temporary_Arrays()
{
  delete [] _changed;
  delete [] _atype;
  delete [] _tmp;

  return;
}

void
RXNFP_Temporary_Arrays::dump(std::ostream& output) const
{
  output << "\n#Atoms: " << _natoms<< "\n";
  output << "#\tatype\tchanged\ttmp\n";
    
  for(int i=0; i < _natoms ; ++i)
  {
      output  << i << "\t " << _atype[i] << "\t" << _changed[i]<<  "\t" << _tmp[i]<<"\n";
  }  
}

/*
  When doing extra fingerprint bits for the changing atoms we need to save the values
*/

void
RXNFP_Temporary_Arrays::save_and_change_atom_types_of_changed_atoms(const int offset)
{
  for (int i = 0; i < _natoms; ++i)
  {
    if (0 == _changed[i])
      continue;

    _atsave[i] = _atype[i];
    _atype[i] += offset;
  }

  return;
}

void
RXNFP_Temporary_Arrays::restore_changed_atom_types()
{
  for (int i = 0; i < _natoms; ++i)
  {
    if (0 == _changed[i])
      continue;

    _atype[i] = _atsave[i];
  }

  return;
}

static int
do_apply_isotopes(Molecule & m, 
                  const int * changed,
                  const int isotope)
{
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (changed[i])
      m.set_isotope(i, isotope);
  }

  return 1;
}

void 
transfer_to_sparse(const int * b,
                   const int n,
                   Sparse_Fingerprint_Creator & sfc)
{
  for (int i = 0; i < n; ++i)
  {
    if (0 == b[i])
      continue;

    sfc.hit_bit(i, b[i]);
  }
}

static int
identify_atoms_in_fingerprint(Molecule & m,
                              const int rad,
                              const int * changed,
                              int * in_fingerprint)
{
  const int matoms = m.natoms();

  std::copy_n(changed, matoms, in_fingerprint);    // by default, all changing atoms in fp

  if (rad > 0)
  {
    for (int i = 0; i < matoms; ++i)
    {
      if (0 == changed[i])    // we look from each of the changing atoms
        continue;

      const int ifrag = m.fragment_membership(i);

      for (int j = 0; j < matoms; ++j)   // check all the others
      {
        if (in_fingerprint[j])
          continue;

        if (ifrag != m.fragment_membership(j))
          continue;

        if (m.bonds_between(i, j) <= rad)
          in_fingerprint[j] = 1;
      }
    }
  }

  if (expand_across_double_bonds)
  {
    m.compute_aromaticity_if_needed();

    Set_of_Atoms extras;

    for (int i = 0; i < matoms; ++i)
    {
      if (0 == in_fingerprint[i])
        continue;

      const Atom * a = m.atomi(i);

      const int acon = a->ncon();

      if (acon == m.nbonds(i))    // fully saturated, no double bonds here
        continue;

      for (int j = 0; j < acon; ++j)
      {
        const Bond * b = a->item(j);

        if (b->is_single_bond() || b->is_aromatic())
          continue;

        const int k = b->other(i);

        if (in_fingerprint[k])
          continue;

        extras.add(k);
      }
    }

//  cerr << "Identified " << extras.size() << " doubly bonded extra in " << m.name() << ", radius " << rad << endl;
    extras.set_vector(in_fingerprint, 1);
  }

  return 1;
}

static int
produce_product_fingerprint(ISIS_RXN_FILE_Molecule & p,
                            const int rad,
                            int * changed,
                            int * tmp,
                            int * atype,
                            Sparse_Fingerprint_Creator & sfc)
{
  const int matoms = p.natoms();

  std::copy_n(changed, matoms, tmp);

  identify_atoms_in_fingerprint(p, rad, changed, tmp);

#ifdef DEBUG_PRODUCE_PRODUCT_FP
  for (int i = 0; i < matoms; ++i)
  {
    if (tmp[i])
      cerr << "Atom " << i << " in fp " << p.smarts_equivalent_for_atom(i) << " atype " << atype[i] << " rad " << rad << endl;
  }
#endif

  set_iwmfingerprint_nbits(default_nbits_sparse);

  IWMFingerprint fp;

  fp.construct_fingerprint(p, atype, tmp);

  transfer_to_sparse(fp.vector(), iwmfingerprint_nbits(), sfc);

  return 1;
}


static void
regenerate_with_changed_atom_offsets(ISIS_RXN_FILE_Molecule & r,
                                     EC_Fingerprint_Generator & ecfp,
                                     RXNFP_Temporary_Arrays & rxnfpta,
                                     Sparse_Fingerprint_Creator & sfc)
{
  rxnfpta.save_and_change_atom_types_of_changed_atoms(1702);

  ecfp.generate_fingerprint(r, rxnfpta.atype(), rxnfpta.tmp(), 1, sfc);

  rxnfpta.restore_changed_atom_types();

  return;
}

static int
produce_ec_fingerprint(ISIS_RXN_FILE_Molecule & r,
                       const int rad,
                       const IWString & tag,
                       RXNFP_Temporary_Arrays & rxnfpta,
                       IWString_and_File_Descriptor & output)
{
//cerr << "EC FP's, md5 " << columnar_output_specifications.produce_md5_fingerprints() << endl;
//cerr << " bit_separator " << columnar_output_specifications.bit_separator() << endl;

  EC_Fingerprint_Generator ecfp;
  ecfp.set_differentiate_ring_and_chain_bonds(1);

  Sparse_Fingerprint_Creator sfc;
  ecfp.generate_fingerprint(r, rxnfpta.atype(), rxnfpta.tmp(), 1, sfc);

  if (rad > 0 && set_extra_bits_for_changing_atoms)
    regenerate_with_changed_atom_offsets(r, ecfp, rxnfpta, sfc);

  if (columnar_output_specifications.active())
  {
    output << columnar_output_specifications.output_separator();
    if (columnar_output_specifications.produce_md5_fingerprints())
      sfc.write_as_md5_sum(output);
    else
      sfc.write_as_feature_count(columnar_output_specifications.bit_separator()[0], output);
  }
  else
  {
    IWString tmp;

    sfc.daylight_ascii_form_with_counts_encoded(tmp);

    output << tag << tmp << ">\n";
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

template <typename G, typename B>
void
regenerate_with_changed_atom_offsets(ISIS_RXN_FILE_Molecule & r,
                                     G & lfpc,
                                     RXNFP_Temporary_Arrays & rxnfpta,
                                     B & bits)
{
  rxnfpta.save_and_change_atom_types_of_changed_atoms(992);

  lfpc.process(r, rxnfpta.atype(), rxnfpta.tmp(), bits);

  rxnfpta.restore_changed_atom_types();

  return;
}

/*
  Linear is more complicated because it can do either fixed width or sparse
*/

static int
produce_linear_fingerprint(ISIS_RXN_FILE_Molecule & r,
                           const int rad,
                           const IWString & tag,
                           RXNFP_Temporary_Arrays & rxnfpta,
                           IWString_and_File_Descriptor & output)
{
  LFP::Linear_Fingerprint_Creator lfpc;
  lfpc.set_defaults(&lfpd);

  if (tag.starts_with("FP"))
  {
    LFP::Fixed_Width_Fingerprint_Bits fwfpb(default_nbits_dense);
    lfpc.process(r, rxnfpta.atype(), rxnfpta.tmp(), fwfpb);

    if (set_extra_bits_for_changing_atoms)
      regenerate_with_changed_atom_offsets(r, lfpc, rxnfpta, fwfpb);
      
//cerr << "Fingerprint contains " << fwfpb.nset() << " bits " << r.smiles() << endl;
    if (columnar_output_specifications.active())
    {
      output << columnar_output_specifications.output_separator();
      if (columnar_output_specifications.produce_md5_fingerprints())
        fwfpb.write_as_md5_sum(output);
      else
        fwfpb.write_as_feature_count(columnar_output_specifications.bit_separator(), output);
    }
    else
      fwfpb.append_fingerprint(tag, output);
  }
  else
  {
    LFP::Sparse_Fingerprint_Bits sfpb;
    lfpc.process(r, rxnfpta.atype(), rxnfpta.tmp(), sfpb);

    if (set_extra_bits_for_changing_atoms)
      regenerate_with_changed_atom_offsets(r, lfpc, rxnfpta, sfpb);
      
    if (columnar_output_specifications.active())
    {
      output << columnar_output_specifications.output_separator();
      if (columnar_output_specifications.produce_md5_fingerprints())
        sfpb.write_as_md5_sum(output);
      else
        sfpb.write_as_feature_count(columnar_output_specifications.bit_separator(), output);
    }
    else
      sfpb.append_fingerprint(tag, output);
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
do_generate_difference_fingerprint(RXN_File & r,
                                   RXNFP_Temporary_Arrays & rxnfpta,
                                   const IWString & tag,
                                   IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfc;

  Product_Atom_Types pat;

  pat.initialise(r, atom_typing_specification);

  const int * changed = rxnfpta.changed();

  r.reaction_fingerprint(changed, rxnfpta.atype(), pat, 0, sfc);

  if (columnar_output_specifications.active())
  {
    output << columnar_output_specifications.output_separator();

    if (columnar_output_specifications.produce_md5_fingerprints())
      unordered_map_to_md5(sfc.bits_found(), output);
    else
      sfc.write_as_feature_count(columnar_output_specifications.bit_separator()[0], output);
  }
  else
  {
    IWString tmp;

    sfc.daylight_ascii_form_with_counts_encoded(tmp);

    output << tag << tmp << ">\n";
  }

  return 1;
}

static int
produce_fingerprint(ISIS_RXN_FILE_Molecule & r,  
                    const int rad, 
                    const IWString & tag,
                    RXNFP_Temporary_Arrays & rxnfpta,
                    IWString_and_File_Descriptor & output)
{
  const int matoms = r.natoms();

  int * changed = rxnfpta.changed();
  int * tmp = rxnfpta.tmp();

  std::copy_n(changed, matoms, tmp);

  identify_atoms_in_fingerprint(r, rad, changed, tmp);

//#define DEBUG_PRODUCE_REAGENT_FP
#ifdef DEBUG_PRODUCE_REAGENT_FP
  Molecule x(r);
  x.reset_all_atom_map_numbers();
  for (int i = 0; i < matoms; ++i)
  {
    if (0 == tmp[i])
      continue;

    cerr << "Atom " << i << " in fp " << r.smarts_equivalent_for_atom(i) << " atype " << rxnfpta.atype()[i] << " rad " << rad << endl;

    if (changed[i])
      x.set_isotope(i, 1);
    else
      x.set_isotope(i, 2);
  }
  cerr << x.smiles() << " RAD " << rad <<  " tag " << tag << endl;
#endif

//cerr << "NC " << columnar_output_specifications.produce_NC_fingerprints() << endl;

  if (tag.starts_with("NCEC") || columnar_output_specifications.produce_NC_fingerprints())
    return produce_ec_fingerprint(r, rad, tag, rxnfpta, output);
  else     // assume linear
    return produce_linear_fingerprint(r, rad, tag, rxnfpta, output);
}

/*
  this has really not bee used, not sure how useful it is...
*/

static int
do_fingerprint_changing_product_atoms(RXN_File & rxn,
                                      IWString_and_File_Descriptor & output)
{
  const int np = rxn.number_products();

  if (0 == np)   // hard to imagine
    return 0;

  int max_atoms = rxn.product(0).natoms();

  for (int i = 1; i < np; ++i)
  {
    const int t = rxn.product(i).natoms();
    if (t > max_atoms)
      max_atoms = t;
  }

  if (0 == max_atoms)    // hard to imagine
    return 0;

  int * changed = new_int(max_atoms + max_atoms + max_atoms); std::unique_ptr<int[]> free_changed(changed);
  int * tmp     = changed + max_atoms;

  int * atype;
  if (atom_typing_specification.active())
    atype = changed + max_atoms + max_atoms;
  else
    atype = nullptr;

  Changing_Atom_Conditions cac;
  if (any_changing_bond_means_a_changing_atoms)
    cac.set_is_changing_if_different_neighbours(1);
  if (only_consider_largest_reagent_fragment)
    cac.set_only_consider_largest_reagent_fragment(1);
  if (discern_changing_atoms_only_in_first_fragment)
    cac.set_discern_changing_atoms_only_in_first_fragment(1);
  cac.set_consider_aromatic_bonds(consider_aromatic_bonds);
  if (consider_aromatic_bonds)
  {
    rxn.set_preserve_kekule_forms(0);
    rxn.set_aromatic_bonds_lose_kekule_identity(1);
  }

  cac.set_ignore_lost_atom_if_isolated(1);

  Sparse_Fingerprint_Creator * sfc = new Sparse_Fingerprint_Creator[nr]; std::unique_ptr<Sparse_Fingerprint_Creator[]> free_sfc(sfc);

  for (int p = 0; p < np; ++p)
  {
    output << product_smiles_tag << rxn.product(p).smiles() << ">\n";

    std::fill_n(changed, max_atoms, 0);

    int c;
    if (atom_typing_specification.active())
      c = rxn.identify_atoms_changing_product(p, atom_typing_specification, changed, cac);
    else
      c = rxn.identify_atoms_changing_product(p, changed, cac);
    (void) c;
//  cerr << " found " << c << " changing atoms in product\n";

    if (atom_typing_specification.active())
    {
      std::fill_n(atype, max_atoms, 0);          // is this necessary?
      if (! atom_typing_specification.assign_atom_types(rxn.product(p), atype))
      {
        cerr << "Cannot assign atom types '" << rxn.name() << "'\n";
        return 0;
      }

      for (int j = 0; j < nr; ++j)
      {
        if (! produce_product_fingerprint(rxn.product(p), radius[j], changed, tmp, atype, sfc[j]))
          return 0;
      }
    }
    else
    {
      for (int j = 0; j < nr; ++j)
      {
        if (! produce_product_fingerprint(rxn.product(p), radius[j], changed, tmp, atype, sfc[j]))
          return 0;
      }
    }
  }

  for (int r = 0; r < nr; ++r)
  {
    IWString tmp;

    sfc[r].daylight_ascii_form_with_counts_encoded(tmp);

    output << product_tag[r] << tmp << ">\n";
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static void
translate_strange_characters(IWString & s)
{
  const int n = s.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const char c = s[i];
    if (isalnum(c))
      continue;

    if ('_' == c)
      continue;
    if ('.' == c)
      continue;
    if ('-' == c)
      continue;

    if (use_first_token_of_name)
    {
//    cerr << "use_first_token_of_name:truncating at '" << c << "', i = " << i << endl;
      s.iwtruncate(i);
      return;
    }

    s[i] = '_';
  }

  return;
}

static int
generate_externally_visible_reaction_name(const RXN_File & rxn,
                                          IWString & s,
                                          bool do_tranlate = false)
{
  if (reaction_name_is_file_name)
  {
    s = rxn.fname();
    const int last_dir_separator = s.rindex('/');
    if (last_dir_separator > 0)
      s.remove_leading_chars(last_dir_separator + 1);
    
  }
  else if (reaction_name_is_path_name)
    s = rxn.fname();
  else
    s = rxn.name();

  if (strip_rxn_from_file_name && s.ends_with(".rxn"))
    s.chop(4);

  if (do_tranlate)
    translate_strange_characters(s);

  return 1;
}

static void
add_changing_atom_smarts_to_hash(ISIS_RXN_FILE_Molecule & r,
                                 const int * changed,
                                 IW_STL_Hash_Map_int & h)
{
  const int matoms = r.natoms();

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == changed[i])
      continue;

    IWString s = r.smarts_equivalent_for_atom(i);

    h[s]++;
  }

  return;
}

static int
all_changing_atoms_in_embedding(const int matoms,
                                const int * changed,
                                const Set_of_Atoms * e)
{
  for (int i = 0; i < matoms; ++i)
  {
    if (0 == changed[i])
      continue;

    if (! e->contains(i))
      return 0;
  }

  return 1;
}

/*
  We have the changing atoms. We also have some queries that say that the changing atoms
  must consist of certain motifs.
  Make sure that the changed atoms match a query match
*/

static int
changing_atoms_match(ISIS_RXN_FILE_Molecule & r,
                     const int c,
                     const int * changed,
                     resizable_array_p<Substructure_Query> & changing_atoms_must_match)
{
  const int matoms = r.natoms();

  Molecule_to_Match target(&r);

  for (int i = 0; i < changing_atoms_must_match.number_elements(); ++i)
  {
    Substructure_Results sresults;

    const int nhits = changing_atoms_must_match[i]->substructure_search(target, sresults);
//  cerr << nhits << " matches to changing atoms must match query " << i << endl;
    if (0 == nhits)
      continue;

    for (int j = 0; j < nhits; ++j)
    {
      const Set_of_Atoms * e = sresults.embedding(j);

      if (c != e->number_elements())   // cannot match
        continue;

      if (all_changing_atoms_in_embedding(matoms, changed, e))
        return 1;
    }
  }

  return 0;
}

static int
append_changing_atom_count_fingerprint(const int c,
                                       const IWString & changing_atom_tag,
                                       IWString_and_File_Descriptor & output)
{
  Sparse_Fingerprint_Creator sfp;
  sfp.hit_bit(78, c);     // some arbitrary number

  IWString tmp;

  sfp.daylight_ascii_form_with_counts_encoded(tmp);

  output << changing_atom_tag << tmp << ">\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
write_reagent_count(RXN_File & rxn,
                    IWString_and_File_Descriptor & output)
{
  const int nreagents = rxn.number_reagents();

  IWString s(rxn.name());
  translate_strange_characters(s);
  output << s << ' ' << nreagents;

  if (nreagents > 1)
  {
    for (int i = 0; i < nreagents; ++i)
    {
      const auto & ri = rxn.reagent(i);
      output << ' ' << ri.natoms();
    }
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
write_changing_atom_count(RXN_File & rxn,
                          const int c,
                          IWString_and_File_Descriptor & output)
{
  IWString s;
  generate_externally_visible_reaction_name(rxn, s, false);

  output << s;
  output << ',' << c << '\n';   // Note csv, in case reaction name is multiple tokens.
   
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
write_reaction_smiles(RXN_File & rxn,
                      IWString_and_File_Descriptor & output)
{
  Reaction_Smiles_Options opts;

  if (! write_agent)
    opts.set_write_agent(0);

  rxn.write_rxn_smiles(opts, output);

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
write_reaction_and_maybe_changing_atom_count(RXN_File & rxn,
                                             ISIS_RXN_FILE_Molecule & r,
                                             const int c,
                                             IWString_and_File_Descriptor & output)
{
  if (columnar_output_specifications.active())
  {
    Reaction_Smiles_Options opts;
    opts.set_write_reaction_name(0);
    if (! write_agent)
      opts.set_write_agent(0);

    rxn.write_rxn_smiles(opts, output);

    if (append_changing_atoms_count_to_name)
      output << columnar_output_specifications.output_separator() << c;

    IWString s;
    generate_externally_visible_reaction_name(rxn, s, false);
    output << columnar_output_specifications.output_separator() << s;

    if (changing_atom_count_tag.length())
      output << columnar_output_specifications.output_separator() << c;
  }
  else                // tdt
  {
    output << smiles_tag << r.smiles() << ">\n";

    output << identifier_tag ;
    IWString s;
    generate_externally_visible_reaction_name(rxn, s, false);
    output << s;

    if (append_changing_atoms_count_to_name)
      output << ' ' << c;

    output << ">\n";

    if (changing_atom_count_tag.length())
      append_changing_atom_count_fingerprint(c, changing_atom_count_tag, output);
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

static int
ensure_locator_arrays_filled(RXN_File & rxn,
                             int * tmp)
{
  Molecule_to_Query_Specifications mqs;
  mqs.set_make_embedding(1);
  mqs.set_interpret_atom_alias_as_smarts(0);

  IWReaction notused;
  RXN_File_Create_Reaction_Options rxcfcro;
  rxcfcro.set_only_create_query_from_first_reagent(only_create_query_from_first_reagent);

  if (! rxn.create_reaction(notused, rxcfcro,  mqs, tmp))
  {
    cerr << "Cannot create reaction\n";
    return 0;
  }

  return 1;
}

/*
  Generate reaction fingerprints.
*/

static int
rxn_fingerprint(RXN_File & rxn,
                IWString_and_File_Descriptor & output)
{
  const int initial_nr = rxn.number_reagents();

  if (0 == initial_nr)
  {
    cerr << "Skipping reaction with no reagents " << rxn.name() << endl;
    if (stream_for_reagent_count.is_open())
      write_changing_atom_count(rxn, 0, stream_for_reagent_count);
    return 1;
  }

#define STUFF_NOW_DONE_BY_RXN2RXNSMI
#ifdef STUFF_NOW_DONE_BY_RXN2RXNSMI
  if (rxn.remove_fragments_not_participating())
    cerr << "Removed fragments not reacting from '" << rxn.name() << ", now " << rxn.number_reagents() << " reagents\n";

  if (rxn.eliminate_reagents_not_participating())
    cerr << "Removed reagents not reacting from '" << rxn.name() << "', now " << rxn.number_reagents() << " reagents\n";

  if (rxn.all_reagents_the_same())
    cerr << "Removed duplicate reagents from '" << rxn.name() << "' now " << rxn.number_reagents() << " reagents\n";

  if (remove_duplicate_reagents_atom_maps_scrambled && rxn.remove_duplicate_reagents_atom_maps_scrambled())
    cerr << "Removed duplicate reagents, but with different atom maps from '" << rxn.name() << "' now " << rxn.number_reagents() << " reagents\n";
//rxn.print_atom_map_into(cerr);

//const int nreagents = rxn.number_reagents();

//if (initial_nr != nreagents)
//  cerr << "Reagent count changed " << initial_nr << " to " << nreagents << " in " << rxn.name() << endl;
#endif

  const int nreagents = rxn.number_reagents();

  if (stream_for_reagent_count.is_open())
    write_reagent_count(rxn, stream_for_reagent_count);

  if (1 == nreagents)
    ;
  else if (take_first_of_multiple_reagents)
    took_first_of_multiple_reagents++;
  else if (skip_reactions_with_multiple_ragents)
  {
    if (verbose)
      cerr << "Skipping multi reagent " << nreagents << " reaction '" << rxn.name() << "'\n";

    if (stream_for_multiple_reagents.is_open())
    {
      stream_for_multiple_reagents << rxn.name() << '\n';
      stream_for_multiple_reagents.write_if_buffer_holds_more_than(4096);
    }

    reagents_with_multiple_reagents_skipped++;

    return 1;
  }
  else
  {
    cerr << "multiple reagents " << nreagents << " in '" << rxn.name() << "', cannot process\n";
    return 0;
  }

//if (any_changing_bond_means_a_changing_atoms)
//  rxn.set_is_changing_if_different_neighbours(1);

  auto & r = rxn.reagent(0);

  preprocess(r);

  const int matoms = r.natoms();

  int xx = rxn.highest_atom_map_number();
  if (xx < matoms)
    xx = matoms;

  RXNFP_Temporary_Arrays rxnfpta(matoms, xx + 1);

  int * changed = rxnfpta.changed();
  int * tmp = rxnfpta.tmp();

// We need to create a reaction in order for the reagent and product locator arrays to be filled

  if (consider_aromatic_bonds)
  {
    rxn.set_preserve_kekule_forms(0);
    rxn.set_aromatic_bonds_lose_kekule_identity(1);
  }
   
  ensure_locator_arrays_filled(rxn, tmp);

  if (discard_existing_atom_map && discard_reaction_unless_all_atoms_mapped)
  {
    const int mapped = r.number_mapped_atoms();
    if (mapped != r.natoms())
    {
      reactions_discarded_for_incomplete_atom_map++;
      if (verbose)
        cerr << rxn.name() << " discarded for unmapped atoms. Have " << r.natoms() << " atoms, mapped " << mapped << endl;
      unmapped_atom_count[matoms - mapped]++;
      acc_unmapped_atom_fraction.extra(static_cast<double>(matoms - mapped) / static_cast<double>(matoms));
      if (write_not_full_mapped_reactins_to_reaction_stream)
        rxn.do_write(stream_for_reactions);

      return 1;
    }
  }

  // Check for duplicate atom maps before anything else.
  // If left to later, this reaction caused problems.
  // CC1N=CC2C(C=1)=C([N+]([O-])=O)C=CC=2.[Cl:15][C:16]1[C:25]2[C:20](=[CH:21][CH:22]=[CH:23][CH:24]=2)[CH:19]=[CH:18][N:17]=1>>[ClH:15].[Cl:15][C:16]1[C:25]2[C:20](=[CH:21][CH:22]=[CH:23][CH:24]=2)[CH:19]=[CH:18][N:17]=1 |f:2.3|	US03930837		1976		

  if (rxn.contains_duplicate_atom_map_numbers())
  {
    reactions_containing_duplicate_atom_map_numbers++;
    return 1;
  }

  Changing_Atom_Conditions cac;
  if (any_changing_bond_means_a_changing_atoms)
    cac.set_is_changing_if_different_neighbours(1);
  if (only_consider_largest_reagent_fragment)
    cac.set_only_consider_largest_reagent_fragment(1);
  cac.set_ignore_lost_atom_if_isolated(1);
  cac.set_include_changing_bonds_in_changing_atom_count(include_changing_bonds_in_changing_atom_count);
  if (discern_changing_atoms_only_in_first_fragment)
    cac.set_discern_changing_atoms_only_in_first_fragment(1);
  cac.set_consider_aromatic_bonds(consider_aromatic_bonds);

//rxn.do_write("foo.rxn");

  int c;
  if (atom_typing_specification.active())
    c = rxn.identify_atoms_changing_reagent(0, atom_typing_specification, changed, cac);
  else
    c = rxn.identify_atoms_changing_reagent(0, changed, cac);

  if (verbose > 1)
    cerr << rxn.name() << " changing atom count " << c << " with " << rxn.number_reagents() << " reagents\n";

  if (changing_atoms_should_be.number_elements() > 0 &&
      ! changing_atoms_should_be.contains(c))
  {
    reactions_with_unexpected_changing_atom_counts++;

    rxn.discard_atom_map();
    IWReaction notused2;
    std::fill_n(tmp, matoms, 1);
    RXN_File_Create_Reaction_Options rxnfcro;
    rxnfcro.set_only_create_query_from_first_reagent(only_create_query_from_first_reagent);
    rxn.create_reaction(notused2, rxnfcro, tmp);
    if (atom_typing_specification.active())
      c = rxn.identify_atoms_changing_reagent(0, atom_typing_specification, changed, cac);
    else
      c = rxn.identify_atoms_changing_reagent(0, changed, cac);

    if (changing_atoms_should_be.contains(c))
    {
      reactions_with_corrected_atom_counts++;
      if (verbose > 1)
        cerr << rxn.name() << " corrected changing atom count\n";
    }
    else 
    {
      if (verbose > 1)
        cerr << rxn.name() << " cannot correct changing atom count, now " << c << endl;
      if (stream_for_reactions.is_open() && write_failed_changing_atom_counts_to_wfile)
        rxn.do_write(stream_for_reactions);
      return 1;    // ignore
    }
  }

  if (! rxn.at_least_some_mapped_atoms_common_btw_reagents_and_products())
  {
    reactions_with_no_reagent_atoms_in_products++;
    return 1;
  }

  if (rxn.reagent(0).number_fragments() > 1)
  {
    cerr << "Caution, first reagent has multiple fragments\n";
//  return 1;
  }

  if (stream_for_changing_atoms.is_open())
    write_changing_atom_count(rxn, c, stream_for_changing_atoms);

  acc_changing_atoms[c]++;

  if (0 == c && discard_reactions_with_no_changing_atoms)
  {
    reactions_with_no_changing_atoms_discarded++;
    if (verbose > 2)
      cerr << "No changing atoms\n";
    return 1;
  }

  if (changing_atoms_must_match.number_elements() && ! changing_atoms_match(r, c, changed, changing_atoms_must_match))
  {
    reactions_rejected_for_changing_atoms_must_match++;
    if (stream_for_reactions.is_open() && write_failed_changing_atom_counts_to_wfile)
      rxn.do_write(stream_for_reactions);
    if (verbose)
      cerr << "Reaction " << r.name() << " rejected for no match to changing atoms query, " << c << " changing atoms\n";
    return 1;
  }

  if (remove_non_participating_fragments)
    rxn.remove_non_participating_fragments();

  if (! atom_typing_specification.assign_atom_types(r, rxnfpta.atype()))
  {
    cerr << "Cannot assign atom types '" << rxn.name() << "'\n";
    return 0;
  }

// Now we can do what we came here for

  if (stream_for_reactions.is_open())
    rxn.do_write(stream_for_reactions);

  if (isotope > 0)
    do_apply_isotopes(r, changed, isotope);

  if (stream_for_reaction_smiles.is_open())
    write_reaction_smiles(rxn, stream_for_reaction_smiles);

  if (c < 100)    // Fixed size array, horrible.
    add_changing_atom_smarts_to_hash(r, changed, changing_atom_smarts[c]);

  write_reaction_and_maybe_changing_atom_count(rxn, r, c, output);

  for (int i = 0; i < nr; ++i)
  {
    if (! produce_fingerprint(r, radius[i], reagent_tag[i], rxnfpta, output))
      return 0;
  }

  if (fingerprint_changing_product_atoms)
    do_fingerprint_changing_product_atoms(rxn, output);

  if (difference_fingerprint_tag.length())
    do_generate_difference_fingerprint(rxn, rxnfpta, difference_fingerprint_tag, output);

  if (columnar_output_specifications.active())
    output << '\n';
  else
    output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return output.good();
}

static int
echo_bad_data(iwstring_data_source & input,
              const off_t initial_offset,
              IWString_and_File_Descriptor & output)
{
  const auto current_offset = input.tellg();
  input.seekg(initial_offset);
//cerr << "Bad reaction begin " << initial_offset << " now " << current_offset << endl;
  input.echo(output, (current_offset - initial_offset));
  input.seekg(current_offset, 0);
  assert (input.tellg() == current_offset);
  output.write_if_buffer_holds_more_than(4096);

  return input.tellg() == current_offset;
}

static int
next_reaction_smiles(iwstring_data_source & input,
                     RXN_File & rxn)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
    return 0;

//cerr << "SMILES INPUT '" << buffer << "'\n";

  if (! rxn.build_from_reaction_smiles(buffer, 1))
  {
    cerr << "Cannot interpret " << buffer << endl;
    return 0;
  }

  return 1;
}

static int
get_next_reaction(iwstring_data_source & input,
                  const char * fname,
                  RXN_File & rxn)
{
  if (input_is_reaction_smiles)
    return next_reaction_smiles(input, rxn);

  const auto initial_offset = input.tellg();

  if (rxn.do_read(input))
    ;
  else if (input.eof())
    return 0;
  else if (ignore_bad_reactions)
  {
    bad_reactions_ignored++;

    if (stream_for_bad_reactions.is_open())
      echo_bad_data(input, initial_offset, stream_for_bad_reactions);
  }
  else
  {
    cerr << "Fatal error reading reaction, now at line " << input.lines_read() << endl;
    return 0;
  }

  rxn.set_fname(fname);

  return 1;
}

static int
rxn_fingerprint(iwstring_data_source & input,
                const char * fname,
                IWString_and_File_Descriptor & output)
{
  input.set_translate_tabs(1);

  while (1)
  {
    RXN_File rxn;
    rxn.set_do_automatic_atom_mapping(0);

    const auto initial_offset = input.tellg();

    if (! get_next_reaction(input, fname, rxn))
    {
      if (input.eof())
        return 1;

      return 0;
    }

    if (discard_reactions_with_isotopes && rxn.contains_isotopic_reagent_atoms())
    {
      if (verbose)
        cerr << rxn.name() << " contains isotopic atoms\n";
      reactions_discarded_for_isotopic_atoms++;
      continue;
    }

    if (remove_cis_trans_bonding)
      rxn.remove_cis_trans_bonding();

    if (discard_existing_atom_map)
      rxn.discard_atom_map();

    reactions_read++;

    if (verbose > 1)
      cerr << "Processing '" << rxn.name() << "'\n";

    if (! rxn_fingerprint(rxn, output))
    {
      cerr << "rxn_fingerprint:fatal error processing '" << rxn.name() << "' now at line " << input.lines_read() << endl;

      if (! ignore_bad_reactions)
        return 0;

      if (stream_for_bad_reactions.is_open())
        echo_bad_data(input, initial_offset, stream_for_bad_reactions);

      bad_reactions_ignored++;
    }

//  if (reactions_read >= 3)
//    return 1;
  }

  return 1;
}

static int rxn_fingerprint(const char * fname, IWString_and_File_Descriptor & output);

static int
rxn_fingerprint_file_containing_reaction_files(iwstring_data_source & input,
                                         IWString_and_File_Descriptor & output)
{
  IWString buffer;

  while (input.next_record(buffer))
  {
    if (! rxn_fingerprint(buffer.null_terminated_chars(), output))
      return 0;
  }

  return 1;
}

static int
rxn_fingerprint_file_containing_reaction_files(const char * fname,
                                               IWString_and_File_Descriptor & output)
{
  IWString tmp(fname);
  assert (tmp.starts_with("F:"));
  tmp.remove_leading_chars(2);

  iwstring_data_source input(tmp.null_terminated_chars());

  if (! input.good())
  {
    cerr << "Cannot open file of reactions '" << tmp << "'\n";
    return 0;
  }

  return rxn_fingerprint_file_containing_reaction_files(input, output);
}

static int
rxn_fingerprint(const char * fname,
                IWString_and_File_Descriptor & output)
{
  assert(nullptr != fname);

  const_IWSubstring tmp(fname);
  if (tmp.starts_with("F:"))
    return rxn_fingerprint_file_containing_reaction_files(fname, output);

  iwstring_data_source input(fname);
  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  return rxn_fingerprint(input, fname, output);
}

void
WriteSortedByValueHashMap(const IW_STL_Hash_Map_int& hash,
                          std::ostream& output)
{
  const auto n = hash.size();
  int * tmp = new int[n + n]; std::unique_ptr<int[]> free_tmp(tmp);
  const_IWSubstring * keys = new const_IWSubstring[n]; std::unique_ptr<const_IWSubstring[]> free_keys(keys);

  // Load smarts and counts into arrays.
  int ndx = 0;
  for (const auto& kv : hash)
  {
    keys[ndx] = kv.first;
    tmp[ndx] = kv.second;
    ndx++;
  }

  int * permutations = tmp + n;
  std::iota(permutations, permutations + n, 0);

  std::sort(permutations, permutations + n, [&tmp](const int x1, const int x2) {
    return tmp[x1] > tmp[x2];
  });

  for (size_t i = 0; i < n; ++i)
  {
    output << ' ' << tmp[permutations[i]] << " occurrences of " << keys[permutations[i]] << "\n";
  }

  return;
}

static int
get_radii(const Command_Line & cl,
          const char flag,
          resizable_array<int> & radii)
{
  const_IWSubstring r;
  for (int i = 0; cl.value(flag, r, i); ++i)
  {
    int j = 0;
    const_IWSubstring token;
    while (r.nextword(token, j, ','))
    {
      int x;
      if (! token.numeric_value(x) || x < 0)
      {
        cerr << "Invalid radius '" << token << "'\n";
        return 0;
      }

//    cerr << "interpreted as " << x << endl;

      radii.add_if_not_already_present(x);
    }
  }

  radii.iwqsort_lambda([](const int r1, const int r2) {    // not really necessary to sort them
                           if (r1 < r2)
                             return -1;
                           else if (r1 > r2)
                             return 1;
                           else
                             return 0;         // will never happen
                            });

  return radii.number_elements();
}

static void
display_misc_options(std::ostream & output)
{
  output << " -M naqm         suppress 'no atoms in query' message\n";
  output << " -M CTAG=<tag>   add a fingerprint that contains the changing atom count\n";
  output << " -M noiso        discard reactions containing isotopic atoms\n";
  output << " -M cacff        discern changing atom count on the first fragment only\n";
  output << " -M cbc          changes in bonding will be included with the changing atom count\n";
  output << " -M skip0c       skip reactions with no changing atoms\n";
  output << " -M rmnp         remove non participating fragments from reactions\n";
  output << " -M rmagent      when writing reaction smiles, do NOT write the agent(s)\n";
  output << " -M ebch         generate extra bits for the changing atoms\n";
  output << " -M xctb         remove cis trans bonding from input\n";
  output << " -M anyf         changing atoms can be in any reagent\n";
  output << " -M aromb        when looking at changing bonds, do not consider aromatic Kekule forms\n";

  exit(1);
}

static int
rxn_fingerprint(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vA:E:g:m:r:J:P:I:xduC:W:ahc:ws:q:j:fbD:U:ByeX:F:M:RiK:V:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else 
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }
  else
    set_auto_create_new_elements(1);

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will remove small fragments in reagent molecules\n";
  }

  if (cl.option_present('P'))
  {
    const_IWSubstring p = cl.string_value('P');

    if (! atom_typing_specification.build(p))
    {
      cerr << "INvalid atom typing specification '" << p << "'\n";
      return 1;
    }
  }
  else    // too many difficulties if the atom typing ends up unset
  {
    atom_typing_specification.build("UST:ACHSY");
  }

  if (cl.option_present('M'))
  {
    const_IWSubstring m;
    for (int i = 0; cl.value('M', m, i); ++i)
    {
      if ("naqm" == m)
      {
        set_iwreaction_display_no_atoms_in_query_message(0);
      }
      else if (m.starts_with("CTAG="))
      {
        m.remove_leading_chars(5);
        changing_atom_count_tag = m;
        if (verbose)
          cerr << "changing atom count written to fingerprint with tag '" << changing_atom_count_tag << "'\n";
        if (! changing_atom_count_tag.ends_with('<'))
          changing_atom_count_tag << '<';
      }
      else if ("noiso" == m)
      {
        discard_reactions_with_isotopes = 1;
        if (verbose)
          cerr << "Will discard reactions containing isotopic atoms\n";
      }
      else if ("cacff" == m)
      {
        discern_changing_atoms_only_in_first_fragment = 1;
        if (verbose)
          cerr << "Will only discern changing atom counts in the first reagent fragment\n";
      }
      else if ("cbc" == m)
      {
        include_changing_bonds_in_changing_atom_count = 1;
        if (verbose)
          cerr << "Will discern changing bonds when computing changing atom count\n";
      }
      else if ("skip0c" == m)
      {
        discard_reactions_with_no_changing_atoms = 1;
        if (verbose)
          cerr << "Will skip reactions with no changing atoms\n";
      }
      else if ("rmnp" == m)
      {
        remove_non_participating_fragments = 1;
        if (verbose)
          cerr << "Will remove non participating fragments\n";
      }
      else if ("rmagent" == m)
      {
        write_agent = 0;
        if (verbose)
          cerr << "Will not write agents\n";
      }
      else if ("ebch" == m)
      {
        set_extra_bits_for_changing_atoms = 1;
        if (verbose)
          cerr << "Will set extra bits for the changing atoms\n";
      }
      else if ("xctb" == m)
      {
        remove_cis_trans_bonding = 1;

        if (verbose)
          cerr << "Will remove cis-trans bonding information\n";
      }
      else if ("anyf" == m)
      {
        only_create_query_from_first_reagent = 0;
        if (verbose)
          cerr << "Changing atoms can be in any reagent\n";
      }
      else if ("aromb" == m)
      {
        consider_aromatic_bonds = 1;
        if (verbose)
          cerr << "Will look at bond aromaticity in determining changing bonds\n";
      }
      else if ("help" == m)
      {
        display_misc_options(cerr);
      }
      else
      {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_misc_options(cerr);
      }
    }
  }

  if (cl.option_present('y'))
  {
    lfpd.set_place_bits_for_rings(0);

    if (verbose)
      cerr <<"Will exclude ring information from fingerprints\n";
  }

  lfpd.set_differentiate_ring_and_chain_bonds(1);

  if (cl.option_present('e'))
  {
    only_consider_largest_reagent_fragment = 1;
    if (verbose)
      cerr << "Will ignore small reagent fragments when determining chaning atom count\n";
  }

  if (cl.option_present('i'))
  {
    input_is_reaction_smiles = 1;

    if (verbose)
      cerr << "Input is reaction smiles\n";
  }

  if (cl.option_present('m'))
  {
    const_IWSubstring m;

    for (int i = 0; cl.value('m', m, i); ++i)    // does not need to be a loop, these are mutually exclusive...
    {
      if ("first" == m)
      {
        take_first_of_multiple_reagents = 1;
        if (verbose)
          cerr << "Will process the first of multiple reagents\n";
      }
      else if ("skip" == m)
      {
        skip_reactions_with_multiple_ragents = 1;
        if (verbose)
          cerr << "Will skip reactions with multiple reagents\n";
      }
      else
      {
        cerr << "Unrecognised -m qualifier '" << m << "'\n";
        usage(1);
      }
    }
  }

  if (cl.option_present('s'))
  {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i)
    {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s))
      {
        cerr << "Invalid matched atoms should be smarts '" << s << "'\n";
        return 1;
      }

      changing_atoms_must_match.add(q.release());
    }
  }

  if (cl.option_present('q'))
  {
    if (! process_queries(cl, changing_atoms_must_match, 'q', verbose))
    {
      cerr << "Cannot process changing atoms must match specification (-q)\n";
      return 1;
    }
  }

  if (cl.option_present('a'))
  {
    append_changing_atoms_count_to_name = 1;

    if (verbose)
      cerr << "Will append changing atom count to name\n";
  }

  if (cl.option_present('h'))
  {
    any_changing_bond_means_a_changing_atoms = 1;

    if (verbose)
      cerr << "If any different mapped atoms are attached, then an atom is marked changing\n";
  }

  if (cl.option_present('c'))
  {
    int c;
    for (int i = 0; cl.value('c', c, i); ++i)
    {
      changing_atoms_should_be.add(c);
      if (verbose)
        cerr << "If the number of changing atoms is not " << c << " will recompute the atom map\n";
    }
  }

  if (cl.option_present('x'))
  {
    discard_existing_atom_map = 1;

    if (verbose)
      cerr << "Will discard existing atom map data\n";

    if (cl.option_present('j'))
    {
      discard_reaction_unless_all_atoms_mapped = 1;
      if (verbose)
        cerr << "Unless all reagent atoms are mapped, reaction will be discarded\n";
      const_IWSubstring j = cl.string_value('j');
      if ("write" == j)
        write_not_full_mapped_reactins_to_reaction_stream = 1;
    }
  }

  if (cl.option_present('f'))
  {
    use_first_token_of_name = 1;

    if (verbose)
      cerr << "Will truncate reaction names to the first token\n";
  }

  if (cl.option_present('U'))
  {
    if (input_is_reaction_smiles)
    {
      cerr << "The -U option is incompatible with reaction smiles input\n";
      usage(1);
    }

    const_IWSubstring q = cl.string_value('U');

    if ('f' == q)
    {
      reaction_name_is_file_name = 1;
      if (verbose)
        cerr << "Will use the reaction file name as the reaction name\n";
    }
    else if ("nrxn" == q)
    {
      reaction_name_is_file_name = 1;
      strip_rxn_from_file_name = 1;
      if (verbose)
        cerr << "Will use the reaction file name (without suffix) as the reaction name\n";
    }
    else if ('p' == q)
    {
      reaction_name_is_path_name = 1;
      if (verbose)
        cerr << "Will use the reaction path name as the reaction name\n";
    }
    else
    {
      cerr << "Unrecognised reaction name modifier '" << q << "'\n";
      usage(1);
    }
  }

  if (cl.option_present('B'))
  {
    remove_duplicate_reagents_atom_maps_scrambled = 1;
    if (verbose)
      cerr << "Will remove duplicate reagents even if atom maps scrambled\n";
  }

  if (cl.option_present('b'))
  {
    convert_products_to_agents_when_possible = 1;
    if (verbose)
      cerr << "Where possible, will convert product molecules to agents\n";
  }

  if (cl.option_present('d'))
  {
    expand_across_double_bonds = 1;

    if (verbose)
      cerr << "Will include doubly bonded atoms in shell expansion\n";
  }

  if (cl.option_present('I'))
  {
    if (! cl.value('I', isotope) || isotope < 1)
    {
      cerr << "The isotope flag (-I) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will label changing atoms with isotope " << isotope << endl;
  }

  if (cl.option_present('u'))
  {
    fingerprint_changing_product_atoms = 1;

    if (verbose)
      cerr << "Will product fingerprints of changing product atoms\n";
  }

  resizable_array<int> tmpradii;

  nr = get_radii(cl, 'r', tmpradii);

  if (0 == nr)
  {
    cerr << "Must specify radius from changing atoms via the -r option\n";
    usage(1);
  }

  if (cl.option_present('K'))
  {
    if (cl.option_present('J'))
    {
      cerr << "When requesting columnar output (-K) cannot also specify fingerprint tag (-J)\n";
      usage(1);
    }

    if (! columnar_output_specifications.build(cl, 'K', verbose))
    {
      cerr << "Cannot initialise columnar output specifications (-K)\n";
      return 1;
    }
  }
  else if (cl.option_count('J') == nr)
    ;
  else if (1 == cl.option_count('J'))
    ;
  else
  {
    cerr << "Must have a fingerprint tag (-J) for each radius specified\n";
    usage(1);
  }

  radius = new_int(nr);std::unique_ptr<int[]> free_radius(radius);
  reagent_tag = new IWString[nr];std::unique_ptr<IWString[]> free_reagent_tag(reagent_tag);
  product_tag = new IWString[nr];std::unique_ptr<IWString[]> free_product_tag(product_tag);   // even if not used

  if (cl.option_present('J'))
  {
    IWString common_j;
    if (1 == cl.option_count('J'))
      cl.value('J', common_j);

    for (int i = 0; i < nr; ++i)
    {
      radius[i] = tmpradii[i];

      IWString j;
      if (common_j.length())
      {
        j = common_j;
        j << radius[i];
      }
      else
        cl.value('J', j, i);

      if (j.ends_with('<'))
        j.chop();         // yes, we might just add it back just below

      if (fingerprint_changing_product_atoms && ! j.starts_with("NC"))
      {
        cerr << "Sorry, when fingerprinting product atoms can only use sparse fingerprints (NC)\n";
        return 1;
      }

      if (j.starts_with("FP") || j.starts_with("NC"))
      {
        if (fingerprint_changing_product_atoms)
        {
          reagent_tag[i] = j;
          product_tag[i] = j;
          reagent_tag[i] << "R<";
          product_tag[i] << "P<";
        }
        else
        {
          reagent_tag[i] = j;
          reagent_tag[i] << '<';
        }
      }
      else
      {
        cerr << "Unrecognised tag form '" << j << "'\n";
        return 1;
      }
    }
  }
  else if (columnar_output_specifications.active())
  {
    for (int i = 0; i < nr; ++i)
    {
      radius[i] = tmpradii[i];
    }
  }

  if (verbose)
  {
    cerr << "radii";
    for (int i = 0; i < nr; ++i)
    {
      cerr << ' ' << radius[i];
    }
    cerr << endl;
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('C'))
  {
    const char * c = cl.option_value('C');

    if (! stream_for_changing_atoms.open(c))
    {
      cerr << "Cannot open stream for changing atom counts '" << c << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Counts of changing atoms written to '" << c << "'\n";

    stream_for_changing_atoms << "ID RGNT_Change\n";
  }

  if (cl.option_present('D'))
  {
    const char * d = cl.option_value('D');

    if (! stream_for_reagent_count.open(d))
    {
      cerr << "cannot open stream for reagent counts '" << d << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Reagent counts written to '" << d << "'\n";

    stream_for_reagent_count << "ID NReagents\n";
  }

  if (cl.option_present('W'))
  {
    const char * w = cl.option_value('W');

    stream_for_reactions.open(w, std::ios::out);
    if (! stream_for_reactions.good())
    {
      cerr << "Cannot open stream for reactions '" << w << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Reactions written to '" << w << "'\n";

    if (cl.option_present('w'))
    {
      write_failed_changing_atom_counts_to_wfile = 1;
      if (verbose)
        cerr << "Reactions that cannot be re-mapped will also be written to the -W file\n";
    }
  }

  if (cl.option_present('X'))
  {
    const char * x = cl.option_value('X');

    if (! stream_for_reaction_smiles.open(x))
    {
      cerr << "Cannot open stream for reaction smiles '" << x << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Reaction smiles written to '" << x << "'\n";
  }

  if (cl.option_present('F'))
  {
    IWString f = cl.string_value('F');

    ignore_bad_reactions = 1;

    if ('.' == f)
    {
      if (verbose)
        cerr << "Will ignore failed reactions\n";
    }
    else
    {
      if (! stream_for_bad_reactions.open(f.null_terminated_chars()))
      {
        cerr << "Cannot open stream for bad reactions '" << f << "'\n";
        return 1;
      }

      if (verbose)
        cerr << "Will write failed reactions to '" << f << "'\n";
    }
  }

  if (cl.option_present('V'))
  {
    cl.value('V', difference_fingerprint_tag);

    if (! difference_fingerprint_tag.starts_with("NC"))
    {
      cerr << "Difference fingerprint tags must be non colliding form '" << difference_fingerprint_tag << "' invalid\n";
      return 1;
    }

    if (! difference_fingerprint_tag.ends_with("<"))
      difference_fingerprint_tag << '<';

    if (verbose)
      cerr << "Will form reaction difference fingerprint with tag " << difference_fingerprint_tag << endl;
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! rxn_fingerprint(cl[i], output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  if (verbose)
  {
    cerr << "Read " << reactions_read << " reactions\n";
    Accumulator_Int<int> acc;

    int most_common_count = 0;
    for (int i = 0; i < acc_changing_atoms.number_elements(); ++i)
    {
      const int c = acc_changing_atoms[i];
      if (0 == c)
        continue;

      acc.extra(i, c);
      cerr << c << " reactions had " << i << " changing atoms\n";

      if (c > most_common_count)
        most_common_count = c;
    }

    if (acc.n() > 0)
    {
      cerr << "Changing atom counts btw " << acc.minval() << " and " << acc.maxval() << " ave " << static_cast<float>(acc.average()) << endl;
      cerr << "Purity " << static_cast<float>(most_common_count) / static_cast<float>(acc.n()) << endl;
    }

    if (changing_atoms_should_be.number_elements() > 0)
    {
      cerr << reactions_with_unexpected_changing_atom_counts << " reactions initially had changing atom counts different from";
      for (int i = 0; i < changing_atoms_should_be.number_elements(); ++i)
      {
        cerr << ' ' << changing_atoms_should_be[i];
      }
      cerr << endl;
      cerr << reactions_with_corrected_atom_counts << " of those reactions re-mapped to yield correct changing atom counts. Leaves " << (reactions_with_unexpected_changing_atom_counts - reactions_with_corrected_atom_counts) << " invalid, init " << reactions_read << endl;
    }

    if (discard_reaction_unless_all_atoms_mapped)
    {
      cerr << "discarded " << reactions_discarded_for_incomplete_atom_map << " reactions for incomplete atom mapping\n";
      for (int i = 0; i < unmapped_atom_count.number_elements(); ++i)
      {
        if (unmapped_atom_count[i])
          cerr << unmapped_atom_count[i] << " reactions had " << i << " unmapped atoms\n";
      }

      if (acc_unmapped_atom_fraction.n())
        cerr << "Unmapped atom fraction btw " << static_cast<float>(acc_unmapped_atom_fraction.minval()) << " and " << static_cast<float>(acc_unmapped_atom_fraction.maxval()) << " ave " << static_cast<float>(acc_unmapped_atom_fraction.average()) << endl;
    }

    if (changing_atoms_must_match.number_elements())
      cerr << reactions_rejected_for_changing_atoms_must_match << " reactions rejected for mismatch to changing atoms query\n";

    for (int i = 0; i < 100; ++i)
    {
      const IW_STL_Hash_Map_int & hi = changing_atom_smarts[i];

      if (0 == hi.size())
        continue;

      cerr << "When " << i << " changing atoms\n";
      WriteSortedByValueHashMap(hi, cerr);
    }

    cerr << reagents_with_multiple_reagents_skipped << " reactions skipped for multiple reagents\n";

    if (took_first_of_multiple_reagents)
      cerr << " took the first of multiple reagents " << took_first_of_multiple_reagents << " times\n";
  }

  if (bad_reactions_ignored)
    cerr << "SKIPPED " << bad_reactions_ignored << " bad reactions\n";

  if (reactions_discarded_for_isotopic_atoms)
    cerr << "SKIPPED " << reactions_discarded_for_isotopic_atoms << " reactions containing isotopic atoms\n";

  if (reactions_with_no_changing_atoms_discarded)
    cerr << "SKIPPED " << reactions_with_no_changing_atoms_discarded << " reactions with no changing atoms\n";

  if (reactions_with_no_reagent_atoms_in_products)
    cerr << "SKIPPED " << reactions_with_no_reagent_atoms_in_products << " reactions_with_no_reagent_atoms_in_products\n";

  if (reactions_containing_duplicate_atom_map_numbers)
    cerr << "SKIPPED " << reactions_containing_duplicate_atom_map_numbers << " reactions with duplicate atom map numbers\n";

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = rxn_fingerprint(argc, argv);

  return rc;
}
