// Descriptor computation.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>

#include <math.h>
#include <signal.h>
#include <sys/signal.h>
#include <unistd.h>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwmisc/set_or_unset.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/donor_acceptor.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/path_around_ring.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/qry_wcharge.h"
#include "Molecule_Tools/partial_symmetry.h"
#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/xlogp.h"

#include "Molecule_Tools/iwdescr.pb.h"

using std::cerr;
using iwmisc::Fraction;

const char * prog_name = nullptr;

static int reduce_to_largest_fragment = 0;

/*
  In order to deal with undefined values, we designate a special
  character for undefined values
*/

static IWString undefined_value = '.';

// Prepended to the name of each descriptor.
static IWString descriptor_prefix = "w_";

static Charge_Assigner charge_assigner;

static Donor_Acceptor_Assigner donor_acceptor_assigner;

static Chemical_Standardisation chemical_standardisation;

static int min_hbond_feature_separation = 0;
static int max_hbond_feature_separation = std::numeric_limits<int>::max();

// We can read a mapping from legacy names to externally specified names.
static IW_STL_Hash_Map_String name_translation;

static int flush_after_each_molecule = 0;

static alogp::ALogP alogp_engine;

// Normally we read a smiles file and produce a descriptor file.
// Optionally we can read a descriptor file with a smiles as the first
// column, and optionally, write the same thing.
static int read_descriptor_file_pipeline = 0;
static int write_descriptor_file_pipeline = 0;

/*
  Since the strongly fused ring system descriptors are sparse, they are
  optional too
*/

// Various features can be computed optionally. Put them
// all into a struct that can be initialised with the -O option.
// Some of these have default values of 1, some 0
struct DescriptorsToCompute {
  int adjacent_ring_fusion_descriptors = 1;
  int bonds_between_rings = 1;
  int charge_descriptors = 1;
  int complexity_descriptors = 1;
  int crowding_descriptors = 1;
  int distance_matrix_descriptors = 1;
  int donor_acceptor = 1;
  int hbond_descriptors = 1;
  int simple_hbond_descriptors = 1;
  int ncon_descriptors = 1;
  int perform_expensive_chirality_perception = 1;
  int partial_symmetry_descriptors = 1;
  int polar_bond_descriptors = 1;
  int psa = 1;
  int ramey_descriptors = 1;
  int ring_chain_descriptors = 1;
  int ring_fusion_descriptors = 1;
  int ring_substitution_descriptors = 1;
  int compute_alogp = 1;
  int compute_xlogp = 1;
  int ring_substitution_ratio_descriptors = 1;
  int specific_groups = 1;
  int spinach_descriptors = 1;
  int symmetry_descriptors = 1;
#ifdef MCGOWAN
  int mcgowan = 1;
#endif

  // As new descriptors are added, make sure you add them to the "all" and "none"
  // directives.

  // Set all known optional descriptor types to `s`.
  // As new descriptors are added, make sure you add them to this function.
  void SetAll(int s);

  // Display usage message and return 0.
  int DisplayUsage(char flag) const;

  int ReadDescriptorsToCompute(const const_IWSubstring& fname);
  int ReadDescriptorsToCompute(iwstring_data_source& input);

  public:
    int Initialise(Command_Line& cl);
};

void
DescriptorsToCompute::SetAll(int s) {
  adjacent_ring_fusion_descriptors = s;
  bonds_between_rings = s;
  charge_descriptors = s;
  complexity_descriptors = s;
  crowding_descriptors = s;
  distance_matrix_descriptors = s;
  donor_acceptor = s;
  hbond_descriptors = s;
  ncon_descriptors = s;
  perform_expensive_chirality_perception = s;
  partial_symmetry_descriptors = s;
  polar_bond_descriptors = s;
  ring_chain_descriptors = s;
  ring_substitution_descriptors = s;
  compute_alogp = s;
  compute_xlogp = s;
  ring_substitution_ratio_descriptors = s;
  simple_hbond_descriptors = s;
  specific_groups = s;
  spinach_descriptors = s;
  symmetry_descriptors = s;
  psa = s;
#ifdef MCGOWAN
  mcgowan = s;
#endif  // MCGOWAN
  ring_fusion_descriptors = s;
  ramey_descriptors = s;
}

int
DescriptorsToCompute::DisplayUsage(char flag) const {
  // clang-format off
  IWString dash_o;
  dash_o << " -" << flag << ' ';
  cerr << dash_o << "none          turn off all optional descriptors\n";
  cerr << dash_o << "all           turn on  all optional descriptors\n";
  cerr << dash_o << "adjring       enable adjacent to ring fusion descriptors\n";
  cerr << dash_o << "bbr           enable bonds between rings descriptors\n";
  cerr << dash_o << "charge        enable all formal charge descriptors\n";
  cerr << dash_o << "chiral        enable all expensive chirality perception descriptors\n";
  cerr << dash_o << "complex       enable planar fused ring descriptors\n";
  cerr << dash_o << "crowd         enable atomic crowding descriptors\n";
  cerr << dash_o << "dm            enable descriptors based on the distance matrix\n";
  cerr << dash_o << "donacc        enable donor acceptor derived descriptors\n";
  cerr << dash_o << "hbond         enable donor/acceptor derived descriptors\n";
  cerr << dash_o << "shbond        enable simplistic donor/acceptor derived descriptors\n";
  cerr << dash_o << "ncon          enable ncon and fncon descriptors\n";
  cerr << dash_o << "pbond         enable polar bond derived descriptors\n";
  cerr << dash_o << "psa           enable Novartis polar surface area descriptors\n";
#ifdef MCGOWAN
  cerr << dash_o << "mcgowan       enable McGowan volume descriptors\n";
#endif  // MCGOWAN
  cerr << dash_o << "psymm         enable partial symmetry derived descriptors\n";
  cerr << dash_o << "ramey         enable Ramey (element count) descriptors\n";
  cerr << dash_o << "rcj           enable ring chain join descriptors\n";
  cerr << dash_o << "rfuse         enable ring fusion descriptors\n";
  cerr << dash_o << "rss           enable ring substitution descriptors\n";
  cerr << dash_o << "rssr          enable ring substitution ratio descriptors\n";
  cerr << dash_o << "spch          enable spinach related descriptors\n";
  cerr << dash_o << "spcgrp        enable specific group descriptors\n";
  cerr << dash_o << "symm          enable symmetry related descriptors\n";
  cerr << dash_o << "xlogp         enable xlogp\n";
  // clang-format on

  return 0;
}

int
DescriptorsToCompute::Initialise(Command_Line& cl) {
  static constexpr char kFlag = 'O';
  const int verbose = cl.option_present('v');
  const_IWSubstring o;
  for (int i = 0; cl.value(kFlag, o, i); ++i) {
    if (o == "none") {
      SetAll(0);
      if (verbose) {
        cerr << "All optional descriptors turned off\n";
      }
    } else if (o == "all") {
      SetAll(1);
      if (verbose) {
        cerr << "All optional descriptors enabled\n";
      }
    } else if (o == "adjring") {
      adjacent_ring_fusion_descriptors = 1;
      if (verbose) {
        cerr << "bonding adjacent to ring fusion descriptors enabled\n";
      }
    } else if (o == "bbr") {
      bonds_between_rings = 1;
      if (verbose) {
        cerr << "bonds between rings descriptors enabled\n";
      }
    } else if (o == "charge") {
      charge_descriptors = 1;
      if (verbose) {
        cerr << "charge derived descriptors enabled\n";
      }
    } else if (o == "complex") {
      complexity_descriptors = 1;
      if (verbose) {
        cerr << "complexity derived descriptors enabled\n";
      }
    } else if (o == "crowd") {
      crowding_descriptors = 1;
      if (verbose) {
        cerr << "atomic crowding descriptors enabled\n";
      }
    } else if (o == "dm") {
      distance_matrix_descriptors = 1;
      if (verbose) {
        cerr << "distance matrix derived descriptors enabled\n";
      }
    } else if (o == "donacc") {
      donor_acceptor = 1;
      if (verbose) {
        cerr << "donor acceptor derived descriptors enabled\n";
      }
    } else if (o == "hbond") {
      hbond_descriptors = 1;
      if (verbose) {
        cerr << "hydrogen bonding derived descriptors enabled\n";
      }
    } else if (o == "shbond") {
      simple_hbond_descriptors = 1;
      if (verbose) {
        cerr << "simplistic hydrogen bonding derived descriptors enabled\n";
      }
    } else if (o == "chiral" || o == "expchiral") {
      perform_expensive_chirality_perception = 1;;
      if (verbose) {
        cerr << "Will check all unlabelled atoms for possible chirality\n";
      }
    } else if (o == "ncon") {
      ncon_descriptors = 1;;
      if (verbose) {
        cerr << "ncon and fncon descriptors enabled\n";
      }
    } else if (o == "pbond") {
      polar_bond_descriptors = 1;;
      if (verbose) {
        cerr << "polar bond derived descriptors enabled\n";
      }
    } else if (o == "psa") {
      psa = 1;;
      if (verbose) {
        cerr << "Novartis polar surface area descriptor enabled\n";
      }
#ifdef MCGOWAN
    } else if (o == "mcgowan") {
      mcgowan = 1;;
      if (verbose) {
        cerr << "McGowan volumes descriptor enabled\n";
      }
#endif
    } else if (o == "psymm") {
      partial_symmetry_descriptors = 1;;
      if (verbose) {
        cerr << "partial symmetry derived descriptors enabled\n";
      }
    } else if (o == "symm") {
      symmetry_descriptors = 1;;
      if (verbose) {
        cerr << "symmetry derived descriptors enabled\n";
      }
    } else if (o == "ramey") {
      ramey_descriptors = 1;;
      if (verbose) {
        cerr << "element count descriptors enabled\n";
      }
    } else if (o == "rcj") {
      ring_chain_descriptors = 1;;
      if (verbose) {
        cerr << "ring chain join descriptors enabled\n";
      }
    } else if (o == "rfuse") {
      ring_fusion_descriptors = 1;;
      if (verbose) {
        cerr << "ring fusion descriptors enabled\n";
      }
    } else if (o == "rss") {
      ring_substitution_descriptors = 1;;
      if (verbose) {
        cerr << "ring substitution descriptors enabled\n";
      }
    } else if (o == "rssr") {
      ring_substitution_ratio_descriptors = 1;;
      if (verbose) {
        cerr << "ring substitution ratio descriptors enabled\n";
      }
    } else if (o == "spch") {
      spinach_descriptors = 1;;
      if (verbose) {
        cerr << "spinach derived descriptors enabled\n";
      }
    } else if (o == "spcgrp") {
      specific_groups = 1;;
      if (verbose) {
        cerr << "specific substructure descriptors enabled\n";
      }
    } else if (o == "alogp") {
      compute_alogp = 1;
      if (verbose) {
        cerr << "alogp will be computed\n";
      }
    } else if (o == "xlogp") {
      compute_xlogp = 1;
      if (verbose) {
        cerr << "xlogp will be computed\n";
      }
    } else if (o.starts_with("F:")) {
      o.remove_leading_chars(2);
      if (! ReadDescriptorsToCompute(o)) {
        cerr << "DescriptorsToCompute::Initialise:cannot process '" << o << "'\n";
        return 0;
      }
    } else if (o == "help") {
      return DisplayUsage(kFlag);
    } else {
      cerr << "DescriptorsToCompute::Initialise:unrecognised '" << o << "'\n";
      return DisplayUsage(kFlag);
    }
  }

  return 1;
}

int
DescriptorsToCompute::ReadDescriptorsToCompute(const const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "DescriptorsToCompute::ReadDescriptorsToCompute:cannot open '" << fname << "'\n";
    return 0;
  }

  return DescriptorsToCompute::ReadDescriptorsToCompute(input);
}

int
DescriptorsToCompute::ReadDescriptorsToCompute(iwstring_data_source& input) {
  input.set_strip_trailing_blanks(1);

  // Turn off everything to start
  SetAll(0);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with("#")) {
      continue;
    }

    // This code is copied from above, so keep the same variable name.
    // Should set up an Activate() function...
    const const_IWSubstring o(buffer);

    if (o == "none") {
      SetAll(0);
    } else if (o == "all") {
      SetAll(1);
    } else if (o == "adjring") {
      adjacent_ring_fusion_descriptors = 1;
    } else if (o == "bbr") {
      bonds_between_rings = 1;
    } else if (o == "charge") {
      charge_descriptors = 1;
    } else if (o == "complex") {
      complexity_descriptors = 1;
    } else if (o == "crowd") {
      crowding_descriptors = 1;
    } else if (o == "dm") {
      distance_matrix_descriptors = 1;
    } else if (o == "donacc") {
      donor_acceptor = 1;
    } else if (o == "hbond") {
      hbond_descriptors = 1;
    } else if (o == "shbond") {
      simple_hbond_descriptors = 1;
    } else if (o == "chiral" || o == "expchiral") {
      perform_expensive_chirality_perception = 1;;
    } else if (o == "ncon") {
      ncon_descriptors = 1;;
    } else if (o == "pbond") {
      polar_bond_descriptors = 1;;
    } else if (o == "psa") {
      psa = 1;;
    } else if (o == "psymm") {
      partial_symmetry_descriptors = 1;;
    } else if (o == "symm") {
      symmetry_descriptors = 1;;
    } else if (o == "ramey") {
      ramey_descriptors = 1;;
    } else if (o == "rcj") {
      ring_chain_descriptors = 1;;
    } else if (o == "rfuse") {
      ring_fusion_descriptors = 1;;
    } else if (o == "rss") {
      ring_substitution_descriptors = 1;;
    } else if (o == "rssr") {
      ring_substitution_ratio_descriptors = 1;;
    } else if (o == "spch") {
      spinach_descriptors = 1;;
    } else if (o == "spcgrp") {
      specific_groups = 1;;
    } else if (o == "alogp") {
      compute_alogp = 1;
    } else if (o == "xlogp") {
      compute_xlogp = 1;
    } else {
      cerr << "DescriptorsToCompute::ReadDescriptorsToCompute:unrecognised '" << o << "'\n";
      return 0;
    }
  }

  return 1;
}

static DescriptorsToCompute descriptors_to_compute;

/*
  Strongly fused is intended to apply to molecules that will have strained systems.
  But if the pair of rings that is share > 1 bond are of very different sizes, then
  that doesn't apply
*/

static int max_difference_in_ring_size_for_strongly_fused = 0;

/*
  partial charge descriptors are based on substructure queries - these
  aren't used
*/

static resizable_array_p<Query_and_Charge_Stats> charge_queries;

template class resizable_array_p<Query_and_Charge_Stats>;

/*
  We can test random smiles permutations
*/

static int ntest = 0;

static int test_failures = 0;

static int keep_going_after_test_failure = 0;

static int output_precision = 4;

static int compute_spiro_fusions = 1;

static int verbose = 0;
static uint64_t molecules_read = 0;
static int molecules_with_no_rings = 0;

static unsigned int alarm_time = 0;

static int molecules_skipped_by_timer = 0;

static IWString tag;           // if producing fingerprints
static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int work_as_tdt_filter = 0;

/*
  Apr 2015 can be handy to have the smiles in the output
*/

static int include_smiles_as_descriptor = 0;

/*
  We have an optional behaviour for the fsdrng* descriptors.
  Do they just consider ring system sizes of two, or any context
*/

static int fsdrng_descriptors_consider_just_two_ring_systems = 1;

static Molecular_Weight_Control mwc;

static int ignore_molecules_with_no_atoms = 0;

// Many of the features are int forms.
static IWDigits iwdigits;

enum IWDescr_Enum
{
  iwdescr_natoms,
  iwdescr_nrings,
  iwdescr_nelem,
  iwdescr_amw,
  iwdescr_ncon1,
  iwdescr_fncon1,
  iwdescr_ncon2,
  iwdescr_fncon2,
  iwdescr_ncon3,
  iwdescr_fncon3,
  iwdescr_ncon4,
  iwdescr_fncon4,
  iwdescr_frhc,
  iwdescr_mltbd,
  iwdescr_fmltbd,
  iwdescr_chmltbd,
  iwdescr_fchmltbd,
  iwdescr_rgmltbd,
  iwdescr_frgmltbd,
  iwdescr_dcca,
  iwdescr_fdcca,
  iwdescr_mxdst,
  iwdescr_fmxdst,
  iwdescr_mxsdlp,
  iwdescr_avsdlp,
  iwdescr_mxsdlprl,
  iwdescr_mdallp,
  iwdescr_fmdallp,
  iwdescr_fdiffallp,
  iwdescr_harary,
  iwdescr_rotbond,
  iwdescr_frotbond,
  iwdescr_ringatom,
  iwdescr_rhacnt,
  iwdescr_rhaf,
  iwdescr_frafus,
  iwdescr_rngatmf,
  iwdescr_aroma,
  iwdescr_aromha,
  iwdescr_fraromha,
  iwdescr_aromdens,
  iwdescr_ch2,
  iwdescr_ch,
  iwdescr_htroatom,
  iwdescr_htroaf,
  iwdescr_nrgnhlht,
  iwdescr_ohsh,
  iwdescr_co2h,
  iwdescr_amine,
  iwdescr_pyridine,
  iwdescr_pyrrole,
  iwdescr_hacts,
  iwdescr_hdons,
  iwdescr_hduals,
  iwdescr_mhr,
  iwdescr_mxhrf,
  iwdescr_mnhrf,
  iwdescr_lrsysz,
  iwdescr_srsz,
  iwdescr_lrsz,
  iwdescr_rng7atoms,
  iwdescr_nrsyscmr,
  iwdescr_mars,
  iwdescr_frspch,
  iwdescr_spchtro,
  iwdescr_rbfrspch,
  iwdescr_satspcha,
  iwdescr_unsatspcha,
  iwdescr_fsatspcha,
  iwdescr_scaffoldbranches,
  iwdescr_nrnspch,
  iwdescr_fnrnspc,
  iwdescr_trmnlrng,
  iwdescr_intrnlrng,
  iwdescr_rng2spch,
  iwdescr_rng2bridge,
  iwdescr_rcj,
  iwdescr_rchj,
  iwdescr_amrcj,
  iwdescr_alrcj,
  iwdescr_pbcount,
  iwdescr_frpbond,
  iwdescr_nonpbond,
  iwdescr_pbarom,
  iwdescr_npbarom,
  iwdescr_pbunset,
  iwdescr_dvinylb,
  iwdescr_ringsys,
  iwdescr_arring,
  iwdescr_alring,
  iwdescr_excybond,
  iwdescr_excydbond,
  iwdescr_excydscon,
  iwdescr_excydsconh,
  iwdescr_excydscondon,
//iwdescr_scra,
//iwdescr_scrha,
//iwdescr_scrd,
  iwdescr_atmpiele,
  iwdescr_fratmpie,
  iwdescr_unsatura,
  iwdescr_funsatura,
  iwdescr_ringisol,
  iwdescr_isolrc,
  iwdescr_isolhtrc,
  iwdescr_erichsct,
  iwdescr_aiercsct,
  iwdescr_lercsct,
  iwdescr_faiercst,
  iwdescr_avcon,
  iwdescr_avchcon,
  iwdescr_avalcon,
  iwdescr_platt,
  iwdescr_weiner,
  iwdescr_crowding,
  iwdescr_fcrowdng,
  iwdescr_halogen,
  iwdescr_halogena,
  iwdescr_bigatom,
  iwdescr_fbigatom,
  iwdescr_csp3,
  iwdescr_fcsp3,
  iwdescr_fccsp3,
  iwdescr_csp3_chain,
  iwdescr_aromc,
  iwdescr_aliphc,
  iwdescr_numcdb,
  iwdescr_totdbsub,
  iwdescr_avcdbsub,
  iwdescr_nflxchn,
  iwdescr_atflxchn,
  iwdescr_faflxchn,
  iwdescr_fnflxchn,
  iwdescr_lflxchn,
  iwdescr_avflxchn,
  iwdescr_rkentrpy,
  iwdescr_nconjgsc,
  iwdescr_atincnjs,
  iwdescr_mxcnjscz,
  iwdescr_cinconjs,
  iwdescr_brunsneg,
  iwdescr_brunspos,
  iwdescr_formal_charge,
  iwdescr_brunsacc,
  iwdescr_brnsdual,
  iwdescr_brunsdon,
  iwdescr_brunshbdsum,
  iwdescr_nplus,
  iwdescr_nminus,
  iwdescr_muldiam,
  iwdescr_rad,
  iwdescr_mulrad,
  iwdescr_tm,
  iwdescr_tg3,
  iwdescr_ishape,
  iwdescr_maxdrng,
  iwdescr_maxdarom,
  iwdescr_maxdhtro, 
  iwdescr_maxdons, 
  iwdescr_avebbtwn,
  iwdescr_normbbtwn,
  iwdescr_compact,
  iwdescr_nolp,
  iwdescr_avdcentre,
  iwdescr_stddcentre,
  iwdescr_centre3,
  iwdescr_centre3h,
  iwdescr_mh3b,
  iwdescr_cntrdgncy,
  iwdescr_cntrdshell1,
  iwdescr_cntrdshell2,
  iwdescr_cntrdshell3,
  iwdescr_aveshell1,
  iwdescr_aveshell2,
  iwdescr_aveshell3,
  iwdescr_maxshell3,
  iwdescr_nnsssrng,
  iwdescr_nrings3,
  iwdescr_nrings4,
  iwdescr_nrings5,
  iwdescr_nrings6,
  iwdescr_nrings7,
  iwdescr_nrings8,       // must keep this in sync with MAX_RING_SIZE
  iwdescr_rsarom1,
  iwdescr_rsarom2,
  iwdescr_rsarom3,
  iwdescr_rsaliph1,
  iwdescr_rsaliph2,
  iwdescr_rsaliph3,
  iwdescr_rsaliph4,
  iwdescr_rssys1,
  iwdescr_rssys2,
  iwdescr_rssys3,
  iwdescr_rssys4,
  iwdescr_rssys5,
  iwdescr_rssys6,
  iwdescr_rssys7,
  iwdescr_rssys8,
  iwdescr_rssys9,
  iwdescr_ar5,
  iwdescr_ar6,
  iwdescr_al5,
  iwdescr_al6,
  iwdescr_fsdrng5l5l,
  iwdescr_fsdrng5l5r,
  iwdescr_fsdrng5r5r,
  iwdescr_fsdrng5r6l,
  iwdescr_fsdrng5l6r,
  iwdescr_fsdrng5l6l,
  iwdescr_fsdrng5r6r,
  iwdescr_fsdrng6r6r,
  iwdescr_fsdrng6l6r,
  iwdescr_fsdrng6l6l,
  iwdescr_fsdrngarar,
  iwdescr_fsdrngalar,
  iwdescr_fsdrngalal,
  iwdescr_nspiro,
  iwdescr_nchiral,
  iwdescr_nvrtspsa,
#ifdef MCGOWAN
  iwdescr_mcgowan,
#endif
  iwdescr_acmbe,
  iwdescr_cmr,
  iwdescr_cd4ring,
  iwdescr_cd4chain,
  iwdescr_frsub,
  iwdescr_frssub,
  iwdescr_arorthoring,
  iwdescr_alorthoring,
  iwdescr_bbr1,
  iwdescr_bbr2,
  iwdescr_bbr3,
  iwdescr_bbr4,
  iwdescr_bbr5,
  iwdescr_bbr6,
  iwdescr_sboradjf,
  iwdescr_dboradjf,
  iwdescr_hcount,
  iwdescr_hperatom,
  iwdescr_ro5_ohnh,
  iwdescr_ro5_on,        // the last one which will always be computed

// complexity descriptors

  iwdescr_nsfsdsys,
  iwdescr_rnginsfs,
  iwdescr_lgstrfsy,
  iwdescr_htrcsfsy,
  iwdescr_mxhtsfsy,

  iwdescr_npfsdsys,
  iwdescr_rnginpfs,
  iwdescr_lgplnfsy,
  iwdescr_htrcpfsy,
  iwdescr_mxhtpfsy,

  iwdescr_aamind,
  iwdescr_aa2mind,
  iwdescr_aaave,
  iwdescr_admind,
  iwdescr_ad2mind,
  iwdescr_adave,
  iwdescr_ddmind,
  iwdescr_dd2mind,
  iwdescr_ddave,

// symmetry related descriptors

  iwdescr_symmatom,
  iwdescr_fsymmatom,
  iwdescr_lsepsymatom,
  iwdescr_flsepsymatom,
  iwdescr_maxsymmclass,

  iwdescr_maxpsymd,
  iwdescr_fmaxpsymd,
  iwdescr_maxpsymdmean,
  iwdescr_psymdnumzero,

  iwdescr_alogp,
  iwdescr_xlogp,

// Ramey descriptors

  iwdescr_obalance,
  iwdescr_rmync,
  iwdescr_rmynn,
  iwdescr_rmyno,
  iwdescr_rmynf,
  iwdescr_rmyns,
  iwdescr_rmyncl,
  iwdescr_rmynbr,
  iwdescr_rmyni,
  iwdescr_rmy_heavy_halogen   // remember to update NUMBER_DESCRIPTORS below.
};

/*
  In order to keep the array large enough, make sure you use the last
  member of the enumeration here
*/

#define NUMBER_DESCRIPTORS (iwdescr_rmy_heavy_halogen +  1)

/*
  Descriptor names
*/

class Descriptor : public Set_or_Unset<float>, public Accumulator<float>
{
  private:
    IWString _name;

    int _active;

//  When accumulating statistics, we keep track of the number of times we are zero

    int _zero_value_count;

//  Some descriptors have default values - the value they take when there are no
//  rings for example

    Set_or_Unset<float> _default_value;

//  If we are generating fingerprints, we need to know how to bucketise
//  these values, and the number of replicates

    int   _fingerprint_replicates;
    float _min;
    float _max;
    float _dy;

//  In calibration we find that some fingerprints do really well as individuals

    int _best_fingerprint;

  public:
    Descriptor();

    void set_name(const char * newname);
    const IWString& name() const {
      return _name;
    }

    void set_best_fingerprint(int s) { _best_fingerprint = s;}
    int  best_fingerprint() const { return _best_fingerprint;}

    void set_min_max_dy(float v1, float v2, float v3) { _min = v1; _max = v2; _dy = v3;}
    void set_min_max_resolution(float v1, float v2, int r);
    int set_range(float dmin, float dmax);

    int produce_fingerprint() const { return _fingerprint_replicates;}
    void set_produce_fingerprint(int s) { _fingerprint_replicates = s;}

    int bit_replicates() const {return _fingerprint_replicates;}

    void produce_fingerprint(int bitnum, Sparse_Fingerprint_Creator &) const;

    const IWString & descriptor_name() const { return _name;}

    int active() const { return _active;}

    int report_statistics(std::ostream &) const;

    void update_statistics();

    void set_default_value(float d) { _default_value.set(d);}

    void reset();

    void set(int s);
    void set(float s) { Set_or_Unset<float>::set(s);}
    void set(double s) { Set_or_Unset<float>::set(static_cast<float>(s));}

    Descriptor operator++ (int);

    Descriptor& operator= (int s) {
      set(s);
      return *this;
    }
    Descriptor& operator= (float s) {
      set(s);
      return *this;
    }
    Descriptor& operator= (double s) {
      set(s);
      return *this;
    }
};

Descriptor::Descriptor()
{
  _active = 0;

  _zero_value_count = 0;

  _fingerprint_replicates = 0;

  _best_fingerprint = 0;

  return;
}

void
Descriptor::set_name(const char * newname) {
  _active = 1;

  if (descriptor_prefix.empty()) {
    _name = newname;
    return;
  }

  _name = descriptor_prefix;
  _name << newname;
  return;
}

int
Descriptor::set_range(float dmin, float dmax) {
  if (dmax < dmin) {
    cerr << "Descriptor::set_range:invalid range " << dmin << ',' << dmax << '\n';
    return 0;
  }

  _min = dmin;
  _max = dmax;

  return 1;
}

int
Descriptor::report_statistics (std::ostream & os) const
{
  const int nsamples = Accumulator<float>::n();
  os << "Descriptor '" << _name << "' " << nsamples << " values sampled\n";
  if (nsamples == 0) {
    return 1;
  }

  os << " between " << Accumulator<float>::minval() << " and " << Accumulator<float>::maxval();
  if (nsamples > 1)
    os << " ave " << Accumulator<float>::average() << " sd " << Accumulator<float>::variance();

  os << ' ' << _zero_value_count << " instances of zero\n";

  return os.good();
}

void
Descriptor::produce_fingerprint (int bitnum, Sparse_Fingerprint_Creator & sfc) const
{
  if (0 == _fingerprint_replicates)
  {
    cerr << "Descriptor::produce_fingerprint:not initialised '" << _name << "'\n";
    return;
  }

  float v;

  if (! Set_or_Unset<float>::value(v))
    return;

  int c;

  if (v <= _min)
    c = 1;
  else
  {
    if (v >= _max)
      v = _max;

    c = static_cast<int>((v - _min) / _dy + 0.5f);

    c++;   // ensure non zero
  }
  
//cerr << "Descriptor::produce_fingerprint: descriptor " << _name << " value " << v << " bit " << bitnum << " value " << c << " rep " << _fingerprint_replicates << '\n';

  if (1 == _fingerprint_replicates)    // presumably a common case
  {
    sfc.hit_bit(bitnum, c);
    return;
  }

  for (int i = 0; i < _fingerprint_replicates; ++i)
  {
    sfc.hit_bit(bitnum + i, c);
  }

  return;
}

void
Descriptor::update_statistics()
{
  float v;
  if (! Set_or_Unset<float>::value(v))
  {
    if (verbose > 2)
      cerr << "Descriptor::update_statistics: '" << _name << "' not set\n";
    return;
  }

  Accumulator<float>::extra(v);

  if (0.0 == v)
    _zero_value_count++;

  return;
}

void
Descriptor::reset()
{
  float f;
  if (_default_value.value(f))
  {
    Set_or_Unset<float>::set(f);
  }
  else
    Set_or_Unset<float>::unset();

  return;
}

void
Descriptor::set_min_max_resolution (float v1, float v2, int r)
{
  assert (v1 < v2);

  _min = v1;
  _max = v2;

  _dy = (v2 - v1) / static_cast<float>(r);

  return;
}

void
Descriptor::set (int s)
{
  Set_or_Unset<float>::set(static_cast<float>(s));

  return;
}

static void
fill_descriptor_extremeties (Descriptor * d,
                             const int resolution)
{
  d[iwdescr_natoms].set_min_max_resolution(3.0f, 50.0f, resolution);
  d[iwdescr_nrings].set_min_max_resolution(0.0f, 9.35196f, resolution);
  d[iwdescr_nrings3].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_nrings4].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_nrings5].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_nrings6].set_min_max_resolution(0.0f, 7.76418f, resolution);
  d[iwdescr_nrings7].set_min_max_resolution(0.0f, 3.0f, resolution);
  d[iwdescr_nrings8].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_nelem].set_min_max_resolution(2.0f, 8.0f, resolution);
  d[iwdescr_amw].set_min_max_resolution(100.0f, 600.0f, resolution);           // manually overwritten
  d[iwdescr_ncon1].set_min_max_resolution(0.0f, 14.95436f, resolution);
  d[iwdescr_fncon1].set_min_max_resolution(0.0f, 0.8f, resolution);
  d[iwdescr_ncon2].set_min_max_resolution(0.0f, 36.7907f, resolution);
  d[iwdescr_fncon2].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_ncon3].set_min_max_resolution(0.0f, 21.12081f, resolution);
  d[iwdescr_fncon3].set_min_max_resolution(0.0f, 0.7778f, resolution);
  d[iwdescr_ncon4].set_min_max_resolution(0.0f, 3.904072f, resolution);
  d[iwdescr_fncon4].set_min_max_resolution(0.0f, 0.4545f, resolution);
  d[iwdescr_frhc].set_min_max_resolution(0.2f, 1.0f, resolution);
  d[iwdescr_mltbd].set_min_max_resolution(0.0f, 45.67025f, resolution);
  d[iwdescr_fmltbd].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_chmltbd].set_min_max_resolution(0.0f, 7.63907f, resolution);
  d[iwdescr_fchmltbd].set_min_max_resolution(0.0f, 0.75f, resolution);
  d[iwdescr_rgmltbd].set_min_max_resolution(0.0f, 44.0818f, resolution);
  d[iwdescr_frgmltbd].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_dcca].set_min_max_resolution(0.0f, 14.28554f, resolution);
  d[iwdescr_fdcca].set_min_max_resolution(0.0f, 0.96f, resolution);
  d[iwdescr_mxdst].set_min_max_resolution(1.0f, 32.1864f, resolution);
  d[iwdescr_fmxdst].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_mxsdlp].set_min_max_resolution(0.0f, 8.80433f, resolution);
  d[iwdescr_avsdlp].set_min_max_resolution(0.0f, 6.5f, resolution);
  d[iwdescr_mxsdlprl].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_mdallp].set_min_max_resolution(0.0f, 15.0f, resolution);
  d[iwdescr_fmdallp].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_fdiffallp].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_harary].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_rotbond].set_min_max_resolution(0.0f, 19.3958f, resolution);
  d[iwdescr_frotbond].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_ringatom].set_min_max_resolution(0.0f, 19.3958f, resolution);
  d[iwdescr_rhacnt].set_min_max_resolution(0.0f, 10.70095f, resolution);
  d[iwdescr_rhaf].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_frafus].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_rngatmf].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_aroma].set_min_max_resolution(0.0f, 42.5879f, resolution);
  d[iwdescr_aromha].set_min_max_resolution(0.0f, 9.12088f, resolution);
  d[iwdescr_fraromha].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_aromdens].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_ch2].set_min_max_resolution(0.0f, 19.30857f, resolution);
  d[iwdescr_chmltbd].set_min_max_resolution(0.0f, 7.63907f, resolution);
  d[iwdescr_ch2].set_min_max_resolution(0.0f, 19.30857f, resolution);
  d[iwdescr_ch].set_min_max_resolution(0.0f, 36.0265f, resolution);
  d[iwdescr_htroatom].set_min_max_resolution(0.0f, 18.7231f, resolution);
  d[iwdescr_htroaf].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_nrgnhlht].set_min_max_resolution(0.0f, 13.73698f, resolution);
  d[iwdescr_ohsh].set_min_max_resolution(0.0f, 8.10969f, resolution);
  d[iwdescr_co2h].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_amine].set_min_max_resolution(0.0f, 5.253233f, resolution);
  d[iwdescr_pyridine].set_min_max_resolution(0.0f, 6.031777f, resolution);
  d[iwdescr_pyrrole].set_min_max_resolution(0.0f, 5.0f, resolution);
  d[iwdescr_hacts].set_min_max_resolution(0.0f, 10.11104f, resolution);
  d[iwdescr_hdons].set_min_max_resolution(0.0f, 5.309234f, resolution);
  d[iwdescr_hduals].set_min_max_resolution(0.0f, 2.995694f, resolution);
  d[iwdescr_mhr].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_mxhrf].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_mnhrf].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_lrsysz].set_min_max_resolution(0.0f, 5.19944f, resolution);
  d[iwdescr_srsz].set_min_max_resolution(3.0f, 8.0f, resolution);    // fix some time
  d[iwdescr_lrsz].set_min_max_resolution(4.0f, 8.0f, resolution);    // fix some time
  d[iwdescr_rng7atoms].set_min_max_resolution(4.0f, 16.0f, resolution);    // fix some time
  d[iwdescr_nrsyscmr].set_min_max_resolution(0.0f, 2.33944f, resolution);
  d[iwdescr_mars].set_min_max_resolution(0.0f, 21.39167f, resolution);
  d[iwdescr_frspch].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_spchtro].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_rbfrspch].set_min_max_resolution(0.0f, 0.9524f, resolution);
  d[iwdescr_satspcha].set_min_max_resolution(0.0f, 22.73322f, resolution);
  d[iwdescr_unsatspcha].set_min_max_resolution(0.0f, 12.77837f, resolution);
  d[iwdescr_fsatspcha].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_scaffoldbranches].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_nrnspch].set_min_max_resolution(0.0f, 6.053596f, resolution);
  d[iwdescr_fnrnspc].set_min_max_resolution(0.0f, 0.6842f, resolution);
  d[iwdescr_trmnlrng].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_intrnlrng].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_rng2spch].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_rng2bridge].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_rcj].set_min_max_resolution(0.0f, 13.8975f, resolution);
  d[iwdescr_rchj].set_min_max_resolution(0.0f, 9.34531f, resolution);
  d[iwdescr_amrcj].set_min_max_resolution(0.0f, 12.52716f, resolution);
  d[iwdescr_alrcj].set_min_max_resolution(0.0f, 8.06624f, resolution);
  d[iwdescr_pbcount].set_min_max_resolution(0.0f, 30.0879f, resolution);
  d[iwdescr_frpbond].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_nonpbond].set_min_max_resolution(0.0f, 30.0879f, resolution);
  d[iwdescr_pbarom].set_min_max_resolution(0.0f, 15.23014f, resolution);
  d[iwdescr_npbarom].set_min_max_resolution(0.0f, 37.5792f, resolution);
  d[iwdescr_pbunset].set_min_max_resolution(0.0f, 7.42857f, resolution);
  d[iwdescr_dvinylb].set_min_max_resolution(0.0f, 6.20917f, resolution);
  d[iwdescr_ringsys].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_arring].set_min_max_resolution(0.0f, 8.00442f, resolution);
  d[iwdescr_alring].set_min_max_resolution(0.0f, 5.147917f, resolution);
  d[iwdescr_excybond].set_min_max_resolution(0.0f, 16.91774f, resolution);
  d[iwdescr_excydbond].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_excydscon].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_excydsconh].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_excydscondon].set_min_max_resolution(0.0f, 10.0f, resolution);
//d[iwdescr_scra].set_min_max_resolution(0.0f, 20.0f, resolution);
//d[iwdescr_scrha].set_min_max_resolution(0.0f, 20.0f, resolution);
//d[iwdescr_scrd].set_min_max_resolution(0.0f, 10.0f, resolution);
  d[iwdescr_atmpiele].set_min_max_resolution(0.0f, 47.3891f, resolution);
  d[iwdescr_fratmpie].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_unsatura].set_min_max_resolution(0.0f, 14.46304f, resolution);
  d[iwdescr_funsatura].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_ringisol].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_isolrc].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_isolhtrc].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_erichsct].set_min_max_resolution(0.0f, 10.75268f, resolution);
  d[iwdescr_aiercsct].set_min_max_resolution(0.0f, 50.0f, resolution);
  d[iwdescr_lercsct].set_min_max_resolution(0.0f, 43.18605f, resolution);
  d[iwdescr_faiercst].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_avcon].set_min_max_resolution(1.333f, 2.783f, resolution);
  d[iwdescr_avchcon].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_avalcon].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_platt].set_min_max_resolution(0.0f, 30.0f, resolution);
  d[iwdescr_weiner].set_min_max_resolution(0.0f, 100.0f, resolution);
  d[iwdescr_crowding].set_min_max_resolution(0.0f, 35.21212f, resolution);
  d[iwdescr_fcrowdng].set_min_max_resolution(0.0f, 2.364f, resolution);
  d[iwdescr_halogen].set_min_max_resolution(0.0f, 7.139123f, resolution);
  d[iwdescr_halogena].set_min_max_resolution(0.0f, 4.72754f, resolution);
  d[iwdescr_halogena].set_min_max_resolution(0.0f, 4.72754f, resolution);
  d[iwdescr_bigatom].set_min_max_resolution(0.0f, 5.248915f, resolution);
  d[iwdescr_fbigatom].set_min_max_resolution(0.0f, 0.875f, resolution);
  d[iwdescr_csp3].set_min_max_resolution(0.0f, 26.57935f, resolution);
  d[iwdescr_csp3_chain].set_min_max_resolution(0.0f, 17.05735f, resolution);
  d[iwdescr_fcsp3].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_fccsp3].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_csp3_chain].set_min_max_resolution(0.0f, 17.05735f, resolution);
  d[iwdescr_aromc].set_min_max_resolution(0.0f, 38.29105f, resolution);
  d[iwdescr_aliphc].set_min_max_resolution(0.0f, 29.31654f, resolution);
  d[iwdescr_numcdb].set_min_max_resolution(0.0f, 7.57898f, resolution);
  d[iwdescr_totdbsub].set_min_max_resolution(0.0f, 13.82252f, resolution);
  d[iwdescr_avcdbsub].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_nflxchn].set_min_max_resolution(0.0f, 8.05117f, resolution);
  d[iwdescr_atflxchn].set_min_max_resolution(0.0f, 18.79481f, resolution);
  d[iwdescr_faflxchn].set_min_max_resolution(0.0f, 0.94f, resolution);
  d[iwdescr_fnflxchn].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_lflxchn].set_min_max_resolution(0.0f, 10.29103f, resolution);
  d[iwdescr_avflxchn].set_min_max_resolution(0.0f, 8.4062f, resolution);
  d[iwdescr_rkentrpy].set_min_max_resolution(0.0f, 11.74898f, resolution);
  d[iwdescr_nconjgsc].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_atincnjs].set_min_max_resolution(0.0f, 7.040911f, resolution);
  d[iwdescr_mxcnjscz].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_cinconjs].set_min_max_resolution(0.0f, 27.99497f, resolution);
  d[iwdescr_brunsneg].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_brunspos].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_formal_charge].set_min_max_resolution(0.0f, 10.19615f, resolution);
  d[iwdescr_brunsacc].set_min_max_resolution(0.0f, 10.19615f, resolution);
  d[iwdescr_brnsdual].set_min_max_resolution(0.0f, 2.761233f, resolution);
  d[iwdescr_brunsdon].set_min_max_resolution(0.0f, 6.24665f, resolution);
  d[iwdescr_brunshbdsum].set_min_max_resolution(0.0f, 9.24665f, resolution);
  d[iwdescr_nplus].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_nminus].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_muldiam].set_min_max_resolution(2.0f, 9.4635f, resolution);
  d[iwdescr_rad].set_min_max_resolution(1.0f, 16.37685f, resolution);
  d[iwdescr_mulrad].set_min_max_resolution(1.0f, 6.10493f, resolution);
  d[iwdescr_tm].set_min_max_resolution(0.0f, 8.65778f, resolution);
  d[iwdescr_tg3].set_min_max_resolution(0.0f, 5.955678f, resolution);
  d[iwdescr_ishape].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_maxdrng].set_min_max_resolution(0.0f, 32.90875f, resolution);
  d[iwdescr_maxdarom].set_min_max_resolution(0.0f, 33.0f, resolution);
  d[iwdescr_maxdhtro].set_min_max_resolution(0.0f, 33.0f, resolution);
  d[iwdescr_maxdons].set_min_max_resolution(0.0f, 33.0f, resolution);
  d[iwdescr_avebbtwn].set_min_max_resolution(1.0f, 11.81314f, resolution);
  d[iwdescr_normbbtwn].set_min_max_resolution(0.1089f, 0.4444f, resolution);
  d[iwdescr_compact].set_min_max_resolution(0.02f, 0.7955f, resolution);
  d[iwdescr_nolp].set_min_max_resolution(0.0f, 32.4906f, resolution);
  d[iwdescr_avdcentre].set_min_max_resolution(0.0f, 6.0, resolution);
  d[iwdescr_stddcentre].set_min_max_resolution(0.0f, 32.4906f, resolution);
  d[iwdescr_centre3].set_min_max_resolution(0.0f, 20.0f, resolution);
  d[iwdescr_centre3h].set_min_max_resolution(0.0f, 12.0f, resolution);
  d[iwdescr_mh3b].set_min_max_resolution(0.0f, 12.0f, resolution);
  d[iwdescr_cntrdgncy].set_min_max_resolution(1.0f, 6.0f, resolution);
  d[iwdescr_cntrdshell1].set_min_max_resolution(1.0f, 8.0f, resolution);
  d[iwdescr_cntrdshell2].set_min_max_resolution(1.0f, 14.0f, resolution);
  d[iwdescr_cntrdshell3].set_min_max_resolution(1.0f, 14.0f, resolution);
  d[iwdescr_aveshell1].set_min_max_resolution(1.0f, 3.0f, resolution);
  d[iwdescr_aveshell2].set_min_max_resolution(2.0f, 7.0f, resolution);
  d[iwdescr_aveshell3].set_min_max_resolution(3.0f, 12.0f, resolution);
  d[iwdescr_maxshell3].set_min_max_resolution(3.0f, 12.0f, resolution);
  d[iwdescr_nnsssrng].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_nrings3].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_nrings4].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_nrings5].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_nrings6].set_min_max_resolution(0.0f, 7.76418f, resolution);
  d[iwdescr_nrings7].set_min_max_resolution(0.0f, 3.0f, resolution);
  d[iwdescr_nrings8].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_rsarom1].set_min_max_resolution(0.0f, 5.833233f, resolution);
  d[iwdescr_rsarom2].set_min_max_resolution(0.0f, 6.561037f, resolution);
  d[iwdescr_rsarom3].set_min_max_resolution(0.0f, 9.0f, resolution);
  d[iwdescr_rsaliph1].set_min_max_resolution(0.0f, 3.628616f, resolution);
  d[iwdescr_rsaliph2].set_min_max_resolution(0.0f, 3.79946f, resolution);
  d[iwdescr_rsaliph3].set_min_max_resolution(0.0f, 2.328854f, resolution);
  d[iwdescr_rsaliph4].set_min_max_resolution(0.0f, 2.0307791f, resolution);
  d[iwdescr_rssys1].set_min_max_resolution(0.0f, 3.223665f, resolution);
  d[iwdescr_rssys2].set_min_max_resolution(0.0f, 2.76466f, resolution);
  d[iwdescr_rssys3].set_min_max_resolution(0.0f, 2.636825f, resolution);
  d[iwdescr_rssys4].set_min_max_resolution(0.0f, 3.131377f, resolution);
  d[iwdescr_rssys5].set_min_max_resolution(0.0f, 1.5407535f, resolution);
  d[iwdescr_rssys6].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_rssys7].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_rssys8].set_min_max_resolution(0.0f, 9.0f, resolution);
  d[iwdescr_rssys9].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_ar5].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_ar6].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_al5].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_al6].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_fsdrng5l5l].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_fsdrng5l5r].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_fsdrng5r5r].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_fsdrng5r6l].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_fsdrng5l6r].set_min_max_resolution(0.0f, 3.0f, resolution);
  d[iwdescr_fsdrng5l6l].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_fsdrng5r6r].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_fsdrng6r6r].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_fsdrng6l6r].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_fsdrng6l6l].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_fsdrngarar].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_fsdrngalar].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_fsdrngalal].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_nchiral].set_min_max_resolution(0.0f, 5.821545f, resolution);
  d[iwdescr_nvrtspsa].set_min_max_resolution(0.0f, 220.2642f, resolution);
#ifdef MCGOWAN
  d[iwdescr_mcgowan].set_min_max_resolution(50.0f, 900.0f, resolution);
#endif
  d[iwdescr_acmbe].set_min_max_resolution(3.2f, 113.3805f, resolution);
  d[iwdescr_cmr].set_min_max_resolution(1.32f, 23.31713f, resolution);
  d[iwdescr_cd4ring].set_min_max_resolution(0.0f, 2.134513f, resolution);
  d[iwdescr_cd4chain].set_min_max_resolution(0.0f, 2.508958f, resolution);
  d[iwdescr_frsub].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_frssub].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_arorthoring].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_alorthoring].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_bbr1].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_bbr2].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_bbr3].set_min_max_resolution(0.0f, 5.0f, resolution);
  d[iwdescr_bbr4].set_min_max_resolution(0.0f, 4.0f, resolution);
  d[iwdescr_bbr5].set_min_max_resolution(0.0f, 3.0f, resolution);
  d[iwdescr_bbr6].set_min_max_resolution(0.0f, 3.0f, resolution);
  d[iwdescr_sboradjf].set_min_max_resolution(0.0f, 4.618648f, resolution);
  d[iwdescr_dboradjf].set_min_max_resolution(0.0f, 2.224352f, resolution);
  d[iwdescr_hcount].set_min_max_resolution(0.0f, 20.0f, resolution);    // just a guess
  d[iwdescr_hperatom].set_min_max_resolution(0.0f, 1.0f, resolution);    // just a guess
  d[iwdescr_ro5_ohnh].set_min_max_resolution(0.0f, 6.058745f, resolution);
  d[iwdescr_ro5_on].set_min_max_resolution(0.0f, 15.53816f, resolution);

  d[iwdescr_symmatom].set_min_max_resolution(0.0f, 20.0f, resolution);
  d[iwdescr_fsymmatom].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_lsepsymatom].set_min_max_resolution(0.0f, 30.0f, resolution);
  d[iwdescr_flsepsymatom].set_min_max_resolution(0.0f, 1.0f, resolution);
  d[iwdescr_maxsymmclass].set_min_max_resolution(0.0f, 12.0f, resolution);

  d[iwdescr_maxpsymd].set_min_max_resolution(0, 30.0f, resolution);
  d[iwdescr_fmaxpsymd].set_min_max_resolution(0, 30.0f, resolution);
  d[iwdescr_maxpsymdmean].set_min_max_resolution(0, 10.0f, resolution);
  d[iwdescr_psymdnumzero].set_min_max_resolution(0, 10.0f, resolution);

  d[iwdescr_rmync].set_min_max_resolution(0.0f, 48.2944f, resolution);
  d[iwdescr_rmyncl].set_min_max_resolution(0.0f, 3.300221f, resolution);
  d[iwdescr_rmynn].set_min_max_resolution(0.0f, 10.8482f, resolution);
  d[iwdescr_rmyno].set_min_max_resolution(0.0f, 10.46232f, resolution);
  d[iwdescr_rmynf].set_min_max_resolution(0.0f, 6.060531f, resolution);
  d[iwdescr_rmyns].set_min_max_resolution(0.0f, 3.578964f, resolution);
  d[iwdescr_rmyncl].set_min_max_resolution(0.0f, 3.300221f, resolution);
  d[iwdescr_rmynbr].set_min_max_resolution(0.0f, 9.0f, resolution);
  d[iwdescr_rmyni].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_rmy_heavy_halogen].set_min_max_resolution(0.0f, 10.0f, resolution);

  d[iwdescr_alogp].set_min_max_resolution(-4.0f, 10.0f, resolution);
  d[iwdescr_xlogp].set_min_max_resolution(-4.0f, 10.0f, resolution);

  return;
}

/*Descriptor
Descriptor::operator++ (int)
{
  Accumulator<float>::extra(static_cast<float>(1.0));

  return *this;
}*/

static Descriptor * descriptor = nullptr;

// Linear search through the global descriptors array for a match to `dname`.
static std::optional<Descriptor*>
GetDescriptor(const Descriptor* descriptors,
              int number_descriptors,
              const std::string& dname) {
  IWString name(dname);
  for (int i = 0; i < number_descriptors; ++i) {
    if (descriptor[i].name() == name) {
      return descriptor + i;
    }
  }

  // Try with a prefix.
  name = descriptor_prefix << dname;
  for (int i = 0; i < number_descriptors; ++i) {
    if (descriptor[i].name() == name) {
      return descriptor + i;
    }
  }

  return std::nullopt;
}

// Return the index of the Descriptor having name `dname`.
static int
descriptor_name_2_index(const Descriptor * d,
                        int n,
                        const const_IWSubstring & dname)
{
  if (descriptor_prefix.empty()) {
  } else if (! dname.starts_with(descriptor_prefix)) {
    IWString s;
    s << descriptor_prefix << dname;
    return descriptor_name_2_index(d, n, s);
  }

  for (int i = 0; i < n; i++) {
    if (dname == descriptor[i].descriptor_name())
      return i;
  }

  return -1;
}

class Descriptor_Filter
{
  private:
    IWString _descriptor_name;

    int _descriptor_number;

    Min_Max_Specifier<float> _cond;

    int _items_rejected;

  public:
    Descriptor_Filter();

    int build (const const_IWSubstring &);

    int report (std::ostream &) const;

    int items_rejected () const { return _items_rejected;}

    int satisfied (const Descriptor *);
};

Descriptor_Filter::Descriptor_Filter()
{
  _items_rejected = 0;

  return;
}

static int
identify_descriptor (const IWString & dname,
                     const Descriptor * descriptor,
                     int number_descriptors)
{
  IWString mydname;
  if (dname.starts_with(descriptor_prefix))
    mydname = dname;
  else
    mydname << descriptor_prefix << dname;

  for (int i = 0; i < number_descriptors; i++)
  {
    if (mydname == descriptor[i].descriptor_name())
      return i;
  }

  return -1;
}

/*
  Syntax will be

  dname=3 or dname.eq.3
  dname>4 or dname.gt.4
  dname<8 or dname.lt.4
  dname:n1,n2,n3
*/

#define DFILTER_OP_EQUALS 1
#define DFILTER_OP_LT 2
#define DFILTER_OP_GT 3
#define DFILTER_OP_LIST 4

static int
look_for_dot_operators (const const_IWSubstring & s,
                        IWString & dname,
                        int & op,
                        int & begin_numeric)

{
  int n = s.length() - 4;  // min would be dname.op.1

  for (int i = 1; i < n; i++)
  {
    if ('.' != s[i])
      continue;

    if ('.' != s[i + 3])
      return 0;

    if ('e' == s[i+1] && 'q' == s[i+2])
    {
      s.from_to(0, i - 1, dname);
      op = DFILTER_OP_EQUALS;
      begin_numeric = i + 4;
      return 1;
    }
    if ('l' == s[i+1] && 't' == s[i+2])
    {
      s.from_to(0, i - 1, dname);
      op = DFILTER_OP_LT;
      begin_numeric = i + 4;
      return 1;
    }
    if ('g' == s[i+1] && 't' == s[i+2])
    {
      s.from_to(0, i - 1, dname);
      op = DFILTER_OP_GT;
      begin_numeric = i + 4;
      return 1;
    }

    return 0;
  }

  return 0;
}

int
Descriptor_Filter::build (const const_IWSubstring & s)
{
  int op = 0;
  int begin_numeric = -1;

  int n = s.length() - 1;    // min form is dname=3

  for (int i = 1; i < n; i++)
  {
    char c = s[i];

    if ('=' == c)
      op = DFILTER_OP_EQUALS;
    else if ('>' == c)
      op = DFILTER_OP_GT;
    else if ('<' == c)
      op = DFILTER_OP_LT;
    else if (':' == c)
      op = DFILTER_OP_LIST;
    else
      continue;

    s.from_to(0, i - 1, _descriptor_name);
    begin_numeric = i + 1;
    break;
  }

  if (0 == _descriptor_name.length())
    look_for_dot_operators(s, _descriptor_name, op, begin_numeric);

  if (0 == _descriptor_name.length())
  {
    cerr << "Descriptor_Filter::build:no operator found '" << s << "'\n";
    return 0;
  }

  const_IWSubstring numeric_stuff;

  s.from_to(begin_numeric, s.length() - 1, numeric_stuff);

  _descriptor_number = identify_descriptor(_descriptor_name, descriptor, NUMBER_DESCRIPTORS);

  if (_descriptor_number < 0)
  {
    cerr << "Descriptor_Filter::build:cannot find '" << _descriptor_name << "'\n";
    return 0;
  }

  if (DFILTER_OP_LIST != op)
  {
    float v;
    if (! numeric_stuff.numeric_value(v))
    {
      cerr << "Descriptor_Filter::build:invalid numeric '" << numeric_stuff << "'\n";
      return 0;
    }

    if (DFILTER_OP_EQUALS == op) {
      _cond.add(v);
    } else if (DFILTER_OP_LT == op) {
      if (v == 0.0) {
        _cond.set_max(- std::numeric_limits<float>::min());
      } else {
        _cond.set_max(0.999999 * v);
      }
    } else if (DFILTER_OP_GT == op) {
      if (v == 0.0) {
        _cond.set_min(std::numeric_limits<float>::min());
      } else {
        _cond.set_min(1.000001 * v);
      }
    }

    return 1;
  }

  const_IWSubstring token;
  int i = 0;
  while (numeric_stuff.nextword(token, i, ','))
  {
    float v;
    if (! token.numeric_value(v))
    {
      cerr << "Descriptor_Filter::build:non numeric in list specification '" << s << "'\n";
      return 0;
    }

    _cond.add(v);
  }

  return 1;
}

int
Descriptor_Filter::satisfied (const Descriptor * descriptors)
{
  const Descriptor & myd = descriptors[_descriptor_number];

  float v;

  if (! myd.value(v))    // no value available
    return 0;

  if (_cond.matches(v))
    return 1;

  _items_rejected++;

  return 0;
}

int
Descriptor_Filter::report (std::ostream & os) const
{
  os << "Filter on '" << _descriptor_name << "', rejected " << _items_rejected << " items\n";

  return 1;
}

static Descriptor_Filter * filter = nullptr;
static int number_filters = 0;
static IWString smiles_as_input;

static int rejected_by_filters = 0;

static int
allocate_descriptors()
{
  assert (nullptr == descriptor);

  descriptor = new Descriptor[NUMBER_DESCRIPTORS];

  if (nullptr == descriptor)
  {
    cerr << "Memory failure, cannot allocate " << NUMBER_DESCRIPTORS << " descriptors\n";
    return 0;
  }

  descriptor[iwdescr_natoms].set_name("natoms");
  descriptor[iwdescr_nrings].set_name("nrings");
  descriptor[iwdescr_nelem].set_name("nelem");
  descriptor[iwdescr_amw].set_name("amw");
  if (descriptors_to_compute.ncon_descriptors) {
    descriptor[iwdescr_ncon1].set_name("ncon1");
    descriptor[iwdescr_fncon1].set_name("fncon1");
    descriptor[iwdescr_ncon2].set_name("ncon2");
    descriptor[iwdescr_fncon2].set_name("fncon2");
    descriptor[iwdescr_ncon3].set_name("ncon3");
    descriptor[iwdescr_fncon3].set_name("fncon3");
    descriptor[iwdescr_ncon4].set_name("ncon4");
    descriptor[iwdescr_fncon4].set_name("fncon4");
    descriptor[iwdescr_platt].set_name("platt");
  }
  descriptor[iwdescr_frhc].set_name("frhc");
  descriptor[iwdescr_mltbd].set_name("mltbd");
  descriptor[iwdescr_fmltbd].set_name("fmltbd");
  descriptor[iwdescr_chmltbd].set_name("chmltbd");
  descriptor[iwdescr_fchmltbd].set_name("fchmltbd");
  descriptor[iwdescr_rgmltbd].set_name("rgmltbd");
  descriptor[iwdescr_frgmltbd].set_name("frgmltbd");
  descriptor[iwdescr_dcca].set_name("dcca");
  descriptor[iwdescr_fdcca].set_name("fdcca");
  if (descriptors_to_compute.distance_matrix_descriptors) {
    descriptor[iwdescr_mxdst].set_name("mxdst");
    descriptor[iwdescr_fmxdst].set_name("fmxdst");
    descriptor[iwdescr_mxsdlp].set_name("mxsdlp");
    descriptor[iwdescr_avsdlp].set_name("avsdlp");
    descriptor[iwdescr_mxsdlprl].set_name("mxsdlprl");
    descriptor[iwdescr_mdallp].set_name("mdallp");
    descriptor[iwdescr_fmdallp].set_name("fmdallp");
    descriptor[iwdescr_fdiffallp].set_name("fdiffallp");
    descriptor[iwdescr_harary].set_name("harary");
  }
  descriptor[iwdescr_rotbond].set_name("rotbond");
  descriptor[iwdescr_frotbond].set_name("frotbond");
  descriptor[iwdescr_ringatom].set_name("ringatom");
  descriptor[iwdescr_rhacnt].set_name("rhacnt");
  descriptor[iwdescr_rhaf].set_name("rhaf");
  descriptor[iwdescr_frafus].set_name("frafus");
  descriptor[iwdescr_rngatmf].set_name("rngatmf");
  descriptor[iwdescr_aroma].set_name("aroma");
  descriptor[iwdescr_aromha].set_name("aromha");
  descriptor[iwdescr_fraromha].set_name("fraromha");
  descriptor[iwdescr_aromdens].set_name("aromdens");
  descriptor[iwdescr_ch2].set_name("ch2");
  descriptor[iwdescr_ch].set_name("ch");
  descriptor[iwdescr_htroatom].set_name("htroatom");
  descriptor[iwdescr_htroaf].set_name("htroaf");
  descriptor[iwdescr_nrgnhlht].set_name("nrgnhlht");
  descriptor[iwdescr_ohsh].set_name("ohsh");
  descriptor[iwdescr_co2h].set_name("co2h");

  if (descriptors_to_compute.specific_groups) {
    descriptor[iwdescr_amine].set_name("amine");
    descriptor[iwdescr_pyridine].set_name("pyridine");
    descriptor[iwdescr_pyrrole].set_name("pyrrole");
  }

  if (descriptors_to_compute.simple_hbond_descriptors) {
    descriptor[iwdescr_hacts].set_name("hacts");
    descriptor[iwdescr_hdons].set_name("hdons");
    descriptor[iwdescr_hduals].set_name("hduals");
  }
  descriptor[iwdescr_mhr].set_name("mhr");
  descriptor[iwdescr_mxhrf].set_name("mxhrf");
  descriptor[iwdescr_mnhrf].set_name("mnhrf");
  descriptor[iwdescr_lrsysz].set_name("lrsysz");
  descriptor[iwdescr_srsz].set_name("srsz");
  descriptor[iwdescr_lrsz].set_name("lrsz");
  descriptor[iwdescr_rng7atoms].set_name("rng7atoms");
  descriptor[iwdescr_nrsyscmr].set_name("nrsyscmr");
  descriptor[iwdescr_mars].set_name("mars");
  if (descriptors_to_compute.spinach_descriptors) {
    descriptor[iwdescr_frspch].set_name("frspch");
    descriptor[iwdescr_spchtro].set_name("spchtro");
    descriptor[iwdescr_rbfrspch].set_name("rbfrspch");
    descriptor[iwdescr_satspcha].set_name("satspcha");
    descriptor[iwdescr_unsatspcha].set_name("unsatspcha");
    descriptor[iwdescr_fsatspcha].set_name("fsatspcha");
    descriptor[iwdescr_scaffoldbranches].set_name("scaffoldbranches");
    descriptor[iwdescr_nrnspch].set_name("nrnspch");
    descriptor[iwdescr_fnrnspc].set_name("fnrnspc");

    descriptor[iwdescr_trmnlrng].set_name("trmnlrng");
    descriptor[iwdescr_intrnlrng].set_name("intrnlrng");
    descriptor[iwdescr_rng2spch].set_name("rng2spch");
    descriptor[iwdescr_rng2bridge].set_name("rng2bridge");
  }

  if (descriptors_to_compute.ring_chain_descriptors) {
    descriptor[iwdescr_rcj].set_name("rcj");
    descriptor[iwdescr_rchj].set_name("rchj");
    descriptor[iwdescr_amrcj].set_name("amrcj");
    descriptor[iwdescr_alrcj].set_name("alrcj");
  }

  if (descriptors_to_compute.polar_bond_descriptors) {
    descriptor[iwdescr_pbcount].set_name("pbcount");
    descriptor[iwdescr_frpbond].set_name("frpbond");
    descriptor[iwdescr_nonpbond].set_name("nonpbond");
    descriptor[iwdescr_pbarom].set_name("pbarom");
    descriptor[iwdescr_npbarom].set_name("npbarom");
    descriptor[iwdescr_pbunset].set_name("pbunset");
    descriptor[iwdescr_dvinylb].set_name("dvinylb");
  }
  descriptor[iwdescr_ringsys].set_name("ringsys");
  descriptor[iwdescr_arring].set_name("arring");
  descriptor[iwdescr_alring].set_name("alring");
  descriptor[iwdescr_excybond].set_name("excybond");
  descriptor[iwdescr_excydbond].set_name("excydbond");
  descriptor[iwdescr_excydscon].set_name("excydscon");
  descriptor[iwdescr_excydsconh].set_name("excydsconh");
  descriptor[iwdescr_excydscondon].set_name("excydscondon");
//descriptor[iwdescr_scra].set_name("scra");
//descriptor[iwdescr_scrha].set_name("scrha");
//descriptor[iwdescr_scrd].set_name("scrd");
  descriptor[iwdescr_atmpiele].set_name("atmpiele");
  descriptor[iwdescr_fratmpie].set_name("fratmpie");
  descriptor[iwdescr_unsatura].set_name("unsatura");
  descriptor[iwdescr_funsatura].set_name("funsatura");
  descriptor[iwdescr_ringisol].set_name("ringisol");
  descriptor[iwdescr_isolrc].set_name("isolrc");
  descriptor[iwdescr_isolhtrc].set_name("isolhtrc");
  descriptor[iwdescr_erichsct].set_name("erichsct");
  descriptor[iwdescr_aiercsct].set_name("aiercsct");
  descriptor[iwdescr_lercsct].set_name("lercsct");
  descriptor[iwdescr_faiercst].set_name("faiercst");
  descriptor[iwdescr_avcon].set_name("avcon");
  descriptor[iwdescr_avchcon].set_name("avchcon");
  descriptor[iwdescr_avalcon].set_name("avalcon");
  descriptor[iwdescr_platt].set_name("platt");
  descriptor[iwdescr_weiner].set_name("weiner");
  if (descriptors_to_compute.crowding_descriptors) {
    descriptor[iwdescr_crowding].set_name("crowding");
    descriptor[iwdescr_fcrowdng].set_name("fcrowdng");
  }
  descriptor[iwdescr_halogen].set_name("halogen");
  descriptor[iwdescr_halogena].set_name("halogena");
  descriptor[iwdescr_bigatom].set_name("bigatom");
  descriptor[iwdescr_fbigatom].set_name("fbigatom");
  descriptor[iwdescr_csp3].set_name("csp3");
  descriptor[iwdescr_fcsp3].set_name("fcsp3");
  descriptor[iwdescr_fccsp3].set_name("fccsp3");
  descriptor[iwdescr_csp3_chain].set_name("csp3_chain");
  descriptor[iwdescr_aromc].set_name("aromc");
  descriptor[iwdescr_aliphc].set_name("aliphc");
  descriptor[iwdescr_numcdb].set_name("numcdb");
  descriptor[iwdescr_totdbsub].set_name("totdbsub");
  descriptor[iwdescr_avcdbsub].set_name("avcdbsub");
  descriptor[iwdescr_nflxchn].set_name("nflxchn");
  descriptor[iwdescr_atflxchn].set_name("atflxchn");
  descriptor[iwdescr_faflxchn].set_name("faflxchn");
  descriptor[iwdescr_fnflxchn].set_name("fnflxchn");
  descriptor[iwdescr_lflxchn].set_name("lflxchn");
  descriptor[iwdescr_avflxchn].set_name("avflxchn");
  descriptor[iwdescr_rkentrpy].set_name("rkentrpy");
  descriptor[iwdescr_nconjgsc].set_name("nconjgsc");
  descriptor[iwdescr_atincnjs].set_name("atincnjs");
  descriptor[iwdescr_mxcnjscz].set_name("mxcnjscz");
  descriptor[iwdescr_cinconjs].set_name("cinconjs");
  if (descriptors_to_compute.charge_descriptors) {
    descriptor[iwdescr_brunsneg].set_name("brunsneg");
    descriptor[iwdescr_brunspos].set_name("brunspos");
    descriptor[iwdescr_formal_charge].set_name("formal_charge");
  }
  if (descriptors_to_compute.donor_acceptor) {
    descriptor[iwdescr_brunsacc].set_name("brunsacc");
    descriptor[iwdescr_brnsdual].set_name("brnsdual");
    descriptor[iwdescr_brunsdon].set_name("brunsdon");
    descriptor[iwdescr_brunshbdsum].set_name("brunshbdsum");
    descriptor[iwdescr_nplus].set_name("nplus");
      descriptor[iwdescr_nminus].set_name("nminus");
  }

  if (descriptors_to_compute.distance_matrix_descriptors) {
    descriptor[iwdescr_muldiam].set_name("muldiam");
    descriptor[iwdescr_rad].set_name("rad");
    descriptor[iwdescr_mulrad].set_name("mulrad");
    descriptor[iwdescr_tm].set_name("tm");
    descriptor[iwdescr_tg3].set_name("tg3");
    descriptor[iwdescr_ishape].set_name("ishape");
    descriptor[iwdescr_maxdrng].set_name("maxdrng");
    descriptor[iwdescr_maxdarom].set_name("maxdarom");
    descriptor[iwdescr_maxdhtro].set_name("maxdhtro");
    descriptor[iwdescr_maxdons].set_name("maxdons");
    descriptor[iwdescr_avebbtwn].set_name("avebbtwn");
    descriptor[iwdescr_normbbtwn].set_name("normbbtwn");
    descriptor[iwdescr_compact].set_name("compact");
    descriptor[iwdescr_nolp].set_name("nolp");
    descriptor[iwdescr_avdcentre].set_name("avdcentre");
    descriptor[iwdescr_stddcentre].set_name("stddcentre");
    descriptor[iwdescr_centre3].set_name("centre3");
    descriptor[iwdescr_centre3h].set_name("centre3h");
    descriptor[iwdescr_mh3b].set_name("mh3b");
    descriptor[iwdescr_cntrdgncy].set_name("cntrdgncy");
    descriptor[iwdescr_cntrdshell1].set_name("cntrdshell1");
    descriptor[iwdescr_cntrdshell2].set_name("cntrdshell2");
    descriptor[iwdescr_cntrdshell3].set_name("cntrdshell3");

    descriptor[iwdescr_aveshell1].set_name("aveshell1");
    descriptor[iwdescr_aveshell2].set_name("aveshell2");
    descriptor[iwdescr_aveshell3].set_name("aveshell3");
    descriptor[iwdescr_maxshell3].set_name("maxshell3");
  }

  descriptor[iwdescr_nnsssrng].set_name("nnsssrng");
  descriptor[iwdescr_nrings3].set_name("nrings3");
  descriptor[iwdescr_nrings4].set_name("nrings4");
  descriptor[iwdescr_nrings5].set_name("nrings5");
  descriptor[iwdescr_nrings6].set_name("nrings6");
  descriptor[iwdescr_nrings7].set_name("nrings7");
  descriptor[iwdescr_nrings8].set_name("nrings8");

  if (descriptors_to_compute.ring_substitution_descriptors)  {
    descriptor[iwdescr_rsarom1].set_name("rsarom1");
    descriptor[iwdescr_rsarom2].set_name("rsarom2");
    descriptor[iwdescr_rsarom3].set_name("rsarom3");

    descriptor[iwdescr_rsaliph1].set_name("rsaliph1");
    descriptor[iwdescr_rsaliph2].set_name("rsaliph2");
    descriptor[iwdescr_rsaliph3].set_name("rsaliph3");
    descriptor[iwdescr_rsaliph4].set_name("rsaliph4");
  
    descriptor[iwdescr_rssys1].set_name("rssys1");
    descriptor[iwdescr_rssys2].set_name("rssys2");
    descriptor[iwdescr_rssys3].set_name("rssys3");
    descriptor[iwdescr_rssys4].set_name("rssys4");
    descriptor[iwdescr_rssys5].set_name("rssys5");
    descriptor[iwdescr_rssys6].set_name("rssys6");
    descriptor[iwdescr_rssys7].set_name("rssys7");
    descriptor[iwdescr_rssys8].set_name("rssys8");
    descriptor[iwdescr_rssys9].set_name("rssys9");
  }

  if (descriptors_to_compute.ring_fusion_descriptors) {
    descriptor[iwdescr_ar5].set_name("ar5");
    descriptor[iwdescr_ar6].set_name("ar6");
    descriptor[iwdescr_al5].set_name("al5");
    descriptor[iwdescr_al6].set_name("al6");

    descriptor[iwdescr_fsdrng5l5l].set_name("fsdrng5l5l");
    descriptor[iwdescr_fsdrng5l5r].set_name("fsdrng5l5r");
    descriptor[iwdescr_fsdrng5r5r].set_name("fsdrng5r5r");
    descriptor[iwdescr_fsdrng5l6l].set_name("fsdrng5l6l");
    descriptor[iwdescr_fsdrng5l6r].set_name("fsdrng5l6r");
    descriptor[iwdescr_fsdrng5r6l].set_name("fsdrng5r6l");
    descriptor[iwdescr_fsdrng5r6r].set_name("fsdrng5r6r");
    descriptor[iwdescr_fsdrng6r6r].set_name("fsdrng6r6r");
    descriptor[iwdescr_fsdrng6l6r].set_name("fsdrng6l6r");
    descriptor[iwdescr_fsdrng6l6l].set_name("fsdrng6l6l");

    descriptor[iwdescr_fsdrngarar].set_name("fsdrngarar");
    descriptor[iwdescr_fsdrngalar].set_name("fsdrngalar");
    descriptor[iwdescr_fsdrngalal].set_name("fsdrngalal");
  }

  descriptor[iwdescr_nchiral].set_name("nchiral");
  if (descriptors_to_compute.psa) {
    descriptor[iwdescr_nvrtspsa].set_name("nvrtspsa");
  }
#ifdef MCGOWAN
  if (descriptors_to_compute.mcgowan) {
    descriptor[iwdescr_mcgowan].set_name("mcgowan");
  }
#endif // MCGOWAN
  if (descriptors_to_compute.charge_descriptors) {
    descriptor[iwdescr_acmbe].set_name("acmbe");
  }
  descriptor[iwdescr_cmr].set_name("cmr");
  descriptor[iwdescr_cd4ring].set_name("cd4ring");
  descriptor[iwdescr_cd4chain].set_name("cd4chain");
  if (descriptors_to_compute.ring_substitution_ratio_descriptors) {
    descriptor[iwdescr_frsub].set_name("frsub");
    descriptor[iwdescr_frssub].set_name("frssub");
    descriptor[iwdescr_arorthoring].set_name("arorthoring");
    descriptor[iwdescr_alorthoring].set_name("alorthoring");
  }
  if (descriptors_to_compute.bonds_between_rings) {
    descriptor[iwdescr_bbr1].set_name("bbr1");
    descriptor[iwdescr_bbr2].set_name("bbr2");
    descriptor[iwdescr_bbr3].set_name("bbr3");
    descriptor[iwdescr_bbr4].set_name("bbr4");
    descriptor[iwdescr_bbr5].set_name("bbr5");
    descriptor[iwdescr_bbr6].set_name("bbr6");
  }
  if (descriptors_to_compute.adjacent_ring_fusion_descriptors) {
    descriptor[iwdescr_sboradjf].set_name("sboradjf");
    descriptor[iwdescr_dboradjf].set_name("dboradjf");
  }
  descriptor[iwdescr_hcount].set_name("hcount");
  descriptor[iwdescr_hperatom].set_name("hperatom");
  descriptor[iwdescr_ro5_ohnh].set_name("ro5_ohnh");
  descriptor[iwdescr_ro5_on].set_name("ro5_on");    // the last one which will always be computed

  if (descriptors_to_compute.donor_acceptor && min_hbond_feature_separation > 0)
  {
    descriptor[iwdescr_aamind].set_name("aamind");
    descriptor[iwdescr_aa2mind].set_name("aa2mind");
    descriptor[iwdescr_aaave].set_name("aaave");
    descriptor[iwdescr_admind].set_name("admind");
    descriptor[iwdescr_ad2mind].set_name("ad2mind");
    descriptor[iwdescr_adave].set_name("adave");
    descriptor[iwdescr_ddmind].set_name("ddmind");
    descriptor[iwdescr_dd2mind].set_name("dd2mind");
    descriptor[iwdescr_ddave].set_name("ddave");
  }

  if (descriptors_to_compute.complexity_descriptors)
  {
    descriptor[iwdescr_nspiro].set_name("nspiro");
    descriptor[iwdescr_nsfsdsys].set_name("nsfsdsys");
    descriptor[iwdescr_rnginsfs].set_name("rnginsfs");
    descriptor[iwdescr_lgstrfsy].set_name("lgstrfsy");
    descriptor[iwdescr_htrcsfsy].set_name("htrcsfsy");
    descriptor[iwdescr_mxhtsfsy].set_name("mxhtsfsy");
    descriptor[iwdescr_npfsdsys].set_name("npfsdsys");
    descriptor[iwdescr_rnginpfs].set_name("rnginpfs");
    descriptor[iwdescr_lgplnfsy].set_name("lgplnfsy");
    descriptor[iwdescr_htrcpfsy].set_name("htrcpfsy");
    descriptor[iwdescr_mxhtpfsy].set_name("mxhtpfsy");
  }

  if (descriptors_to_compute.symmetry_descriptors) {
    descriptor[iwdescr_symmatom].set_name("symmatom");
    descriptor[iwdescr_fsymmatom].set_name("fsymmatom");
    descriptor[iwdescr_lsepsymatom].set_name("lsepsymatom");
    descriptor[iwdescr_flsepsymatom].set_name("flsepsymatom");
    descriptor[iwdescr_maxsymmclass].set_name("maxsymmclass");
  }

  if (descriptors_to_compute.partial_symmetry_descriptors) {
    descriptor[iwdescr_maxpsymd].set_name("maxpsymd");
    descriptor[iwdescr_fmaxpsymd].set_name("fmaxpsymd");
    descriptor[iwdescr_maxpsymdmean].set_name("maxpsymdmean");
    descriptor[iwdescr_psymdnumzero].set_name("psymdnzero");
  }

  if (descriptors_to_compute.ramey_descriptors)
  {
    descriptor[iwdescr_obalance].set_name("obalance");
    descriptor[iwdescr_rmync].set_name("rmync");
    descriptor[iwdescr_rmynn].set_name("rmynn");
    descriptor[iwdescr_rmyno].set_name("rmyno");
    descriptor[iwdescr_rmynf].set_name("rmynf");
    descriptor[iwdescr_rmyns].set_name("rmyns");
    descriptor[iwdescr_rmyncl].set_name("rmyncl");
    descriptor[iwdescr_rmynbr].set_name("rmynbr");
    descriptor[iwdescr_rmyni].set_name("rmyni");
    descriptor[iwdescr_rmy_heavy_halogen].set_name("heavy_halogen");
  }
  if (descriptors_to_compute.compute_alogp) {
    descriptor[iwdescr_alogp].set_name("alogp");
  }
  if (descriptors_to_compute.compute_xlogp) {
    descriptor[iwdescr_xlogp].set_name("xlogp");
  }

  if (descriptors_to_compute.distance_matrix_descriptors) {
    descriptor[iwdescr_maxdarom].set_best_fingerprint(1);
  }
  if (descriptors_to_compute.psa) {
    descriptor[iwdescr_nvrtspsa].set_best_fingerprint(1);
  }
  descriptor[iwdescr_natoms].set_best_fingerprint(1);
  descriptor[iwdescr_frafus].set_best_fingerprint(1);
  if (descriptors_to_compute.distance_matrix_descriptors) {
    descriptor[iwdescr_maxdrng].set_best_fingerprint(1);
  }
  descriptor[iwdescr_aromc].set_best_fingerprint(1);
  descriptor[iwdescr_ro5_ohnh].set_best_fingerprint(1);
  descriptor[iwdescr_rmync].set_best_fingerprint(1);
  descriptor[iwdescr_fraromha].set_best_fingerprint(1);
  descriptor[iwdescr_ringatom].set_best_fingerprint(1);
  descriptor[iwdescr_rhaf].set_best_fingerprint(1);
  descriptor[iwdescr_rgmltbd].set_best_fingerprint(1);
  descriptor[iwdescr_ncon3].set_best_fingerprint(1);
  descriptor[iwdescr_brunsdon].set_best_fingerprint(1);
  descriptor[iwdescr_brunshbdsum].set_best_fingerprint(1);
  descriptor[iwdescr_avebbtwn].set_best_fingerprint(1);
  descriptor[iwdescr_amine].set_best_fingerprint(1);
  descriptor[iwdescr_npbarom].set_best_fingerprint(1);
  descriptor[iwdescr_hacts].set_best_fingerprint(1);
  descriptor[iwdescr_fchmltbd].set_best_fingerprint(1);
  descriptor[iwdescr_fnrnspc].set_best_fingerprint(1);
  descriptor[iwdescr_atmpiele].set_best_fingerprint(1);
  descriptor[iwdescr_ncon2].set_best_fingerprint(1);
  descriptor[iwdescr_funsatura].set_best_fingerprint(1);
  descriptor[iwdescr_csp3].set_best_fingerprint(1);
  descriptor[iwdescr_ch].set_best_fingerprint(1);
  descriptor[iwdescr_avsdlp].set_best_fingerprint(1);

  int rc = 1;

#ifdef TEST_FOR_NAMES

// All descriptors must have a name!

  for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
  {
    if (0 == descriptor[i].descriptor_name().length())
    {
      cerr << "Yipes, no name for descriptor " << i << '\n';
      rc = 0;
    }
  }
#endif

  return rc;
}

/*
  Some descriptors have default values - most of these are ring based
  descriptors that default to 0
*/

static void
initialise_descriptor_defaults()
{
  descriptor[iwdescr_rcj].set_default_value(0.0);
  descriptor[iwdescr_rchj].set_default_value(0.0);
  descriptor[iwdescr_amrcj].set_default_value(0.0);
  descriptor[iwdescr_alrcj].set_default_value(0.0);

  return;
}

/*
  When running tests on random smiles, we need to be able to save the
  previously computed results
*/

static Set_or_Unset<float> * saved_result = nullptr;

static void
MaybeFlush(IWString_and_File_Descriptor& output) {
  if (flush_after_each_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(4096);
  }
}

static int
write_fingerprint (Molecule & m,
                   const Descriptor * descriptor,
                   int n,
                   IWString_and_File_Descriptor & output)
{
  if (! work_as_tdt_filter)    // reading a smiles file
  {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  Sparse_Fingerprint_Creator sfc;

  int bstart = 0;

  for (int i = 0; i < n; i++)
  {
    if (! descriptor[i].produce_fingerprint())
      continue;

    descriptor[i].produce_fingerprint(bstart, sfc);

    bstart += descriptor[i].bit_replicates();
  }

  IWString tmp;
  sfc.daylight_ascii_form_with_counts_encoded(tag, tmp);
  output << tmp << '\n';

  if (! work_as_tdt_filter) {
    output << "|\n";
  }

  MaybeFlush(output);

  return 1;
}

static int
create_header(const Descriptor* descriptor,
              const IW_STL_Hash_Map_String& name_translation,
              const char output_separator,
              IWString & header)
{
  header.resize(500);

  int not_found_in_name_translation = 0;

  int rc = 0;
  for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
  {
    if (descriptor[i].active())
    {
      if (header.length())
        header << output_separator;

      if (name_translation.empty()) {
        header << descriptor[i].descriptor_name();
        ++rc;
        continue;
      }

      auto iter = name_translation.find(descriptor[i].descriptor_name());
      if (iter == name_translation.end()) {
        header << descriptor[i].descriptor_name();
        ++not_found_in_name_translation;
        if (verbose > 1) {
          cerr << "Feature " << descriptor[i].descriptor_name() << " not in cross reference\n";
        }
      } else {
        header << iter->second;
      }
      rc++;
    }
  }

  if (not_found_in_name_translation) {
    cerr << "Warning " << not_found_in_name_translation << " features not found in name translation table\n";
  }

  return rc;
}


static int
WriteHeader(const Descriptor* descriptor,
             const IW_STL_Hash_Map_String& name_translation,
             const char output_separator,
             IWString_and_File_Descriptor & output) {
  IWString header;

  int tokens = create_header(descriptor, name_translation, output_separator, header);

  if (include_smiles_as_descriptor) {
    output << "smiles" << output_separator;
  }

  output << header << '\n';

  if (verbose) {
    cerr << "Output will contain " << tokens << " descriptors\n";
  }

  return 1;
}

static int
do_filtering(Molecule & m,
             IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < number_filters; i++)
  {
    if (filter[i].satisfied(descriptor))
      continue;

    rejected_by_filters++;
    return 1;
  }

  output << smiles_as_input << ' ' << m.name() << '\n';

  return 1;
}

// If `v` is close to being a whole number, write as an integer to `output`.
int
OutputAsWholeNumber(float v,
                    IWString& output) {
  if (fabs(static_cast<int>(v) - v) > 1.0e-05) {
    return 0;
  }

  iwdigits.append_number(output, static_cast<int>(v));

  return 1;
}

static int
write_the_output(Molecule & m,
                 const char output_separator,
                 IWString_and_File_Descriptor & output)
{
  if (ntest) {    // no output in test mode
    return 1;
  }

  if (number_filters) {
    return do_filtering(m, output);
  }

  if (tag.length()) {
    return write_fingerprint(m, descriptor, NUMBER_DESCRIPTORS, output);
  }

  if (read_descriptor_file_pipeline && write_descriptor_file_pipeline) {
    output << m.smiles() << ' ';
    output << m.name();  // includes all previously calculated descriptors.
  } else if (read_descriptor_file_pipeline) {
    output << m.name();  // includes all previously calculated descriptors.
  } else if (write_descriptor_file_pipeline) {
    output << m.smiles() << ' ';
    append_first_token_of_name(m.name(), output);
  } else {
    append_first_token_of_name(m.name(), output);
  }

  if (include_smiles_as_descriptor) {
    output << output_separator << smiles_as_input;
  }

  for (int i = 0; i < NUMBER_DESCRIPTORS; i++) {
    if (! descriptor[i].active()) {
      continue;
    }

    output << output_separator;
    cerr << "writing " << i << std::endl;

    float v;
    if (descriptor[i].value(v)) {
      if (OutputAsWholeNumber(v, output)) {
      } else if (v > 100.0f)
        output.append_number(v, 5);
      else
        output << v;
    } else {
      output << undefined_value;
    }
  }

  output << '\n';

  MaybeFlush(output);

  if (verbose)
  {
    for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
    {
      if (descriptor[i].active())
        descriptor[i].update_statistics();
    }
  }

  return output.good();
}

//#define DEBUG_RING_SUBSTITUTION_STUFF
#ifdef DEBUG_RING_SUBSTITUTION_STUFF
static int
print_array (const extending_resizable_array<int> & z,
             const char * c,
             IWString & output)
{
  if (z.empty())
    return output.good();

  output << c << '\n';

  for (int i = 0; i < z.number_elements(); i++)
  {
    if (z[i] > 0)
      output << i << " bonds, count " << z[i] << '\n';
  }

  return output.good();
}
#endif

#ifdef MCGOWAN
static float
McGowanVolume(Molecule& m) {
  float rc = 0.0;

  // Arbitrarily stop copying data once we have the organic set covered.
  static float mcgowan[57] = {
    0.0,    // 0
    8.71,   // 1
    0.0,        // 2
    22.23,      // 3
    20.27,      // 4
    18.31,      // 5
    16.35,      // 6
    14.39,      // 7
    12.43,      // 8
    10.47,      // 9
    8.51,       // 10
    32.71,      // 11
    30.75,      // 12
    28.79,      // 13
    26.83,      // 14
    24.87,      // 15
    22.91,      // 16
    20.95,      // 17
    18.99,      // 18
    51.89,      // 19
    50.28,      // 20
    48.68,      // 21
    47.07,      // 22
    45.47,      // 23
    43.86,      // 24
    42.26,      // 25
    40.65,      // 26
    39.05,      // 27
    37.44,      // 28
    35.84,      // 29
    34.23,      // 30
    32.63,      // 31
    31.02,      // 32
    29.42,      // 33
    27.81,      // 34
    26.21,      // 35
    24.60,      // 36
    60.22,      // 37
    58.61,      // 38
    57.01,      // 39
    55.40,      // 40
    53.80,      // 41
    52.19,      // 42
    50.59,      // 43
    48.98,      // 44
    47.38,      // 45
    45.77,      // 46
    44.17,      // 47
    42.56,      // 48
    40.96,      // 49
    39.35,      // 50
    37.75,      // 51
    36.14,      // 52
    34.54,      // 53
    32.93,      // 54
    77.25,      // 55
    76.00       // 56
  };

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    rc += mcgowan[m.atomic_number(i)];
    rc += m.implicit_hydrogens(i);
  }

  return rc;
}
#endif

/*
  In a ring path of size N, we have a difference of D bonds. If this is larger than
  half the ring size (N2), then we need to go the other way around the ring to get
  the shortest distance
*/

static void
compute_ring_substitution_descriptors (extending_resizable_array<int> & diffs,
                                       int n,
                                       int n2,
                                       int d)
{
  assert (d > 0);

//cerr << "N = " << n << " N2 = " << n2 << " diff " << d << '\n';

  if (d > n2)
    d = n - d;

//cerr << "d = " << d << '\n';

  diffs[d]++;

  return;
}

static int
compute_ring_substitution_descriptors(Molecule & m,
                                      const int * ncon,
                                      const Set_of_Atoms & par,
                                      extending_resizable_array<int> & diffs)
{
//#define DEBUG_COMPUTE_RING_SUBSTITUTION_DESCRIPTORS
#ifdef DEBUG_COMPUTE_RING_SUBSTITUTION_DESCRIPTORS
  cerr << "Processing " << par << "\n";
#endif

  int n = par.number_elements();

  resizable_array<int> offset_of_branch;

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = par[i];

#ifdef DEBUG_COMPUTE_RING_SUBSTITUTION_DESCRIPTORS
    cerr << "Atom " << j << " is part of the ring edge\n";
#endif

    if (2 == ncon[j])
      continue;

    if (m.nrings(j) > 1)
      continue;

    offset_of_branch.add(i);
  }

  int noffset = offset_of_branch.number_elements();

  if (noffset < 2)     // ring with 0 or just 1 attachment
    return 1;

#ifdef DEBUG_COMPUTE_RING_SUBSTITUTION_DESCRIPTORS
  cerr << "Found " << noffset << " substitution points from the ring\n";
  if (2 == noffset)
    cerr << "Atoms " << par[offset_of_branch[1]] << " and " << par[offset_of_branch[0]] << '\n';
#endif

  int n2 = n / 2;     // half the ring size. Works with odd sized rings too

  if (2 == noffset)    // a very common case
  {
    compute_ring_substitution_descriptors(diffs, n, n2, offset_of_branch[1] - offset_of_branch[0]);

    return 1;
  }

  for (int i = 0; i < noffset; i++)
  {
    int di = offset_of_branch[i];

    for (int j = i + 1; j < noffset; j++)
    {
      int dj = offset_of_branch[j];

      compute_ring_substitution_descriptors(diffs, n, n2, dj - di);
    }
  }

  return 1;
}

int
BondedToAnotherRing(Molecule& m,
                    const int* ring_membership,
                    atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    if (b->nrings()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (ring_membership[o]) {
      return 1;
    }
  }

  return 0;
}

// Compute the number of ortho rings attached to `ring`.
static int
OrthoRingCount(Molecule& m,
               const int* ncon,
               const int* ring_membership,
               const Ring& ring) {
  int rc = 0;
  const int ring_size = ring.number_elements();

  for (int i = 0; i < ring_size; ++i) {
    atom_number_t a1 = ring[i];
    if (ncon[a1] == 2) {
      continue;
    }
    atom_number_t a2 = ring[(i + 1) % ring_size];
    if (ncon[a2] == 2) {
      continue;
    }
    if (! BondedToAnotherRing(m, ring_membership, a1)) {
      continue;
    }

    if (BondedToAnotherRing(m, ring_membership, a2)) {
      ++rc;
    }
  }

  return rc;
}

/*
  Measures of how crowded the rings are
*/

static int
compute_ring_substitution_ratios(Molecule & m,
                                 const int * ncon,
                                 const Atom * const * atom,
                                 const int * ring_membership)
{
  int ring_atoms = 0;
  int substituted_outside_ring = 0;
  int simple_substituted_outside_ring = 0;

  if (0 == m.nrings())
  {
    descriptor[iwdescr_frsub].set(static_cast<float>(0.0));
    descriptor[iwdescr_frssub].set(static_cast<float>(0.0));
    descriptor[iwdescr_arorthoring].set(static_cast<float>(0.0));
    descriptor[iwdescr_alorthoring].set(static_cast<float>(0.0));

    return 1;
  }

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == ring_membership[i])
      continue;

    ring_atoms++;

    if (2 == ncon[i])    // definitely not substituted outside the ring
      continue;

    const Atom * a = atom[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = a->other(i, j);

      if (ring_membership[k])
        continue;

      substituted_outside_ring++;
      if (1 == ncon[k])
        simple_substituted_outside_ring++;
    }
  }

  descriptor[iwdescr_frsub].set(static_cast<float>(substituted_outside_ring) / static_cast<float>(ring_atoms));
  descriptor[iwdescr_frssub].set(static_cast<float>(simple_substituted_outside_ring) / static_cast<float>(ring_atoms));

  m.compute_aromaticity_if_needed();

  int aromatic_ortho_ring_count = 0;
  int aliphatic_ortho_ring_count = 0;
  for (const Ring* r : m.sssr_rings()) {
    if (r->is_aromatic()) {
      aromatic_ortho_ring_count += OrthoRingCount(m, ncon, ring_membership, *r);
    } else {
      aliphatic_ortho_ring_count += OrthoRingCount(m, ncon, ring_membership, *r);
    }
  }

  descriptor[iwdescr_arorthoring].set(aromatic_ortho_ring_count);
  descriptor[iwdescr_alorthoring].set(aliphatic_ortho_ring_count);

  return 1;
}

static int
bonds_to_nearest_ring(const Molecule & m,
                      atom_number_t previous_atom,
                      atom_number_t current_atom,
                      const Atom * const * atom,
                      const int * fsid)
{
  const Atom * a = atom[current_atom];

  int acon = a->ncon();

  int matoms = m.natoms();

  int rc = matoms;    // we look for the shortest path

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other(current_atom, i);

    if (previous_atom == j)
      continue;

    if (fsid[j] > 0)
      return 1;

    int tmp = bonds_to_nearest_ring(m, current_atom, j, atom, fsid);

    if (tmp <= 0)
      continue;

    tmp++;

    if (tmp < rc)
      rc = tmp;
  }

  if (matoms == rc)    // didn't find anything
    return -1;

  return rc;
}

static int
compute_between_ring_descriptors(Molecule & m,
                                 const int * ncon,
                                 const Atom * const * atom)
{
  int nr = m.nrings();

  if (nr < 2)    // we deal with between ring information
  {
    descriptor[iwdescr_bbr1].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr2].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr3].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr4].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr5].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr6].set(static_cast<float>(0.0));

    return 1;
  }

  int matoms = m.natoms();

  int * fsid = new int[matoms]; std::unique_ptr<int[]> free_fsid(fsid);

  if (1 == m.label_atoms_by_ring_system(fsid))   // more than 1 ring, but 1 ring system
  {
    descriptor[iwdescr_bbr1].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr2].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr3].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr4].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr5].set(static_cast<float>(0.0));
    descriptor[iwdescr_bbr6].set(static_cast<float>(0.0));

    return 1;
  }

  extending_resizable_array<int> bonds_between_rings;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == fsid[i])    // not in a ring
      continue;

    if (2 == ncon[i])
      continue;

    const Atom * ai = atom[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      const Bond * b = ai->item(j);

      if (b->nrings())
        continue;

      atom_number_t k = b->other(i);

      if (fsid[k])
        bonds_between_rings[1]++;
      else
      {
        int tmp = bonds_to_nearest_ring(m, i, k, atom, fsid);
        if (tmp > 0)
          bonds_between_rings[tmp + 1]++;
      }
    }
  }

  descriptor[iwdescr_bbr1].set(static_cast<float>(bonds_between_rings[1] / 2));
  descriptor[iwdescr_bbr2].set(static_cast<float>(bonds_between_rings[2] / 2));
  descriptor[iwdescr_bbr3].set(static_cast<float>(bonds_between_rings[3] / 2));
  descriptor[iwdescr_bbr4].set(static_cast<float>(bonds_between_rings[4] / 2));
  descriptor[iwdescr_bbr5].set(static_cast<float>(bonds_between_rings[5] / 2));
  descriptor[iwdescr_bbr6].set(static_cast<float>(bonds_between_rings[6] / 2));

  return 1;
}

static int
count_ring_bonds(const Atom & a)
{
  int rc = 0;

  for (const Bond*b : a) {
    if (b->nrings() > 0) {
      ++rc;
    }
  }

  return rc;
}

/*
  Look for ring atoms that are adjacent to a fusion and are connected
  outside the ring
*/

static int
compute_adjacent_ring_fusion_descriptors(Molecule & m,
                                         const int * ncon,
                                         const int * ring_membership)
{
  int single_bonds_found = 0;
  int double_bonds_found = 0;

  (void) m.ring_membership();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (ring_membership[i] < 2)
      continue;

    if (3 != ncon[i])
      continue;

    const Atom * ai = m.atomi(i);

//  cerr << "Atom " << i << " 3 connections, rbc " << count_ring_bonds(*ai) << '\n';

    // very important, otherwise C1CC2C1C(C2CCCNC)C1=CC=CC=C1 fails
    if (3 != count_ring_bonds(*ai)) {
      continue;
    }

    for (const Bond* b : *ai) {
//    cerr << "Checking bond to " << b->other(i) << " nrings " << b->nrings() << '\n';

      atom_number_t k = b->other(i);

      const Atom * ak = m.atomi(k);

      if (2 != count_ring_bonds(*ak))
        continue;

      // No exocyclic bonds here.
      if (ak->ncon() == 2) {
        continue;
      }
      
      for (const Bond* b: *ak) {
//      cerr << "From " << k << " to " << b->other(k) << " bond in " << b->nrings() << " rings\n";
        if (b->nrings())
          ;
        else if (b->is_single_bond())
          single_bonds_found++;
        else
          double_bonds_found++;
      }
    }
  }

  descriptor[iwdescr_sboradjf].set(static_cast<float>(single_bonds_found));
  descriptor[iwdescr_dboradjf].set(static_cast<float>(double_bonds_found));

  return 1;
}

static int
compute_xlogp(Molecule& m) {
  std::optional<double> x = xlogp::XLogP(m);
  if (x) {
    descriptor[iwdescr_xlogp] = *x;
    return 1;
  }

  return 0;
}

static int
compute_alogp(Molecule& m) {
  std::optional<double> x = alogp_engine.LogP(m);
  if (x) {
    descriptor[iwdescr_alogp] = *x;
    return 1;
  }

  return 0;
}

static int
compute_ring_substitution_descriptors(Molecule & m,
                                      const int * ncon)
{
  int nr = m.nrings();

  if (0 == nr)
    return 1;

  int matoms = m.natoms();

  extending_resizable_array<int> bonds_between_ring_substitutions_aromatic_ring;
  extending_resizable_array<int> bonds_between_ring_substitutions_aliphatic_ring;

// Do all the unfused rings first

  int fused_rings_present = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_fused())
    {
      fused_rings_present++;
      continue;
    }

    if (ri->is_aromatic())
      compute_ring_substitution_descriptors(m, ncon, *ri, bonds_between_ring_substitutions_aromatic_ring);
    else
      compute_ring_substitution_descriptors(m, ncon, *ri, bonds_between_ring_substitutions_aliphatic_ring);
  }

  extending_resizable_array<int> bonds_between_ring_substitutions_ring_system;

// Now any ring systems

  if (fused_rings_present)
  {
    int * fsid = new_int(matoms, -1); std::unique_ptr<int[]> free_fsid(fsid);
    int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

    for (int i = 0; i < nr; i++)
    {
      if (ring_already_done[i])
        continue;

      const Ring * ri = m.ringi(i);

      int f = ri->fused_system_identifier();

      if (f < 0)    // not fused
        continue;

      int strongly_fused_ring_found;

      if (ri->largest_number_of_bonds_shared_with_another_ring() > 1)
        strongly_fused_ring_found = 1;
      else
        strongly_fused_ring_found = 0;

      ri->set_vector(fsid, f + 1);

      for (int j = i + 1; j < nr; j++)
      {
        if (ring_already_done[j])
          continue;

        const Ring * rj = m.ringi(j);
        if (f != rj->fused_system_identifier())
          continue;

        ring_already_done[j] = 1;

//      cerr << "Ring " << j << " has " << rj->largest_number_of_bonds_shared_with_another_ring() << " bonds shared\n";
        if (strongly_fused_ring_found)
          ;
        else if (rj->largest_number_of_bonds_shared_with_another_ring() > 1)
          strongly_fused_ring_found = 1;
        else
          rj->set_vector(fsid, f + 1);
      }

//    cerr << "Any strongly fused rings? " << strongly_fused_ring_found << '\n';
      if (strongly_fused_ring_found)   // can't do these yet
        continue;

      Set_of_Atoms s;

//    cerr << "Processing atoms with flag " << (f + 1) << ", i = " << i << '\n';
      if (! path_around_edge_of_ring_system(m, fsid, f + 1, s))
        continue;

//    cerr << "Atoms around ring " << s << '\n';
//    cerr << "S contains " << s.number_elements() << " atoms\n";
      compute_ring_substitution_descriptors(m, ncon, s, bonds_between_ring_substitutions_ring_system);
    }

#ifdef ECHO_FSID
    for (int i = 0; i < matoms; i++)
    {
      cerr << "Atom " << i << " fsid " << fsid[i] << '\n';
    }
#endif
  }

#ifdef DEBUG_RING_SUBSTITUTION_STUFF
  cerr << m.name() << '\n';
  print_array(bonds_between_ring_substitutions_aromatic_ring, "aromatic ring", cerr);
  print_array(bonds_between_ring_substitutions_aliphatic_ring, "aliphatic ring", cerr);
  print_array(bonds_between_ring_substitutions_ring_system, "ring system", cerr);
#endif

  descriptor[iwdescr_rsarom1].set(bonds_between_ring_substitutions_aromatic_ring[1]);
  descriptor[iwdescr_rsarom2].set(bonds_between_ring_substitutions_aromatic_ring[2]);
  descriptor[iwdescr_rsarom3].set(bonds_between_ring_substitutions_aromatic_ring[3]);

  descriptor[iwdescr_rsaliph1].set(bonds_between_ring_substitutions_aliphatic_ring[1]);
  descriptor[iwdescr_rsaliph2].set(bonds_between_ring_substitutions_aliphatic_ring[2]);
  descriptor[iwdescr_rsaliph3].set(bonds_between_ring_substitutions_aliphatic_ring[3]);

  if (bonds_between_ring_substitutions_aliphatic_ring.number_elements() < 4)
    descriptor[iwdescr_rsaliph4].set(0);
  else
  {
    int n4 = 0;
    for (int i = 4; i < bonds_between_ring_substitutions_aliphatic_ring.number_elements(); i++)
    {
      n4 += bonds_between_ring_substitutions_aliphatic_ring[i];
    }
    descriptor[iwdescr_rsaliph4].set(n4);
  }
  
  int istop = bonds_between_ring_substitutions_ring_system.number_elements();
  if (istop > 8)
    istop = 8;

  descriptor[iwdescr_rssys1].set(bonds_between_ring_substitutions_ring_system[1]);
  descriptor[iwdescr_rssys2].set(bonds_between_ring_substitutions_ring_system[2]);
  descriptor[iwdescr_rssys3].set(bonds_between_ring_substitutions_ring_system[3]);
  descriptor[iwdescr_rssys4].set(bonds_between_ring_substitutions_ring_system[4]);
  descriptor[iwdescr_rssys5].set(bonds_between_ring_substitutions_ring_system[5]);
  descriptor[iwdescr_rssys6].set(bonds_between_ring_substitutions_ring_system[6]);
  descriptor[iwdescr_rssys7].set(bonds_between_ring_substitutions_ring_system[7]);
  descriptor[iwdescr_rssys8].set(bonds_between_ring_substitutions_ring_system[8]);
  descriptor[iwdescr_rssys9].set(bonds_between_ring_substitutions_ring_system[9]);

  return 1;
}

static void
compute_ramey_descriptors(Molecule & m,
                          const atomic_number_t * z)
{
  if (m.contains_non_periodic_table_elements())
  {
    descriptor[iwdescr_obalance].set(0.0f);
    return;
  }

  int matoms = m.natoms();

  int nc = 0;
  int nn = 0;
  int no = 0;
  int nf = 0;
  int ns = 0;
  int ncl = 0;
  int nbr = 0;
  int ni = 0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t zi = z[i];

    if (6 == zi)
      nc++;
    else if (7 == zi)
      nn++;
    else if (8 == zi)
      no++;
    else if (9 == zi)
      nf++;
    else if (16 == zi)
      ns++;
    else if (17 == zi)
      ncl++;
    else if (35 == zi)
      nbr++;
    else if (53 == zi)
      ni++;
  }

  descriptor[iwdescr_rmync].set(static_cast<float>(nc));
  descriptor[iwdescr_rmynn].set(static_cast<float>(nn));
  descriptor[iwdescr_rmyno].set(static_cast<float>(no));
  descriptor[iwdescr_rmynf].set(static_cast<float>(nf));
  descriptor[iwdescr_rmyns].set(static_cast<float>(ns));
  descriptor[iwdescr_rmyncl].set(static_cast<float>(ncl));
  descriptor[iwdescr_rmynbr].set(static_cast<float>(nbr));
  descriptor[iwdescr_rmyni].set(static_cast<float>(ni));
  descriptor[iwdescr_rmy_heavy_halogen].set(static_cast<float>(ncl + nbr + ni));

  const int nh = m.implicit_hydrogens();

  const auto oxygen_balance = - 1600.0 * (2 * nc + static_cast<double>(nh) / 2.0 - no) / m.molecular_weight_ignore_isotopes();

  descriptor[iwdescr_obalance].set(static_cast<float>(oxygen_balance));

  return;
}

/*
  Is an atom part of a spiro ring fusion.

  note that we exclude the case of a spiro fusion within a ring system:

    -----------
    |         |
     \       /
      \     /
       spiro

  Too wierd to worry about...

  We look for atoms that have 4 connections and are in 2 rings. If either of
  the rings that contain the atom are non-fused, then we immediately have a 
  spiro fusion
*/

static int
is_spiro_fused(Molecule & m,
               atom_number_t a)
{
  int nr = m.nrings();
  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (! ri->contains(a))
      continue;

    if (! ri->is_fused())    // not fused, must be spiro
      return 1;

//  RI is fused. Look for the other ring that also contains A

    for (int j = i + 1; j < nr; j++)
    {
      const Ring * rj = m.ringi(j);

      if (! rj->contains(a))
        continue;

      if (! rj->is_fused())     // got it
        return 1;

//    Both rings are fused. If they are in different ring systems, we have it

      if (ri->fused_system_identifier() != rj->fused_system_identifier())
        return 1;

      return 0;       // we have checked 2 rings and no spiro fusion
    }
  }

  return 0;
}

static int
do_compute_spiro_fusions(Molecule & m,
                         const int * ncon,
                         const int * ring_membership)
{
  int nr = m.nrings();

  if (! compute_spiro_fusions || nr < 2) {
    descriptor[iwdescr_nspiro].set(0.0);
    return 1;
  }

  int matoms = m.natoms();

  int spiro_fusions = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (4 != ncon[i])
      continue;

    if (ring_membership[i] < 2)
      continue;

    if (is_spiro_fused(m, i))
      spiro_fusions++;
  }

  descriptor[iwdescr_nspiro].set(static_cast<float>(spiro_fusions));

  return 1;
}

static int
do_compute_chirality_descriptors(Molecule & m,
                                 const atomic_number_t * z,
                                 const int * ncon)
{
  int matoms = m.natoms();

  int chiral_centres = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (6 == z[i])
    {
      if (4 == ncon[i])
        ;
      else if (3 == ncon[i] && 1 == m.implicit_hydrogens(i))
        ;
      else
        continue;
    }
    else if (16 == z[i])
    {
      if (3 == ncon[i] && 4 == m.nbonds(i))
        ;
      else
        continue;
    }
    else
      continue;

    if (nullptr != m.chiral_centre_at_atom(i))
    {
      chiral_centres++;
      continue;
    }

    if (! descriptors_to_compute.perform_expensive_chirality_perception)
      ;
    else if (m.is_ring_atom(i))    // C1C2C3C1N1C2C31 has 5 chiral centers otherwise.
      ;
    else if (is_actually_chiral(m, i))  {   // Too expensive to compute by default.
      chiral_centres++;
    }
  }

  descriptor[iwdescr_nchiral].set(static_cast<float>(chiral_centres));

  return 1;
}

// Partial symmetry descriptors.
void
NearSymmetricDescriptors(Molecule& m) {
  partial_symmetry::PartialSymmetry psim(m);
  const int * symmetry = psim.SymmetricAtRadius();
#ifdef DEBUG_PARTIAL_SYMMETRY
  Molecule tmp(m);
  write_isotopically_labelled_smiles(tmp, false, cerr);
  cerr << '\n';
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << " " << tmp.smarts_equivalent_for_atom(i) << " value " << symmetry[i] << '\n';
  }
#endif

  Accumulator_Int<int> acc;
  const int matoms = m.natoms();
  acc.extra(symmetry, matoms);
  descriptor[iwdescr_maxpsymd] = acc.maxval();
  descriptor[iwdescr_fmaxpsymd] = Fraction<float>(acc.maxval(), matoms);
  descriptor[iwdescr_maxpsymdmean] = acc.average();

  const int nzero = std::count(symmetry, symmetry + matoms, 0);
  descriptor[iwdescr_psymdnumzero] = nzero;

  return;
}

// Given separated atoms `a1` and `a2` compute the mean distance of
// all other atoms to these.
void
compute_atom_distribution_along_longest_path(const Molecule & m,
                        const int * dm,
                        atom_number_t a1,
                        atom_number_t a2) {

  Accumulator_Int<int> to_a1, to_a2;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == a1 || i == a2) {
      continue;
    }
    int d = dm[a1 * matoms + i];
    to_a1.extra(d);
    d = dm[a2 * matoms + i];
    to_a2.extra(d);
  }

//cerr << "Between atoms " << a1 << ' ' << to_a1 << '\n';
//cerr << "Between atoms " << a2 << ' ' << to_a2 << '\n';

  if (to_a1.n() == 0) {
    descriptor[iwdescr_mdallp]  = 0.0f;
    descriptor[iwdescr_fmdallp]  = 0.0f;
    descriptor[iwdescr_fdiffallp]  = 0.0f;
    return;
  }

  const float mean1 = to_a1.average();
  const float mean2 = to_a2.average();
  float minval = std::min(mean1, mean2);
  descriptor[iwdescr_mdallp] = minval;
  descriptor[iwdescr_fmdallp] = minval / static_cast<float>(dm[a1 * matoms + a2]);
  const float diff = abs(mean1 - mean2);
  descriptor[iwdescr_fdiffallp] = diff / static_cast<float>(dm[a1 * matoms + a2]);

  return;
}

void
AppendShell(const Molecule& m,
            atom_number_t zatom,
            int radius,
            const int * dm,
            Set_of_Atoms& shell) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == zatom) {
      continue;
    }
    if (dm[zatom * matoms + i] == radius) {
      shell << i;
    }
  }
}

static int
do_compute_mean_shell_occupancies(const Molecule& m,
                                  const int * dm) {
  const int matoms = m.natoms();
  const int max_radius = 3;

  std::unique_ptr<Accumulator_Int<int>[]> acc(new Accumulator_Int<int>[max_radius + 1]);

  Set_of_Atoms shell;  // Scope here for efficiency.
  for (int i = 0; i < matoms; ++i) {
    shell.resize_keep_storage(0);
    for (int r = 1; r <= max_radius; ++r) {
      AppendShell(m, i, r, dm, shell);
      acc[r].extra(shell.number_elements());
    }
  }

  descriptor[iwdescr_aveshell1].set(acc[1].average());
  descriptor[iwdescr_aveshell2].set(acc[2].average());
  descriptor[iwdescr_aveshell3].set(acc[3].average());
  descriptor[iwdescr_maxshell3].set(acc[3].maxval());

  return 1;
}

#define CURRENT_SHELL 30
#define NEXT_SHELL 99
#define SHELL_COMPLETED -5

/*
  We determine the atoms that have the smallest total distance to
  other atoms in the molecule.
  We count them.

  Then we do a couple of shell expansions to count the number of atoms
  encountered.
*/

static int
do_compute_descriptors_related_to_centroid_atom(Molecule & m,
                                                const int * dm,
                                                const int longest_path,
                                                const atomic_number_t * z,
                                                const int * ncon,
                                                int * in_shell)
{
  const auto matoms = m.natoms();

  Set_of_Atoms atoms_at_min;
  int mind = longest_path * longest_path;

  for (auto i = 0; i < matoms; ++i)
  {
    int tot_dist = std::accumulate(dm + i * matoms, dm + i * matoms + matoms, 0);

    if (tot_dist < mind) {
      atoms_at_min.resize_keep_storage(0);
      atoms_at_min.add(i);
      mind = tot_dist;
    }
    else if (tot_dist == mind) {
      atoms_at_min.add(i);
    }
  }

  std::fill_n(in_shell, matoms, 0);

  int atoms_in_next_shell = 0;

  for (auto i = 0; i < atoms_at_min.number_elements(); ++i)
  {
    const auto j = atoms_at_min[i];
    // m.set_isotope(j, 1);
    in_shell[j] = CURRENT_SHELL;

    const auto a = m.atomi(j);
    for (auto k = 0; k < ncon[j]; ++k)
    {
      const atom_number_t s = a->other(j, k);
      if (atoms_at_min.contains(s))
        continue;

      atoms_in_next_shell++;
      in_shell[s] = NEXT_SHELL;
    }
  }

  descriptor[iwdescr_cntrdgncy].set(atoms_at_min.number_elements());
  descriptor[iwdescr_cntrdshell1].set(atoms_in_next_shell);

  const int max_centroid_bond_radius = 2;    // change if ever needed

  for (auto r = 1; r <= max_centroid_bond_radius; ++r)
  {
    for (auto i = 0; i < matoms; ++i)
    {
      if (0 == in_shell[i])
        continue;

      if (CURRENT_SHELL == in_shell[i])
        in_shell[i] = SHELL_COMPLETED;
      else if (NEXT_SHELL == in_shell[i])
        in_shell[i] = CURRENT_SHELL;
    }

    atoms_in_next_shell = 0;

    for (auto i = 0; i < matoms; ++i)
    {
      if (CURRENT_SHELL != in_shell[i])
        continue;

      const auto a = m.atomi(i);

      //cerr << "Processing atom " << i << " ncon " << a->ncon() << " array " << ncon[i] << '\n';

      for (auto j = 0; j < ncon[i]; ++j)
      {
//      cerr << "from atom " << i << " ncon " << ncon[i] << " connection " << j << " type " << m.smarts_equivalent_for_atom(i) << '\n';
        atom_number_t k = a->other(i, j);

        if (0 != in_shell[k])
          continue;

        atoms_in_next_shell++;
        in_shell[k] = NEXT_SHELL;
      }
    }
    descriptor[iwdescr_cntrdshell1 + r].set(atoms_in_next_shell);
  }

//cerr << m.smiles() << ' ' << m.name() << '\n';

  return 1;
}

static void
compute_shortest_distance_from_longest_path(Molecule & m,
                                            const int * dm,
                                            const int * in_path,
                                            int * shortest_distance_from_longest_path)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)    // atom I is outside the longest path
  {
    if (in_path[i])
      continue;

    for (int j = 0; j < matoms; j++)   // atom J is in the longest path
    {
      if (! in_path[j])
        continue;

      int b = dm[i * matoms + j];

      if (b < shortest_distance_from_longest_path[i])
        shortest_distance_from_longest_path[i] = b;
    }
  }

  return;
}

static void
do_compute_distances_from_longest_path_descriptors(Molecule & m,
                                                   const int * dm,
                                                   atom_number_t astart,
                                                   atom_number_t astop,
                                                   int * in_path,
                                                   int * shortest_distance_from_longest_path)
{
  in_path[astart] = 1;

  const int matoms = m.natoms();

  const Atom * a = m.atomi(astart);

//const int current_distance = m.bonds_between(astart, astop);
  const int current_distance = dm[astart * matoms + astop];

  for (const Bond* b : *a) {
    const atom_number_t j = b->other(astart);

    if (j == astop)
      compute_shortest_distance_from_longest_path(m, dm, in_path, shortest_distance_from_longest_path);
    else if (in_path[j])    // no doubling back
      ;
    else if (current_distance - 1 == dm[j * matoms + astop])
      do_compute_distances_from_longest_path_descriptors(m, dm, j, astop, in_path, shortest_distance_from_longest_path);
  }

  in_path[astart] = 0;

  return;
}

static int
do_compute_distances_from_longest_path_descriptors(Molecule & m,
                                                   const int * dm,
                                                   atom_number_t a1,
                                                   atom_number_t a2,
                                                   int * in_path,
                                                   const int longest_path)
{
  compute_atom_distribution_along_longest_path(m, dm, a1, a2);

  int matoms = m.natoms();
  
  set_vector(in_path, matoms, 0);
  in_path[a2] = 1;

  int * shortest_distance_from_longest_path = new_int(matoms, matoms); std::unique_ptr<int[]> free_sdlp(shortest_distance_from_longest_path);

  shortest_distance_from_longest_path[a2] = 0;
  do_compute_distances_from_longest_path_descriptors(m, dm, a1, a2, in_path, shortest_distance_from_longest_path);

  Accumulator_Int<int> acc;
  for (int i = 0; i < matoms; i++)
  {
    if (matoms != shortest_distance_from_longest_path[i])
      acc.extra(shortest_distance_from_longest_path[i]);
  }

  if (0 == acc.n())    // linear chain type molecule
  {
    descriptor[iwdescr_mxsdlp].set(0.0);
    descriptor[iwdescr_avsdlp].set(0.0);
    descriptor[iwdescr_mxsdlprl].set(0.0);
  }
  else
  {
    descriptor[iwdescr_mxsdlp].set(acc.maxval());
    descriptor[iwdescr_avsdlp].set(acc.average());
    descriptor[iwdescr_mxsdlprl].set(static_cast<float>(acc.maxval()) / static_cast<float>(longest_path));
  }

  return 1;
}

/*
  This is complicated by the fact that there can be many paths of
  maximum length. To break ties, we use the molecular canonicalisation
  function
*/

static int
do_compute_distances_from_longest_path_descriptors(Molecule & m,
                                                   const int * dm,
                                                   const int longest_path,
                                                   int * in_path)
{
  assert (nullptr != dm);

  const int matoms = m.natoms();

  Set_of_Atoms a1, a2;

  for (int i = 0; i < matoms; i++)
  {
    for (int j = i + 1; j < matoms; j++)
    {
      if (longest_path != dm[i * matoms + j])
        continue;

      a1.add(i);
      a2.add(j);
    }
  }

  // Single atom molecule, or all atoms disconnected.
  if (a1.empty()) {
    return 1;
  }

  int n = a1.number_elements();

  if (1 == n)
    return do_compute_distances_from_longest_path_descriptors(m, dm, a1[0], a2[0], in_path, longest_path);

  int highest_score = -1;
  int pair_with_highest_score = -1;

  for (int i = 0; i < n; i++)
  {
    int c1 = m.canonical_rank(a1[i]);
    int c2 = m.canonical_rank(a2[i]);

    int canonical_sum;

    if (c1 > c2)
      canonical_sum = c1 * (matoms + matoms) + c2;
    else
      canonical_sum = c2 * (matoms + matoms) + c1;

    if (canonical_sum > highest_score)
    {
      highest_score = canonical_sum;
      pair_with_highest_score = i;
    }
  }

  return do_compute_distances_from_longest_path_descriptors(m, dm, a1[pair_with_highest_score], a2[pair_with_highest_score], in_path, longest_path);
}

// We store an atom status into an int.
static constexpr uint32_t kIsRing = 1;
static constexpr uint32_t kIsAromatic = 2;
static constexpr uint32_t kIsHeteroatom = 4;
static constexpr uint32_t kIsONS = 8;

static std::unique_ptr<uint32_t[]>
AssignDistanceMatrixAtomTypes(Molecule& m,
                              const atomic_number_t* z) {
  const int matoms = m.natoms();

  // To be returned
  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);

  for (int i = 0; i < matoms; ++i){
    atype[i] = 0;
    if (m.ring_bond_count(i) > 0) {
      atype[i] |= kIsRing;
      if (m.is_aromatic(i)) {
        atype[i] |= kIsAromatic;
      }
    }
    if (z[i] == 6) {
      continue;
    }
    atype[i] |= kIsHeteroatom;
    if (z[i] == 7 || z[i] == 8 || z[i] == 16) {
      atype[i] |= kIsONS;
    }
  }

  return atype;
}

static int
do_compute_distance_matrix_descriptors(Molecule & m,
                                       const atomic_number_t * z,
                                       const int * ncon,
                                       int * eccentricity,
                                       const int * dm)
{
  const auto matoms = m.natoms();
  if (1 == matoms)    // these parameters don't make sense
    return 1;

  int * totd = eccentricity + matoms;

  int tm = 0;                 // terminal methyl groups
  int tg3 = 0;                // terminal groups separated by 3 bonds
  int max_distance_between_ring_atoms = 0;
  int max_distance_between_aromatic_atoms = 0;
  int max_distance_between_heteratoms = 0;
  int max_distance_between_ons = 0;
  double harary = 0.0;

  Accumulator_Int<int> bonds_between_stats;

  std::unique_ptr<uint32_t[]> atype = AssignDistanceMatrixAtomTypes(m, z);

  for (int i = 0; i < matoms; i++) {
    const int tg = (1 == ncon[i]);     // is this a terminal group

    const int i_is_ring = atype[i] & kIsRing;
    const int i_is_aromatic = atype[i] & kIsAromatic;
    const int i_is_heteroatom = atype[i] & kIsHeteroatom;
    const int i_is_ons = atype[i] & kIsONS;

    if (tg && 6 == z[i] && 1 == m.nbonds(i)) {
      tm++;
    }

    for (int j = i + 1; j < matoms; j++) {
      int d = dm[i * matoms + j];

      bonds_between_stats.extra(d);

      harary += 1.0 / static_cast<float>(d);

      totd[i] += d;
      totd[j] += d;

      if (tg && 1 == ncon[j] && 3 == d)
        tg3++;
        
      if (d > eccentricity[i])
        eccentricity[i] = d;
      if (d > eccentricity[j])
        eccentricity[j] = d;

      if (i_is_ring && d > max_distance_between_ring_atoms && (atype[j] & kIsRing)) {
        max_distance_between_ring_atoms = d;
      }
      if (i_is_aromatic && d > max_distance_between_aromatic_atoms && (atype[j] & kIsAromatic)) {
        max_distance_between_aromatic_atoms = d;
      }
      if (i_is_heteroatom && d > max_distance_between_heteratoms && (atype[j] & kIsHeteroatom)) {
        max_distance_between_heteratoms = d;
      }
      if (i_is_ons && d > max_distance_between_ons && (atype[j] & kIsONS)) {
        max_distance_between_ons = d;
      }
    }
  }

  descriptor[iwdescr_harary].set(harary);
  descriptor[iwdescr_weiner].set(static_cast<float>(bonds_between_stats.sum()));

  int max_eccentricity = 0;     // same as longest path
  int muldiam = 0;
  int min_eccentricity = matoms;
  int mulrad = 0;

  int min_tot_d = totd[0];
  Set_of_Atoms atoms_at_min_tot_d;

  int max_tot_d = totd[0];
//atom_number_t atom_max_tot_d  = 0;

  for (int i = 0; i < matoms; i++)
  {
    int ecc = eccentricity[i];

    if (ecc > max_eccentricity)
    {
      max_eccentricity = ecc;
      muldiam = 1;
    } else if (ecc == max_eccentricity) {
      muldiam++;
    }

    if (ecc < min_eccentricity)
    {
      min_eccentricity = ecc;
      mulrad = 1;
    }
    else if (ecc == min_eccentricity)
      mulrad++;

    if (totd[i] > max_tot_d)
    {
      max_tot_d = totd[i];
//    atom_max_tot_d = i;
    }
    else if (totd[i] < min_tot_d)
    {
      min_tot_d = totd[i];
      atoms_at_min_tot_d.resize_keep_storage(0);
      atoms_at_min_tot_d.add(i);
    }
    else if (totd[i] == min_tot_d)
      atoms_at_min_tot_d.add(i);
  }

  descriptor[iwdescr_muldiam].set(static_cast<float>(muldiam));
  descriptor[iwdescr_rad].set    (static_cast<float>(min_eccentricity));
  descriptor[iwdescr_mulrad].set (static_cast<float>(mulrad));

  descriptor[iwdescr_tm].set     (static_cast<float>(tm));
  descriptor[iwdescr_tg3].set    (static_cast<float>(tg3));

  descriptor[iwdescr_mxdst].set(static_cast<float>(max_eccentricity));
  descriptor[iwdescr_fmxdst].set(iwmisc::Fraction<float>(max_eccentricity, matoms));

  if (min_eccentricity > 0) {
    descriptor[iwdescr_ishape].set(static_cast<float>(max_eccentricity - min_eccentricity) / static_cast<float>(min_eccentricity));
  }

  descriptor[iwdescr_maxdrng].set(max_distance_between_ring_atoms);
  descriptor[iwdescr_maxdarom].set(max_distance_between_aromatic_atoms);
  descriptor[iwdescr_maxdhtro].set(max_distance_between_heteratoms);
  descriptor[iwdescr_maxdons].set(max_distance_between_ons);
  descriptor[iwdescr_avebbtwn].set(bonds_between_stats.average_if_available_minval_if_not());
  descriptor[iwdescr_normbbtwn].set(bonds_between_stats.average_if_available_minval_if_not() / static_cast<float>(matoms));

  descriptor[iwdescr_compact].set(1.0f - static_cast<float>(max_eccentricity) / static_cast<float>(matoms));
  descriptor[iwdescr_nolp].set(matoms - max_eccentricity -1);

  descriptor[iwdescr_avdcentre].set(static_cast<float>(min_tot_d) / static_cast<float>(matoms - 1));

  Accumulator_Int<int> from_middle;

  int atoms_within_3_bonds = 0;
  int heteroatoms_within_3_bonds = 0;

  for (int i = 0; i < atoms_at_min_tot_d.number_elements(); ++i)
  {
    const auto j = atoms_at_min_tot_d[i];
    for (int k = 0; k < matoms; ++k)
    {
      if (k == j)
        continue;

      int d = dm[matoms * j + k];

      from_middle.extra(d);

      if (d > 3)
        continue;

      atoms_within_3_bonds++;

      if (6 != z[k])
        heteroatoms_within_3_bonds++;
    }
  }

  descriptor[iwdescr_stddcentre].set(sqrt(from_middle.variance()));
  descriptor[iwdescr_centre3].set(static_cast<float>(atoms_within_3_bonds));
  descriptor[iwdescr_centre3h].set(static_cast<float>(heteroatoms_within_3_bonds));

//set_vector(totd, matoms, 0);

  // Actually this is wrong, it is measuring heteroatoms further
  // than 3 bonds. Turns out this is useful this way.
  int largest_heteroatoms_within_three_bonds = 0;

  for (auto i = 0; i < matoms; ++i)
  {
    int heteroatoms = (6 != z[i]);

    for (auto j = 0; j < matoms; ++j)
    {
      if (j == i)
        continue;

//    int d = m.bonds_between(i, j);
      int d = dm[i * matoms + j];

      // Wrong, but useful as is.
      if (d < 3)
        continue;

      if (6 != z[j])
        heteroatoms++;
    }

    if (heteroatoms > largest_heteroatoms_within_three_bonds)
      largest_heteroatoms_within_three_bonds = heteroatoms;
  }

  descriptor[iwdescr_mh3b].set(static_cast<float>(largest_heteroatoms_within_three_bonds));

  return max_eccentricity;
}

static int
do_compute_distance_matrix_descriptors(Molecule & m,
                                       const atomic_number_t * z,
                                       const int * ncon)
{
  const auto matoms = m.natoms();

  m.recompute_distance_matrix();

  int * eccentricity = new_int(matoms + matoms);std::unique_ptr<int[]> free_eccentricity(eccentricity);

  int * dm = new int[matoms * matoms]; std::unique_ptr<int[]> free_dm(dm);
  assert (nullptr != dm);

  std::copy_n(m.distance_matrix_warning_may_change(), matoms * matoms, dm);

  int longest_path = do_compute_distance_matrix_descriptors(m, z, ncon, eccentricity, dm);

  do_compute_distances_from_longest_path_descriptors(m, dm, longest_path, eccentricity);   // re-use eccentricity array

  do_compute_descriptors_related_to_centroid_atom(m, dm, longest_path, z, ncon, eccentricity);   // eccentricity just used as a temporary array, existing contents destroyed

  do_compute_mean_shell_occupancies(m, dm);

  return longest_path;
}
                                    
static int
compute_polar_bond_descriptors(Molecule & m,
                               const atomic_number_t * z,
                               const int * ncon)
{
  int nb = m.nedges();

  int polar_bond_count = 0;
  int aromatic_polar_bond_count = 0;
  int aromatic_non_polar_bond_count = 0;
  int polar_multiple_bonds = 0;
  int di_vinyl_bonds = 0;      // matches the single bond in  *=*-*=*, not in an aromatic ring

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (z[a1] != z[a2])
      polar_bond_count++;

    int isar = m.in_same_aromatic_ring(a1, a2);

    if (isar)     // in same aromatic ring(s)
    {
      if (z[a1] == z[a2])
        aromatic_non_polar_bond_count++;
      else
        aromatic_polar_bond_count++;
    }
    else if (z[a1] == z[a2])     // a non-polar bond
      ;
    else if (! b->is_single_bond())
      polar_multiple_bonds++;

    if (! isar && b->is_single_bond() &&
        (ncon[a1] < m.nbonds(a1)) && (ncon[a2] < m.nbonds(a2)))
      di_vinyl_bonds++;
  }

  descriptor[iwdescr_pbcount].set( static_cast<float>(polar_bond_count));
  if (nb > 0)
    descriptor[iwdescr_frpbond].set( static_cast<float>(polar_bond_count) / static_cast<float>(nb));
  descriptor[iwdescr_nonpbond].set( static_cast<float>(nb - polar_bond_count));
  descriptor[iwdescr_pbarom].set( static_cast<float>(aromatic_polar_bond_count));
  descriptor[iwdescr_npbarom].set( static_cast<float>(aromatic_non_polar_bond_count));
  descriptor[iwdescr_pbunset].set( static_cast<float>(polar_multiple_bonds));
  descriptor[iwdescr_dvinylb].set( static_cast<float>(di_vinyl_bonds));

  return 1;
}

/*
  Two >3 connected atoms connected to each other contributes 1.
  If there is an atom with 2 connections between them, the contribution is 0.5
*/

static int
compute_crowding_descriptors(Molecule & m,
                             const int * ncon,
                             Atom ** atom)
{
  const int matoms = m.natoms();

  float rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    if (ncon[i] < 3)
      continue;

    const Atom * ai = atom[i];
    for (int j = 0; j < ncon[i]; j++)
    {
      const atom_number_t k = ai->other(i, j);

      if (ncon[k] > 2)      // two crowded atoms directly connected
        rc += 1.0;
      else if (2 == ncon[k])   // look for a 3 connected atom attached to K
      {
        const Atom * ak = atom[k];

        for (int l = 0; l < 2; l++)
        {
          const atom_number_t n = ak->other(k, l);

          if (n == i)
            continue;

          if (ncon[n] > 2)
            rc += 0.5;

          break;     // there are only two connections to check
        }
      }
    }
  }

  rc = rc * 0.5;    // the loop above double counts

  descriptor[iwdescr_crowding].set(rc);
  descriptor[iwdescr_fcrowdng].set(rc / static_cast<float>(matoms));

  return 1;
}

static int
compute_crowding_descriptors(Molecule & m,
                             const int * ncon)
{
  const Atom ** a = new const Atom *[m.natoms()]; std::unique_ptr<const Atom *[]> free_a(a);

  m.atoms(a);

  return compute_crowding_descriptors(m, ncon, (Atom **) a);
}

/*
  Various function supporting the inter-feature hydrogen bonding descriptors
*/

/*
  We have computed distances between Hydrogen bonding features and perhaps
  there were none found - or just one. We need to set the reported values
  for these missing distances to 0
*/

static void
check_for_no_values(int & shortest, int & nshortest, const int matoms)
{
  if (matoms == shortest)
    shortest = nshortest = 0;
  else if (matoms == nshortest)
    nshortest = 0;

  return;
}

static void
compute_average(Accumulator_Int<int> & acc, double & zresult)
{
  if (0 == acc.n())
  {
    zresult = 0.0;
    return;
  }

  if (1 == acc.n())
    zresult = acc.sum();
  else
    zresult = acc.average();

  return;
}

/*
  We are trying to keep track of the shortest and next-to-shortest distances
*/

static void
check_against_two(int tocheck, int & shortest, int & nshortest)
{
  if (tocheck < shortest)
  {
    nshortest = shortest;
    shortest = tocheck;
  }
  else if (tocheck < nshortest)
    nshortest = tocheck;

  return;
}

static int
compute_bond_separation_parameters(Molecule & m,
                                   const int * isotope)
{
  int matoms = m.natoms();

  int shortest_aa = matoms;
  int shortest_ad = matoms;
  int shortest_dd = matoms;

  int nshortest_aa = matoms;
  int nshortest_ad = matoms;
  int nshortest_dd = matoms;

  Accumulator_Int<int> aa, dd, ad;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == isotope[i])
      continue;

    boolean i_is_acceptor = (isotope[i] >= 2);

    for (int j = i + 1; j < matoms; j++)
    {
      if (0 == isotope[j])
        continue;

      int d = m.bonds_between(i, j);

      if (d < min_hbond_feature_separation || d > max_hbond_feature_separation)
        continue;

      boolean j_is_acceptor = (isotope[j] >= 2);

      if (i_is_acceptor && j_is_acceptor)
      {
        aa.extra(d);
        check_against_two(d, shortest_aa, nshortest_aa);
      }
      if (! i_is_acceptor && ! j_is_acceptor)
      {
        dd.extra(d);
        check_against_two(d, shortest_dd, nshortest_dd);
      }
      else
      {
        ad.extra(d);
        check_against_two(d, shortest_ad, nshortest_ad);
      }
    }
  }

  check_for_no_values(shortest_aa, nshortest_aa, matoms);
  check_for_no_values(shortest_ad, nshortest_ad, matoms);
  check_for_no_values(shortest_dd, nshortest_dd, matoms);

  double ave_aa, ave_ad, ave_dd;
  compute_average(aa, ave_aa);
  compute_average(ad, ave_ad);
  compute_average(dd, ave_dd);

  descriptor[iwdescr_aamind].set(static_cast<float>(shortest_aa));
  descriptor[iwdescr_aa2mind].set(static_cast<float>(nshortest_aa));
  descriptor[iwdescr_aaave].set(static_cast<float>(ave_aa));

  descriptor[iwdescr_admind].set(static_cast<float>(shortest_ad));
  descriptor[iwdescr_ad2mind].set(static_cast<float>(nshortest_ad));
  descriptor[iwdescr_adave].set(static_cast<float>(ave_ad));

  descriptor[iwdescr_ddmind].set(static_cast<float>(shortest_dd));
  descriptor[iwdescr_dd2mind].set(static_cast<float>(nshortest_dd));
  descriptor[iwdescr_ddave].set(static_cast<float>(ave_dd));

  return 1;
}

#define MAX_RING_SIZE 9

/*
  We count the number of types of ring fusions.
  These are laid out in the descriptor array as:

  5l5l    
  5l5r
  5r5r
  5l6l
  5l6r
  5r6l
  5r6r
  6l6l
  6l6r
  6r6r

  Where l means aLiphatic and r means aRomatic
*/

static int
compute_ring_fusion_descriptors(Molecule & m,
                                const int * ncon,
                                const int * ring_membership)
{
  int ar5 = 0;
  int ar6 = 0;
  int al5 = 0;
  int al6 = 0;

  int fsdrng5l5l = 0;
  int fsdrng5l5r = 0;
  int fsdrng5r5r = 0;
  int fsdrng5l6l = 0;
  int fsdrng5l6r = 0;
  int fsdrng5r6l = 0;
  int fsdrng5r6r = 0;
  int fsdrng6l6l = 0;
  int fsdrng6l6r = 0;
  int fsdrng6r6r = 0;

  m.compute_aromaticity_if_needed();

  int nr = m.nrings();

  if (nr < 2) {
    return 1;
  }

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    int rsi = ri->number_elements();

    if (rsi < 5 || rsi > 6)
      continue;

    int n = ri->fused_ring_neighbours();

    if (0 == n)    // isolated ring
    {
      if (ri->is_aromatic())
      {
        if (6 == rsi)
          ar6++;
        else if (5 == rsi)
          ar5++;
      }
      else
      {
        if (6 == rsi)
          al6++;
        else if (5 == rsi)
          al5++;
      }
      continue;
    }

    if (fsdrng_descriptors_consider_just_two_ring_systems && n > 1)
      continue;

    for (int j = 0; j < n; j++)
    {
      const Ring * rj = ri->fused_neighbour(j);

      if (rj->ring_number() < ri->ring_number())
        continue;

      int rsj = rj->number_elements();

      if (rsj < 5 || rsj > 6)
        continue;

//    cerr << "lbswar " << ri->largest_number_of_bonds_shared_with_another_ring() << '\n';
//    cerr << "Shared " << ri->compute_bonds_shared_with(*rj) << '\n';

      if (ri->largest_number_of_bonds_shared_with_another_ring() > 1 &&   // strongly fused
          ri->compute_bonds_shared_with(*rj) > 1)
        continue;

      if (fsdrng_descriptors_consider_just_two_ring_systems && rj->fused_ring_neighbours() > 1)
        continue;

      if (ri->is_aromatic() && rj->is_aromatic())
      {
        if (5 == rsi && 5 == rsj)
          fsdrng5r5r++;
        else if (6 == rsi && 6 == rsj)
          fsdrng6r6r++;
        else
          fsdrng5r6r++;
      }
      else if (! ri->is_aromatic() && ! rj->is_aromatic())   // both aliphatic
      {
        if (5 == rsi && 5 == rsj)
          fsdrng5l5l++;
        else if (6 == rsi && 6 == rsj)
          fsdrng6l6l++;
        else
          fsdrng5l6l++;
      }
      else    // rings have different aromaticity
      {
        if (5 == rsi && 5 == rsj)               // both of size 5
          fsdrng5l5r++;
        else if (6 == rsi && 6 == rsj)          // both of size 6
          fsdrng6l6r++;
        else if (ri->is_aromatic())   // therefore rj is aliphatic
        {
          if (5 == rsi)
            fsdrng5r6l++;
          else
            fsdrng5l6r++;
        }
        else     // rj is aromatic and ri is aliphatic
        {
          if (5 == rsi)
            fsdrng5l6r++;
          else
            fsdrng5r6l++;
        }
      }
    }
  }

  descriptor[iwdescr_ar5].set(ar5);
  descriptor[iwdescr_ar6].set(ar6);
  descriptor[iwdescr_al5].set(al5);
  descriptor[iwdescr_al6].set(al6);

  descriptor[iwdescr_fsdrng5l5l].set(fsdrng5l5l);
  descriptor[iwdescr_fsdrng5l5r].set(fsdrng5l5r);
  descriptor[iwdescr_fsdrng5r5r].set(fsdrng5r5r);
  descriptor[iwdescr_fsdrng5l6l].set(fsdrng5l6l);
  descriptor[iwdescr_fsdrng5l6r].set(fsdrng5l6r);
  descriptor[iwdescr_fsdrng5r6l].set(fsdrng5r6l);
  descriptor[iwdescr_fsdrng5r6r].set(fsdrng5r6r);
  descriptor[iwdescr_fsdrng6l6l].set(fsdrng6l6l);
  descriptor[iwdescr_fsdrng6l6r].set(fsdrng6l6r);
  descriptor[iwdescr_fsdrng6r6r].set(fsdrng6r6r);

  descriptor[iwdescr_fsdrngarar].set(fsdrng5r5r + fsdrng5r6r + fsdrng6r6r);
  descriptor[iwdescr_fsdrngalar].set(fsdrng5l5r + fsdrng5l6r + fsdrng6l6r);
  descriptor[iwdescr_fsdrngalal].set(fsdrng5l5l + fsdrng5l6l + fsdrng6l6l);

  return 1;
}

//#define DEBUG_ATOMS_NOT_ALREADY_MARKED

/*
  We are trying to determine the largest number of atoms in a ring system.
  We have a new ring. 
*/

static int
atoms_not_already_marked(const Ring & r,
                         int * atom_already_done)
{
  int rc = 0;       // the number of items not already set

  int ring_size = r.number_elements();

#ifdef DEBUG_ATOMS_NOT_ALREADY_MARKED
  cerr << "Checking ring of size " << ring_size << '\n';
#endif

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];

#ifdef DEBUG_ATOMS_NOT_ALREADY_MARKED
    if (atom_already_done[j])
      cerr << "Atom " << j << " already included\n";
    else
      cerr << "Atom " << j << " being added\n";
#endif

    if (atom_already_done[j])
      continue;

    atom_already_done[j] = 1;
    rc++;
  }

  return rc;
}

#ifdef OLD_VERSION_NO_LONGER_USED
static int
compute_exocyclic_bonds(Molecule & m,
                        const Ring & r,
                        const int *ncon)
{
  int rc = 0;

  int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];

    if (2 == ncon[j])
      continue;

    const Atom * a = m.atomi(j);

    for (int k = 0; k < ncon[j]; k++)
    {
      const Bond * b = a->item(k);

      if (0 == b->nrings())
        rc++;
    }
  }

  return rc;
}
#endif

static int
compute_exocyclic_bonds(Molecule & m,
                        const atomic_number_t * z,
                        const Ring & r,
                        const int *ncon,
                        int & double_bond_attachments,
                        int & singly_connected_attachments,
                        int & singly_connected_heteroatoms,
                        int & singly_connected_donors)
{
  int rc = 0;

  int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];

    if (2 == ncon[j])
      continue;

    const Atom * a = m.atomi(j);

    for (int k = 0; k < ncon[j]; k++)
    {
      const Bond * b = a->item(k);

      if (b->nrings())
        continue;

      rc++;

      const atom_number_t o = b->other(j);

      if (1 != ncon[o])
        continue;

      if (b->is_double_bond())
        double_bond_attachments++;

      singly_connected_attachments++;
      if (6 == z[o])
        continue;

      singly_connected_heteroatoms++;
      if (m.hcount(o))
        singly_connected_donors++;
    }
  }

  return rc;
}

/*
  Atom A is a member of ring R. Scan the atoms attached to A and determine
  whether or not they are all terminal groups
*/

static int
just_terminal_groups_outside_ring(const Molecule & m,
                                  const int * ncon,
                                  const Ring & r,
                                  atom_number_t a)
{
  const Atom * ai = m.atomi(a);

  for (int i = 0; i < ncon[a]; i++)
  {
    atom_number_t j = ai->other(a, i);

    if (1 == ncon[j])     // must be terminal group outside the ring
      continue;

    if (r.contains(j))
      continue;

    if (ncon[j] > 1)     // this connection is going somewhere
      return 0;
  }

  return 1;     // Yep, just singly connected things branching off
}

static float
compute_ring_isolation(Molecule & m,
                       const int * ncon,
                       const Ring & r)
{
  int ring_size = r.number_elements();

  int branches = 0;

  int terminal_groups_found_here = 0;
  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];
    if (2 == ncon[j])         // must be in the ring
      continue;

    branches++;

    if (just_terminal_groups_outside_ring(m, ncon, r, j))
      terminal_groups_found_here++;
  }

  branches -= terminal_groups_found_here;

  if (0 == branches)     // this ring must be the whole molecule
    return 2.0;

  return 1.0 / static_cast<float>(branches);
}

static int
heteroatoms_in_ring(const Ring * r,
                    const atomic_number_t * z)
{
  int rc = 0;

  int rs = r->number_elements();
  for (int i = 0; i < rs; i++)
  {
    atom_number_t j = r->item(i);
    if (6 != z[j])
      rc++;
  }

  return rc;
}

/*
  We compute the number of each ring size
*/

static int
compute_ring_descriptors(Molecule & m,
               const atomic_number_t * z,
               const int * ncon,
               int * ring_already_done,
               int * atom_already_done)
{
  int matoms = m.natoms();

  int nrings[MAX_RING_SIZE];
  set_vector(nrings, MAX_RING_SIZE, 0);

  int rings_in_largest_system = 1;
  int atoms_in_largest_system = 0;
  int ring_systems_containing_multiple_rings = 0;

  int isolated_rings = 0;
  int isolated_heterocycles = 0;

  int max_heteroatoms_in_ring = 0;
  float max_ring_heteroatom_fraction = 0.0;
  float min_ring_heteroatom_fraction = 1.0;

  float ring_isolation_score = 0.0;

  int aromatic_rings = 0;
  int aliphatic_rings = 0;

  int smallest_ring_size = 0;
  int largest_ring_size = 0;

  int exocyclic_bonds = 0;
  int double_bond_attachments = 0;
  int singly_connected_attachments = 0;
  int singly_connected_heteroatoms = 0;
  int singly_connected_donors = 0;

  int more_than_7_atoms = 0;

  int nr = m.nrings();

  m.compute_aromaticity_if_needed();

  int ring_systems = nr;

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (ri->is_aromatic())
      aromatic_rings++;
    else
      aliphatic_rings++;

    if (! ri->is_fused())
    {
      float tmp = compute_ring_isolation(m, ncon, *ri);
      ring_isolation_score += tmp;
    }

    int hac = heteroatoms_in_ring(ri, z);
    if (hac > max_heteroatoms_in_ring)
      max_heteroatoms_in_ring = hac;

    exocyclic_bonds += compute_exocyclic_bonds(m, z, *ri, ncon, double_bond_attachments, 
                                               singly_connected_attachments, singly_connected_heteroatoms,
                                               singly_connected_donors);

    if (! ri->is_fused())
    {
      isolated_rings++;
      if (hac > 0)
        isolated_heterocycles++;
    }

    int rs = ri->number_elements();

    if (smallest_ring_size == 0) {
      smallest_ring_size = rs;
    }
    if (rs > largest_ring_size) {
      largest_ring_size = rs;
    }
    if (rs > 7) {
      ++more_than_7_atoms;
    }

    float tmp = static_cast<float>(hac) / static_cast<float>(rs);
    if (tmp > max_ring_heteroatom_fraction)
      max_ring_heteroatom_fraction = tmp;

    if (tmp < min_ring_heteroatom_fraction)
      min_ring_heteroatom_fraction = tmp;

    if (rs < MAX_RING_SIZE)
      nrings[rs]++;

    if (ring_already_done[i])
      continue;

    int system_size = 1;
    int atoms_in_system = ri->number_elements();

    if (ri->is_fused())
    {
      set_vector(atom_already_done, matoms, 0);

      ri->set_vector(atom_already_done, 1);

      for (int j = i + 1; j < nr; j++)
      {
        const Ring * rj = m.ringi(j);
        if (rj->fused_system_identifier() == ri->fused_system_identifier())
        {
          system_size++;
          atoms_in_system += atoms_not_already_marked(*rj, atom_already_done);
          ring_already_done[j] = 1;
          ring_systems--;
        }
      }
    }

#ifdef DEBUG_ATOMS_NOT_ALREADY_MARKED
    cerr << "System contains " << system_size << " rings and " << atoms_in_system << " atoms\n";
#endif

    if (system_size > rings_in_largest_system)
      rings_in_largest_system = system_size;

    if (system_size > 1)
      ring_systems_containing_multiple_rings++;

    if (atoms_in_system > atoms_in_largest_system)
      atoms_in_largest_system = atoms_in_system;
  }

// for some descriptors, examine the non sssr rings too

  for (int i = 0; i < m.non_sssr_rings(); i++)
  {
    const Ring * ri = m.non_sssr_ring(i);

    int hac = heteroatoms_in_ring(ri, z);
    if (hac > max_heteroatoms_in_ring)
      max_heteroatoms_in_ring = hac;

    int rs = ri->number_elements();

    float tmp = static_cast<float>(hac) / static_cast<float>(rs);
    if (tmp > max_ring_heteroatom_fraction)
      max_ring_heteroatom_fraction = tmp;

//  if (tmp < min_ring_heteroatom_fraction)
//    min_ring_heteroatom_fraction = tmp;
  }

  descriptor[iwdescr_mhr].set(static_cast<float>(max_heteroatoms_in_ring));
  descriptor[iwdescr_mxhrf].set(static_cast<float>(max_ring_heteroatom_fraction));
  descriptor[iwdescr_mnhrf].set(static_cast<float>(min_ring_heteroatom_fraction));
  descriptor[iwdescr_lrsysz].set(static_cast<float>(rings_in_largest_system));
  descriptor[iwdescr_srsz].set(static_cast<float>(smallest_ring_size));
  descriptor[iwdescr_lrsz].set(static_cast<float>(largest_ring_size));
  descriptor[iwdescr_rng7atoms].set(static_cast<float>(more_than_7_atoms));

  descriptor[iwdescr_nrsyscmr].set(static_cast<float>(ring_systems_containing_multiple_rings));
  descriptor[iwdescr_mars].set(static_cast<float>(atoms_in_largest_system));
  descriptor[iwdescr_ringsys].set(static_cast<float>(ring_systems));
  descriptor[iwdescr_ringisol].set(ring_isolation_score);
  descriptor[iwdescr_isolrc].set(isolated_rings);
  descriptor[iwdescr_isolhtrc].set(isolated_heterocycles);

  descriptor[iwdescr_arring].set(aromatic_rings);
  descriptor[iwdescr_alring].set(aliphatic_rings);

//descriptor[iwdescr_scra].set(static_cast<float>(singly_connected_attachments));
//descriptor[iwdescr_scrha].set(static_cast<float>(singly_connected_heteroatoms));
//descriptor[iwdescr_scrd].set(static_cast<float>(singly_connected_donors));

  descriptor[iwdescr_excybond].set(static_cast<float>(exocyclic_bonds));
  descriptor[iwdescr_excydbond].set(static_cast<float>(double_bond_attachments));
  descriptor[iwdescr_excydscon].set(static_cast<float>(singly_connected_attachments));
  descriptor[iwdescr_excydsconh].set(static_cast<float>(singly_connected_heteroatoms));
  descriptor[iwdescr_excydscondon].set(static_cast<float>(singly_connected_donors));

  int offset = iwdescr_nrings3;    // dangerous case from enumeration to int

  for (int i = 3; i < MAX_RING_SIZE; i++)
  {
    descriptor[offset].set(static_cast<float>(nrings[i]));
    offset++;
  }

  return 1;
}

static int
compute_ring_descriptors(Molecule & m,
               const atomic_number_t * z,
               const int * ncon)
{
  assert (m.ok());

  int nr = m.nrings();

  // Initialise ring descriptors, missing values do not make sense.
  // Note that some of these may not get computed depending on
  // what has been requested...
  descriptor[iwdescr_trmnlrng].set(0.0);
  descriptor[iwdescr_intrnlrng].set(0.0);
  descriptor[iwdescr_rng2spch].set(0.0);
  descriptor[iwdescr_rng2bridge].set(0.0);

  // Various fused ring descriptors do not get computed
  // if there is only 1 ring in the molecule.
  if (nr < 2) {
    descriptor[iwdescr_fsdrng5l5l].set(0.0);
    descriptor[iwdescr_fsdrng5l5r].set(0.0);
    descriptor[iwdescr_fsdrng5r5r].set(0.0);
    descriptor[iwdescr_fsdrng5l6l].set(0.0);
    descriptor[iwdescr_fsdrng5l6r].set(0.0);
    descriptor[iwdescr_fsdrng5r6l].set(0.0);
    descriptor[iwdescr_fsdrng5r6r].set(0.0);
    descriptor[iwdescr_fsdrng6r6r].set(0.0);
    descriptor[iwdescr_fsdrng6l6r].set(0.0);
    descriptor[iwdescr_fsdrng6l6l].set(0.0);

    descriptor[iwdescr_fsdrngarar].set(0.0);
    descriptor[iwdescr_fsdrngalar].set(0.0);
    descriptor[iwdescr_fsdrngalal].set(0.0);

    descriptor[iwdescr_nsfsdsys].set(0.0);
    descriptor[iwdescr_rnginsfs].set(0.0);
    descriptor[iwdescr_lgstrfsy].set(0.0);
    descriptor[iwdescr_htrcsfsy].set(0.0);
    descriptor[iwdescr_mxhtsfsy].set(0.0);
  }


  if (0 == nr) {
    descriptor[iwdescr_mhr].set(0.0);
    descriptor[iwdescr_mxhrf].set(0.0);
    descriptor[iwdescr_mnhrf].set(0.0);
    descriptor[iwdescr_lrsysz].set(0.0);
    descriptor[iwdescr_lrsz].set(0.0);
    descriptor[iwdescr_rng7atoms].set(0.0);
    descriptor[iwdescr_mars].set(0.0);
    descriptor[iwdescr_ringsys].set(0.0);
    descriptor[iwdescr_ringisol].set(0.0);
    descriptor[iwdescr_isolrc].set(0.0);
    descriptor[iwdescr_isolhtrc].set(0.0);
    descriptor[iwdescr_arring].set(0.0);
    descriptor[iwdescr_alring].set(0.0);
    descriptor[iwdescr_excybond].set(0.0);
    descriptor[iwdescr_excydbond].set(0.0);
    descriptor[iwdescr_excydscon].set(0.0);
    descriptor[iwdescr_excydsconh].set(0.0);
    descriptor[iwdescr_excydscondon].set(0.0);

    descriptor[iwdescr_nrings3].set(0.0);
    descriptor[iwdescr_nrings4].set(0.0);
    descriptor[iwdescr_nrings5].set(0.0);
    descriptor[iwdescr_nrings6].set(0.0);
    descriptor[iwdescr_nrings7].set(0.0);
    descriptor[iwdescr_nrings8].set(0.0);
  
    descriptor[iwdescr_rsarom1].set(0.0);
    descriptor[iwdescr_rsarom2].set(0.0);
    descriptor[iwdescr_rsarom3].set(0.0);

    descriptor[iwdescr_rsaliph1].set(0.0);
    descriptor[iwdescr_rsaliph2].set(0.0);
    descriptor[iwdescr_rsaliph3].set(0.0);
    descriptor[iwdescr_rsaliph4].set(0.0);


    descriptor[iwdescr_rssys1].set(0.0);
    descriptor[iwdescr_rssys2].set(0.0);
    descriptor[iwdescr_rssys3].set(0.0);
    descriptor[iwdescr_rssys4].set(0.0);
    descriptor[iwdescr_rssys5].set(0.0);
    descriptor[iwdescr_rssys6].set(0.0);
    descriptor[iwdescr_rssys7].set(0.0);
    descriptor[iwdescr_rssys8].set(0.0);
    descriptor[iwdescr_rssys9].set(0.0);

    descriptor[iwdescr_ar5].set(0.0);
    descriptor[iwdescr_ar6].set(0.0);
    descriptor[iwdescr_al5].set(0.0);
    descriptor[iwdescr_al6].set(0.0);
  
    // If the molecule has no rings, then the spinach is not defined.

    descriptor[iwdescr_rhaf].set(0.0);
    descriptor[iwdescr_frafus].set(0.0);
    descriptor[iwdescr_fraromha].set(0.0);
    descriptor[iwdescr_srsz].set(0.0);
    descriptor[iwdescr_nrsyscmr].set(0.0);

    descriptor[iwdescr_spchtro].set(0.0);
    descriptor[iwdescr_rbfrspch].set(0.0);
    descriptor[iwdescr_satspcha].set(0.0);
    descriptor[iwdescr_unsatspcha].set(0.0);
    descriptor[iwdescr_fsatspcha].set(0.0);
    descriptor[iwdescr_scaffoldbranches].set(0.0);
    descriptor[iwdescr_nrnspch].set(0.0);
    descriptor[iwdescr_fnrnspc].set(0.0);

    descriptor[iwdescr_avcdbsub].set(0.0);

    descriptor[iwdescr_npfsdsys].set(0.0);
    descriptor[iwdescr_rnginpfs].set(0.0);
    descriptor[iwdescr_lgplnfsy].set(0.0);
    descriptor[iwdescr_htrcpfsy].set(0.0);
    descriptor[iwdescr_mxhtpfsy].set(0.0);


    return 1;
  }

  int * ring_already_done = new_int(nr); std::unique_ptr<int[]> free_ring_already_done(ring_already_done);
  int * atom_already_done = new_int(m.natoms()); std::unique_ptr<int[]> free_atom_already_done(atom_already_done);

  return compute_ring_descriptors(m, z, ncon, ring_already_done, atom_already_done);
}

static int
compute_hbond_descriptors (Molecule & m,
               const atomic_number_t * z,
               const int * ncon,
               const Atom ** atom)
{
  int matoms = m.natoms();
  double acceptor_score = 0.0;
  double donor_score    = 0.0;
  double dual_score     = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    Atom * ai = const_cast<Atom *>(atom[i]);

    int nbi = ai->nbonds();
    int hci = ai->implicit_hydrogens();

    if (8 == z[i] && 1 == ncon[i] && 2 == nbi)    // carbonyl
      acceptor_score += 1.0;
    else if (8 == z[i] && 1 == ncon[i] && 1 == nbi)   // hydroxy
      dual_score += 1.0;
    else if (8 == z[i] && 2 == ncon[i])      // ether
      acceptor_score += 0.5;

    else if (7 == z[i] && 2 == hci)      // primary amine
      donor_score += 1.0;
    else if (7 == z[i] && 1 == hci)      // secondary amine
      donor_score += 1.0;
    else if (7 == z[i] && 2 == ncon[i] && 3 == nbi)     // =N-
      acceptor_score += 1.0;
    else if (7 == z[i] && 0 == hci)        // tertiary amine
      acceptor_score += 0.2;

    else if (16 == z[i] && 1 == hci)
      donor_score += 1.0;
    else if (16 == z[i] && 2 == ncon[i])
      acceptor_score += 0.2;
  }

  descriptor[iwdescr_hacts].set(static_cast<float>(acceptor_score));
  descriptor[iwdescr_hdons].set(static_cast<float>(donor_score));
  descriptor[iwdescr_hduals].set(static_cast<float>(dual_score));

  return 1;
}

/*
  For each substructure, we compute both the charge on each atom,
  together with the dipole components for the embedding.
*/

static int
compute_substructure_based_charge_descriptors(Molecule & m,
                                              IWString & output,
                                              Query_and_Charge_Stats * query,
                                              int ii)
{
  Substructure_Results sresults;

  int nmatches = query->substructure_search(&m, sresults);
  if (verbose)
    cerr << "query '" << query->comment() << "' matched " << nmatches << " times\n";

  if (0 == nmatches)
    return 1;

  query->tally_embedding(m, sresults);

  for (int i = 0; i < nmatches; i++)
  {
    const Set_of_Atoms * s = sresults.embedding(i);
    int atoms_matched = s->number_elements();
    charge_t sum_charge_on_matched_atoms = 0.0;

//  To compute the dipole, we keep track of the sum of the coordinates,
//  the sum of the charges, and the sum of charge * coord

    Coordinates xsum (0.0, 0.0, 0.0);     // sum of coordinates
    Coordinates qxsum(0.0, 0.0, 0.0);     // coords * charge
    Coordinates xxsum(0.0, 0.0, 0.0);     // coords * coords

    for (int j = 0; j < atoms_matched; j++)
    {
      atom_number_t a = s->item(j);
      output << ii << "_q_" << j;
      if (i > 0)
        output << "_" << i;
      output << " " << m.charge_on_atom(a) << '\n';

      Coordinates coords = * (m.atomi(a));

      xsum += coords;      // the sum of the coordinates
      charge_t q = m.charge_on_atom(a);
      coord_t tx = coords.x();
      coord_t ty = coords.y();
      coord_t tz = coords.z();
      qxsum.add(tx * q, ty * q, tz * q);  
      xxsum.add(tx * tx, ty * ty, tz * tz);
      sum_charge_on_matched_atoms += q;
    }

    xsum *= 1.0 / static_cast<float>(atoms_matched);   // now the average of the coordinates

//  Compute the "geometric" dipole

    Coordinates gdipole = xxsum - xsum * xsum * atoms_matched;

//  Compute the charge dipole
    
    Coordinates dipole = qxsum - xsum * sum_charge_on_matched_atoms;

    output << ii << "_dipx_" << i << " " << dipole.x() << '\n';
    double sum = dipole.x();

    output << ii << "_dipy_" << i << " " << dipole.y() << '\n';
    sum += dipole.y();

    output << ii << "_dipz_" << i << " " << dipole.z() << '\n';
    sum += dipole.z();

    output << ii << "_totdip_" << i << " " << sum << '\n';

    dipole.normalise();
    gdipole.normalise();

    output << ii << "_dip_mas_" << i << " " << dipole.angle_between_unit_vectors(gdipole) << '\n';
  }

  return 0;
}

static int
compute_radha_entropy_descriptors(Molecule & m,
                                  atom_number_t zatom,
                                  const int * ncon,
                                  const int * ring_membership,
                                  const Atom ** atom,
                                  int * already_done)
{
  const Atom * ai = atom[zatom];

  already_done[zatom] = 1;

  int rc = 0;

  for (int i = 0; i < ncon[zatom]; i++)
  {
    const Bond * b = ai->item(i);

    const atom_number_t j = b->other(zatom);

    if (already_done[j])
      continue;

    if (ring_membership[j])
      ;
    else if (ncon[j] > 2)    // don't follow this
      rc++;
    else if (1 == ncon[j])
      continue;
    else
    {
      rc++;   // atom J is part of the chain
      rc += compute_radha_entropy_descriptors(m, j, ncon, ring_membership, atom, already_done);
    }
  }

  return rc;
}

/*
  Although called ENTROPY, this focusses on flexible chains


*/

static int
compute_radha_entropy_descriptors(Molecule & m,
                                  const int * ncon,
                                  const int * ring_membership,
                                  const Atom ** atom,
                                  int * already_done)
{
  int matoms = m.natoms();

  int longest_chain = 0;
  int total_chain_atoms = 0;
  int nchains = 0;

  int non_ring_atoms = 0;

  double radha_entropy = 1.0;

  for (int i = 0; i < matoms; i++)
  {
    if (ring_membership[i])
      continue;

    non_ring_atoms++;

    if (already_done[i] || 2 != ncon[i])
      continue;

    int chain_length = compute_radha_entropy_descriptors(m, i, ncon, ring_membership, atom, already_done);

    if (0 == chain_length)
      continue;

    if (1 == chain_length)
      radha_entropy += 3.0;
    else if (2 == chain_length)
      radha_entropy += 9.0;
    else
      radha_entropy += pow(3.0, static_cast<double>(chain_length));

    nchains++;
    total_chain_atoms += chain_length;
    if (chain_length > longest_chain)
      longest_chain = chain_length;
  }

  descriptor[iwdescr_nflxchn].set(static_cast<float>(nchains));
  descriptor[iwdescr_lflxchn].set(static_cast<float>(longest_chain));
  descriptor[iwdescr_atflxchn].set(static_cast<float>(total_chain_atoms));

  if (nchains)
    descriptor[iwdescr_avflxchn].set(static_cast<float>(total_chain_atoms) / static_cast<float>(nchains));
  else
    descriptor[iwdescr_avflxchn].set(static_cast<float>(0.0));

  descriptor[iwdescr_faflxchn].set(static_cast<float>(total_chain_atoms) / static_cast<float>(matoms));

  if (non_ring_atoms)
    descriptor[iwdescr_fnflxchn].set(static_cast<float>(total_chain_atoms) / static_cast<float>(non_ring_atoms));
  else
    descriptor[iwdescr_fnflxchn].set(static_cast<float>(0.0));

  descriptor[iwdescr_rkentrpy].set(log(radha_entropy));

  return 1;
}

static int
compute_radha_entropy_descriptors(Molecule & m,
                                  const int * ncon,
                                  const int * ring_membership,
                                  const Atom ** atom)
{
  int * tmp = new_int(m.natoms()); std::unique_ptr<int[]> free_tmp(tmp);

  return compute_radha_entropy_descriptors(m, ncon, ring_membership, atom, tmp);
}

static int
compute_substructure_based_charge_descriptors (Molecule & m,
                             IWString & output)
{
  int rc = 0;

  for (int i = 0; i < charge_queries.number_elements(); i++)
  {
    rc += compute_substructure_based_charge_descriptors(m, output, charge_queries[i], i);
  }

  return rc;
}

static int
look_for_inter_ring_atoms(const Molecule & m,
                          atom_number_t avoid,
                          atom_number_t zatom, 
                          int * inter_ring,
                          const int * ring_membership,
                          const Atom ** atom)
{
  const Atom * a = atom[zatom];

  if (1 == a->ncon())    // did not find another ring atom
    return 0;

  int rc = 0;

  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);

    if (j == avoid)
      continue;

    if (ring_membership[j])   // found ring, great
    {
      inter_ring[zatom] = 1;
      rc = 1;
    }
    else if (look_for_inter_ring_atoms(m, zatom, j, inter_ring, ring_membership, atom))
    {
      inter_ring[zatom] = 1;
      rc = 1;
    }
  }

  return rc;
}

/*
  Identify all the atoms that bridge the rings
*/

static void
id_inter_ring_atoms(const Molecule & m,
                    int * inter_ring,
                    const int * ring_membership,
                    const Atom ** atom)
{
  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == ring_membership[i])
      continue;

    const Atom * a = atom[i];
    if (a->ncon() == 2) {
      continue;
    }

    for (const Bond* b : *a) {
      if (b->nrings())     // we are looking out from the ring
        continue;

      atom_number_t k = b->other(i);

      if (ring_membership[k])    // biphenyl type linkage
        continue;

      if (inter_ring[k])    // already processed
        continue;

      look_for_inter_ring_atoms(m, i, k, inter_ring, ring_membership, atom);
    }
  }

  return;
}

static int
id_spinach(const Molecule & m,
           int * spinach,
           const int * ring_membership,
           const Atom ** atom)
{
  id_inter_ring_atoms(m, spinach, ring_membership, atom);

  const int matoms = m.natoms();

// First task is to add doubly bonded, singly connected atoms to the ring and scaffold atoms
  for (int i = 0; i < matoms; i++)
  {
    if (ring_membership[i] || spinach[i])    // ring atoms or bridging atoms, add doubly bonded atoms
    {
      const Atom * a = atom[i];
      const int acon = a->ncon();
      if (acon > 2 && acon < a->nbonds())
      {
        for (const Bond* b: *a) {
          if (! b->is_double_bond()) {
            continue;
          }

          atom_number_t k = b->other(i);

          if (1 != atom[k]->ncon()) {
            continue;
          }

          spinach[k] = 1;     // will be included with scaffold
        }
      }
    }
  }

// Now that the doubly bonded atoms are taken care of, invert things

  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (ring_membership[i] || spinach[i])
      spinach[i] = 0;
    else
    {
      spinach[i] = 1;
      rc++;
    }
  }

  return rc;
}

/*
  We count the number of times a ring atom is bonded to a non-ring atom
*/

static int
compute_ring_chain_descriptors(Molecule & m,
                               const atomic_number_t * z,
                               const int * ncon,
                               const Atom ** atom,
                               const int * ring_membership)
{
  int ring_chain_join_count = 0;
  int ring_chain_heteroatom_join_count = 0;
  int aromatic_ring_chain_join_count = 0;
  int aliphatic_ring_chain_join_count = 0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
#ifdef DEBUG_RING_CHAIN_DESCRIPTORS
    if (0 == ring_membership[i])
      cerr << "Atom " << i << " not part of a ring\n";
#endif

    if (0 == ring_membership[i])     // atom I must be in a ring
      continue;

    if (2 == ncon[i])     // 2 connections in a ring, cannot have a chain joined
      continue;

    aromaticity_type_t aromi;
    (void) m.aromaticity(i, aromi);

#ifdef DEBUG_RING_CHAIN_DESCRIPTORS
    cerr << "Atom " << i << " aromaticity " << aromi << '\n';
#endif

    const Atom * ai = atom[i];

    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (ring_membership[k])
        continue;

#ifdef DEBUG_RING_CHAIN_DESCRIPTORS
      cerr << "    atom " << k << " is a non-ring attached atom\n";
#endif
      ring_chain_join_count++;
      if (6 != z[k])
        ring_chain_heteroatom_join_count++;

      if (is_aromatic_atom(aromi))
        aromatic_ring_chain_join_count++;
      else
        aliphatic_ring_chain_join_count++;
    }
  }

#ifdef DEBUG_RING_CHAIN_DESCRIPTORS
  cerr << " rcj = " << ring_chain_join_count << '\n';
  cerr << " amrcj = " << aromatic_ring_chain_join_count << '\n';
  cerr << " alrcj = " << aliphatic_ring_chain_join_count << '\n';
#endif

  descriptor[iwdescr_rcj].set(static_cast<float>(ring_chain_join_count));
  descriptor[iwdescr_rchj].set(static_cast<float>(ring_chain_heteroatom_join_count));
  descriptor[iwdescr_amrcj].set(static_cast<float>(aromatic_ring_chain_join_count));
  descriptor[iwdescr_alrcj].set(static_cast<float>(aliphatic_ring_chain_join_count));

  return 1;
}

// Is `zatom` within `m` a branch point within the scaffold.
// It must be a non ring atom. If it is joined to a doubly bonded
// singly attached atom, it is not a join.
static int
WithinScaffoldBranched(const Molecule& m, 
                       atom_number_t zatom) {
  const Atom& a = m[zatom];
  if (a.ncon() <= 2) {
    return 0;
  }

  int branches = 0;
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);
    // cerr << "From atom " << zatom << " to atom " << j << '\n';
    if (b->is_double_bond() && m.ncon(j) == 1) {
      continue;
    }
    ++branches;
  }

  // cerr << "Branches from " << zatom << " " << branches << '\n';

  assert (branches >= 2);

  // If just two attachments, no branching of the scaffold here.
  return branches - 2;
}

static void
count_spinach_connections(const Molecule & m,
                          const Ring & r,
                          const int * spinach,
                          int & spinach_connections, 
                          int & non_spinach_connections)
{
  spinach_connections = 0;
  non_spinach_connections = 0;

  int ring_size = r.number_elements();

  for (int i = 0; i < ring_size; i++)
  {
    atom_number_t j = r[i];

    const Atom * a = m.atomi(j);

    // no branching outside the ring here
    if (2 == a->ncon()) {
      continue;
    }

    for (const Bond* b : *a) {
      if (b->nrings())
        continue;

      atom_number_t l = b->other(j);

//    cerr << "Atom " << j << " in ring, to atom " << l << " ?spinach " << spinach[k] << '\n';

      if (spinach[l])
        spinach_connections++;
      else
        non_spinach_connections++;
    }
  }

  return;
}

//#define DEBUG_SPINACH_STUFF

static int
compute_spinach_descriptors(Molecule & m, int matoms,
                            int * spinach,
                            const atomic_number_t * z,
                            const int * ncon,
                            const int * ring_membership,
                            const Atom ** atom)
{
  assert (matoms > 0);
  assert (nullptr != spinach);

  if (0 == m.nrings()) {
    molecules_with_no_rings++;
    descriptor[iwdescr_frspch].set(static_cast<float>(1.0));

    return 1;
  }

  int spinach_atoms = id_spinach(m, spinach, ring_membership, atom);

#ifdef DEBUG_SPINACH_STUFF
  m.set_isotopes(spinach);
  cerr << m.smiles() << ' ' << m.name() << "\n";
  m.transform_to_non_isotopic_form();

//for (int i = 0; i < matoms; i++)
//{
//  cerr << "Atom " << i << " spinach " << spinach[i] << '\n';
//}
#endif

  descriptor[iwdescr_frspch].set(static_cast<float>(spinach_atoms) / static_cast<float>(matoms));

  int non_ring_non_spinach_atoms = 0;    // atoms between rings
  int heteroatoms_in_spinach = 0;
  int rotatable_bonds_in_spinach = 0;
  int rotatable_bonds_in_scaffold = 0;
  int bonds_in_spinach = 0;
  int saturated_spinach_atoms = 0;
  int unsaturated_spinach_atoms = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == spinach[i])
    {
      if (0 == ring_membership[i])
        non_ring_non_spinach_atoms++;
    }
    else if (6 != z[i])
      heteroatoms_in_spinach++;

    const Atom * ai = atom[i];

    int icon = ai->ncon();

    if (spinach[i])
    {
      if (icon == ai->nbonds())
        saturated_spinach_atoms++;
      else
        unsaturated_spinach_atoms++;
    }

    if (spinach[i])
    {
      for (int j = 0; j < icon; j++)
      {
        const Bond * b = ai->item(j);

        atom_number_t k = b->other(i);

        if (k < i || 0 == spinach[k])    // only count bonds once. Want both atoms in spinach
          continue;

        bonds_in_spinach++;
  
        if (b->is_single_bond() && ncon[i] > 1 && ncon[k] > 1)
          rotatable_bonds_in_spinach++;
      }
    }
    else   // atom I is NOT in spinach, in scaffold
    {
      for (int j = 0; j < icon; j++)
      {
        const Bond * b = ai->item(j);

        if (b->nrings())
          continue;

        atom_number_t k = b->other(i);

        if (b->is_single_bond() && ncon[i] > 1 && ncon[k] > 1 && 0 == ring_membership[i] && 0 == ring_membership[k])
          rotatable_bonds_in_scaffold++;
      }
    }
  }

  if (spinach_atoms > 0) {
    descriptor[iwdescr_spchtro].set(iwmisc::Fraction<float>(heteroatoms_in_spinach, spinach_atoms));
  }

  m.ring_membership();

  descriptor[iwdescr_nrnspch].set(static_cast<float>(non_ring_non_spinach_atoms));
  if (matoms == spinach_atoms)
    descriptor[iwdescr_fnrnspc].set(static_cast<float>(1.0));
  else
    descriptor[iwdescr_fnrnspc].set(static_cast<float>(non_ring_non_spinach_atoms) / static_cast<float>(matoms - spinach_atoms));

  if (bonds_in_spinach)
    descriptor[iwdescr_rbfrspch].set(static_cast<float>(rotatable_bonds_in_spinach) / static_cast<float>(bonds_in_spinach));

  descriptor[iwdescr_satspcha].set   (static_cast<float>(saturated_spinach_atoms));
  descriptor[iwdescr_unsatspcha].set (static_cast<float>(unsaturated_spinach_atoms));
  if (spinach_atoms > 0) {
    descriptor[iwdescr_fsatspcha].set(iwmisc::Fraction<float>(saturated_spinach_atoms, spinach_atoms));
  }

  // Work out the number of branches within the scaffold.
  // Examine the non-ring, scaffold atoms - the linker atoms.
  int branches_in_scaffold = 0;
  for (int i = 0; i < matoms; ++i) {
    if (spinach[i] || ring_membership[i]) {
      continue;
    }

    branches_in_scaffold += WithinScaffoldBranched(m, i);
  }

  descriptor[iwdescr_scaffoldbranches].set(static_cast<float>(branches_in_scaffold));

  int terminal_rings = 0;
  int internal_rings = 0;
  int spinach_connections = 0;
  int non_spinach_connections = 0;

  int nr = m.nrings();

//cerr << "Molecule has " << nr << " rings\n";

  resizable_array<int> strongly_fused;
  for (int i = 0; i < nr; ++i)
  {
    const Ring * ri = m.ringi(i);
    if (ri->strongly_fused_ring_neighbours())
      strongly_fused.add_if_not_already_present(ri->fused_system_identifier());
  }

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    if (strongly_fused.contains(ri->fused_system_identifier()))
      continue;

    int s, ns;     // spinach and non-spinach
    count_spinach_connections(m, *ri, spinach, s, ns);

//  cerr << "s " << s << " ns " << ns << '\n';

    if (1 == ns)
      terminal_rings++;
    else
      internal_rings++;

    spinach_connections += s;
    non_spinach_connections += ns;
  }

//cerr << "spinach_connections " << spinach_connections << " non_spinach_connections " << non_spinach_connections << " terminal_rings " << terminal_rings << '\n';

  descriptor[iwdescr_trmnlrng].set(static_cast<float>(terminal_rings));
  descriptor[iwdescr_intrnlrng].set(static_cast<float>(internal_rings));
  descriptor[iwdescr_rng2spch].set(static_cast<float>(spinach_connections));
  descriptor[iwdescr_rng2bridge].set(static_cast<float>(non_spinach_connections));

  return 1;
}

static int
compute_spinach_descriptors(Molecule & m,
               const atomic_number_t * z,
               const int * ncon,
               const int * ring_membership,
               const Atom ** atom)
{
  int matoms = m.natoms();

  int * spinach = new_int(matoms); std::unique_ptr<int[]> free_spinach(spinach);

  return compute_spinach_descriptors(m, matoms, spinach, z, ncon, ring_membership, atom);
}

static int
compute_rule_of_five_stuff(Molecule & m,
                           const atomic_number_t * z)
{
  int matoms = m.natoms();

  int ohnh = 0;
  int on = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 == z[i] || 8 == z[i])
    {
      on++;

      const int h = m.hcount(i);

      if (0 == h)    // acceptor
        continue;

      if (7 == z[i] && 2 == h)
        ohnh += 2;
      else
        ohnh++;
    }
  }

  descriptor[iwdescr_ro5_ohnh].set(static_cast<float>(ohnh));
  descriptor[iwdescr_ro5_on].set(static_cast<float>(on));

  return 1;
}


static double
identify_mr_halogen(const Molecule & m,
                    const atomic_number_t * z,
                    int * already_done)
{
  int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    if (9 == z[i])
    {
      rc += 1.0632;
      already_done[i] = 15;
    }
    else if (17 == z[i])
    {
      rc += 5.6105;
      already_done[i] = 16;
    }
    else if (35 == z[i])
    {
      rc += 8.6782;
      already_done[i] = 17;
    }
    else if (53 == z[i])
    {
      rc += 13.8741;
      already_done[i] = 18;
    }
  }

  return rc;
}

static double
identify_mr_oxygen(const Molecule & m,
                   const atomic_number_t * z,
                   const Atom * const * atom,
                   int * already_done)
{
  int matoms = m.natoms();

#ifdef DEBUG_IDENTIFY_MR_OXYGEN
  cerr << "On entry to identify_mr_oxygen\n";
  for (int i = 0; i < matoms; i++)
  {
    cerr << "atom " << i << " already_done " << already_done[i] << '\n';
  }
#endif

  double rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i] || already_done[i])
      continue;

    const Atom * ai = atom[i];

    if (2 == ai->ncon())
    {
      rc += 1.6351;
      already_done[i] = 8;
      continue;
    }

    if (0 == ai->ncon())
      continue;

    const Bond * b = ai->item(0);

    if (b->is_single_bond())
    {
      rc += 1.6351;
      already_done[i] = 7;
      continue;
    }

    atom_number_t n = b->other(i);

    if (7 == z[n])
    {
      rc += 2.1407;
      already_done[i] = 9;
    }
    else
    {
      rc += 1.7956;
      already_done[i] = 8;
    }
  }

  return rc;
}

static double
identify_mr_sulphur(const Molecule & m,
                    const atomic_number_t * z,
                    const Atom * const * atom,
                    int * already_done)
{
  double rc = 0.0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (16 != z[i] || already_done[i])
      continue;

    const Atom * ai = atom[i];

    int icon = ai->ncon();

    if (1 == icon)
    {
      if (1 == ai->nbonds())
      {
        rc += 7.3190;
        already_done[i] = 19;
      }
      else
      {
        rc += 9.1680;
        already_done[i] = 20;
      }

      continue;
    }

    if (icon == ai->nbonds())
    {
      rc += 7.3190;
      already_done[i] = 19;
      continue;
    }

    atom_number_t o1 = INVALID_ATOM_NUMBER;
    atom_number_t o2 = INVALID_ATOM_NUMBER;

    for (int j = 0; j < icon; j++)
    {
      const Bond * b = ai->item(j);

      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other(i);

      if (8 != z[k] || already_done[k])
        continue;

      if (INVALID_ATOM_NUMBER == o1)
        o1 = k;
      else
        o2 = k;
    }

    if (INVALID_ATOM_NUMBER != o2)
    {
      rc += 5.3321;
      already_done[i] = 22;
    }
    else if (INVALID_ATOM_NUMBER != o1)
    {
      rc += 6.0762;
      already_done[i] = 21;
    }
    else
    {
      rc += 9.1680;
      already_done[i] = 20;
    }
  }

  return rc;
}

/*
  I'm not really sure about Phosphorus, so I'll just guess
*/

static double
identify_mr_phosphorus(const Molecule & m,
                       const atomic_number_t * z,
                       const Atom * const * atom,
                       int * already_done)
{
  double rc = 0.0;

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    if (15 != z[i])
      continue;

    already_done[i] = 24;
    rc++;
  }

  return 5.3 * rc;
}

static double
identify_mr_nitro(const Molecule & m,
                  const atomic_number_t * z,
                  const Atom * const * atom,
                  int * already_done)
{
  int matoms = m.natoms();

  float rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    if (8 != z[i])
      continue;

    if (already_done[i])
      continue;

    const Atom * ai = atom[i];

    if (1 != ai->ncon())
      continue;

    const Bond * b = ai->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t n = b->other(i);

    if (7 != z[n] || already_done[n])
      continue;

    const Atom * an = atom[n];

    if (3 != an->ncon())
      continue;

    if (5 != an->nbonds())
      continue;

    atom_number_t other_oxygen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = an->item(j);

      if (! b->is_double_bond())
        continue;

      atom_number_t k = b->other(n);

      if (8 != z[k] || already_done[k])
        continue;

      other_oxygen = k;
      break;
    }

    if (INVALID_ATOM_NUMBER == other_oxygen)
      continue;

//  already_done[i] = 13;
    already_done[n] = 13;
//  already_done[other_oxygen] = 13;
    rc++;
  }

  return rc * 3.5054;
}

static double
identify_mr_nitrogen(Molecule & m,
                     const atomic_number_t * z,
                     const Atom * const * atom,
                     int * already_done)
{
  int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i] || already_done[i])
      continue;

    if (m.is_aromatic(i))
    {
      rc += 2.7662;
      already_done[i] = 12;
      continue;
    }

    const Atom * n = atom[i];

    int ncon = n->ncon();

    if (n->nbonds() == ncon)
    {
      rc += 3.0100;
      already_done[i] = 10;
      continue;
    }

    if (1 == ncon)
    {
      rc += 3.2009;
      already_done[i] = 11;
      continue;
    }

//  Check if doubly bonded to a heteroatom, or if part of =N-O group
    
    atom_number_t aromatic_neighbour = INVALID_ATOM_NUMBER;
    atom_number_t doubly_bonded_heteroatom = INVALID_ATOM_NUMBER;
    atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;
    atom_number_t singly_bonded_nitrogen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < ncon; j++)
    {
      const Bond * b = n->item(j);

      atom_number_t k = b->other(i);

      if (b->is_single_bond())
      {
        if (m.is_aromatic(k))
          aromatic_neighbour = k;
        else if (8 == z[k])
          singly_bonded_oxygen = k;
        else if (7 == z[k])
          singly_bonded_nitrogen = k;

        continue;
      }
      else if (b->is_double_bond() && 6 != z[k])
        doubly_bonded_heteroatom = k;
    }

    if (INVALID_ATOM_NUMBER != doubly_bonded_heteroatom && INVALID_ATOM_NUMBER != aromatic_neighbour)
    {
      already_done[i] = 14;
      rc += 3.8095;
    }
    else if (INVALID_ATOM_NUMBER != singly_bonded_oxygen)
    {
      already_done[i] = 23;
      rc += 4.2109;
    }
    else if (INVALID_ATOM_NUMBER != singly_bonded_nitrogen)
    {
      already_done[i] = 23;
      rc += 3.9809;
    }
    else
    {
      already_done[i] = 11;
      rc += 3.2009;
    }
  }

  return rc;
}


static double
identify_mr_carbon(Molecule & m,
                   const atomic_number_t * z,
                   const Atom * const * atom,
                   int * already_done)
{
  int matoms = m.natoms();

  double rc = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    if (6 != z[i] || already_done[i])
      continue;

    if (m.is_aromatic(i))
    {
      already_done[i] = 4;
      rc += 3.5090;
      continue;
    }

    const Atom * c = atom[i];

    int ncon = c->ncon();

    if (ncon == c->nbonds())
    {
      already_done[i] = 1;
      rc += 2.8158;
      continue;
    }

    atom_number_t doubly_bonded_heteroatom = INVALID_ATOM_NUMBER;
    int found_triple_bond = 0;

    for (int j = 0; j < ncon; j++)
    {
      const Bond * b = c->item(j);

      if (b->is_single_bond())
        continue;

      if (b->is_triple_bond())
      {
        found_triple_bond = 1;
        break;
      }

      atom_number_t k = b->other(i);

      if (6 != z[k])
        doubly_bonded_heteroatom = k;
    }

    if (INVALID_ATOM_NUMBER != doubly_bonded_heteroatom)
    {
      rc += 3.0887;
      already_done[i] = 5;
    }
    else if (found_triple_bond)
    {
      already_done[i] = 3;
      rc += 3.8974;
    }
    else
    {
      rc += 3.8278;
      already_done[i] = 2;
    }
  }

  return rc;
}

static int
compute_molar_refractivity(Molecule & m,
                    const atomic_number_t * z,
                    const Atom * const * atom,
                    const int * is_aromatic)
{
  int matoms = m.natoms();

  int * already_done = new_int(matoms); std::unique_ptr<int[]> free_already_done(already_done);

  double rc = 0.0;

  rc += identify_mr_halogen(m, z, already_done);
  rc += identify_mr_sulphur(m, z, atom, already_done);
  rc += identify_mr_phosphorus(m, z, atom, already_done);
  rc += identify_mr_nitro(m, z, atom, already_done);
  rc += identify_mr_oxygen(m, z, atom, already_done);
  rc += identify_mr_nitrogen(m, z, atom, already_done);
  rc += identify_mr_carbon(m, z, atom, already_done);

  for (int i = 0; i < matoms; i++)
  {
    rc += m.implicit_hydrogens(i) * 0.9155;
  }

  descriptor[iwdescr_cmr].set(rc / 10.0);

//#define DEBUG_COMPUTE_MOLAR_REFRACTIVITY
#ifdef DEBUG_COMPUTE_MOLAR_REFRACTIVITY
 
  int write_smiles = 0;
  for (int i = 0;i < matoms; i++)
  {
    if (0 == already_done[i])
    {
      if (report_unclassified_atoms)
        cerr << "Unclassified CMR atom, molecule '" << m.name() << "', atom " << i << " '" << m.smarts_equivalent_for_atom(i) << '\n';
      write_smiles = 1;
    }
  }


  if (write_smiles)
  {
    m.set_isotopes(already_done);
    cerr << m.smiles() << ' ' << m.name() << '\n';
    m.transform_to_non_isotopic_form();
  }

#endif

  return 1;
}


static float
identify_acm_n(const Molecule & m,
                int matoms,
                const atomic_number_t * z,
                const Atom * const * atom,
                int * already_done)
{
  int nplus = 0;
  int n = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (7 != z[i])
      continue;

    if (1 == atom[i]->formal_charge())
      nplus++;
    else
      n++;

    already_done[i] = 1;
  }

  return nplus * 11.5 + n * 1.2;
}

static float
identify_acm_c(const Molecule & m,
               int matoms,
               const atomic_number_t * z,
               const Atom * const * atom,
               int * already_done)
{
  int csp2 = 0;
  int csp3 = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (6 != z[i])
      continue;

    const Atom * ac = atom[i];

    if (ac->ncon() < ac->nbonds())
      csp2++;
    else
      csp3++;

    already_done[i] = 1;
  }

  return csp3 * 0.8 + csp2 * 0.7;
}

static float
identify_acm_halogen(const Molecule & m,
                 int matoms,
                 const atomic_number_t * z,
                 const Atom * const * atom,
                 int * already_done)
{
  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (9 == z[i])
      ;
    else if (17 == z[i])
      ;
    else if (35 == z[i])
      ;
    else if (53 == z[i])
      ;
    else
      continue;

    already_done[i] = 1;
    rc++;
  }

  return 1.3 * rc;
}

static float
identify_acm_os(const Molecule & m,
                int matoms,
                const atomic_number_t * z,
                const Atom * const * atom,
                int * already_done)
{
  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (8 == z[i] || 16 == z[i])
    {
      already_done[i] = 1;
      rc++;
    }
  }

  return 10.0 * rc;
}

static float
identify_acm_co(const Molecule & m,
                int matoms,
                const atomic_number_t * z,
                const Atom * const * atom,
                int * already_done)
{
  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (8 != z[i])
      continue;

    if (1 != atom[i]->ncon())
      continue;

    const Bond * b = atom[i]->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t c = b->other(i);

    if (already_done[c])   // strange
      continue;

    already_done[c] = 1;
    already_done[i] = 1;
    rc++;
  }

  return rc * 3.5;
}

static float
identify_acm_oh(const Molecule & m,
                int matoms,
                const atomic_number_t * z,
                const Atom * const * atom,
                int * already_done)
{
  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (8 != z[i])
      continue;

    if (1 != atom[i]->ncon())
      continue;

    if (1 == const_cast<Atom *>(atom[i])->implicit_hydrogens())
      rc++;
  }

  return rc * 2.5;
}

static float
identify_acm_acids(const Molecule & m,
                   int matoms,
                   const atomic_number_t * z,
                   const Atom * const * atom,
                   int * already_done)
{
  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    if (8 != z[i])
      continue;

    const Atom * ai = atom[i];

    if (1 != ai->ncon())
      continue;

    const Bond * b = ai->item(0);

    if (! b->is_double_bond())
      continue;

    atom_number_t c = b->other(i);

    if (already_done[c])
      continue;

    if (6 != z[c])
      continue;

    const Atom * ac = atom[c];

    if (3 != ac->ncon())
      continue;

    atom_number_t singly_bonded_oxygen = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 3; j++)
    {
      const Bond * b = ac->item(j);

      if (b->is_double_bond())
        continue;

      atom_number_t o = b->other(c);

      if (8 != z[o])
        continue;

      if (already_done[o])
        continue;

      if (1 != atom[o]->ncon())
        continue;

      singly_bonded_oxygen = o;
      break;
    }

    if (INVALID_ATOM_NUMBER == singly_bonded_oxygen)
      continue;

    already_done[i] = 1;
    already_done[c] = 1;
    already_done[singly_bonded_oxygen] = 1;
    rc++;
  }

  return rc * 8.2;
}

static float
identify_acm_phosphoric_acids(const Molecule & m,
                    int matoms,
                    const atomic_number_t * z,
                    const Atom * const * atom,
                    int * already_done)
{
  int rc = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (15 != z[i])
      continue;

    const Atom * ai = atom[i];

    if (4 != ai->ncon())
      continue;

    atom_number_t o1 = INVALID_ATOM_NUMBER;
    atom_number_t o2 = INVALID_ATOM_NUMBER;
    atom_number_t o3 = INVALID_ATOM_NUMBER;

    for (int j = 0; j < 4; j++)
    {
      atom_number_t k = ai->other(i, j);

      if (8 != z[k])
        continue;

      if (1 != atom[k]->ncon())
        continue;

      if (INVALID_ATOM_NUMBER == o1)
        o1 = k;
      else if (INVALID_ATOM_NUMBER == o2)
        o2 = k;
      else if (INVALID_ATOM_NUMBER == o3)
        o3 = k;
    }

    if (INVALID_ATOM_NUMBER == o3)
      continue;

    already_done[i] = 1;
    already_done[o1] = 1;
    already_done[o2] = 1;
    already_done[o3] = 1;

    rc++;
  }

  return rc * 10.0;
}

static int
andrews_craik_martin(Molecule & m,
                     const atomic_number_t * z,
                     const Atom * const * atom)
{
  int matoms = m.natoms();

  int * already_done = new_int(matoms); std::unique_ptr<int[]> free_already_done(already_done);

  float rc = static_cast<float>(0.0);

//cerr << "Using '" << m.smiles() << "'\n";

  rc += identify_acm_phosphoric_acids(m, matoms, z, atom, already_done);
//cerr << "After P, rc " << rc << '\n';
  rc += identify_acm_acids(m, matoms, z, atom, already_done);
//cerr << "After acids " << rc << '\n';
  rc += identify_acm_oh(m, matoms, z, atom, already_done);
//cerr << "After oh " << rc << '\n';
  rc += identify_acm_co(m, matoms, z, atom, already_done);
//cerr << "After co " << rc << '\n';
  rc += identify_acm_os(m, matoms, z, atom, already_done);
//cerr << "After os " << rc << '\n';
  rc += identify_acm_halogen(m, matoms, z, atom, already_done);
//cerr << "After halogen " << rc << '\n';
  rc += identify_acm_c(m, matoms, z, atom, already_done);
//cerr << "After c " << rc << '\n';
  rc += identify_acm_n(m, matoms, z, atom, already_done);
//cerr << "After n " << rc << '\n';

  descriptor[iwdescr_acmbe].set(rc);

  return 1;
}

/*
  The number of atoms which are attached to a halogen
*/

static int
compute_halogen_attachments(Molecule & m,
                            int matoms,
                            const atomic_number_t * z,
                            const int * ncon,
                            const Atom ** atom)
{
  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (9 == z[i] || 17 == z[i] || 35 == z[i] || 53 == z[i])
      continue;

    if (1 == ncon[i])
      continue;

    const Atom * ai = atom[i];
    for (int j = 0; j < ncon[i]; j++)
    {
      atom_number_t k = ai->other(i, j);
      if (9 == z[k] || 17 == z[k] || 35 == z[k] || 53 == z[k])
      {
        rc++;
        break;
      }
    }
  }

  descriptor[iwdescr_halogena].set(static_cast<float>(rc));

  return rc;
}

/*
  An electron rich section of a molecule is one in which adjacent atoms
  have pi electrons
*/

static int
compute_electron_rich(const Molecule & m,
                      atom_number_t a,
                      const Atom ** atom,
                      int * already_done)
{
  already_done[a] = 1;

  int rc = 1;

  const Atom * ai = atom[a];

  for (const Bond* b : *ai) {
    atom_number_t j = b->other(a);
    if (already_done[j])
      continue;

    Atom * aj = const_cast<Atom *>(atom[j]);

    boolean electron_rich = false;
    int lp = 0;

    if (aj->ncon() < aj->nbonds())
      electron_rich = true;
    else if (aj->lone_pair_count(lp) && lp)
     electron_rich = true;

    if (electron_rich)
      rc += compute_electron_rich(m, j, atom, already_done);
  }

  return rc;
}

/*
  A doubly bonded N has been identified. Is it a guanidine
*/
static int
guanidine(Molecule & m,
          atom_number_t doubly_bonded_nitrogen,
          int * nitrogens,
          const Atom ** atom,
          const atomic_number_t * z,
          const int * ncon,
          const int * ring_membership)
{
  atom_number_t c = atom[doubly_bonded_nitrogen]->other (doubly_bonded_nitrogen, 0);

  if (3 != ncon[c] || ring_membership[c])
    return 0;

  const Atom * carbon = atom[c];

  if (4 != carbon->nbonds())
    return 0;

// * --- [linker_atom] ---- C ( ==== doubly_bonded_nitrogen) --- [terminal_atom]

  atom_number_t terminal_atom = INVALID_ATOM_NUMBER;
  atom_number_t linker_atom   = INVALID_ATOM_NUMBER;

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = carbon->item(i);

    if (b->is_double_bond())     // must be the bond to doubly_bonded_nitrogen
      continue;

    atom_number_t j = b->other(c);

    if (7 != z[j] || ring_membership[j])
      return 0;

    if (1 == ncon[j] && INVALID_ATOM_NUMBER == terminal_atom)
      terminal_atom = j;
    else if (2 == ncon[j] && INVALID_ATOM_NUMBER == linker_atom)
      linker_atom = j;
    else
      return 0;
  }

  if (INVALID_ATOM_NUMBER == terminal_atom || INVALID_ATOM_NUMBER == linker_atom)
    return 0;
   
  nitrogens[doubly_bonded_nitrogen] = 1;
  nitrogens[terminal_atom] = 1;
  nitrogens[linker_atom] = 1;

  return 1;
}

static int
compute_symmetry_related_descriptors(Molecule & m) {
  const int * dm = m.distance_matrix_warning_may_change();

  const int matoms = m.natoms();

  const int * s = m.symmetry_classes();

  int asymmetric = 0;

  int furthest_separated_symmetry_related_atoms = 0;

  int max_equivalent_atoms = 0;

  int longest_path = 0;

  extending_resizable_array<int> class_done;

  for (int i = 0; i < matoms; ++i) {
    const int si = s[i];
    if (class_done[si] > 0) {
      continue;
    }

    class_done[si] = 1;

    // Counters for the atoms in this class.
    int ns = 1;
    int maxsep = 0;
    for (int j = 0; j < matoms; ++j) {
      if (i == j) {
        continue;
      }
      const int d = dm[i * matoms + j];
      if (d > longest_path) {
        longest_path = d;
      }
      if (si == s[j]) {
        ns++;
        if (d > maxsep)
          maxsep = d;
      }
    }

    if (1 == ns)         // not related to anything else
    {
      asymmetric++;
      continue;
    }

    if (ns > max_equivalent_atoms)
      max_equivalent_atoms = ns;

    if (maxsep > furthest_separated_symmetry_related_atoms)
      furthest_separated_symmetry_related_atoms = maxsep;
  }

  descriptor[iwdescr_symmatom].set(matoms - asymmetric);
  descriptor[iwdescr_fsymmatom].set(Fraction<float>(matoms - asymmetric, matoms));
  descriptor[iwdescr_lsepsymatom].set(furthest_separated_symmetry_related_atoms);
  if (longest_path > 0)
    descriptor[iwdescr_flsepsymatom].set(Fraction<float>(furthest_separated_symmetry_related_atoms, longest_path));
  else
    descriptor[iwdescr_flsepsymatom].set(0.0f);
  descriptor[iwdescr_maxsymmclass].set(max_equivalent_atoms);

  return 1;
}

/*
  Counts the number of 4 connected carbons that are in chains and rings
  asymc is the number of 4 connected atoms with 4 different atom types connected.
*/

static void
compute_4_connected_carbon_stuff(const Molecule & m,
                                 const atomic_number_t * z,
                                 const int * ncon,
                                 const int * ring_membership)
{
  int cd4_ring = 0;
  int cd4_chain = 0;

  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (4 != ncon[i])
      continue;

    if (6 != z[i])
      continue;

    if (ring_membership[i])
      cd4_ring++;
    else
      cd4_chain++;
  }

  descriptor[iwdescr_cd4ring].set(static_cast<float>(cd4_ring));
  descriptor[iwdescr_cd4chain].set(static_cast<float>(cd4_chain));

  return;
}

/*
  Is atom ZATOM a nitro group

  Thought I was going to need this for extended conjugation
*/

/*static int
is_nitro (Molecule & m,
          atom_number_t zatom,
          const Atom ** atom)
{
  const Atom * a = atom[zatom];

  int doubly_bonded_oxygens = 0;

  int acon = a->ncon();
  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item (i);

    if (! b->is_double_bond())
      continue;

    atom_number_t j = b->other (zatom);

    const Atom * aj = atom[j];

    if (8 == aj->atomic_number())
      doubly_bonded_oxygens++;
  }

  return 2 == doubly_bonded_oxygens;
}*/

/*
  Do two rings share 3 or more consecutive atoms in common
*/

static int
strongly_fused(const Ring & r1,
               const Ring & r2,
               int matoms,
               int * tmp)
{

  set_vector(tmp, matoms, 0);
  r1.set_vector(tmp, 1);
  r2.increment_vector(tmp);

// If there are only two atoms in common, then these cannot be planar fused

  if (2 == count_occurrences_of_item_in_array(2, matoms, tmp))
    return 0;

// Well, there could be a case where two rings had 3 atoms in common but not
// two bonds. Imagine the case of 2 rings that are fused somewhere, and then
// join somewhere else in a spiro fusion. Let's assume that's too rare to worry about

  if (0 == max_difference_in_ring_size_for_strongly_fused)  // no checking sizes
    ;
  else if (r1.number_elements() - r2.number_elements() > max_difference_in_ring_size_for_strongly_fused)
    return 0;
  else if (r2.number_elements() - r1.number_elements() > max_difference_in_ring_size_for_strongly_fused)
    return 0;

  return 1;
}

/*
  Turns out that the code for growing a strongly or weakly fused ring system
  is very similar. This function can do either one
*/

static int
grow_fused_system(Molecule & m,
                  int i,
                  int * ring_already_done,
                  int flag,
                  int * atmp,
                  int growing_strongly_fused_system)
{
  const Ring * ri = m.ringi(i);

  int max_bonds_shared = ri->largest_number_of_bonds_shared_with_another_ring();
   
  ring_already_done[i] = flag;

  int rc = 1;

  int matoms = m.natoms();

  int number_neighbours = ri->fused_ring_neighbours();

  for (int j = 0; j < number_neighbours; j++)
  {
    const Ring * rj = ri->fused_neighbour(j);

    int k = rj->ring_number();

//  Decide whether we are growing a strong or weakly fused system

    if (1 == max_bonds_shared)    // growing a weakly fused system
    {
      if (rj->largest_number_of_bonds_shared_with_another_ring() > 1)
        continue;
    }
    else if (rj->largest_number_of_bonds_shared_with_another_ring() < 2)   // pair of rings not strongly fused
      continue;

    if (ring_already_done[k])
      continue;

    int sf = strongly_fused(*ri, *rj, matoms, atmp);

    if (growing_strongly_fused_system && sf)
      ;
    else if (! growing_strongly_fused_system && ! sf)
      ;
    else
      continue;

    rc += grow_fused_system(m, k, ring_already_done, flag, atmp, growing_strongly_fused_system);
  }

  return rc;
}

static int
compute_planar_fused_rings(Molecule & m,
                           const atomic_number_t * z,
                           int * ring_already_done,
                           int * atmp)
{
  int matoms = m.natoms();

  int nr = m.nrings();

  int number_planar_fused_systems_found = 0;
  int rings_in_planar_fused_systems = 0;
  int largest_planar_fused_system_size = 0;
  int max_heteroatoms_in_system = 0;
  int heterocycles_in_planar_fused_systems = 0;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    if (! ri->is_fused())
      continue;

    if (ri->largest_number_of_bonds_shared_with_another_ring() > 1)    // in a strongly fused system
      continue;

    int flag = matoms + 2 * i;     // just a unique number that won't collide with a fused system identifier

    int planar_system_size = grow_fused_system(m, i, ring_already_done, flag, atmp, 0);

    set_vector(atmp, matoms, 0);
    ri->set_vector(atmp, 1);

    number_planar_fused_systems_found++;
    rings_in_planar_fused_systems += planar_system_size;

    if (planar_system_size > largest_planar_fused_system_size)
      largest_planar_fused_system_size = planar_system_size;

    if (heteroatoms_in_ring(ri, z))
      heterocycles_in_planar_fused_systems++;

    for (int j = i + 1; j < nr; j++)
    {
      if (flag != ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);
      rj->set_vector(atmp, 1);

      int h = heteroatoms_in_ring(m.ringi(j), z);

//    cerr << "For ring " << (*rj) << " h = " << h << '\n';

      if (h)
        heterocycles_in_planar_fused_systems++;
    }

    int heteroatoms_in_system = 0;
    for (int j = 0; j < matoms; j++)
    {
      if (0 == atmp[j])    // not in system
        continue;

      if (6 != z[j])
        heteroatoms_in_system++;
    }

    if (heteroatoms_in_system > max_heteroatoms_in_system)
      max_heteroatoms_in_system = heteroatoms_in_system;
  }

  descriptor[iwdescr_npfsdsys].set(static_cast<float>(number_planar_fused_systems_found));
  descriptor[iwdescr_rnginpfs].set(static_cast<float>(rings_in_planar_fused_systems));
  descriptor[iwdescr_lgplnfsy].set(static_cast<float>(largest_planar_fused_system_size));
  descriptor[iwdescr_htrcpfsy].set(static_cast<float>(heterocycles_in_planar_fused_systems));
  descriptor[iwdescr_mxhtpfsy].set(static_cast<float>(max_heteroatoms_in_system));

  return 1;
}

/*
  Just check to see that there is at least one ring that is strongly fused and satisfies the
  size constraint
*/

static int
joined_rings_too_different_in_size(Molecule & m,
                                   const Ring & r,
                                   int * tmp)
{
  int matoms = m.natoms();

  for (int i = 0; i < r.fused_ring_neighbours(); i++)
  {
    const Ring * rj = r.fused_neighbour(i);

    if (strongly_fused(r, *rj, matoms, tmp))
      return 0;
  }

  return 1;
}

/*
*/

static int
compute_non_planar_fused_rings(Molecule & m,
                               const atomic_number_t * z,
                               int * ring_already_done,
                               int * atmp)
{
  int matoms = m.natoms();

  int nr = m.nrings();

  int number_strongly_fused_systems_found = 0;
  int rings_in_strongly_fused_systems = 0;
  int largest_strongly_fused_system_size = 0;     // the number of rings in it
  int max_heteroatoms_in_system = 0;
  int heterocycles_in_strongly_fused_systems = 0;

  for (int i = 0; i < nr; i++)
  {
    if (ring_already_done[i])
      continue;

    const Ring * ri = m.ringi(i);

    if (! ri->is_fused())
      continue;

    if (ri->largest_number_of_bonds_shared_with_another_ring() < 2)
      continue;

    if (0 == max_difference_in_ring_size_for_strongly_fused)
      ;
    else if (joined_rings_too_different_in_size(m, *ri, atmp))
      continue;

    int flag = ri->fused_system_identifier() + 1;    // a unique identifier

    int strongly_fused_size = grow_fused_system(m, i, ring_already_done, flag, atmp, 1);

    number_strongly_fused_systems_found++;
    rings_in_strongly_fused_systems += strongly_fused_size;

    if (strongly_fused_size > largest_strongly_fused_system_size)
      largest_strongly_fused_system_size = strongly_fused_size;

    set_vector(atmp, matoms, 0);

    for (int j = 0; j < nr; j++)
    {
      if (flag != ring_already_done[j])
        continue;

      const Ring * rj = m.ringi(j);
      rj->set_vector(atmp, 1);

      int h = heteroatoms_in_ring(rj, z);

      if (h)
        heterocycles_in_strongly_fused_systems++;
    }

    int heteroatoms_in_system = 0;
    for (int j = 0; j < matoms; j++)
    {
      if (atmp[j] && 6 != z[j])
        heteroatoms_in_system++;
    }

    if (heteroatoms_in_system > max_heteroatoms_in_system)
      max_heteroatoms_in_system = heteroatoms_in_system;
  }

  descriptor[iwdescr_nsfsdsys].set(static_cast<float>(number_strongly_fused_systems_found));
  descriptor[iwdescr_rnginsfs].set(static_cast<float>(rings_in_strongly_fused_systems));
  descriptor[iwdescr_lgstrfsy].set(static_cast<float>(largest_strongly_fused_system_size));
  descriptor[iwdescr_htrcsfsy].set(static_cast<float>(heterocycles_in_strongly_fused_systems));
  descriptor[iwdescr_mxhtsfsy].set(static_cast<float>(max_heteroatoms_in_system));

  return 1;
}

static int
compute_fused_rings(Molecule & m,
                    const atomic_number_t * z)
{
  int nr = m.nrings();

  if (nr < 2)
    return 1;

  int * rtmp = new_int(nr); std::unique_ptr<int[]> free_rtmp(rtmp);

  int * atmp = new int[m.natoms()]; std::unique_ptr<int[]> free_atmp(atmp);

  (void) compute_non_planar_fused_rings(m, z, rtmp, atmp);

  set_vector(rtmp, nr, 0);

  return compute_planar_fused_rings(m, z, rtmp, atmp);
}

template <typename T>
void
back_to_zero(T * v,
             int n,
             T vfrom, 
             T vto)
{
  for (int i = 0; i < n; ++i)
  {
    if (vfrom == v[i])
      v[i] = vto;
  }

  return;
}

static int
discern_conjugated_section(Molecule & m,
                           atom_number_t zatom,
                           const Atom ** atom,
                           const int * is_aromatic,
                           int * already_done,
                           int flag)
{
  already_done[zatom] = flag;

  int rc = 1;

  const Atom * a = atom[zatom];

  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);

    if (already_done[j])
      continue;

    int pe;

    if (b->is_double_bond() || b->is_triple_bond())    // will definitely extend across B
      ;
    if (const_cast<Atom *>(atom[j])->pi_electrons(pe) && pe)
      ;
    else
      continue;

#ifdef DEBUG_COMPUTE_EXTENDED_CONJUGATION
    cerr << "From atom " << zatom << " to " << j << " amide? " << is_amide <<'\n';
#endif

#ifdef DEBUG_COMPUTE_EXTENDED_CONJUGATION
    cerr << "From atom " << zatom << " extend to atom " << j << '\n';
#endif
    rc += discern_conjugated_section(m, j, atom, is_aromatic, already_done, flag);
  }

  return rc;
}

/*
  this turned out to be really hard, so we get a very simplified version.
  An extended conjugation section is just a set of contiguous atoms
  all of which have pi electrons
*/

static int
compute_extended_conjugation(Molecule & m,
                             const Atom ** atom,
                             const int * is_aromatic,
                             int * already_done)
{
  int matoms = m.natoms();

  int number_conjugated_sections = 0;
  int max_conjugated_section_size = 0;
  int atoms_in_conjugated_sections = 0;

  m.ring_membership();

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i])
      continue;

    int pe;
    if (! const_cast<Atom *>(atom[i])->pi_electrons(pe) || 0 == pe)
      continue;

    int conjugated_section_size = discern_conjugated_section(m, i, atom, is_aromatic, already_done, number_conjugated_sections + 1);

#ifdef DEBUG_COMPUTE_EXTENDED_CONJUGATION
    cerr << "In '" << m.name() << "' found section of size " << conjugated_section_size << ", ndx " << (number_conjugated_sections+1) << " starting with atom " << i << '\n';
#endif

    if (conjugated_section_size <= 2)    // Just an isolated double bond, or maybe a section that didn't extend
    {
      back_to_zero(already_done, matoms, number_conjugated_sections+1, 0);
      continue;
    }

    atoms_in_conjugated_sections += conjugated_section_size;
    number_conjugated_sections++;

    if (conjugated_section_size > max_conjugated_section_size)
      max_conjugated_section_size = conjugated_section_size;
  }

  int carbons_in_conjugated_sections = 0;

  if (number_conjugated_sections)
  {
    for (int i = 0; i < matoms; i++)
    {
      if (already_done[i] > 0 && 6 == atom[i]->atomic_number())
        carbons_in_conjugated_sections++;
    }
  }

//#define TEST_EXTENDED_CONJUGATED_SECTIONS
#ifdef TEST_EXTENDED_CONJUGATED_SECTIONS
  Molecule mcopy = m;

  for (int i = 0; i < matoms; i++)
  {
    if (already_done[i] > 0)
      mcopy.set_isotope(i, already_done[i]);
  }

  cerr << mcopy.smiles() << ' ' << m.name() << " extended " << number_conjugated_sections << ", " << atoms_in_conjugated_sections << " atoms\n";

#endif

  descriptor[iwdescr_nconjgsc].set(static_cast<float>(number_conjugated_sections));
  descriptor[iwdescr_atincnjs].set(static_cast<float>(atoms_in_conjugated_sections));
  descriptor[iwdescr_mxcnjscz].set(static_cast<float>(max_conjugated_section_size));
  descriptor[iwdescr_cinconjs].set(static_cast<float>(carbons_in_conjugated_sections));

//m.transform_to_non_isotopic_form ();

  return 1;
}

static int
compute_extended_conjugation(Molecule & m,
                             const Atom ** atom,
                             const int * is_aromatic)
{
  int * already_done = new_int(m.natoms()); std::unique_ptr<int[]> free_already_done(already_done);

  return compute_extended_conjugation(m, atom, is_aromatic, already_done);
}

static int
compute_double_bond_substitution(Molecule & m,
                                 const atomic_number_t * z,
                                 const int * ncon)
{
  int nb = m.nedges();

  int number_doubly_bonded_carbons = 0;

  int total_double_bond_substituents = 0;

  for (int i = 0; i < nb; i++)
  {
    const Bond * b = m.bondi(i);

    if (! b->is_double_bond())
      continue;

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (m.in_same_aromatic_ring(a1, a2))
      continue;

    if (6 == z[a1])
    {
      number_doubly_bonded_carbons++;
      total_double_bond_substituents += ncon[a1] - 1;
    }

    if (6 == z[a2])
    {
      number_doubly_bonded_carbons++;
      total_double_bond_substituents += ncon[a2] - 1;
    }
  }

  descriptor[iwdescr_numcdb].set(static_cast<float>(number_doubly_bonded_carbons));
  descriptor[iwdescr_totdbsub].set(static_cast<float>(total_double_bond_substituents));

  if (number_doubly_bonded_carbons)
    descriptor[iwdescr_avcdbsub].set(static_cast<float>(total_double_bond_substituents) / static_cast<float>(number_doubly_bonded_carbons));

  return 1;
}
static int
compute_amine_count(Molecule & m,
                    const Atom ** atom,
                    const atomic_number_t * z,
                    const int * ncon,
                    const int * ring_membership,
                    int * already_done)
{
// First perceive any guanadines  CNC(=N)N

  int amine_count = 0;
  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    if (7 == z[i] && 1 == ncon[i] && 2 == atom[i]->nbonds())
      amine_count += guanidine(m, i, already_done, atom, z, ncon, ring_membership);
  }

// Now process all remaining N's

  for (int i = 0; i < matoms; i++)
  {
    if (7 == z[i] && 0 == already_done[i] && m.hcount(i))
      amine_count++;
  }

  descriptor[iwdescr_amine].set(static_cast<float>(amine_count));

  return 1;
}

static int
compute_pyridine_pyrrole(Molecule & m,
                         const Atom ** atom,
                         const atomic_number_t * z,
                         const int * ncon,
                         const int * ring_membership)
{
  int matoms = m.natoms();

  int pyridine = 0;
  int pyrrole = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (7 != z[i])
      continue;

    if (0 == ring_membership[i])
      continue;

    if (! m.is_aromatic(i))
      continue;

    if (3 == ncon[i])   // not interested in these
      continue;

    if (0 == m.implicit_hydrogens(i))
      pyridine++;
    else
      pyrrole++;
  }

  descriptor[iwdescr_pyridine].set(static_cast<float>(pyridine));
  descriptor[iwdescr_pyrrole].set(static_cast<float>(pyrrole));

  return 1;
}

static int
carboxylic_acid(const Atom ** atom, atom_number_t doubly_bonded_oxygen,
               const atomic_number_t * z,
               const int * ncon)
{
  atom_number_t c = atom[doubly_bonded_oxygen]->other(doubly_bonded_oxygen, 0);

  if (3 != ncon[c] || 6 != z[c])
    return 0;

  const Atom * ac = atom[c];

  int found_singly_bonded_oxygen = 0;

  for (int i = 0; i < 3; i++)
  {
    const Bond * b = ac->item(i);

    if (! b->is_single_bond())
      continue;

    atom_number_t j = b->other(c);

    if (8 == z[j] && 1 == ncon[j])
      found_singly_bonded_oxygen++;
  }

  return (1 == found_singly_bonded_oxygen);
}

static int
compute_topological_descriptors(Molecule & m,
               const atomic_number_t * z,
               const int * ncon,
               const int * ring_membership,
               const Atom ** atom,
               const int * is_aromatic,
               int * already_done)
{
  int matoms = m.natoms();

  int heavy_atom_count = 0;
  int two_connected_chain_atom = 0;
  int rotatable_bonds  = 0;
  int ring_atom_count = 0;
  int heteroatom_count = 0;
  int ring_heteroatom_count = 0;
  int fused_ring_atom_count = 0;
  int singly_connected_oxygen_or_sulphur_count = 0;
  int carboxylic_acid_count = 0;
  int aromatic_atom_count = 0;
  int aromatic_heteroatom_count = 0;
  int chain_multiple_bonds = 0;
  int ring_multiple_bonds = 0;
  int carbon_hydrogen_count = 0;
  int ch2 = 0;
  int unsaturation = 0;       // non aromatic
  int atoms_with_pi_electrons = 0;      // aromatic, unsaturated and lone pair electrons - this last type aren't really pi electrons...
  int electron_rich_sections = 0;       // contiguous parts of the molecule with lots of electrons
  int atoms_in_electron_rich_sections = 0;   // the number of atoms in electron rich sections
  int largest_electron_rich_section = 0;
  int halogen_count = 0;
  int bigatom_count = 0;
  int aliphatic_carbon_count = 0;
  int non_ring_non_halogen_heteroatoms = 0;
  int csp3 = 0;
  int csp3_chain = 0;
  int hcount = 0;

// We count the number of atoms with a given connectivity

  int connected[6];
  set_vector(connected, 6, 0);

  int total_connectivity = 0;
  int total_aliphatic_connectivity = 0;
  int total_chain_connectivity = 0;

  m.ring_membership();    // force SSSR if needed

  for (int i = 0; i < matoms; i++)
  {
    total_connectivity += ncon[i];

    if (ncon[i] <= 4)
      connected[ncon[i]]++;
    else
      connected[5]++;

    atomic_number_t zi = z[i];

    if (1 == zi)
      hcount++;
    else
    {
      heavy_atom_count++;
      hcount += m.hcount(i);
    }

    if (6 != zi)
      heteroatom_count++;

    if (zi > 10)
      bigatom_count++;

    if ((8 == zi || 16 == zi) && 1 == ncon[i])
      singly_connected_oxygen_or_sulphur_count++;

    Atom * ai = const_cast<Atom *>(atom[i]);

    if (8 == zi && 1 == ncon[i] && 2 == ai->nbonds())
      carboxylic_acid_count += carboxylic_acid(atom, i, z, ncon);

    if (6 == zi && ncon[i] == ai->nbonds())
    {
      csp3++;
      if (0 == ring_membership[i])
        csp3_chain++;
    }

    int is_halogen;
    if (9 == zi || 17 == zi || 35 == zi || 53 == zi)
    {
      is_halogen = 1;
      halogen_count++;
    }
    else
      is_halogen = 0;

    if (6 != z[i] && 0 == ring_membership[i] && ! is_halogen)
      non_ring_non_halogen_heteroatoms++;

    boolean aromatic = false;
    if (ring_membership[i])
    {
      ring_atom_count++;

      if (is_aromatic[i])
      {
        aromatic_atom_count++;
        aromatic = true;
      }

      if (6 != zi)
      {
        ring_heteroatom_count++;
        if (aromatic)
          aromatic_heteroatom_count++;
      }
      if (ring_membership[i] > 1)
        fused_ring_atom_count++;
    }
    else
    {
      total_chain_connectivity += ncon[i];

      if (2 == ncon[i])
        two_connected_chain_atom++;
    }

    if (6 == zi && ! aromatic)
      aliphatic_carbon_count++;
    if (6 == zi && ai->implicit_hydrogens())
    {
      carbon_hydrogen_count++;
      if (2 == ai->implicit_hydrogens())
        ch2++;
    }

//  Determine electron rich sections. Skip if this atom was already included
//  in an electron rich contiguous section

    int lp = 0;
    boolean electron_rich = false;

    if (aromatic)
    {
      atoms_with_pi_electrons++;
      electron_rich = true;
    }
    else
    {
      total_aliphatic_connectivity += ncon[i];

      if (ncon[i] < ai->nbonds())
      {
        atoms_with_pi_electrons++;
        unsaturation++;
        electron_rich = true;
      }
      else if (ai->lone_pair_count(lp) && lp)
      {
        atoms_with_pi_electrons++;
        electron_rich = true;
      }
    }

    if (electron_rich && 0 == already_done[i])
    {
      electron_rich_sections++;
      int tmp = compute_electron_rich(m, i, atom, already_done);
      if (tmp > largest_electron_rich_section)
        largest_electron_rich_section = tmp;
      atoms_in_electron_rich_sections += tmp;
//    cerr << electron_rich_sections << " electon rich sections, na = " << atoms_in_electron_rich_sections << '\n';
    }

    int icon = ncon[i];
    for (int j = 0; j < icon; j++)
    {
      const Bond * b = ai->item(j);

      atom_number_t k = b->other(i);

      if (k < i)        // don't want to count bonds twice.
        continue;

      if (b->nrings())
      {
        if (b->is_aromatic())      // Sept 2013, used to count these, but final result was too correlated with aromatic atom count
          ;
        else if (! b->is_single_bond())
          ring_multiple_bonds++;
        continue;
      }

      if (b->is_single_bond() && ncon[i] > 1 && ncon[k] > 1 && ! triple_bond_at_either_end(m, b) && ! part_of_otherwise_non_rotabable_entity(m, i, k))
        rotatable_bonds++;
      else if (! b->is_single_bond())
        chain_multiple_bonds++;
    }
  }

  descriptor[iwdescr_natoms].set(static_cast<float>(matoms));

  Molecular_Weight_Calculation_Result mwcr;
  if (m.molecular_weight(mwc, mwcr))
  {
//  cerr << "AMW " << mwcr.amw() << ' ' << m.name() << '\n';
    descriptor[iwdescr_amw].set(mwcr.amw());
//  float q;
//  descriptor[iwdescr_amw].value(q);
//  cerr << "Value being used " << q << " hydrogens " << m.implicit_hydrogens() << '\n';
  } else {
    descriptor[iwdescr_amw].set(0.0f);
  }

  descriptor[iwdescr_bigatom].set(static_cast<float>(bigatom_count));
  descriptor[iwdescr_fbigatom].set(static_cast<float>(bigatom_count) / static_cast<float>(matoms));
  descriptor[iwdescr_halogen].set(static_cast<float>(halogen_count));
  descriptor[iwdescr_nelem].set(static_cast<float>(m.number_different_elements()));

  descriptor[iwdescr_hcount].set(static_cast<float>(hcount));
  descriptor[iwdescr_hperatom].set(static_cast<float>(hcount) / static_cast<float>(matoms));

  descriptor[iwdescr_csp3].set(static_cast<float>(csp3));
  descriptor[iwdescr_fcsp3].set(static_cast<float>(csp3)/static_cast<float>(matoms));

  if (heteroatom_count == matoms)
    descriptor[iwdescr_fccsp3].set(static_cast<float>(0.0));   // hard to know what this should be
  else
    descriptor[iwdescr_fccsp3].set(static_cast<float>(csp3)/static_cast<float>(matoms - heteroatom_count));

  descriptor[iwdescr_csp3_chain].set(static_cast<float>(csp3_chain));

  if (descriptors_to_compute.ncon_descriptors) {
    const int offset = iwdescr_ncon1;    // possibly dangerous cast from enumeration to int

    for (int i = 1; i < 5; i++)
    {
      descriptor[offset + 2 * (i - 1)].set(static_cast<float>(connected[i]));
      descriptor[offset + 2 * (i - 1) + 1].set(static_cast<float>(connected[i]) / static_cast<float>(matoms));
    }
  }

  int highly_connected = 0;
  for (int i = 2; i < 5; i++) {
    highly_connected += connected[i];
  }

  descriptor[iwdescr_frhc].set(static_cast<float>(highly_connected) / static_cast<float>(matoms));

  int mbonds = m.nedges();
  if (0 == mbonds)           // let's hope no-one does computations on single atom molecules
  {
    cerr << "Warning, no bonds in molecule\n";
    mbonds = 1;
  }

  descriptor[iwdescr_mltbd].set(static_cast<float>(chain_multiple_bonds + ring_multiple_bonds));
  descriptor[iwdescr_fmltbd].set(static_cast<float>(chain_multiple_bonds + ring_multiple_bonds) / static_cast<float>(mbonds));
  descriptor[iwdescr_chmltbd].set(static_cast<float>(chain_multiple_bonds));
  descriptor[iwdescr_fchmltbd].set(static_cast<float>(chain_multiple_bonds) / static_cast<float>(mbonds));
  descriptor[iwdescr_rgmltbd].set(static_cast<float>(ring_multiple_bonds));
  descriptor[iwdescr_frgmltbd].set(static_cast<float>(ring_multiple_bonds) / static_cast<float>(mbonds));

//int chain_atoms = matoms - ring_atom_count;
//output << "frmcc " << static_cast<float>(singly_connected_atom) / static_cast<float>(matoms)  << '\n';
//output << "mcc " << singly_connected_atom << '\n';

  descriptor[iwdescr_dcca].set(static_cast<float>(two_connected_chain_atom));
  descriptor[iwdescr_fdcca].set(static_cast<float>(two_connected_chain_atom) / static_cast<float>(matoms));

  descriptor[iwdescr_rotbond].set(static_cast<float>(rotatable_bonds));
  descriptor[iwdescr_frotbond].set(iwmisc::Fraction<float>(rotatable_bonds, m.nedges()));
  descriptor[iwdescr_ringatom].set(static_cast<float>(ring_atom_count));
  if (ring_atom_count)
  {
    descriptor[iwdescr_nrings].set(static_cast<float>(m.nedges() - matoms + 1));
    descriptor[iwdescr_nnsssrng].set(static_cast<float>(m.non_sssr_rings()));
    descriptor[iwdescr_rhacnt].set(static_cast<float>(ring_heteroatom_count));
    descriptor[iwdescr_rhaf].set(static_cast<float>(ring_heteroatom_count) / static_cast<float>(ring_atom_count));

//   The fraction of rings atoms which are fused;

    descriptor[iwdescr_frafus].set(static_cast<float>(fused_ring_atom_count) / static_cast<float>(ring_atom_count));
  }
  else     // some of the rings descriptors are truly zero's, others are missing or zero by definition
  {
    descriptor[iwdescr_nrings].set(static_cast<float>(0.0));
    descriptor[iwdescr_nnsssrng].set(static_cast<float>(0.0));
    descriptor[iwdescr_rhacnt].set(static_cast<float>(0.0));
  }

// The fraction of ring atoms.

  if (heavy_atom_count > 0)
    descriptor[iwdescr_rngatmf].set(static_cast<float>(ring_atom_count) / static_cast<float>(heavy_atom_count));
  else
    descriptor[iwdescr_rngatmf].set(0.0f);

  descriptor[iwdescr_aroma].set(static_cast<float>(aromatic_atom_count));
  descriptor[iwdescr_aromha].set(static_cast<float>(aromatic_heteroatom_count));
  descriptor[iwdescr_aromc].set(static_cast<float>(aromatic_atom_count - aromatic_heteroatom_count));
  descriptor[iwdescr_aliphc].set(static_cast<float>(aliphatic_carbon_count));
  if (ring_atom_count)
  {
    descriptor[iwdescr_fraromha].set(static_cast<float>(aromatic_heteroatom_count) / static_cast<float>(ring_atom_count));

    descriptor[iwdescr_aromdens].set(static_cast<float>(aromatic_atom_count) / static_cast<float>(matoms));
  }
  else
  {
    descriptor[iwdescr_aromdens].set(static_cast<float>(0.0));
  }

  descriptor[iwdescr_ch2].set(static_cast<float>(ch2));

  descriptor[iwdescr_ch].set(static_cast<float>(carbon_hydrogen_count));

// The fraction of heavy atoms which are heteroatoms

  descriptor[iwdescr_htroatom].set(static_cast<float>(heteroatom_count));
  if (heavy_atom_count > 0)
    descriptor[iwdescr_htroaf].set(static_cast<float>(heteroatom_count) / static_cast<float>(heavy_atom_count));
  else
    descriptor[iwdescr_htroaf].set(0.0f);

  descriptor[iwdescr_nrgnhlht].set(static_cast<float>(non_ring_non_halogen_heteroatoms));

  descriptor[iwdescr_ohsh].set(static_cast<float>(singly_connected_oxygen_or_sulphur_count));

  descriptor[iwdescr_co2h].set(static_cast<float>(carboxylic_acid_count));

  descriptor[iwdescr_atmpiele].set(atoms_with_pi_electrons);
  descriptor[iwdescr_fratmpie].set(static_cast<float>(atoms_with_pi_electrons) / static_cast<float>(matoms));

  descriptor[iwdescr_unsatura].set(static_cast<float>(unsaturation));
  descriptor[iwdescr_funsatura].set(static_cast<float>(unsaturation)/static_cast<float>(matoms));
  descriptor[iwdescr_erichsct].set(static_cast<float>(electron_rich_sections));
  descriptor[iwdescr_aiercsct].set(static_cast<float>(atoms_in_electron_rich_sections));
  descriptor[iwdescr_lercsct].set(static_cast<float>(largest_electron_rich_section));
  descriptor[iwdescr_faiercst].set(static_cast<float>(atoms_in_electron_rich_sections) / static_cast<float>(matoms));

  descriptor[iwdescr_platt].set(static_cast<float>(total_connectivity));
  descriptor[iwdescr_avcon].set(static_cast<float>(total_connectivity) / static_cast<float>(matoms));
  if (matoms - aromatic_atom_count > 0) {
    descriptor[iwdescr_avalcon].set(static_cast<float>(total_aliphatic_connectivity) / static_cast<float>(matoms - aromatic_atom_count));
  }

  if (matoms - ring_atom_count > 0) {
    descriptor[iwdescr_avchcon].set(static_cast<float>(total_chain_connectivity) / static_cast<float>(matoms - ring_atom_count));
  }

  if (halogen_count <= 1)
    descriptor[iwdescr_halogena].set(static_cast<float>(halogen_count));
  else
    compute_halogen_attachments(m, matoms, z, ncon, atom);

  // Harder to individualise, TODO...
  compute_radha_entropy_descriptors(m, ncon, ring_membership, atom);

  set_vector(already_done, matoms, 0);

  if (descriptors_to_compute.specific_groups) {
    compute_amine_count(m, atom, z, ncon, ring_membership, already_done);
    compute_pyridine_pyrrole(m, atom, z, ncon, ring_membership);
  }

  if (descriptors_to_compute.simple_hbond_descriptors) {
    compute_hbond_descriptors(m, z, ncon, atom);
  }

  compute_ring_descriptors(m, z, ncon);

  if (descriptors_to_compute.complexity_descriptors) {
    do_compute_spiro_fusions(m, ncon, ring_membership);
  }

  if (descriptors_to_compute.ramey_descriptors) {
    compute_ramey_descriptors(m, z);
  }

  if (descriptors_to_compute.compute_xlogp) {
    compute_xlogp(m);
  }
  if (descriptors_to_compute.compute_alogp) {
    compute_alogp(m);
  }


  if (descriptors_to_compute.ring_substitution_descriptors) {
    compute_ring_substitution_descriptors(m, ncon);
  }

  if (descriptors_to_compute.adjacent_ring_fusion_descriptors) {
    compute_adjacent_ring_fusion_descriptors(m, ncon, ring_membership);
  }

  if (descriptors_to_compute.ring_substitution_ratio_descriptors) {
    compute_ring_substitution_ratios(m, ncon, atom, ring_membership);
  }

  compute_4_connected_carbon_stuff(m, z, ncon, ring_membership);

  if (descriptors_to_compute.bonds_between_rings) {
    compute_between_ring_descriptors(m, ncon, atom);
  }

  if (descriptors_to_compute.ring_fusion_descriptors) {
    compute_ring_fusion_descriptors(m, ncon, ring_membership);
  }

  if (descriptors_to_compute.polar_bond_descriptors) {
    compute_polar_bond_descriptors(m, z, ncon);
  }

  compute_double_bond_substitution(m, z, ncon);

  compute_extended_conjugation(m, atom, is_aromatic);

  if (descriptors_to_compute.partial_symmetry_descriptors) {
    NearSymmetricDescriptors(m);
  }

  return 1;
}

static int
compute_distance_charge_descriptors(Molecule & m, 
                                    IWString & output)
{
  const int maxdist = 5;
  int matoms = m.natoms();

  for (int dist = 1; dist < maxdist; dist++)
  {
    charge_t max_charge_separation = 0.0;
    charge_t sum_charge_sepatation = 0.0;
    int number_atoms_this_separation = 0;
    for (int i = 0; i < matoms; i++)
    {
      for (int j = i + 1; j < matoms; j++)
      {
        if (dist != m.bonds_between(i, j))
          continue;
        number_atoms_this_separation++;
        charge_t diff = fabs(m.charge_on_atom(i) - m.charge_on_atom(j));
        if (diff > max_charge_separation)
          max_charge_separation = diff;
        sum_charge_sepatation += diff;
      }
    }

    if (number_atoms_this_separation)
    {
      output << "mxqdif" << dist << " " << max_charge_separation << '\n';
      output << "avqdif" << dist << " " << sum_charge_sepatation / static_cast<float>(number_atoms_this_separation) << '\n';
    }
    else
      break;
  }

  return 0;
}

static int
compute_distance_4_charge_heteroatom_descriptors(Molecule & mol, 
               IWString & output,
               const atomic_number_t * z,
               const int * ncon,
               atom_number_t * conn)
{

  int matoms = mol.natoms();
  int count = 0;
  double sum_charge_separation = 0.0;
  double max_charge_separation = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t zi = z[i];
    if (6 == zi || 1 == zi)
      continue;

    int icon = ncon[i];
    for (int j = 0; j < icon; j++)
    {
      atom_number_t k = mol.other(i, j);
      if (k < i)
        continue;

//    Atoms i and k now define the middle of the dihedral
//    Now find atoms bonded to i and k

      int kcon = mol.connections(k, conn);
      for (int l = 0; l < kcon; l++)
      {
        atom_number_t m = conn[l];
        if (m == i)
          continue;

        charge_t mcharge = mol.charge_on_atom(m);

        for (int n = 0; n < icon; n++)
        {
          atom_number_t o = mol.other(i, n);
          if (o == k)
            continue;

          double diff = fabs(mcharge - mol.charge_on_atom(o));
          count++;
          sum_charge_separation += diff;
          if (diff > max_charge_separation)
            max_charge_separation = diff; 
        }
      }
    }
  }


  if (count)
  {
    output << "hmxqdif3 " << max_charge_separation << '\n';
    output << "havqdif3 " << sum_charge_separation / count << '\n';
  }

  return 1;
}

static int
compute_distance_3_charge_heteroatom_descriptors(Molecule & mol, 
               IWString & output,
               const atomic_number_t * z,
               const int * ncon,
               atom_number_t * conn)
{

  (void) ncon;

  int matoms = mol.natoms();
  int count = 0;
  double sum_charge_separation = 0.0;
  double max_charge_separation = 0.0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t zi = z[i];
    if (6 == zi || 1 == zi)
      continue;

    int icon = mol.connections(i, conn);
    for (int j = 0; j < icon; j++)
    {
      atom_number_t k = conn[j];
      charge_t kcharge = mol.charge_on_atom(k);
      for (int l = j + 1; l < icon; l++)
      {
        atom_number_t m = conn[l];
        double diff = fabs(kcharge - mol.charge_on_atom(m));
        count++;
        sum_charge_separation += diff;
        if (diff > max_charge_separation)
          max_charge_separation = diff; 
//      Coordinates dipole;
//      compute_dipole(m, i, k, m, dipole);
      }
    }
  }


  if (count)
  {
    output << "hmxqdif2 " << max_charge_separation << '\n';
    output << "havqdif2 " << sum_charge_separation / count << '\n';
  }

  return 1;
}

static int
compute_distance_2_charge_heteroatom_descriptors(Molecule & m, 
               IWString & output,
               const atomic_number_t * z,
               atom_number_t * conn)
{

  int matoms = m.natoms();

  double sum_charge_separation = 0.0;
  double max_charge_separation   = 0.0;
  int count = 0;

  for (int i = 0; i < matoms; i++)
  {
    atomic_number_t zi = z[i];
    if (6 == zi || 1 == zi)
      continue;

    charge_t icharge = m.charge_on_atom(i);

    int icon = m.connections(i, conn);
    for (int j = 0; j < icon; j++)
    {
      atom_number_t k = conn[j];
      charge_t diff = fabs(icharge - m.charge_on_atom(k));
      sum_charge_separation += diff;
      if (diff > max_charge_separation)
        max_charge_separation = diff;
      count++;
    }
  }

  if (count)
  {
    output << "hmxqdif1 " << max_charge_separation << '\n';
    output << "havqdif1 " << sum_charge_separation / count << '\n';
  }

  return 1;
}

static int
compute_distance_charge_heteroatom_descriptors(Molecule & m, 
               IWString & output,
               const atomic_number_t * z,
               const int * ncon,
               atom_number_t * conn)
{
  compute_distance_2_charge_heteroatom_descriptors(m, output, z, conn);
  compute_distance_3_charge_heteroatom_descriptors(m, output, z, ncon, conn);
  compute_distance_4_charge_heteroatom_descriptors(m, output, z, ncon, conn);

  return 1;
}

static const charge_t LARGE_POSITIVE = 19.0;
static const charge_t LARGE_NEGATIVE = -19.0;

static int
compute_charge_descriptors(Molecule & m, 
               IWString & output,
               const atomic_number_t * z)
{
  charge_t sum_charge = 0.0;
  charge_t sum_negative_charges = 0.0;
  charge_t sum_positive_charges = 0.0;

  charge_t sum_charge_heavy_atoms = 0.0;
  charge_t max_charge_heavy_atom  = 0.0;
  int      heavy_atom_count = 0;

  charge_t sum_charge_heteroatoms = 0.0;
  charge_t max_charge_heteroatom  = 0.0;
  int      heteroatom_count = 0;

  charge_t sum_charge_ring_atoms = 0.0;
  charge_t max_charge_ring_atom  = 0.0;
  int      ring_atom_count = 0;

  charge_t sum_charge_ring_heteroatoms = 0.0;
  charge_t max_charge_ring_heteroatom  = 0.0;
  int      ring_heteroatom_count = 0;

  int fused_ring_atom_count = 0;

  charge_t most_positive = LARGE_NEGATIVE;
  atom_number_t most_positive_atom = INVALID_ATOM_NUMBER;
  charge_t most_negative = LARGE_POSITIVE;
  atom_number_t most_negative_atom = INVALID_ATOM_NUMBER;

  charge_t most_positive_ring = LARGE_NEGATIVE;
  atom_number_t most_positive_ring_atom = INVALID_ATOM_NUMBER;

  charge_t most_negative_ring = LARGE_POSITIVE;
  atom_number_t most_negative_ring_atom = INVALID_ATOM_NUMBER;

  charge_t most_positive_non_ring = LARGE_NEGATIVE;
  atom_number_t most_positive_non_ring_atom = INVALID_ATOM_NUMBER;

  charge_t most_negative_non_ring = LARGE_POSITIVE;
  atom_number_t most_negative_non_ring_atom = INVALID_ATOM_NUMBER;

  int matoms = m.natoms();
  for (int i = 0; i < matoms; i++)
  {
    charge_t q = m.charge_on_atom(i);
    double fabsq = fabs(q);

    sum_charge += fabsq;
    if (q < 0.0)
      sum_negative_charges += q;
    else if (q > 0.0)
      sum_positive_charges += q;

    if (q < most_negative)
    {
      most_negative = q;
      most_negative_atom = i;
    }
    if (q > most_positive)
    {
      most_positive = q;
      most_positive_atom = i;
    }

    int zi = z[i];
    if (1 == zi)    // not interested in hydrogens beyond here
      continue;

    if (m.is_non_ring_atom(i))
    {
      if (q > most_positive_non_ring)
      {
        most_positive_non_ring_atom = i;
        most_positive_non_ring = q;
      }
      if (q < most_negative_non_ring)
      {
        most_negative_non_ring = q;
        most_negative_non_ring_atom = i;
      }
    }
    else
    {
      if (q > most_positive_ring)
      {
        most_positive_ring = q;
        most_positive_ring_atom = i;
      }
      if (q < most_negative_ring)
      {
        most_negative_ring = q;
        most_negative_ring_atom = i;
      }
    }

    heavy_atom_count++;
    sum_charge_heavy_atoms += q;
    if (fabsq > max_charge_heavy_atom)
      max_charge_heavy_atom = fabsq;

    if (6 != zi)
    {
      heteroatom_count++;
      sum_charge_heteroatoms += fabsq;
      if (fabsq > max_charge_heteroatom)
        max_charge_heteroatom = fabsq;
    }
    if (m.is_ring_atom(i))
    {
      ring_atom_count++;
      sum_charge_ring_atoms += fabsq;
      if (m.nrings(i) > 1)
        fused_ring_atom_count++;
      if (fabsq > max_charge_ring_atom)
        max_charge_ring_atom = fabsq;
      if (6 != zi)
      {
        ring_heteroatom_count++;
        sum_charge_ring_heteroatoms += fabsq;
        if (fabsq > max_charge_ring_heteroatom)
          max_charge_ring_heteroatom = fabsq;
      }
    }
  }

// The average charge in the molecule
// If there is more negative charge in the molecule than positive, we reverse its
// sign.

  if ((- sum_negative_charges) > sum_positive_charges)
    sum_charge = -sum_charge;

  output << "avch " << sum_charge / static_cast<float>(matoms) << '\n';
  output << "tnch " << sum_negative_charges << '\n';
  output << "tpch " << sum_positive_charges << '\n';

// The average charge on heavy atoms.

  output << "avchhvy " << sum_charge_heavy_atoms / static_cast<float>(heavy_atom_count) << '\n';

// Average charge on heteroatoms

  if (heteroatom_count)
    output << "avchhet " << sum_charge_heteroatoms / static_cast<float>(heteroatom_count) << '\n';
  else
    output << "avchhet " << undefined_value << '\n';

// most positive and most negative

  output << "mpc " << most_positive << '\n';
  output << "mnc " << most_negative << '\n';

  if (m.ok_2_atoms(most_positive_atom, most_negative_atom))
  {
    int tmp = m.bonds_between(most_positive_atom, most_negative_atom);
    output << "dist " << tmp << '\n';
    output << "pd " << (most_positive - most_negative) / static_cast<float>(tmp * tmp) << '\n';
  }

// Most positive and most negative non ring atoms

  if (INVALID_ATOM_NUMBER != most_negative_non_ring_atom)
    output << "mnnr " << most_negative_non_ring << '\n';
  if (INVALID_ATOM_NUMBER != most_positive_non_ring_atom)
    output << "mpnr " << most_positive_non_ring << '\n';

// The distance between them

  if (m.ok_2_atoms(most_negative_non_ring_atom, most_positive_non_ring_atom))
    output << "mnnrmpnr " << m.bonds_between(most_negative_non_ring_atom, most_positive_non_ring_atom) << '\n';

// All the various ring things are computed only if rings are present.

  if (ring_atom_count)
  {
//  Average charge on ring atoms

    output << "avchrng " << sum_charge_ring_atoms / static_cast<float>(ring_atom_count) << '\n';

//  Most positive and most negative ring atoms

    if (INVALID_ATOM_NUMBER != most_positive_ring_atom)
      output << "mpr "  << most_positive_ring << '\n';
    if (INVALID_ATOM_NUMBER != most_negative_ring_atom)
      output << "mnr "  << most_negative_ring << '\n';

//  Various distances

    if (m.ok_2_atoms(most_negative_ring_atom, most_positive_ring_atom))
      output << "mnrmpr " << m.bonds_between(most_negative_ring_atom, most_positive_ring_atom) << '\n';
    if (m.ok_2_atoms(most_negative_ring_atom, most_positive_non_ring_atom))
      output << "mnrmpnr " << m.bonds_between(most_negative_ring_atom, most_positive_non_ring_atom) << '\n';
    if (m.ok_2_atoms(most_negative_ring_atom, most_negative_non_ring_atom))
      output << "mnrmnnr " << m.bonds_between(most_negative_ring_atom, most_negative_non_ring_atom) << '\n';
    if (m.ok_2_atoms(most_positive_ring_atom, most_negative_ring_atom))
      output << "mprmnr "  << m.bonds_between(most_positive_ring_atom, most_negative_ring_atom) << '\n';
    if (m.ok_2_atoms(most_positive_ring_atom, most_positive_non_ring_atom))
      output << "mprmpnr " << m.bonds_between(most_positive_ring_atom, most_positive_non_ring_atom) << '\n';
    if (m.ok_2_atoms(most_positive_ring_atom, most_negative_non_ring_atom))
      output << "mprmnnr " << m.bonds_between(most_positive_ring_atom, most_negative_non_ring_atom) << '\n';

//  Average charge on ring heteroatoms

    if (ring_heteroatom_count)
      output << "avchrnghet " << sum_charge_ring_heteroatoms / static_cast<float>(ring_atom_count) << '\n';
  }

  return 0;
}

static void
store_charge_assigner_results(const Molecule & m,
                              const atomic_number_t * z,
                              const Atom ** atom)
{
  int npos = 0;
  int nneg = 0;

  int matoms = m.natoms();

  int positive_nitrogen = 0;
  int negative_nitrogen = 0;

  for (int i = 0; i < matoms; i++)
  {
    formal_charge_t fc = atom[i]->formal_charge();
    if (0 == fc)
      continue;

    if (fc < 0)
    {
      nneg++;
      if (7 == z[i])
        negative_nitrogen++;
    }
    else
    {
      npos++;
      if (7 == z[i])
        positive_nitrogen++;
    }
  }

  descriptor[iwdescr_brunsneg].set(static_cast<float>(nneg));
  descriptor[iwdescr_brunspos].set(static_cast<float>(npos));
  descriptor[iwdescr_formal_charge].set(static_cast<int>(nneg + npos));
  descriptor[iwdescr_nplus].set   (static_cast<float>(positive_nitrogen));
  descriptor[iwdescr_nminus].set  (static_cast<float>(negative_nitrogen));

  return;
}

static void
store_donor_acceptor_assigner_results(int matoms, const int * isotope)
{
  int acc = 0;
  int dual = 0;
  int don = 0;

  for (int i = 0; i < matoms; i++)
  {
    if (0 == isotope[i])
      ;
    else if (1 == isotope[i])
      acc++;
    else if (2 == isotope[i])
      dual++;
    else if (3 == isotope[i])
      don++;
  }

  acc += dual;
  don += dual;

  descriptor[iwdescr_brunsacc].set(static_cast<float>(acc));
  descriptor[iwdescr_brnsdual].set(static_cast<float>(dual));
  descriptor[iwdescr_brunsdon].set(static_cast<float>(don));
  descriptor[iwdescr_brunshbdsum].set(static_cast<float>(acc + don - dual));

  return;
}

static int
iwdescriptors(Molecule & m, 
              IWString_and_File_Descriptor & output, 
              const atomic_number_t * z,
              int * ncon,  // not const because of charge assigner.
              const int * ring_membership,
              atom_number_t * conn,
              const Atom ** atom,
              const int * is_aromatic,
              const char output_separator)
{
  assert (m.ok());
  assert (output.good());

  const int matoms = m.natoms();

  do_compute_chirality_descriptors(m, z, ncon);

  // Once chirality descriptors are computed, remove them. This will cause tests to
  // fail, so tests can only really be done on a-chiral molecules.
  m.remove_all_chiral_centres();

  int * already_done = new_int(matoms); std::unique_ptr<int[]> free_already_done(already_done);
  compute_topological_descriptors(m, z, ncon, ring_membership, atom, is_aromatic, already_done);

  if (descriptors_to_compute.crowding_descriptors) {
    compute_crowding_descriptors(m, ncon);
  }

  if (descriptors_to_compute.spinach_descriptors) {
    compute_spinach_descriptors(m, z, ncon, ring_membership, atom);
  }

  if (descriptors_to_compute.ring_chain_descriptors) {
    compute_ring_chain_descriptors(m, z, ncon, atom, ring_membership);
  }

  if (descriptors_to_compute.complexity_descriptors) {
    compute_fused_rings(m, z);
  }

// It looks as if the Novartis PSA descriptor should be done before charges are assigned

  if (descriptors_to_compute.psa) {
    descriptor[iwdescr_nvrtspsa].set(static_cast<float>(novartis_polar_surface_area(m, z, atom, is_aromatic)));
  }

#ifdef MCGOWAN
  if (descriptors_to_compute.mcgowan) {
    descriptor[iwdescr_mcgowan].set(McGowanVolume(m));
  }
#endif

  compute_molar_refractivity(m, z, atom, is_aromatic);

  compute_rule_of_five_stuff(m, z);

  if (descriptors_to_compute.distance_matrix_descriptors) {
    do_compute_distance_matrix_descriptors(m, z, ncon);
  }

  if (descriptors_to_compute.symmetry_descriptors) {
    compute_symmetry_related_descriptors(m);
  }

  if (descriptors_to_compute.charge_descriptors && charge_assigner.active())
  {
    (void) charge_assigner.process(m);
    store_charge_assigner_results(m, z, atom);
    andrews_craik_martin(m, z, atom);
    for (int i = 0; i < matoms; ++i) {
      ncon[i] = m.ncon(i);
    }
  }

  if (descriptors_to_compute.donor_acceptor && donor_acceptor_assigner.active())
  {
    int * da_results = new int[matoms]; std::unique_ptr<int[]> free_da_results(da_results);

    donor_acceptor_assigner.process(m, da_results);
    store_donor_acceptor_assigner_results(matoms, da_results);

    if (min_hbond_feature_separation > 0) {
      compute_bond_separation_parameters(m, da_results);
    }
  }

  // Not implemented, we are not dealing with partial charges here.
  if (0 && m.has_charges())
  {
    compute_charge_descriptors(m, output, z);
    compute_distance_charge_descriptors(m, output);
    compute_distance_charge_heteroatom_descriptors(m, output, z, ncon, conn);
    if (charge_queries.number_elements())
      compute_substructure_based_charge_descriptors(m, output);
  }

  return write_the_output(m, output_separator, output);
}

static int
iwdescriptors(Molecule & m,
              const char output_separator,
              IWString_and_File_Descriptor & output)
{
  for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
  {
    descriptor[i].reset();
  }

// For efficiency, we put many atom properties into arrays

  int matoms = m.natoms();

  atomic_number_t * z = new atomic_number_t[matoms]; std::unique_ptr<atomic_number_t[]> free_z(z);
  m.atomic_numbers(z);

  int * nc = new int[matoms]; std::unique_ptr<int[]> free_nc(nc);
  int maxcon = m.ncon(nc);

  int * ring_membership = new int[matoms]; std::unique_ptr<int[]> free_ring_membership(ring_membership);
//m.ring_membership(ring_membership);
  m.ring_membership_including_non_sssr_rings(ring_membership);

  atom_number_t * conn = new atom_number_t[maxcon]; std::unique_ptr<int[]> free_conn(conn);

  const Atom ** atom = new const Atom *[matoms]; std::unique_ptr<const Atom *[]> free_atom(atom);
  m.atoms(atom);

  int * is_aromatic = new int[matoms]; std::unique_ptr<int[]> free_is_aromatic(is_aromatic);

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < matoms; i++) {
    is_aromatic[i] = m.is_aromatic(i);
  }

  return iwdescriptors(m, output, z, nc, ring_membership, conn, atom, is_aromatic, output_separator);
}

/*
  Two values have been reported as different. Are they close enough?
*/

static int
values_are_equivalent(const Set_or_Unset<float> & s1,
                       const Set_or_Unset<float> & s2)
{
  float v1, v2;

  if (! s1.value(v1) || ! s2.value(v2))    // they should both be set
    return 0;

  if (fabs(v1 - v2) < 1.0e-05)
    return 1;

  return 0;
}

static int
results_are_different (const IWString & mname)
{
  assert (nullptr != descriptor);
  assert (nullptr != saved_result);

  int rc = 0;

  for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
  {
    if (! descriptor[i].active())
      continue;

    if (descriptor[i] == saved_result[i])
      continue;

//  If they are close enough, that's fine

    if (values_are_equivalent(descriptor[i], saved_result[i]))
      continue;

    if (0 == rc)
      cerr << "Mismatch '" << mname << "'\n";

    cerr << "Result mismatch, descriptor " << i << " '" << descriptor[i].descriptor_name() << "'\n";
    float s{};
    (void) saved_result[i].value(s);
    float d{};
    (void) descriptor[i].value(d);

    cerr << saved_result[i] << " vs " << descriptor[i] << " difference " << (s - d) << '\n';
    rc++;
  }

  return rc;
}

/*
  Doing tests is complicated by the fact that the underlying functions change
  the molecule. 
*/

static int
test_iwdescriptors(const IWString * rsmi, const IWString & mname)
{
  assert (ntest > 0);

// Save the results

  for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
  {
    if (! descriptor[i].active())
      continue;

    saved_result[i] = descriptor[i];
  }

  IWString_and_File_Descriptor notused;

  for (int i = 0; i < ntest; i++)
  {
    Molecule tmp;
    if (! tmp.build_from_smiles(rsmi[i])) {
      cerr << "Yipes, cannot parse random smiles '" << rsmi[i] << "'\n";
      return 0;
    }

    tmp.set_name(mname);

    // no output is performed when ntest > 0
    if (! iwdescriptors(tmp, ' ', notused)) {
      return 0;
    }

    if (results_are_different(mname)) {
      cerr << "Yipes, one or more result mismatches\n";
      cerr << "Random '" << rsmi[i] << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
test_iwdescriptors(Molecule & m)
{
  IWString * rsmi = new IWString[ntest]; std::unique_ptr<IWString[]> free_rsmi(rsmi);
  cerr << "Parent smiles " << m.smiles() << '\n';
  for (int i = 0; i < ntest; i++) {
    rsmi[i] = m.random_smiles();
    cerr << rsmi[i] << ' ' << m.name() << '\n';
  }

  int rc = test_iwdescriptors(rsmi, m.name());

  if (0 == rc) {
    cerr << "One or more test failures '" << m.name() << "', molecule " << molecules_read << '\n';
    cerr << m.smiles() << '\n';
    test_failures++;
  }

  if (keep_going_after_test_failure)
    return 1;

  return rc;
}

static void
preprocess(Molecule & m)
{
  if (1 == m.natoms() && 1 == m.atomic_number(0)) {   // we want to be able to compute a result for Hydrogen "molecule"
    return;
  }

  // Revert any hydrogen isotopes. Various other functions may be
  // reluctant to remove them otherwise.
  // Neutralize any formal charge on a Hydrogen atom.
  // https://dot-jira.lilly.com/browse/GC3TK-652?jql=assignee%20in%20(RX87690)

  int explicit_hydrogen_found = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    if (1 != a->atomic_number())
      continue;

    ++explicit_hydrogen_found;

    if (0 != a->isotope())
      m.set_isotope(i, 0);

    if (a->formal_charge()) {
      m.set_formal_charge(i, 0);
    }
  }

  if (chemical_standardisation.active())
    (void) chemical_standardisation.process(m);

  if (m.number_fragments() > 1)
  {
    if (0 == reduce_to_largest_fragment)
    {
      cerr << "Fatal, '" << m.name() << " has " << m.number_fragments() << " components\n";
      //iwabort();
    }

    (void) m.reduce_to_largest_fragment_carefully();
  }

  // The charge assigner calls move_hydrogens_to_end_of_connection_table,
  // which messes up cached atomic properties. Make that call now to
  // hopefully make the subsequent call a no-op.
  if (explicit_hydrogen_found) {
    m.MoveToEndOfConnectionTable(1);
  }

  const IWString & mname = m.name();

  if (0 == mname.length())
  {
    IWString newname;
    newname << "IWD" << molecules_read;

    cerr << "No name, set to '" << newname << "'\n";

    m.set_name(newname);
  }
  else if (number_filters)
    ;
  else
  {
    IWString newname(mname);
    newname.truncate_at_first(' ');

    m.set_name(newname);
  }

  return;
}

static int
handle_zero_atom_molecule(const Molecule & m,
                          IWString_and_File_Descriptor & output)
{
  if (verbose)
    cerr << "No atoms in " << m.name() << '\n';

  if (ignore_molecules_with_no_atoms)
    return 1;

  if (write_descriptor_file_pipeline) {
    output << m.name();
  } else {
    append_first_token_of_name(m.name(), output);
  }

  if (include_smiles_as_descriptor) {
    output << ' ' << smiles_as_input;
  }

  for (int i = 0; i < NUMBER_DESCRIPTORS; ++i)
  {
    output << " 0";
  }

  output << "\n";

  return 1;
}

class Alarm_Tripped
{
  private:
  public:
    Alarm_Tripped() {};
};

static void
alarm_handler(int sig)
{
  throw Alarm_Tripped();
}

static int
iwdescriptors_alarm(data_source_and_type<Molecule> & input,
                    const char output_separator,
                    IWString_and_File_Descriptor & output)
{
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);

    if (number_filters || include_smiles_as_descriptor) {
      smiles_as_input = m->smiles();
    }

    molecules_read++;

    preprocess(*m);

    if (m->empty()) {
      handle_zero_atom_molecule(*m, output);
      continue;
    }

    alarm(alarm_time);
    signal(SIGALRM, alarm_handler);
    
    int initial_size = output.length();

    int rc = 0;
    try{
      rc = iwdescriptors(*m, output_separator, output);
    }
    catch (Alarm_Tripped  at)
    {
      output.resize_keep_storage(initial_size);
      molecules_skipped_by_timer++;
      if (verbose)
        cerr << "Skipped '" << m->name() << "' by timer " << alarm_time << '\n';

      continue;
    }
    catch (...)
    {
      cerr << "Caught generic exception\n";
    }

    MaybeFlush(output);

    if (rc && ntest) {
      rc = test_iwdescriptors(*m);
    }

    if (0 == rc) {
      return 0;
    }
  }

  return 1;
}

static int
iwdescriptors(data_source_and_type<Molecule> & input,
              const char output_separator,
              IWString_and_File_Descriptor & output)
{
  assert (input.ok());

  Molecule * m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (number_filters || include_smiles_as_descriptor)
      smiles_as_input = m->smiles();

    molecules_read++;

    preprocess(*m);

    if (m->empty()) {
      handle_zero_atom_molecule(*m, output);
      continue;
    }

    if (! iwdescriptors(*m, output_separator, output)) {
      return 0;
    }

    MaybeFlush(output);

    if (ntest) {
      if (! test_iwdescriptors(*m)) {
        cerr << "Tests fail " << m->name() << '\n';
        return 0;
      }
    }
  }

  return 1;
}

static int
iwdescriptors_filter(iwstring_data_source & input,
                     const char output_separator,
                     IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    output << buffer << '\n';

    if (! buffer.starts_with(smiles_tag))
      continue;

    buffer.chop();
    buffer.remove_up_to_first('<');

    Molecule m;

    if (! m.build_from_smiles(buffer))
    {
      cerr << "Cannot parse smiles '" << buffer << "'\n";
      return 0;
    }

    if (! iwdescriptors(m, output_separator, output))
      return 0;
  }

  return 1;
}

static int
iwdescriptors_filter(const char * fname,
                     const char output_separator,
                     IWString_and_File_Descriptor & output)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open filter input stream\n";
    return 0;
  }

  return iwdescriptors_filter(input, output_separator, output);
}

static int
iwdescriptors_rpipe_line(const const_IWSubstring& buffer,
                    char output_separator,
                    IWString_and_File_Descriptor& output) {
  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "iwdescriptors_rpipe_line:invalid smiles\n";
    return 0;
  }

  return iwdescriptors(m, output_separator, output);
}

static int
iwdescriptors_rpipe(iwstring_data_source& input,
                    char output_separator,
                    IWString_and_File_Descriptor& output) {
  static int first_call = 1;
  const_IWSubstring buffer;

  if (first_call) {
    if (! input.next_record(buffer)) {
      cerr << "iwdescriptors_rpipe:cannot read header\n";
      return 0;
    }

    if (write_descriptor_file_pipeline) {
      output << buffer << output_separator;
      WriteHeader(descriptor, name_translation, output_separator, output);
    } else {
      buffer.remove_leading_words(1);
      output << buffer << output_separator;
      WriteHeader(descriptor, name_translation, output_separator, output);
    }

    first_call = 0;
  }

  while (input.next_record(buffer)) {
    if (! iwdescriptors_rpipe_line(buffer, output_separator, output)) {
      cerr << "Fatal error\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

static int
iwdescriptors_rpipe(const char* fname,
                    char output_separator,
                    IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "iwdescriptors_rpipe:cannot open '" << fname << "'\n";
    return 0;
  }

  return iwdescriptors_rpipe(input, output_separator, output);
}

static int
iwdescriptors(const char * fname,
              FileType input_type,
              const char output_separator,
              IWString_and_File_Descriptor & output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 1;
  }

  if (verbose > 1)
    input.set_verbose(1);

  if (alarm_time > 0)
    return iwdescriptors_alarm(input, output_separator, output);
  else
    return iwdescriptors(input, output_separator, output);
}

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
// clang-format off
  cerr << "  -u <char>      specify char to be used for undefined values\n";
  cerr << "  -q <file>      specify substructure query file\n";
  cerr << "  -b <dist>      extra (more expensive) Hydrogen bonding descriptors\n";
  cerr << "  -N <qualifier> charge assigner options, enter '-N help' for info\n";
  cerr << "  -H <qualifier> donor/accepter  options, enter '-H help' for info\n";
  cerr << "  -T <ntest>     make <ntest> random permutations for testing\n";
  cerr << "  -T kg          keep going after a test failure\n";
  cerr << "  -O ...         control of features to generate, enter '-O help'\n";
  cerr << "  -a <seconds>   abandon any computation after <seconds> seconds\n";
  cerr << "  -F <filter>    filter specification\n";
  cerr << "  -G ...         fingerprint specification\n";    
  cerr << "  -B ...         other optional features\n";
  cerr << "  -s .           include the smiles as a descriptor - for now . means first descriptor\n";
  cerr << "  -S             zero value for SD2 in Polar Surface Area computations\n";
  cerr << "  -d             write normally integer descriptors as floats\n";
  cerr << "  -z             write empty descriptors for zero atom molecules\n";
  cerr << "  -i <type>      specify input file type\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -E <symbol>    standard element related options\n";
  display_standard_chemical_standardisation_options(cerr, 'g');
  display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
parse_replicates_specification(const_IWSubstring & dname,
                               int & replicates)
{
  if (! dname.contains(':'))
    return 1;

  const_IWSubstring tmp;
  if (! dname.split_into_directive_and_value(tmp, ':', replicates))
    return 0;

  dname.truncate_at_first(':');

  return 1;
}

static void
DisplayDashBOptions(std::ostream& output) {
// clang-format off
  output << " -B prefix=...   prefix for descriptor names\n";
  output << " -B sep=.        output token separator (def space)\n";
  output << " -B mxsfsdf=nn   max difference in ring size for strongly fused rings\n";
  output << " -B quiet        do not report unclassified atoms\n";
  output << " -B namexref=<fname> file containing w.NameTranslation proto for name translations\n";
  output << " -B ranges=<fname>  file containing w.Ranges text proto with descriptor ranges\n";
  output << " -B rpipe        input is a smiles pipeline 'smiles id descriptors...'\n";
  output << " -B wpipe        write a smiles pipeline 'smiles id descriptors...'\n";
  output << " -B flush        flush output after each molecule\n";
// clang-format on
}

static int
ReadNameTranslation(const const_IWSubstring& fname,
                    IW_STL_Hash_Map_String& name_translation) {
  IWString tmp(fname);
  std::optional<w::Features> maybe_proto = 
        iwmisc::ReadTextProto<w::Features>(tmp);
  if (! maybe_proto) {
    cerr << "ReadNameTranslation:cannot read '" << fname << "'\n";
    return 0;
  }

  for (const auto& feature : (*maybe_proto).feature()) {
    const IWString old_name = feature.computed_name();
    const IWString new_name = feature.name();
    name_translation[old_name] = new_name;
  }

  return name_translation.size();
}
static int
ReadDescriptorRanges(const const_IWSubstring& fname) {
  IWString tmp(fname);
  std::optional<w::Ranges> maybe_proto = 
        iwmisc::ReadTextProto<w::Ranges>(tmp);
  if (! maybe_proto) {
    cerr << "ReadDescriptorRanges:cannot read '" << fname << "'\n";
    return 0;
  }

  // Ignore missing descriptors, those may not be turned on. Too complicated
  // otherwise. But what this means is that truly bad input will be silently
  // ignored. What we need is a hash containing all known descriptor names,
  // but that does not exist.
  for (const auto& range : (*maybe_proto).range()) {
    std::optional<Descriptor*> maybe_descriptor = 
                        GetDescriptor(descriptor, NUMBER_DESCRIPTORS, range.name());
    if (maybe_descriptor) {
      (*maybe_descriptor)->set_range(range.min(), range.max());
    } else if (verbose) {
      cerr << "ReadDescriptorRanges:no match for '" << range.name() << "'\n";
    }
  }

  return 1;
}

static void
DisplayFingerprintOptions(std::ostream& output) {
  output << " -G FILTER         work as a TDT filter\n";
  output << " -G RE=<n>         the number of buckets used in discretising the values\n";
  output << " -G ALL            use all descriptors to generate the fingperint\n";
  output << " -G BEST           from calibration runs, certain descriptors have been designated\n";
  output << "                   as the 'best'. Use these designated to generate the fingperint\n";
  output << " -G d1,d2,d3...    specify individual descriptors to be fingerprinted\n";
  output << " <n>               generate <n> replicates for each bit.\n";
  output << " If any descriptor name is followed by a ':n', that feature will include 'n'\n";
  output << " replicates of that bit in the output. By default, all features get the same number\n";
}

static void
DisplayFilterOptions(std::ostream& output) {
  output <<  R"(
Use the -F option to filter molecules. The syntax uses fortran like comparison operators.
-F w_natoms.lt.50
only writes molecules that have fewer than 50 heavy atoms. Note that when using property filters
output is a smiles file.

Multiple filters are compiled as and conditions,
iwdescr ... -F w_natoms.gt.12 -F w_natoms.lt.50 -F w_nrings.gt.1 -F w_nrings.lt.5 file.smi > filtered.smi
would only write molecules that satisfy all conditions.
)";
}

int
iwdescr(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "g:N:u:A:fq:E:vi:lH:b:T:O:F:a:G:s:zSB:d");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('A'))
    set_global_aromaticity_type(Daylight);
  else if (! process_standard_aromaticity_options(cl, verbose)) {
    usage(7);
  }

  if (verbose)
    set_display_psa_unclassified_atom_mesages(1);
  else
    set_display_psa_unclassified_atom_mesages(0);

  if (verbose == 0) {
    alogp_engine.set_display_error_messages(0);
  }

  if (cl.option_present('S')) {
    set_non_zero_constribution_for_SD2(0);
    if (verbose) {
      cerr << "SD2 atoms given zero weight in PSA calculations\n";
    }
  }

  if (cl.option_present ('O')) {
    if (! descriptors_to_compute.Initialise(cl)) {
      cerr << "Cannot initialise descriptors to compute (-O)\n";
      return 1;
    }
  }

  // We have an initialisation order problem with descriptor range
  // specifications. That cannot be done until the Descriptor array
  // has been allocated, but it is not allocated yet. If there is
  // a range specification given via the -B option, store the file
  // name and use it after descriptors are allocated.
  IWString range_proto_file_name;

  char output_separator = ' ';

  if (cl.option_present('B')) {
    IWString b;
    for (int i = 0; cl.value('B', b, i); ++i) {
      if (b.starts_with("prefix=")) {
        b.remove_leading_chars(7);
        constexpr int kMessage = 0;
        char_name_to_char(b, kMessage);
        descriptor_prefix = b;
        if (verbose)
          cerr << "Descriptors generated with prefix '" << descriptor_prefix << "'\n";
      } else if (b == "quiet") {
        set_display_psa_unclassified_atom_mesages(0);

        if (verbose) {
          cerr << "Will not report unclassified atoms\n";
        }
      } else if (b.starts_with("sep=")) {
        b.remove_leading_chars(4);
        char_name_to_char(b);  // Should check this...
        output_separator = b[0];
        if (verbose) {
          cerr << "Output separator set to '" << output_separator << "'\n";
        }
      } else if (b.starts_with("mxsfsdf=")) {
        b.remove_leading_chars(8);
        if (! b.numeric_value(max_difference_in_ring_size_for_strongly_fused) || max_difference_in_ring_size_for_strongly_fused < 1)
        {
          cerr << "The max difference in ring size in a strongly fused system must be positive\n";
          usage(5);
        }

        if (verbose) {
          cerr << "Strongly fused ring systems will only grow if the ring sizes differ by " << max_difference_in_ring_size_for_strongly_fused << " atoms or less\n";
        }
      } else if (b.starts_with("namexref=")) {
        b.remove_leading_chars(9);
        if (! ReadNameTranslation(b, name_translation)) {
          cerr << "Cannot read name translation '" << b << "'\n";
          return 1;
        }
      } else if (b.starts_with("ranges=")) {
        b.remove_leading_chars(7);
        range_proto_file_name = b;
        if (verbose) {
          cerr << "Descriptor range info read from " << range_proto_file_name << '\n';
        }
      } else if (b == "rpipe") {
        read_descriptor_file_pipeline = 1;
        if (verbose) {
          cerr << "Input assumed to be a descriptor file pipeline\n";
        }
      } else if (b == "wpipe") {
        write_descriptor_file_pipeline = 1;
        if (verbose) {
          cerr << "Will write a descriptor file pipeline\n";
        }
      } else if (b == "flush") {
        flush_after_each_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after each molecule\n";
        }
      } else if (b == "help") {
        DisplayDashBOptions(cerr);
        return 0;
      } else {
        cerr << "Unrecognized -B qualifier '" << b << "'\n";
        DisplayDashBOptions(cerr);
        return 1;
      }
    }
  }

  xlogp::SetIssueUnclassifiedAtomMessages(0);

  FileType input_type = FILE_TYPE_INVALID;
  if (read_descriptor_file_pipeline) {
    // don't care about input type.
  } else if (! cl.option_present('i'))
  {
    if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))
      input_type = FILE_TYPE_SMI;
    else if (! all_files_recognised_by_suffix(cl))
    {
      cerr << "Cannot auto-detect file types\n";
      return 8;
    }
  }
  else if (! process_input_type(cl, input_type))
  {
    cerr << prog_name << ": cannot discern input type\n";
    usage(3);
  }

  process_elements(cl);

  if (cl.option_present('g'))
  {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot initialise chemical standardisation object from command line\n";
      usage(7);
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;
    if (verbose)
      cerr << "Will reduce multi-fragment molecules to largest fragment\n";
  }

  if (cl.option_present('b'))
  {
    if (! cl.value('b', min_hbond_feature_separation) || min_hbond_feature_separation < 1)
    {
      cerr << "The minimum Hydrogen bond feature separation (-b) option must be followed by a whole number\n";
      usage(82);
    }

    if (verbose) {
      cerr << "Will perceive H-Bond features separated by more than " << min_hbond_feature_separation << " bonds or more\n";
    }

    descriptors_to_compute.hbond_descriptors = 1;
  }

  iwdigits.initialise(1024);
  if (cl.option_present('d')) {
    iwdigits.append_to_each_stored_string(".");
    if (verbose) {
      cerr << "Whole numbers written as floats\n";
    }
  }

  set_aromatic_bonds_lose_kekule_identity(0);

  if (cl.option_present('q') && (! process_queries(cl, charge_queries, verbose, 'q')))
  {
    cerr << prog_name << ": cannot process command line queries\n";
    usage(7);
  }

  if (cl.option_present('f'))
  {
    for (int i = 0; i < charge_queries.number_elements(); i++)
    {
      charge_queries[i]->set_max_matches_to_find(1);
    }

    if (verbose)
      cerr << "Queries will match only one embedding\n";
  }

  if (cl.option_present('N'))
  {
    if (! charge_assigner.construct_from_command_line(cl, (verbose > 1), 'N'))
    {
      return 7;
    }
  }

  if (cl.option_present('H'))
  {
    if (! donor_acceptor_assigner.construct_from_command_line(cl, 'H', verbose))
    {
      cerr << "Cannot initialise donor/acceptor assignment object\n";
      usage(4);
    }
  }

  if (! allocate_descriptors())
  {
    cerr << "Yipes, internal error allocating names\n";
    return 15;
  }

  initialise_descriptor_defaults();

  if (range_proto_file_name.size() > 0) {
    if (! ReadDescriptorRanges(range_proto_file_name)) {
      cerr << "Cannot read descriptor range text proto '" << range_proto_file_name << "'\n";
      return 1;
    }
    if (verbose) {
      cerr << "Descriptor ranges read from '" << range_proto_file_name << '\n';
    }
  }

// Aug 99. I built a shell wrapper with '-u 0' as a default. Some want to
// try otherwise. So, scan all instances of the -u option and take the last

  if (cl.option_present('u')) {
    int i = 0;
    while (cl.value('u', undefined_value, i++)) {
      if (verbose)
        cerr << "Undefined values will be written as '" << undefined_value << "'\n";
    }
    if (undefined_value == "none") {
      undefined_value = "";
    }
  }

// Ran into problems trying to do test with a charge assigner present

  if (cl.option_present('T')) {
    if (cl.option_present('N') || charge_assigner.active())
    {
      cerr << "Test mode does not work with a charge assigner\n";
      return 63;
    }

    const_IWSubstring t;
    int i = 0;
    while (cl.value('T', t, i++))
    {
      if ("kg" == t)
      {
        keep_going_after_test_failure = 1;
      }
      else if (! t.numeric_value(ntest) || ntest < 1)
      {
        cerr << "The test option (-T) must be followed by a whole positive number\n";
        usage(4);
      }
    }

    if (0 == ntest)
    {
      cerr << "Must specify the number of test to perform with the -T option\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will perform tests on " << ntest << " random smiles permutations\n";

    saved_result = new Set_or_Unset<float>[NUMBER_DESCRIPTORS];
  }

  number_filters = cl.option_count('F');

  if (number_filters)
  {
    filter = new Descriptor_Filter[number_filters];

    const_IWSubstring f;
    for (int i = 0; i < number_filters; i++)
    {
      cl.value('F', f, i);
      if (f == "help") {
        DisplayFilterOptions(cerr);
        return 1;
      }
      if (! filter[i].build(f))
      {
        cerr << "INvalid filter specification '" << f << "'\n";
        return 3;
      }
    }

    if (verbose)
      cerr << "Defined " << number_filters << " filters\n";
  }

  set_default_iwstring_float_concatenation_precision(output_precision);

  if (cl.option_present('a'))
  {
    if (! cl.value('a', alarm_time) || 0 == alarm_time)
    {
      cerr << "The alarm time (-a) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will time out molecules that take more than " << alarm_time << " seconds\n";
  }

  if (cl.option_present('G'))
  {
    int replicates = 1;
    int resolution = 10;

    int i = 0;
    const_IWSubstring s;

    while (cl.value('G', s, i++))
    {
      if (s == "help") {
        DisplayFingerprintOptions(cerr);
        return 0;
      }

      if ("FILTER" == s) {
        work_as_tdt_filter = 1;
        continue;
      }

      if (s.starts_with("R=")) {
        s.remove_leading_chars(2);
        if (! s.numeric_value(resolution) || resolution < 2) {
          cerr << "The resolution on fingerprints (R=) must be a whole +ve number\n";
          return 1;
        }
        continue;
      }

      int tmp;
      if (s.numeric_value(tmp) && tmp > 0)
      {
        replicates = tmp;
        continue;
      }

      if (s.starts_with("ALL")) {
        int tmpr = replicates;
        if (! parse_replicates_specification(s, tmpr)) {
          cerr << "Invalid ALL,replicates specification '" << s << "'\n";
          return 2;
        }

        for (int i = 0; i < NUMBER_DESCRIPTORS; ++i)
        {
          descriptor[i].set_produce_fingerprint(tmpr);
        }
        continue;
      }

      if (s.starts_with("BEST"))     // leave open possibility of adding a number later
      {
        int tmpr = replicates;
        if (! parse_replicates_specification(s, tmpr))
        {
          cerr << "Invalid BEST,replicates specification '" << s << "'\n";
          return 2;
        }

        for (int i = 0; i < NUMBER_DESCRIPTORS; ++i)
        {
          if (descriptor[i].best_fingerprint())
            descriptor[i].set_produce_fingerprint(tmpr);
        }
        continue;
      }

      const_IWSubstring dname;
      for (int j = 0; s.nextword(dname, j, ',');)
      {
        int tmpr = replicates;
        if (! parse_replicates_specification(dname, tmpr))
        {
          cerr << "Invalid descriptor:replicate specification '" << dname << "'\n";
          return 3;
        }

        int k = descriptor_name_2_index(descriptor, NUMBER_DESCRIPTORS, dname);

        if (k < 0)
        {
          cerr << "Unrecognised descriptor name '" << dname << "' for fingerprint\n";
          return 3;
        }

        cerr << "descriptor " << descriptor[k].descriptor_name() << " getting " << tmpr << " replicates " << replicates << '\n';
        descriptor[k].set_produce_fingerprint(tmpr);
      }
    }

    tag = "NCIWD<";

    fill_descriptor_extremeties(descriptor, resolution);
  }

  if (cl.option_present('s'))
  {
    if (write_descriptor_file_pipeline) {
      cerr << "Writing a descriptor file pipeline is incompatible with the -s option\n";
      return 1;
    }

    include_smiles_as_descriptor = 1;

    if (verbose)
      cerr << "Will include the smiles as an output descriptor\n";
  }

  if (cl.option_present('z'))
  {
    ignore_molecules_with_no_atoms = 0;

    if (verbose)
      cerr << "Will write zero descriptors for empty molecules\n";
  }

  mwc.set_ignore_isotopes(1);

  // Even if alogp is not being computed, set some useful default values.
  alogp_engine.set_use_alcohol_for_acid(1);
  alogp_engine.set_rdkit_phoshoric_acid_hydrogen(1);

  if (cl.empty()) {
    cerr << prog_name << ": no files specified\n";
    usage(6);
  }

  IWString_and_File_Descriptor output(1);

// write the descriptor file header if needed.

  if (ntest)     // no header needed
    ;
  else if (number_filters > 0)   // no header
    ;
  else if (tag.length() > 0)   // fingerprints, no header
    ;
  else if (read_descriptor_file_pipeline)   // handled elsewhere.
    ;
  else if (write_descriptor_file_pipeline) {
    output << "smiles" << output_separator << "Id" << output_separator;
    if (!WriteHeader(descriptor, name_translation, output_separator, output)) {
      cerr << "Cannot write descriptor file header\n";
      return 1;
    }
  }
  else {
    output << "Name" << output_separator;
    if (!WriteHeader(descriptor, name_translation, output_separator, output)) {
      cerr << "Cannot write descriptor file header\n";
      return 1;
    }
  }

//cerr << "HEader contains " << output.size() << " bytes, filters " << number_filters << '\n';

  if (read_descriptor_file_pipeline) {
    for (const char* fname : cl) {
      if (! iwdescriptors_rpipe(fname, output_separator, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  } else {
    for (const char* fname : cl) {
      int rc;
      if (work_as_tdt_filter) {
        rc = iwdescriptors_filter(fname, output_separator, output);
      } else {
        rc = iwdescriptors(fname, input_type, output_separator, output);
      }

      if (! rc) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << molecules_read << " molecules read\n";
    if (molecules_with_no_rings)
      cerr << molecules_with_no_rings << " molecules have no rings\n";

    for (int i = 0; i < NUMBER_DESCRIPTORS; i++)
    {
      if (descriptor[i].active())
        descriptor[i].report_statistics(cerr);
    }

    if (number_filters)
    {
      cerr << rejected_by_filters << " molecules rejected by " << number_filters << " filters\n";
      for (int i = 0; i < number_filters; i++)
      {
        filter[i].report(cerr);
      }
    }
  }

  if (test_failures)
    cerr << test_failures << " test failures\n";
  else if (ntest && verbose)
    cerr << "NO test failures\n";

  if (molecules_skipped_by_timer) {
    cerr << " Skipped " << molecules_skipped_by_timer << " molecules because of timer " << alarm_time << '\n';
  }

  for (int i = 0; i < charge_queries.number_elements(); i++) {
    charge_queries[i]->report(i, cerr);
  }

  delete [] descriptor;

  return 0;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = iwdescr(argc, argv);

  return rc;
}
