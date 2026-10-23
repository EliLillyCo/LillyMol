#ifndef MOLECULE_TOOLS_IWDESCR_INTERNAL_H_
#define MOLECULE_TOOLS_IWDESCR_INTERNAL_H_

#include <iosfwd>
#include <iostream>

#include "Foundational/iwstring/iwstring.h"
class iwstring_data_source;

class const_IWSubstring;

// Private implementation details for iwdescr.
//
// This header is intentionally not part of the stable public API. It exists so
// iwdescr_lib.cc can make the IWDescr-owned state concrete without exposing the
// Descriptor implementation through iwdescr_lib.h.
//
// Descriptor is entirely internal. DescriptorsToCompute is currently internal to
// the C++ implementation, but IWDescr exposes const/mutable accessors so it can
// later be bound to Python or moved to a narrower public configuration API.

enum IWDescr_Enum {
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
  iwdescr_d2sp3,
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
  iwdescr_scafrotb,
  iwdescr_spchrotb,
  iwdescr_maxrotbgrp,
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
  iwdescr_biphenyl,
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
  // iwdescr_scra,
  // iwdescr_scrha,
  // iwdescr_scrd,
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
  iwdescr_internalhbd,
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
  iwdescr_nrings8,  // must keep this in sync with MAX_RING_SIZE
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
  iwdescr_mxlencchain2,
  iwdescr_mxlencchain3,
  // Saturated chain features.
  iwdescr_nsatchain,
  iwdescr_mxsatchain,
  iwdescr_fsatchain,
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
  iwdescr_ro5_on,  // the last one which will always be computed

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
  iwdescr_rmy_heavy_halogen  // remember to update NUMBER_DESCRIPTORS below.
};

constexpr int kNumberDescriptors = iwdescr_rmy_heavy_halogen + 1;


// Various features can be computed optionally. Put them all into a struct that
// can be initialised with the -O option. Some of these have default values of 1,
// some 0.
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
  int long_carbon_chains = 1;
  int saturated_chains = 1;

  void SetAll(int s);
  int DisplayUsage(std::ostream& output) const;
  int ReadDescriptorsToCompute(const const_IWSubstring& fname);
  int ReadDescriptorsToCompute(iwstring_data_source& input);
  int Initialise(Command_Line& cl);
};

#endif  // MOLECULE_TOOLS_IWDESCR_INTERNAL_H_
