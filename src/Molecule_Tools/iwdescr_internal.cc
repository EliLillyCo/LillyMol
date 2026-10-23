
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Tools/iwdescr_internal.h"
#include "Molecule_Tools/iwdescr_lib.h"

using std::cerr;

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
  long_carbon_chains = s;
  saturated_chains = s;
  ring_fusion_descriptors = s;
  ramey_descriptors = s;
}

int
DescriptorsToCompute::DisplayUsage(std::ostream& output) const {
  // clang-format off
  output << R"(Specify descriptor(s) to be disabled or enabled.
If the feature name starts with '-', the feature is disabled.
 -O none          turn off all optional descriptors
 -O all           turn on  all optional descriptors
 -O adjring       enable adjacent to ring fusion descriptors
 -O bbr           enable bonds between rings descriptors
 -O charge        enable all formal charge descriptors
 -O chiral        enable all expensive chirality perception descriptors
 -O complex       enable planar fused ring descriptors
 -O crowd         enable atomic crowding descriptors
 -O lcc           enable atomic crowding descriptors
 -O dm            enable descriptors based on the distance matrix
 -O donacc        enable donor acceptor derived descriptors
 -O hbond         enable donor/acceptor derived descriptors
 -O shbond        enable simplistic donor/acceptor derived descriptors
 -O ncon          enable ncon and fncon descriptors
 -O pbond         enable polar bond derived descriptors
 -O psa           enable Novartis polar surface area descriptors
 -O psymm         enable partial symmetry derived descriptors
 -O ramey         enable Ramey (element count) descriptors
 -O rcj           enable ring chain join descriptors
 -O rfuse         enable ring fusion descriptors
 -O rss           enable ring substitution descriptors
 -O rssr          enable ring substitution ratio descriptors
 -O spch          enable spinach related descriptors
 -O spcgrp        enable specific group descriptors
 -O satchain      enable computation of saturated chain descriptors
 -O symm          enable symmetry related descriptors
 -O xlogp         enable xlogp
)";
  // clang-format on

  return 0;
}

int
DescriptorsToCompute::Initialise(Command_Line& cl) {
  static constexpr char kFlag = 'O';
  const int verbose = cl.option_present('v');
  const_IWSubstring o;
  for (int i = 0; cl.value(kFlag, o, i); ++i) {
    int s = 1;
    if (o[0] == '-') {
      s = 0;
      o++;
    }
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
      adjacent_ring_fusion_descriptors = s;
      if (verbose) {
        cerr << "bonding adjacent to ring fusion descriptors enabled\n";
      }
    } else if (o == "bbr") {
      bonds_between_rings = s;
      if (verbose) {
        cerr << "bonds between rings descriptors enabled\n";
      }
    } else if (o == "charge") {
      charge_descriptors = s;
      if (verbose) {
        cerr << "charge derived descriptors enabled\n";
      }
    } else if (o == "complex") {
      complexity_descriptors = s;
      if (verbose) {
        cerr << "complexity derived descriptors enabled\n";
      }
    } else if (o == "crowd") {
      crowding_descriptors = s;
      if (verbose) {
        cerr << "atomic crowding descriptors enabled\n";
      }
    } else if (o == "lcc") {
      long_carbon_chains = s;
      if (verbose) {
        cerr << "long carbon chain enabled\n";
      }
    } else if (o == "dm") {
      distance_matrix_descriptors = s;
      if (verbose) {
        cerr << "distance matrix derived descriptors enabled\n";
      }
    } else if (o == "donacc") {
      donor_acceptor = s;
      if (verbose) {
        cerr << "donor acceptor derived descriptors enabled\n";
      }
    } else if (o == "hbond") {
      hbond_descriptors = s;
      if (verbose) {
        cerr << "hydrogen bonding derived descriptors enabled\n";
      }
    } else if (o == "shbond") {
      simple_hbond_descriptors = s;
      if (verbose) {
        cerr << "simplistic hydrogen bonding derived descriptors enabled\n";
      }
    } else if (o == "chiral" || o == "expchiral") {
      perform_expensive_chirality_perception = s;
      if (verbose) {
        cerr << "Will check all unlabelled atoms for possible chirality\n";
      }
    } else if (o == "ncon") {
      ncon_descriptors = s;
      if (verbose) {
        cerr << "ncon and fncon descriptors enabled\n";
      }
    } else if (o == "pbond") {
      polar_bond_descriptors = s;
      if (verbose) {
        cerr << "polar bond derived descriptors enabled\n";
      }
    } else if (o == "psa") {
      psa = s;
      if (verbose) {
        cerr << "Novartis polar surface area descriptor enabled\n";
      }
    } else if (o == "psymm") {
      partial_symmetry_descriptors = s;
      if (verbose) {
        cerr << "partial symmetry derived descriptors enabled\n";
      }
    } else if (o == "symm") {
      symmetry_descriptors = s;
      if (verbose) {
        cerr << "symmetry derived descriptors enabled\n";
      }
    } else if (o == "ramey") {
      ramey_descriptors = s;
      if (verbose) {
        cerr << "element count descriptors enabled\n";
      }
    } else if (o == "rcj") {
      ring_chain_descriptors = s;
      if (verbose) {
        cerr << "ring chain join descriptors enabled\n";
      }
    } else if (o == "rfuse") {
      ring_fusion_descriptors = s;
      if (verbose) {
        cerr << "ring fusion descriptors enabled\n";
      }
    } else if (o == "rss") {
      ring_substitution_descriptors = s;
      if (verbose) {
        cerr << "ring substitution descriptors enabled\n";
      }
    } else if (o == "rssr") {
      ring_substitution_ratio_descriptors = s;
      if (verbose) {
        cerr << "ring substitution ratio descriptors enabled\n";
      }
    } else if (o == "spch") {
      spinach_descriptors = s;
      if (verbose) {
        cerr << "spinach derived descriptors enabled\n";
      }
    } else if (o == "spcgrp") {
      specific_groups = s;
      if (verbose) {
        cerr << "specific substructure descriptors enabled\n";
      }
    } else if (o == "alogp") {
      compute_alogp = s;
      if (verbose) {
        cerr << "alogp will be computed\n";
      }
    } else if (o == "xlogp") {
      compute_xlogp = s;
      if (verbose) {
        cerr << "xlogp will be computed\n";
      }
    } else if (o == "satchain") {
      saturated_chains = s;
      if (verbose) {
        cerr << "Saturated chain descriptors computed\n";
      }
    } else if (o.starts_with("F:")) {
      o.remove_leading_chars(2);
      if (!ReadDescriptorsToCompute(o)) {
        cerr << "DescriptorsToCompute::Initialise:cannot process '" << o << "'\n";
        return 0;
      }
    } else if (o == "help") {
      return DisplayUsage(cerr);
    } else {
      cerr << "DescriptorsToCompute::Initialise:unrecognised '" << o << "'\n";
      return DisplayUsage(cerr);
    }
  }

  return 1;
}

int
DescriptorsToCompute::ReadDescriptorsToCompute(const const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "DescriptorsToCompute::ReadDescriptorsToCompute:cannot open '" << fname
         << "'\n";
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
      perform_expensive_chirality_perception = 1;
    } else if (o == "ncon") {
      ncon_descriptors = 1;
    } else if (o == "pbond") {
      polar_bond_descriptors = 1;
    } else if (o == "psa") {
      psa = 1;
    } else if (o == "psymm") {
      partial_symmetry_descriptors = 1;
    } else if (o == "symm") {
      symmetry_descriptors = 1;
    } else if (o == "ramey") {
      ramey_descriptors = 1;
    } else if (o == "rcj") {
      ring_chain_descriptors = 1;
    } else if (o == "rfuse") {
      ring_fusion_descriptors = 1;
    } else if (o == "rss") {
      ring_substitution_descriptors = 1;
    } else if (o == "rssr") {
      ring_substitution_ratio_descriptors = 1;
    } else if (o == "spch") {
      spinach_descriptors = 1;
    } else if (o == "spcgrp") {
      specific_groups = 1;
    } else if (o == "alogp") {
      compute_alogp = 1;
    } else if (o == "xlogp") {
      compute_xlogp = 1;
    } else {
      cerr << "DescriptorsToCompute::ReadDescriptorsToCompute:unrecognised '" << o
           << "'\n";
      return 0;
    }
  }

  return 1;
}

using iwmisc::Fraction;
using std::cerr;

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
look_for_dot_operators(const const_IWSubstring& s, IWString& dname, int& op,
                       int& begin_numeric) {
  int n = s.length() - 4;  // min would be dname.op.1

  for (int i = 1; i < n; i++) {
    if ('.' != s[i]) {
      continue;
    }

    if ('.' != s[i + 3]) {
      return 0;
    }

    if ('e' == s[i + 1] && 'q' == s[i + 2]) {
      s.from_to(0, i - 1, dname);
      op = DFILTER_OP_EQUALS;
      begin_numeric = i + 4;
      return 1;
    }
    if ('l' == s[i + 1] && 't' == s[i + 2]) {
      s.from_to(0, i - 1, dname);
      op = DFILTER_OP_LT;
      begin_numeric = i + 4;
      return 1;
    }
    if ('g' == s[i + 1] && 't' == s[i + 2]) {
      s.from_to(0, i - 1, dname);
      op = DFILTER_OP_GT;
      begin_numeric = i + 4;
      return 1;
    }

    return 0;
  }

  return 0;
}

Descriptor_Filter::Descriptor_Filter() {
  return;
}

// Find the descriptor in `descriptors` that has name `name`.
int
DescriptorNumber(const IWString& prefix, const IWString& dname,
        const std::span<const Descriptor>& descriptors) {
  IWString mydname(dname);
  if (! prefix.empty() && dname.starts_with(prefix)) {
    mydname.remove_leading_chars(prefix.length());
  }

#ifdef LILLYMOL_HAS_CPP23
  for (auto [i, d] : std::views::enumerate(descriptors)) {
#else
  for (uint32_t i = 0; i < descriptors.size(); ++i) {
    const Descriptor& d = descriptors[i];
#endif
    if (mydname == d.descriptor_name()) {
      return i;
    }
  }

  return -1;
}

int
Descriptor_Filter::build(const IWString& prefix, 
                         const const_IWSubstring& s,
                         const std::span<const Descriptor>& descriptors) {
  int op = 0;
  int begin_numeric = -1;

  int n = s.length() - 1;  // min form is dname=3

  for (int i = 1; i < n; i++) {
    char c = s[i];

    if ('=' == c) {
      op = DFILTER_OP_EQUALS;
    } else if ('>' == c) {
      op = DFILTER_OP_GT;
    } else if ('<' == c) {
      op = DFILTER_OP_LT;
    } else if (':' == c) {
      op = DFILTER_OP_LIST;
    } else {
      continue;
    }

    s.from_to(0, i - 1, _descriptor_name);
    begin_numeric = i + 1;
    break;
  }

  if (_descriptor_name.empty()) {
    look_for_dot_operators(s, _descriptor_name, op, begin_numeric);
  }

  if (_descriptor_name.empty()) {
    cerr << "Descriptor_Filter::build:no operator found '" << s << "'\n";
    return 0;
  }

  const_IWSubstring numeric_stuff;

  s.from_to(begin_numeric, s.length() - 1, numeric_stuff);

  _descriptor_number = DescriptorNumber(prefix, _descriptor_name, descriptors);

  if (_descriptor_number < 0) {
    cerr << "Descriptor_Filter::build:cannot find '" << _descriptor_name << "'\n";
    return 0;
  }

  if (DFILTER_OP_LIST != op) {
    float v;
    if (!numeric_stuff.numeric_value(v)) {
      cerr << "Descriptor_Filter::build:invalid numeric '" << numeric_stuff << "'\n";
      return 0;
    }

    if (DFILTER_OP_EQUALS == op) {
      _cond.add(v);
    } else if (DFILTER_OP_LT == op) {
      if (v == 0.0) {
        _cond.set_max(-std::numeric_limits<float>::min());
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
  while (numeric_stuff.nextword(token, i, ',')) {
    float v;
    if (!token.numeric_value(v)) {
      cerr << "Descriptor_Filter::build:non numeric in list specification '" << s
           << "'\n";
      return 0;
    }

    _cond.add(v);
  }

  return 1;
}

int
Descriptor_Filter::satisfied(const Descriptor* descriptors) {
  const Descriptor& myd = descriptors[_descriptor_number];

  float v;

  if (!myd.value(v)) {  // no value available
    return 0;
  }

  if (_cond.matches(v)) {
    return 1;
  }

  _items_rejected++;

  return 0;
}

int
Descriptor_Filter::Report(std::ostream& os) const {
  os << "Filter on '" << _descriptor_name << "', rejected " << _items_rejected
     << " items\n";

  return 1;
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
Descriptor::report_statistics(std::ostream& os) const {
  const uint32_t nsamples = _stats.n();
  os << "Descriptor '" << _name << "' " << nsamples << " values sampled\n";
  if (nsamples == 0) {
    return 1;
  }

  os << " between " << _stats.minval() << " and " << _stats.maxval();
  if (nsamples > 1) {
    os << " ave " << static_cast<float>(_stats.average()) << " sd " <<
          static_cast<float>(_stats.variance());
  }

  os << ' ' << _zero_value_count << " instances of zero\n";

  return os.good();
}

void
Descriptor::produce_fingerprint(int bitnum, Sparse_Fingerprint_Creator& sfc) const {
  if (0 == _fingerprint_replicates) {
    cerr << "Descriptor::produce_fingerprint:not initialised '" << _name << "'\n";
    return;
  }

  float v;

  if (!Set_or_Unset<float>::value(v)) {
    return;
  }

  int c;

  if (v <= _min) {
    c = 1;
  } else {
    if (v >= _max) {
      v = _max;
    }

    c = static_cast<int>((v - _min) / _dy + 0.5f);

    c++;  // ensure non zero
  }


  // cerr << "Descriptor::produce_fingerprint: descriptor " << _name << " value " << v <<
  // " bit " << bitnum << " count " << c << " rep " << _fingerprint_replicates << '\n';

  if (1 == _fingerprint_replicates) {  // presumably a common case
    sfc.hit_bit(bitnum, c);
    return;
  }

  for (int i = 0; i < _fingerprint_replicates; ++i) {
    sfc.hit_bit(bitnum + i, c);
  }

  return;
}

void
Descriptor::reset() {
  if (_default_value) {
    Set_or_Unset<float>::set(_default_value.value());
  } else {
    Set_or_Unset<float>::unset();
  }

  return;
}

void
Descriptor::set_min_max_resolution(float v1, float v2, int r) {
  assert(v1 < v2);

  _min = v1;
  _max = v2;

  _dy = (v2 - v1) / static_cast<float>(r);

  return;
}

void
Descriptor::set(int s) {
  Set_or_Unset<float>::set(static_cast<float>(s));

  return;
}

void
FillDescriptorExtremeties(std::span<Descriptor>&  d, const int resolution) {
  d[iwdescr_natoms].set_min_max_resolution(3.0f, 50.0f, resolution);
  d[iwdescr_nrings].set_min_max_resolution(0.0f, 9.35196f, resolution);
  d[iwdescr_nrings3].set_min_max_resolution(0.0f, 6.0f, resolution);
  d[iwdescr_nrings4].set_min_max_resolution(0.0f, 8.0f, resolution);
  d[iwdescr_nrings5].set_min_max_resolution(0.0f, 7.0f, resolution);
  d[iwdescr_nrings6].set_min_max_resolution(0.0f, 7.76418f, resolution);
  d[iwdescr_nrings7].set_min_max_resolution(0.0f, 3.0f, resolution);
  d[iwdescr_nrings8].set_min_max_resolution(0.0f, 2.0f, resolution);
  d[iwdescr_nelem].set_min_max_resolution(2.0f, 8.0f, resolution);
  d[iwdescr_amw].set_min_max_resolution(100.0f, 600.0f, resolution); 
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
  d[iwdescr_d2sp3].set_min_max_resolution(0.0f, 19.30857f, resolution);
  d[iwdescr_chmltbd].set_min_max_resolution(0.0f, 7.63907f, resolution);
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
  d[iwdescr_srsz].set_min_max_resolution(3.0f, 8.0f, resolution);        // fix some time
  d[iwdescr_lrsz].set_min_max_resolution(4.0f, 8.0f, resolution);        // fix some time
  d[iwdescr_rng7atoms].set_min_max_resolution(4.0f, 16.0f, resolution);  // fix some time
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
  // d[iwdescr_scra].set_min_max_resolution(0.0f, 20.0f, resolution);
  // d[iwdescr_scrha].set_min_max_resolution(0.0f, 20.0f, resolution);
  // d[iwdescr_scrd].set_min_max_resolution(0.0f, 10.0f, resolution);
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
  d[iwdescr_internalhbd].set_min_max_resolution(0.0f, 10.0f, resolution);
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
  d[iwdescr_mxlencchain2].set_min_max_resolution(0.0f, 20.0f, resolution);
  d[iwdescr_mxlencchain3].set_min_max_resolution(0.0f, 20.0f, resolution);
  d[iwdescr_nsatchain].set_min_max_resolution(0.0f, 20.0f, resolution);
  d[iwdescr_mxsatchain].set_min_max_resolution(0.0f, 20.0f, resolution);
  d[iwdescr_fsatchain].set_min_max_resolution(0.0f, 1.0f, resolution);
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
  d[iwdescr_hcount].set_min_max_resolution(0.0f, 20.0f, resolution);   // just a guess
  d[iwdescr_hperatom].set_min_max_resolution(0.0f, 1.0f, resolution);  // just a guess
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

void
Descriptor::set_name(const char* newname) {
  _active = 1;

  _name = newname;

  return;
}

Descriptor::Descriptor() {
  _active = 0;

  _zero_value_count = 0;

  _fingerprint_replicates = 0;

  _best_fingerprint = 0;

  return;
}

void
Descriptor::update_statistics(int verbose) {
  float v;
  if (!Set_or_Unset<float>::value(v)) {
    if (verbose > 2) {
      cerr << "Descriptor::update_statistics: '" << _name << "' not set\n";
    }
    return;
  }

  _stats.extra(v);

  if (v == 0.0f) {
    _zero_value_count++;
  }

  return;
}
