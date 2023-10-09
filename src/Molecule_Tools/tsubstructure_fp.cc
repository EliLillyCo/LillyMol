#include <iostream>

#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "tsubstructure_fp.h"

using std::cerr;
using std::endl;

static IWString smiles_tag("$SMI<");
static IWString pcn_tag("PCN<");


TSubstructure_FP::TSubstructure_FP() {
   _work_as_filter = 0;
   _bit_replicates = 1;
   _default_fingerprint_nbits = 1;
}

template <typename OUTPUT> int
TSubstructure_FP::do_fingerprint_output(Molecule & m, const int nq, const int * hits, OUTPUT & output) const {
  if (0 == _tag.length()) {
    cerr << "TSubstructure_FP::do_fingerprint_output:no tag\n";
    return 0;
  }

  if (! _work_as_filter) {
    write_smiles_and_pcn(m, output);
  }
  
  if (_tag.starts_with("FP"))
    return _do_fingerprint_output(nq, hits, output);
  if (_tag.starts_with("NC"))
    return _do_sparse_fingerprint_output(nq, hits, output);

  cerr << "TSubstructure_FP::do_fingerprint_output:unknown tag type " << _tag << endl;
  return 0;
}

template int TSubstructure_FP::do_fingerprint_output(Molecule& m, int nq, const int * hits, std::ostream& output) const;

template <typename OUTPUT> int
TSubstructure_FP::_do_sparse_fingerprint_output(const int nq, const int * hits, OUTPUT& output) const
{
  Sparse_Fingerprint_Creator sfc;

  for (int i = 0; i < nq; ++i) {
    if (0 == hits[i])
      continue;

    sfc.hit_bit(i * _bit_replicates, hits[i]);  // always.
    for (int j = 1; j < _bit_replicates; ++j) {
      sfc.hit_bit(i * _bit_replicates + j, hits[i]);
    }
  }

  sfc.write_fingerprint(_tag, output);
  if (! _work_as_filter) {
    output << "|\n";
  }

  return 1;
}

template <typename OUTPUT> int
TSubstructure_FP::_do_fingerprint_output(const int nq, const int * hits, OUTPUT& output) const
{
  int bits_needed = 0;
  if (_default_fingerprint_nbits) {
    bits_needed = _default_fingerprint_nbits * _bit_replicates;
  } else if (_bit_replicates > 0) {
    bits_needed = nq * _bit_replicates;
  } else {
    bits_needed = nq;
  }

  if (0 != bits_needed % 8) {
    bits_needed = (bits_needed / 8 + 1) * 8;
  }

  IW_Bits_Base fp;
  fp.allocate_space_for_bits(bits_needed);

  for (int i = 0; i < nq; ++i) {
    if (0 == hits[i])
      continue;

    fp.set(i * _bit_replicates);
    for (int j = 1; j < _bit_replicates; ++j) {
      fp.set(i * _bit_replicates + j);
    }
  }

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);
  output << _tag << tmp << ">\n";
  if (! _work_as_filter) {
    output << "|\n";
  }

  return 1;
}

template <typename OUTPUT> int
write_smiles_and_pcn(Molecule& m, OUTPUT& output) {
  output << smiles_tag << m.smiles() << ">\n";
  output << pcn_tag << m.name() << ">\n";

  return output.good();
}

template int write_smiles_and_pcn(Molecule&, std::ostream&);
