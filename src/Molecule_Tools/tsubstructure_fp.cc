#include <iostream>
#include <limits>

#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "tsubstructure_fp.h"

using std::cerr;

static IWString smiles_tag("$SMI<");
static IWString pcn_tag("PCN<");

TSubstructure_FP::TSubstructure_FP() {
   _work_as_filter = 0;
   _bit_replicates = 1;
   _default_fingerprint_nbits = 1;
   _extra_bit_total_hits = 0;
}

template <typename OUTPUT> int
TSubstructure_FP::do_fingerprint_output(Molecule & m, const int nq, const int * hits, OUTPUT & output) const {
  if (_tag.empty()) {
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

  cerr << "TSubstructure_FP::do_fingerprint_output:unknown tag type " << _tag << '\n';
  return 0;
}

template int TSubstructure_FP::do_fingerprint_output(Molecule& m, int nq, const int * hits, std::ostream& output) const;

template <typename OUTPUT> int
TSubstructure_FP::_do_sparse_fingerprint_output(const int nq, const int * hits, OUTPUT& output) const
{
  Sparse_Fingerprint_Creator sfc;

  int total_hits = 0;
  for (int i = 0; i < nq; ++i) {
    if (0 == hits[i]) {
      continue;
    }

    total_hits += hits[i];

    sfc.hit_bit(i * _bit_replicates, hits[i]);  // always.
    for (int j = 1; j < _bit_replicates; ++j) {
      sfc.hit_bit(i * _bit_replicates + j, hits[i]);
    }
  }

  if (_extra_bit_total_hits && total_hits > 0) {
    if (total_hits > std::numeric_limits<uint8_t>::max()) {
      total_hits = std::numeric_limits<uint8_t>::max();
    }

    for (int i = 0; i < _extra_bit_total_hits; ++i) {
      sfc.hit_bit(_bit_replicates * nq, total_hits);
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

  const int extra_bit_bstart = bits_needed;

  // In a binary fingerprint this just indicates whether any of the queries matched.
  if (_extra_bit_total_hits) {
    bits_needed += _extra_bit_total_hits;
  }

  if (0 != bits_needed % 8) {
    bits_needed = (bits_needed / 8 + 1) * 8;
  }

  IW_Bits_Base fp;
  fp.allocate_space_for_bits(bits_needed);

  // If we are writing an extra bit for the total number of hits.
  int total_hits = 0;
  for (int i = 0; i < nq; ++i) {
    if (0 == hits[i]) {
      continue;
    }

    ++total_hits;

    fp.set(i * _bit_replicates);
    for (int j = 1; j < _bit_replicates; ++j) {
      fp.set(i * _bit_replicates + j);
    }
  }

  if (total_hits > 0 && _extra_bit_total_hits) {
    //int bstart = nq * _bit_replicates;
    for (int i = 0; i < _bit_replicates; ++i) {
      fp.set(extra_bit_bstart + i);
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
