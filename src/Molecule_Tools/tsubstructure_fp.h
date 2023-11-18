#ifndef TSUBSTRUCTURE_FP_H
#define TSUBSTRUCTURE_FP_H

#include <ostream>

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/molecule.h"

class TSubstructure_FP
{
  private:
    int _work_as_filter;
    int _bit_replicates;
    int _default_fingerprint_nbits;
    IWString _tag;

    // We can optionally set a bit for the total number of hits
    // across all the queries. Note that the value of the variable controls
    // the number of bit replicates associated with the total hits bit.
    int _extra_bit_total_hits;

//  private functions

    template <typename OUTPUT> int _do_fingerprint_output(int nq, const int * hits, OUTPUT& output) const ;
    template <typename OUTPUT> int _do_sparse_fingerprint_output (int nq, const int * hits, OUTPUT& output) const ;

  public:
    TSubstructure_FP();

    void set_work_as_filter(int s) { _work_as_filter = s;}
    void set_bit_replicates(int s) { _bit_replicates = s;}
    void set_fingerprint_tag(const const_IWSubstring & s) { _tag = s;}
    void set_default_fingerprint_nbits(int s) { _default_fingerprint_nbits = s;}
    void set_extra_bit_total_hits(int s) {
      _extra_bit_total_hits = s;
    }

    template <typename OUTPUT>
      int do_fingerprint_output(Molecule &, int nq, const int * hits, OUTPUT& output) const;

    int active() const { return _tag.length();}
};

template <typename OUTPUT> int write_smiles_and_pcn(Molecule& m, OUTPUT& output);

#endif
