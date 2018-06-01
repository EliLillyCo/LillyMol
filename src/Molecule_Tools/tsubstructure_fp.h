#ifndef TSUBSTRUCTURE_FP_H
#define TSUBSTRUCTURE_FP_H

class TSubstructure_FP
{
  private:
    int _work_as_filter;
    int _bit_replicates;
    int _default_fingerprint_nbits;
    IWString _tag;

//  private functions

    int _do_fingerprint_output(int nq, const int * hits, std::ostream & output);
    int _do_sparse_fingerprint_output (int nq, const int * hits, std::ostream & output);

  public:
#ifdef INTERNAL_DISTRIBUTION
    TSubstructure_FP();
#else
    TSubstructure_FP() { return; }
#endif // INTERNAL_DISTRIBUTION

    void set_work_as_filter(int s) { _work_as_filter = s;}
    void set_bit_replicates(int s) { _bit_replicates = s;}
    void set_fingerprint_tag(const const_IWSubstring & s) { _tag = s;}
    void set_default_fingerprint_nbits(int s) { _default_fingerprint_nbits = s;}

#ifdef INTERNAL_DISTRIBUTION
    int do_fingerprint_output(Molecule &, int nq, const int * hits, std::ostream & output);
#else
    int do_fingerprint_output(Molecule &, int nq, const int * hits, std::ostream & output) { return 1; }
#endif // INTERNAL_DISTRIBUTION

    int active() const { return _tag.length();}
};
#ifdef INTERNAL_DISTRIBUTION
extern int write_smiles_and_pcn(Molecule & m, std::ostream & output);
#else
static int write_smiles_and_pcn(Molecule & m, std::ostream & output) { return 1; }
#endif // INTERNAL_DISTRIBUTION


#endif
