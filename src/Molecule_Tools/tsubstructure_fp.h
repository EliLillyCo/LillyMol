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
    TSubstructure_FP();
    void set_work_as_filter(int s) { _work_as_filter = s;}
    void set_bit_replicates(int s) { _bit_replicates = s;}
    void set_fingerprint_tag(const const_IWSubstring & s) { _tag = s;}
    void set_default_fingerprint_nbits(int s) { _default_fingerprint_nbits = s;}
    int do_fingerprint_output(Molecule &, int nq, const int * hits, std::ostream & output);

    int active() const { return _tag.length();}
};
extern int write_smiles_and_pcn(Molecule & m, std::ostream & output);

#endif
