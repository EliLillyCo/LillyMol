#include "sparse_fp_creator.h"
#include "iwbits.h"

#include "molecule.h"
#include "tsubstructure_fp.h"

TSubstructure_FP::TSubstructure_FP()
{
  _work_as_filter = 0;
  _bit_replicates = 1;
  _default_fingerprint_nbits = 0;

  return;
}

/*
  We have two scenarios for fingerprint output. 
    We are reading a file of molecules
    We are reading a TDT file as a filter
*/

int
TSubstructure_FP::_do_fingerprint_output(int nq,
                      const int * hits,
                      std::ostream & output)
{
  IW_Bits_Base fp;
  if (_default_fingerprint_nbits)
    fp.allocate_space_for_bits(_default_fingerprint_nbits * _bit_replicates);

  if (1 == _bit_replicates)
    fp.construct_from_array_of_ints(hits, _default_fingerprint_nbits);   // 2nd arg default_fingerprint_nbits rather than nq to make sure the bit vector is filled all the way
  else
  {
    for (int i = 0; i < nq; i++)
    {
      if (0 == hits[i])
        continue;

      for (int j = 0; j < _bit_replicates; j++)
      {
        if (hits[j])
          fp.set(j * _default_fingerprint_nbits + i);
      }
    }
  }

  fp.write_daylight_ascii_representation(output, _tag);

  return output.good();
}

int
TSubstructure_FP::_do_sparse_fingerprint_output (int nq,
                              const int * hits,
                              std::ostream & output)
{
  Sparse_Fingerprint_Creator sfp;

  if (1 == _bit_replicates)
    sfp.create_from_array_of_ints(hits, nq);
  else
  {
    for (int i = 0; i < nq; i++)
    {
      if (0 == hits[i])
        continue;

      for (int j = 0; j < _bit_replicates; j++)
      {
        sfp.hit_bit(j * _bit_replicates + i);
      }
    }
  }

  IWString ascii_rep;
  sfp.daylight_ascii_form_with_counts_encoded(ascii_rep);

  output << _tag << ascii_rep << ">\n";

  return output.good();
}

int
write_smiles_and_pcn(Molecule & m,
                     std::ostream & output)
{
  output << "$SMI<" << m.smiles() << ">\n";
  output << "PCN<" << m.molecule_name() << ">\n";

  return output.good();
}

int
TSubstructure_FP::do_fingerprint_output(Molecule & m,
                      int nq,
                      const int * hits,
                      std::ostream & output)
{
  if (! _work_as_filter)
   write_smiles_and_pcn(m, output);

  if (_tag.starts_with("NC"))
    (void) _do_sparse_fingerprint_output(nq, hits, output);
  else
    (void) _do_fingerprint_output(nq, hits, output);

  if (! _work_as_filter)
    output << "|\n";

  return output.good();
}
