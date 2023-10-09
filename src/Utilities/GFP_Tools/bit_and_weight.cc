#include "bit_and_weight.h"

void
Bit_and_Weight::increment_weight(unsigned int b,
                                 double & rc) const
{
//IW_Hash_Map<unsigned int, double>::const_iterator f = _weight.find(b);
  auto f = _weight.find(b);

  if (f == _weight.end())
    return;

  rc += (*f).second;

  return;
}

Fixed_and_Sparse_Bit_and_Weight::Fixed_and_Sparse_Bit_and_Weight()
{
  _nfp = NBC_FIXED_ARRAY_SIZE;
  _fp = new Bit_and_Weight[_nfp];

  _nsfp = NBC_FIXED_ARRAY_SIZE;
  _sparsefp = new Bit_and_Weight[_nsfp];

  return;
}

Fixed_and_Sparse_Bit_and_Weight::~Fixed_and_Sparse_Bit_and_Weight()
{
  if (nullptr != _fp)
    delete [] _fp;

  if (nullptr != _sparsefp)
    delete [] _sparsefp;

  return;
}

int
Fixed_and_Sparse_Bit_and_Weight::status(std::ostream & os) const
{
  os << "Fixed_and_Sparse_Bit_and_Weight::status: ";

  if (0 == _nfp && 0 == _nsfp)
  {
    os << "inactive\n";
    return 1;
  }

  os << _nfp << " fixed, " << _nsfp << " sparse fingerprints\n";

  for (int i = 0; i < _nfp; i++)
  {
    if (_fp[i].size() > 0)
      os << "  fp " << i << " has " << _fp[i].size() << " bits\n";
  }

  for (int i = 0; i < _nsfp; i++)
  {
    if (_sparsefp[i].size() > 0)
      os << " sfp " << i << " has " << _sparsefp[i].size() << " bits\n";
  }

  return 1;
}

int
Fixed_and_Sparse_Bit_and_Weight::initialise(const IW_General_Fingerprint & fp)
{
  _nfp = fp.nfingerprints();
  if (_nfp)
    _fp = new Bit_and_Weight[_nfp];

  _nsfp = fp.number_sparse_fingerprints();

  if (_nsfp)
    _sparsefp = new Bit_and_Weight[_nsfp];

  return 1;
}

/*
  In THIS we have bit numbers and probability values in some large collection,
  likely the Lilly collection. RHS contains the distributions for what is
  presumed to be a set of active molecules. Change that into something that
  can be used for Naive Bayesian scoring
*/

int
Fixed_and_Sparse_Bit_and_Weight::build_scoring_model (Fixed_and_Sparse_Bit_and_Weight & rhs,
                                        int population) const
{
  for (int i = 0; i < _nfp; i++)
  {
    if (_fp[i].size() > 0)
      _fp[i].build_scoring_model(rhs._fp[i], population);
  }

  for (int i = 0; i < _nsfp; i++)
  {
    if (_sparsefp[i].size() > 0)
      _sparsefp[i].build_scoring_model(rhs._sparsefp[i], population);
  }

  return 1;
}

/*
  On entry, RHS contains bit numbers and counts. POPULATION is
  the number of molecules used to derive the profile held in RHS.

  We compare the probability implied by RHS with the probabilities
  in THIS, and convert RHS to a scoring model
*/

int
Bit_and_Weight::build_scoring_model(Bit_and_Weight & rhs, int population) const
{
  double reciprocal_p = 1.0 / static_cast<double>(population);

//for (BWHT::iterator i = rhs._weight.begin(); i != rhs._weight.end(); ++i)
  for (auto i = rhs._weight.begin(); i != rhs._weight.end(); ++i)
  {
    unsigned int b = (*i).first;

    BWHT::const_iterator f = _weight.find(b);

    if (f == _weight.end())   // bit found in RHS missing from THIS
    {
      (*i).second = 1.0;   // zero contribution
      continue;
    }

    double baseline = (*f).second;    // fraction hit in the large population

    double rhs_pop = (*i).second * reciprocal_p;   // population in RHS

    (*i).second = rhs_pop / baseline;
  }

  return 1;
}

int
Fixed_and_Sparse_Bit_and_Weight::write_scoring_model (IWString_and_File_Descriptor & output) const
{
  for (int i = 0; i < _nfp; i++)
  {
    if (0 == _fp[i].size())
      continue;

    IWString prefix;

    prefix << "FP " << i;

    _fp[i].write_scoring_model(prefix, output);

    if (output.length() > 32768)
      output.write_whole_blocks_shift_unwritten();
  }

  for (int i = 0; i < _nsfp; i++)
  {
    if (0 == _sparsefp[i].size())
      continue;

    IWString prefix;

    prefix << "SFP " << i;

    _sparsefp[i].write_scoring_model(prefix, output);

    if (output.length() > 32768)
      output.write_whole_blocks_shift_unwritten();
  }

  return 1;
}

int
Bit_and_Weight::write_scoring_model(const IWString & prefix,
                                    IWString_and_File_Descriptor & output) const
{
//for (BWHT::const_iterator i = _weight.begin(); i != _weight.end(); ++i)
  for (auto i = _weight.cbegin(); i != _weight.cend(); ++i)
  {
    unsigned int b = (*i).first;
    double w = (*i).second;

    output << prefix << " bit " << b << " NB " << w << '\n';
  }

  return 1;
}
