#ifndef NAIVE_BAYESIAN_CLASSIFIER_H
#define NAIVE_BAYESIAN_CLASSIFIER_H

#include <math.h>

#include "iw_stl_hash_map.h"

#include "gfp.h"
#include "bit_and_weight.h"

class Naive_Bayesian_Classifier
{
  private:
    Fixed_and_Sparse_Bit_and_Weight _fsbw;

    IWString _zclass;

//  private functions

    int _build(const const_IWSubstring & buffer);

    int _classify(int, const Sparse_Fingerprint &, double &) const;
    int _classify(int, const IWDYFP &, double &) const;

  public:
    Naive_Bayesian_Classifier();

    int build(const const_IWSubstring & fname);
    int build(iwstring_data_source &);

    int set_sfp_weight(int f, unsigned int b, double);
    int set_fp_weight (int f, unsigned int b, double);

    double classify(const IW_General_Fingerprint & fp) const;
};

#endif
