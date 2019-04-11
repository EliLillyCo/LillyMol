#ifndef FB_BITSNWEIGHTS_H
#define FB_BITSNWEIGHTS_H

class const_IWSubstring;

/*
  We want to have a weight for each bit in a set of bits
*/

class FB_Bits_and_Weights_Fixed_Width
{
  private:
    int _nbits;
    double * _wmatch;
    double * _wnmatch;

    int _include_non_matching_bits_in_similarity;

//  private functions

    int  _parse_weight_record (const const_IWSubstring & buffer, int ndx);

    double _tanimoto_only_consider_matching_bits (const int * fp1, const int * fp2) const;
    double _tanimoto_include_non_matching_bits (const int * fp1, const int * fp2) const;

  public:
    FB_Bits_and_Weights_Fixed_Width();
    ~FB_Bits_and_Weights_Fixed_Width();

    int build_from_file (const const_IWSubstring &);
    int build_from_file (iwstring_data_source &);
    int construct_from_command_line (Command_Line & cl,
                                                char flag,
                                                int verbose);

    int active () const { return _nbits;}

    double tanimoto (const int * fp1, const int * fp2) const;
};

#endif
