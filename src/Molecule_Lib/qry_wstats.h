#ifndef MOLECULE_LIB_QRY_WSTATS_H_
#define MOLECULE_LIB_QRY_WSTATS_H_

#include <iostream>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

#include "ostream_and_type.h"
#include "substructure.h"
#include "Molecule_Lib/substructure.pb.h"

/*
  This extension of a substructure query keeps track of the number of
  hits, and can write matches (or non matches) to a file
*/

class Substructure_Hit_Statistics : public Substructure_Query
{
  private:
    ofstream_and_type _stream_for_matches;
    ofstream_and_type _stream_for_non_matches;
    int _append_match_details_to_molecule_name;
    int _append_non_match_details_to_molecule_name;

//  Sept 2010. The format of the matching info is awful. Allow a new format

    int _use_vertical_bars_for_query_details;

    // Feb 2022. CSR format for query matches. Queries that have names
    // with multiple tokens will be quoted.

    int _csr_query_details;
    // Defaults to space.
    IWString _csr_separator;

  protected:
    int _verbose;
    int _molecules_which_match;
    int _molecules_which_do_not_match;
    extending_resizable_array<int> _molecules_which_match_n_times;

//  private functions

  private:

    void _default_values ();
    int  _set_stream (FileType, const char *, ofstream_and_type &);
    int  _update_matches (int, Molecule *, const Substructure_Results &);
    int  _update_name_if_needed (int, Molecule *);

  public:
    Substructure_Hit_Statistics ();
    Substructure_Hit_Statistics (const const_IWSubstring &);
    ~Substructure_Hit_Statistics ();

    int ok () const;

    int set_verbose (int i) { return _verbose = i;}

    void set_use_vertical_bars_for_query_details (int s) { _use_vertical_bars_for_query_details = s;}
    void set_csr_query_details(int s) {
      _csr_query_details = s;
    }
    void set_csr_separator(const const_IWSubstring& s) {
      _csr_separator = s;
      _csr_query_details = 1;
    }

    int substructure_search (Molecule *);
    int substructure_search (Molecule & m) { return substructure_search(&m);}
    int substructure_search (Molecule *, Substructure_Results &);
    int substructure_search (Molecule &, Substructure_Results &);
    int substructure_search (Molecule_to_Match &);
    int substructure_search (Molecule_to_Match &, Substructure_Results &);

    int set_stream_for_matches     (FileType, const char *);
    int set_stream_for_non_matches (FileType, const char *);

    int molecules_which_match () const { return _molecules_which_match;}
    int molecules_which_do_not_match () const { return _molecules_which_do_not_match;}

    void set_append_match_details_to_molecule_name (int ii)
      {_append_match_details_to_molecule_name = ii;}

    void set_append_non_match_details_to_molecule_name (int ii)
      {_append_non_match_details_to_molecule_name = ii;}

    int report (std::ostream & os, int verbose) const;

    int ConstructFromProto(const SubstructureSearch::SubstructureQuery& proto);
};

extern std::ostream & operator << (std::ostream &, const Substructure_Hit_Statistics &);


#endif  // MOLECULE_LIB_QRY_WSTATS_H_
