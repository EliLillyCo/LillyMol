/*
  When doing lots of substructure searches, we often know that some
  of the queries are dependent on other queries.

  We introduce the notion of NEED and AVAILABLE.

  Some queries in the list do nothing but make AVAILABLE a value.
  Others NEED a given value

  Note that the individual queries make available IWString values.
  It is up to the array object to convert these to array items
*/

#include "Molecule_Lib/qry_wstats.h"

class NA_Substructure_Query : public Substructure_Hit_Statistics
{
  private:
    IWString _makes_available;

    resizable_array_p<IWString> _string_needs;
    resizable_array<int> _needs;

  public:
    NA_Substructure_Query ();

    int read (const const_IWSubstring &);
    int read (iwstring_data_source &);

    const IWString & makes_available () const { return _makes_available;}
    void set_makes_available (const IWString & s) { _makes_available = s;}

    resizable_array_p<IWString> & string_needs () { return _string_needs;}
    const resizable_array_p<IWString> & string_needs () const { return _string_needs;}
    const resizable_array<int> & needs () const { return _needs;}
    resizable_array<int> & needs () { return _needs;}

    int also_needs (const const_IWSubstring &);

    int has_preconditions () const { return _string_needs.number_elements ();}
};

/*
  Leading space stuff is significant
*/

#define NASS_AVAIL_TOKEN " AVAIL="

#define NASS_NEED_TOKEN "NEED="

class Set_of_NA_Substructure_Query: public resizable_array_p<NA_Substructure_Query>
{
  private:

//  During evaluation we evaluate any queries with no NEED values first. Do a one
//  time count of the number of such queries.

    int _queries_not_needing_anything;

//  When filtering for bad motifs, we can break once we've found a match to a
//  query. We ignore matches to queries which make AVAILABLE features.

    int _break_at_first_match;

//  Conversely for non matches

    int _break_at_first_non_match;

//  If we are using preconditions, we will probably want to not count
//  a match to a precondition as a match

    int _ignore_precondition_matches_for_break;

//  Everything below here is thread unsafe

//  We save the substructure results in these

    Substructure_Results * _sresults;

//  We save the level of optimisation we achieve - basically a count of the
//  number of searches we don't do

    int _searches_performed;
    int _optimisation_level;

// private functions

    void _count_queries_not_needing_anything ();

    int _build_from_query_file (const const_IWSubstring & c, int verbose);
    int _build_from_query_file (iwstring_data_source & input, const IWString & directory, int verbose);
    int _build_from_smarts_file (const const_IWSubstring & c, int verbose);
    int _build_from_smarts_file (iwstring_data_source & input, int verbose);
    int _build_from_smarts (const const_IWSubstring &, int);

    int _get_need_avail (NA_Substructure_Query * q, const const_IWSubstring & buffer);
    int _parse_avail_token (NA_Substructure_Query * q, const const_IWSubstring & buffer);
    int _parse_need_token (NA_Substructure_Query * q, const const_IWSubstring & buffer);

    int _find_circular_dependency (NA_Substructure_Query * q, int * already_hit);
    int _find_circular_dependencies (int * already_hit);

    int _substructure_search (Molecule & m, int * result);

    int _establish_preconditions (NA_Substructure_Query * q, Molecule_to_Match &, int * result);

  public:
    Set_of_NA_Substructure_Query ();
    ~Set_of_NA_Substructure_Query ();

    int  debug_print (std::ostream &) const;

    int  build_from_command_line (Command_Line & cl, char cflag, int verbose);

    int  build_from_option (const const_IWSubstring &, int);

    void set_ignore_precondition_matches_for_break (int s) { _ignore_precondition_matches_for_break = s;}

//  Call establish_cross_references once the queries have been read

    int  establish_cross_references ();

    void set_break_at_first_match (int s)     { _break_at_first_match = s;}
    void set_break_at_first_non_match (int s) { _break_at_first_non_match = s;}

    int  substructure_search (Molecule & m, int * result);

    const Substructure_Results * sresults () const { return _sresults;}

    const Substructure_Results & sresults (int i) { return _sresults[i];}

    int  optimisation_level () const { return _optimisation_level;}
    int  searches_performed () const { return _searches_performed;}
};
