#include <stdlib.h>
#include <assert.h>

#if defined (IWSGI) || (__GNUC_MINOR__ == 95)
#include <strstream>
#else
#include <sstream>
#endif

#include <iostream>

#include "qry_wstats.h"
#include "target.h"

void

Substructure_Hit_Statistics::_default_values ()
{
  _verbose = 0;

  _molecules_which_match = 0;
  _molecules_which_do_not_match = 0;

   _append_match_details_to_molecule_name = 0;
   _append_non_match_details_to_molecule_name = 0;

   _use_vertical_bars_for_query_details = 0;

  return;
}

Substructure_Hit_Statistics::Substructure_Hit_Statistics ()
{
  _default_values ();

  return;
}

/*Substructure_Hit_Statistics::Substructure_Hit_Statistics (const IWString & comment) :
        Substructure_Query (comment)
{
  _default_values ();

  return;
}*/

Substructure_Hit_Statistics::Substructure_Hit_Statistics (const const_IWSubstring & comment) :
        Substructure_Query (comment)
{
  _default_values ();

  return;
}

Substructure_Hit_Statistics::~Substructure_Hit_Statistics ()
{

  if (-9023 == _molecules_which_match)
    cerr << "Substructure_Hit_Statistics: freeing already freed object\n";

  _molecules_which_match = -9023;

  return;
}

int
Substructure_Hit_Statistics::ok () const
{
  return Substructure_Query::ok ();
}

int
Substructure_Hit_Statistics::report (std::ostream & os, int verbose) const
{
  assert (os.good ());

  os << "Details on hits for query '" << Substructure_Query::comment () << "'";
  if (_verbose)
    os << " (verbose)";
  os << endl;

  if (0 == _molecules_which_match && 0 == _molecules_which_do_not_match)
  {
    os << "No match attempts made\n";
    return 1;
  }

  float fraction = static_cast<float>(_molecules_which_match) / static_cast<float>(_molecules_which_match + _molecules_which_do_not_match);

  os << "Tested " <<  _molecules_which_match + _molecules_which_do_not_match <<
        " molecules, " << _molecules_which_match << " matched, " <<
        _molecules_which_do_not_match << " did not. Fraction " << fraction << endl;

  int sites = 0;
  int number_different_contributions = 0;
  for (int i = 0; i < _molecules_which_match_n_times.number_elements (); i++)
  {
    if (_molecules_which_match_n_times.item (i))
    {
      os << " " << _molecules_which_match_n_times[i] << " molecules had " << i << " hits\n";
      sites += i * _molecules_which_match_n_times[i];
      number_different_contributions++;
    }
  }

  if (number_different_contributions)
    os << sites << " reactive sites\n";
 
  if (_molecules_which_match && _stream_for_matches.valid ())
    os << " Matches written to '" << _stream_for_matches.fname () << "'\n";

  if (_molecules_which_do_not_match && _stream_for_non_matches.valid ())
    os << " Non matches written to '" << _stream_for_non_matches.fname () << "'\n";

  if (verbose)
    print_environment_matches (os);

  return os.good ();
}

int
Substructure_Hit_Statistics::_set_stream (int output_type, const char * fname,
                                          ofstream_and_type & stream_for)
{
  if (! valid_file_type (output_type))
  {
    cerr << "Substructure_Hit_Statistics::_set_stream: invalid type " << output_type << endl;
    return 0;
  }

  if (! stream_for.open (fname))
  {
    cerr << "Substructure_Hit_Statistics::_set_stream: cannot open '" << fname << "'\n";
    return 0;
  }

  stream_for.set_type (output_type);

  return 1;
}

int
Substructure_Hit_Statistics::set_stream_for_matches (int output_type, const char * fname)
{
  assert (_stream_for_matches.good ());

  if (_stream_for_non_matches.valid () && 
      _stream_for_non_matches.fname () == fname)
  {
    cerr << "Substructure_Hit_Statistics::set_stream_for_matches: name collision '" << fname << "'\n";
    return 0;
  }

  if (! _set_stream (output_type, fname, _stream_for_matches))
  {
    cerr << "Substructure_Hit_Statistics::set_stream_for_matches: cannot use '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Substructure_Hit_Statistics::set_stream_for_non_matches (int output_type, const char * fname)
{
  assert (_stream_for_matches.good ());

  if (_stream_for_matches.valid () && 
      _stream_for_matches.fname () == fname)
  {
    cerr << "Substructure_Hit_Statistics::set_stream_for_non_matches: name collision '" << fname << "'\n";
    return 0;
  }

  if (! _set_stream (output_type, fname, _stream_for_non_matches))
  {
    cerr << "Substructure_Hit_Statistics::set_stream_for_non_matches: cannot use '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Substructure_Hit_Statistics::_update_name_if_needed (int nmatches, Molecule * m)
{
  if (nmatches && _append_match_details_to_molecule_name && nmatches >= _append_match_details_to_molecule_name)
  {
    IWString new_name;
    new_name.resize (512);    // void excessive mallocing

    new_name = m->molecule_name ();

    if (_use_vertical_bars_for_query_details)
    {
      if (! new_name.ends_with('|'))
        new_name << " |";
      new_name << nmatches << " matches " << Substructure_Query::comment () << "|";
    }
    else
    {
      new_name << " (" << nmatches << " matches to '" << Substructure_Query::comment () << "')";
    }

    m->set_name (new_name);
  }

  if (0 == nmatches && _append_non_match_details_to_molecule_name)
  {
    IWString new_name = m->molecule_name ();

    new_name << " (0 matches to '" << Substructure_Query::comment () << "')";

    m->set_name (new_name);
  }

  return 1;
}

int
Substructure_Hit_Statistics::substructure_search (Molecule * m)
{
  Substructure_Results sresults;

  int nmatches = Substructure_Query::substructure_search (m, sresults);

  _update_matches (nmatches, m, sresults);

  _update_name_if_needed (nmatches, m);

  return nmatches;
}

int 
Substructure_Hit_Statistics::substructure_search (Molecule * m,
                                                  Substructure_Results & sresults)
{
  return substructure_search (*m, sresults);
}

int
Substructure_Hit_Statistics::substructure_search (Molecule & m,
                                                  Substructure_Results & results)
{
  Molecule_to_Match target (&m);

  return substructure_search (target, results);
}

int
Substructure_Hit_Statistics::substructure_search (Molecule_to_Match & m)
{
  Substructure_Results not_used;

  return substructure_search (m, not_used);
}

int
Substructure_Hit_Statistics::substructure_search (Molecule_to_Match & m,
                                                  Substructure_Results & results)
{
  int nmatches = Substructure_Query::substructure_search (m, results);

  _update_matches (nmatches, m.molecule (), results);

  _update_name_if_needed (nmatches, m.molecule ());

  return nmatches;
}

/*
  By one of our substructure_search interfaces, we have NMATCHES hits for
  molecule M
*/

int
Substructure_Hit_Statistics::_update_matches (int nmatches,
                                              Molecule * m,
                                              const Substructure_Results & sresults)
{

  if (0 == nmatches)
  {
    _molecules_which_do_not_match++;
    if (_stream_for_non_matches.valid ())
    {
      if (1 == _molecules_which_do_not_match)
        _stream_for_non_matches << "# Non matches to query '" << comment () << "'\n";
      _stream_for_non_matches.write_molecule (m);
    }

    if (_verbose > 1)
    {
      cerr << "no match to '" << _comment << "' ";
      for (int i = 0; i < _number_elements; i++)
      {
//      const Single_Substructure_Query * qi = _things[i];
//      assert (qi->ok ());

        if (_number_elements > 1)
          cerr << "query " << i;
        cerr << " only matched " << sresults.max_query_atoms_matched_in_search () <<
                " atoms in the query\n";
      }
    }

    return 0;
  }

  _molecules_which_match++;

  if (_stream_for_matches.valid ())
  {
//  if (1 == _molecules_which_match)
//    _stream_for_matches << "# Matches to query '" << comment () << "'\n";

    _stream_for_matches.write_molecule (m);
  }

  _molecules_which_match_n_times[nmatches]++;

  return nmatches;
}


