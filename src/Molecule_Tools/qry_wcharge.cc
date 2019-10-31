
#include "iwaray.h"
#include "qry_wcharge.h"
#include "path.h"

/*
  We must have the same number of charge accumulators as there are query atoms
  in the query.
*/

void
Query_and_Charge_Stats::_default_values ()
{
  _embeddings = 0;
  _max_atoms = 0;

  return;
}

Query_and_Charge_Stats::Query_and_Charge_Stats (const const_IWSubstring & comment) :
                               Substructure_Query (comment)
{
  _default_values ();
}

Query_and_Charge_Stats::Query_and_Charge_Stats ()
{
  _default_values ();
}

Query_and_Charge_Stats::~Query_and_Charge_Stats ()
{
  assert (ok ());

  if (-6 == _embeddings)
    cerr << "freeing previously freed Query_and_Charge_Stats\n";

  _embeddings = -6;

  return;
}

int
Query_and_Charge_Stats::_determine_max_atoms ()
{
  for (int i = 0; i < _number_elements; i++)
  {
    Single_Substructure_Query * qi = _things[i];
    int tmp = qi->max_atoms_in_query ();
    if (tmp > _max_atoms)
      _max_atoms = tmp;
  }
  
  return 1;
}

/*
  This was initially developed when a Substructure_Query was a single
  query. When Substructure_Queries became composite objects, I just
  made simple changes to get it to compile. I'm really not sure how
  usable this would be with multiple queries.
*/

int
Query_and_Charge_Stats::_allocate_charge_distributions ()
{
  if (0 == _max_atoms)
    _determine_max_atoms ();

  if (_charge_distributions.number_elements () < _max_atoms)
    _charge_distributions.resize (_max_atoms);

  for (int i = _charge_distributions.number_elements (); i < _max_atoms; i++)
  {
    Charge_Distribution * tmp = new Charge_Distribution;
    _charge_distributions.add (tmp);
  }

  return 1;
}

int
Query_and_Charge_Stats::tally_embedding (const Molecule & m,
                                         Substructure_Results & sresults)
{
  if (0 == _max_atoms)
    _determine_max_atoms ();

  _embeddings++;

  if (_charge_distributions.number_elements () < _max_atoms)
    _allocate_charge_distributions ();

  const Set_of_Atoms * s = sresults.embedding (0);

  int atoms_matched = s->number_elements ();

  for (int i = 0; i < atoms_matched; i++)
  {
    atom_number_t a = s->item (i);

    charge_t q = m.charge_on_atom (a);
    _charge_distributions[i]->extra (q);
  }

  return 1;
}

#include <iomanip>
using namespace std;

int
Query_and_Charge_Stats::report (int query_number, std::ostream & os) const
{
  assert (os.good ());

  assert (ok ());
  assert (_charge_distributions.ok ());

  double largest_range = 0.0;
  int    atom_with_largest_range = -6;
  double largest_variance = 0.0;
  int    atom_with_largest_variance = -6;

  assert (_charge_distributions.number_elements ());

  os << "Report on Query Charge Distribution number " << query_number << ", "
     << _embeddings << " embeddings\n";

  for (int i = 0; i < _charge_distributions.number_elements (); i++)
  {
    Charge_Distribution * dist = _charge_distributions[i];
    assert (_embeddings == dist->n ());

    os << "Atom " << setw (3) << i;
    os << setprecision (3) << " min " << setw (8) << dist->minval () << 
                              " max " << setw (8) << dist->maxval ();
    if (_embeddings > 1)
    {
      double tmp = dist->range ();
      if (tmp > largest_range)
      {
        largest_range = tmp;
        atom_with_largest_range = i;
      }
      os << setprecision (3) << " range " << setw (8) << tmp;

      os << setprecision (3) << " ave " << setw (8) << dist->average ();
    
      tmp = dist->variance ();
      if (tmp > largest_variance)
      {
        largest_variance = tmp;
        atom_with_largest_variance = i;
      }
      os << setprecision (3) << " var " << setw (8) << tmp;
    }
    os << endl;
  }

  if (_embeddings > 1)
  {
    os << "Atom " << setw(3) << atom_with_largest_range <<    " has the largest range     " << 
          setprecision (3) << setw (8) << largest_range << endl;
    os << "Atom " <<setw (3) << atom_with_largest_variance << " has the largest variiance " << 
          setprecision (3) << setw (8) << largest_variance << endl;
  }

  return 1;
}
