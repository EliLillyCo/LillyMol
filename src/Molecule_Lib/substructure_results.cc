#include <stdlib.h>
#include <memory>
#include <algorithm>

#include "substructure.h"
#include "target.h"
#include "misc.h"
#include "misc2.h"

/*
  June 2010. We get undesirable behaviour when there are multiple fragments each
  with the same atom count and we are asked to remove matches not in the largest 
  fragment. The results may change depending on the fragment ordering.
  As an optional behaviour, we can have it keep track of all fragments with the
  largest atom count.
*/

static int remove_hits_not_in_largest_fragment_behaviour = 0;

void
set_remove_hits_not_in_largest_fragment_behaviour(int s)
{
  remove_hits_not_in_largest_fragment_behaviour = s;
}

Substructure_Results::Substructure_Results ()
{
  _save_matched_atoms = 1;

  _symmetry_class = NULL;

  _already_matched = NULL;

  _just_matched = NULL;

  _save_query_atoms_matched = 1;

  _remove_invalid_atom_numbers_from_new_embeddings = 0;

  _max_query_atoms_matched_in_search = 0;

  return;
}

Substructure_Results::~Substructure_Results ()
{
  if (NULL != _already_matched)
    delete [] _already_matched;

  if (NULL != _just_matched)
    delete [] _just_matched;

  return;
}

int
Substructure_Results::ok () const
{
  return 1;
}

int
Substructure_Results::return_code() const
{
  if (0 == _hits_found)
    return 0;
    
  if (! _save_matched_atoms)
    return _hits_found;

  return _embedding.number_elements();
}

int
Substructure_Results::copy_embeddings (const Substructure_Results & rhs)
{
  _symmetry_class = NULL;

  _already_matched = NULL;    // should we check to see if it needs to be deleted??

  _just_matched = NULL;

  _save_matched_atoms = rhs._save_matched_atoms;

  _save_query_atoms_matched = rhs._save_query_atoms_matched;

  _remove_invalid_atom_numbers_from_new_embeddings = rhs._remove_invalid_atom_numbers_from_new_embeddings;

  int ne = rhs._embedding.number_elements();

  _embedding.resize_keep_storage(0);

  _embedding.resize(ne);

  _query_atoms_matched.resize(0);

  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms * e = rhs._embedding[i];

    Set_of_Atoms * tmp = new Set_of_Atoms(*e);

    _embedding.add(tmp);
  }

  _hits_found = ne;

  return ne;
}

/*
  Note that we aren't saving the Query_Atoms_Matched. Add that feature
  if it ever becomes an issue
*/

int
Substructure_Results::add_embeddings (const Substructure_Results & rhs)
{
  int ne = rhs._embedding.number_elements();

  if (0 == ne)
    return 0;

  if (! _save_matched_atoms)   // huh!!
  {
    cerr << "Substructure_Results::add_embeddings: not saving embeddings\n";
    return 0;
  }

  _embedding.make_room_for_extra_items(ne);

  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms * e = rhs._embedding[i];

    Set_of_Atoms * tmp = new Set_of_Atoms(*e);

    _embedding.add(tmp);
  }

  _hits_found += ne;

  return 1;
}

int
Substructure_Results::set_save_matched_atoms (int s)
{
  _save_matched_atoms = s;

  return 1;
}

int
Substructure_Results::set_save_query_atoms_matched (int s)
{
  _save_query_atoms_matched = s;

  return 1;
}

/*
  Initialise is called when we are just about to start a search
*/

void
Substructure_Results::initialise(int m)
{
  reset();

  _atoms_in_target_molecule = m;

  if (_save_matched_atoms && _embedding.elements_allocated() < m)
  {
    _embedding.resize(m + 3);
    if (_save_query_atoms_matched)
      _query_atoms_matched.resize(m + 3);
  }

  _complete = 0;

  _max_query_atoms_matched_in_search = 0;

  return;
}

/*
  Reset is designed to remove the results of a previous search.
*/

void
Substructure_Results::reset ()
{
  _hits_found = 0;
  _embedding.resize_keep_storage(0);

  if (_save_query_atoms_matched)
    _query_atoms_matched.resize_keep_storage(0);

  _symmetry_class = NULL;

  if (NULL != _already_matched)
  {
    delete [] _already_matched;
    _already_matched = NULL;
  }

  if (NULL != _just_matched)
  {
    delete [] _just_matched;
    _just_matched = NULL;
  }

  _hits_per_fragment.resize_keep_storage(0);

  _embeddings_violating_distance_constraints = 0;

  _max_query_atoms_matched_in_search = 0;

  return;
}


int
Substructure_Results::embedding_is_unique(const Set_of_Atoms & new_embedding)
{
  if (NULL == _just_matched)
    _just_matched = new int[_atoms_in_target_molecule];

  std::fill_n(_just_matched, _atoms_in_target_molecule, 0);
  new_embedding.set_vector(_just_matched, 1);

  const int n = _embedding.number_elements();

  for (int i = 0; i < n; ++i)
  {
    if (_embedding[i]->all_members_non_zero_in_array(_just_matched))
      return 0;
  }

  return 1;
}

int 
Substructure_Results::embedding_overlaps_previous_embeddings (const Set_of_Atoms & new_embedding)
{
  if (NULL == _already_matched)
  {
    _already_matched = new_int(_atoms_in_target_molecule);
    new_embedding.set_vector(_already_matched, 1);
    return 0;
  }

  if (new_embedding.any_members_set_in_array(_already_matched))
    return 1;

  new_embedding.set_vector(_already_matched, 1);

  return 0;
}

int
Substructure_Results::add_embedding(Set_of_Atoms * e,
                                    const Query_Atoms_Matched & qam)
{
  _hits_found++;

  if (_hits_found > 9999 && 0 == _hits_found % 10000)
    cerr << "Query got " << _hits_found << " embeddings so far\n";

//cerr << "_save_matched_atoms " << _save_matched_atoms << " _save_query_atoms_matched " << _save_query_atoms_matched << endl;

  if (! _save_matched_atoms)
    return 1;

  if (_remove_invalid_atom_numbers_from_new_embeddings)
    e->remove_all(INVALID_ATOM_NUMBER);

  _embedding.add(e);

  if (_save_query_atoms_matched)
  {
    Query_Atoms_Matched * tmp = new Query_Atoms_Matched;

    int qam_size = qam.number_elements();

    tmp->resize(qam_size);
    for (int i = 0; i < qam_size; i++)
    {
      Substructure_Atom * a = qam[i];
      tmp->add(a);
      tmp->increment_preference_value(a->preference_value_current_match());
    }

    _query_atoms_matched.add(tmp);
  }

  return 1;
}

//#define DEBUG_EMBEDDING_IS_SYMMETRY_RELATED

int
Substructure_Results::_are_symmetry_related (const Set_of_Atoms & e1, const Set_of_Atoms & e2) const
{

  const int n = e1.number_elements();

  for (int i = 0; i < n; i++)
  {
    atom_number_t a1 = e1[i];
    atom_number_t a2 = e2[i];

    if (a1 < 0 || a2 < 0)
      continue;

#ifdef DEBUG_EMBEDDING_IS_SYMMETRY_RELATED
    cerr << " checking atom " << a1 << " (" << _symmetry_class[a1] << ") and " << a2 << " (" << _symmetry_class[a2] << ")\n";

#endif
    if (_symmetry_class[a1] != _symmetry_class[a2])
      return 0;
  }

#ifdef DEBUG_EMBEDDING_IS_SYMMETRY_RELATED
  cerr << "Embeddings are symmetry related\n";
#endif

  return 1;
}

/*
  Is a proposed new embedding symmetrically equivalent to an existing embedding
*/

int
Substructure_Results::embedding_is_symmetry_related (const Set_of_Atoms & new_embedding) const
{
  assert (_symmetry_class);

  int ne = _embedding.number_elements();

#ifdef DEBUG_EMBEDDING_IS_SYMMETRY_RELATED
  cerr << "Checking " << new_embedding << " against " << ne << " embeddings\n";
#endif

  if (0 == ne)
    return 0;

  int np = new_embedding.number_elements();

  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms * p = _embedding[i];

    if (np != p->number_elements())      // pretty unlikely to happen
      continue;

    if (_are_symmetry_related(new_embedding, *p))
      return 1;
  }

#ifdef DEBUG_EMBEDDING_IS_SYMMETRY_RELATED
  cerr << "New embedding not symmetry related\n";
#endif

  return 0;
}

/*
  The sorting process is complicated by the fact that two arrays must be
  sorted, the array of hits, and _query_atoms_matched.
  We create a little class to bond them together
  Note that this uses the implicit copy operator in iwqsort.
  The number of embeddings is expected to be small, so we
  don't worry too much about efficiency here.
*/

class QamSa
{
  private:
    Query_Atoms_Matched * _qam;
    Set_of_Atoms        * _sa;
    int                   _sort_value;
  public:
    QamSa (Query_Atoms_Matched * qam, Set_of_Atoms * s) : _qam (qam), _sa (s) {}
    QamSa () { _qam = NULL; _sa = NULL; _sort_value = 0;}

    void set_sort_value (int s) { _sort_value = s;}
    int  sort_value () const { return _sort_value;}

    Query_Atoms_Matched * qam () const { return _qam;}
    Set_of_Atoms * sa () const { return _sa;}

    void set_qam (Query_Atoms_Matched * q) { _qam = q;}
    void set_sa  (Set_of_Atoms * s)        { _sa  = s;}
};

static int
QamSa_comparitor (const QamSa * q1, const QamSa * q2)
{
  if (q1->sort_value() < q2->sort_value())
    return 1;
  else if (q1->sort_value() > q2->sort_value())
    return -1;
  else
    return 0;
}


//#define DEBUG_DO_SORT_BY_PREFERENCE_VALUE

int
Substructure_Results::sort_by_preference_value ()
{
  assert (_save_query_atoms_matched);

  int ne = _embedding.number_elements();

#ifdef DEBUG_DO_SORT_BY_PREFERENCE_VALUE
  cerr << "Before sorting " << ne << " embeddings\n";
  for (int i = 0; i < ne; i++)
  {
    Query_Atoms_Matched * qam = _query_atoms_matched[i];
    Set_of_Atoms * s = _embedding[i];

    cerr << "Embedding " << i << " preference value " << qam->preference_value() << 
            " atoms " << (*s) << endl;
  }
#endif

  QamSa * qamsa = new QamSa[ne]; std::unique_ptr<QamSa[]> free_qamsa(qamsa);

  for (int i = 0; i < ne; i++)
  {
    QamSa & t = qamsa[i];
    t.set_qam(_query_atoms_matched[i]);
    t.set_sa (_embedding[i]);
    t.set_sort_value(_query_atoms_matched[i]->preference_value());
//  cerr << "Preference " << t.sort_value() << endl;
  }

  ::qsort(qamsa, ne, sizeof(QamSa), (int (*) (const void *, const void *) ) QamSa_comparitor);

  for (int i = 0; i < ne; i++)
  {
    QamSa & t = qamsa[i];

    _embedding[i] = t.sa();
    _query_atoms_matched[i] = t.qam();
  }

#ifdef DEBUG_DO_SORT_BY_PREFERENCE_VALUE
  cerr << "After sorting\n";
  for (int i = 0; i < ne; i++)
  {
    Query_Atoms_Matched * qam = _query_atoms_matched[i];
    Set_of_Atoms * s = _embedding[i];

    cerr << "Embedding " << i << " preference value " << qam->preference_value() << 
            " atoms " << (*s) << endl;
  }
#endif

  return ne;
}

//#define DEBUG_REMOVE_HITS_VIOLATING_DISTANCE

/*
  Might look strange that I have a special case for 1 == ncheck, but this is 
  legacy of a bug. The original implementation just checked the first matched
  atom, and so rather than throwing that code away, I just made it a special
  case. Maybe a bad idea....
*/

int
Substructure_Results::remove_hits_violating_distance (Molecule_to_Match & target,
                              const Min_Max_Specifier<int> & distance_between_hits,
                              int ncheck)
{
  Molecule * m = target.molecule();

  if (1 != ncheck)
    return _remove_hits_violating_distance_check_all_atoms(*m, distance_between_hits, ncheck);

  int ne = _embedding.number_elements();

#ifdef DEBUG_REMOVE_HITS_VIOLATING_DISTANCE
  cerr << "Checking " << ne << " embeddings for distance violations\n";
#endif

  int rc = 0;    // the number of embeddings we remove

  for (int i = 0; i < ne; i++)
  {
    atom_number_t ai = _embedding[i]->item(0);

    for (int j = ne - 1; j > i; j--)
    {
      atom_number_t aj = _embedding[j]->item(0);

      if (m->fragment_membership(ai) != m->fragment_membership(aj))
        continue;

      int d = m->bonds_between(ai, aj);

#ifdef DEBUG_REMOVE_HITS_VIOLATING_DISTANCE
      cerr << "Hit " << i << " atom " << ai << " and hit " << j << " atom " << aj << 
              " are separated by " << d << " bonds, matches = " << distance_between_hits.matches(d) << endl;
#endif

      if (! distance_between_hits.matches(d))
      {
        _embedding.remove_item(j);
        if (_save_query_atoms_matched)
          _query_atoms_matched.remove_item(j);
        _embeddings_violating_distance_constraints++;
        j--;
        ne--;
        rc++;
      }
    }
  }

#ifdef DEBUG_REMOVE_HITS_VIOLATING_DISTANCE
  cerr << "Returning with " << ne << " valid matches\n";
#endif

  assert (ne == _embedding.number_elements());

  return rc;
}

static int
embeddings_too_close (Molecule & m,
                      const Set_of_Atoms & e1,
                      const Set_of_Atoms & e2,
                      const Min_Max_Specifier<int> & distance_between_hits,
                      int ncheck)
{
  int n1 = e1.number_elements();
  int n2 = e2.number_elements();

  if (n1 > ncheck)
    n1 = ncheck;

  if (n2 > ncheck)
    n2 = ncheck;

  for (int i = 0; i < n1; i++)
  {
    atom_number_t j = e1[i];

    for (int k = 0; k < n2; k++)
    {
      atom_number_t l = e2[k];

      int d;

      if (j == l)
        d = 0;
      else if (m.fragment_membership(j) != m.fragment_membership(l))
        continue;
      else
        d = m.bonds_between(j, l);

      if (! distance_between_hits.matches(d))
        return 1;
    }
  }

  return 0;
}

int
Substructure_Results::_remove_hits_violating_distance_check_all_atoms (Molecule & m,
                              const Min_Max_Specifier<int> & distance_between_hits,
                              int ncheck)
{
  int ne = _embedding.number_elements();

#ifdef DEBUG_REMOVE_HITS_VIOLATING_DISTANCE
  cerr << "Checking " << ne << " embeddings for distance violations\n";
#endif

  if (0 == ncheck)    // set to largest known embedding
  {
    for (int i = 0; i < ne; i++)
    {
      int s = _embedding[i]->number_elements();
      if (s > ncheck)
        ncheck = s;
    }
  }

  int rc = 0;    // the number of embeddings we remove

  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms * ei = _embedding[i];

    for (int j = ne - 1; j > i; j--)
    {
      const Set_of_Atoms * ej = _embedding[j];

//    cerr << " i = " << i << " ai " << ei->item(0) << " j = " << j << " aj " << ej->item(0) << endl;

      if (! embeddings_too_close(m, *ei, *ej, distance_between_hits, ncheck))
        continue;

      _embedding.remove_item(j);
      if (_save_query_atoms_matched)
        _query_atoms_matched.remove_item(j);
      _embeddings_violating_distance_constraints++;
      ne--;
      rc++;
    }
  }

#ifdef DEBUG_REMOVE_HITS_VIOLATING_DISTANCE
  cerr << "Returning with " << ne << " valid matches\n";
#endif

  assert (ne == _embedding.number_elements());

  return rc;
}

int
Substructure_Results::remove_hits_not_in_largest_fragment (Molecule_to_Match & target)
{
  Molecule * m = target.molecule();

  if (1 == remove_hits_not_in_largest_fragment_behaviour)
    return _remove_hits_not_in_largest_fragment_multiple_largest(*m);

  int nf = m->number_fragments();
  assert (nf > 1);

  int largest_fragment = 0;
  int atoms_in_largest_fragment = m->atoms_in_fragment(0);

  for (int i = 1; i < nf; i++)
  {
    int a = m->atoms_in_fragment(i);
    if (a > atoms_in_largest_fragment)
    {
      atoms_in_largest_fragment = a;
      largest_fragment = i;
    }
  }

//cerr << "Target has " << nf << " fragments, largest " << atoms_in_largest_fragment << ", checking " << _embedding.number_elements() << " embeddings\n";

  int rc = 0;

  for (int i = _embedding.number_elements() - 1; i >= 0; i--)
  {
    const Set_of_Atoms * e = _embedding[i];
    atom_number_t j = e->item(0);
    if (m->fragment_membership(j) != largest_fragment)
    {
      _embedding.remove_item(i);
      rc++;
    }
  }

  return rc;
}

int
Substructure_Results::_remove_hits_not_in_largest_fragment_multiple_largest (Molecule & m)
{
  int nf = m.number_fragments();
  assert (nf > 1);

  int * is_largest = new_int(nf); std::unique_ptr<int[]> free_is_largest(is_largest);

  int atoms_in_largest_fragment = m.atoms_in_fragment(0);
  is_largest[0] = 1;

  for (int i = 1; i < nf; i++)
  {
    int a = m.atoms_in_fragment(i);

    if (a < atoms_in_largest_fragment)
      continue;

    if (a > atoms_in_largest_fragment)
    {
      set_vector(is_largest, nf, 0);
      is_largest[i] = 1;
      atoms_in_largest_fragment = a;
    }
    else
      is_largest[i] = 1;
  }

#ifdef DEBUG_REMOVE_HITS_NOT_IN_LARGEST_FRAGMENT_MULTIPLE_LARGEST
  cerr << atoms_in_largest_fragment << " atoms in largest fragment\n";
  for (int i = 0; i < nf; i++)
  {
    cerr << "frag " << i << " largest " << is_largest[i] << endl;
  }
#endif

  int rc = 0;

  for (int i = _embedding.number_elements() - 1; i >= 0; i--)
  {
    const Set_of_Atoms * e = _embedding[i];

    atom_number_t j = e->item(0);

#ifdef DEBUG_REMOVE_HITS_NOT_IN_LARGEST_FRAGMENT_MULTIPLE_LARGEST
    cerr << "Embedding " << i << " in frag " << m.fragment_membership(j) << endl;
#endif

    if (0 == is_largest[m.fragment_membership(j)])
    {
      _embedding.remove_item(i);
      rc++;
    }
  }

#ifdef DEBUG_REMOVE_HITS_NOT_IN_LARGEST_FRAGMENT_MULTIPLE_LARGEST
  cerr << "Embedding contains " << _embedding.number_elements() << endl;
#endif

  return rc;
}

int
Substructure_Results::print_embeddings (std::ostream & os, int print_query_atoms) const
{
  assert (ok());
  assert (os.good ());

  int ne = _embedding.number_elements();
  os << "Query contains " << ne << " embeddings\n";
  for (int i = 0; i < ne; i++)
  {
    const Set_of_Atoms & p = *(_embedding[i]);
    os << "Embedding " << i << ": " << p;
    os << endl;
    if (print_query_atoms && _save_query_atoms_matched)
    {
      os << "Query atoms:      ";
      const Query_Atoms_Matched & q = *(_query_atoms_matched[i]);
      for (int j = 0; j < q.number_elements(); j++)
      {
        os << " " << q[j]->unique_id();
      }

      os << endl;
    }
  }

  return os.good();
}

/*const Substructure_Atom *
Substructure_Results::query_atom_matching  (int i, int j) const
{
  const Query_Atoms_Matched * m = query_atoms_matching(i);

  return m->item(j);
}*/

/*int
Substructure_Results::query_atoms_in_match (int i) const
{
  const Set_of_Atoms * p = embedding(i);

  return p->number_elements();
}*/

//#define DEBUG_REMOVE_LOW_PREFERENCE_HITS

int
Substructure_Results::remove_low_preference_hits ()
{
  assert (_save_query_atoms_matched);

  int min_preference = 0;
  int max_preference = 0;

  int ne = _embedding.number_elements();
  for (int i = 0; i < ne; i++)
  {
    Query_Atoms_Matched * qam = _query_atoms_matched[i];

    int p = qam->preference_value();

#ifdef DEBUG_REMOVE_LOW_PREFERENCE_HITS
    cerr << "Embedding " << i << " has preference value " << p << ' ' << *(_embedding[i]) << endl;
#endif

    if (0 == i)
    {
      min_preference = p;
      max_preference = p;
    }
    else
    {
      if (p < min_preference)
        min_preference = p;
      if (p > max_preference)
        max_preference = p;
    }
  }

#ifdef DEBUG_REMOVE_LOW_PREFERENCE_HITS
  if (min_preference == max_preference)
    cerr << "All preference values the same, not reduced\n";
#endif

  if (min_preference == max_preference)
    return _embedding.number_elements();

  for (int i = ne - 1; i >= 0; i--)
  {
    Query_Atoms_Matched * qam = _query_atoms_matched[i];
    if (qam->preference_value() < max_preference)
    {
      _embedding.remove_item(i);
      _query_atoms_matched.remove_item(i);
    }
  }

#ifdef DEBUG_REMOVE_LOW_PREFERENCE_HITS
  cerr << "Reduced to " << _embedding.number_elements() << " embeddings\n";
#endif

  return _query_atoms_matched.number_elements();;
}

int
Substructure_Results::remove_embedding (int e)
{
  assert (_hits_found> 0);

  assert (_embedding.ok_index(e));

  _embedding.remove_item(e);
  if (_query_atoms_matched.number_elements())
    _query_atoms_matched.remove_item(e);

  _hits_found--;

  return 1;
}

void
Substructure_Results::size_hits_per_fragment_array (int s)
{
  if (_hits_per_fragment.number_elements())
    _hits_per_fragment.resize(0);

  _hits_per_fragment.extend(s);

  return;
}

void
Substructure_Results::each_embedding_set_vector (int * v,
                                                 int s) const
{
  for (int i = 0; i < _embedding.number_elements(); i++)
  {
    _embedding[i]->set_vector(v, s);
  }

  return;
}

/*
  Max distance from a matched atom to an unmatched atom.
  Would be much faster with an array
*/

static int
compute_maxd (Molecule & m,
              const Set_of_Atoms & s)
{
  int n = s.number_elements();

  int matoms = m.natoms();

  int maxd = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t a = s[i];

    for (int j = 0; j < matoms; j++)
    {
      if (a == j || s.contains(j))
        continue;

      int d = m.bonds_between(a, j);

      if (d > maxd)
        maxd = d;
    }
  }

  return maxd;
}

static int
compute_ncon (const Molecule & m,
              const Set_of_Atoms & s)
{
  int n = s.number_elements();

  int rc = 0;

  for (int i = 0; i < n; i++)
  {
    atom_number_t a = s[i];

    const Atom * ai = m.atomi(a);

    int acon = ai->ncon();

    for (int j = 0; j < acon; j++)
    {
      atom_number_t k = ai->other(a, j);

      if (! s.contains(k))
        rc++;
    }
  }

  return rc;
}

int
Substructure_Results::sort_matches (Molecule_to_Match & target,
                                    int sort_specification)
{
  int ne = _embedding.number_elements();

//cerr << "Sorting " << ne << " embeddings\n";

  QamSa * qamsa = new QamSa[ne]; std::unique_ptr<QamSa[]> free_qamsa(qamsa);  // the sorter is set up to sort high to low, so we need to be careful how we load things here

  for (int i = 0; i < ne; i++)
  {
    QamSa & t = qamsa[i];
    t.set_qam(_query_atoms_matched[i]);
    t.set_sa (_embedding[i]);

    int s = 0;

    if ((SORT_MATCHES_BY_INCREASING_NCON & sort_specification) ||
        (SORT_MATCHES_BY_DECREASING_NCON & sort_specification))
    {
      int ncon = compute_ncon(*(target.molecule()), *(_embedding[i]));
      if (SORT_MATCHES_BY_INCREASING_NCON & sort_specification)
        s = -ncon;
      else
        s = ncon;
    }

    if ((SORT_MATCHES_BY_INCREASING_MAXD & sort_specification) ||
        (SORT_MATCHES_BY_DECREASING_MAXD & sort_specification))
    {
      int maxd = compute_maxd(*(target.molecule()), *(_embedding[i]));
      if (SORT_MATCHES_BY_INCREASING_MAXD & sort_specification)
        s -= 1000 * maxd;    // might break with very large molecules, 1000 is arbitrary
      else
        s += 1000 * maxd;

//    cerr << *(_embedding[i]) << " maxd " << maxd << endl;
    }

//  cerr << "Sort value " << s << endl;
    t.set_sort_value(s);
  }

  ::qsort(qamsa, ne, sizeof(QamSa), (int (*) (const void *, const void *) ) QamSa_comparitor);

  for (int i = 0; i < ne; i++)
  {
    QamSa & t = qamsa[i];

    _embedding[i] = t.sa();
//  cerr << "Now " << *(_embedding[i]) << " sort value " << t.sort_value() << endl;
    _query_atoms_matched[i] = t.qam();
  }

  return 1;
}

int
Substructure_Results::set_embeddings(Set_of_Atoms const** e, const int s)
{
  assert (_embedding.number_elements() == s);

  std::copy_n(const_cast<Set_of_Atoms**>(e), s, _embedding.rawdata());    // very dangerous const case since we are giving those Set_of_Atoms objects to _embedding, which will delete them!

  return 1;
}

int
Substructure_Results::overlaps_with_existing_embedding(const Set_of_Atoms & s)
{
  cerr << "Substructure_Results::overlaps_with_existing_embedding:checking " << s << ", atoms_in_target_molecule " << _atoms_in_target_molecule << endl;

  if (NULL == _already_matched)
  {
    _already_matched = new_int(_atoms_in_target_molecule);
    s.set_vector(_already_matched, 1);
    return 0;
  }

  if (s.any_members_set_in_array(_already_matched))
  {
    cerr << "Hits previously matched atoms\n";
    for (int i = 0; i < _atoms_in_target_molecule; ++i)
    {
      cerr << " matched atom " << i << " " << _already_matched[i] << endl;
    }
  }

  if (s.any_members_set_in_array(_already_matched))
    return 1;

  s.set_vector(_already_matched, 1);

  return 0;
}
