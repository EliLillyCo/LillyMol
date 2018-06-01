#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "misc.h"

#include "substructure.h"
#include "target.h"
#include "path.h"

Substructure_Ring_Specification::Substructure_Ring_Specification()
{
  _aromatic = SUBSTRUCTURE_NOT_SPECIFIED;

  return;
}

Substructure_Ring_Specification::~Substructure_Ring_Specification()
{
  if (-5 == _match_as_match_or_rejection)
    cerr << "Deleting already deleted Substructure_Ring_Specification\n";

  _match_as_match_or_rejection = -5;

  return;
}

int
Substructure_Ring_Specification::ok() const
{
  if (_match_as_match_or_rejection < 0)
    return 0;

  return 1;
}

int
Substructure_Ring_Specification::debug_print (std::ostream & os,
                                              const IWString & indentation) const
{
  IWString ind = indentation;
  ind += "  ";

  os << ind << "Details on Substructure_Ring_Specification\n";

  if (0 == _match_as_match_or_rejection)
    os << ", REJECTION";

  if (SUBSTRUCTURE_NOT_SPECIFIED == _aromatic)
    ;
  else if (_aromatic)
    os << ind << "  Aromatic";
  else
    os << ind << "  Non aromatic";

  os << endl;

  if (_comment.length())
    os << ind << "  " << _comment << endl;

  os << ind << "  Hits needed " << _hits_needed << endl;
  os << ind << "  Ring size   " << _ring_size << endl;
  os << ind << "  Ncon        " << _ncon << endl;
  os << ind << "  Unsat       " << _within_ring_unsaturation << endl;
  os << ind << "  AtwPi e's   " << _atoms_with_pi_electrons << endl;
  os << ind << "  Fused       " << _fused << endl;
  os << ind << "  Heteroatoms " << _heteroatom_count << endl;
  os << ind << "  Arom Nbrs   " << _fused_aromatic_neighbours << endl;
  os << ind << "  Aliph Nbrs  " << _fused_non_aromatic_neighbours << endl;
  os << ind << "  Max Fused N " << _largest_number_of_bonds_shared_with_another_ring << endl;
  os << ind << "  Strongly Fs " << _strongly_fused_ring_neighbours << endl;

  return os.good();
}

int
Substructure_Ring_Specification::terse_details (std::ostream & os,
                                      const IWString & indentation) const
{
  IWString ind = indentation;
  ind += "  ";

  os << ind << "Details on Substructure_Ring_Specification\n";

  if (0 == _match_as_match_or_rejection)
    os << ", REJECTION";

  if (SUBSTRUCTURE_NOT_SPECIFIED == _aromatic)
    ;
  else if (_aromatic)
    os << ind << "  Aromatic";
  else
    os << ind << "  Non aromatic";

  os << endl;

  if (_hits_needed.is_set())
    os << ind << "  Hits needed " << _hits_needed << endl;
  if (_ring_size.is_set())
    os << ind << "  Ring size   " << _ring_size << endl;
  if (_ncon.is_set())
    os << ind << "  Ncon        " << _ncon << endl;
  if (_within_ring_unsaturation.is_set())
    os << ind << "  Unsat       " << _within_ring_unsaturation << endl;
  if (_atoms_with_pi_electrons.is_set())
    os << ind << "  AtwPie's    " << _atoms_with_pi_electrons << endl;
  if (_fused.is_set())
    os << ind << "  Fused       " << _fused << endl;
  if (_heteroatom_count.is_set())
    os << ind << "  Heteroatoms " << _heteroatom_count << endl;
  if (_fused_aromatic_neighbours.is_set())
    os << ind << "  Arom Nbrs   " << _fused_aromatic_neighbours << endl;
  if (_fused_non_aromatic_neighbours.is_set())
    os << ind << "  Aliph Nbrs  " << _fused_aromatic_neighbours << endl;

  return os.good();
}

static int
compute_atoms_with_pi_electrons (const Ring * r, Molecule_to_Match & target)
{
  int rc = 0;

  int rsize = r->number_elements();

  for (int i = 0; i < rsize; i++)
  {
    atom_number_t j = r->item(i);
  
    Target_Atom & a = target[j];
  
    int acon = a.ncon();

    if (a.nbonds() > acon)     // unsaturation so we assume pi electrons
      rc++;
    else
    {
      Atom * a1 = const_cast<Atom *>(a.atom());

      int pi;
      if (a1->pi_electrons(pi) && pi > 0)
        rc++;
    }
  }

  cerr << "Ring has " << rc << " of " << rsize << " atoms with pi electrons\n";

  return rc;
}

static int
compute_within_ring_unsaturation (const Ring * r, Molecule_to_Match & target)
{
  int rc = 0;

  int rsize = r->number_elements();

  for (int i = 0; i < rsize; i++)
  {
    atom_number_t j = r->item(i);
  
    Target_Atom & a = target[j];
  
    int acon = a.ncon();

    if (a.nbonds() == acon)     // no unsaturation here
      continue;

//  We must make sure that the unsaturation is within the ring

    for (int k = 0; k < acon; k++)
    {
      Bond_and_Target_Atom & bata = a.other(k);
      if (bata.bond()->is_single_bond())
        continue;

      atom_number_t o = bata.other()->atom_number();
      if (o < j)     
        continue;

      if (r->contains(o))     // multiple bond to another ring atom
        rc++;
    }
  }

  return rc;
}

//#define DEBUG_SS_RING_MATCHES

int
Substructure_Ring_Specification::matches (Molecule_to_Match & target)
{
  if (SUBSTRUCTURE_NOT_SPECIFIED != _aromatic)
    target.molecule()->compute_aromaticity();

  int nr = target.nrings();

  extending_resizable_array<int> hits_in_fragment;    // used if _all_hits_in_same_fragment specified
  Molecule * m;
  if (_all_hits_in_same_fragment)
    m = target.molecule();
  else
    m = NULL;

  int nhits = 0;

  for (int i = 0; i < nr; i++)
  {
    const Ring * r = target.ringi(i);

    int rsize = r->number_elements();

#ifdef DEBUG_SS_RING_MATCHES
    cerr << "Testing ring " << i << " of size " << rsize << " matches is " <<
          _ring_size.matches(rsize) << endl;
#endif

    if (! _ring_size.matches(rsize))
      continue;

    if (! _fused.matches(r->fused_ring_neighbours()))
      continue;

    if (! _largest_number_of_bonds_shared_with_another_ring.matches(r->largest_number_of_bonds_shared_with_another_ring()))
      continue;

    if (! _strongly_fused_ring_neighbours.matches(r->strongly_fused_ring_neighbours()))
      continue;

    if (SUBSTRUCTURE_NOT_SPECIFIED == _aromatic)
      ;
    else if (0 == _aromatic && r->is_aromatic())
      continue;
    else if (_aromatic && r->is_non_aromatic())
      continue;

    if (_ncon.is_set() || _attached_heteroatom_count.is_set() || _heteroatom_count.is_set())
    {
      int rcon   = 0;
      int rh     = 0;
      int ahc    = 0;
      for (int j = 0; j < rsize; j++)
      {
        atom_number_t k = r->item(j);
  
        Target_Atom & a = target[k];
  
        if (_is_heteroatom[a.atomic_number()])
          rh++;
  
        int acon = a.ncon();

        if (acon > 2)        // if only 2 connections, all nbrs in ring
          rcon += acon - 2;    // two of its neighbours must be in the ring

        if (_attached_heteroatom_count.is_set())    // only compute if needed
        {
          for (int l = 0; l < acon; l++)
          {
            const Bond_and_Target_Atom & bata = a.other(l);

            const Target_Atom * n = bata.other();

            if (r->contains(n->atom_number()))
              continue;

            if (_is_heteroatom[n->atomic_number()])
              ahc++;
          }
        }
      }

      if (! _ncon.matches(rcon))
        continue;
  
      if (! _attached_heteroatom_count.matches(ahc))
        continue;
  
#ifdef DEBUG_SS_RING_MATCHES
      cerr << "rh = " << rh << " match " << _heteroatom_count.matches(rh) << endl;
#endif

      if (! _heteroatom_count.matches(rh))
        continue;
    }

    if (_fused_aromatic_neighbours.is_set())
    {
      int arfsn = 0;
      for (int j = 0; j < r->fused_ring_neighbours(); j++)
      {
        const Ring * rnj = r->fused_neighbour(j);
        if (rnj->is_aromatic())
          arfsn++;
      }

      if (! _fused_aromatic_neighbours.matches(arfsn))
        continue;
    }

    if (_fused_non_aromatic_neighbours.is_set())
    {
      int narfsn = 0;
      for (int j = 0; j < r->fused_ring_neighbours(); j++)
      {
        const Ring * rnj = r->fused_neighbour(j);
        if (rnj->is_non_aromatic())
          narfsn++;
      }

      if (! _fused_non_aromatic_neighbours.matches(narfsn))
        continue;
    }

    if (_all_hits_in_same_fragment)
    {
      atom_number_t a = r->item(0);
      hits_in_fragment[m->fragment_membership(a)]++;
    }

    if (_within_ring_unsaturation.is_set())
    {
      int ring_unsaturation = compute_within_ring_unsaturation(r, target);

      if (! _within_ring_unsaturation.matches(ring_unsaturation))
        continue;
    }

    if (_atoms_with_pi_electrons.is_set())
    {
      int awpe = compute_atoms_with_pi_electrons(r, target);
      if (! _atoms_with_pi_electrons.matches(awpe))
        continue;
    }

    if (_environment_atom.number_elements() && ! _environment_matches(target, *r))
      continue;

    nhits++;
  }

  if (0 == nhits)
    return ! _match_as_match_or_rejection;

  if (_all_hits_in_same_fragment)
  {
    for (int i = 0; i < hits_in_fragment.number_elements(); i++)
    {
      if (_hits_needed.matches(hits_in_fragment[i]))
        return _match_as_match_or_rejection;
    }

    return ! _match_as_match_or_rejection;
  }

  if (_hits_needed.matches(nhits))
    return _match_as_match_or_rejection;

  if (_match_as_match_or_rejection)
    return 0;
  else
    return nhits;
}

int
Substructure_Ring_Specification::_environment_matches (Molecule_to_Match & target,
                                   const Ring & r)
{
  int * tmp = new_int(target.natoms()); std::unique_ptr<int[]> free_tmp(tmp);

  r.set_vector(tmp, 1);

  int nhits = Substructure_Ring_Base::_environment_matches(target, tmp);

  if (0 == nhits)
    return ! _match_as_match_or_rejection;

  return _match_as_match_or_rejection;
}
