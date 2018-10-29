#include <stdlib.h>

#include "path.h"

#include "ring_ext_rep.h"


Ring_Extraction_Replacement_Conditions::Ring_Extraction_Replacement_Conditions()
{
  _ring_size_needed = 0;

  _ring_aromaticity_needed = 0;

  _isotope_for_ring_fusion = 0;

  _ring_fusion_element = NULL;

  _isotope_for_substitution_points = 0;

  _include_substituents = 0;

  _remove_chirality = 0;

  _only_process_fused_rings = 0;

  _only_process_unfused_rings = 0;

  _fused_neighbours_allowed = 0;

  return;
}

Ring_Extraction_Replacement_Conditions::~ Ring_Extraction_Replacement_Conditions ()
{
  return;
}

int
Ring_Extraction_Replacement_Conditions::initialise (Command_Line & cl,
                                                int verbose)
{
  if (cl.option_present('r'))
  {
    if (! cl.value('r', _ring_size_needed) || _ring_size_needed < 3)
    {
      cerr << "The ring size value (-r) must be a whole +ve number > 2\n";
      return 0;
    }

    if (verbose)
      cerr << "Will only extract rings of size " << _ring_size_needed << endl;
  }

  if (cl.option_present('f'))
  {
    const_IWSubstring f = cl.string_value('f');

    if (f.numeric_value(_isotope_for_ring_fusion))
    {
      if (_isotope_for_ring_fusion <= 0)
      {
        cerr << "Invalid isotope for ring fusion '" << f << "'\n";
        return 3;
      }

      if (verbose)
        cerr << "Will label ring fusion points with isotope '" << _isotope_for_ring_fusion << endl;
    }
    else
    {
      _ring_fusion_element = get_element_from_symbol_no_case_conversion(f);
      if (NULL == _ring_fusion_element)
        _ring_fusion_element = create_element_with_symbol(f);

      if (NULL == _ring_fusion_element)
      {
        cerr << "INvalid ring fusion element '" << f << "'\n";
        return 4;
      }

      if (verbose)
        cerr << "Ring fusions designated with '" << _ring_fusion_element->symbol() << "'\n";
    }

  }

  if (cl.option_present('F'))
  { 
    int i = 0;
    const_IWSubstring f;
    while (cl.value('F', f, i++))
    {
      if ("only" == f || f.starts_with("fuse"))
      {
        _only_process_fused_rings = 1;
        if (verbose)
          cerr << "Only fused rings will be extracted\n";
        continue;
      }
      else if ("isol" == f)
      {
        _only_process_unfused_rings = 1;
        if (verbose)
          cerr << "Only isolated rings will be extracted\n";
        continue;
      }

      if (! f.numeric_value(_fused_neighbours_allowed) || _fused_neighbours_allowed < 1)
      {
        cerr << "Invalid degree of fusion specifier '" << f << "'\n";
        return 0;
      }

      if (verbose)
        cerr << "Will only process fused rings with exactly " << _fused_neighbours_allowed << " fused neighbours\n";
    }
    
    if (_only_process_fused_rings && _only_process_unfused_rings)
    {
      cerr << "Cannot specify both '-f isol' and '-f only'\n";
      return 0;
    }
  }

  if (cl.option_present('c'))
  {
    _remove_chirality = 1;
    if (verbose)
      cerr << "All chiral information removed\n";
  }

  if (cl.option_present ('a'))
  {
    const_IWSubstring a = cl.string_value('a');
    if ("ar" == a)
    {
      _ring_aromaticity_needed = 2;
      if (verbose)
        cerr << "Will only fetch for aromatic rings\n";
    }
    else if ("al" == a)
    {
      _ring_aromaticity_needed = 1;
      if (verbose)
        cerr << "Will only fetch for aliphatic rings\n";
    }
    else
    {
      cerr << "Unrecognised -a qualifier '" << a << "'\n";
      return 4;
    }
  }

  if (cl.option_present('b'))
  {
    if (! cl.value('b', _isotope_for_substitution_points) || _isotope_for_substitution_points <= 0)
    {
      cerr << "The isotope for subsitution points (-b) must be a whole +ve number\n";
      return 0;
    }

    if (verbose)
      cerr << "Substitution points labelled with isotope " << _isotope_for_substitution_points << endl;
  }

  return 1;
}


int
display_standard_ring_ext_rep_options (std::ostream & os)
{
  os << "  -r <size>     ring size to fetch\n";
  os << "  -a ar         fetch only aromatic rings\n";
  os << "  -a al         fetch only aliphatic rings\n";
  os << "  -f <el>       use element <el> to designate ring fusion points\n";
  os << "  -f <iso>      use isotope <iso> to designate ring fusion points\n";
  os << "  -F fuse       only extract fused rings\n";
  os << "  -F isol       only extract unfused rings\n";
  os << "  -F <num>      only extract fused rings that have <num> fused neighbours\n";
  os << "  -b <iso>      isotopically label all points of substitution\n";
  os << "  -c            remove all chirality\n";

  return 1;
}


void
Ring_Extraction_Replacement_Conditions::_add_substituents (Molecule & m, 
                  atom_number_t ring_atom,
                  atom_number_t outside_ring,
                  int * include_atom) const
{
  const Atom * a = m.atomi(outside_ring);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++)
  {
    const Bond * b = a->item (i);

    atom_number_t j = b->other(outside_ring);

    if (ring_atom == j)
      continue;

    if (b->is_double_bond())
      include_atom[j] = 1;
    else if (6 != m.atomic_number(j) && 1 == m.ncon(j))
      include_atom[j] = 1;
  }

  return;
}

int
Ring_Extraction_Replacement_Conditions::identify_atoms_associated_with_ring (Molecule & m,
                                              const Ring & r,
                                              const int * in_same_ring,
                                              int * include_atom) const
{
  int matoms = m.natoms();

  int n = r.number_elements ();

  for (int i = 0; i < n; i++)
  {
    atom_number_t j = r[i];

    include_atom[j] = 1;

//  cerr << "Added atom " << j << " to the system\n";

    const Atom * aj = m.atomi(j);

    int jcon = aj->ncon ();

    for (int k = 0; k < jcon; k++)
    {
      const Bond * b = aj->item (k);

      atom_number_t l = b->other(j);

      if (r.contains(l))
        continue;

      if (in_same_ring[j * matoms + l])   // part of a fusion
      {
        if (NULL != _ring_fusion_element)
        {
          m.set_element (l, _ring_fusion_element);
          include_atom[l] = 1;
        }
        else if (_isotope_for_ring_fusion)
          m.set_isotope (j, _isotope_for_ring_fusion);
      }
      else if (b->is_double_bond())   // double bond fused to ring is always part of ring
        include_atom[l] = 1;
      else   // something outside the ring
      {
        if (_isotope_for_substitution_points)
          m.set_isotope(j, _isotope_for_substitution_points);
        if (_include_substituents)
          _add_substituents (m, j, l, include_atom);
      }
    }
  }

  return 1;
}

int
Ring_Extraction_Replacement_Conditions::identify_atoms_associated_with_ring_system (Molecule & m,
                                              int * include_atom) const
{
  int matoms = m.natoms();

  Set_of_Atoms to_be_added;

  int join_to_ring_atom = 0;

  for (int i = 0; i < matoms; i++)
  {
//  cerr << "Atom " << i << " type " << m.smarts_equivalent_for_atom(i) << " inc? " << include_atom[i] << endl;

    if (! include_atom[i])
      continue;

    const Atom * ai = m.atomi(i);

    const auto icon = ai->ncon ();

    if (2 == icon)     // nothing outside the ring system here
      continue;

    for (auto j = 0; j < icon; j++)
    {
      const Bond * b = ai->item (j);

      atom_number_t k = b->other(i);

      if (include_atom[k])
        continue;

      if (b->is_double_bond())   // double bond fused to ring is always part of ring
        to_be_added.add(k);
      else   // something outside the ring
      {
        join_to_ring_atom++;
        if (_isotope_for_substitution_points)
          m.set_isotope(i, _isotope_for_substitution_points);
//      if (_include_substituents)     does not work because it alters the include atom array
//        _add_substituents (m, i, k, include_atom);
      }
    }
  }

//cerr << "to_be_added " << to_be_added << endl;

  if (to_be_added.size() > 0)
    to_be_added.set_vector(include_atom, 1);

  return 1;
}

int 
Ring_Extraction_Replacement_Conditions::can_be_processed (Molecule & m,
                                                          const Ring & r) const
{
  if (_ring_size_needed > 0 &&  r.number_elements() != _ring_size_needed)
    return 0;

  if (r.strongly_fused_ring_neighbours())
    return 0;

  if (_only_process_unfused_rings && r.is_fused())
    return 0;

  if (_only_process_fused_rings && ! r.is_fused())
    return 0;

  if (0 == _fused_neighbours_allowed)
    ;
  else if (_fused_neighbours_allowed != r.fused_ring_neighbours())
    return 0;

  if (0 == _ring_aromaticity_needed)
    ;
  else if (1 == _ring_aromaticity_needed)
  {
    if (r.is_aromatic())
      return 0;
  }
  else if (2 == _ring_aromaticity_needed)
  {
    if (! r.is_aromatic())
      return 0;
  }

  return 1;
}

static int
count_included_connections (const Molecule & m,
                            atom_number_t zatom,
                            const int * include_atom)
{
  const Atom * a = m.atomi(zatom);

  int acon = a->ncon();

  int rc = 0;

  for (int i = 0; i < acon; i++)
  {
    atom_number_t j = a->other (zatom, i);

    if (include_atom[j])
      rc++;
  }

  return rc;
}

int 
Ring_Extraction_Replacement_Conditions::append_connectivity_smarts (Molecule & m,
                                                        atom_number_t zatom,
                                                        const int * include_atom,
                                                        int aromatic,
                                                        IWString & smarts) const
{
//cerr << "append_connectivity_smarts, aromatic? " << aromatic << endl;

  smarts << '[';
  if (! m.is_ring_atom(zatom))
    ;
  else if (aromatic)
    smarts << 'a';
  else
    smarts << 'A';

  int iso = m.isotope(zatom);

  if (0 == iso)
//  smarts << "R1D" << count_included_connections(m, zatom, include_atom) << ']';
    smarts << "R1D" << m.ncon(zatom) << ']';
  else if (_isotope_for_ring_fusion == iso)
    smarts << "R2D" << m.ncon(zatom) << "]";
  else if (_isotope_for_substitution_points == iso)
    smarts << "R1D" << m.ncon(zatom) << "]";

  return 1;
}

/*
  This variant is called when processing ring systems
*/

void
Ring_Extraction_Replacement_Conditions::append_connectivity_smarts (Molecule & m,
                                                                    atom_number_t zatom,
                                                                    int aromatic,
                                                                    IWString & smarts) const
{
  smarts << '[';

  if (! m.is_ring_atom(zatom))
    ;
  else if (aromatic)
    smarts << 'a';
  else if (2 == m.ncon(zatom))
    smarts << 'A';

  smarts << 'R' << m.nrings(zatom);

  const auto iso = m.isotope(zatom);

  if (0 == iso)
    smarts << 'D' << m.ncon(zatom);
  else if (_isotope_for_substitution_points == iso)
  {
    if (1 == m.formal_charge(zatom) && 7 == m.atomic_number(zatom))
      smarts << "D4";
    else
      smarts << "D>2";
  }

  smarts << ']';

  return;
}

int
initialise_in_same_ring_array (Molecule & m,
                               int * in_same_ring)
{
  int matoms = m.natoms ();

  int nr = m.nrings ();

  for (int i = 0; i < nr; i++)
  {
    const Ring * ri = m.ringi(i);

    int ring_size = ri->number_elements();
    for (int j = 0; j < ring_size; j++)
    {
      atom_number_t k = ri->item(j);
      for (int l = j + 1; l < ring_size; l++)
      {
        atom_number_t n = ri->item(l);
        in_same_ring[n * matoms + k] = 1;
        in_same_ring[k * matoms + n] = 1;
      }
    }
  }

  return 1;
}
