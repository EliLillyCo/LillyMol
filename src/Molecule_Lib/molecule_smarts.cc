#include <stdlib.h>

#define COMPILING_MOLECULE_SMARTS

#include "aromatic.h"
#include "molecule.h"
#include "path.h"
#include "smiles.h"

/*
  Are we making a smarts to describe the atom, or to allow subsequent substructure
  searches?
*/

static int make_smarts_embedding_value = 0;

void
set_make_smarts_embedding(int s)
{
  make_smarts_embedding_value = s;

  return;
}

int
make_smarts_embedding()
{
  return make_smarts_embedding_value;
}

static int include_unsaturation_in_embeddings_value = 0;

int
include_unsaturation_in_smarts_embeddings()
{
  return include_unsaturation_in_embeddings_value;
}

void
set_include_unsaturation_in_smarts_embeddings(int s)
{
  include_unsaturation_in_embeddings_value = s;
}

static int include_whole_molecule_ring_membership_in_subset_smarts_value = 0;

int
include_whole_molecule_ring_membership_in_subset_smarts()
{
  return include_whole_molecule_ring_membership_in_subset_smarts_value;
}

void
set_include_whole_molecule_ring_membership_in_subset_smarts(int s)
{
  include_whole_molecule_ring_membership_in_subset_smarts_value = s;
}

static int include_D_in_smarts_value = 1;

void 
set_include_D_in_smarts(int s)
{
  include_D_in_smarts_value = s;
}

int
include_D_in_smarts()
{
  return include_D_in_smarts_value;
}

IWString
Molecule::smarts_equivalent_for_atom(atom_number_t zatom)
{
  IWString s;

  append_smarts_equivalent_for_atom(zatom, s);

  return s;
}

IWString
Molecule::smarts_equivalent_for_atom(atom_number_t zatom) const
{
  IWString s;

  append_smarts_equivalent_for_atom(zatom, s);

  return s;
}

IWString
Molecule::const_smarts_equivalent_for_atom(atom_number_t zatom) const
{
  IWString s;

  append_smarts_equivalent_for_atom(zatom, s);

  return IWString(s);
}

#ifdef NOT_NEEDED_NOW_MAYBE_SOMETIME
static void
append_D_smarts(IWString & s,
                int ncon,
                int eh)
{
  if (eh)
  {
    ncon = ncon - eh;
    assert (ncon >= 0);
  }

  if (ncon > 0)
    s << "D>" << (ncon - 1);

  return;
}
#endif


static void
append_formal_charge(IWString & s,
                     formal_charge_t fc)
{
  if (0 == fc)
    ;
  else if (1 == fc)
    s += '+';
  else if (-1 == fc)
    s += '-';

  return;
}

int
Molecule::_append_smarts_equivalent_for_atom(atom_number_t zatom,
                                             int ncon,
                                             int rm,
                                             IWString & s) const
{
  const Atom * a = _things[zatom];

  s << '[';

  if (a->is_isotope())
    s << a->isotope();

  if (1 == a->atomic_number())
    s << "#1";
  else if (nullptr != _aromaticity && is_aromatic_atom(_aromaticity[zatom]))
    s << a->element()->aromatic_symbol();
  else if (a->element()->symbol().length() > 2)
  {
    s << "#{" << a->element()->symbol() << '}';
  }
  else
    s << a->atomic_symbol();

//cerr << "Atom " << zatom << " value " << make_smarts_embedding_value << endl;

  if (! include_D_in_smarts_value)
    ;
  else if (! make_smarts_embedding_value)
    s << 'D' << ncon;
  else if (ncon > 1)
    s << "D>" << (ncon - 1);

  append_formal_charge(s, a->formal_charge());

  if (0 == rm)
    ;
  else if (1 == rm)
    s << 'R';
  else
    s << 'R' << rm;

  if (! include_unsaturation_in_embeddings_value)
    ;
  else if (a->ncon() == a->nbonds())
    ;
  else if (nullptr != _aromaticity && is_aromatic_atom(_aromaticity[zatom]))
    ;
  else
    s << "G>0";

  s << ']';

  return 1;
}

int
Molecule::append_smarts_equivalent_for_atom(atom_number_t zatom,
                                            IWString & s) const
{
  assert (ok_atom_number(zatom));
  
  const Atom * a = _things[zatom];

  return _append_smarts_equivalent_for_atom(zatom, a->ncon(), 0, s);
}

void
Molecule::_compute_ncon_and_explicit_hydrogens(atom_number_t zatom,
                                               int & ncon,
                                               int & eh,
                                               const int * include_atom) const
{
  ncon = 0;
  eh = 0;

  const Atom * a = _things[zatom];

  for (int i = 0; i < a->ncon(); i++)
  {
    atom_number_t j = a->other(zatom, i);

    if (1 == _things[j]->atomic_number())    // hmmm, should this go after the check for inclusion?
      eh++;

    if (include_atom[j])
      ncon++;
  }

  return;
}

void
Molecule::_append_isotope_and_atomic_symbol(atom_number_t zatom,
                                            IWString & s)
{
  const Atom * a = _things[zatom];

  if (a->is_isotope())
    s << a->isotope();

  if (1 == a->atomic_number())
    s << "#1";
  else if (is_aromatic(zatom))
    s << a->element()->aromatic_symbol();
  else if (a->atomic_symbol().length() > 2)
  {
    s << "#{" << a->atomic_symbol() << '}';
  }
  else
    s << a->atomic_symbol();

  return;
}

int
Molecule::append_smarts_equivalent_for_atom(atom_number_t zatom,
                                            IWString & s,
                                            const int * include_atom)
{
  int ncon, eh;

  _compute_ncon_and_explicit_hydrogens(zatom, ncon, eh, include_atom);

  int nr = nrings();

  if (0 == nr || 0 == nrings(zatom))
    return _append_smarts_equivalent_for_atom(zatom, ncon - eh, 0, s);

// For this to be a ring atom, zatom must be part of a ring that is fully included in include_atom

  int ring_membership = 0;

  if (include_whole_molecule_ring_membership_in_subset_smarts_value)
    ring_membership = nrings(zatom);
  else
  {
    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = ringi(i);

      if (ri->all_members_non_zero_in_array(include_atom) && ri->contains(zatom))
        ring_membership++;
    }
  }

  return _append_smarts_equivalent_for_atom(zatom, ncon - eh, ring_membership, s);
}

/*
  Non const version can do rings and such
*/

int
Molecule::append_smarts_equivalent_for_atom(atom_number_t zatom,
                                            IWString & s)
{
  assert (ok_atom_number(zatom));
  
  s << '[';

  _append_isotope_and_atomic_symbol(zatom, s);

  const Atom * a = _things[zatom];

// Apr 2005. If making an embedding, including D>xxx directives just isn't useful

  if (! make_smarts_embedding_value)
  {
    if (include_D_in_smarts_value)
      s << 'D' << a->ncon();

    s << 'H' << hcount(zatom);

    if (a->ncon() > 0)
      s << 'v' << a->nbonds();
  }

  append_formal_charge(s, a->formal_charge());
  
  if (! make_smarts_embedding_value && is_ring_atom(zatom))
  {
    int nr = nrings();

    for (int i = 0; i < nr; i++)
    {
      const Ring * ri = ringi(i);

      if (ri->contains(zatom))
        s << ";r" << ri->number_elements();
    }
  }

  s << ']';

  return 1;
}

/*

  Seems like a good idea, flesh this out sometime

int
Molecule::append_smarts_equivalent_for_atom(atom_number_t zatom,
                                            IWString & s,
                                            const Smiles_Formation_Info & sfi)
{
  s << '[';

  _append_isotope_and_atomic_symbol(zatom, s);

  const Atom * a = _things[zatom];

  const int * include_atom = sfi.include_atom();

  int ncon, eh;

  if (nullptr == include_atom)
  {
    ncon = a->ncon();
    eh = explicit_hydrogens(zatom);
  }
  else
     _compute_ncon_and_explicit_hydrogens(zatom, ncon, eh, include_atom);

  append_D_smarts(s, ncon, eh);

  append_formal_charge(s, a->formal_charge());

  if (is_ring_atom(zatom))
    s << 'R';

   s << ']';

  return 1;
}*/

const IWString &
Molecule::smarts()
{
  return smarts(_smiles_information, nullptr);
}

const IWString &
Molecule::smarts(Smiles_Information & smi_info,
                 const int * include_atom)
{
  smi_info.set_smiles_is_smarts(1);

  smi_info.prepare_to_build_ordering(_number_elements);

  smi_info.set_create_smarts_embedding(1);

  (void) ring_membership();

  int ss = write_smiles_with_smarts_atoms();

  int sa = get_include_aromaticity_in_smiles();

  set_write_smiles_with_smarts_atoms(1);

  set_include_aromaticity_in_smiles(1);

  (void) smiles(smi_info, include_atom);

  set_write_smiles_with_smarts_atoms(ss);

  set_include_aromaticity_in_smiles(sa);

  smi_info.set_smiles_is_smarts(1);

  return smi_info.smiles();
}

const IWString &
Molecule::smarts_starting_with_atom(atom_number_t a,
                                    Smiles_Information & smi_info,
                                    const int * include_atom)
{
  (void) ring_membership();    // lots of problems without this - because the smiles subset building does not properly populate the fragment info

  int ss = write_smiles_with_smarts_atoms();

  int sa = get_include_aromaticity_in_smiles();

  set_write_smiles_with_smarts_atoms(1);

  smi_info.set_smiles_is_smarts(1);

  compute_aromaticity_if_needed();

  (void) smiles_starting_with_atom(a, smi_info, include_atom);

  set_write_smiles_with_smarts_atoms(ss);

  set_include_aromaticity_in_smiles(sa);

  smi_info.set_smiles_is_smarts(1);

  return smi_info.smiles();
}

int
Molecule::write_molecule_smarts(std::ostream & os)
{
  invalidate_smiles();

  int tmp = make_smarts_embedding_value;

  make_smarts_embedding_value = 1;

  os << smarts() << ' ' << _molecule_name << '\n';

  make_smarts_embedding_value = tmp;

  return os.good();
}
