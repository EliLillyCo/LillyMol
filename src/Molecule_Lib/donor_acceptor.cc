#include <stdlib.h>

#include "iwbits.h"
#include "cmdline.h"

#include "rwsubstructure.h"
#include "donor_acceptor.h"
#include "qry_wstats.h"
#include "target.h"
#include "misc.h"

void
display_standard_donor_acceptor_assigner_options (std::ostream & os, char cflag)
{
  os << "  -" << cflag << " a=<query>   specify query(s) for acceptors\n";
  os << "  -" << cflag << " d=<query>   specify query(s) for donors\n";
  os << "  -" << cflag << " label       apply isotopic labels to the matched atoms\n";
  os << "  -" << cflag << " add         add isotopic label to any existing isotopic label\n";
  os << "  -" << cflag << " atype       set entries in the atom type set\n";
  os << "  -" << cflag << " nodetach    don't check for explicit Hydrogens\n";
  os << "  -" << cflag << " write=<file> write labelled molecules to <file>\n";
  os << "  -" << cflag << " help        this message\n";

  return;
}


Donor_Acceptor_Assigner::Donor_Acceptor_Assigner ()
{
  _apply_isotopic_labels = 0;

  _apply_atom_type_labels = 0;

  _add_to_existing_isotopic_label = 0;

  return;
}

int
Donor_Acceptor_Assigner::ok () const
{
  return 1;
}

int
Donor_Acceptor_Assigner::construct_from_command_line (Command_Line & cl,
                                   char cflag,
                                   int verbose)
{
  if (! cl.option_present(cflag))
  {
    return 0;
  }

  IWString fname_for_labelled_molecules;

  int i = 0;
  const_IWSubstring c;
  while (cl.value(cflag, c, i++))
  {
    if (_temp_detach_hydrogens.recognised_directive(c))
      continue;

    if (c.starts_with("a="))
    {
      if (! _fetch_queries(c, _acceptor_queries))
      {
        cerr << "Cannot process acceptor queries '" << c << "'\n";
        return 0;
      }
    }
    else if (c.starts_with("d="))
    {
      if (! _fetch_queries(c, _donor_queries))
      {
        cerr << "Cannot process donor queries '" << c << "'\n";
        return 0;
      }
    }
    else if ("label" == c)
    {
      _apply_isotopic_labels = 1;
      if (verbose)
        cerr << "Isotopic labels will be applied\n";
    }
    else if ("add" == c)
    {
      _apply_isotopic_labels = 1;
      _add_to_existing_isotopic_label = 1;
      if (verbose)
        cerr << "Will append to any existing isotopic label\n";
    }
    else if ("atype" == c)
    {
      _apply_atom_type_labels = 1;
      if (verbose)
        cerr << "Will update the atom type of each atom\n";
    }
    else if (c.starts_with("write=") || c.starts_with("stream="))
    {
      fname_for_labelled_molecules = c.after('=');
    }
    else if ("help" == c)
    {
      display_standard_donor_acceptor_assigner_options(cerr, cflag);
      exit(0);
    }
    else
    {
      cerr << "Donor_Acceptor_Assigner::construct_from_command_line: unrecognised qualifier '" << c << "'\n";
      display_standard_donor_acceptor_assigner_options(cerr, cflag);
      exit(1);
    }
  }

  if (0 == fname_for_labelled_molecules.length())
    return 1;

  if (0 == _apply_isotopic_labels)
  {
    cerr << "Donor_Acceptor_Assigner::construct_from_command_line: label implied\n";
    _apply_isotopic_labels = 1;
  }

  if (! _open_stream_for_labelled_molecules(fname_for_labelled_molecules))
  {
    cerr << "Donor_Acceptor_Assigner::construct_from_command_line: cannot open '" << fname_for_labelled_molecules << "'\n";
    return 0;
  }

  if (verbose)
    cerr << "Isotopically labelled molecules written to '" << fname_for_labelled_molecules << "'\n";

  return 1;
}

int
Donor_Acceptor_Assigner::_assign_acceptors (Molecule_to_Match & target, int * isotope)
{
  int na = _acceptor_queries.number_elements();

  int rc = 0;
  for (int i = 0; i < na; i++)
  {
    Substructure_Results sresults;

    int nhits = _acceptor_queries[i]->substructure_search(target, sresults);
    if (0 == nhits)
      continue;

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * embedding = sresults.embedding(j);
      if (1 != embedding->number_elements())
      {
        cerr << "Yipes, donor embedding contains " << embedding->number_elements() << " atoms\n";
        return 0;
      }

      atom_number_t k = embedding->item(0);

      isotope[k] = DEFAULT_ACCEPTOR_ISOTOPIC_LABEL;
    }

    rc += nhits;
  }

  return rc;
}

//#define DEBUG_ASSIGN_DONORS

int
Donor_Acceptor_Assigner::_assign_donors (Molecule_to_Match & target, int * isotope)
{
  int na = _donor_queries.number_elements();

#ifdef DEBUG_ASSIGN_DONORS
  cerr << "Processing " << na << " donor queries\n";
#endif

  int rc = 0;
  for (int i = 0; i < na; i++)
  {
    Substructure_Results sresults;
    int nhits = _donor_queries[i]->substructure_search(target, sresults);
    if (0 == nhits)
      continue;

#ifdef DEBUG_ASSIGN_DONORS
    cerr << nhits << " hits to query " << i << endl;
#endif

    for (int j = 0; j < nhits; j++)
    {
      const Set_of_Atoms * embedding = sresults.embedding(j);
      if (1 != embedding->number_elements())
      {
        cerr << "Yipes, donor embedding contains " << embedding->number_elements() << " atoms\n";
        return 0;
      }

      atom_number_t k = embedding->item(0);

#ifdef DEBUG_ASSIGN_DONORS
      cerr << "Hits atom " << k << " '" << target.molecule()->atomic_symbol(k) << "'\n";
#endif

      if (0 == isotope[k])
        isotope[k] = DEFAULT_DONOR_ISOTOPIC_LABEL;
      else if (1 == isotope[k])
      {
        isotope[k] = DEFAULT_DUAL_ISOTOPIC_LABEL;           // dual nature, donor and acceptor
//      dual++;
      }
    }

    rc += nhits;
  }

  return rc;
}

int
Donor_Acceptor_Assigner::_process (Molecule_to_Match & target,
                                  int * isotope)
{
  int rc = _assign_acceptors(target, isotope);
  rc    += _assign_donors(target, isotope);

  return rc;
}

int
Donor_Acceptor_Assigner::__process (Molecule & m, 
                                    int * isotope)
{
  Molecule_to_Match target(&m);

  int rc = _process(target, isotope);

//cerr << "Donor_Acceptor_Assigner::__process:finished " << m.smiles() << " _apply_isotopic_labels " << _apply_isotopic_labels << endl;
  if (_apply_isotopic_labels)
  {
    _do_apply_isotopic_labels(m, isotope);

    if (_stream_for_labelled_molecules.active())
      _stream_for_labelled_molecules.write(m);
  }

  if (_apply_atom_type_labels)
    _do_apply_atom_type_labels(m, isotope);

  return rc;
}

int
Donor_Acceptor_Assigner::_process (Molecule & m,
                                  int * isotope)
{
  _temp_detach_hydrogens.detach_atoms(m);

  int rc = __process(m, isotope);

  _temp_detach_hydrogens.reattach_atoms(m);

  return rc;
}

int
Donor_Acceptor_Assigner::process (Molecule & m, int * isotope)
{
  int i_own_the_vector;

  const int matoms = m.natoms();

  if (NULL == isotope)
  {
    isotope = new int[matoms];
    i_own_the_vector = 1;
  }
  else
  {
    i_own_the_vector = 0;
  }

  set_vector(isotope, matoms, 0);

  int rc = _process(m, isotope);

  if (i_own_the_vector)
    delete [] isotope;

  return rc;
}

int
Donor_Acceptor_Assigner::_do_apply_isotopic_labels (Molecule & m, const int * isotope) const
{
  int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; i++)
  {
    if (0 == isotope[i])
      continue;

    if (_add_to_existing_isotopic_label)
      m.set_isotope(i, isotope[i] + m.isotope(i));
    else
      m.set_isotope(i, isotope[i]);

    rc++;
  }

  return rc;
}

/*
  The default label for acceptors is bit 31, while donors get bit 30
*/

void
Donor_Acceptor_Assigner::_do_apply_atom_type_labels (Molecule & m,
                                                     const int * isotope) const
{
  int matoms = m.natoms();

  Atom_Types & atom_type = m.atom_types();

  for (int i = 0; i < matoms; i++)
  {
    if (0 == isotope[i])
      continue;

    if (DEFAULT_ACCEPTOR_ISOTOPIC_LABEL == isotope[i] || DEFAULT_DUAL_ISOTOPIC_LABEL == isotope[i])
      atom_type[i] |= one_bit_32[31];

    if (DEFAULT_DUAL_ISOTOPIC_LABEL == isotope[i] || DEFAULT_DONOR_ISOTOPIC_LABEL == isotope[i])
      atom_type[i] |= one_bit_32[30];
  }

  return;
}

int
Donor_Acceptor_Assigner::_fetch_queries (const_IWSubstring & c,
                                         resizable_array_p<Substructure_Hit_Statistics> & queries)
{
  assert ('=' == c[1]);

  c.remove_leading_chars(2);     // srtip off 'a=' or 'd='

  if (c.starts_with("F:"))
  {
    c.remove_leading_chars(2);
    if (! queries_from_file(c, queries, 1, 0))
    {
      cerr << "Cannot read query file 'F:" << c << "'\n";
      return 0;
    }
  }
  else
  {
    Substructure_Hit_Statistics * q = new Substructure_Hit_Statistics;
    if (! q->read(c))
    {
      cerr << "Donor_Acceptor_Assigner::_fetch_queries: cannot construct query from '" << c << "'\n";
      delete q;
      return 0;
    }

    queries.add(q);
  }

  return queries.number_elements();
}

int
Donor_Acceptor_Assigner::_open_stream_for_labelled_molecules (const IWString & stem)
{
  _stream_for_labelled_molecules.add_output_type(SMI);

  return _stream_for_labelled_molecules.new_stem(stem, 1);
}

void
Donor_Acceptor_Assigner::deactivate()
{
 _donor_queries.resize(0);

 _acceptor_queries.resize(0);

 _apply_isotopic_labels = 0;

 _apply_atom_type_labels = 0;

 return;
}

int
Donor_Acceptor_Assigner::build (const const_IWSubstring & s)
{
  const_IWSubstring c;
  for (auto i = 0; s.nextword(c, i); )
  {
    if (c.starts_with("a="))
    {
      if (! _fetch_queries(c, _acceptor_queries))
      {
        cerr << "Cannot process acceptor queries '" << c << "'\n";
        return 0;
      }
    }
    else if (c.starts_with("d="))
    {
      if (! _fetch_queries(c, _donor_queries))
      {
        cerr << "Cannot process donor queries '" << c << "'\n";
        return 0;
      }
    }
    else
    {
      cerr << "Donor_Acceptor_Assigner::build:unrecognised directive '" << c << "'\n";
      return 0;
    }
  }

  if (0 == _donor_queries.number_elements() || 0 == _acceptor_queries.number_elements())
  {
    cerr << "Donor_Acceptor_Assigner::build:incomplete specification\n";
    return 0;
  }

  _apply_isotopic_labels = 1;

  return 1;
}
