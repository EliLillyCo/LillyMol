/*
  Look up fragments in a dicer database
*/

#include <stdlib.h>
#include <memory>
#include <limits>
using namespace std;

#include "db_cxx.h"

#include "cmdline.h"
#include "iwstring_data_source.h"
#include "accumulator.h"
#include "misc.h"

#include "molecule.h"
#include "molecule_to_query.h"
#include "target.h"
#include "aromatic.h"

const char * prog_name = NULL;

static int verbose = 0;

static int molecules_read = 0;

static int reduce_to_largest_fragment = 0;

static int lower_atom_count_cutoff = 0;

static char output_separator = ' ';

static int molecules_with_no_fragments = 0;

static int molecules_with_zero_fragments_found = 0;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Builds a dicer fragment database - takes output from dicer\n";
  cerr << "  -d <dbname>   database name\n";
  cerr << "  -c <natoms>   ignore fragments with fewer than <natoms> atoms\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";

  exit(rc);
}

static void
preprocess (Molecule & m)
{
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  return;
}

static int
do_no_fragments_found (Molecule & parent,
                       IWString_and_File_Descriptor & output)
{
  molecules_with_zero_fragments_found++;

  output << output_separator << '1';     // fraction zero

  output << output_separator << '0' << output_separator << '0' << output_separator << '0';

  output << output_separator << '0';

  output << output_separator << '0';

  return 1;
}

static int
do_output (Molecule & parent,
           const int * t,
           IWString_and_File_Descriptor & output)
{
  const int matoms = parent.natoms();

  const int nedges = parent.nedges();

  const int * max_size_hit       = t;
  const int * max_size_hit_count = t + matoms;
  const int * max_count          = t + 2 * matoms;
  const int * max_count_size     = t + 3 * matoms;
  const int * total_hits         = t + 4 * matoms;
  const int * hits_to_bonds      = t + 5 * matoms;

  int zero_hits = 0;

  int largest_fragment = 0;
  int largest_fragment_count = 0;

  int smallest_fragment = matoms;
  int smallest_fragment_count = 0;

  int highest_count = 0;
  int highest_count_size = 0;

  Accumulator_Int<int> tot;

  cerr << "Atom contains " << matoms << " atoms\n";
  for (int i = 0; i < matoms; ++i)
  {
    cerr << total_hits[i] << " hits to atom " << i << endl;

    if (0 == total_hits[i])
    {
      zero_hits++;
      continue;
    }

    if (max_size_hit[i] < smallest_fragment)
    {
      smallest_fragment = max_size_hit[i];
      smallest_fragment_count = max_size_hit_count[i];
    }

    if (max_size_hit[i] > largest_fragment)
    {
      largest_fragment = max_size_hit[i];
      largest_fragment_count = max_size_hit_count[i];
    }

    if (max_count[i] > highest_count)
    {
      highest_count = max_count[i];
      highest_count_size = max_count_size[i];
    }

    tot.extra(total_hits[i]);
  }

  int zero_hits_bonds = 0;

  for (int i = 0; i < nedges; ++i)
  {
    if (0 == hits_to_bonds[i])
      zero_hits_bonds++;
  }

  output << output_separator;

  if (0 == zero_hits)     // fraction zero
    output << '0';
  else
    output << static_cast<float>(zero_hits) / static_cast<float>(matoms);

  cerr << "tot contains " << tot.n() << " values\n";
  output << output_separator << tot.minval() << output_separator << tot.maxval() << output_separator << static_cast<float>(tot.average());

  output << output_separator << static_cast<float>(largest_fragment) / static_cast<float>(matoms);

  output << output_separator << static_cast<float>(zero_hits_bonds) / static_cast<float>(nedges);

  return 1;
}

static int
append_count (const const_IWSubstring & buffer,
              resizable_array<int> & fragment_counts)
{
  int i = buffer.rindex(' ');

  if (i < 0)   // will never happen
    return 0;

  const_IWSubstring s;
  buffer.from_to(i + 1, buffer.length() - 1, s);

  int c;

  if (! s.numeric_value(c) || c < 1)
  {
    cerr << "Invalid count '" << buffer << "'\n";
    return 0;
  }

  fragment_counts.add(c);

  return 1;
}

#ifdef NOTUSEDQ
static int
dicer_bdb_lookup_record (const const_IWSubstring & buffer,
                     const IWString & parent,
                     Db & database)
{
  IWString usmi;
  int i = 0;

  buffer.nextword(usmi, i);

  Dbt dbkey;

  dbkey.set_data(const_cast<char *> (usmi.rawchars ()));    // loss of const OK
  dbkey.set_size(usmi.length ());

  Dbt already_in_database;

  if (0 != database.get(NULL, &dbkey, &already_in_database, 0))  // not already in database
  {
#ifdef DEBUG_dicer_bdb_lookup_RECORD
    cerr << "Storing parent is '" << parent << "'\n";
#endif
    new_fragments_stored++;

    return do_store (database, dbkey, parent);
  }

// An entry with this key already in the database.

  IWString tmp;
  tmp.set (reinterpret_cast<char *>(already_in_database.get_data()), static_cast<int>(already_in_database.get_size()));

  if (3 != tmp.nwords())
  {
    cerr << "dicer_bdb_lookup_record:unexpected database contents '" << tmp << "' should contain 3 tokens\n";
    return 0;
  }

  int r = tmp.rindex(' ');

  const_IWSubstring zdigits;

  tmp.from_to(r + 1, tmp.length() - 1, zdigits);

#ifdef DEBUG_dicer_bdb_lookup_RECORD
  cerr << "From '" << tmp << "' extract string count '" << zdigits << "'\n";
#endif

  int c;
  if (! zdigits.numeric_value(c) || c < 1)
  {
    cerr << "dicer_bdb_lookup_record:the count stored in the database must be a whole +ve number '" << tmp << "' invalid\n";
    return 0;
  }

#ifdef DEBUG_dicer_bdb_lookup_RECORD
  cerr << "Count determined to be " << c << endl;
#endif

  tmp.resize_keep_storage(r + 1);

  tmp << (++c);

#ifdef DEBUG_dicer_bdb_lookup_RECORD
  cerr << "Storing updated db contents '" << tmp << "'\n";
#endif

  fragments_appended_to_existing_db_entries++;

  return do_store(database, dbkey, tmp);
}
#endif

static int
update_coverage (Molecule & parent,
                 Molecule_to_Match & target,
                 Molecule & f,
                 const int c,
                 int * t)
{
  Molecule_to_Query_Specifications mqs;

  mqs.set_atoms_conserve_ring_membership(1);
  mqs.set_substituents_only_at_isotopic_atoms(1);

  Substructure_Query q;

  if (! q.create_from_molecule(f, mqs))
  {
    cerr << "Cannot convert molecule to substructure query '" << f.smiles() << "'\n";
    return 0;
  }

  q.set_find_one_embedding_per_atom(1);
  q.set_find_unique_embeddings_only(1);

  Substructure_Results sresults;

  int nhits = q.substructure_search(target, sresults);

  if (0 == nhits)
  {
    cerr << "Very strange, no hits in parent " << parent.smiles() << " from query generated by " << f.smiles() << endl;
    return 0;
  }

  const int matoms = parent.natoms();

  const int nedges = parent.nedges();

  int * max_size_hit       = t;
  int * max_size_hit_count = t + matoms;
  int * max_count          = t + 2 * matoms;
  int * max_count_size     = t + 3 * matoms;
  int * total_hits         = t + 4 * matoms;
  int * hits_to_bonds      = t + 5 * matoms;

  cerr << "Processing fragment found " << c << " times in database\n";
  for (int i = 0; i < nhits; ++i)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    const int esize = e->number_elements();

    for (int j = 0; j < esize; ++j)
    {
      const auto k = e->item(j);
      if (esize > max_size_hit[k])
      {
        max_size_hit[k] = esize;
        max_size_hit_count[k] = c;
      }

      if (c > max_count[k])
      {
        max_count[k] = c;
        max_count_size[k] = esize;
      }
    }

    e->each(total_hits, [c] (int & j) { j += c;});

    for (int j = 0; j < nedges; ++j)
    {
      const Bond * b = parent.bondi(i);
      if (e->contains_atoms (b->a1(), b->a2()))
      {
        hits_to_bonds[j] += c;
      }
    }
  }

  return 1;
}

static int
update_coverage (Molecule & parent,
                 Molecule_to_Match & target,
                 const IWString & fsmiles,
                 const int c,
                 int * t)
{
  if (count_atoms_in_smiles(fsmiles) < lower_atom_count_cutoff)
    return 1;

  Molecule f;

  if (! f.build_from_smiles(fsmiles))
  {
    cerr << "Cannot interpret fragment smiles " << fsmiles << endl;
    return 0;
  }

  if (f.natoms() < lower_atom_count_cutoff) // redundant
    return 1;

    cerr << "LIne " << __LINE__ << " fragment smiles " << fsmiles << endl;
  return update_coverage (parent, target, f, c, t);
}

static int
extract_count (IWString & buffer,
               int & c)
{
  int i = buffer.rindex(' ');

  if (i < 0)
    return 0;

  const_IWSubstring s;

  buffer.from_to(i + 1, buffer.length() - 1, s);

  if (! s.numeric_value(c) || c < 1)
  {
    cerr << "Invalid numeric '" << buffer << "'\n";
    return 0;
  }

  buffer.iwtruncate(i - 1);

  return 1;
}

static int
update_coverage (Molecule & parent,
                 Molecule_to_Match & target,
                 const IWString & fsmiles,
                 const Dbt & indb,
                 int * t)
{
  IWString f;

  f.set (reinterpret_cast<char *>(indb.get_data()), static_cast<int>(indb.get_size()));

  int c;

  if (! extract_count(f, c))
    return 0;

  cerr << "Fragment found " << c << " times in db\n";
  return update_coverage (parent, target, fsmiles, c, t);
}

static int
build_molecule (const const_IWSubstring & buffer,
                Molecule & m)
{
  m.resize(0);

  int i = 0;
  const_IWSubstring token;

  if (! buffer.nextword(token, i))
  {
    cerr << "Cannot extract smiles from '" << buffer << "'\n";
    return 0;
  }

  if (! m.build_from_smiles(token))
  {
    cerr << "Invalid smiles '" << buffer << "'\n";
    return 0;
  }

  if (! buffer.nextword(token, i))
  {
    cerr << "No hame available '" << buffer << "'\n";
    return 0;
  }

  m.set_name(token);

  return 1;
}

static int
dicer_bdb_lookup (Molecule & parent,
                  resizable_array_p<IWString> & fragments,
                  Db ** database,
                  const int ndb,
                  IWString_and_File_Descriptor & output)
{
  const int matoms = parent.natoms();

  int * t = new_int(5 * matoms + parent.nedges());  unique_ptr<int[]> free_t(t);

  Molecule_to_Match target(&parent);

  int dbfound = 0;

  const int nf = fragments.number_elements();

  for (int i = 0; i < nf; ++i)
  {
    Dbt dbkey;

    dbkey.set_data((void *)fragments[i]->rawchars());
    dbkey.set_size(fragments[i]->length());

    for (int j = 0; j < ndb; ++j)
    {
//    cerr << "Looking up " << *fragments[i] << endl;

      Dbt indb;

      if (0 != database[j]->get(NULL, &dbkey, &indb, 0))
        continue;

      cerr << "Found " << *fragments[i] << endl;

      update_coverage(parent, target, *fragments[i], indb, t);

      dbfound++;
    }
  }

  cerr << parent.name() << " found " << dbfound << " of " << nf << " fragments\n";

  output << parent.smiles() << output_separator << parent.name();

  if (0 == dbfound)
    do_no_fragments_found(parent, output);
  else
    do_output (parent, t, output);

  output << '\n';

  output.write_if_buffer_holds_more_than(4196);

  return 1;
}

static int
dicer_bdb_lookup (iwstring_data_source & input,
                  Db ** database,
                  const int ndb,
                  IWString_and_File_Descriptor & output)
{
  const_IWSubstring buffer;

  Molecule parent;
  resizable_array_p<IWString> fragments;

  while (input.next_record (buffer))
  {
    if (buffer.contains (" B="))
    {
      if (fragments.size () > 0)
      {
        if (! dicer_bdb_lookup(parent, fragments, database, ndb, output))
          return 0;
      }
      else
        molecules_with_no_fragments++;

      if (! build_molecule (buffer, parent))
      {
        cerr << "dicer_bdb_lookup:invalid parent molecule record '" << buffer << "'\n";
        return 0;
      }

      fragments.resize_keep_storage(0);

      molecules_read++;
    }
    else
    {
      IWString * f = new IWString;
      int i = 0;
      buffer.nextword(*f, i);       // first token on the line is the unique smiles
      fragments.add(f);
    }
  }

  if (fragments.size() > 0)
    return dicer_bdb_lookup(parent, fragments, database, ndb, output);

  return 1;
}

static int
dicer_bdb_lookup (const char * fname,
                  Db ** database,
                  const int ndb,
                  IWString_and_File_Descriptor & output)
{
  assert (NULL != fname);

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  return dicer_bdb_lookup(input, database, ndb, output);
}

static int
dicer_bdb_lookup (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:E:ld:c:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

//set_substituents_only_at_isotopic_atoms(1);    moved to mqs object

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  }
  else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('l'))
  {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', lower_atom_count_cutoff) || lower_atom_count_cutoff < 1)
    {
      cerr << "The fragment lower atom count limit (-c) must be a whole +ve number\n";
      usage(1);
    }

    if (verbose)
      cerr << "Will discard fragments having fewer than " << lower_atom_count_cutoff << " atoms\n";
  }

  if (! cl.option_present('d'))
  {
    cerr << "Must specify dicer database via the -d option\n";
    usage(1);
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  int ndb = cl.option_count('d');

  Db ** database = new Db *[ndb]; unique_ptr<Db*[]> free_database(database);

  if (cl.option_present('d'))
  {
    IWString dbname;
    for (int i = 0; cl.value('d', dbname, i); ++i)
    {
      int flags;
      DBTYPE dbtype;
      int mode;
    
      if (dash_s(dbname.null_terminated_chars()))
      {
        dbtype = DB_UNKNOWN;
        flags = 0;
        mode = 0;
      }
      else
      {
        dbtype = DB_BTREE;
        flags = DB_CREATE;
        mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;
      }

      database[i] = new Db(NULL, DB_CXX_NO_EXCEPTIONS);

//    cerr << "Just about to open database " << i << " name '" << dbname << "'\n";
      int rc = database[i]->open (NULL, dbname.null_terminated_chars(), NULL, dbtype, flags, mode);

      if (0 != rc)
      {
        cerr << "Cannot open database '" << dbname << "'\n";
        database[i]->err(rc, "");
        return 2;
      }

      if (verbose)
        cerr << "Smiles will be read from database '" << dbname << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = 0; i < cl.number_elements(); i++)
  {
    if (! dicer_bdb_lookup(cl[i], database, ndb, output))
    {
      rc = i + 1;
      break;
    }
  }

  output.flush();

  for (int i = 0; i < ndb; ++i)
  {
    database[i]->close(0);
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules\n";
    cerr << molecules_with_no_fragments << " molecules came in with no fragments\n";
    cerr << molecules_with_zero_fragments_found << " molecules found zero fragments in db\n";
  }

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = dicer_bdb_lookup (argc, argv);

  return rc;
}
