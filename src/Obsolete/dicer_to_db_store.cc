/*
  Scans a descriptor file for similarity to a given vector
*/

#include <stdlib.h>
#include <sys/stat.h>

#include "gdbm.h"

#include "cmdline.h"
#include "iw_stl_hash_map.h"
#include "iwstring_data_source.h"

/*
  Bit of a kludge here. I don't want to include molecule.h in this programme
*/

extern int count_atoms_in_smiles(const const_IWSubstring & smiles);

const char * prog_name = NULL;

static int verbose = 0;

static int min_atoms_in_fragment = 0;
static int max_atoms_in_fragment = 0;

static int fragments_with_too_few_atoms = 0;
static int fragments_with_too_many_atoms = 0;

static int max_examples_to_store = 1;
static int fragments_with_already_enough_examples = 0;

static int molecules_read = 0;
static int fragments_read = 0;

static int new_fragments_encountered = 0;

static GDBM_FILE database = NULL;

static int input_is_sorted = 0;

static extending_resizable_array<int> examples_per_fragment;

/*
  We can speed things up with a hash
*/

static IW_STL_Hash_Map_int hash;

static int maximum_hash_size = 0;

static IWString * most_recently_used = NULL;

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Take dicer output and builds a gdbm database\n";
  cerr << " -d <dbname>    name of database to create\n";
  cerr << " -n <number>    how many ID's (examples) to store for each fragment\n";
  cerr << " -c <number>    minimum number of atoms in fragment\n";
  cerr << " -C <number>    maximum number of atoms in fragment\n";
  cerr << " -s             input is sorted by smiles - much faster processing\n";
  cerr << "                consider also sorting by gdbm_sort_file_by_gdbm_hash_function\n";
  cerr << " -x <size>      maximum hash size\n";
  cerr << " -v             verbose output\n";

  exit (rc);
}

static int
do_store(datum & dkey,
         const const_IWSubstring & smiles,
         const IWString & new_data_to_store)
{

  datum to_store;
  to_store.dptr = const_cast<char *>(new_data_to_store.rawchars());
  to_store.dsize = new_data_to_store.length();

  if (0 == gdbm_store (database, dkey, to_store, GDBM_REPLACE))
    return 1;

  cerr << "Yipes, cannot store key '" << smiles << "', " << gdbm_strerror (gdbm_errno) << endl;

  return 0;
}

static int
store_new_record(datum & dkey,
                 const const_IWSubstring & smiles,
                 const const_IWSubstring & id,
                 int nsamples)
{
//cerr << "Storing '" << smiles << "' id '" << id << "' nsamples " << nsamples << endl;

  IWString new_data_to_store;

  new_data_to_store << nsamples << ' ' << id;

  new_fragments_encountered++;

  return do_store(dkey, smiles, new_data_to_store);
}

static int
dicer_to_db(const const_IWSubstring & smiles,
            const const_IWSubstring & id,
            int nsamples = 1)
{
//cerr << "Storing '" << smiles << "', id '" << id << "', count " << nsamples << endl;

  datum dkey;
  dkey.dptr = const_cast<char *>(smiles.rawchars());
  dkey.dsize = smiles.length();

  datum fromdb = gdbm_fetch(database, dkey);

  if (NULL == fromdb.dptr)
    return store_new_record(dkey, smiles, id, nsamples);

  IWString existing_data;

  existing_data.set_and_assume_ownership(fromdb.dptr, fromdb.dsize);

  int nw = existing_data.nwords();

  if (nw < 2)
  {
    cerr << "Huh, existing database data must be at least 2 tokens\n";
    cerr << "'" << existing_data << "' cannot be right, key '" << smiles << "'\n";
    return 0;
  }

  const_IWSubstring string_count, examples;

  existing_data.split(string_count, ' ', examples);

  int c;
  if (! string_count.numeric_value(c) || c <= 0)
  {
    cerr << "First token of stored data must be count '" << existing_data << "' is invalid\n";
    cerr << "Key '" << smiles << "'\n";
    return 0;
  }

  c += nsamples;

  IWString new_data_to_store;
  new_data_to_store.resize(existing_data.size() + id.length() + 4);

  new_data_to_store << c << ' ' << examples;

  if (nw - 1 > max_examples_to_store)
    fragments_with_already_enough_examples++;
  else
    new_data_to_store << ' ' << id;

  return do_store(dkey, smiles, new_data_to_store);
}

/*
  We don't want to do anything chemically sensible, so we do a very rough
  text based atom counting
*/

static int
atoms_in_smiles_within_range(const const_IWSubstring & smiles)
{
  if (0 == min_atoms_in_fragment && 0 == max_atoms_in_fragment)
    return 1;

  int n = smiles.length();

  if (min_atoms_in_fragment > n)   // smiles with N chars cannot represent more than N atoms
    return 0;

  int rc = count_atoms_in_smiles(smiles);

  if (min_atoms_in_fragment > 0 && rc < min_atoms_in_fragment)
  {
    fragments_with_too_few_atoms++;
    return 0;
  }

  if (max_atoms_in_fragment > 0 && rc > max_atoms_in_fragment)
  {
    fragments_with_too_many_atoms++;
    return 0;
  }

  return 1;
}

static int
extract_first_two_tokens (const const_IWSubstring & buffer,
                          const_IWSubstring & smiles,
                          const_IWSubstring & id)
{
  int i = 0;
  if (! buffer.nextword(smiles, i))
    return 0;

  return buffer.nextword(id, i);
}

static int
dicer_to_db_record(const const_IWSubstring & buffer)
{
  if (buffer.nwords() < 2)
  {
    cerr << "Dicer output must have at least two tokens\n";
    return 0;
  }

  const_IWSubstring smiles, id;

  if (! extract_first_two_tokens(buffer, smiles, id))
    return 0;

  if (! atoms_in_smiles_within_range(smiles))
    return 1;

  return dicer_to_db(smiles, id);
}

static int
dicer_to_db_sorted (iwstring_data_source & input)
{
  IWString previous_smiles;
  IWString ids;
  int count = 0;

  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.contains(" B="))
    {
      molecules_read++;
      continue;
    }

    const_IWSubstring smiles, id;
    if (! extract_first_two_tokens(buffer, smiles, id))
    {
      cerr << "Cannot extract smiles and id '" << buffer << "'\n";
      return 0;
    }

    fragments_read++;

    if (smiles == previous_smiles)
    {
//    cerr << "Another example of " << smiles << "'\n";
      if (count < max_examples_to_store)
        ids.append_with_spacer(id);
      count++;
      continue;
    }

    if (0 == previous_smiles.length())
      ;
    else if (! atoms_in_smiles_within_range(previous_smiles))
      ;
    else
    {
//    cerr << count << " examples from '" << ids << "'\n";

      if (verbose)
        examples_per_fragment[count]++;

      if (! dicer_to_db (previous_smiles, ids, count))
        return 0;
    }

    previous_smiles = smiles;
    ids.resize_keep_storage(0);
    ids = id;
    count = 1;
  }

  if (count > 0)
  {
    if (verbose)
      examples_per_fragment[count]++;

    return dicer_to_db (previous_smiles, ids, count);
  }

  return 1;
}

static int
dicer_to_db(iwstring_data_source & input)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.contains(" B="))
    {
      molecules_read++;
      continue;
    }

    fragments_read++;

    if (! dicer_to_db_record(buffer))
    {
      cerr << "Fatal error reading '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
dicer_to_db (const char * fname)
{
  iwstring_data_source input (fname);

  if (! input.good ())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (input_is_sorted)
    return dicer_to_db_sorted (input);

  return dicer_to_db (input);
}


static int
dicer_to_db (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vd:n:c:C:x:s");

  if (cl.unrecognised_options_encountered ())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count('v');

  if (! cl.option_present('d'))
  {
    cerr << "Must specify name of database via the -d option\n";
    usage(4);
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', max_examples_to_store) || max_examples_to_store < 1)
    {
      cerr << "The max number of examples to store option(-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will store a max of " << max_examples_to_store << " examples of each fragment\n";
  }

  if (cl.option_present('c'))
  {
    if (! cl.value('c', min_atoms_in_fragment) || min_atoms_in_fragment < 1)
    {
      cerr << "The minumum atoms in a fragment value (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will only store fragments with " << min_atoms_in_fragment << " or more atoms\n";
  }

  if (cl.option_present('C'))
  {
    if (! cl.value('C', max_atoms_in_fragment) || max_atoms_in_fragment < min_atoms_in_fragment)
    {
      cerr << "The maxumum atoms in a fragment value (-C) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Will only store fragments with " << max_atoms_in_fragment << " or fewer atoms\n";
  }

  if (cl.option_present('x'))
  {
    if (! cl.value('x', maximum_hash_size) || maximum_hash_size < 1)
    {
      cerr << "The maximum hash size (-x) must be a whole +ve number\n";
      usage(23);
    }

    if (verbose)
      cerr << "Will maintain a cache of " << maximum_hash_size << " items\n";

    most_recently_used = new IWString[maximum_hash_size];
  }

  if (cl.option_present('s'))
  {
    input_is_sorted = 1;

    if (verbose)
      cerr << "Input will be assumed to be sorted\n";
  }

  if (0 == cl.number_elements ())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (cl.option_present('d'))
  {
    const char * dbname = cl.option_value('d');

    int flags = GDBM_WRCREAT;

    int mode = (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);

    database = gdbm_open (const_cast<char *> (dbname),
                          0,
                          flags,
                          mode,
                          NULL);

    if (NULL == database)
    {
      cerr << "Cannot open database '" << dbname << "' : " << gdbm_strerror (gdbm_errno) << endl;
      return 1;
    }

    if (verbose)
      cerr << "Opened database '" << dbname << "'\n";
  }

  int rc = 0;
  for (int i = 0; i < cl.number_elements (); i++)
  {
    if (! dicer_to_db (cl[i]))
    {
      rc = i + 1;
      break;
    }
  }

  if (verbose)
  {
    cerr << "Read " << molecules_read << " molecules, and " << fragments_read << " fragments\n";
    cerr << new_fragments_encountered << " new fragments encountered\n";
    cerr << fragments_with_already_enough_examples << " fragments already had " << max_examples_to_store << " or more examples stored\n";
    if (min_atoms_in_fragment)
      cerr << fragments_with_too_few_atoms << " fragments had fewer than " << fragments_with_too_few_atoms << " atoms\n";
    if (max_atoms_in_fragment)
      cerr << fragments_with_too_many_atoms << " fragments had more than " << max_atoms_in_fragment << " atoms\n";

    if (input_is_sorted)
    {
      for (int i = 0; i < examples_per_fragment.number_elements(); i++)
      {
        if (examples_per_fragment[i])
          cerr << examples_per_fragment[i] << " fragments with " << i << " examples\n";
      }
    }
  }

  cerr << "Closing database\n";
  if (NULL != database)
    gdbm_close(database);

  cerr << "Database closed\n";

  if (NULL != most_recently_used)
    delete [] most_recently_used;

  return rc;
}

int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = dicer_to_db (argc, argv);

  return rc;
}
