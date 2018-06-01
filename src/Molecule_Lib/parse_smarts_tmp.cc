#include <stdlib.h>

#ifdef PARSE_SMARTS_TMP_MALLOC_CHECK
#include "iwmalloc.h"
#endif

#include "substructure.h"
#include "parse_smarts_tmp.h"

Parse_Smarts_Tmp::Parse_Smarts_Tmp ()
{
  _last_query_atom_created = -1;

  return;
}

Parse_Smarts_Tmp::~Parse_Smarts_Tmp ()
{
#ifdef PARSE_SMARTS_TMP_MALLOC_CHECK
  check_malloc_magic ();
#endif

  return;
}

int
Parse_Smarts_Tmp::set_natoms (int n)
{
  assert (n > 0);

//cerr << "Parse_Smarts_Tmp:set_natoms: natoms " << n << endl;

  return 1;
}

#ifdef PARSE_SMARTS_TMP_MALLOC_CHECK

int
Parse_Smarts_Tmp::check_malloc_magic () const
{
  if (_root.number_elements ())
    iwmalloc_check_malloc_magic (_root.rawdata ());
  if (_no_matched_atoms_between.number_elements ())
    iwmalloc_check_malloc_magic (_no_matched_atoms_between.rawdata ());
  if (_link_atom.number_elements ())
    iwmalloc_check_malloc_magic (_link_atom.rawdata ());

  return 1;
}

#endif

// arch-tag: 3fa138c0-8e33-4af1-b8fa-4074ffbc7818
