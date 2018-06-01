#ifndef RWMOLECULE_H
#define RWMOLECULE_H

#define EXTRA_STRING_RECORD(ds, b, c) \
  if (! (ds).next_record ((b)))\
  {\
    cerr << (c) << " eof\n";\
    return 0;\
  }

template <typename T> int rwmolecule_error (const char *, T &);

extern int write_coordinates (std::ostream &, const Atom * a, int = 0);

extern int append_sybyl_atom_type (std::ostream & os, int atype, const IWString & asymbol);

extern void reset_rwmolecule_file_scope_variables ();

#endif
