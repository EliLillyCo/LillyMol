#ifndef MOLECULE_LIB_OSTREAM_AND_TYPE_H_
#define MOLECULE_LIB_OSTREAM_AND_TYPE_H_

#include <fstream>

#include "Foundational/iwstring/iwstring.h"
#include "iwmtypes.h"

class Molecule;

/*
  This class consists of an ofstream which knows which kind of
  structure file to write.
*/

class ofstream_and_type : public std::ofstream
{
  private:
    FileType _output_type;
    IWString _fname;
    int _valid;
    int _molecules_written;
    int _verbose;

//  private functions

    int _default_values ();

  public:
    ofstream_and_type ();
    ofstream_and_type (FileType);
    ofstream_and_type (FileType, const char *);
    ofstream_and_type (FileType, IWString &);
    ~ofstream_and_type ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int valid () const { return _valid;}

    int open (const char *);
    int open (IWString &);

//  int  set_type (int);
    int  set_type (FileType);
    void set_verbose (int verbose) {_verbose = verbose;}

    const IWString & fname () const { return _fname;}

    int molecules_written () const { return _molecules_written;}

    int write_molecule (Molecule *);
    int write_molecules (const resizable_array_p<Molecule> &);
};

#endif  // MOLECULE_LIB_OSTREAM_AND_TYPE_H_
