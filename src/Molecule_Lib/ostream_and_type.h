#ifndef IW_STREAMTYPE_H
#define IW_STREAMTYPE_H

#include <fstream>

#include "iwstring.h"

class Molecule;

#include "iwaray.h"

/*
  This class consists of an ofstream which knows which kind of
  structure file to write.
*/

class ofstream_and_type : public std::ofstream
{
  private:
    int _output_type;
    IWString _fname;
    int _valid;
    int _molecules_written;
    int _verbose;

//  private functions

    int _default_values ();

  public:
    ofstream_and_type ();
    ofstream_and_type (int);
    ofstream_and_type (int, const char *);
    ofstream_and_type (int, IWString &);
    ~ofstream_and_type ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int valid () const { return _valid;}

    int open (const char *);
    int open (IWString &);

    int  set_type (int);
    void set_verbose (int verbose) {_verbose = verbose;}

    const IWString & fname () const { return _fname;}

    int molecules_written () const { return _molecules_written;}

    int write_molecule (Molecule *);
    int write_molecules (const resizable_array_p<Molecule> &);
};

#endif
