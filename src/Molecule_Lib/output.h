#ifndef IW_OUTPUT_H
#define IW_OUTPUT_H

/*
  This is an attempt at generalised molecule output.

  The intention is that it will take care of all details of writing
  a molecule.

  For each ofstream_and_type, the molecule will be written to that stream.

  For each _output_types, the molecule will be written to a file of that
  type.
*/

#include "iwstring.h"

#include "ostream_and_type.h"
#include "iwaray.h"

class Command_Line;
class Molecule;
class MDL_File_Supporting_Material;


class Molecule_Output_Object : public resizable_array_p<ofstream_and_type>
{
  private:
    int _single_file_per_molecule;
    int _molecules_per_file;

    resizable_array<int> _output_types;

    int _verbose;

    int    _molecules_written;

//  When doing a fixed number of molecules per file, we need a file name stem

    IWString _stem;

//  We also need a counter on the number of files created

    int _files_created;

//  And we can also use a given token from the name field for the file name step

    int _name_token_for_fname;

//  We can write to stdout if required

    int _use_stdout;

//  private functions

    void _default_values ();

    int _open_new_stems (const const_IWSubstring & nstem, int = 0);

    int _suffix_already_used (const char * suffix);

    int _add_output_type (int, const const_IWSubstring &, MDL_File_Supporting_Material &);

  public:
    Molecule_Output_Object ();
    ~Molecule_Output_Object ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int active () const;
    int is_open () const { return active();}

    int good () const;

    int set_verbose (int i) { return _verbose = i;}

    int set_molecules_per_file (int);
    int molecules_per_file () const { return _molecules_per_file;}

    int name_token_for_file_name () const { return _name_token_for_fname;}
    void set_name_token_for_file_name (int t) { _name_token_for_fname = t;}

    int molecules_written () const { return _molecules_written;}

    int new_stem (const const_IWSubstring &, int = 0);
    const IWString & stem () const {return _stem;}

    int would_use_name (const char *) const;
    int would_use_name (const const_IWSubstring &, const char *) const;

    int would_overwrite_input_files (const Command_Line &, const const_IWSubstring &) const;

    int determine_output_types (const Command_Line &, char = 'o');
    int add_output_type (int);
    int number_output_types () const { return _output_types.number_elements ();}

    int write (Molecule *);
    int write (Molecule &);

//  very dangerous method. I did this for tsubstructure which needed to write atom and bond lists

    std::ostream & stream_for_type (int) const;

    int do_close();
    int do_flush();

    int first_output_type () const { return _output_types[0];}
};

#endif
