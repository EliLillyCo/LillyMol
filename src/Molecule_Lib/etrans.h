#ifndef MOLECULE_LIB_ETRANS_H_
#define MOLECULE_LIB_ETRANS_H_

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

#ifdef BUILD_BAZEL
#include "Molecule_Lib/etrans.pb.h"
#else
#include "etrans.pb.h"
#endif

class Element;
class Molecule;
class Molecule_to_Match;

#include "ematch.h"

class Element_Transformation
{
  private:
    int   _transform_every_atom_type;
    Element_Matcher _from;
    const Element * _to;
    int _isotope;

    // When changing from I to Cl, we can get a bad valence if we are
    // translating a 3 connected Iodone.
    int _reverse_bad_valence;

    int _molecules_processed;
    int _molecules_changed;
    int _atoms_changed;

//  private functions

    void _default_values ();
    int MaybeChangeElement(Molecule& m, atom_number_t zatom);

  public:
    Element_Transformation ();

    int ok () const;
    int debug_print (std::ostream &) const;

//  int build(const char *);
    int build(const IWString& value);

    int Build(const element_transformation::ElementTransformation& proto);

    int process(Molecule &);

    int process(Molecule_to_Match &);
};

class Element_Transformations : public resizable_array_p<Element_Transformation>
{
  private:
  public:

    int ok() const;
    int debug_print(std::ostream &) const;

    int active() const { return _number_elements;}

    int construct_from_command_line(Command_Line &, int verbose= 0, char = 't');

    // Add a transformation directive 'Br=Cl' for example.
    int Add(const IWString& token);

    int Build(const element_transformation::ElementTransformations& proto);

    int process(Molecule *);

    int process(Molecule &);

    int process(Molecule_to_Match &);
};

extern int display_standard_etrans_options(std::ostream &, char = 't');

extern int process_element_transformations(Command_Line &,
                                           Element_Transformations &,
                                           int = 0,
                                           char = 't');
#endif  // MOLECULE_LIB_ETRANS_H_
