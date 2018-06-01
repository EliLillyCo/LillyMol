#ifndef IW_MOLECULAR_PROPERTIES_H
#define IW_MOLECULAR_PROPERTIES_H

#define NPROPERTIES 8

class Molecule;

class Molecular_Properties_Generator
{
  private:
    int _unsaturation_includes_aromatic;

// private functions

    template <typename T> int _molecular_properties_generation (Molecule & m, const int * ring_membership, T * properties) const;

  public:
    Molecular_Properties_Generator ();
    
    void set_unsaturation_includes_aromatic (int s) { _unsaturation_includes_aromatic = s;}

//  int operator () (Molecule &, unsigned char[]) const;
    template <typename T> int operator () (Molecule &, T[]) const;
};

#endif
