#ifndef IW_TRIPOS_H

/*
  We need a class to describe the residue information in a .mol2 file
*/

class Molecule;

class Tripos_Residue_Information
{
  private:
    int _number_residues;
    IWString * _residue_name;
    int * _atom_number_to_residue_number;

  public:
    Tripos_Residue_Information (int na);
    ~Tripos_Residue_Information ();
    
    int update_residue_information (int atom, int rnum, const IWString & rname);

    int remove_all_atoms_except (Molecule & m, const IWString & rname) const;
};

#endif
