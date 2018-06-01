#ifndef MOLECULE_SMARTS_H
#define MOLECULE_SMARTS_H

  private:
    int _smarts (atom_number_t astart,
                   int * include_atom,
                   int flag,
                   IWString & s);

    void _compute_ncon_and_explicit_hydrogens (atom_number_t zatom,
                                                int & ncon,
                                                int & eh,
                                                const int * include_atom) const;
    void _append_isotope_and_atomic_symbol (atom_number_t zatom,
                                             IWString & smiles);
    int _append_smarts_equivalent_for_atom (atom_number_t zatom,
                                              int ncon,
                                              int rm,
                                              IWString & s) const;


#endif

/* arch-tag: 12f3775a-c3d9-4d19-a884-b63de1606b98

*/
