#ifndef MOLECULE_LIB_DONOR_ACCEPTOR_H_
#define MOLECULE_LIB_DONOR_ACCEPTOR_H_

#include <iostream>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/donor_acceptor.pb.h"
#include "output.h"
#include "temp_detach_atoms.h"

#include "Molecule_Lib/pharmacophore.pb.h"

class Molecule_to_Match;
class Substructure_Hit_Statistics;

#define DEFAULT_ACCEPTOR_ISOTOPIC_LABEL 1
#define DEFAULT_DONOR_ISOTOPIC_LABEL 3
#define DEFAULT_DUAL_ISOTOPIC_LABEL 2

class Donor_Acceptor_Assigner 
{
  private:
    resizable_array_p<Substructure_Hit_Statistics> _donor_queries;

    resizable_array_p<Substructure_Hit_Statistics> _acceptor_queries;

    int _apply_isotopic_labels;

    int _apply_atom_type_labels;

    int _add_to_existing_isotopic_label;

    Molecule_Output_Object _stream_for_labelled_molecules;

    Temp_Detach_Atoms _temp_detach_hydrogens;

//  private functions

    int  _assign_acceptors(Molecule_to_Match & target, int * isotope);
    int  _assign_donors(Molecule_to_Match & target, int * isotope);

    int  _do_apply_isotopic_labels(Molecule &, const int *) const;
    void _do_apply_atom_type_labels(Molecule &, const int *) const;

    // T will be either IWString or const_IWSubstring.
    template <typename T>
    int  _fetch_queries(T & c, resizable_array_p<Substructure_Hit_Statistics> &);

    int  _open_stream_for_labelled_molecules(const IWString &);

    int  _process(Molecule_to_Match &, int *);
    int  __process(Molecule &, int *);
    int  _process(Molecule &, int *);

  public:
    Donor_Acceptor_Assigner();

    int ok() const;
    int debug_print(std::ostream &) const;

    int active() const { return (_donor_queries.number_elements() || _acceptor_queries.number_elements()) ; }
    void deactivate();

    void set_apply_isotopic_labels(int s) { _apply_isotopic_labels = s;}

    int construct_from_command_line(Command_Line &, char, int = 0);
    int build(const const_IWSubstring &);
    int BuildFromProto(const Pharmacophore::DonorAcceptor& proto);

    // Construct from proto files.
    int BuildFromProto(const IWString& fname);
    int BuildFromProto(const BrunsDonorAcceptor::BrunsDonorAcceptor& proto, const IWString& dirname);

    int process(Molecule &, int * = nullptr);
};

void display_standard_donor_acceptor_assigner_options(std::ostream &, char);

#endif  // MOLECULE_LIB_DONOR_ACCEPTOR_H_
