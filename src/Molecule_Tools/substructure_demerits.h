#ifndef SUBSTRUCTURE_DEMERIT_H
#define SUBSTRUCTURE_DEMERIT_H

class Molecule;
class Demerit;
class Charge_Assigner;

namespace substructure_demerits
{
void set_verbose (int);

void set_keep_going_after_rejection (int);

/*
  Feb 2005. People want to be able to apply just the rejections from here
*/

void set_only_apply_rejection_rules ();

int hard_coded_queries_statistics (std::ostream &);
int initialise_hard_coded_queries_to_do (IWString &);

int hard_coded_queries (Molecule &, Demerit &);

/*
  We need a means of passing the charge assigner in substructure_demerits.cc back
  to the calling programme so it can be initialised from the command line. Awful!

  Or you can pass a string
*/

Charge_Assigner & charge_assigner ();
int initialise_charge_assigner(const char *);

void set_substructure_demerits_too_many_rings (int s);
void set_substructure_demerits_ring_size_too_large (int s);


void set_cx_chain_rejection_length(int s);
void set_all_numeric_demerit_values (int s);
};

#endif
