#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwaray/iwaray.h"

#include "qry_wcharge.h"
#include "rwsubstructure.h"

template class resizable_array_p<Charge_Distribution>;
template class resizable_array_base<Charge_Distribution *>;

template class resizable_array_p<Query_and_Charge_Stats>;
template class resizable_array_base<Query_and_Charge_Stats *>;

template int process_queries (Command_Line &, resizable_array_p<Query_and_Charge_Stats> &, int, char);
template int process_files_of_queries (Command_Line &, resizable_array_p<Query_and_Charge_Stats> &, int, int, char);
template int process_cmdline_token (const char, const const_IWSubstring &, resizable_array_p<Query_and_Charge_Stats> &, int);

template int queries_from_file (iwstring_data_source &, resizable_array_p<Query_and_Charge_Stats> &, const IWString &, int);
template int queries_from_file (IWString const &, resizable_array_p<Query_and_Charge_Stats> &, int, int);
template int queries_from_file (const_IWSubstring const &, resizable_array_p<Query_and_Charge_Stats> &, int, int);

template int file_record_is_file (resizable_array_p<Query_and_Charge_Stats> & queries, const IWString & directory_path, IWString & buffer, int verbose);
template int file_record_is_smarts (resizable_array_p<Query_and_Charge_Stats> & queries, IWString & buffer, int verbose);

template int smarts_from_file (const const_IWSubstring &, resizable_array_p<Query_and_Charge_Stats> &, int);
template int smarts_from_file (iwstring_data_source &, resizable_array_p<Query_and_Charge_Stats> &, int);

template int queries_from_file_of_molecules<Query_and_Charge_Stats>(const_IWSubstring const &,
                                resizable_array_p<Query_and_Charge_Stats> &, int);
template int queries_from_file_of_molecules<Query_and_Charge_Stats>(const_IWSubstring const &,
                                Molecule_to_Query_Specifications & mqs,
                                resizable_array_p<Query_and_Charge_Stats> &, int);
template int queries_from_file_of_molecules<Query_and_Charge_Stats>(data_source_and_type<MDL_Molecule> &,
                                Molecule_to_Query_Specifications & mqs,
                                resizable_array_p<Query_and_Charge_Stats> &, int);
template int queries_from_file_of_molecules<Query_and_Charge_Stats>(MDL_Molecule &,
                                Molecule_to_Query_Specifications & mqs,
                                resizable_array_p<Query_and_Charge_Stats> &, int);

template int read_one_or_more_queries_from_file (resizable_array_p<Query_and_Charge_Stats> & queries, const const_IWSubstring & fname, int verbose);
template int read_one_or_more_queries_from_file (resizable_array_p<Query_and_Charge_Stats> & queries, iwstring_data_source & input, int verbose);

template int queries_from_ISIS_query_file<Query_and_Charge_Stats>(const_IWSubstring const&, Molecule_to_Query_Specifications&, resizable_array_p<Query_and_Charge_Stats>&, int);
template int queries_from_ISIS_query_file<Query_and_Charge_Stats>(data_source_and_type<MDL_Molecule>&, Molecule_to_Query_Specifications&, resizable_array_p<Query_and_Charge_Stats>&, int);
template int queries_from_ISIS_query_file<Query_and_Charge_Stats>(const_IWSubstring const&, resizable_array_p<Query_and_Charge_Stats>&, int);
template int query_from_ISIS_query_file<Query_and_Charge_Stats>(MDL_Molecule&, Molecule_to_Query_Specifications&, resizable_array_p<Query_and_Charge_Stats>&, int);

template int
build_query_from_smiles (const const_IWSubstring & smiles,
                         resizable_array_p<Query_and_Charge_Stats> & queries,
                         int verbose);
