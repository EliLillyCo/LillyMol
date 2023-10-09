#include <iostream>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"

#include "molecule_to_query.h"
#include "qry_wstats.h"
#include "rwsubstructure.h"


template class resizable_array_p<Substructure_Hit_Statistics>;
template class resizable_array_base<Substructure_Hit_Statistics *>;

template int process_queries(Command_Line &, resizable_array_p<Substructure_Hit_Statistics> &, int, char);
template int process_cmdline_token(char, const const_IWSubstring &, resizable_array_p<Substructure_Hit_Statistics> &, int);
template int process_files_of_queries(Command_Line &, resizable_array_p<Substructure_Hit_Statistics> &, int, int, char);
//template int process_cmdline_token(const char, const const_IWSubstring&, resizable_array_p<Substructure_Hit_Statistics>&, int);
//template int process_cmdline_token<Substructure_Hit_Statistics>(char, const_IWSubstring&, resizable_array_p<Substructure_Hit_Statistics>&, int);

template int queries_from_file(iwstring_data_source &, resizable_array_p<Substructure_Hit_Statistics> &, const IWString &, int);
template int queries_from_file(const_IWSubstring const &, resizable_array_p<Substructure_Hit_Statistics> &, int, int);
template int queries_from_file(IWString const &, resizable_array_p<Substructure_Hit_Statistics> &, int, int);

template int smarts_from_file(const const_IWSubstring &, resizable_array_p<Substructure_Hit_Statistics> &, int);
template int smarts_from_file(iwstring_data_source &, resizable_array_p<Substructure_Hit_Statistics> &, int);
template int smarts_or_smiles_from_file(iwstring_data_source &, resizable_array_p<Substructure_Hit_Statistics> &, int);

template int file_record_is_file(resizable_array_p<Substructure_Hit_Statistics> & queries, const IWString & directory_path, IWString & buffer, int verbose);
template int file_record_is_smarts(resizable_array_p<Substructure_Hit_Statistics> & queries, IWString & buffer, int verbose);

template int queries_from_file_of_molecules<Substructure_Hit_Statistics>(const_IWSubstring const &, resizable_array_p<Substructure_Hit_Statistics> &, int);
template int queries_from_file_of_molecules<Substructure_Hit_Statistics>(data_source_and_type<MDL_Molecule> &, Molecule_to_Query_Specifications &, resizable_array_p<Substructure_Hit_Statistics> &, int);
template int queries_from_file_of_molecules<Substructure_Hit_Statistics>(MDL_Molecule &, Molecule_to_Query_Specifications  &, resizable_array_p<Substructure_Hit_Statistics> &, int);
template int queries_from_file_of_molecules<Substructure_Hit_Statistics>(const_IWSubstring const&, Molecule_to_Query_Specifications &, resizable_array_p<Substructure_Hit_Statistics>&, int);

template int read_one_or_more_queries_from_file(resizable_array_p<Substructure_Hit_Statistics> & queries, const const_IWSubstring & fname, int verbose);
template int read_one_or_more_queries_from_file(resizable_array_p<Substructure_Hit_Statistics> & queries, iwstring_data_source & input, int verbose);

template int queries_from_ISIS_query_file<Substructure_Hit_Statistics>(const_IWSubstring const&, resizable_array_p<Substructure_Hit_Statistics>&, int);
template int queries_from_ISIS_query_file<Substructure_Hit_Statistics>(const_IWSubstring const&, Molecule_to_Query_Specifications&, resizable_array_p<Substructure_Hit_Statistics>&, int);
template int queries_from_ISIS_query_file<Substructure_Hit_Statistics>(data_source_and_type<MDL_Molecule>&, Molecule_to_Query_Specifications&, resizable_array_p<Substructure_Hit_Statistics>&, int);
template int query_from_ISIS_query_file<Substructure_Hit_Statistics>(MDL_Molecule&, Molecule_to_Query_Specifications&, resizable_array_p<Substructure_Hit_Statistics>&, int);
template int build_query_from_smiles(const const_IWSubstring & smiles, resizable_array_p<Substructure_Hit_Statistics> & queries, int verbose);

