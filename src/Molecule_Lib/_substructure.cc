#include <stdlib.h>

/*
  Used for instantiating templates associated with substructure_Queries.
*/

#include "Foundational/cmdline/cmdline.h"

#include "substructure.h"
#include "rwsubstructure.h"

template int process_queries(Command_Line &, resizable_array_p<Substructure_Query> &, int, char);
template int process_cmdline_token(const char, const const_IWSubstring &, resizable_array_p<Substructure_Query> &, int);

template int queries_from_file(iwstring_data_source &, resizable_array_p<Substructure_Query> &, const IWString &, int);
template int queries_from_file(const_IWSubstring const &, resizable_array_p<Substructure_Query> &, int, int);

template int smarts_from_file(const const_IWSubstring &, resizable_array_p<Substructure_Query> &, int);
template int smarts_from_file(iwstring_data_source &, resizable_array_p<Substructure_Query> &, int);

template int smiles_from_file(const const_IWSubstring &, resizable_array_p<Substructure_Query> &, int);
template int smiles_from_file(iwstring_data_source &, resizable_array_p<Substructure_Query> &, int);

template int process_files_of_queries(Command_Line &, resizable_array_p<Substructure_Query> &, int, int, char);

template int file_record_is_file(resizable_array_p<Substructure_Query> & queries, const IWString & directory_path, IWString & buffer, int verbose);

template int file_record_is_smarts(resizable_array_p<Substructure_Query> & queries, IWString & buffer, int verbose);

template int queries_from_file_of_molecules<Substructure_Query>(const_IWSubstring const &, resizable_array_p<Substructure_Query> &, int);
template int queries_from_file_of_molecules<Substructure_Query>(data_source_and_type<MDL_Molecule> &, Molecule_to_Query_Specifications  &, resizable_array_p<Substructure_Query> &, int);
template int queries_from_file_of_molecules<Substructure_Query>(MDL_Molecule &, Molecule_to_Query_Specifications  &, resizable_array_p<Substructure_Query> &, int);
template int queries_from_file_of_molecules<Substructure_Query>(const_IWSubstring const&, Molecule_to_Query_Specifications &, resizable_array_p<Substructure_Query>&, int);

template int read_one_or_more_queries_from_file(resizable_array_p<Substructure_Query> & queries, const const_IWSubstring & fname, int verbose);
template int read_one_or_more_queries_from_file(resizable_array_p<Substructure_Query> & queries, iwstring_data_source & input, int verbose);

template int queries_from_ISIS_query_file<Substructure_Query>(const_IWSubstring const&, Molecule_to_Query_Specifications&, resizable_array_p<Substructure_Query>&, int);
template int queries_from_ISIS_query_file<Substructure_Query>(data_source_and_type<MDL_Molecule>&, Molecule_to_Query_Specifications&, resizable_array_p<Substructure_Query>&, int);
template int queries_from_ISIS_query_file<Substructure_Query>(const_IWSubstring const&, resizable_array_p<Substructure_Query>&, int);
template int query_from_ISIS_query_file<Substructure_Query>(MDL_Molecule&, Molecule_to_Query_Specifications&, resizable_array_p<Substructure_Query>&, int);
template int build_query_from_smiles<Substructure_Query>(const const_IWSubstring & smiles, resizable_array_p<Substructure_Query> & queries, int verbose);
