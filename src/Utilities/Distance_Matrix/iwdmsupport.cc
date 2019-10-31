#include <stdlib.h>

#include "iwdmsupport.h"

#include "iwstring_data_source.h"

void
write_index_and_id (int i,
                    const IWString * _id,
                    std::ostream & output)
{
  output << i;

  if (NULL == _id)
    return;

  output << " '" << *(_id + i) << "'";

  return;
}

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");

void
write_smiles_and_id (const IWString & id,
                     const IW_STL_Hash_Map_String & id_to_smiles,
                     std::ostream & output)
{
  IW_STL_Hash_Map_String::const_iterator f = id_to_smiles.find (id);

  if (f == id_to_smiles.end () && id.contains (' '))
  {
    IWString tmp (id);
    tmp.truncate_at_first (' ');
    f = id_to_smiles.find (tmp);
  }

  if (f == id_to_smiles.end ())
    output << smiles_tag << "C>\n";
  else
    output << smiles_tag << (*f).second << ">\n";

  output << identifier_tag << id << ">\n";

  return;
}

template <typename T>
int
parse_directive (const IWString & n,
                 T & v)
{
  const_IWSubstring d, sv;
  if (! n.split (d, '=', sv))
  {
    cerr << "Directives must be of the form 'XXX=YYY', '" << n << "' is invalid\n";
    return 0;
  }

  return sv.numeric_value (v);
}

template int parse_directive (const IWString &, float &);
template int parse_directive (const IWString &, int &);
template int parse_directive (const IWString &, double &);
template int parse_directive (const IWString &, unsigned char &);
