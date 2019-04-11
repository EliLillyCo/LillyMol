#include <stdlib.h>

#include "smiles_id_dist.h"

#include "iwstring_data_source.h"

#define MAX_DISTANCE 1.0

static IWString smiles_tag ("$SMI<");
static IWString identifier_tag ("PCN<");
static IWString distance_tag ("DIST<");

int
extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       IWString & zresult)
{
  assert (tag.ends_with ('<'));

  if (! buffer.starts_with (tag))
  {
    cerr << "Tag mismatch. Expecting '" << tag << "' got '" << buffer << "'\n";
    return 0;
  }

  if (! buffer.ends_with ('>'))
  {
    cerr << "Malformed TDT item '" << buffer << "'\n";
    return 0;
  }

  zresult = buffer;
  zresult.remove_leading_chars (tag.length ());
  zresult.chop ();

  return 1;
}

Smiles_ID_Dist::Smiles_ID_Dist ()
{
  _distance = static_cast<similarity_type_t> (2.0);

  return;
}

Smiles_ID_Dist::Smiles_ID_Dist (const IWString & s, const IWString & i, similarity_type_t d) :
    _smiles (s),
    _id (i),
    _distance (d)

{
}

Smiles_ID_Dist::Smiles_ID_Dist (const const_IWSubstring & s, const const_IWSubstring & i, similarity_type_t d) :
    _smiles (s),
    _id (i),
    _distance (d)

{
}

Smiles_ID_Dist &
Smiles_ID_Dist::operator = (const Smiles_ID_Dist & rhs)
{
  _smiles   = rhs._smiles;
  _id       = rhs._id;
  _distance = rhs._distance;

  return *this;
}

int
Smiles_ID_Dist::build (iwstring_data_source & input, int & fatal)
{
  fatal = 0;

  int got_distance = 0;

  const_IWSubstring buffer;

  while (input.next_record (buffer))
  {
    if (buffer.starts_with (smiles_tag))
    {
      if (_smiles.length ())
      {
        cerr << "Smiles_ID_Dist::build: duplicate smiles on line " << input.lines_read () << endl;
        return 0;
      }

      if (! extract_from_tdt_form (buffer, smiles_tag, _smiles))
      {
        cerr << "Smiles_ID_Dist::build: invalid smiles record on line " << input.lines_read () << endl;
        return 0;
      }
    }
    else if (buffer.starts_with (identifier_tag))
    {
      if (_id.length ())
      {
        cerr << "Smiles_ID_Dist::build: duplicate id on line " << input.lines_read () << endl;
        return 0;
      }
      if (! extract_from_tdt_form (buffer, identifier_tag, _id))
      {
        cerr << "Smiles_ID_Dist::build: invalid identifier record on line " << input.lines_read () << endl;
        return 0;
      }
    }
    else if (buffer.starts_with (distance_tag))
    {
      if (got_distance)
      {
        cerr << "Smiles_ID_Dist::build: duplicate distance on line " << input.lines_read () << endl;
        return 0;
      }

      IWString d;
      if (! extract_from_tdt_form (buffer, distance_tag, d))
      {
        cerr << "NN_Item_Base::build: invalid distance form '" << buffer << "', line " << input.lines_read () << endl;
        return 0;
      }

      if (! d.numeric_value (_distance) || _distance < 0.0 || _distance > MAX_DISTANCE)
      {
        cerr << "NN_Item_Base::build: invalid distance '" << buffer << "'\n";
        return 0;
      }

      got_distance = 1;
    }
    else if ('|' == buffer)
    {
      if (_smiles.length () || _id.length () || got_distance)
      {
        cerr << "Smiles_ID_Dist::build: incomplete specification\n";
        if (_smiles.length ())
          cerr << "smiles " << _smiles << endl;
        if (_id.length ())
          cerr << "id " << _id << endl;
        if (got_distance)
          cerr << "dist " << _distance << endl;
        fatal = 1;

        return 0;
      }

      return 0;     // got to the end of the TDT, no information
    }
    else
      continue;

    if (_smiles.length () && _id.length () && got_distance)
      return 1;
  }

  cerr << "Smiles_ID_Dist::build: unexpected EOF\n";
  return 0;
}
