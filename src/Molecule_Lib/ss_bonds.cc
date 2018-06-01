/*
  Functions for Substructure_Bond objects
*/

#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include "misc.h"

#include "substructure.h"
#include "target.h"


Substructure_Bond_Specifier_Base::Substructure_Bond_Specifier_Base()
{
  _next = NULL;

  return;
}

Substructure_Bond_Specifier_Base::~Substructure_Bond_Specifier_Base()
{
  if (NULL != _next)
    delete _next;

  return;
}

void
Substructure_Bond_Specifier_Base::add_to_chain (Substructure_Bond_Specifier_Base * b)
{
  if (NULL == _next)
    _next = b;
  else
    _next->add_to_chain(b);

  return;
}

/*
  
*/

Substructure_Bond_Specifier_Type::Substructure_Bond_Specifier_Type (bond_type_t b)
{
  assert (SINGLE_BOND == b || DOUBLE_BOND == b || TRIPLE_BOND == b);    // single, double or triple only

  _bond_type = b;

  return;
}

int
Substructure_Bond_Specifier_Type::ok() const
{
  if (SINGLE_BOND == _bond_type)
    return 1;

  if (DOUBLE_BOND == _bond_type)
    return 1;

  if (TRIPLE_BOND == _bond_type)
    return 1;

  return 0;
}

int
Substructure_Bond_Specifier_Type::debug_print (std::ostream & os, const IWString & ind) const
{
  os << ind << "Substructure_Bond_Specifier_Type: type ";

  if (SINGLE_BOND & _bond_type)
    os << "single";
  else if (DOUBLE_BOND & _bond_type)
    os << "double";
  else if (TRIPLE_BOND & _bond_type)
    os << "triple";

  if (AROMATIC_BOND & _bond_type)
    os << " aromatic";

  os << endl;

  return os.good();
}

int
Substructure_Bond_Specifier_Type::smarts (std::ostream & os) const
{
  if (SINGLE_BOND & _bond_type)
    os << '-';
  else if (DOUBLE_BOND & _bond_type)
    os << '=';
  else if (TRIPLE_BOND & _bond_type)
    os << '#';
  else if (AROMATIC_BOND & _bond_type)
    os << ':';
  else
    os << "???";

  return os.good();
}

int
Substructure_Bond_Specifier_Type::involves_aromatic_bond_specification (int & r) const
{
  if (NULL == _next)
    return 0;

  return _next->involves_aromatic_bond_specification(r);
}

int
Substructure_Bond_Specifier_Type::matches (Bond_and_Target_Atom & bata) const
{
//cerr << "Substructure_Bond_Specifier_Type::matches: _bond_type " << _bond_type << " target " << bata.btype() <<  " returning " << (0 != ((BOND_TYPE_ONLY_MASK &bata.btype()) & (BOND_TYPE_ONLY_MASK &_bond_type))) << endl;

  return (0 != ((BOND_TYPE_ONLY_MASK & bata.btype()) & (BOND_TYPE_ONLY_MASK & _bond_type)));   // Oct 2014. Otherwise =: (double and aromatic) does not work
//return bata.btype() == _bond_type;
}

Substructure_Bond_Specifier_Ring::Substructure_Bond_Specifier_Ring (boolean r) : _ring(r)
{
}

int
Substructure_Bond_Specifier_Ring::ok() const
{
  return (_ring >= 0);
}

int
Substructure_Bond_Specifier_Ring::debug_print (std::ostream & os, const IWString & ind) const
{
  os << ind << "Substructure_Bond_Specifier_Ring::ring = " << _ring << endl;

  return os.good();
}

int
Substructure_Bond_Specifier_Ring::smarts (std::ostream & os) const
{
  if (_ring)
    os << '@';
  else
    os << "!@";

  return os.good();
}

int
Substructure_Bond_Specifier_Ring::involves_aromatic_bond_specification (int & r) const
{
  r++;

  if (NULL == _next)
    return 0;

  return _next->involves_aromatic_bond_specification(r);
}

int
Substructure_Bond_Specifier_Ring::matches (Bond_and_Target_Atom & bata) const
{
  boolean r = boolean(bata.nrings() > 0);

//cerr << "Substructure_Bond_Specifier_Ring::matches: rings for bond " << r << " from " << bata.nrings() << " must match " << _ring << " returning " << (_ring == r) << endl;

  return (_ring == r);
}

Substructure_Bond_Specifier_NRings::Substructure_Bond_Specifier_NRings (int r) : _nrings (r)
{
  assert (r >= 0);

  return;
}

int
Substructure_Bond_Specifier_NRings::ok() const
{
  if (_nrings < 0)
    return 0;

  return 1;

}

int
Substructure_Bond_Specifier_NRings::debug_print (std::ostream & os, const IWString & ind) const
{
  os << ind << "Substructure_Bond_Specifier_NRings: nrings = " << _nrings << endl;

  return os.good();
}

int
Substructure_Bond_Specifier_NRings::smarts (std::ostream & os) const
{
  return os.good();
}

int
Substructure_Bond_Specifier_NRings::involves_aromatic_bond_specification (int & r) const
{
  r = 1;
  if (NULL == _next)
    return 0;

  return _next->involves_aromatic_bond_specification(r);
}

int
Substructure_Bond_Specifier_NRings::matches (Bond_and_Target_Atom & bata) const
{
  int r = bata.nrings();

  return (_nrings == r);
}

Substructure_Bond_Specifier_Aromatic::Substructure_Bond_Specifier_Aromatic (boolean a) : _aromatic(a)
{
  return;
}

int
Substructure_Bond_Specifier_Aromatic::ok() const
{
  return 1;
}

int
Substructure_Bond_Specifier_Aromatic::involves_aromatic_bond_specification (int & r) const
{
  return 1;
}

int
Substructure_Bond_Specifier_Aromatic::debug_print (std::ostream & os, const IWString & ind) const
{
  os << ind << "Substructure_Bond_Specifier_Aromatic: matches " << _aromatic << endl;

  return os.good();
}

int
Substructure_Bond_Specifier_Aromatic::smarts (std::ostream & os ) const
{
  if (! _aromatic)
    os << "!:";
  else
    os << ':';

  return os.good();
}

int
Substructure_Bond_Specifier_Aromatic::matches (Bond_and_Target_Atom & bata) const
{
//cerr << "Substructure_Bond_Specifier_Aromatic::matches: _aromatic = " << _aromatic << " target = " << bata.aromatic() << endl;

  boolean ar = bata.aromatic();

  return _aromatic == ar;
}

void
Substructure_Bond::_default_values()
{
  _a1 = NULL;

  _b = NULL;

  _bond_types = 0;

  return;
}

Substructure_Bond::Substructure_Bond()
{
  _default_values();

  return;
}

int
Substructure_Bond::ok() const
{
  if (_a1 && ! _a1->ok())
    return 0;

  if (NULL != _b)
    return _b->ok();

  return 1;
}

int 
Substructure_Bond::debug_print (std::ostream & os, 
                                const IWString & indentation) const
{
  assert (os.good());

  os << indentation << "Bond Descriptor:";

  if (! ok())
    os << indentation << " ** Warning, ok fails **";

  if (0 != _bond_types)
  {
    if (SINGLE_BOND & _bond_types)
      os << "single";
    if (DOUBLE_BOND & _bond_types)
      os << "double";
    if (TRIPLE_BOND & _bond_types)
      os << "triple";
    if (AROMATIC_BOND & _bond_types)
      os << "aromatic";
    os << '\n';
  }

  os << indentation;
  _logexp.debug_print(os);

  Substructure_Bond_Specifier_Base * b = _b;
  os << indentation << "Start bond info " << b << endl;
  while (NULL != b)
  {
    b->debug_print(os, indentation);
    b = b->next();
  }

  return os.good();
}

void
Substructure_Bond::set_atom (Substructure_Atom * a)
{
  assert (ok());

  assert (a && a->ok());

  _a1 = a;

  return;
}

/*
  Aromatic bonds (type 4) must be treated very carefully.
  We set bond type 4 - in case the target object is returning
  type 4 for aromatic bonds, but we must also set a single
  aromatic bond in the alternates
*/

/*int
Substructure_Bond::add_type (int bt)
{
  assert (bt >= 0 && bt <= 4);

  assert (NULL == _b);

  if (0 == bt)
  {
    set_match_any();

    return 1;
  }

  _btype[bt] = 1;

  if (4 == bt)
  {
    assert (NULL == _b);
    _b = new Substructure_Bond_Specifier_Aromatic(1);
  }

  return 1;
}*/

void
Substructure_Bond::set_match_any()
{
  assert (ok());

  assert (NULL == _b);

  _bond_types = (SINGLE_BOND | DOUBLE_BOND | TRIPLE_BOND | AROMATIC_BOND);

  return;
}

int
Substructure_Bond::bond_type_as_string (IWString & zresult) const
{
  assert (ok());

  assert (NULL == _b);    // can only be a simple query

  int rc = 0;

  if (SINGLE_BOND & _bond_types)
  {
    zresult = "single";
    rc = 1;
  }
  
  if (DOUBLE_BOND & _bond_types)
  {
    if (rc)
      zresult << "|double";
    else
      zresult = "double";
    rc++;
  }
  
  if (TRIPLE_BOND & _bond_types)
  {
    if (rc)
      zresult << "|triple";
    else
      zresult = "triple";
    rc++;
  }
  
  if (AROMATIC_BOND & _bond_types)
  {
    if (rc)
      zresult << "|aromatic";
    else
      zresult = "aromatic";
    rc++;
  }
  
  if (0 == rc)
  {
    cerr << "Substructure_Bond::bond_type_as_string: huh, nothing set\n";
    debug_print(cerr, "");
    zresult = "???";
    return 0;
  }

  return rc;
}

int
Substructure_Bond::involves_aromatic_bond_specification (int & need_rings) const
{
  assert (ok());

  if (AROMATIC_BOND & _bond_types)
    return 1;

  if (NULL == _b)    // just matches bond types
    return 0;

  return _b->involves_aromatic_bond_specification(need_rings);
}

/*
  In a smarts, a common thing is to omit the bond, which defaults to
  single or aromatic
*/

int
Substructure_Bond::make_single_or_aromatic()
{
  assert (ok());
  assert (NULL == _b);

  _bond_types = (SINGLE_BOND | AROMATIC_BOND);

  return 1;
}

//#define DEBUG_BOND_MATCH

/*
  Bond matching function. If we are just a bond type, catch that first.
  Otherwise, go an evaluate the logical expression
*/

int
Substructure_Bond::matches (Bond_and_Target_Atom & bata)
{
#ifdef DEBUG_BOND_MATCH
  cerr << "Checking bond to atom " << bata.other()->atom_number() << " type(s) " << bata.btype() << " my type(s) = " << _bond_types << " match " << (_bond_types & bata.btype()) << endl;
#endif

  if (0 == _bond_types)    // no types set, no need to check
    ;
  else if (0 == ((_bond_types) & (bata.btype())))    // mismatch on types specified
    return 0;

#ifdef DEBUG_BOND_MATCH
  if (NULL == _b)
    cerr << "NO components to check\n";
  else
    cerr << "Checking " << _logexp.number_results() << " components\n";
#endif

  if (NULL == _b)
    return 1;

  _logexp.reset();

  int i = 0;     // which result are we generating

  Substructure_Bond_Specifier_Base * b = _b;

  while (NULL != b)
  {
    if (! _logexp.result_needed(i))
    {
      i++;
      b = b->next();
      continue;
    }

    int result = b->matches(bata);

    _logexp.set_result(i, result);

#ifdef DEBUG_BOND_MATCH
    cerr << "Bond component " << i << " matching bond ";
    b->debug_print(cerr);
    cerr << " match is " << result;
#endif

    int rc;
    if (_logexp.evaluate(rc))
    {
#ifdef DEBUG_BOND_MATCH
      cerr << " expression is complete: rc = " << rc << endl;
      _logexp.debug_print(cerr);
#endif

      return rc;
    }

#ifdef DEBUG_BOND_MATCH
    cerr << " expression not complete yet\n";
#endif

    i++;
    b = b->next();
  }

  assert (NULL == "Substructure_Bond_Specifier_Base::matches: should not come to here");

  return 1;
}

static int
count_bond_charactes (const char * smarts,
                      int chars_to_process)
{
  for (int i = 0; i < chars_to_process; i++)
  {
    char s = smarts[i];

    if ('-' == s)
      ;
    else if ('=' == s)
      ;
    else if ('#' == s)
      ;
    else if (':' == s)
      ;
    else if ('~' == s)
      ;
    else if ('&' == s)
      ;
    else if (',' == s)
      ;
    else if ('^' == s)
      ;
    else if (';' == s)
      ;
    else if ('@' == s)
      ;
    else if ('!' == s)
      ;
    else
      return i;
  }

  return chars_to_process;
}

/*
  We make some extensions to smarts.
  @@ means two rings, @@@ means 3 rings, etc...

  We return the number of characters we process
*/

static int
fetch_ring_specs (const char * smarts,
                  int max_chars,
                  int & nrings)
{
  assert ('@' == smarts[0]);

  if (1 == max_chars)
  {
    nrings = 1;
    return 1;
  }

  if ('@' != smarts[1])
  {
    nrings = 1;
    return 1;
  }

  nrings = 2;
  for (int i = 2; i < max_chars; i++)
  {
    if ('@' != smarts[i])
      return i;

    nrings++;
  }

  return max_chars;
}

static int
fetch_ring_specs (const char * smarts,
                  int characters_processed,
                  int chars_to_process,
                  int & nr)
{
  return fetch_ring_specs(smarts + characters_processed, chars_to_process - characters_processed, nr);
}
       
static bond_type_t
char_to_btype (char c)
{
  if ('-' == c)
    return SINGLE_BOND;

  if ('=' == c)
    return DOUBLE_BOND;

  if ('#' == c)
    return TRIPLE_BOND;

  if (':' == c)
    return AROMATIC_BOND;

  return 0;
}

//#define DEBUG_BOND_FROM_SMARTS

/*
  We need to be specially careful with aromatic bonds.
*/

int
Substructure_Bond::_construct_from_smarts (const char * smarts,
                             int chars_to_process,
                             int & characters_processed)
{
  int nchar = count_bond_charactes(smarts, chars_to_process);

#ifdef DEBUG_BOND_FROM_SMARTS
  cerr << "Substructure_Bond::_construct_from_smarts: processing '";
  for (int i = characters_processed; i < chars_to_process; i++)
  {
    cerr << smarts[i];
  }
  cerr << "', nchars = " << nchar << endl;
#endif

  if (0 == nchar)    // take the default of single or aromatic
  {
    make_single_or_aromatic();
    return 1;
  }

// The special case of match any bond

  if (1 == nchar && '~' == smarts[0])
  {
    set_match_any();
    characters_processed = 1;
    return 1;
  }

  if (1 == nchar)
  {
    bond_type_t bt = char_to_btype(smarts[0]);
    if (0 != bt)
    {
      _bond_types = bt;
      characters_processed = 1;

      return 1;
    }
  }

// For efficiency, try to process the special case of 'type,type'


  if (3 == nchar && ',' == smarts[1])
  {
    bond_type_t bt1 = char_to_btype(smarts[0]);
    bond_type_t bt2 = char_to_btype(smarts[2]);

    if (0 != bt1 && 0 != bt2)
    {
      _bond_types = bt1 | bt2;

      characters_processed = 3;

      return 1;
    }
  }

// The special case of '!type' is also easy sometimes. We don't do ! aromatic here

  if (2 == nchar && '!' == smarts[0])
  {
    bond_type_t bt = char_to_btype(smarts[1]);
    if (0 != bt && AROMATIC_BOND != bt)
    {
      if (SINGLE_BOND == bt)
        _bond_types = DOUBLE_BOND | TRIPLE_BOND;
      else if (DOUBLE_BOND == bt)
        _bond_types = SINGLE_BOND | TRIPLE_BOND;
      else if (TRIPLE_BOND == bt)
        _bond_types = SINGLE_BOND | DOUBLE_BOND;
      else
      {
        cerr << "Substructure_Bond::_construct_from_smarts: what is this '" << smarts[0] << smarts[1] << "'\n";
        return 0;
      }

      characters_processed = 2;

      return 1;
    }
  }

  if (';' == smarts[0] || '&' == smarts[0] || ',' == smarts[0] || '^' == smarts[0])
  {
    cerr << "Substructure_Bond::_construct_from_smarts:cannot start bond smarts with '" << smarts[0] << "'\n";
    return 0;
  }

#ifdef DEBUG_BOND_FROM_SMARTS
  cerr << "Substructure_Bond::_construct_from_smarts: building operators\n";
#endif

// All other cases are handled the hard way

  int logop = IW_LOGEXP_UNDEFINED;

  int natt = 0;    // the number of components of the logical expression

  boolean previous_token_was_operator = false;

  while (characters_processed < chars_to_process)
  {
    int unary_op = 1;

    boolean token_is_operator = false;

    if ('!' == smarts[characters_processed])
    {
      unary_op = 0;
      characters_processed++;

      if (characters_processed >= chars_to_process)
      {
        cerr << "Substructure_Bond::_construct_from_smarts: smarts cannot end with '!'\n";
        return 0;
      }
    }

    char s = smarts[characters_processed];

    Substructure_Bond_Specifier_Base * b = NULL;

    unsigned int bt = char_to_btype(s);

#ifdef DEBUG_BOND_FROM_SMARTS
    cerr << "Substructure_Bond_Specifier_Base::_construct_from_smarts: char '" << s << "' btype " << bt << ", characters_processed " << characters_processed << endl;
#endif

    if (SINGLE_BOND == bt || DOUBLE_BOND == bt || TRIPLE_BOND == bt)
    {
      b = new Substructure_Bond_Specifier_Type(bt);     // can be a memory leak at times, fix sometime...
      characters_processed++;
    }
    else if (':' == s)
    {
      if (0 == unary_op)
        b = new Substructure_Bond_Specifier_Aromatic(0);
      else
        b = new Substructure_Bond_Specifier_Aromatic(1);

       unary_op = 1;
       characters_processed++;
    }
    else if ('@' == s)
    {
      int nr;
      int nch = fetch_ring_specs(smarts, characters_processed, chars_to_process, nr);

//    Deal with some common special cases

      if (1 == nch)
      {
        if (0 == unary_op)
        {
          b = new Substructure_Bond_Specifier_Ring(0);
          unary_op = 1;
        }
        else
          b = new Substructure_Bond_Specifier_Ring(1);
      }
      else    // must be number of rings specifier
      {
        b = new Substructure_Bond_Specifier_NRings(nr);
      }

      characters_processed += nch;
    }
    else if ('&' == s)
    {
      logop = IW_LOGEXP_AND;
      characters_processed++;
      token_is_operator = true;
    }
    else if (',' == s)
    {
      logop = IW_LOGEXP_OR;
      characters_processed++;
      token_is_operator = true;
    }
    else if (';' == s)
    {
      logop = IW_LOGEXP_LOW_PRIORITY_AND;
      characters_processed++;
      token_is_operator = true;
    }
    else if ('^' == s)
    {
      logop = IW_LOGEXP_XOR;
      characters_processed++;
      token_is_operator = true;
    }
    else if ('~' == s)   // does not restrict anything, so just skip. This is wrong, the logical expression gets out of sync with the bonds, fix sometime
    {
      cerr << "Substructure_Bond::_construct_from_smarts:~ specification as part of composite makes no sense, ignored. Beware!\n";
      characters_processed++;
      continue;
    }
    else if (0 == unary_op)
    {
      cerr << "Substructure_Bond::_construct_from_smarts: unary operator '!' not followed by bond specifier\n";
      return 0;
    }
    else      // must be done with this bond
    {
      if (NULL != b)
        delete b;
      break;
    }

    if (token_is_operator)
    {
      assert (IW_LOGEXP_UNDEFINED != logop);
      if (0 == unary_op)
      {
        cerr << "Substructure_Bond::_construct_from_smarts: operator cannot follow unary operator\n";
        return 0;
      }
      _logexp.add_operator(logop);
      previous_token_was_operator = true;
    }
    else
    {
      if (NULL == _b)
        _b = b;
      else
      {
        _b->add_to_chain(b);
        if (! previous_token_was_operator)
          _logexp.add_operator(IW_LOGEXP_AND);
      }
      _logexp.set_unary_operator(natt, unary_op);
      natt++;
      previous_token_was_operator = false;
    }
  }

#ifdef DEBUG_BOND_FROM_SMARTS
  cerr << "After parsing bond, characters processed = " << characters_processed << endl;
#endif

  return 1;
}

int
Substructure_Bond::construct_from_smarts (const char * smarts,
                             int chars_to_process,
                             int & characters_processed)
{
  assert (NULL == _b);

  characters_processed = 0;

  if (chars_to_process <= 0)
  {
    cerr << "Substructure_Bond::construct_from_smarts: no chars to process\n";
    return 0;
  }

  int rc = _construct_from_smarts (smarts, chars_to_process, characters_processed);

  return rc;
}

/*
  Get bond info from a Molecule's bond.
  Note that ring info is optional
*/

int
Substructure_Bond::copy (const Bond * b,
                         int copy_bond_attributes)
{
  if (b->is_aromatic())
    _bond_types = AROMATIC_BOND;
  else if (b->is_single_bond())
    _bond_types = SINGLE_BOND;
  else if (b->is_double_bond())
    _bond_types = DOUBLE_BOND;
  else if (b->is_triple_bond())
    _bond_types = TRIPLE_BOND;

  if (copy_bond_attributes)
  {
    cerr << "Substructure_Bond::copy: sorry, copy_bond_attributes not implemented\n";
    return 0;
  }

  return 1;
}

int
Substructure_Bond::set_must_be_in_a_ring (int m)
{
  Substructure_Bond_Specifier_Ring * r = new Substructure_Bond_Specifier_Ring (m);

  if (NULL == _b)
    _b = r;
  else
  {
    _b->add_to_chain(r);
    _logexp.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
  }

  return 1;
}
