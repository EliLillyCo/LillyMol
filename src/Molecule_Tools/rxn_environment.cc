#include <iostream>
#include <memory>
using std::cerr;
using std::endl;

#include "Foundational/iwmisc/misc.h"

#include "rxn_environment.h"

#define MAX_RADIUS_KEY "_MAX_RADIUS"
#define NUMBER_REACTIONS_KEY "_NREACTIONS"

RXN_Environment::RXN_Environment()
{
  _max_radius = 0;
  _nreactions = 0;
  _found = nullptr;

  return;
}

RXN_Environment::~RXN_Environment()
{
  if (nullptr != _found)
    delete [] _found;

  return;
}

int
RXN_Environment::initialise(const int s)
{
  _max_radius = s;
  
  assert (nullptr == _found);

  _found = new std::unordered_map<uint64_t,int>[_max_radius + 1];

  return 1;
}

static int
shortest_distance_to_changing_atom(Molecule & m,
                                   const int * changed,
                                   const atom_number_t zatom)
{
  const int matoms = m.natoms();

  int mindist = matoms;

  for (int i = 0; i < matoms; ++i)
  {
    if (0 == changed[i])
      continue;

    if (i == zatom)
    {
      return 0;
    }

    const int d = m.bonds_between(zatom, i);

    if (d < mindist)
      mindist = d;
  }

  return mindist;
}

static int
bond_numeric(const Bond * b)
{
  int rc = 0;
  if (b->is_aromatic())
    rc = 97;
  else if (b->is_single_bond())
    rc = 11;
  else if (b->is_double_bond())
    rc = 24;
  else if (b->is_triple_bond())
    rc = 33;

  if (b->nrings())
    rc += 25;

  return rc;
}


uint64_t
RXN_Environment::_form_hash(Molecule & m,
                            const Bond * b,
                            const int * atype,
                            const int * changed,
                            int & d)  const
{
  const atom_number_t a1 = b->a1();
  const atom_number_t a2 = b->a2();

  const int d1 = shortest_distance_to_changing_atom(m, changed, a1);
  if (d1 > (_max_radius + 1))
    return 0;

  const int d2 = shortest_distance_to_changing_atom(m, changed, a2);

  if (d1 > _max_radius && d2 > _max_radius)
    return 0;

  uint64_t x = static_cast<uint64_t>(atype[a1]) * static_cast<uint64_t>(atype[a2]);

  x += 9 * bond_numeric(b);

  d = std::min(d1, d2);

  return x;
}

int
RXN_Environment::gather(ISIS_RXN_FILE_Molecule & m,
                        const int * atype,
                        const int * changed)
{
  m.recompute_distance_matrix();
  m.ring_membership();

  assert (nullptr != _found);

  const int nb = m.nedges();

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = m.bondi(i);

    int d;
    const uint64_t x = _form_hash(m, b, atype, changed, d);
    if (0 == x)
      continue;

    _found[d][x]++;
  }

  _nreactions++;

  return 1;
}


template <typename T>
int
RXN_Environment::_do_store_value(Db & db,
                                 const char * key,
                                 const T & v) const
{
  Dbt dbkey((char *) key, ::strlen(key));
  Dbt dbdata((char *) &v, sizeof(v));

  const int rc = db.put(nullptr, &dbkey, &dbdata, 0);

  if (0 == rc)
    return 1;

  cerr << "RXN_Environment::_do_store_value:cannot store, " << key << ' ';
  db.err(rc, "");
  return 0;
}

int
RXN_Environment::do_store(Db & db) const
{
  for (int i = 0; i <= _max_radius; ++i)
  {
    for (auto j : _found[i])
    {
      _do_store(db, i, j);
    }
  }

  if (! _do_store_value(db, MAX_RADIUS_KEY, _max_radius))
    return 0;

  if (! _do_store_value(db, NUMBER_REACTIONS_KEY, _nreactions))
    return 0;

  Dbt dbkey((char *) MAX_RADIUS_KEY, ::strlen(MAX_RADIUS_KEY));
  Dbt dbdata((char *) &_max_radius, sizeof(_max_radius));

  const int rc = db.put(nullptr, &dbkey, &dbdata, 0);

  if (0 == rc)
    return 1;

  cerr << "RXN_Environment::_do_store:cannot store, max radius ";
  db.err(rc, "");
  return 0;
}

/*
  Store a 9 byte item. Radius int the first byte, the atom types and bond info in the other 8 bytes.
  The value is the count
*/

int
RXN_Environment::_do_store(Db & database,
                           const int radius,
                           const std::pair<uint64_t, int> & f) const
{
  unsigned char buffer[1 + sizeof(uint64_t)];

  buffer[0] = static_cast<unsigned char>(radius);

//cerr << "set first byte " << static_cast<int>(buffer[0]) << endl;

  uint64_t x = f.first;

  ::memcpy(buffer + 1, reinterpret_cast<void *>(&x), sizeof(uint64_t));

  Dbt dbkey((char *) buffer, sizeof(buffer));

  int c = f.second;
  Dbt dbdata((char *)&c, sizeof(c));

//cerr << "JUst about to call put\n";

  const int rc = database.put(nullptr, &dbkey, &dbdata, 0);

  if (0 == rc)
    return 1;

//cerr << "rc " << rc << endl;

  cerr << "RXN_Environment::_do_store:cannot store, radius " << radius << " ";
  database.err(rc, "");
  cerr << endl;
  return 0;
}

int
RXN_Environment::do_read (Db & database)
{
  Dbt dbkey((char *)(MAX_RADIUS_KEY), ::strlen(MAX_RADIUS_KEY));

  Dbt dbvalue;

  auto rc = database.get(nullptr, &dbkey, &dbvalue, 0);
  if (0 != rc)
  {
    cerr << "RXN_Environment::do_read:cannot retrieve radius " << MAX_RADIUS_KEY << "\n";
    return 0;
  }

  ::memcpy(&_max_radius, dbvalue.get_data(), sizeof(_max_radius));

  Dbc * cursor = nullptr;

  rc = database.cursor(nullptr, &cursor, 0);
  if (0 != rc)
  {
    database.err(rc, "cannot acquire cursor");
    return 0;
  }

  const int expected_key_size  = 1 + sizeof(uint64_t);
  const int expected_data_size = sizeof(int);

  Dbt zkey, zdata;

  while (0 == (rc = cursor->get(&zkey, &zdata, DB_NEXT)))
  {
    if (expected_key_size != zdata.get_size() || expected_data_size != zdata.get_size())
    {
      cerr << "RXN_Environment::do_read:skipping invalid size " << zkey.get_size() << " " << zdata.get_size() << "\n";
      continue;
    }

    unsigned char * buffer = reinterpret_cast<unsigned char *>(zkey.get_data());

    const int r = static_cast<int>(buffer[0]);
    if (r < 0 || r > _max_radius)
    {
      cerr << "RXN_Environment::do_read:invalid radius stored " << r << " max is " << _max_radius << endl;
      return 0;
    }

    uint64_t x;
    ::memcpy(&x, buffer + 1, sizeof(x));

    buffer = reinterpret_cast<unsigned char *>(zdata.get_data());
    int c;

    ::memcpy(&c, buffer + 1 + sizeof(uint64_t), sizeof(c));

    _found[r][x] = c;
  }

  return 1;
}

int
RXN_Environment::summary(std::ostream & output) const
{
  output << "RXN_Environment::summary: _nreactions " << _nreactions << " reactions,_max_radius " << _max_radius << endl;

  for (int i = 0; i <= _max_radius; ++i)
  {
    output << " r " << i << " " << _found[i].size() << " features,";

    int max_count = 0;
    int tot = 0;
    for (auto x : _found[i])
    {
      tot += x.second;

      if (x.second > max_count)
        max_count = x.second;
    }

    output << " max count " << max_count << ' ';
    if (_nreactions > 0)
      output << (static_cast<float>(max_count) / static_cast<float>(_nreactions)) << " per reaction\n";
    else
      output << '\n';
  }

  return 1;
}

int
RXN_Environment::assess(Molecule & m,
                        const int * atype,
                        const int * changed,
                        RXN_Environment_Assessment & rxnenva)  const
{
  rxnenva.reset();

  if (nullptr == atype)
    return 0.0f;

  m.recompute_distance_matrix();
  m.ring_membership();

  assert (nullptr != _found);

  const int nb = m.nedges();

  int notfound = 0;
  int found = 0;

  float max_ratio = 0.0f;

  for (int i = 0; i < nb; ++i)
  {
    const Bond * b = m.bondi(i);

    int d;    // not used here
    const uint64_t x = _form_hash(m, b, atype, changed, d);
    if (0 == x)
    {
      notfound++;
      continue;
    }

    int r;
    const int f = _fetch_count(x, r);

    if (0 == f)
    {
      notfound++;
      continue;
    }

    found++;

    const float y = static_cast<float>(f) / static_cast<float>(_nreactions);

    rxnenva.extra_fraction(y);

    if (y > max_ratio)
      max_ratio = y;
  }

  rxnenva.set_found(found, notfound);

  if (0 == found)     // very unlikely
    return 0.0f;

  return max_ratio * static_cast<float>(found) / static_cast<float>(notfound + found);   // heuristic
}

int
RXN_Environment::_fetch_count(const uint64_t x,
                              int & r) const
{
  for (int i = 0; i <= _max_radius; ++i)
  {
    const auto f = _found[i].find(x);
    if (f == _found[i].end())
      continue;

    r = i;
    return f->second;
  }

  return 0;
}

int
RXN_Environment::assess(Molecule & m,
                        const int * atype,
                        const Set_of_Atoms * e,
                        RXN_Environment_Assessment & rxnenva)  const
{
  const int matoms = m.natoms();

  int * changed = new_int(matoms);std::unique_ptr<int[]> free_changed(changed);

  e->set_vector(changed, 1);

  return assess(m, atype, changed, rxnenva);
}

void
RXN_Environment_Assessment::_default_values()
{
  _fraction_found = 0.0f;

  _found = 0;
  _notfound = 0;

  return;
}

RXN_Environment_Assessment::RXN_Environment_Assessment()
{
  _default_values();
}

void
RXN_Environment_Assessment::reset()
{
  _default_values();

  _acc_fraction.reset();

  return;
}

template <typename T>
int
RXN_Environment_Assessment::do_print(T &output) const
{
  static const char output_separator = ' ';

  if (0 == _acc_fraction.n())
  {
    output << output_separator << '0' << output_separator << '0';
    return 1;
  }

  output << output_separator << _acc_fraction.maxval() << output_separator << static_cast<float>(_acc_fraction.average()) << 
            output_separator << (static_cast<float>(_found) / static_cast<float>(_found + _notfound));

  return 1;
}
