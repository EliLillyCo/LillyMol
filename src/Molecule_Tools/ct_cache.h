#ifndef CT_CACHE_H
#define CT_CACHE_H

class Molecule;

template <typename A, typename B>
class Bond_Cache_Template
{
  private:
    B _bt;
    B _ring;
    A _other;

  public:
    Bond_Cache_Template();
};

typedef Bond_Cache_Template<unsigned char, unsigned char> Bond_Cache;

template <typename T>
class AT_Cache
{
  private:
    T _atomic_number;
    T _hcount;
    T _ncon;
    T _arom;


  public:
};

template <typename T>
class CT_Cache
{
  private:
    const int _matoms;
    const int _nrings;

    int _items_per_atom;

    T * _zdata;

  public:
    CT_Cache(Molecule & m);
    ~CT_Cache();
};

#ifdef CT_CACHE_IMPLEMENTATION

template <typename T>
CT_Cache<T>::CT_Cache(Molecule & m) : _matoms(m.natoms()), _nrings(m.nrings())
{
  if (0 == _matoms)
  {
    _zdata = nullptr;
    return;
  }

  const int items_per_atomic_info =  (
                                      1 +       // atomic number
                                      1 +       // ncon
                                      1 +       // hcount
                                      1 +       // fragment
                                      1 +       // ring bond count
                                      1);       // arom

  const int items_per_bond = 3;    // bond type, rings, other

  _items_per_atom = items_per_atomic_info + 4 * items_per_bond;

  const int items_needed = _matoms * _items_per_atom;

  _zdata = new T[items_needed];

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < _matoms; ++i)
  {
    const Atom * a = m.atomi(i);

    int offset = i * _items_per_atom;

    const Element * e = a->element();

    if (! e->organic())
    {
      cerr << "CT_Cache::CT_Cache:cannot handle non periodic table elements " << e->symbol() << " in " << m.name() << endl;
      return;
    }

    const int acon = a->ncon();
    if (acon > 4)
    {
      cerr << "CT_Cache::CT_Cache:cannot handle highly connected molecules " << acon << " in " << m.name() << endl;
    }

    _zdata[offset++] = e->atomic_number();
    _zdata[offset++] = acon;
    _zdata[offset++] = m.hcount(i);
    _zdata[offset++] = m.fragment_membership(i);

    const int rbc_offset = offset++;

    _zdata[offset++] = m.is_aromatic(i);

    int rbc = 0;

    for (int j = 0; j < acon; ++j)
    {
      const Bond * b = a->item(j);
      if (b->nrings())
        rbc++;

      if (b->is_aromatic())
        _zdata[offset++] = AROMATIC_BOND;
      else if (b->is_single_bond())
        _zdata[offset++] = SINGLE_BOND;
      else if (b->is_double_bond())
        _zdata[offset++] = DOUBLE_BOND;
      else if (b->is_triple_bond())
        _zdata[offset++] = TRIPLE_BOND;
    }

    _zdata[rbc_offset] = rbc;

  }

  return;
}

template <typename T>
CT_Cache<T>::~CT_Cache()
{
  if (nullptr != _zdata)
    delete [] _zdata;
}

#endif

#endif
