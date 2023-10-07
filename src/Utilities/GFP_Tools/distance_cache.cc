#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string.h>

#include "distance_cache.h"

IW_Distance_Cache_Line::IW_Distance_Cache_Line ()
{
  _id = -1;
  _cache = nullptr;

  return;
}

IW_Distance_Cache::IW_Distance_Cache ()
{
  _cache = nullptr;
  _ncols = 0;
  _nrows = 0;
  _next_open_slot = nullptr;

  return;
}

IW_Distance_Cache::~IW_Distance_Cache ()
{
  if (nullptr != _cache)
    delete [] _cache;

  return;
}

int
IW_Distance_Cache::_compute_width_of_cache_line (int ncols) const
{
  return (sizeof(IW_Distance_Cache_Line) + ncols * sizeof(iw_cache_t));
}

int
IW_Distance_Cache::resize (int ncols, size_t bytes)
{
  _width_of_cache_line = _compute_width_of_cache_line(ncols);

  std::cerr << "Size for each struct " << _width_of_cache_line << std::endl;

  _nrows = bytes / _width_of_cache_line;

  if (_nrows <= 0)
  {
    std::cerr << "IW_Distance_Cache::initialise:ncols " << ncols << " but only " << bytes << " bytes, not big enough\n";
    return 0;
  }

  if (nullptr != _cache)
    delete [] _cache;

  _cache = new iw_cache_t[bytes];
  _next_open_slot = reinterpret_cast<IW_Distance_Cache_Line *> (_cache);

  _ncols = ncols;

  unsigned char * v = reinterpret_cast<unsigned char *>(_cache);

  size_t step = _width_of_cache_line;

  for (int i = 0; i < _nrows; i++, v += _width_of_cache_line)
  {
    IW_Distance_Cache_Line * l = reinterpret_cast<IW_Distance_Cache_Line *>(v);

    l->set_id(-1);

    l->set_cache(v + sizeof(IW_Distance_Cache_Line));
  }

  return 1;
}

/*
  Mark all lines empty
*/

int
IW_Distance_Cache::clear ()
{
  assert (nullptr != _cache);

  unsigned char * v = reinterpret_cast<unsigned char *>(_cache);

  size_t step = _width_of_cache_line;

  for (int i = 0; i < _nrows; i++, v += step)
  {
    IW_Distance_Cache_Line * l = reinterpret_cast<IW_Distance_Cache_Line *>(v);
    l->set_id(-1);
  }

  _next_open_slot = reinterpret_cast<IW_Distance_Cache_Line *>(_cache);

  return 1;
}

template <typename T>
int
IW_Distance_Cache::store_line (int id, const T * to_store)
{
  unsigned char * v = reinterpret_cast<unsigned char *>(_next_open_slot);
  unsigned char * c = reinterpret_cast<unsigned char *>(_cache);

  int nstored = (v - c) / _width_of_cache_line;

  if (nstored == _nrows)
    return 0;

  _next_open_slot->store (id, _ncols, to_store);

  v += _width_of_cache_line;

  _next_open_slot = reinterpret_cast<IW_Distance_Cache_Line *>(v);

  return 1;
}

template int IW_Distance_Cache::store_line(int, const float *);

template <typename T>
int
IW_Distance_Cache::retrieve_cache_line (int id, T * f) const
{
  unsigned char * v = reinterpret_cast<unsigned char *>(_cache);

  for (int i = 0; i < _nrows; i++, v += _width_of_cache_line)
  {
    const IW_Distance_Cache_Line * l = reinterpret_cast<IW_Distance_Cache_Line *>(v);

    if (l->id() != id)
      continue;

    l->retrieve(_ncols, f);
    return 1;
  }

  return 0;   // not cached
}

template int IW_Distance_Cache::retrieve_cache_line<float>(int, float*) const;

int
IW_Distance_Cache::cached (int id) const
{
  unsigned char * v = reinterpret_cast<unsigned char *>(_cache);

  for (int i = 0; i < _nrows; i++, v += _width_of_cache_line)
  {
    const IW_Distance_Cache_Line * l = reinterpret_cast<IW_Distance_Cache_Line *>(v);

    if (l->id() == id)
      return 1;
  }

  return 0;
}

int
IW_Distance_Cache::active_lines() const
{
  int rc = 0;

  unsigned char * v = reinterpret_cast<unsigned char *>(_cache);

  for (int i = 0; i < _nrows; i++, v += _width_of_cache_line)
  {
    const IW_Distance_Cache_Line * l = reinterpret_cast<IW_Distance_Cache_Line *>(v);

    if (l->id() >= 0)
      rc++;
  }

  return rc;
}

void
IW_Distance_Cache_Line::store (int id, int n, const iw_cache_t * c)
{
  memcpy(_cache, c, n);
  _id = id;

  return;
}

static unsigned char
convert_to_byte (float f)
{
  return static_cast<unsigned char>(f * 255.49 + 0.4999);
}

void
IW_Distance_Cache_Line::store (int id,
                               int n,
                               const float * f)
{
  for (int i = 0; i < n; i++)
  {
    _cache[i] = convert_to_byte(f[i]);
  }

  _id = id;

  return;
}

void
IW_Distance_Cache_Line::retrieve (int n, iw_cache_t * c) const
{
  memcpy(c, _cache, n);
}

class Just_used_For_initialisation
{
  private:
  public:
    Just_used_For_initialisation();
};

static Just_used_For_initialisation jufi;

static float byte_to_float[256];

Just_used_For_initialisation::Just_used_For_initialisation()
{
  for (int i = 0; i < 256; i++)
  {
    byte_to_float[i] = static_cast<float>(i / 256.0);
  }

  return;
}

void
IW_Distance_Cache_Line::retrieve (int n, float * f) const
{
  for (int i = 0; i < n; i++)
  {
    f[i] = byte_to_float[_cache[i]];
  }

  return;
}
