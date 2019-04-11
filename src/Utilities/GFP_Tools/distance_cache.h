#ifndef IW_DISTANCE_CACHE_H
#define IW_DISTANCE_CACHE_H

typedef unsigned char iw_cache_t;

/*
  If the ID of the cache line is negative, that means this cache line
  is available for use
*/

class IW_Distance_Cache_Line
{
  private:
    int _id;

    iw_cache_t * _cache;

  public:
    IW_Distance_Cache_Line ();

    int id () const { return _id;}
    void set_id (int s) { _id = s;}

    void set_cache (iw_cache_t * c) { _cache = c;}

    void store (int id, int n, const iw_cache_t *);
    void store (int id, int n, const float *);

    void retrieve (int n, iw_cache_t*) const;
    void retrieve (int n, float*) const;
};

class IW_Distance_Cache
{
  private:
    iw_cache_t * _cache;

    int _nrows;

    int _ncols;

    int _width_of_cache_line;

//  since we only fill this sequentially, we can save time by
//  keeping track of next open slot

    IW_Distance_Cache_Line * _next_open_slot;
    
//  private functions

    int _compute_width_of_cache_line (int ncols) const;

  public:
    IW_Distance_Cache();
    ~IW_Distance_Cache();

    int resize (int ncols, size_t bytes);
    int clear ();

    int  cached (int id) const;

    int  active_lines () const;

    template <typename T> int store_line (int id, const T * c);
    template <typename T> int retrieve_cache_line (int id, T * f) const;
};

#endif

