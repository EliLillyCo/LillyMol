#ifndef IWGFP_H
#define IWGFP_H

#include <sys/types.h>
#include "iwstring.h"
#include "set_or_unset.h"
#include "misc.h"
#include "cmdline.h"

#include "dyfp.h"
#include "sparsefp.h"
#include "multi_conformer.h"

class iwstring_data_source;
class Tversky;

typedef float iwproperty_t;

class IW_TDT;


/*
  All the build_from_contiguous_storage methods have an ambiguity in
  their behaviour. 
  Should they allocate their own arrays for the data, or should they
  just continue to point at the storage given? All these methods
  have an int argument that indicates whether or not they should
  allocate and copy the data.
  Note that if they just point to the data given, make sure
  their destructors are never called!
*/

//#define FB_ENTROPY_WEIGHTED_FPS
#ifdef FB_ENTROPY_WEIGHTED_FPS
#endif

template <typename T>
class Molecular_Properties
{
  protected:
    int _nproperties;

    T * _property;

  public:
    ~Molecular_Properties ();

    int active () const { return NULL != _property;} 

    Molecular_Properties<T> & operator = (const Molecular_Properties<T> &);

    int nproperties() const { return _nproperties;}

    const T * rawdata() const { return _property;}

    void create_subset (const int *);

    T sum () const { return sum_vector(_property, _nproperties);}

    const void * build_from_contiguous_storage (const void *, int);
    void * copy_to_contiguous_storage (void *) const;
    void * copy_to_contiguous_storage_gpu (void *) const;
};

template <typename T>
void *
Molecular_Properties<T>::copy_to_contiguous_storage (void * p) const
{
  memcpy(p, this, sizeof(*this));

  p = reinterpret_cast<unsigned char *>(p) + sizeof(*this);

  copy_vector(reinterpret_cast<T *>(p), _property, _nproperties);

  p = reinterpret_cast<T *>(p) + _nproperties;

  return p;
}

template <typename T>
void *
Molecular_Properties<T>::copy_to_contiguous_storage_gpu (void * p) const
{
  memcpy(p, &_nproperties, sizeof(int));

  p = reinterpret_cast<unsigned char *>(p) + sizeof(int);

  copy_vector(reinterpret_cast<T *>(p), _property, _nproperties);

  p = reinterpret_cast<T *>(p) + _nproperties;

  return p;
}

template <typename T>
const void *
Molecular_Properties<T>::build_from_contiguous_storage (const void * p,
                                int allocate_array)
{
  if (allocate_array && NULL != _property)
    delete [] _property;

  memcpy(this, p, sizeof(*this));

  p = reinterpret_cast<const unsigned char *>(p) + sizeof(*this);

  if (allocate_array)
  {
    _property = new T[_nproperties];
    copy_vector(_property, reinterpret_cast<const T *>(p), _nproperties);
  }
  else
    _property = (T *) (p);   // dangerous cast

  p = reinterpret_cast<const T *>(p) + _nproperties;

  return p;
}

extern int set_number_integer_molecular_properties (int);
extern int set_number_continuous_molecular_properties (int);
extern int number_integer_molecular_properites();
extern int initialise_properties_ratios();

class Molecular_Properties_Integer : public Molecular_Properties<int>
{
  private:
  public:
    Molecular_Properties_Integer ();

    int construct_from_tdt_fp_record (const const_IWSubstring &);

    int natoms () const;
    int nrings () const;
    int aromatic_atoms() const;

    similarity_type_t similarity (const Molecular_Properties_Integer &) const;

    similarity_type_t dice_coefficient (const Molecular_Properties_Integer &) const;

    int bits_in_common (const Molecular_Properties_Integer &) const;
#ifdef FB_ENTROPY_WEIGHTED_FPS
    int construct_from_bit_vector (const IWDYFP &);
#endif
};

class Molecular_Properties_Continuous : public Molecular_Properties<float>
{
  private:
  public:
    Molecular_Properties_Continuous ();

    int construct_from_tdt_record (const const_IWSubstring &);
    int construct_from_descriptor_record (const const_IWSubstring & buffer);

    similarity_type_t similarity  (const Molecular_Properties_Continuous &) const;
    similarity_type_t cartesian_distance (const Molecular_Properties_Continuous &) const;
    similarity_type_t exp1_distance (const Molecular_Properties_Continuous &) const;
    similarity_type_t exp2_distance (const Molecular_Properties_Continuous &) const;
    similarity_type_t dice_coefficient (const Molecular_Properties_Continuous &) const;
};

#include "fixed_size_counted_fingerprint.h"

#include "sparse_collection.h"

class IW_General_Fingerprint
{
  protected:
    IWDYFP * _fingerprint;

    Sparse_Fingerprint * _sparse_fingerprint;

    Molecular_Properties_Integer _molecular_properties_integer;

    Molecular_Properties_Continuous _molecular_properties_continuous;

    typedef Fixed_Size_Counted_Fingerprint_uchar Fixed_Size_Counted_Fingerprint;
//  typedef Fixed_Size_Counted_Fingerprint_uint Fixed_Size_Counted_Fingerprint;

    Fixed_Size_Counted_Fingerprint * _fixed_size_counted_fingerprint;

    Multiconformer_01 * _multiconformer_01;

    Multiconformer_Fixed_Counted * _multiconformer_fixed_size_counted;

    Multiconformer_Sparse * _multiconformer_sparse;

    IWString  _id;

    Set_or_Unset<__off_t> _offset;

//  private functions

    void _default_values ();

    void _allocate_fingerprint_array ();

    int _construct_from_tdt (IW_TDT &, int &);
    int _read_molecular_properties_integer (IW_TDT & tdt, const IWString & tag);
    int _read_molecular_properties_continuous (IW_TDT & tdt, const IWString & tag);

    int _convert_sparse_fingerprints_to_bits (const Set_of_Sparse_Fingerprint_Collection_Profile & sfpcp);
    int _convert_sparse_fingerprints_to_fixed_width_counted (const Set_of_Sparse_Fingerprint_Collection_Profile & sfpcp);

//  similarity_type_t _property_distance (const IW_General_Fingerprint & rhs) const;

  public:
    IW_General_Fingerprint ();
    IW_General_Fingerprint (const IW_General_Fingerprint &);
    ~IW_General_Fingerprint ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int print_all_bits (std::ostream &) const;

    int construct_from_tdt (IW_TDT &, int &);

    void set_id(const const_IWSubstring & s) { _id = s;}

    int nfingerprints () const;

    int number_sparse_fingerprints () const;

//  Both of these functions will die if called when properties haven't been initialised

    int natoms () const { return _molecular_properties_integer.natoms ();}
    int nrings () const { return _molecular_properties_integer.nrings ();}
    int aromatic_atoms() const { return _molecular_properties_integer.aromatic_atoms();}

    const Molecular_Properties_Integer & molecular_properties_integer() const { return _molecular_properties_integer;}
    Molecular_Properties_Integer & molecular_properties_integer() { return _molecular_properties_integer;}

    const IWString & id () const { return _id;}

    // off_t not defined in C++-11, only __off_t is
    int offset (__off_t &) const;
    void set_offset (__off_t o) { _offset.set (o);}

    IWDYFP & operator [] (int) const;
    IWDYFP & item (int) const;

    const Sparse_Fingerprint & sparse_fingerprint (int i) const { return _sparse_fingerprint[i];}
    Sparse_Fingerprint & sparse_fingerprint (int i) { return _sparse_fingerprint[i];}

    IW_General_Fingerprint & operator = (const IW_General_Fingerprint &);

//  Once things are initialised, we can change the fingerprint tags.

    int set_fingerprint_tag (int, const const_IWSubstring &);

    similarity_type_t tanimoto (IW_General_Fingerprint &);
    similarity_type_t tanimoto (IW_General_Fingerprint * rhs) { return tanimoto (*rhs);}
    similarity_type_t tversky  (IW_General_Fingerprint &, const Tversky &);

    similarity_type_t distance  (IW_General_Fingerprint & rhs) { return static_cast<similarity_type_t> (1.0) - tanimoto (rhs);}

//  possibly faster similarity function with a cutoff

    int tanimoto (IW_General_Fingerprint & rhs, similarity_type_t, similarity_type_t &);

    similarity_type_t optimistic_distance (IW_General_Fingerprint & rhs, const Tversky &);
    similarity_type_t optimistic_distance (IW_General_Fingerprint & rhs, tversky_coeff_t, tversky_coeff_t);

//  We can initialise ourselves based on a fingerprint

    int initialise (const IW_TDT &);

//  One member will need to be called to initialise the static members

    int bits_in_fingerprint (int f) const;

    void flip ();     // invert all fingerprints - used by leader algorithm

    void iwor (const IW_General_Fingerprint &);

    int convert_to_non_sparse_forms (const Set_of_Sparse_Fingerprint_Collection_Profile & sfpcp);

//  Used to emulate svm stuff. All bits from all fingerprints are grouped together

    float equal_weight_tanimoto(IW_General_Fingerprint & rhs);

//  For svmfp models, we want a couple of other kernels
  
    float equal_weight_dot_product (IW_General_Fingerprint & rhs);

    int bytes_needed_for_contiguous_storage() const;
    void * copy_to_contiguous_storage (void *) const;
    int bytes_needed_for_contiguous_storage_gpu() const;
    void * copy_to_contiguous_storage_gpu (void *) const;
    const void * build_from_contiguous_storage (const void * p, int);

#ifdef FB_ENTROPY_WEIGHTED_FPS
  int convert_01_fingerprint_to_integer_molecular_properties();
#endif
  int fixed_fingerprints_as_hex (IWString &) const;
};

/*
  We often need a class with a fingerprint and a smiles
*/

class FP_and_Smiles : public IW_General_Fingerprint
{
  protected:
    IWString _smiles;

  public:
    int construct_from_tdt (IW_TDT &, int &);

    const IWString & smiles () const { return _smiles;}
};

/*
  We often need a class which has a distance value associated with it
*/

class IW_GFP_D : public IW_General_Fingerprint
{
  protected:
    similarity_type_t _distance;

  public:
    IW_GFP_D ();

    similarity_type_t   distance () const { return _distance;}
    void set_distance (similarity_type_t d) { _distance = d;}
};

extern void set_identifier_tag (const const_IWSubstring & id);
extern int  display_standard_gfp_options (std::ostream &);

extern int initialise_fingerprints (Command_Line & cl, int);

extern int number_fingerprints ();
extern int number_sparse_fingerprints ();

/*
  Allows us to control what similarity function to use
*/

extern void set_sparse_fingerprint_counts_limited(const int);

/*
  The other way to initialise the fingerprint environment is to
  examine a file and use whatever fingerprints are in there.
*/

extern int initialise_fingerprints (const char *, int);
extern int initialise_fingerprints (const IWString &, int);

extern int build_pool (iwstring_data_source & input,
            IW_General_Fingerprint * pool,
            int max_pool_size,
            int & pool_size);

extern int build_pool (const const_IWSubstring & fname,
            IW_General_Fingerprint * pool,
            int max_pool_size,
            int & pool_size,
            const IWString & identifier_tag);

extern int can_be_compared (const IW_General_Fingerprint & fp1, const IW_General_Fingerprint & fp2);

extern int get_tversky_specifications (Command_Line & cl,
                            char aflag, float & a,
                            char bflag, float & b,
                            int verbose);


extern int need_to_call_initialise_fingerprints (const Command_Line &);

extern void set_convert_sparse_fingerprint_to_bits (int);
extern void set_convert_sparse_fingerprint_to_fixed_width_counted (int);

extern int count_tdts_in_file (iwstring_data_source & input, const IWString & identifier_tag);

extern void set_property_weight_integer (double w);
extern void set_fixed_width_fingerprint_weight (int ndx, double w);
extern void set_sparse_fingerprint_weight (int ndx, double w);

extern void set_report_fingerprint_status (int s);

extern const IWString & fixed_fingerprint_tag (int i);
extern const IWString & sparse_fingerprint_tag (int i);

extern const IWString & property_tag();

extern void delete_gfp_file_scope_static_objects ();

#ifdef FB_ENTROPY_WEIGHTED_FPS
extern int initialise_entropy_weighting_for_fingerprints (Command_Line &, char, int);
#endif

extern int property_based_windows_present();

#define WINDOW_MAX_NATOMS 256
#define WINDOW_MAX_RINGS 20

template <typename F>
class Window_Specification
{
  private:
    int _atom_count_window_present;
    int _lower_atom_count_cutoff[WINDOW_MAX_NATOMS];
    int _upper_atom_count_cutoff[WINDOW_MAX_NATOMS];

    int _ring_count_window_present;
    int _lower_ring_count_cutoff[WINDOW_MAX_RINGS];
    int _upper_ring_count_cutoff[WINDOW_MAX_RINGS];

    int _aromatic_atom_count_window_present;
    int _lower_aromatic_atom_count_cutoff[WINDOW_MAX_NATOMS];
    int _upper_aromatic_atom_count_cutoff[WINDOW_MAX_NATOMS];

//  private functions

    int _set_atom_count_window (const int, const int);
    int _set_ring_count_window (const int, const int);
    int _set_aromatic_atom_count_window (const int, const int);
    int _set_atom_count_delta_window (const int, const int, const int);

  public:
    Window_Specification () {
        _atom_count_window_present = 0;
        _ring_count_window_present = 0;
        _aromatic_atom_count_window_present = 0;
    }

    int build (Command_Line & cl, char, int verbose);

    int atom_count_window_present () const { return _atom_count_window_present;}
    int ring_count_window_present () const { return _ring_count_window_present;}
    int aromatic_atom_count_window_present () const { return _aromatic_atom_count_window_present;}

    void set_atom_count_window_present (int s) { _atom_count_window_present = s;}
    void set_ring_count_window_present (int s) { _ring_count_window_present = s;}
    void set_aromatic_atom_count_window_present (int s) { _aromatic_atom_count_window_present = s;}

    int can_be_compared (const F & fp1, const F & fp2) const;
};

template <typename F>
int
Window_Specification<F>::build (Command_Line & cl,
                                char flag,
                                int verbose)
{
  const_IWSubstring w;
  for (auto i = 0; cl.value(flag, w, i); ++i)
  {
    int minxtra = -1;
    int maxxtra = -1;

    if (w.starts_with("R:"))
    {
      w += 2;
      int window;
      if (! w.numeric_value(window) || window < 0 || window > 100)
      {
        cerr << "The '-" << flag << " R:' option must be followed by a whole percentage\n";
        return 0;
      }

      _set_ring_count_window(window, verbose);
    }
    else if (w.starts_with("A:minextra="))
    {
      w.remove_leading_chars(11);
      if (! w.numeric_value(minxtra) || minxtra < 0)
      {
        cerr << "Invalid min extra atoms 'A:minextra=" << w << "'\n";
        return 0;
      }
      if (verbose)
        cerr << "Will only compare molecules with at least " << minxtra << " extra atoms\n";
    }
    else if (w.starts_with("A:maxextra="))
    {
      w.remove_leading_chars(11);
      if (! w.numeric_value(maxxtra) || maxxtra < 0)
      {
        cerr << "Invalid max extra atoms 'A:maxextra=" << w << "'\n";
        return 0;
      }
      if (verbose)
        cerr << "Will only compare molecules with at most " << maxxtra << " extra atoms\n";
    }
    else if (w.starts_with("A:"))
    {
      w += 2;
      int window;
      if (! w.numeric_value(window) || window < 0 || window > 100)
      {
        cerr << "The '-w A:' option must be followed by a whole percentage\n";
        return 0;
      }

      _set_atom_count_window(window, verbose);
    }
    else if (w.starts_with("a:"))
    {
      w += 2;
      int window;

      if (! w.numeric_value(window) || window < 0)
      {
        cerr << "The '-W a:' option must be followed by a valid number\n";
        return 0;
      }

      if (verbose)
        cerr << "Aromatic atom count window set to " << window << endl;
      
      _set_aromatic_atom_count_window(window, verbose);
    }
    else if ("help" == w)
    {
      cerr << endl;
      return 0;
    }
    else
    {
      cerr << "Unrecognised window specification (-" << flag << ") qualifier '" << w << "'\n";
      return 0;
    }

    if (maxxtra >= 0 || minxtra >= 0)
    {
      if (_atom_count_window_present)
      {
        cerr << "Sorry, cannot have a percentage based atom count window as well as minextra/maxextra\n";
        return 0;
      }
      if (maxxtra >= 0 && minxtra >= 0 && minxtra > maxxtra)
      {
        cerr << "Inconsistent values for max extra atoms " << maxxtra << " and min extr atoms " << minxtra << endl;
        return 0;
      }

      _set_atom_count_delta_window(minxtra, maxxtra, verbose);
    }
  }

  return 1;
}

template <typename F>
int
Window_Specification<F>::_set_atom_count_window (const int w,
                                                 const int verbose)
{
  assert (w > 1 && w <= 100);     // must be a percentage

  for (int i = 1; i < WINDOW_MAX_NATOMS; i++)
  {
    int delta = int (static_cast<float> (i * w) / 100.0);

    if (i - delta < 1)
      _lower_atom_count_cutoff[i] = 1;
    else
      _lower_atom_count_cutoff[i] = i - delta;

    _upper_atom_count_cutoff[i] = i + delta;

    if (verbose > 1)
      cerr << "Molecules with " << i << " atoms compared with " << _lower_atom_count_cutoff[i] << " to " << _upper_atom_count_cutoff[i] << " atoms\n";
  }

  _atom_count_window_present = 1;

  return 1;
}

template <typename F>
int
Window_Specification<F>::can_be_compared (const F & fp1,
                                          const F & fp2) const
{
  if (_atom_count_window_present)
  {
    int na1 = fp1.natoms();
    int na2 = fp2.natoms();

#ifdef DEBUG_CAN_BE_COMPARED
    cerr << "Can we compare '" << fp1.id() << " " << na1 << " (" << lower_atom_count_cutoff[na1] << ',' << upper_atom_count_cutoff[na1] << ") and '" << fp2.id() << " " << na2 << endl; //" (" << lower_atom_count_cutoff[na2] << ',' << upper_atom_count_cutoff[na2] << ")\n";
#endif

    if (na2 < _lower_atom_count_cutoff[na1])
      return 0;

    if (na2 > _upper_atom_count_cutoff[na1])
      return 0;
  }

  if (_ring_count_window_present)
  {
    int nr1 = fp1.nrings();
    int nr2 = fp2.nrings();

#ifdef DEBUG_CAN_BE_COMPARED
    cerr << "Rings?? " << nr1 << " (" << lower_ring_count_cutoff[nr1] << ',' << upper_ring_count_cutoff[nr1] << ") and " << nr2 << " (" << lower_ring_count_cutoff[nr2] << ',' << upper_ring_count_cutoff[nr2] << ")\n";
#endif
    if (nr2 < _lower_ring_count_cutoff[nr1])
      return 0;

    if (nr2 > _upper_ring_count_cutoff[nr1])
      return 0;
  }

  if (_aromatic_atom_count_window_present)
  {
    int na1 = fp1.aromatic_atoms();
    int na2 = fp2.aromatic_atoms();

    if (na2 < _lower_aromatic_atom_count_cutoff[na1])
      return 0;

    if (na2 > _upper_aromatic_atom_count_cutoff[na1])
      return 0;
  }

#ifdef DEBUG_CAN_BE_COMPARED
  cerr << "Yes, can be compared\n";
#endif

  return 1;      // these fingerprints can be compared
}

template <typename F>
int
Window_Specification<F>::_set_atom_count_delta_window (const int min_extra_atoms, 
                                                     const int max_extra_atoms,
                                                     const int verbose)
{
  for (int i = 1; i < WINDOW_MAX_NATOMS; i++)
  {
    if (max_extra_atoms >= 0)
      _upper_atom_count_cutoff[i] = i + max_extra_atoms;
    else
      _upper_atom_count_cutoff[i] = WINDOW_MAX_NATOMS;

    if (min_extra_atoms >= 0)
      _lower_atom_count_cutoff[i] = i + min_extra_atoms;
    else
      _lower_atom_count_cutoff[i] = 1;
  }

  _atom_count_window_present = 1;

  return 1;
}

template <typename F>
int
Window_Specification<F>::_set_aromatic_atom_count_window (const int w,
                                        const int verbose)
{
  for (int i = 0; i < WINDOW_MAX_NATOMS; i++)
  {
    _lower_aromatic_atom_count_cutoff[i] = i - w;
    if (_lower_aromatic_atom_count_cutoff[i] < 0)
      _lower_aromatic_atom_count_cutoff[i] = 0;

    _upper_aromatic_atom_count_cutoff[i] = i + w;
  }

  _aromatic_atom_count_window_present = 1;

  return 1;
}

template <typename F>
int
Window_Specification<F>::_set_ring_count_window (const int w, const int verbose)
{
  assert (w > 1 && w <= 100);     // must be a percentage

  for (int i = 1; i < WINDOW_MAX_RINGS; i++)
  {
    int delta = int (static_cast<float> (i * w) / 100.0);

    if (i - delta < 0)
      _lower_ring_count_cutoff[i] = 0;
    else
      _lower_ring_count_cutoff[i] = i - delta;

    _upper_ring_count_cutoff[i] = i + delta;

    if (verbose > 1)
      cerr << "Molecules with " << i << " rings compared with " << _lower_ring_count_cutoff[i] << " to " << _upper_ring_count_cutoff[i] << " rings\n";
  }

  _ring_count_window_present = 1;

  return 1;
}


#endif
