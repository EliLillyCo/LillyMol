#ifndef IWHISTOGRAM_H
#define IWHISTOGRAM_H

class const_IWSubstring;

#include <iostream>
using std::cerr;
using std::endl;

class IWString;

class IWHistogram
{
  protected:
    double _min;
    double _max;
    double _dx;

    int _nbuckets;

    unsigned int _nsamples;

    unsigned int * _count;

//  Normally out-of-range values in extra () are ignored. We can
//  optionally put them in the first or last bucket

    int _put_out_of_range_values_in_first_or_last_bucket;

//  private functions

    int _initialise ();

  public:
    IWHistogram ();
    ~IWHistogram ();

    int ok () const;
    int debug_print (std::ostream &) const;

    void set_put_out_of_range_values_in_first_or_last_bucket (int s) { _put_out_of_range_values_in_first_or_last_bucket = s;}

    int initialise (const const_IWSubstring &);
    int initialise (double mn, double mx, double dx);

    int active () const { return NULL != _count;}

    int nbuckets () const { return _nbuckets;}
    double minval() const { return _min;}
    double maxval() const { return _max;}
    double delta () const { return _dx;}

    int nsamples () const { return _nsamples;}

    const unsigned int * raw_counts () const { return _count;}

    int extra (double);

    int write (std::ostream &) const;
    int write_terse (std::ostream &, int = 1) const;

    int reset ();

    IWHistogram & add (const IWHistogram &);

    int write_as_daylight_fingerprint (const IWString &, std::ostream &) const;

    int which_bucket (double, int &) const;

    int number_samples_less_than (double) const;

//  if someone wants to process only the filled buckets

    int highest_filled_bucket() const;
};

class Resizable_Histogram : public IWHistogram
{
  private:
    double _hard_min;
    double _hard_max;

  public:
    Resizable_Histogram ();

    int ok () const;

    int set_min (double);
    int set_max (double);

    int extra (double);     // overload so we can alter the ranges if needed
};

/*
  In a Fuzzy Histogram, every time a value is encountered, we set its bucket,
  together with _bucket_plus and _bucket_minus either size
*/

class Fuzzy_Histogram : public IWHistogram
{
  private:
    int _bucket_minus;
    int _bucket_plus;

  public:
    Fuzzy_Histogram ();

    void set_bucket_minus (int b) { _bucket_minus = b;}
    void set_bucket_plus (int b) { _bucket_plus = b;}

    int extra (double);
};

class Histogram_with_Overlapping_Buckets : public IWHistogram
{
  private:
    double _delta;

  public:
    Histogram_with_Overlapping_Buckets ();

    void set_delta (double d) { _delta = d;}

    int extra (double);
};

#endif
