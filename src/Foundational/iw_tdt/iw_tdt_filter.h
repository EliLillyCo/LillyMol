#ifndef IW_TDT_FILTER_H
#define IW_TDT_FILTER_H

#include "iwcrex.h"
#include "logical_expression.h"

class IW_TDT;

/*
  A TDT_FILTER object yields a yes/no match given a TDT.

  The overall object is a logical expression of various other conditions
*/

class IW_TDT_Filter_Base : public IW_Regular_Expression
{
  protected:

//  If there are multiple dataitems in a TDT, we can possibly
//  process a given one

    int _item_to_process;

    int _column_to_process;

//  When processing clogp data we'd like to truncate perception at the
//  semicolon character

    char _truncate_before_first;

//  protected functions

    int _extract_data (const IW_TDT &, const_IWSubstring &, int which_item_to_retrieve, int = 1);
    int _extract_data (const IW_TDT &, double &, int &, int which_item_to_retrieve);
    int _extract_numeric (const const_IWSubstring & s, const char * op, double & zresult);
    int _extract_numeric (const const_IWSubstring & s, char op, double & zresult);

  public:
    IW_TDT_Filter_Base ();
    virtual ~IW_TDT_Filter_Base ();

    virtual int ok () const = 0;
    virtual int debug_print (std::ostream &) const = 0;

    virtual int construct (const const_IWSubstring &) = 0;

    int set_item_to_process (int s) { _item_to_process = s; return 1;}

//  Second argument is the return code if the item is absent from the TDT

    virtual int matches (const IW_TDT &, int) = 0;
};

class IW_TDT_Filter_lt : public IW_TDT_Filter_Base
{
  private:
    double _value;

  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    int matches (const IW_TDT &, int);
};

class IW_TDT_Filter_gt : public IW_TDT_Filter_Base
{
  private:
    double _value;

  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    int matches (const IW_TDT &, int);
};

class IW_TDT_Filter_eq : public IW_TDT_Filter_Base
{
  private:
    double _value;

  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    int matches (const IW_TDT &, int);
};

class IW_TDT_Filter_present : public IW_TDT_Filter_Base
{
  private:

  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    int matches (const IW_TDT &, int);
};

class IW_TDT_Filter_regexp : public IW_TDT_Filter_Base
{
  private:
    IW_Regular_Expression _regexp;

  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    int matches (const IW_TDT &, int);
};

class IW_TDT_Filter_random : public IW_TDT_Filter_Base
{
  private:
    float _ratio;

  public:

    int ok () const;
    int debug_print (std::ostream &) const;

    int construct (const const_IWSubstring &);

    int matches (const IW_TDT &, int);
};

#if defined(__INTEL_COMPILER)
class IW_TDT_Filter
{
  private:
  public:
    int active () const { return 0;}

    int build_from_string (const const_IWSubstring & s) { return 1;}

    int debug_print (std::ostream &) const { return 1;}
    
    int report (std::ostream & os) const { return 1;}

    int matches (const IW_TDT & tdt) { return 1;}
};
#else
class IW_TDT_Filter
{
  private:

    IW_Logical_Expression _logexp;

    resizable_array_p<IW_TDT_Filter_Base> _f;

    int _tdts_examined;

    int _tdts_passing;

    int _tdts_with_tag_present;     // not implemented

//  Mar 2004. I want the ability to say that if a tag is absent from the TDT
//  that constitutes a match

    int _tag_absent_means_match;

//  private functions

    void _default_values ();

  public:
    IW_TDT_Filter ();
    ~IW_TDT_Filter ();

    int ok () const;
    int debug_print (std::ostream &) const;

    int active () const { return _f.number_elements ();}

    int build_from_string (const const_IWSubstring & s);

    void set_tag_absent_means_match (int s) { _tag_absent_means_match = s;}

    int report (std::ostream &) const;

    int matches (const IW_TDT &);
};

#endif


extern int display_tdt_filter_syntax (std::ostream & os);

#endif
