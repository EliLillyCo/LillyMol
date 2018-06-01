#ifndef IW_MOPAC_H                     
#define IW_MOPAC_H                     

class Mopac_Output_Control_Object 
{
  private:
    atom_number_t _start_atom;
    int           _native_ordering;
    IWString        _keywords;
    IWString        _comments;

  public:
    Mopac_Output_Control_Object (atom_number_t = INVALID_ATOM_NUMBER);

    int            add_keyword (const char *);
    void           set_keywords (const char * c) { _keywords = c;}
    const IWString & keywords () { return _keywords;}

    const IWString & comments () { return _comments;}
    void           set_comment (const char * c) {_comments = c;}

    atom_number_t start_atom () const { return _start_atom;}
    int native_ordering () const { return _native_ordering;}

};

#endif
