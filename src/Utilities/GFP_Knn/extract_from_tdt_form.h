#ifndef IW_EXTRTDTFORM_H
#define IW_EXTRTDTFORM_H

class IWString;
class const_IWSubstring;

extern int extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       IWString & zresult);

extern int extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       int & zresult);
extern int extract_from_tdt_form (const const_IWSubstring & buffer,
                       const IWString & tag,
                       float & zresult);
#endif
