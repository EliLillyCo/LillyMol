#ifndef IW_SHOULD_MATCH_H
#define IW_SHOULD_MATCH_H

extern int should_match (const IWString &, const char *, const char *);
extern int should_match (const const_IWSubstring &, const char *, const char *);
extern int should_match (const char *, const char *, const char *);

#endif
