#include <stdlib.h>
//#include <stdio.h>

#include "All.h"

extern int pmcio_(char *ttypin, char *spew, ftnlen ttypin_len, ftnlen spew_len);
extern int error0_(void);
extern int audit_(char *, ftnlen);
extern int dbinit_(char *type, char *name, ftnlen type_len, ftnlen name_len);

/*int  SMIFILE  = 0;
int  CLOGP    = 1;
int  tabout   = 0;
int  outlev   = 0;*/

int
initialise_BioByte_Arrays()
{
/* smiles variables */
  if (! (      ncon = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    atnumb = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    charge = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    chiral = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    hcount = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (      pmap = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (       con = (integer *)malloc( 2048*sizeof(integer))) ) { exit(1); }
  if (! (      bond = (integer *)malloc( 2048*sizeof(integer))) ) { exit(1); }
  if (! (     dbond = (integer *)malloc( 2048*sizeof(integer))) ) { exit(1); }
  if (! (     vbond = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    atomwt = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (       cid = (integer *)malloc(  512*sizeof(integer))) ) { exit(1); }
  if (! (    aromat = (logical *)malloc(  256*sizeof(logical))) ) { exit(1); }
  if (! (    hknown = (logical *)malloc(  256*sizeof(logical))) ) { exit(1); }
  if (! (    multat = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    valenc = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (     cylen = (integer *)malloc(   32*sizeof(integer))) ) { exit(1); }
  if (! (     cycle = (integer *)malloc( 2048*sizeof(integer))) ) { exit(1); }
  if (! (      cmap = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    cyarom = (logical *)malloc(   32*sizeof(logical))) ) { exit(1); }
  if (! (    cymemb = (integer *)malloc( 8192*sizeof(integer))) ) { exit(1); }
  if (! (      orig = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    redist = (integer *)malloc( 1024*sizeof(integer))) ) { exit(1); }
  if (! (    cybond = (integer *)malloc( 2048*sizeof(integer))) ) { exit(1); }
  if (! (    zclass = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    zorder = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    uniqat = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }

/* thor variables */
if (! (      mtag = (   char *)malloc( 2000*sizeof(char))) ) { exit(1); }
if (! (     mtext = (   char *)malloc( 6400*sizeof(char))) ) { exit(1); }
if (! (    tfllen = (integer *)malloc(  512*sizeof(integer))) ) { exit(1); }
if (! (    bucket = (   char *)malloc(10240*sizeof(char))) ) { exit(1); }
if (! (    record = (   char *)malloc( 2048*sizeof(char))) ) { exit(1); }
if (! (       fmt = (   char *)malloc( 2400*sizeof(char))) ) { exit(1); }
if (! (      atag = (   char *)malloc(  300*sizeof(char))) ) { exit(1); }
if (! (      vtag = (   char *)malloc( 1200*sizeof(char))) ) { exit(1); }
if (! (    picked = (logical *)malloc(   30*sizeof(logical))) ) { exit(1); }
if (! (    bfname = (   char *)malloc(  320*sizeof(char))) ) { exit(1); }
if (! (    bfspec = (   char *)malloc(  320*sizeof(char))) ) { exit(1); }

/* clogp variables */
  if (! (     order = (integer *)malloc( 1050*sizeof(integer))) ) { exit(1); }
  if (! (     elist = (integer *)malloc( 1050*sizeof(integer))) ) { exit(1); }
  if (! (   olist__ = (integer *)malloc( 1050*sizeof(integer))) ) { exit(1); }
  if (! (    branch = (   real *)malloc( 1050*sizeof(real))) ) { exit(1); }
  if (! (     value = (   real *)malloc( 2000*sizeof(real))) ) { exit(1); }
  if (! (    reliab = (integer *)malloc( 2000*sizeof(integer))) ) { exit(1); }
  if (! (    zwiter = (integer *)malloc( 2000*sizeof(integer))) ) { exit(1); }
  if (! (     elink = (integer *)malloc( 2000*sizeof(integer))) ) { exit(1); }
  if (! (       ori = (integer *)malloc(  600*sizeof(integer))) ) { exit(1); }
  if (! (     olink = (integer *)malloc(  600*sizeof(integer))) ) { exit(1); }
  if (! (    oclass = (integer *)malloc(  600*sizeof(integer))) ) { exit(1); }
  if (! (       rho = (   real *)malloc(  600*sizeof(real))) ) { exit(1); }
  if (! (     sigma = (   real *)malloc(  600*sizeof(real))) ) { exit(1); }
  if (! (     gname = (integer *)malloc( 1050*sizeof(integer))) ) { exit(1); }
  if (! (    englsh = (integer *)malloc( 1050*sizeof(integer))) ) { exit(1); }
  if (! (     first = (logical *)malloc( 2600*sizeof(logical))) ) { exit(1); }
  if (! (     slink = (integer *)malloc( 2600*sizeof(integer))) ) { exit(1); }
  if (! (   const__ = (   real *)malloc(  300*sizeof(real))) ) { exit(1); }
  if (! (     ormat = (   real *)malloc( 2025*sizeof(real))) ) { exit(1); }
  if (! (    prxtyp = (   char *)malloc( 2100*sizeof(char))) ) { exit(1); }
  if (! (    ohcont = (   char *)malloc( 1050*sizeof(char))) ) { exit(1); }
  if (! (      nage = (   char *)malloc( 1050*sizeof(char))) ) { exit(1); }
  if (! (       env = (   char *)malloc(30000*sizeof(char))) ) { exit(1); }
  if (! (    stpool = (   char *)malloc(52000*sizeof(char))) ) { exit(1); }
  if (! (    cnstid = (   char *)malloc( 6000*sizeof(char))) ) { exit(1); }

/* fragment variables */
  if (! (      imap = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (      fmap = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    fatmap = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (    infmap = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }
  if (! (     fvals = (   real *)malloc(  256*sizeof(real))) ) { exit(1); }
  if (! (    fsmile = (   char *)malloc( 6400*sizeof(char))) ) { exit(1); }
  if (! (    gsmile = (   char *)malloc( 6400*sizeof(char))) ) { exit(1); }
  if (! (    fenvir = (   char *)malloc( 1200*sizeof(char))) ) { exit(1); }
  if (! (    isotyp = (   char *)malloc(   16*sizeof(char))) ) { exit(1); }
  if (! (   frgcfrg = (integer *)malloc( 2048*sizeof(integer))) ) { exit(1); }
  if (! (   flagfcf = (integer *)malloc(   80*sizeof(integer))) ) { exit(1); }
  if (! (  naromsmi = (integer *)malloc(  256*sizeof(integer))) ) { exit(1); }

/* explain variables */
if (! (      xtag = (   char *)malloc(  120*sizeof(char))) ) { exit(1); }
if (! (    xplain = (   char *)malloc( 9600*sizeof(char))) ) { exit(1); }

/* genie variables */
if (! (     wnode = (   char *)malloc(32000*sizeof(char))) ) { exit(1); }
if (! (     wedge = (   char *)malloc( 2000*sizeof(char))) ) { exit(1); }
if (! (    vector = (logical *)malloc(51456*sizeof(logical))) ) { exit(1); }
if (! (    vecset = (logical *)malloc(  201*sizeof(logical))) ) { exit(1); }
if (! (    vecmem = (logical *)malloc(  256*sizeof(logical))) ) { exit(1); }
if (! (    vecnam = (   char *)malloc(16080*sizeof(char))) ) { exit(1); }
if (! (    vectar = (   char *)malloc(32160*sizeof(char))) ) { exit(1); }
if (! (    hitlst = (integer *)malloc(50000*sizeof(integer))) ) { exit(1); }

/* init terminal, fragDB, audit control */
/* comment the line below for monolithic/standalone operation */
  pmcio_(" ", "QUIET", 1L, 5L);

  error0_();
  audit_("ON", 2L);

/* kludge to see fragdb statistics (need to load fragDB now) */
  dbinit_("LOAD", " ", 4L, 1L);

//  printf("size long int %d ftnlen %d real %d\n", sizeof(long int), sizeof(ftnlen), sizeof(real));

  return 1;
}

