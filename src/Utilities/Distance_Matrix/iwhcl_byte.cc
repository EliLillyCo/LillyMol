/*
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Parameters:                                               C
C                                                            C
C  N                 number of items
C  DISS(LEN)         dissimilarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, DISNN  vectors of length N, used to store      C 
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the dissimilarity assoc. C
C                    with the latter.                        C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C                                                            C
C------------------------------------------------------------C
      SUBROUTINE HC(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,DISNN,
     X                FLAG,DISS)
      implicit none
      integer n
      integer len
      integer iopt
      REAL MEMBR(N),DISS(LEN)
      INTEGER IA(N),IB(N)
      REAL CRIT(N)
      INTEGER NN(n)
      REAL DISNN(n)
      LOGICAL FLAG(N)
*/


/*
  Map row I and column J of upper half diagonal symmetric matrix 
  onto vector.
*/

#include <stdlib.h>
#include <limits>
using namespace std;

static inline unsigned int
ioffset(unsigned int N,
        unsigned int I,
        unsigned int J)
{
      return J+(I-1)*N-(I*(I+1))/2;
}

template <typename T>
T
iwmin (T f1,
       T f2)
{
  if (f1 < f2)
    return f1;

  return f2;
}

template <typename T>
T
iwmax (T f1,
       T f2)
{
  if (f1 > f2)
    return f1;

  return f2;
}

/*
  Special purpose version of clustering tool to deal with
  byte distances. The
*/

int
iwhcl_nn_byte (unsigned int n,
       unsigned int len,
       int iopt,
       int * ia,
       int * ib,
       float * crit,
       float * membr,
       int * nn,
       unsigned char * disnn,
       int * flag,
       unsigned char * diss)
{

  unsigned int i, j, jm, im, i2, j2, k, jj;
  int ind1, ind2, ind3;
  int ncl;
  int ind;
  float x, xx;

  unsigned char dmin;

  float tmp;

  unsigned char inf = numeric_limits<unsigned char>::max();

  for (i = 0; i < n; i++)
  {
    membr[i] = static_cast<float>(1.0);
    flag[i] = 1;
  }

  ncl = n;

// Carry out an agglomeration - first create list of NNs

  for (i = 0; i < n - 1; i++)
  {
    dmin = inf;
    for (j = i + 1; j < n; j++)
    {
      ind = ioffset(n, i, j);
      if (diss[ind] < dmin)
      {
        dmin = diss[ind];
        jm = j;
      }
    }

    nn[i] = jm;
    disnn[i] = dmin;
   }

   label_400:

// Next, determine least diss. using list of NNs

  dmin=inf;
  for (i = 0; i < n - 1; i++)
  {
    if (0 == flag[i])
      ;
    else if (disnn[i] < dmin)
    {
      im = i;
      jm = nn[i];
    }
  }

  ncl=ncl-1;

//This allows an agglomeration to be carried out.

  if (im < jm)
  {
    i2 = im;
    j2 = jm;
  }
  else
  {
    i2 = jm;
    j2 = im;
  }
  ia[n-ncl]=i2;
  ib[n-ncl]=j2;
  crit[n-ncl]=dmin;
  ind3=ioffset(n,i2,j2);

// Update dissimilarities from new cluster.

  flag[j2] = 0;
  dmin=inf;
  for (k = 0; k < n; k++)
  {
    if (0 == flag[k])
      continue;

    if (k == i2)
      continue;

     if (i2 < k)
       ind1=ioffset(n,i2,k);
     else
       ind1=ioffset(n,k,i2);

     if (j2 < k)
       ind2=ioffset(n,j2,k);
     else
       ind2=ioffset(n,k,j2);

     xx = static_cast<float>(diss[ind3]);

// ward's minimum variance method - iopt=1.

    if (iopt == 1)
    {
      tmp = (membr[i2]+membr[k])*diss[ind1]+
            (membr[j2]+membr[k])*diss[ind2]-
             membr[k]*xx;
      x = membr[i2] + membr[j2] + membr[k];
      diss[ind1]=static_cast<unsigned int>(tmp/x);
    }

// single link method - iopt=2.

   else if (iopt == 2)
     diss[ind1]=iwmin(diss[ind1],diss[ind2]);

// complete link method - iopt=3.

   else if (iopt == 3)
     diss[ind1]=iwmax(diss[ind1],diss[ind2]);

// average link (or group average) method - iopt=4.

   else if (iopt == 4)
     diss[ind1]=(membr[i2]*diss[ind1]+membr[j2]*diss[ind2])/
                       (membr[i2]+membr[j2]);

// mcquitty's method - iopt=5.

    else if (iopt == 5)
      diss[ind1] = (diss[ind1] + diss[ind2]) * 0.5;

// median (gower's) method - iopt=6.

    else if (iopt == 6)
      diss[ind1] = 0.5*diss[ind1] + 0.5*diss[ind2] - 0.25*xx;

// centroid method - iopt=7.

    else if (iopt == 7)
      diss[ind1]=(membr[i2]*diss[ind1]+membr[j2]*diss[ind2]-
                membr[i2]*membr[j2]*xx/(membr[i2]+membr[j2]))/
                (membr[i2]+membr[j2]);

    if (i2 > k)
      ;
    else if (diss[ind1] >= dmin)
      ;
    else
    {
      dmin=diss[ind1];
      jj=k;
    }
  }

  membr[i2]=membr[i2]+membr[j2];
  disnn[i2]=dmin;
  nn[i2]=jj;

// update list of nns insofar as this is required.

  for (i = 0; i < n - 1; i++)
  {
    if (0 == flag[i])
      continue;

    if (nn[i] == i2 || nn[2] == j2)
    {
//       (redetermine nn of i:)
      dmin=inf;
      for (j = i + 1; j < n; j++)
      {
        if (0 == flag[j])
          continue;

        ind=ioffset(n,i,j);
        if (diss[ind] < dmin)
        {
          dmin=diss[ind];
          jj=j;
        }
      }
      nn[i]=jj;
      disnn[i]=dmin;
    }
  }

// repeat previous steps until n-1 agglomerations carried out.

 if (ncl > 1)
   goto label_400;

  return 1;
}

