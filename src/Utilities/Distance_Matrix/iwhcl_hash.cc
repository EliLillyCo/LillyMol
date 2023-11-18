#include <stdlib.h>
#include <iostream>

#include "iwhcl_hash.h"

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

Lots of problems with this, finish it sometime
Leaks memory, I'm not sure of the data structure, etc....
Rubbish right now...


int
iwhcl_hash (const unsigned int n,
       unsigned int multiplier,
       const int iopt,
       int * ia,
       int * ib,
       float * crit,
       float * membr,
       unsigned int * nn,
       unsigned char * disnn,
       int * flag,
       hash_map<Pair_of_Uints, unsigned char, Pair_of_Uints_Hash_Fn> & diss)
{
  typedef hash_map<Pair_of_Uints, unsigned char, Pair_of_Uints_Hash_Fn> MAP;

  int jm, im, jj;
  unsigned int i2, j2;

  Pair_of_Uints ind(multiplier);
  Pair_of_Uints ind1(multiplier);
  Pair_of_Uints ind2(multiplier);
  Pair_of_Uints ind3(multiplier);
  Pair_of_Uints ndx(multiplier);

  unsigned char dmin;

  float diss_ind1, diss_ind2, x, xx, tmp, new_value;   // scope here for efficiency

  int ncl = n;

  for (unsigned int i = 0; i < n; i++)
  {
    membr[i] = static_cast<float>(1.0);
    flag[i] = 1;
  }

// Carry out an agglomeration - first create list of NNs

  unsigned char inf = numeric_limits<unsigned char>::max();

  for (unsigned int i = 0; i < n; i++)
  {
    disnn[i] = inf;
    nn[i] = n + n;    // something out of range
  }

  unsigned int n1, n2;
  unsigned char d;
  unsigned char longest_distance_present = static_cast<unsigned char>(0);

  for (MAP::const_iterator i = diss.begin(); i != diss.end(); i++)
  {
    n1 = (*i).first.n1();
    n2 = (*i).first.n2();

    d = (*i).second;

    if (d > longest_distance_present)
      longest_distance_present = d;

    if (d < disnn[n1])
    {
      disnn[n1] = d;
      nn[n1] = n2;
    }

    if (d < disnn[n2])
    {
      disnn[n2] = d;
      nn[n2] = n1;
    }
  }

  float longest_distance_present_float = static_cast<float>(longest_distance_present);

  label_400:

// Next, determine least diss. using list of NNs

  dmin=inf;
  for (unsigned i = 0; i < n; i++)
  {
    if (0 == flag[i])
      continue;

    if (disnn[i] < dmin)
    {
      dmin = disnn[i];
      im = i;
      jm = nn[i];
    }
  }

//cerr << "Closest pair " << im << " to " << jm << " dist " << static_cast<int>(dmin) << endl;

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

//cerr << "Grouping " << i2 << " (" << membr[i2] << ") with " << j2 << " (" << membr[j2] << ")\n";
  ia[n-ncl]=i2;
  ib[n-ncl]=j2;  
  crit[n-ncl]=dmin;

  ind3.set(i2, j2);     // ind3 will always be set in the hash

  ncl=ncl-1;
  cerr << "ncl = " << ncl << ", hash size " << diss.size() << endl;

  if (1 == ncl)
    return 1;

// Update dissimilarities from new cluster.

  flag[j2] = 0;
  dmin=inf;
  for (unsigned int k = 0; k < n; k++)
  {
    if (0 == flag[k])
      continue;

    if (k == i2)
      continue;

    if (i2 < k)
      ind1.set(i2, k);
    else
      ind1.set(k, i2);

    MAP::const_iterator f1 = diss.find(ind1);

    if (f1 == diss.end())
      diss_ind1 = longest_distance_present_float;
    else
      diss_ind1 = (*f1).second;

    if (j2 < k)
      ind2.set(j2, k);
    else
      ind2.set(k, j2);

    MAP::const_iterator f2 = diss.find(ind2);
    if (f2 == diss.end())
      diss_ind2 = longest_distance_present_float;
    else
      diss_ind2 = (*f2).second;

    if (iopt == 1)           // ward's minimum variance method - iopt=1.
    {
      xx = static_cast<float>(diss[ind3]);

      tmp=(membr[i2]+membr[k])*diss_ind1+
          (membr[j2]+membr[k])*diss_ind2-
              membr[k]*xx;
      x=membr[i2]+membr[j2]+membr[k];
      new_value = static_cast<unsigned char>(tmp/x);
    }

    else if (iopt == 2)            // single link method - iopt=2.
      new_value =iwmin(diss_ind1,diss_ind2);

    else if (iopt == 3)          // complete link method - iopt=3.
      new_value =iwmax(diss_ind1,diss_ind2);

    else if (iopt == 4)           // average link (or group average) method - iopt=4.
    {
      new_value = static_cast<unsigned char>((membr[i2]*diss_ind1+membr[j2]*diss_ind2)/
                        (membr[i2]+membr[j2]));
    }

    else if (iopt == 5)          // mcquitty's method - iopt=5.
    {
      new_value  = static_cast<unsigned char>((diss_ind1+diss_ind2)*0.5);
    }

    else if (iopt == 6)          // median (gower's) method - iopt=6.
    {
      xx = static_cast<float>(diss[ind3]);

      new_value = static_cast<unsigned char>(0.5*diss_ind1+0.5*diss_ind2-0.25*xx);
    }

    else if (iopt == 7)          // centroid method - iopt=7.
    {
      xx = static_cast<float>(diss[ind3]);

      new_value = static_cast<unsigned char>((membr[i2]*diss_ind1+membr[j2]*diss_ind2-
                membr[i2]*membr[j2]*xx/(membr[i2]+membr[j2]))/
                (membr[i2]+membr[j2]));
    }

    diss[ind1] = static_cast<unsigned char>(new_value + 0.4999);

    if (new_value < dmin)
    {
      dmin = static_cast<unsigned char>(new_value + 0.4999);
      jj=k;
    }
  }

  membr[i2]=membr[i2]+membr[j2];
  disnn[i2]=dmin;
  if (jj < 0)
    cerr << "No nearest neighbour found when joining " << i2 << " and " << j2 << endl;
  else
    nn[i2]=jj;

// update list of nns insofar as this is required.

  for (unsigned int i = 0; i < n - 1; i++)
  {
    if (0 == flag[i])
      continue;

    if (nn[i] == i2 || nn[i] == j2)    // redetermine nn of I
    {
      dmin=inf;
      jj = -1;
      for (unsigned int j = 0; j < n; j++)
      {
        if (j == i || 0 == flag[j])
          continue;

        if (i < j)
          ind.set(i, j);
        else
          ind.set(j, i);

        MAP::const_iterator f = diss.find(ind);
        if (f == diss.end())
          continue;

        if ((*f).second < dmin)
        {
          dmin=(*f).second;
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


size_t
IW64bithash::operator () (const iwuint64 & s) const
{
  unsigned int t[2];

  memcpy(t, &s, sizeof(s));

  return (t[0] + 1) * t[1];
}

template unsigned char iwmin(unsigned char, unsigned char);
template unsigned char iwmax(unsigned char, unsigned char);

int
Pair_of_Uints::operator== (const Pair_of_Uints & rhs) const
{
  if (_n1 != rhs._n1)
    return 0;

  return _n2 == rhs._n2;
}

size_t
Pair_of_Uints_Hash_Fn::operator() (const Pair_of_Uints & p) const
{
  return p.n1() * p.multiplier() + p.n2();
}
