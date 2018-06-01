#ifndef KAHANSUM_J
#define KAHANSUM_J

class KahanSum
{
  private:
    double _sum;
    double _c;

  public:
    KahanSum ();

    KahanSum & operator += (double);
    KahanSum & operator += (const KahanSum &);
    KahanSum & operator = (double);

    operator double () const {return _sum;}

    double sum() const { return _sum;}
};

#endif
