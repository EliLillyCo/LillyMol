#ifndef NNORM_H
#define NNORM_H

class NNorm
{
  private:
    int _n;

  public:
    NNorm();

    double operator() (float, float) const;

    double final() (float, int, int) const;
};

#endif
