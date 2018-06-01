#ifndef IW_BF_H
#define IW_BF_H

/*
  Various convenient extensions of library functors

  The order of the variables in the comparison is very important
*/

template <typename T, typename C>
class IW_Binary_Comparison
{
  protected:
    T _v;
    C _c;

    typedef C Comparitor;

  public:
    IW_Binary_Comparison (const T & s) : _v (s) {};

    int operator () (const T & rhs) const { return _c (rhs, _v);}
};

template <typename T>
class IW_Less_Equal : public IW_Binary_Comparison<T, std::less_equal<T> >
{
  private:

  public:
    IW_Less_Equal (T f) : IW_Binary_Comparison<T, std::less_equal<T> > (f) {};
};

template <typename T>
class IW_Less_Than : public IW_Binary_Comparison<T, std::less<T> >
{
  private:

  public:
    IW_Less_Than (T f) : IW_Binary_Comparison<T, std::less<T> > (f) {};
};

template <typename T>
class IW_Greater_Than : public IW_Binary_Comparison<T, std::greater<T> >
{
  private:
  public:
    IW_Greater_Than (const T & f) : IW_Binary_Comparison<T, std::greater<T> > (f) {}
};

template <typename T>
class IW_Greater_Equal : public IW_Binary_Comparison<T, std::greater_equal<T> >
{
  private:
  public:
    IW_Greater_Equal (const T & f) : IW_Binary_Comparison<T, std::greater_equal<T> > (f) {}
};

template <typename T>
class IW_Not_Equal : public IW_Binary_Comparison<T, std::not_equal_to<T> >
{
  private:
    typedef std::not_equal_to<T> Comparitor;

  public:
    IW_Not_Equal (const T & f) : IW_Binary_Comparison<T, Comparitor> (f) {};
};

template <typename T>
class IW_Equal : public IW_Binary_Comparison<T, std::equal_to<T> >
{
  private:
  public:
    IW_Equal (const T & f) : IW_Binary_Comparison<T, std::equal_to<T> > (f) {};
};

class IW_Less_Equal_Int : public IW_Less_Equal<int>
{
  private:
  public:
    IW_Less_Equal_Int (int f) : IW_Less_Equal<int> (f) {};
};

class IW_Not_Equal_Int : public IW_Not_Equal<int>
{
  private:
  public:
    IW_Not_Equal_Int (int f) : IW_Not_Equal<int> (f) {};
};

class IW_Equal_Int : public IW_Equal<int>
{
  private:
  public:
    IW_Equal_Int (int f) : IW_Equal<int> (f) {};
};

class IW_Greater_Than_Int : public IW_Greater_Than<int>
{
  private:
  public:
    IW_Greater_Than_Int (int f) : IW_Greater_Than<int> (f) {};
};

#endif
