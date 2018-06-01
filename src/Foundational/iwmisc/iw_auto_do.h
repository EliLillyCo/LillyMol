#ifndef IW_AUTO_DO_H
#define IW_AUTO_DO_H

/*
  Extending the concept of auto_ptr. It does something in its descructor.
  Why not have something that does something arbitrary on exit
*/

template <typename T, typename C>
class IW_Auto_Do
{
  private:
    T & _t;
    C & _op;

  public:
    IW_Auto_Do (T & t, C & op) : _t (t), _op (op) {}
    ~IW_Auto_Do () { _op (_t);}
};

/*
  A common usage for these is incrementing a variable
*/

template <typename T>
class IWIncrement
{
  private:
    T _inc;

  public:
    IWIncrement (const T & i) : _inc (i) {};

    void operator () (T & i) { i = i + _inc;}
};

/*
  Of course we can create a new type that automatically does increments.
  But, this will be less efficient because we will need to create an
  IWIncrement<T> object within the loop every time
*/

template <typename T, typename I>
class IW_Auto_Increment : public IW_Auto_Do<T, IWIncrement<T> >
{
  private:
    IWIncrement<T> _inc;

  public:
    IW_Auto_Increment (T & t, const T & i) : IW_Auto_Do<T, IWIncrement<T> > (t, _inc), _inc (i) {};

    void operator () (T & i) { _inc (i);}
};

#endif
