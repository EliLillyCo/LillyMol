# Correlate

[![Build Status](https://travis-ci.org/iwatsonlilly/Correlate.jl.svg?branch=master)](https://travis-ci.org/iwatsonlilly/Correlate.jl)

[![Coverage Status](https://coveralls.io/repos/iwatsonlilly/Correlate.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/iwatsonlilly/Correlate.jl?branch=master)

[![codecov.io](http://codecov.io/github/iwatsonlilly/Correlate.jl/coverage.svg?branch=master)](http://codecov.io/github/iwatsonlilly/Correlate.jl?branch=master)

Enforces a given correlation structure on a matrix.

Especially during simulation studies, different sets of numbers may be
generated independently.  But there may be a known correlation
structure needed between these separate sets of numbers.  Correlate
shuffles the rows in order to produce a final dataset that exhibits
the requested correlation structure.

**Installation**: `Pkg.add("Correlate")`

**Usage**:
Correlate starts with a matrix of data lacking the desired correlation
structure.  Use `cor()` to get the existing correlation structure in a matrix.
`correlate!` is given a target correlation structure and a tolerance.

It works by exchanging column data between rows in ways that move the
correlation structure of the data closer to the target.  It can do
this efficiently by keeping track of partial results of the
correlation computation, and then incrementing and decrementing these
sums as a result of the changes. It never does a full matrix correlation
computation.

```julia
    correlate!(data::Array{Float64,2}, corm::Array{Float64,2}, tolm::Array{Float64,2};
               maxiter=100*nrows, report=0, maxsec=0, verbose=false)
```
CORM is the desired correlation matrix. While it should have 1.0 on the diagonal,
that is never checked. In fact, only the upper triangular elements are examined.

TOLM is a matrix (or scalar) of tolerances. The algorithm keeps iterating until the
difference between each correlation and the target is below tolerance. Again, only
the upper triangular part is examined.

It is possible to provide a matrix of correlation values that will be
impossible to achieve.  `correlate!` does not check this.  You
may wish to attempt a Cholesky decomposition, `chol(corm)`, to see if
a solution is possible.  Even though `correlate!` does not check the
diagonals in the correlation matrix, you will need them when
invoking `chol`.

By default, `correlate!` will attempt `100*nrows` iterations. You can control 
how many iterations with the `maxiter` keyword argument.  If it cannot
make any improvement in any `maxiter/10` consecutive attempts, it will give up.
If you are not satisfied with the result at that time, you could invoke it
again, but it may, or may not, be able to improve.

For large matrices, `correlate!` can take a long time to converge. 
You might specify a number of (wall clock) seconds for it to devote to
the problem before it abandons the attempt: the `maxsec` keyword.  In
all cases, a "partial" result will be in the matrix.

The `report=n` keyword argument instructs `correlate!` to provide
information periodically as it works - every `n` iterations.

When the `verbose` keyword argument is set to `true`, `correlate!` will report
results of the computation when done.

This report might look like

    Converged. Performed 41170 iterations, made 3200 switches
      1.0000  0.7990 -0.0009 -0.2008       1.0000  0.8000  0.0000 -0.2000
      0.7990  1.0000 -0.3990 -0.4990       0.8000  1.0000 -0.4000 -0.5000
     -0.0009 -0.3990  1.0000  0.1004       0.0000 -0.4000  1.0000  0.1000
     -0.2008 -0.4990  0.1004  1.0000      -0.2000 -0.5000  0.1000  1.0000
    Sum of errors 0.005128, 0 items not converged
    Switched column 1 1916 times
    Switched column 2 1036 times
    Switched column 3 130 times
    Switched column 4 118 times
    Calculation took 0.03 seconds
  
First the number of attempts (41170) and then the number of times data was
switched (the resulting arrangement was closer to the target).

Then, side by side, are the correlation matrix (from `cor`), and then the
target correlation matrix.

Lastly there is information about how many times data in each column was
switched. This can vary considerably depending on whether the input data
has an existing correlation structure, or how different the correlation
targets might be.

`correlate!` will return true if convergence is achieved, false otherwise.

Whereas correlate was built to enforce a desired correlation structure on a
set of data, it can also be used to remove an existing correlation structure
by submitting a correlation matrix with all off diagonals set to zero `ones(n)`.

## Standalone Tool

In the src directory you will find standalone Linux shell tool, correlate.jl. This provides
an interface to the Julia package via a command line tool.

Usage:

```
Enforces correlation structure in data (using Julia package Correlate)
correlate.jl <cor> <tol> -S <output> <input>
 -n <n>        maximum number of iterations to try (def 1000*nrows)
 -C <fname>    read correlation matrix from <fname>
 -c <cor>      two columns in input, just one correlation value needed
 -T <fname>    read tolerance matrix from <fname>
 -t <tol>      use a uniform tolerance for all pairs
 -x <seconds>  abandon computation after <seconds> elapsed time
 -i <sep>      input separator (default ,) use 'tab' or 'space'
 -nohdr        input file does NOT have a header record
 -r <n>        report progress and check convergence and timeout every <n> steps
 -B <n>        input consists of multiple chunks of data, <n> items in each chunk
 -S <fname>    output file name
 -U <suffix>   if processing multiple files, suffix for files created
 -v            verbose output
```

The only mandatory options are

1. Desired correlation
2. Tolerance
3. Name of output file
4. Name of input file

As an example, in the src directory try

    ./correlate.jl -v -C corm.dat -i space -t 0.01 rand.1000.3.dat > out.dat

You should see a message similar to what was displayed above.

**Maintenance**: Correlate is maintained by i.watson@lilly.com (aka ianiwatson@gmail.com)

