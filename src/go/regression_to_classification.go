package main
import (
  "os"
  "sort"
  "io/ioutil"
  "strings"
  "fmt"
  "strconv"
  "flag"
)

func usage (rc int) {
  fmt.Fprintln(os.Stderr, "Converts a continuous input to discrete outputs")
  fmt.Fprintln(os.Stderr, " -c <col>             column to process (default 2)")
  fmt.Fprintln(os.Stderr, " -C <c1,c2,c3>        class cutoffs (ascending order)")
  fmt.Fprintln(os.Stderr, " -L <l1,l2,l3>        class labels (same order as -C)")
  fmt.Fprintln(os.Stderr, " -min <frac>          ensure a minimum of <frac> items in the smallest class (2 class only)")
  fmt.Fprintln(os.Stderr, " -s <n>               number of header records to skip (0)")
  fmt.Fprintln(os.Stderr, " -r                   replace the input column (default is append)")
  fmt.Fprintln(os.Stderr, " -e                   append at end of record (default is after column)")
  fmt.Fprintln(os.Stderr, " -d <dname>           descriptor to process")
  fmt.Fprintln(os.Stderr, " -D <dname>           name for newly added column")
  fmt.Fprintln(os.Stderr, " -i <sep>             input separator (default ' '), 'tab' and 'space' recognised")

  os.Exit(rc)
}

type regression_to_classification_args struct {
  col int
  input_separator string
  nclasses int
  replace_input bool
  append_to_end bool
  prepend_before_col bool
  cutoff [] float64
  min_fraction float64
  label [] string
  hdr int
  col_name string
  dname string
  colname string
  verbose bool
}

func lowest_assigned_class(assigned []int, n int) (int,float64) {
  smallest_class_count := assigned[0]
  smallest_class := 0

  for i := 1; i < len(assigned); i++ {
    if assigned[i] < smallest_class {
      smallest_class_count = assigned[i]
      smallest_class = i
    }
  }

  return smallest_class, float64(smallest_class_count) / float64(n)
}

type Ndx_Value struct {
  ndx int
  v float64
}

type ByValue [] Ndx_Value

func (a ByValue) Len() int                { return len(a)}
func (a ByValue) Swap(i, j int)           { a[i], a[j] = a[j], a[i]}
func (a ByValue) Less(i, j int) bool      { return a[i].v > a[j].v}

func regression_to_classification(zdata []byte, r2c regression_to_classification_args,
                                  output *os.File) bool {
  col := r2c.col
  input_separator := r2c.input_separator
  cutoff := r2c.cutoff
  nclasses := r2c.nclasses
  verbose := r2c.verbose

  lines := strings.Split(string(zdata), "\n")

//fmt.Fprintf(os.Stderr, "Will process %d header records\n", r2c.hdr)

  for i := 0; i < r2c.hdr; i++ {
    line := string(lines[i])
    if 1 == len(line) {
      fmt.Fprintf(os.Stderr, "Cannot read header record'\n")
      return false;
    }

    f := strings.Split(strings.TrimRight(line,"\n"), input_separator)

    var app string

    for c := 0; c < len(f); c++ {
      if col < 0 && r2c.colname == f[c] {
        col = c
      }

      if c > 0 {
         fmt.Fprintf(output, "%s", input_separator)
      }

      if c != col {
        fmt.Fprintf(output, "%s", f[c])
      } else {
        if r2c.replace_input {
          fmt.Fprintf(output, "%s.C", f[c])
        } else if r2c.append_to_end {
          fmt.Fprintf(output, "%s", f[c])
          app = fmt.Sprintf("%s.C", f[c])
        } else if r2c.prepend_before_col {
          fmt.Fprintf(output, "%s.C%s%s", f[c], input_separator, r2c.dname)
        } else {
          fmt.Fprintf(output, "%s%s%s", f[c], input_separator, r2c.dname)
        }
      }
    }
    if r2c.append_to_end {
      fmt.Fprintf(output, "%s%s", input_separator, app)
    }
    fmt.Fprintf(output, "\n")
  }

  if col < 0 {
    fmt.Fprintf(os.Stderr, "Column '%s' not established\n", r2c.colname)
    return false
  }

  assigned := make([] int, nclasses + 1)    // overall statistics of how many assigned to each class

  individual_assignment := make([]string, 0)    // for each item in the input, what is the class assignment

  numeric_values := make([]Ndx_Value, 0)         // in case we need to re-assign some classes

// First scan through the input and determine the class assignments

  for i := r2c.hdr; i < len(lines) - 1; i++ {
    line := string(lines[i])
//  fmt.Fprintf(os.Stderr, "Processing %s\n", line)

    f := strings.Split(line, input_separator)

    if len(f) < col {
      fmt.Fprintf(os.Stderr, "Too few columns in input '%s'\n", line)
      return false
    }

    v,err := strconv.ParseFloat(strings.TrimRight(f[col],"\n"), 64)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Invalid numeric '%s'\n", f[col])
      return false
    }

    if r2c.min_fraction > 0.0 {
      x := Ndx_Value{i - r2c.hdr, v}
      numeric_values = append(numeric_values, x)
    }

    if v > cutoff[nclasses-1] {
      individual_assignment = append(individual_assignment, r2c.label[nclasses])
      assigned[nclasses] ++
    } else {
      for i := 0; i < nclasses; i++ {
        if v <= cutoff[i] {
          individual_assignment = append(individual_assignment, r2c.label[i])
          assigned[i] ++
          break;
        }
      }
    }
  }

  if verbose {
    for c := 0; c <= nclasses; c++ {
      fmt.Fprintf(os.Stderr, "%d items assigned %s %.2f\n", assigned[c], r2c.label[c], float64(assigned[c])/float64(len(individual_assignment)))
    }
  }

  if (r2c.min_fraction > 0.0) {
    if 2 != len(assigned) {
      fmt.Fprintf(os.Stderr, "Sorry, the minimum class fraction option only works with two class problems, found %d classes\n", len(assigned))
      return false
    }
    smallest_class, smallest_class_fraction := lowest_assigned_class(assigned, len(lines))
    if smallest_class_fraction < r2c.min_fraction {
      needed := int(r2c.min_fraction * float64(len(individual_assignment)) + 0.4999)
      fmt.Fprintf(os.Stderr, "Class %d only assigned %f times, need %f (%d items)\n", smallest_class, smallest_class_fraction, r2c.min_fraction, needed)
      sort.Sort(ByValue(numeric_values))
      to_class := individual_assignment[numeric_values[0].ndx]
      for i := 0; i < needed; i++ {
        j := numeric_values[i].ndx
//      fmt.Fprintf(os.Stderr, "i = %d %d %v current %s\n", i, j, numeric_values[i].v, individual_assignment[j])
        if individual_assignment[j] != to_class{
          individual_assignment[j] = to_class
        }
      }
    }
  }

  n := 0

  for i := r2c.hdr; i < len(lines) - 1; i, n = i + 1, n + 1 {
    line := string(lines[i])
//  fmt.Fprintf(os.Stderr, "Looking at %s\n", line)

    f := strings.Split(line, input_separator)

    var app string

    for i := 0; i < col; i++ {
      if i > 0 {
        fmt.Fprintf(output, "%s", input_separator)
      }
      fmt.Fprintf(output, "%s", f[i])
    }

    if r2c.prepend_before_col {
    } else if ! r2c.replace_input{
      fmt.Fprintf(output, "%s%s", input_separator, strings.TrimRight(f[col],"\n"))
    }

    if r2c.append_to_end {
      app = individual_assignment[n]
    } else {
      fmt.Fprintf(output, "%s%s", input_separator, individual_assignment[n])
    }

    if r2c.prepend_before_col {
      fmt.Fprintf(output, "%s%s", input_separator, strings.TrimRight(f[col], "\n"))
    }

    if col + 1 == len(f) {
      if len(app) > 0 {
        fmt.Fprintf(output, "%s%s", input_separator, app)
      }
      fmt.Fprintf(output, "\n")
    } else {
      for c := col+1; c < len(f); c++ {
        if len(app) > 0 && c == len(f) - 1 {
          fmt.Fprintf(output, "%s%s%s%s\n", input_separator, strings.TrimRight(f[c], "\n"), input_separator, app)
        } else {
          fmt.Fprintf(output, "%s%s", input_separator, f[c])
        }
      }
    }
  }

  return true
}


func main () {
  var str_cut string
  var labels string
  var input_separator string
  var col int
  var hdr int
  var dname string
  var colname string
  var replace_input bool
  var append_to_end bool
  var prepend_before_col bool
  var min_classf float64
  var verbose bool

  flag.StringVar (&str_cut,            "C",       "", "comma separated list of cut points")
  flag.StringVar (&labels,             "L",       "", "comma separated list of class labels")
  flag.IntVar    (&col,                "c",       2, "column")
  flag.IntVar    (&hdr,                "s",       0, "header records to skip")
  flag.Float64Var(&min_classf,         "min",     0.0, "minimum class fraction")
  flag.StringVar (&dname,              "D",       "R2C", "header for newly created column")
  flag.StringVar (&colname,            "d",       "", "descriptor to process")
  flag.StringVar (&input_separator,    "i",       " ",          "input separator (default ' '")
  flag.BoolVar   (&replace_input,      "r",       false, "replace input column")
  flag.BoolVar   (&append_to_end,      "e",       false, "append to end of input record")
  flag.BoolVar   (&prepend_before_col, "p",       false, "prepend class label before numeric column")
  flag.BoolVar   (&verbose,            "v",       false, "verbose output")

  flag.Parse()

  nfiles := len(flag.Args())
  if 0 == nfiles {
    fmt.Fprintln(os.Stderr, "Must specify file(s) to process")
    usage(1)
  }

  if nfiles > 1 {
    fmt.Fprintln(os.Stderr, "Not really designed for handling multiple files\n")
  }

  if 0 == len(str_cut) || 0 == len(labels) {
    fmt.Fprintln(os.Stderr, "Must specify both cut points (-C) and labels (-L)")
    usage(1)
  }

  if 1 != len(input_separator) {
    if "space" == input_separator {
      input_separator = " "
    } else if "tab" == input_separator {
      input_separator = "\t"
    } else if "comma" == input_separator {
      input_separator = ","
    } else {
      fmt.Fprintln(os.Stderr, "Unrecognised input delimiter specification '%s'", input_separator)
      usage(1)
    }
  }

  fc := strings.Split(str_cut, ",")
  fl := strings.Split(labels, ",")

  if len(fc) + 1 != len(fl) {
    fmt.Fprintf(os.Stderr, "Must specify N class cutoffs %d and N+1 labels %d\n", len(fc), len(fl))
    usage(1)
  }

  cutoff := make([]float64, len(fc))
  label := make([]string, len(fc)+1)

  for i := 0; i < len(fc); i++ {
    label[i] = fl[i]
    t, err := strconv.ParseFloat(fc[i], 64)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Invalid numeric '%s'\n", fc[i])
      os.Exit(1)
    }
    cutoff[i] = t
  }
  label[len(fc)] = fl[len(fc)]

  if 0 == hdr && len(colname) > 0 {
    hdr = 1
  }

  if len(colname) > 0 {
    col = -1
  }

  positions := 0
  if replace_input {
    positions += 1
  }
  if append_to_end {
    positions += 1
  }
  if prepend_before_col {
    positions += 1
  }
  if positions > 1 {
    fmt.Fprintf(os.Stderr, "Cannot both replace the current column and append\n")
    usage(1)
  }

  var r2c regression_to_classification_args
  r2c.col = col - 1
  r2c.replace_input = replace_input
  r2c.append_to_end = append_to_end
  r2c.prepend_before_col = prepend_before_col
  r2c.input_separator = input_separator
  r2c.nclasses = len(fc)
  r2c.min_fraction = min_classf
  r2c.cutoff = cutoff
  r2c.label = label
  r2c.hdr = hdr
  r2c.dname = dname
  r2c.colname = colname
  r2c.verbose = verbose

  var rc bool
  for _,fname := range flag.Args() {
    if "-" == fname {
      zdata,err := ioutil.ReadAll(os.Stdin)
      if nil != err {
        fmt.Fprintf(os.Stderr, "Cannot read stdin %v\n", err)
        os.Exit(1)
      }

      rc = regression_to_classification(zdata, r2c, os.Stdout)
    } else {
      zdata,err := ioutil.ReadFile(fname)
      if nil != err {
        fmt.Fprintf(os.Stderr, "Cannot read '%s' %v\n", fname, err)
        os.Exit(1)
      }
      rc = regression_to_classification(zdata, r2c, os.Stdout)
    }
  }

  if ! rc {
    os.Exit(1)
  }

  os.Exit(0)
}

