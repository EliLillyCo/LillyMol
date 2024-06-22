package main

import (
  "os"
  "unicode"   
  "io"
//"regexp"
//"path"
//"path/filepath"
  "strings"
//"os/user"
  "bytes"
  "bufio"
  "fmt"
//"runtime"
  "flag"
)

func usage (rc int) {
  fmt.Fprintln(os.Stderr, "Reverses a reaction file - changes products to reagents")
  fmt.Fprintln(os.Stderr, " -S <prefix>          produce new files with prefix <prefix>")
  fmt.Fprintln(os.Stderr, " -R <fname>           write reagent molecules to <fname>")
  fmt.Fprintln(os.Stderr, " -P <fname>           write product molecules to <fname>")
  fmt.Fprintln(os.Stderr, " -r                   write each reagent set to separate file")
  fmt.Fprintln(os.Stderr, " -p                   write each product set to separate file")
  fmt.Fprintln(os.Stderr, " -A                   write agents (default is discard)\n")
  fmt.Fprintln(os.Stderr, " -N                   file name is reaction name\n")
  fmt.Fprintln(os.Stderr, " -nspace              remove whitespace from names when using -N option\n")
  fmt.Fprintln(os.Stderr, " -ftn                 use the first token of the name when using the -N option\n")
  fmt.Fprintln(os.Stderr, " -v                   verbose output\n")
  os.Exit(rc)
}

var write_agents bool = false

var stream_for_reagents *os.File = nil
var stream_for_products *os.File = nil

var separate_file_each_reagent bool = false
var separate_file_each_product bool = false

var file_name_is_reagent_name bool = false

var remove_name_whitespace bool = false
var first_token_of_name bool = false

func remove_non_alpha(s []byte) {
  for i := 0; i < len(s); i++ {
    c := s[i];
    if unicode.IsLetter(rune(c)) || unicode.IsNumber(rune(c)) || '_' == c || '-' == c {
    } else {
      s[i] = '_'
    }
  }
}

func write_reaction(rxn bytes.Buffer,
                    prefix string,
                    ndx * int) bool {

//fmt.Fprintf(os.Stderr, "In write_reaction len %d\n", len(rxn.String()))
  f := strings.Split(rxn.String(), "\n")
//fmt.Fprintf(os.Stderr, "Split into %d\n", len(f))
  if len(f) < 5 {
    fmt.Fprintf(os.Stderr, "Reaction too short %d first '%s'\n", len(f), f[0])
    return false
  }
  count_record := f[4]
//fmt.Fprintf(os.Stderr, "Count '%s'\n", count_record)
  var nr int
  var np int
  na := 0
  c,err := fmt.Sscanf(count_record, "%3d%3d%3d", &nr,&np,&na)
  if 3 == c && nil == err {
  } else {
    c,err := fmt.Sscanf(count_record, "%3d%3d", &nr,&np)
    if 2 == c && nil == err {
    } else {
      fmt.Fprintf(os.Stderr, "Invalid counts record '%s'\n", strings.Trim(count_record, "\n"))
      return false
    }
    na = 0
  }

  tot := nr + np + na + 1
  offsets := make([]int, tot)

  j := 0
  for i,line := range(f) {
    if "$MOL" == line {
      offsets[j] = i
//    fmt.Fprintf(os.Stderr, "  offset %d at %d\n", j, i)
      j++
    } else if "M  END" == line {
      offsets[j] = i+1         // relying on the fact that J got incremented above
    }
  }

  *ndx++
  var fname string
  if file_name_is_reagent_name {
    tmp := strings.TrimSpace(f[1])
    if first_token_of_name {
      tmp = strings.Split(tmp, " ")[0]
    } else if remove_name_whitespace {
      b := []byte(tmp)
      remove_non_alpha(b)
      tmp = string(b)
//    tmp = strings.Replace(tmp, " ", "_", -1)
    }
    fname = fmt.Sprintf("%s.rxn", tmp)
  } else {
    fname = fmt.Sprintf("%s%d.rxn", prefix, *ndx)
  }

//fmt.Fprintf(os.Stderr, "%d reagents %d products %d agents file '%s'\n", nr, np, na, fname)

  outp,err := os.Create(fname)
  if nil != err {
    fmt.Fprintf(os.Stderr, "Cannot create '%s'\n", fname)
    return false
  }

  for i := 0; i < 4; i++ {
    fmt.Fprintln(outp, f[i])
  }

  if ! write_agents {
    na = 0
  }

  fmt.Fprintf(outp, "%3d%3d%3d\n", np, nr, na)
//fmt.Fprintf(os.Stderr, "%3d%3d%3d\n", np, nr, na)
//fmt.Fprintf(os.Stderr, "tot %d, len(offset) %d\n", tot, len(offsets))

  for i := offsets[nr]; i < offsets[nr+np]; i++ {
    fmt.Fprintln(outp, f[i])
  }
  for i := 5; i < offsets[nr]; i++ {
    fmt.Fprintln(outp, f[i])
  }
  if write_agents {
    for i := offsets[nr + np]; i < len(f); i++ {
      fmt.Fprintln(outp, f[i])
    }
  }

  if nil != stream_for_products && len(stream_for_products.Name()) > 0 {
    if np > 1 {
      fmt.Fprintf(os.Stderr, "Multiple product molecules %d '%s'\n", np, strings.Trim(f[offsets[nr]+1], "\n"))
    }
    for i := 6; i < offsets[nr]; i++ {
      fmt.Fprintln(stream_for_products, f[i])
    }
    fmt.Fprintln(stream_for_products, "$$$$")
  } else if separate_file_each_product {
    fname := fmt.Sprintf("P%s%d.sdf", prefix, *ndx)
    outp,err := os.Create(fname)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Cannot open product file '%s' '%v'\n", fname, err)
      return false
    }
    for i := 6; i < offsets[nr]; i++ {
      if ("$MOL" == f[i]) {
        fmt.Fprintln(outp, "$$$$")
      } else {
        fmt.Fprintln(outp, f[i])
      }
    }
    fmt.Fprintln(outp, "$$$$")
    outp.Close()
  }

  if nil != stream_for_reagents && len(stream_for_reagents.Name()) > 0 {
    if nr > 1 {
      fmt.Fprintf(os.Stderr, "Multiple reagent molecules %d '%s'\n", np, strings.Trim(f[offsets[nr]+1], "\n"))
      nr = 1
    }
    for i := offsets[nr]; i < offsets[nr + np]; i++ {
      if "MOL$" == f[i] {
        fmt.Fprintln(stream_for_reagents, "$$$$")
      } else {
        fmt.Fprintln(stream_for_reagents, f[i])
      }
    }
    fmt.Fprintln(stream_for_reagents, "$$$$")
  } else if separate_file_each_reagent {  
    fname := fmt.Sprintf("R%s%d.sdf", prefix, *ndx)
    outp,err := os.Create(fname)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Cannot open product file '%s' '%v'\n", fname, err)
      return false
    }
    for i := offsets[nr] + 1; i < offsets[nr + np]; i++ {
      if ("$MOL" == f[i]) {
        fmt.Fprintln(outp, "$$$$")
      } else {
        fmt.Fprintln(outp, f[i])
      }
    }
    fmt.Fprintln(outp, "$$$$")
    outp.Close()
  }

  return true
}

func rxn_reverse_2(scanner *bufio.Reader,
                   prefix string,
                   ndx * int) bool {
  var rxn bytes.Buffer

  var lines_read = 0

  for ; ; {
    line, err := scanner.ReadString('\n');
    lines_read++
    if io.EOF == err {
      break
    }

//  fmt.Fprintf(os.Stderr, "Processing '%s'\n", line)
    if ("$RXN\n" == line) {
      if 4 == lines_read {
        rxn.Reset()
      } else if rxn.Len() > 0 {
        write_reaction(rxn, prefix, ndx)
        rxn.Reset()
      }
      rxn.WriteString(line)
    } else {
      rxn.WriteString(line)
    }
  }

  if rxn.Len() > 0 {
    write_reaction(rxn, prefix, ndx)
  }

  return true
}

func rxn_reverse(input_fname string, prefix string, ndx * int) bool {
  input,err := os.Open(input_fname)
  if nil != err {
    fmt.Fprintf(os.Stderr, "Cannot open file '%s' '%v'\n", input_fname, err)
    return false
  }

  scanner := bufio.NewReader(input)

  return rxn_reverse_2(scanner, prefix, ndx)
}

func main() {
  var verbose bool
  var prefix string
  var reagent_fname string
  var product_fname string

  flag.BoolVar  (&write_agents,      "A",      false, "write agents")
  flag.StringVar (&prefix,           "S",      "", "create files with prefix <s>")
  flag.StringVar (&reagent_fname,    "R",      "", "file name for reagents")
  flag.StringVar (&product_fname,    "P",      "", "file name for products")
  flag.BoolVar   (&separate_file_each_reagent, "r", false, "separate file each product")
  flag.BoolVar   (&separate_file_each_product, "p", false, "separate file each product")
  flag.BoolVar   (&file_name_is_reagent_name,  "N", false, "file name is reaction name")
  flag.BoolVar   (&first_token_of_name,        "ftn", false, "first token of reaction name as fname")
  flag.BoolVar   (&remove_name_whitespace,     "nspace", false, "remove whitespace from rxn name for fname")
  flag.BoolVar  (&verbose,           "v",      false, "verbose output")

  flag.Parse()

  nfiles := len(flag.Args())
  if 0 == nfiles {
    fmt.Fprintln(os.Stderr, "Must specify file(s) to process")
    usage(1)
  }

  if len(product_fname) > 0 {
    var err error
    stream_for_products,err = os.Create(product_fname)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Cannot open product file name '%s' %v'\n", product_fname, err)
      os.Exit(1)
    }
  }

  ndx := 0
  for i,fname := range(flag.Args()) {
    if ! rxn_reverse(fname, prefix, &ndx) {
      fmt.Fprintln(os.Stderr, "Failure processing '%s'\n", fname)
      os.Exit(i+1)
    }
  }

  if verbose {
    fmt.Fprintf(os.Stderr, "Processed %d reactions\n", ndx)
  }
}
