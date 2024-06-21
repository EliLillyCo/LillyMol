package main
import (
  "fmt"
  "regexp"
  "strings"
  "os"
  "bufio"
  "flag"
  )

func usage (rc int) {
  fmt.Fprintln(os.Stderr, "grep on the text of a .sdf file")
  fmt.Fprintln(os.Stderr, " -rx <regexp>         regular expression search (default is string match)")
  fmt.Fprintln(os.Stderr, " -tdt                 works witn a .tdt (gfp) file, separator is |\n")
  fmt.Fprintln(os.Stderr, " -v                   verbose output")

  os.Exit(rc)
}

var entity_separator string

func grep_sdf_regexp_match(pattern *regexp.Regexp,
                           reader *bufio.Reader,
                           output *os.File) int {

  sdf := make([] byte, 4096);

  rc := 0;

  for {
    line, err := reader.ReadBytes('\n')
    if nil != err {
      break
    }

    sdf = append(sdf, line...)

//  if "$$$$\n" == string(line) {
    if entity_separator == string(line) {
      found := pattern.Find(sdf)
      if nil != found {
        output.Write(sdf)
//      sdf.WriteTo(output)
        rc += 1
      }
      sdf = nil;
    }
  }
  return rc
}

func grep_sdf_string_match( pattern string,
                           reader *bufio.Reader,
                           output *os.File) int {

  sdf := make([] byte, 4096);

  rc := 0;

  for {
    line, err := reader.ReadBytes('\n')
    if nil != err {
      break
    }

    sdf = append(sdf, line...)

    if entity_separator == string(line) {
      if strings.Index(string(sdf), pattern) >= 0 {
        output.Write(sdf)
//      sdf.WriteTo(output)
        rc += 1
      }
      sdf = nil;
    }
  }
  return rc
}

func main () {
  var rxstring string
  var tdt bool
  var verbose bool
  flag.StringVar   (&rxstring,      "rx",     "",    "regular expression search")
  flag.BoolVar     (&tdt,           "tdt",    false, "processing a .tdt (gfp) file")
  flag.BoolVar     (&verbose,       "v",      false, "verbose output")

  flag.Parse()

  nfiles := len(flag.Args())
  if 0 == nfiles {
    fmt.Fprintln(os.Stderr, "Must specify file(s) to process")
    usage(1)
  }

  var pattern string

  if tdt {
    entity_separator = "|\n"
  } else {
    entity_separator = "$$$$\n"
  }

  istart := 0


  var rx *regexp.Regexp = nil

  if len(rxstring) > 0 {
    var err error
    rx,err = regexp.Compile(rxstring)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Cannot compile '%s'\n", rxstring)
      os.Exit(1)
    }
  } else {
    if 1 == nfiles {
      fmt.Fprintln(os.Stderr, "Must specify pattern and file(s) to process")
      usage(1)
    }

    pattern = flag.Args()[0]
    istart = 1
  }

  files_with_matches := 0;
  matches_found := 0

  for i := istart; i < len(flag.Args()); i++ {
    fname := flag.Args()[i]
    file, err := os.Open(fname)
    if err != nil {
      fmt.Fprintf(os.Stderr, "Cannot open %s %v\n", fname, err)
      os.Exit(1)
    }
    defer file.Close()

    reader := bufio.NewReader(file)
    var tmp int;
    if len(pattern) > 0 {
      tmp = grep_sdf_string_match(pattern, reader, os.Stdout)
    } else {
      tmp = grep_sdf_regexp_match(rx, reader, os.Stdout)
    }
    if tmp > 0 {
      files_with_matches += 1
      matches_found += tmp;
    }
  }

  if 0 == matches_found {
    os.Exit(1)
  }

  if verbose {
    fmt.Fprintf(os.Stderr, "%d of %d files matched. Found %d occurrences\n", files_with_matches, (nfiles-1), matches_found)
  }

  os.Exit(0)
}
