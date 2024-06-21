package main

import (
  "os"
  "unicode"   
  "io/ioutil"
//"io"
  "path"
//"bytes"
  "strings"
  "fmt"
//"runtime"
//"strconv"
  "flag"
)

func usage (rc int) {
  fmt.Fprintln(os.Stderr, "changes file names to get rid of spaces and special characters")
  fmt.Fprintln(os.Stderr, " -t <char>            replacement character (default _)")
  fmt.Fprintln(os.Stderr, " -ok <xxx>            specify characters that will NOT be translated")
  fmt.Fprintln(os.Stderr, " -dry                 dry run, do not actually change anything")
  fmt.Fprintln(os.Stderr, " -v                   verbose output")
  os.Exit(rc)
}

func already_created(formed map[string]bool, s string) bool {
  _,ok := formed[s]
  return ok
}

func file_exists(fname string) bool {
  _,err := os.Stat(fname)
  if nil != err {
    return false
  } else {
    return true
  }
}

func main () {

  var substitution_character string
  var no_translate string
  var alnum bool
  var help bool
  var dryrun bool
  var verbose bool

  flag.StringVar (&substitution_character,  "t", "_", "replacement character (default underscore _)")
  flag.StringVar (&no_translate,     "ok",      "",        "specify characters that will not be translated")
  flag.BoolVar   (&alnum,            "alnum",   false,     "alphanumeric characters only")
  flag.BoolVar   (&dryrun,              "dry",   false,     "dry run")
  flag.BoolVar   (&help,             "help",     false,     "display help message")
  flag.BoolVar   (&verbose,          "v",     false, "verbose output")

  flag.Parse()

  if 0 == len(flag.Args()) {
    fmt.Fprintln(os.Stderr, "Must specify directory to process")
    usage(1)
  }

  if help {
    usage(1)
  }

  if len(substitution_character) > 1 {
    fmt.Fprintf(os.Stderr, "Substitution character must be length 1 '%s' invalid\n", substitution_character)
    os.Exit(1)
  }

  formed := make(map[string]bool)

  rc := 0
  files_examined := 0;

  for _,dir_to_process := range flag.Args() {

    files,ok := ioutil.ReadDir(dir_to_process)
    if nil != ok {
      fmt.Fprintf(os.Stderr, "Cannot open directory %s\n", dir_to_process)
      os.Exit(1)
    }

    if verbose {
      fmt.Fprintf(os.Stderr, "Processing '%s'\n", dir_to_process)
    }

    for _,f := range files {
      s := f.Name()
      if ".." == s || "." == s {
        continue
      }

      files_examined++;

      initial_name := s
//    fmt.Fprintf(os.Stderr, "Initial name '%s'\n", initial_name)
      b := []byte(initial_name)

      for i := 0; i < len(b); i++ {
        c := b[i];
        if unicode.IsLetter(rune(c)) || unicode.IsNumber(rune(c)) || '_' == c || '.' == c || '-' == c {
        } else if substitution_character[0] == c {
        } else if strings.IndexByte(no_translate, c) >= 0 {
        } else if ' ' == c {
          b[i] = substitution_character[0]
        } else if ! alnum {
          b[i] = substitution_character[0]
        } else {
          b[i] = substitution_character[0]
        }
      }

      s = string(b)
      if initial_name == s {
        continue
      }

      f1 := path.Join(dir_to_process, initial_name)
      f2 := path.Join(dir_to_process, s)

      rc++;

      if already_created(formed, f2) {
        fmt.Fprintf(os.Stderr, "Cannot 'mv %s %s' because target already made. Skipped\n", f1, f2)
      } else if file_exists(f2) {
        fmt.Fprintf(os.Stderr, "Cannot 'mv %s %s' because target file already exists. Skipped\n", f1, f2)
      } else if dryrun {
        fmt.Fprintf(os.Stderr, "mv %s %s\n", f1, f2)
        formed[f2] = true
      } else {
        err := os.Rename(f1, f2)
        if nil != err {
          fmt.Fprintf(os.Stderr, "mv %s to %s failed %v\n", f1, f2, err)
          os.Exit(1)
        }
        formed[f2] = true
      }
    }
  }

  if verbose {
    fmt.Fprintf(os.Stderr, "Scanned %d files, ", files_examined)
    if dryrun {
      fmt.Fprintf(os.Stderr, "would change")
    } else {
      fmt.Fprintf(os.Stderr, "changed")
    }
    fmt.Fprintf(os.Stderr, " %d\n", rc)
  }

  os.Exit(0)
}
