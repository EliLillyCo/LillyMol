package main

import (
  "os"
//"os/exec"   
//"io"
//"regexp"
//"path"
//"path/filepath"
  "strings"
//"time"
//"os/user"
  "bufio"
  "fmt"
  "bytes"
//"io/ioutil"
//"runtime"
  "flag"
)

func usage (rc int) {
  fmt.Fprintln(os.Stderr, "Converts reaction smiles to smiles")
  fmt.Fprintln(os.Stderr, "by default, converts all reagents to a single multi-fragment molecule")
  fmt.Fprintln(os.Stderr, " -rev                 reverse the reaction. THIS IS DONE FIRST! - lhs is then rhs")
  fmt.Fprintln(os.Stderr, " -ronly               write only the reagents")
  fmt.Fprintln(os.Stderr, " -ponly               write only the products")
  fmt.Fprintln(os.Stderr, " -all                 keep all reagents (default is to write only the first")
  fmt.Fprintln(os.Stderr, " -rsm                 concatenate all reagents into a single disconnected molecule")
  fmt.Fprintln(os.Stderr, " -tid                 trim to first token of id")
  fmt.Fprintln(os.Stderr, " -plus2dot            change + signs on either side to . (for Marvin compatibility)")
  fmt.Fprintln(os.Stderr, " -M                   append Chemaxon fragment grouping directives |f...|")
  fmt.Fprintln(os.Stderr, " -rmagent             remove agents from output")
  fmt.Fprintln(os.Stderr, " -v                   verbose output")
  os.Exit(rc)
}

type RXNSMILES_Parms struct {
  trim_id_to_first_token bool
  ronly bool
  ponly bool
  reverse bool
  plus2dot bool
  keep_all_reagents bool
  compress_all_reagents_to_one bool
  write_agents bool
  marvin_frag bool
}



func do_plus_to_dot(s string,
                    parm RXNSMILES_Parms) string {
  var buffer bytes.Buffer

  var f bytes.Buffer
  f.WriteString("|f")

  component_number := 0
  component_grouping_active := false
  in_agent := false

  i := 0    // need to scope here because we use it below the loop

  in_square_bracket := false
  for ; i < len(s); i++ {
    c := s[i]
    if '[' == c {
      in_square_bracket = true
      buffer.WriteString("[")
    } else if ']' == c {
      in_square_bracket = false
      buffer.WriteString("]")
    } else if '.' == c {
      buffer.WriteString(".")
      if in_agent {
      } else if component_grouping_active {
        f.WriteString(fmt.Sprintf(".%d", component_number+1))     // extend current group
      } else {                         // start a new component group
        if len(f.String()) > 2 {          // need a comma if we are not the first
          f.WriteString(",")
        }
        f.WriteString(fmt.Sprintf("%d.%d", component_number, component_number+1))    // start a new one
        component_grouping_active = true
      }
      component_number++
    } else if '+' == c && ! in_square_bracket {
      buffer.WriteString(".")
      component_grouping_active = false
      component_number++
    } else if '>' == c {
      buffer.WriteString(">")
      if ! in_agent {     // first >
        in_agent = true
        component_grouping_active = false
        component_number++
      } else {         // the second one
        if '>' != s[i-1] {
          component_number++
        }
        in_agent = false
      }
    } else if ' ' == c {
      break
    } else {
      buffer.WriteString(string(c))
    }
  }

// we have consumed all the actual reaction

  if parm.marvin_frag && len(f.String()) > 2 {     // starts as '|f'
    buffer.WriteString(fmt.Sprintf(" %s|", f.String()))
  }

  return buffer.String()
}

func write_to_first_plus_sign(s string, outp *os.File) {
  var buffer bytes.Buffer

  in_square_bracket := false
  for i := 0; i < len(s); i++ {
    c := s[i]
    if '[' == c {
      in_square_bracket = true
      buffer.WriteString("[")
    } else if ']' == c {
      in_square_bracket = false
      buffer.WriteString("]")
    } else if '+' == c && ! in_square_bracket {
      break
    } else {
      buffer.WriteString(string(c))
    }
  }

  fmt.Fprintf(outp, "%s", buffer.String())
}

func write_id(f []string,
              trim_id_to_first_token bool,
              outp *os.File) {
  fmt.Fprintf(outp, " %s", f[1])    // always

  if trim_id_to_first_token {    // we are done
     return
  }

  for i := 2; i < len(f); i++ {
    fmt.Fprintf(outp, " %s", f[i])
  }

  return
}

func rxnsmiles2smi(inp *os.File,
                   parm RXNSMILES_Parms,
                   outp  *os.File) bool {
  scanner := bufio.NewScanner(inp)

  for scanner.Scan() {
    f1 := strings.Split(scanner.Text(), " ")     // separate into reaction and name

    rxn := f1[0]    // first token on line is the reaction

    if parm.plus2dot {
      rxn = do_plus_to_dot(rxn, parm)    // do the whole thing, even if we may not write parts of it
    }

    f2 := strings.Split(rxn, ">")

    if 3 != len(f2) {
      fmt.Fprintf(os.Stderr, "Invalid reaction '%s'\n", rxn)
      return false
    }

    if parm.reverse {        // important, this is done first
      f2[0],f2[2] = f2[2],f2[0]
    }

    if parm.ponly {       // only products being written
      fmt.Fprintf(outp, "%s", f2[2])
      write_id(f1, parm.trim_id_to_first_token, outp)
      fmt.Fprintf(outp, "\n")
      continue
    }

//  At this stage, the reagents will be written

//  fmt.Fprintf(os.Stderr, "keep_all_reagents %v f2[0] '%s'\n", parm.keep_all_reagents, f2[0])
//  fmt.Fprintf(os.Stderr, "compress_all_reagents_to_one %v\n", parm.compress_all_reagents_to_one)

    if parm.compress_all_reagents_to_one {
      if ! parm.plus2dot {     // if not already done
        f2[0] = do_plus_to_dot(f2[0], parm)
      }
      fmt.Fprintf(outp, "%s", f2[0])
    } else if parm.keep_all_reagents {
      fmt.Fprintf(outp, "%s", f2[0])
    } else {                                    /// hmmm, breaks if we have already converted
      write_to_first_plus_sign(f2[0], outp)
    }

    if parm.ronly {    // only reagents being written
      fmt.Fprintf(outp, "\n")
      continue
    }

    if parm.write_agents {
      fmt.Fprintf(outp, ">%s", f2[1])
    } else {
      fmt.Fprintf(outp, ">")
    }

    fmt.Fprintf(outp, ">%s", f2[2])

    write_id(f1, parm.trim_id_to_first_token, outp)

    fmt.Fprintf(outp, "\n")
  }

  return true
}
/*
func just_translate_plus_to_dot(inp *os.File, 
                                reverse bool,
                                outp *os.File) bool {
  scanner := bufio.NewScanner(inp)

  if ! reverse {
    for scanner.Scan() {
      do_plus_to_dot(scanner.Text(), outp)    // would accidentially translate a + in the ID
      fmt.Fprintf(outp, "\n")
    }
    return true
  }

  for scanner.Scan() {
    f1 := strings.Split(scanner.Text(), " ")    // separate into reaction and name

    rxn := f1[0]

    f2 := strings.Split(rxn, ">")

    f2[0],f2[2] = f2[2],f2[1]      // reverse them
    do_plus_to_dot(f2[0], outp)
    fmt.Fprintf(outp, ">>")
    do_plus_to_dot(f2[2], outp)
    write_id(f1, trim_id_to_first_token, outp)
    fmt.Fprintf(outp, "\n")
  }

  return true
}
*/

/*func do_translate_plus_to_dot(buffer string,
                              parm RXNSMILES_Parms) string {
  f := strings.Split(buffer, " ")

  s := do_plus_to_dot(f[0], parm)

  var data []byte
  data = append(data, []byte(s)...)
  for i := 1; i < len(f); i++ {
    data = append(data, []byte(f[i])...)
  }

  return string(data);
}*/

func main () {

  var keep_all_reagents bool
  var compress_all_reagents_to_one bool
  var reverse bool
  var rmagent bool
  var plus2dot bool
  var marvin_frag bool
  var verbose bool
  var trim_id_to_first_token bool
  var ronly bool
  var ponly bool

  flag.BoolVar  (&ronly,                  "ronly", false,  "Only write the reagent molecule(s)")
  flag.BoolVar  (&ponly,                  "ponly", false, "Only write the product molecule(s)")
  flag.BoolVar  (&reverse,                "rev",   false, "Reverse the reaction")
  flag.BoolVar  (&keep_all_reagents,      "all",   false, "Retain all reagents (default is only first)")
  flag.BoolVar  (&compress_all_reagents_to_one, "rsm", false, "Concatenate all reagents into single molecule")
  flag.BoolVar  (&rmagent,                "rmagent",   false, "Remove the agent")
  flag.BoolVar  (&plus2dot,               "plus2dot",  false, "Convert + characters to .")
  flag.BoolVar  (&marvin_frag,            "M",      false, "append Marvin fragment grouping")
  flag.BoolVar  (&trim_id_to_first_token, "tid",    false, "Truncate name to first token")
  flag.BoolVar  (&verbose,                "v",      false, "verbose output")

  flag.Parse()

  if 0 == len(flag.Args()) {
    fmt.Fprintf(os.Stderr, "Must specify input file\n")
    usage(1)
  }

  if ! reverse && ! ronly && ! ponly && ! plus2dot && ! rmagent && ! compress_all_reagents_to_one {
    ronly = true
    plus2dot = true
    keep_all_reagents = true
  }

  if compress_all_reagents_to_one {
    keep_all_reagents = true
  }

  var parms RXNSMILES_Parms
  parms.reverse = reverse
  parms.ronly = ronly
  parms.ponly = ponly
  parms.keep_all_reagents = keep_all_reagents
  parms.compress_all_reagents_to_one = compress_all_reagents_to_one
  parms.plus2dot = plus2dot
  parms.marvin_frag = marvin_frag
  parms.trim_id_to_first_token = trim_id_to_first_token
  parms.write_agents = ! rmagent

  fname := flag.Args()[0];

  rc := true

  if 1 == len(flag.Args()) && "-" == fname {
    rc = rxnsmiles2smi(os.Stdin, parms, os.Stdout)
  } else {
    for _,fname := range flag.Args() {
      inp,err := os.Open(fname)
      if nil != err {
        fmt.Fprintf(os.Stderr, "Cannot open '%s'\n", fname)
        rc = false
        break;
      }

      defer inp.Close()
      if ! rxnsmiles2smi(inp, parms, os.Stdout) {
        fmt.Fprintf(os.Stdout, "Error processing '%s'\n", fname)
        rc = false
        break
      }
    }
  }

  if rc {
    os.Exit(0)
  } else {
    os.Exit(1)
  }
}
