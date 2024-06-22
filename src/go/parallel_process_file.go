package main

import (
	"flag"
	"fmt"
	"os"
	"os/exec"
	"runtime"
        "strings"
	"path"
)

func usage(rc int) {
	fmt.Fprintf(os.Stderr, "\n%s [-h <n>] [-s <n>] [-v] -cmd 'command to run' input_file\n\n",
		path.Base(os.Args[0]))
	fmt.Fprintln(os.Stderr,
		"Runs any command in multiple processes on a single file")
	fmt.Fprintln(os.Stderr,
		"Input is assumed to be a SMILES file 'SMILES id exact_mass'")
	fmt.Fprintln(os.Stderr, " -h <n>            number of threads to use (must be at least 2)")
	fmt.Fprintln(os.Stderr, " -cmd <cmd>        command to run - must assume input from stdin")
	fmt.Fprintln(os.Stderr, " -s <n>            header records to skip in input file")
	fmt.Fprintln(os.Stderr, " -v                verbose output\n")
	os.Exit(rc)
}

func next_beginning_of_line(inp *os.File, o int64) int64 {
	b := make([]byte, 512)

	rc := o

	for {
		n, err := inp.Read(b)

		if nil != err {
			fmt.Fprintf(os.Stderr,
				"next_beginning_of_line:error reading file %v\n",
				err)
			return 0
		}

		for j := 0; j < n; j++ {
			if '\n' == b[j] {
				return rc + int64(j) + 1
			}
		}

		rc += int64(len(b))
	}
}

func run_a_chunk_dd(fname string,
	cmd string,
	begin_offset int64,
	end_offset int64,
	waiter chan []byte,
        num uint,
	verbose bool) {

	mycmd := fmt.Sprintf("dd if=%s ibs=1 skip=%d count=%d status=none | %s",
		fname, begin_offset, (end_offset - begin_offset), cmd)

        if (strings.Contains(mycmd, "%d")) {
                mycmd = fmt.Sprintf(mycmd, num)
        }

	if verbose {
		fmt.Fprintf(os.Stderr, "Executing '%s'\n", mycmd)
	}

	out, err := exec.Command("/bin/bash", "-c", mycmd).Output()

	if err != nil {
		fmt.Fprintf(os.Stderr, "Cannot run '%s' %v\n", mycmd, err)
		os.Exit(1)
	}

	waiter <- out
}

func main() {
	var nthreads int
	var verbose bool
	var command_to_run string
	var header_records int

	flag.IntVar(&nthreads, "h", 2, "threads to use")
	flag.StringVar(&command_to_run, "cmd", "", "command to run")
	flag.IntVar(&header_records, "s", 0, "header records to skip")
	flag.BoolVar(&verbose, "v", false, "verbose output")

	flag.Parse()

	if 0 == len(command_to_run) {
		fmt.Fprintf(os.Stderr,
			"Must specify command via the -cmd option\n")
		usage(1)
	}

	nfiles := len(flag.Args())

	if 0 == nfiles {
		fmt.Fprintln(os.Stderr, "Must specify file(s) to process")
		usage(1)
	}

	if nfiles > 1 {
		fmt.Fprintln(os.Stderr,
			"Sorry, only processes one file at a time")
		usage(1)
	}

	if 1 == nthreads {
		fmt.Fprintf(os.Stderr, "only runs parallel\n")
		os.Exit(1)
	}

	fname := flag.Args()[0]

	st, err := os.Stat(fname)

	if nil != err {
		fmt.Fprintf(os.Stderr, "Missing file '%s' stat failed %v\n",
			fname, err)
		os.Exit(1)
	}

	file_size := st.Size()

	inp, err := os.Open(fname)

	if err != nil {
		fmt.Fprintf(os.Stderr, "Cannot open %s %v\n", fname, err)
		os.Exit(1)
	}

	bytes_per_thread := file_size / int64(nthreads)

	if verbose {
		fmt.Fprintf(os.Stderr,
			"File size %v into %d chunks, %v bytes per chunk\n",
			file_size, nthreads, bytes_per_thread)
	}

	if nthreads > 0 {
		runtime.GOMAXPROCS(nthreads)
		if verbose {
			fmt.Fprintf(os.Stderr, "Running %d threads on %d CPU's\n",
				runtime.GOMAXPROCS(0), runtime.NumCPU())
		}
	}

	waiter := make(chan []byte)

	prev, _ := inp.Seek(0, 0) // which should be zero bytes

	for i := 0; i < nthreads; i++ {

		o, err := inp.Seek(bytes_per_thread*int64(i+1), 0)

		if nil != err {
			fmt.Fprintf(os.Stderr, "Cannot seek to start %d, %v\n",
				i, err)
			os.Exit(1)
		}

		if (nthreads - 1) == i {
			o = file_size
		} else {
			o = next_beginning_of_line(inp, o)
		}

		go run_a_chunk_dd(fname, command_to_run, prev, o, waiter,
                                  uint(i+1), verbose)

		prev = o
	}

	for i := 0; i < nthreads; i++ {
		s := <-waiter

		if 0 == i || 0 == header_records {
			os.Stdout.Write(s)
		} else {
			for j := 0; j < len(s); j++ {
				if '\n' != s[j] {
					continue
				}
				os.Stdout.Write(s[j+1:])
				break
			}
		}
	}

	os.Exit(0)
}
