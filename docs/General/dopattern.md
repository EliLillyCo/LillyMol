# dopattern

dopattern is one of the most versatile and useful tools included
with LillyMol.

Everything that dopattern does, could be done at the command
line, but generally with much more complexity.

For example, do write out the numbers from 1 to 5
```
dopattern -a 1 -o 5 echo %
```
The numeric indices to be scanned are specified via the -a and -o
options. For each value, the value replaces the % character.

This is the same as
```
for i in $(seq 5) ; do  echo $i ; done
```
but simpler.

Usually dopattern performs operations on series of files. This is
a common pattern of use when dealing with a large set of molecules
that may have been split into chunks, where the same operation needs
to be performed on each chunk.

For example, if we have split a large set of molecules into chunks
and need to do a substructure search.

For example, split 1M molecules into 100k chunks and search
each one.
```
iwsplit -suffix smi -n 100000 large.smi
dopattern -o 10 tsubstructure -s 'smarts' -m matches% iwsplit%.smi
```
For each value in the range [1,10] the % character is assigned that value, and
tsubstructure is run on the file iwsplit%.smi and matches are written to matches%.smi.

dopattern can also iterate over file names rather than over numbers. For example
if you wanted to process all files with a given suffix in a directory
```
dopattern -suffix smi tsubstructure -m %_matches %.smi
```
So if there is a file 'foo.smi' in the directory, one of the commands executed
would be
```
dopattern -suffix smi tsubstructure -m foo_matches foo.smi
```

Or you can specify the pattens that are iterated
```
dopattern -do foo,bar wc %.smi
```
will execute 'wc foo.smi' and 'wc bar.smi'.

You can also iterate through subdirectories
```
dopattern -subdir 'cd % && wc foo.smi'
```
will cycle through all subdirectories, cd into each one, and count
the file 'foo.smi' in that directory.

For complete flexibility, the tokens to be substituted can be put in
a file, one per line, and processed with the -file option.

LillyMol tool iwsplit is often used to split large smiles files into
smaller chunks, with sequences of names like 'iwsplit1.smi', 'iwsplit2.smi'.
These can be handled via
```
dopattern -stem iwsplit%.smi ls iwsplit%.smi
```
where it identies all the files that are part of that series in the current
directory. Note however that it does this via a sequential scan, upwards
from 1, so if there is a missing file, the sequence will stop there.

## Dry Runs
A tool like dopattern can be both powerful and dangerous. The dry run
'-dry' option does not do anything, just lists the commands that would
be issued without the -dry option.

## Errors
By default, a failed command stops all processing. The `-k` option
suppresses this, while error messages are issued. The `-ks` option
again ignores errors, but does so silently.

## Nested
It is common to have nested dopattern invocations. For example within
each directory, list the first five files in a sequence.
```
dopattern -subdir 'cd % && dopattern -f Q -o 5 ls Q.smi'
```
cd in to each subdirectory. Within each subdirectory, run dopattern,
this time using Q rather than %, and list the files '1.smi',
'2.smi'....

## Quotes
Redirecting stdout and stderr can become problematic, sometimes
just using single and double quotes can solve the problem. But sometimes
it can be too messy. dopattern supports the strings 'squote' and 'dquote'
which are expanded to single and double quote characters as commands
are just about to be executed.

## Parallel
dopattern can use gnu parallel to execute the tasks in parallel on
an SMP machine. In addition, if you have a Grid Engine cluster,
dopattern can be used to launch tasks to a Grid Engine cluster.

# parallel_process_file
This is a Go tool that can in some cases perform the same
functionality as dopattern, but without the need to split an input
file. For example, to search chembl for phenols across 8 CPU's
try
```
parallel_process_file -v -h 8 -cmd 'tsubstructure -s '[OH]-c' -m phenol%d -' chemb.smi
```
which runs in about 6 seconds. 

The 'magic' is the %d in the output file. This will be interpolated with
a unique sequence number, so 'phenol1.smi', 'phenol2.smi'... will be generated.
The other 'magic' is that the command just be reading from standard input. If you
run this with the -v option, you will see why that is the case.

To generate all the cheap descriptors on Chembl.
```
parallel_process_file -v -h 8 -cmd 'iwdescr.sh -O none - > chembl%d.w' chembl.smi
```
which takes about 9 seconds. This time output is to standard output, with each
job writing a uniquely named file, 'chembl1.w', 'chembl2.w' .... Concatenate
these with descriptor_file_cat - which is a cat utility that knows about
header records.
