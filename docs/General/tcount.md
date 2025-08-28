# tcount
`tcount` is like `wc` except that it also counts the number of tokens
on each line. This can be particularly helpful when dealing with tabular
files which should be tabular. Debugging problems in such files can
be difficult.

Usage 
```
Usage: tcount <options> <fname1> <fname2> ...
  -n <number>       process only the first <number> records from each file
  -e <offset>       start processing <offset> bytes before EOF
  -r                report different token count from previous record
  -f                report different token count from first record
  -R <name>         write reported   records to <name>
  -U <name>         write unreported records to <name>
  -t <count>        write records with <count> tokens to -T <fname>
  -T <fname>        file name for corresponding -t
  -o <char>         when writing -T files, use <char> as output delimiter
  -C <fname>        write a descriptor file of 'line.number tcount'
  -b                abort on different token count from previous record
  -w                write reported records to stderr
  -d <delim>        delimiter is any character, or a token word(e.g. 'tab'). see '-d help' for a list (default ' ')
  -c <delim>        record delimiter (default \n, use ^m for DOS)
  -u                consecutive delimiters mean empty words in between
  -g                report non rectangular files
  -q                data is quoted tokens
  -k                data may contain backslashed token separators
  -v                verbose output
```
The `-r` option can be particularly useful, since rather than just accumulating
statics on how many rows have a how many tokens, it will report any occurrence
of differing token counts from one line to the next.

For reasons that are not clear, tcount can also be faster than wc on large files.
I think it might be related to shell variable LANG, but not sure.

Like most LillyMol tools it cannot deal with complex quoted tabular files.
