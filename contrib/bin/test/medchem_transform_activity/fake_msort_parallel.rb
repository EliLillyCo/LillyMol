#!/usr/bin/env ruby
# Deliberately tiny SMILES fixture sorter; production uses LillyMol msort_parallel.
raise 'wrong sort key' unless ARGV[0, 2] == ['-k', 'natoms']
raise 'missing thread flag' unless ARGV[2] == '-h'
raise 'missing input type' unless ARGV[4, 2] == ['-i', 'smi']
lines = File.readlines(ARGV.fetch(6))
threads = Integer(ARGV.fetch(3))
raise 'invalid thread count' unless threads.positive? && threads <= lines.size
File.write(ENV['MCTA_SORT_LOG'], ARGV.join("\n")) if ENV['MCTA_SORT_LOG']
exit 1 if ENV['MCTA_SORT_FAIL']
lines.pop if ENV['MCTA_SORT_TRUNCATE']
lines << lines.first if ENV['MCTA_SORT_EXTRA']
# These test fixtures contain only unbracketed, single-letter heavy atoms.
lines.sort_by { |line| line.split.first.length }.each { |line| puts line }
