#!/usr/bin/env ruby
require 'open3'
exe, = ARGV

def check(exe, args, input, expected = nil, success: true)
  File.write('input.txt', input)
  stdout, stderr, status = Open3.capture3(exe, *args, 'input.txt')
  abort "#{args.inspect}: #{stderr}" unless status.success? == success
  abort "#{args.inspect}: expected #{expected.inspect}, got #{stdout.inspect}" if expected && stdout != expected
end

rows = "a x\na y\nb y\nc y\nd y\ne z\n"
check(exe, ['-column', '1', '-column', '2,n=2', '-D', 'duplicates.txt'], rows,
      "a x\nb y\nc y\ne z\n")
abort 'Wrong duplicate output' unless File.read('duplicates.txt') == "a y\nd y\n"
# Rejection by either column must not consume the other column's allowance.
check(exe, ['-column', '1', '-column', '2'], "a x\nb x\nb y\n", "a x\nb y\n")
check(exe, ['-column', '1', '-column', '2,n=2,subset=(.*)..$'],
      "a abc01\na abc02\nb abc03\nc abc04\nd def01\n", "a abc01\nb abc03\nd def01\n")
# Commas inside a regexp remain part of the regexp, not qualifier separators.
check(exe, ['-column', '1,n=1,subset=(.{1,3})..$'], "abc01\nabc02\ndef01\n",
      "abc01\ndef01\n")
check(exe, ['-column', '-1,n=2'], "a x\nb x\nc x\nd y\n", "a x\nb x\nd y\n")
check(exe, ['-column', '1', '-nc', '-z', '-h', '1'],
      "header\n001Ab\n01aB\n02B\n", "001Ab\n02B\n")
check(exe, ['-column', '1', '-j', '1'], "name\na\na\n", "name\na\n")
check(exe, ['-column', '2', '-csv', '-q'], "a,\"x,y\"\nb,\"x,y\"\nc,z\n",
      "a,\"x,y\"\nc,z\n")
# Counts persist across input files.
File.write('second.txt', "b y\nc z\n")
check(exe, ['-column', '2', 'second.txt'], "a y\nd w\n", "b y\nc z\nd w\n")
%w[0 1,n=0 1,n=x 1,n=2,n=3 1,unknown=3 1, 1,subset=abc 1,subset=( 1,,n=2].each do |spec|
  check(exe, ['-column', spec], "a x\n", success: false)
end
[['-o'], ['-O', 'counts.txt'], ['-c', '1'], ['-x', '1'], ['-r', '2'],
 ['-subset', '(.*)'], ['-p', '0.5'], ['-all']].each do |args|
  check(exe, ['-column', '1', *args], "a x\n", success: false)
end
check(exe, ['-column', '1', '-column', '1'], "a x\n", success: false)
check(exe, ['-column', '3'], "a x\n", success: false)
check(exe, ['-column', '1,subset=(.*)..$'], "x\n", success: false)
# Original composite-key and global-subset behaviour remain available.
check(exe, ['-c', '1', '-c', '2'], "a x\na y\na x\n", "a x\na y\n")
check(exe, ['-c', '2', '-subset', '(.*)..$'], "a abc01\nb abc02\nc def01\n",
      "a abc01\nc def01\n")
# Both paths use the same case/truncation/subset order. Leading delimiters
# remain intact, rather than collapsing distinct values into an empty key.
[['-c', '1', '-subset', '(.*)'], ['-column', '1,subset=(.*)']].each do |args|
  check(exe, [*args, '-nc', '-trunc', '.'],
        "ABC.one\nabc.two\n.def\n.ghi\n.def\n", "ABC.one\n.def\n.ghi\n")
end
# Zero stripping remains per token when constructing a legacy composite key.
check(exe, ['-c', '1', '-c', '2', '-z'], "001 002\n01 02\n01 03\n",
      "001 002\n01 03\n")
puts 'independent column checks passed'
