#!/usr/bin/env ruby

ianhome = "/home/rx87851"

require "#{ianhome}/ruby/lib/iwcmdline.rb"
#require "#{ianhome}/ruby/lib/iwcmdline_v2.rb"
require "#{ianhome}/ruby/lib/bindir.rb"

$expert = false

def usage (rc)
  $stderr.print "What does this programme do?\n"
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

# IWCmdlineV2  if using v2

cl = IWCmdline.new("-v-expert-bindir=dir-keep")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

bindir = Bindir.new(ianhome)

if cl.option_present('bindir')
  cl.values('bindir').each do |d|
    bindir.add_dir(d)
  end
end

if 0 == ARGV.size
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

header = ARGF.gets.split("\t")

fname = Array.new
output = Array.new

istop = 0

(0..header.length).step(2) do |i|
  fname[i] = header[i].gsub(" ", "_") << ".txt"
  
  next if /Reference/.match(fname[i])
  fname[i].gsub!(/\//, "-")
  $stderr << "Col #{i} is #{header[i]}\n"

  output[i] = File.open(fname[i], "w")
  istop = i + 1
end


ARGF.each do |line|
  f = line.chomp.split("\t")
  (0..istop).step(2) do |i|
    rxn = f[i]
    output[i] << f[i] << ' ' << f[i+1] << "\n"
  end
end
