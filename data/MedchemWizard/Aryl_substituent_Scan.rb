#!/usr/bin/env ruby

# [c;H1:1]>>[c:1]F add F
# [c;H1:1]>>[c:1]Cl add Cl
# [c;H1:1]>>[c:1]Br add Br

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

ARGF.each do |line|
  f = line.chomp.split
  next if 0 == f.size
  rxn = f.shift

  name = f.join('_')
  name.gsub!('"', '')
  name.gsub(',', '-')

  j = rxn.split('>>')
  $stderr << "Reagents #{j[0]} products #{j[1]} name #{name}\n"

  fname = "Aryl_substituent_Scan.#{name}.rxn"
  outp = File.open(fname, "w")
  outp << "(0 Reaction\n"
  outp << "  (A C Comment \"Simple_H_Replacement.#{name}\")\n"
  outp << "  (0 Scaffold\n"
  if "[c;H1:1]" == j[0] 
    outp << "    (A C smarts \"[cH1]\")\n"
  else
    $stderr << "Unrecognised lhs '#{j[0]}'\n"
    exit 1
  end
  outp << "  )\n"

  if ! j[1].start_with?("[c:1]")
    $stderr << "RHS not aromatic '#{j[1]}'\n"
    exit 1
  end

  outp << "  (1 Sidechain\n"
  outp << "    (A C smiles \"#{j[1][5,99]} ArylSS.#{name}\")\n"
  outp << "    (A I join (0 0))\n"
  outp << "  )\n"
  outp << ")\n"
end
