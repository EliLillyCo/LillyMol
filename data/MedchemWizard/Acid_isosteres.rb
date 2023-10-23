#!/usr/bin/env ruby

ianhome = "/home/rx87851"

require "#{ianhome}/ruby/lib/iwcmdline.rb"
#require "#{ianhome}/ruby/lib/iwcmdline_v2.rb"
require "#{ianhome}/ruby/lib/bindir.rb"

$expert = false

def usage (rc)
  $stderr.print "Process the Simple H Replacement file\n"
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

# [*:1][C:2]([OH])=O>>[*:1][C:2]=1ON=C(C=1)O 
# [*:1][C:2]([OH])=O>>[*:1][C:2]=1C(=NOC=1)O 
# [*:1][C:2]([OH])=O>>[*:1][C:2]=1C(=NON=1)O 
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

ndx = 1
ARGF.each do |line|
  f = line.chomp.split
  next if 0 == f.size
  rxn = f.shift

  name = f.join('_')
  name.gsub!('"', '')
  name.gsub!(',', '-')

  j = rxn.split('>>')
  $stderr << "Reagents #{j[0]} products #{j[1]} name #{name}\n"

  if "" == name
    name = "#{ndx}"
    ndx += 1
  end

  fname = "Acid_isosteres.#{name}.rxn"
  outp = File.open(fname, "w")
  outp << "(0 Reaction\n"
  outp << "  (A C Comment \"Acid_isosteres.#{name}\")\n"
  outp << "  (0 Scaffold\n"
  if "[*:1][C:2]([OH])=O" == j[0] 
    outp << "    (A C smarts \"[CD3](=O)[OD1]\")\n"
  else
    $stderr << "Unrecognised lhs '#{j[0]}'\n"
    exit 1
  end
  outp << "    (A I remove_atom (1 2))\n"
  outp << "  )\n"

  if ! j[1].start_with?("[*:1]")
    $stderr << "Unrecognised RHS '#{j[1]}'\n"
    exit 1
  end

  outp << "  (1 Sidechain\n"
  outp << "    (A C smiles \"#{j[1][5,99]} AcidIso.#{name}\")\n"
  outp << "    (A I join (0 0))\n"
  outp << "  )\n"
  outp << ")\n"
end
